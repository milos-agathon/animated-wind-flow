# 1) LOAD LIBRARIES
install.packages("pacman")
pacman::p_load(
    ecmwfr, terra, tidyverse, metR,
    classInt, scico, gganimate,
    lubridate, hms
)

# 2) DOWNLOAD ERA5 WIND DATA
my_key <- "****************" # PLEASE REGISTER AT CLIMATE DATA STORE (CDS) AND INSERT YOUR API KEY

ecmwfr::wf_set_key(key = my_key)

year <- 2025
month <- 4
day <- 24

date_time <- seq(
    from = lubridate::ymd_hms(
        paste(
            year, month, day, 00, 00, 00,
            sep = "-"
        )
    ),
    to = lubridate::ymd_hms(
        paste(
            year, month, day, 23, 00, 00,
            sep = "-"
        )
    ), by = "1 hours"
)

time <- hms::as_hms(date_time)

# -9.349365,34.474864,1.115112,38.642618

xmin <- -9.349365
ymin <- 34.474864
xmax <- 1.115112
ymax <- 38.642618
target <- "wind-gibraltar.nc"
path <- getwd()

ecmwfr::wf_request(
    request = list(
        dataset_short_name = "reanalysis-era5-single-levels",
        product_type = "reanalysis",
        format = "netcdf",
        variable = c(
            "10m_u_component_of_wind",
            "10m_v_component_of_wind"
        ),
        year = year, month = month, day = day,
        time = time,
        area = c(ymax, xmin, ymin, xmax),
        grid = c(.25, .25), # arond 5m
        target = target
    ), path = path
)

# 3) READ & TIDY THE DATA
wind_data <- terra::rast(target)
u10 <- terra::rast(target, subds = "u10")
v10 <- terra::rast(target, subds = "v10")

# Rename layers for pivoting
names(u10) <- paste0("u_", date_time)
names(v10) <- paste0("v_", date_time)

# Combine, pivot to long, then widen to get 
# 'u' & 'v' columns
wind_df <- c(u10, v10) %>%
    # 1) One‐shot dump: x, y, plus all u… & v… layers 
    # (assumes layer names carry the time)
    as.data.frame(xy = TRUE, na.rm = TRUE) %>%
    # 2) Pivot longer to get var = {u|v}, 
    # time = layer‐name, speed = value
    tidyr::pivot_longer(
        cols = dplyr::starts_with(
            c("u", "v")
        ), names_to = c("var", "time"),
        # use regex to capture 'u' or 'v'
        names_pattern = "^(u|v)_(.*)$",
        values_to = "speed"
    ) %>%
    # 3) Spread back out so each row has u & v
    tidyr::pivot_wider(
        names_from = var,
        values_from = speed
    ) %>%
    # 4) Parse the time‐string into POSIXct
    dplyr::mutate(
        time = lubridate::ymd_hms(
            time, tz = "UTC"
        )
    ) %>%
    # 5) Group by time & compute spatial means
    dplyr::group_by(x, y, time) %>%
    dplyr::summarise(
        u = mean(u, na.rm = TRUE),
        v = mean(v, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # 6) Return as tibble
    tibble::as_tibble()

head(wind_df)

# 4) PRE-COMPUTE STREAMLINES PER TIME
streamlines <- wind_df %>%
    split(.$time) %>%
    purrr::map_dfr(function(df_hour) {
        tm <- df_hour$time[1]
        p <- ggplot(
            df_hour, aes(
                x = x, y = y, dx = u, dy = v,
                color = sqrt(
                    after_stat(dx)^2 + after_stat(dy)^2
                ),
                alpha = sqrt(
                    after_stat(dx)^2 + after_stat(dy)^2
                ),
                linewidth = sqrt(
                    after_stat(dx)^2 + after_stat(dy)^2
                ),

            )
        ) +
        metR::stat_streamline(
            L = 7, res = 2, geom = "path"
        )
        ld <- layer_data(p)
        ld$speed <- sqrt(ld$dx^2 + ld$dy^2)
        ld$time <- tm
        ld
    })

# 5) DEFINE THE BREAKS
breaks <- classInt::classIntervals(
    streamlines$speed,
    n = 6, style = "equal"
)$brks

# 6) PLOT & ANIMATE
wind_plot <- ggplot(
    streamlines, aes(
        x = x, y = y, group = group,
        color = speed, linewidth = speed
    )
) +
# wind‐lines, thicker where faster
geom_path(
    lineend = "round", show.legend = TRUE
) +
# richer palette, reversed so high speeds pop
scico::scale_color_scico(
    palette = "batlow", direction = -1,
    name = "Wind speed (m/s)",
    breaks = breaks,
    labels = round(breaks, 1)
) +
# tie line size and transparency to speed
scale_alpha(
    range = c(.2, 1),
    breaks = breaks, guide = "none"
) +
scale_linewidth(
    range = c(.1, 1),
    breaks = breaks, guide = "none"
) +
coord_quickmap(expand = FALSE) +
labs(
    title = "Animated Wind Flow Over Gibraltar",
    subtitle = "ERA5 10m winds - April 24, 2025",
    caption = "Data: copernicus ERA5 | Author: Milos makes maps"
) +
theme_minimal() +
theme(
    plot.background = element_rect(
        fill = "white", color = NA
    ),
    panel.background = element_rect(
        fill = "white", color = NA
    ),
    panel.grid.major = element_line(
        size = .25, color = "#EDF6F9"
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    plot.title = element_text(
        size = 20, face = "bold", hjust = .5
    ),
    plot.subtitle = element_text(
        size = 14, hjust = .5, 
        margin = margin(b = 10)
    ),
    plot.caption = element_text(
        size = 8, hjust = 1
    ),
    legend.position.inside = c(.92, .25),
    legend.background = element_rect(
        fill = "white", color = NA
    ),
    legend.key = element_blank()
) +
# smooth, cubic tweening between time steps
gganimate::transition_time(time) +
gganimate::ease_aes("cubic-in-out")

# Render and save the GIF
wind_gif <- gganimate::animate(
    wind_plot,
    nframes = length(
        unique(streamlines$time)
    ) * 2,
    fps = 8, width = 600, height = 400,
    renderer = gganimate::gifski_renderer()
)

gganimate::anim_save(
    "wind_flow_gibraltar_hi.gif", wind_gif
)
