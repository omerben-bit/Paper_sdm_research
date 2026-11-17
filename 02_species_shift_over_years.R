# =============================================================================
# Same-season shifts for top 6 species: mean depth, lon, lat over years
# Dark-mode figures
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

dir.create("outputs_med_3d_pa", showWarnings = FALSE, recursive = TRUE)

# ---- dark theme -------------------------------------------------------------
theme_dark_map <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.background   = element_rect(fill = "#0d1117", color = NA),
      panel.background  = element_rect(fill = "#161b22", color = NA),
      panel.border      = element_rect(fill = NA, color = "#30363d"),
      panel.grid.major  = element_line(color = "#21262d", linewidth = 0.3),
      panel.grid.minor  = element_line(color = "#21262d", linewidth = 0.15),
      strip.background  = element_rect(fill = "#21262d", color = "#30363d"),
      strip.text        = element_text(color = "#e6edf3", face = "bold"),
      axis.text         = element_text(color = "#c9d1d9"),
      axis.title        = element_text(color = "#e6edf3"),
      axis.ticks        = element_line(color = "#30363d"),
      plot.title        = element_text(color = "#f0f6fc", face = "bold", hjust = 0),
      plot.subtitle     = element_text(color = "#8b949e", hjust = 0),
      legend.background = element_rect(fill = "#0d1117", color = NA),
      legend.key        = element_rect(fill = "#161b22", color = NA),
      legend.text       = element_text(color = "#c9d1d9"),
      legend.title      = element_text(color = "#e6edf3"),
      complete = TRUE
    )
}

source("gbif_paper.R")

gbif_paper <- load_gbif_paper_cleaned()
occ <- read_csv("data/processed/med_fish_occurrences_combined.csv", show_col_types = FALSE)

# ---- build working table with year + month (paper-cleaned GBIF only) --------
work <- gbif_paper %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(scientificName)) %>%
  mutate(
    eventDate_chr = as.character(eventDate),
    year = suppressWarnings(as.integer(substr(eventDate_chr, 1, 4))),
    month = suppressWarnings(as.integer(substr(eventDate_chr, 6, 7))),
    depth = as.numeric(depth)
  )

# Meteorological seasons (same temp regime within season)
# DJF = winter, MAM = spring, JJA = summer, SON = autumn
season_of <- function(m) {
  case_when(
    m %in% c(12, 1, 2) ~ "Winter (DJF)",
    m %in% c(3, 4, 5)  ~ "Spring (MAM)",
    m %in% c(6, 7, 8)  ~ "Summer (JJA)",
    m %in% c(9, 10, 11) ~ "Autumn (SON)",
    TRUE ~ NA_character_
  )
}

top6 <- occ %>%
  filter(occurrenceStatus == "present", !str_detect(scientificName, "^Rajidae$")) %>%
  count(scientificName, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(scientificName)

message("Top 6 species:")
print(top6)

dat <- work %>%
  filter(scientificName %in% top6) %>%
  filter(!is.na(year), year >= 1950, year <= as.integer(format(Sys.Date(), "%Y"))) %>%
  mutate(
    season = season_of(month),
    short_name = scientificName %>%
      str_replace(" \\(.*\\)$", "") %>%
      str_replace(" Linnaeus.*$", "") %>%
      str_replace(" Gmelin.*$", "") %>%
      str_replace(" Lesueur.*$", "") %>%
      str_replace(" Geoffroy.*$", "") %>%
      str_trim()
  ) %>%
  filter(!is.na(season))

message("Records by season:")
print(count(dat, season, sort = TRUE))

# Focus season = most records (usually best for Med GBIF); still plot all seasons
focus_season <- dat %>%
  count(season, sort = TRUE) %>%
  slice_head(n = 1) %>%
  pull(season)
message("Focus season (most data): ", focus_season)

# ---- means within year × season ---------------------------------------------
yearly_season <- dat %>%
  group_by(scientificName, short_name, season, year) %>%
  summarise(
    n = n(),
    mean_lon = mean(decimalLongitude, na.rm = TRUE),
    mean_lat = mean(decimalLatitude, na.rm = TRUE),
    mean_depth = mean(depth, na.rm = TRUE),
    n_depth = sum(!is.na(depth)),
    .groups = "drop"
  ) %>%
  mutate(mean_depth = ifelse(n_depth > 0, mean_depth, NA_real_))

write_csv(yearly_season, "outputs_med_3d_pa/shifts_yearly_means_by_season_top6.csv")

species_season_summary <- dat %>%
  group_by(scientificName, short_name, season) %>%
  summarise(
    n = n(),
    year_min = min(year, na.rm = TRUE),
    year_max = max(year, na.rm = TRUE),
    mean_lon = mean(decimalLongitude, na.rm = TRUE),
    mean_lat = mean(decimalLatitude, na.rm = TRUE),
    mean_depth = mean(depth, na.rm = TRUE),
    n_with_depth = sum(!is.na(depth)),
    .groups = "drop"
  ) %>%
  arrange(season, desc(n))
write_csv(species_season_summary, "outputs_med_3d_pa/shifts_species_summary_by_season_top6.csv")
print(species_season_summary)

focus <- yearly_season %>% filter(season == focus_season)

# ---- Plot 1: lon/lat tracks within focus season -----------------------------
p_track <- ggplot(focus, aes(mean_lon, mean_lat, color = year)) +
  geom_path(linewidth = 0.7, alpha = 0.75) +
  geom_point(aes(size = n), alpha = 0.9) +
  scale_color_viridis_c(option = "C") +
  scale_size_continuous(range = c(1.5, 5), name = "n records") +
  facet_wrap(~short_name, scales = "free") +
  coord_fixed(1.2) +
  labs(
    title = paste0("Same-season geographic shift — ", focus_season),
    subtitle = "Yearly mean lon/lat using only records from this season (comparable temperature regime)",
    x = "Mean longitude (°E)", y = "Mean latitude (°N)", color = "Year"
  ) +
  theme_dark_map()
ggsave("outputs_med_3d_pa/shifts_01_lonlat_tracks_SAME_SEASON.png",
       p_track, width = 11, height = 8, dpi = 150, bg = "#0d1117")

# ---- Plot 2: lon & lat time series, focus season ----------------------------
ts_long <- focus %>%
  pivot_longer(c(mean_lon, mean_lat), names_to = "coord", values_to = "value") %>%
  mutate(coord = recode(coord,
                        mean_lon = "Mean longitude (°E)",
                        mean_lat = "Mean latitude (°N)"))

p_ts_xy <- ggplot(ts_long, aes(year, value)) +
  geom_line(color = "#58a6ff", linewidth = 0.75) +
  geom_point(aes(size = n), color = "#79c0ff", alpha = 0.9) +
  facet_grid(coord ~ short_name, scales = "free_y") +
  scale_size_continuous(range = c(1, 3.5), guide = "none") +
  labs(
    title = paste0("Yearly mean lon/lat — ", focus_season, " only"),
    subtitle = "Dark mode · same months each year (see season label)",
    x = "Year", y = NULL
  ) +
  theme_dark_map(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "#c9d1d9"))
ggsave("outputs_med_3d_pa/shifts_02_lonlat_timeseries_SAME_SEASON.png",
       p_ts_xy, width = 12, height = 5.5, dpi = 150, bg = "#0d1117")

# ---- Plot 3: depth within focus season --------------------------------------
focus_depth <- focus %>% filter(!is.na(mean_depth), n_depth >= 1)
if (nrow(focus_depth) > 0) {
  p_depth <- ggplot(focus_depth, aes(year, mean_depth)) +
    geom_line(color = "#3fb950", linewidth = 0.8) +
    geom_point(aes(size = n_depth), color = "#56d364", alpha = 0.9) +
    scale_y_reverse() +
    scale_size_continuous(range = c(1.5, 4), name = "n with depth") +
    facet_wrap(~short_name, scales = "free_y") +
    labs(
      title = paste0("Yearly mean depth — ", focus_season, " only"),
      subtitle = "Y-axis reversed (deeper downward). Same season → more comparable temperature.",
      x = "Year", y = "Mean depth (m)"
    ) +
    theme_dark_map()
  ggsave("outputs_med_3d_pa/shifts_03_depth_timeseries_SAME_SEASON.png",
         p_depth, width = 11, height = 7, dpi = 150, bg = "#0d1117")
} else {
  message("No depth values in focus season for depth plot.")
}

# ---- Plot 4: combo for focus season -----------------------------------------
combo <- focus %>%
  pivot_longer(
    c(mean_lon, mean_lat, mean_depth),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("mean_lon", "mean_lat", "mean_depth"),
      labels = c("Mean lon (°E)", "Mean lat (°N)", "Mean depth (m)")
    )
  )

p_combo <- ggplot(combo, aes(year, value)) +
  geom_line(color = "#58a6ff", linewidth = 0.65, alpha = 0.9) +
  geom_point(aes(size = n), color = "#79c0ff", alpha = 0.85) +
  facet_grid(metric ~ short_name, scales = "free_y") +
  scale_size_continuous(range = c(0.8, 2.8), guide = "none") +
  labs(
    title = paste0("Same-season shifts — ", focus_season),
    subtitle = "Average position & depth by year within one season (dark mode)",
    x = "Year", y = NULL
  ) +
  theme_dark_map(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "#c9d1d9"))
ggsave("outputs_med_3d_pa/shifts_04_combo_SAME_SEASON.png",
       p_combo, width = 13, height = 7, dpi = 150, bg = "#0d1117")

# ---- Plot 5: all seasons comparison (lon over years) ------------------------
p_all_seasons <- ggplot(
  yearly_season,
  aes(year, mean_lat, color = season, group = season)
) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_point(aes(size = n), alpha = 0.85) +
  scale_color_manual(values = c(
    "Winter (DJF)" = "#79c0ff",
    "Spring (MAM)" = "#3fb950",
    "Summer (JJA)" = "#f0883e",
    "Autumn (SON)" = "#d2a8ff"
  )) +
  scale_size_continuous(range = c(1, 3), guide = "none") +
  facet_wrap(~short_name, scales = "free_y") +
  labs(
    title = "Mean latitude over years — separated by season",
    subtitle = "Do not pool seasons: each colour is its own comparable temperature window",
    x = "Year", y = "Mean latitude (°N)", color = "Season"
  ) +
  theme_dark_map()
ggsave("outputs_med_3d_pa/shifts_05_lat_all_seasons.png",
       p_all_seasons, width = 11, height = 8, dpi = 150, bg = "#0d1117")

# Trends within focus season only
trends <- focus %>%
  group_by(short_name) %>%
  group_modify(function(d, ...) {
    out <- list(season = focus_season)
    if (sum(!is.na(d$mean_lon)) >= 3) {
      out$lon_slope_per_year <- coef(lm(mean_lon ~ year, data = d))[["year"]]
    } else out$lon_slope_per_year <- NA_real_
    if (sum(!is.na(d$mean_lat)) >= 3) {
      out$lat_slope_per_year <- coef(lm(mean_lat ~ year, data = d))[["year"]]
    } else out$lat_slope_per_year <- NA_real_
    if (sum(!is.na(d$mean_depth)) >= 3) {
      out$depth_slope_per_year <- coef(lm(mean_depth ~ year, data = d))[["year"]]
    } else out$depth_slope_per_year <- NA_real_
    as.data.frame(out, stringsAsFactors = FALSE)
  }) %>%
  ungroup()
print(trends)
write_csv(trends, "outputs_med_3d_pa/shifts_linear_trends_SAME_SEASON_top6.csv")

message("\nDone (same season + dark mode).")
message("Focus season used for main plots: ", focus_season)
message("Open outputs_med_3d_pa/shifts_*SAME_SEASON*.png")
