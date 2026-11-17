# June–August haul counts by depth band (0–50, 50–100, 100–150 m) and year
# Output: species_selection/jun_aug_observations_by_depth_year.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

OUT <- "species_selection/jun_aug_observations_by_depth_year.csv"
SUMMARY <- "species_selection/jun_aug_observations_summary.csv"

SUMMER_MONTHS <- c(6L, 7L, 8L)
MONTH_LAB <- c(`6` = "June", `7` = "July", `8` = "August")

BBOX_ISRAEL <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 34.5)
BBOX_MEDITS <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 31.0, lat_max = 35.5)

in_bbox <- function(df, box) {
  df %>% filter(
    decimalLongitude >= box$lon_min, decimalLongitude <= box$lon_max,
    decimalLatitude >= box$lat_min, decimalLatitude <= box$lat_max
  )
}

medits_dd <- function(x) {
  x <- as.numeric(x)
  deg <- floor(x / 100)
  deg + (x - deg * 100) / 60
}

medits_haul_key <- function(country, area, vessel, year, haul_number) {
  paste(country, area, vessel, year, haul_number, sep = "|")
}

depth_layer3 <- function(d) {
  d <- as.numeric(d)
  case_when(
    !is.finite(d) ~ "depth_unknown",
    d < 0 ~ "<0 m",
    d < 50 ~ "0–50 m",
    d < 100 ~ "50–100 m",
    d < 150 ~ "100–150 m",
    TRUE ~ "150+ m (outside table)"
  )
}

load_all <- function() {
  parts <- list()
  if (file.exists("Combined_Menachem_Itai_corrected_May_2026_DBG.csv")) {
    parts[[length(parts) + 1]] <- read_csv(
      "Combined_Menachem_Itai_corrected_May_2026_DBG.csv", show_col_types = FALSE
    ) %>%
      transmute(
        dataset = "Israel_trawl_MIT",
        obs_key = as.character(SAMPLE_DESC),
        year = as.integer(YEAR),
        month = as.integer(MONTH),
        depth = as.numeric(DEPTH),
        decimalLatitude = as.numeric(LATITUDE),
        decimalLongitude = as.numeric(LONGITUDE)
      ) %>%
      in_bbox(BBOX_ISRAEL)
  }
  if (file.exists("data/raw/edelist_occurrence.txt")) {
    ed <- read_tsv("data/raw/edelist_occurrence.txt", show_col_types = FALSE, guess_max = 50000)
    id_col <- if ("id" %in% names(ed)) "id" else "occurrenceID"
    parts[[length(parts) + 1]] <- ed %>%
      filter(class == "Actinopterygii" | is.na(class)) %>%
      mutate(meta = str_match(
        .data[[id_col]],
        "^([0-9.]+)_([0-9.]+)_([0-9.]+)_([0-9]+)_([0-9]+)_([0-9]{4})"
      )) %>%
      transmute(
        dataset = "Edelist_OBIS",
        obs_key = str_replace(.data[[id_col]], "_Fish_.*$", ""),
        year = as.integer(meta[, 6]),
        month = as.integer(meta[, 5]),
        depth = as.numeric(meta[, 4]),
        decimalLatitude = as.numeric(meta[, 2]),
        decimalLongitude = as.numeric(meta[, 3])
      ) %>%
      in_bbox(BBOX_ISRAEL)
  }
  if (file.exists("MEDITS/medits_raw_data.csv") && file.exists("MEDITS/hauls_data.csv")) {
    hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE) %>%
      filter(validity == "V") %>%
      mutate(
        obs_key = medits_haul_key(country, area, vessel, year, haul_number),
        haul_month = as.integer(month),
        decimalLatitude = medits_dd(shooting_latitude),
        decimalLongitude = medits_dd(shooting_longitude),
        depth = (as.numeric(shooting_depth) + as.numeric(hauling_depth)) / 2
      )
    parts[[length(parts) + 1]] <- hauls %>%
      transmute(
        dataset = "MEDITS_raw",
        obs_key, year = as.integer(year), month = haul_month,
        depth, decimalLatitude, decimalLongitude, area = as.integer(area)
      ) %>%
      filter(area %in% c(22L, 25L)) %>%
      in_bbox(BBOX_MEDITS) %>%
      select(-area)
  }
  bind_rows(parts)
}

dat <- load_all()
hauls <- dat %>%
  filter(month %in% SUMMER_MONTHS, !is.na(obs_key), obs_key != "") %>%
  group_by(dataset, obs_key, month) %>%
  summarise(
    year = first(year[!is.na(year)]),
    haul_depth = median(depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    depth_band = depth_layer3(haul_depth),
    month_name = MONTH_LAB[as.character(month)]
  )

focus <- hauls %>%
  filter(depth_band %in% c("0–50 m", "50–100 m", "100–150 m"))

by_detail <- focus %>%
  count(year, month, month_name, depth_band, dataset, name = "observations") %>%
  arrange(depth_band, year, month, dataset)

write_csv(by_detail, OUT)

years_per_band <- focus %>%
  group_by(depth_band) %>%
  summarise(
    observations_jun_aug = n(),
    years = paste(sort(unique(year[!is.na(year)])), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(match(depth_band, c("0–50 m", "50–100 m", "100–150 m")))

write_csv(
  focus %>% count(depth_band, month_name, name = "observations"),
  "species_selection/jun_aug_observations_by_month.csv"
)
write_csv(
  focus %>% count(year, depth_band, name = "observations"),
  "species_selection/jun_aug_observations_by_year_depth.csv"
)
write_csv(years_per_band, SUMMARY)

message("June–August observations in 0–50 / 50–100 / 100–150 m:")
print(years_per_band)
message("Wrote ", OUT)
message("Wrote ", SUMMARY)
