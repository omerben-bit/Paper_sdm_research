# =============================================================================
# Summer (JJA) presence rate by year — one bar per summer year
# Y: % of eligible hauls (depth-matched, other fish present) where species caught
# X: year (discrete bars; each bar = June–August pooled)
# Output: primary_sdm/change_in_time/bar_summer_year_<species>.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

OUTPUT_ROOT <- "primary_sdm"
OUT_DIR <- file.path(OUTPUT_ROOT, "change_in_time")
OUT_DEPTH_DIR <- file.path(OUT_DIR, "by_depth")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DEPTH_DIR, showWarnings = FALSE, recursive = TRUE)

BBOX_ISRAEL <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 34.5)
SUMMER_MONTHS <- c(6L, 7L, 8L)

norm_species <- function(x) {
  x %>% as.character() %>% str_trim() %>%
    str_remove("\\s*\\([^)]*\\)$") %>% str_squish()
}

slug_species <- function(x) {
  x %>% str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>% str_sub(1, 60)
}

in_bbox <- function(df, box = BBOX_ISRAEL) {
  df %>%
    filter(
      decimalLongitude >= box$lon_min, decimalLongitude <= box$lon_max,
      decimalLatitude >= box$lat_min, decimalLatitude <= box$lat_max
    )
}

medits_dd <- function(x) {
  x <- as.numeric(x)
  deg <- floor(x / 100)
  deg + (x - deg * 100) / 60
}

theme_dark_bars <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.background   = element_rect(fill = "#0d1117", color = NA),
      panel.background  = element_rect(fill = "#161b22", color = NA),
      panel.border      = element_rect(fill = NA, color = "#30363d"),
      panel.grid.major  = element_line(color = "#21262d", linewidth = 0.3),
      axis.text         = element_text(color = "#c9d1d9"),
      axis.title        = element_text(color = "#e6edf3"),
      plot.title        = element_text(color = "#f0f6fc", face = "bold"),
      plot.subtitle     = element_text(color = "#8b949e"),
      legend.background = element_rect(fill = "#0d1117", color = NA),
      legend.key        = element_rect(fill = "#161b22", color = NA),
      legend.text       = element_text(color = "#c9d1d9"),
      legend.title      = element_text(color = "#e6edf3")
    )
}

load_survey_data <- function() {
  parts <- list()

  mit_path <- "Combined_Menachem_Itai_corrected_May_2026_DBG.csv"
  if (file.exists(mit_path)) {
    parts[["Israel_trawl_MIT"]] <- read_csv(mit_path, show_col_types = FALSE) %>%
      transmute(
        scientificName = as.character(ID_SPECIES),
        decimalLatitude = as.numeric(LATITUDE),
        decimalLongitude = as.numeric(LONGITUDE),
        dataset = "Israel_trawl_MIT",
        month = as.integer(MONTH),
        year = as.integer(YEAR),
        obs_key = as.character(SAMPLE_DESC),
        depth = as.numeric(DEPTH)
      )
  }

  if (file.exists("data/raw/edelist_occurrence.txt")) {
    ed <- read_tsv("data/raw/edelist_occurrence.txt", show_col_types = FALSE, guess_max = 50000)
    id_col <- if ("id" %in% names(ed)) "id" else "occurrenceID"
    parts[["Edelist_OBIS"]] <- ed %>%
      filter(class == "Actinopterygii" | is.na(class)) %>%
      mutate(meta = str_match(.data[[id_col]], "^([0-9.]+)_([0-9.]+)_([0-9.]+)_([0-9]+)_([0-9]+)_([0-9]{4})")) %>%
      transmute(
        scientificName = as.character(scientificName),
        decimalLatitude = as.numeric(meta[, 2]),
        decimalLongitude = as.numeric(meta[, 3]),
        dataset = "Edelist_OBIS",
        month = as.integer(meta[, 5]),
        year = as.integer(meta[, 6]),
        obs_key = str_replace(.data[[id_col]], "_Fish_.*$", ""),
        depth = as.numeric(meta[, 3])
      )
  }

  if (file.exists("MEDITS/medits_raw_data.csv") && file.exists("MEDITS/hauls_data.csv")) {
    catch <- read_csv("MEDITS/medits_raw_data.csv", show_col_types = FALSE,
                      locale = locale(encoding = "Latin1")) %>%
      mutate(
        country = as.character(country), area = as.integer(area),
        vessel = as.character(vessel), year = as.integer(year),
        haul_number = as.numeric(haul_number)
      ) %>%
      filter(nbtot > 0)
    hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE) %>%
      mutate(
        country = as.character(country), area = as.integer(area),
        vessel = as.character(vessel), year = as.integer(year),
        haul_number = as.numeric(haul_number),
        month = as.integer(month)
      ) %>%
      filter(validity == "V") %>%
      transmute(
        country, area, vessel, year, haul_number,
        haul_month = as.integer(month),
        decimalLatitude = medits_dd(shooting_latitude),
        decimalLongitude = medits_dd(shooting_longitude),
        depth = (as.numeric(shooting_depth) + as.numeric(hauling_depth)) / 2,
        obs_key = paste(country, area, vessel, year, haul_number, sep = "|")
      )
    taxa <- NULL
    if (file.exists("MEDITS/taxa_medits.csv")) {
      taxa <- read_csv("MEDITS/taxa_medits.csv", show_col_types = FALSE,
                       locale = locale(encoding = "Latin1"))
      code_col <- names(taxa)[str_detect(names(taxa), "^MEDITS")]
      name_col <- names(taxa)[str_detect(names(taxa), "Scientific Name")]
      taxa <- taxa %>%
        mutate(
          medits_code = .data[[code_col[1]]],
          scientificName = .data[[name_col[1]]]
        ) %>%
        select(medits_code, scientificName)
      catch <- catch %>% mutate(medits_code = paste0(genus, species)) %>%
        left_join(taxa, by = "medits_code")
    } else {
      catch <- catch %>% mutate(scientificName = paste0(genus, species))
    }
    parts[["MEDITS_raw"]] <- catch %>%
      inner_join(hauls, by = c("country", "area", "vessel", "year", "haul_number")) %>%
      transmute(
        scientificName = coalesce(scientificName, medits_code),
        decimalLatitude, decimalLongitude,
        dataset = "MEDITS_raw",
        month = haul_month,
        year = as.integer(year),
        obs_key,
        depth
      )
  }

  bind_rows(parts) %>%
    filter(
      !is.na(scientificName), scientificName != "",
      !is.na(year), is.finite(year),
      is.finite(decimalLatitude), is.finite(decimalLongitude),
      is.finite(depth)
    ) %>%
    in_bbox(BBOX_ISRAEL) %>%
    mutate(sp_key = norm_species(scientificName))
}

species_depth_band <- function(dat, focus_key, q = c(0.05, 0.95)) {
  band <- dat %>%
    filter(sp_key == focus_key, is.finite(depth)) %>%
    summarise(
      depth_lo = unname(quantile(depth, q[1], na.rm = TRUE)),
      depth_hi = unname(quantile(depth, q[2], na.rm = TRUE)),
      .groups = "drop"
    )
  if (nrow(band) == 0 || !is.finite(band$depth_lo)) {
    return(NULL)
  }
  band
}

depth_tertiles <- function(depth_lo, depth_hi) {
  if (!is.finite(depth_lo) || !is.finite(depth_hi)) {
    return(NULL)
  }
  if (depth_hi <= depth_lo) {
    depth_hi <- depth_lo + 1
  }
  breaks <- unique(seq(depth_lo, depth_hi, length.out = 4))
  if (length(breaks) < 4) {
    breaks <- c(depth_lo, depth_lo, depth_hi, depth_hi)
    breaks <- unique(breaks)
    while (length(breaks) < 4) {
      breaks <- c(breaks, depth_hi + 1)
    }
    breaks <- breaks[1:4]
  }
  labels <- sprintf("%.0f–%.0f m", breaks[1:3], breaks[2:4])
  list(breaks = breaks, labels = labels)
}

eligible_summer_hauls <- function(dat, focus_key) {
  band <- species_depth_band(dat, focus_key)
  if (is.null(band)) {
    return(list(hauls = tibble(), band = NULL))
  }

  summer <- dat %>%
    filter(month %in% SUMMER_MONTHS, !is.na(obs_key), obs_key != "", is.finite(depth))

  hauls <- summer %>%
    group_by(dataset, year, obs_key) %>%
    summarise(
      haul_depth = median(depth, na.rm = TRUE),
      has_focus = any(sp_key == focus_key),
      focus_sightings = sum(sp_key == focus_key),
      has_other = any(sp_key != focus_key),
      .groups = "drop"
    ) %>%
    filter(
      haul_depth >= band$depth_lo, haul_depth <= band$depth_hi,
      has_other
    )

  list(hauls = hauls, band = band)
}

# JJA hauls per year at species depth, other fish present
summer_year_pa <- function(dat, focus_key) {
  el <- eligible_summer_hauls(dat, focus_key)
  hauls <- el$hauls
  band <- el$band
  if (is.null(band) || nrow(hauls) == 0) {
    return(list(metrics = tibble(), band = band))
  }

  metrics <- hauls %>%
    group_by(year, dataset) %>%
    summarise(
      presence_hauls = sum(has_focus),
      absence_hauls = sum(!has_focus),
      total_hauls = n(),
      sightings = sum(focus_sightings),
      presence_pct = round(100 * sum(has_focus) / n(), 1),
      .groups = "drop"
    ) %>%
    arrange(dataset, year) %>%
    mutate(
      year_lab = factor(year, levels = sort(unique(year))),
      scientificName = NA_character_
    )

  list(metrics = metrics, band = band)
}

# Sightings on presence hauls, split into 3 depth bins within species band
summer_year_sightings_by_depth <- function(dat, focus_key) {
  el <- eligible_summer_hauls(dat, focus_key)
  hauls <- el$hauls
  band <- el$band
  if (is.null(band) || nrow(hauls) == 0) {
    return(list(metrics = tibble(), band = band, tertiles = NULL))
  }

  tri <- depth_tertiles(band$depth_lo, band$depth_hi)
  if (is.null(tri)) {
    return(list(metrics = tibble(), band = band, tertiles = NULL))
  }

  pres <- hauls %>%
    filter(has_focus, focus_sightings > 0) %>%
    mutate(
      depth_bin = cut(
        haul_depth,
        breaks = tri$breaks,
        include.lowest = TRUE,
        right = TRUE,
        labels = tri$labels
      )
    )

  if (nrow(pres) == 0) {
    return(list(metrics = tibble(), band = band, tertiles = tri))
  }

  metrics <- pres %>%
    group_by(year, dataset, depth_bin) %>%
    summarise(sightings = sum(focus_sightings), presence_hauls = n(), .groups = "drop") %>%
    arrange(dataset, year, depth_bin) %>%
    mutate(
      year_lab = factor(year, levels = sort(unique(pres$year))),
      depth_bin = factor(depth_bin, levels = tri$labels)
    )

  list(metrics = metrics, band = band, tertiles = tri)
}

plot_summer_year_depth_bars <- function(depth_pa, species, outfile) {
  metrics <- depth_pa$metrics
  band <- depth_pa$band
  if (nrow(metrics) == 0) {
    message("Skip (no depth-binned sightings): ", species)
    return(invisible(NULL))
  }

  depth_txt <- sprintf(
    "Species depth band %.0f–%.0f m · 3 equal depth bins · JJA · sightings on presence hauls",
    band$depth_lo, band$depth_hi
  )

  p <- ggplot(metrics, aes(year_lab, sightings, fill = depth_bin)) +
    geom_col(
      position = position_dodge(width = 0.82),
      width = 0.24,
      alpha = 0.92
    ) +
    facet_wrap(~dataset, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = c("#58a6ff", "#3fb950", "#d29922"),
      name = "Haul depth"
    ) +
    labs(
      title = paste(species, "— summer sightings by depth (per year)"),
      subtitle = depth_txt,
      x = "Summer year",
      y = "Sightings (records)"
    ) +
    theme_dark_bars() +
    theme(legend.position = "top")

  n_yr <- n_distinct(metrics$year)
  ggsave(
    outfile, p,
    width = max(7, 0.55 * n_yr + 4),
    height = 3.2 + 2.8 * n_distinct(metrics$dataset),
    dpi = 150, bg = "#0d1117"
  )
  message("Wrote ", outfile)

  metrics %>%
    mutate(scientificName = species) %>%
    select(scientificName, year, dataset, depth_bin, sightings, presence_hauls) %>%
    write_csv(sub("\\.png$", "_metrics.csv", outfile))
}

plot_summer_year_bars <- function(pa, species, outfile) {
  metrics <- pa$metrics
  band <- pa$band
  if (nrow(metrics) == 0) {
    message("Skip (no depth-matched summer-by-year data): ", species)
    return(invisible(NULL))
  }

  depth_txt <- sprintf(
    "Depth %.0f–%.0f m (p5–p95 of species records) · JJA pooled per year · ",
    band$depth_lo, band$depth_hi
  )

  p <- ggplot(metrics, aes(year_lab, presence_pct, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72, alpha = 0.92) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_brewer(palette = "Set2", name = "Survey") +
    labs(
      title = paste(species, "— summer presence rate over time"),
      subtitle = paste0(
        depth_txt,
        "share of eligible hauls (other fish caught) with species present"
      ),
      x = "Summer year",
      y = "Presence (% of hauls)"
    ) +
    theme_dark_bars() +
    theme(legend.position = "top")

  n_yr <- n_distinct(metrics$year)
  ggsave(
    outfile, p,
    width = max(7, 0.45 * n_yr + 4),
    height = 5.5,
    dpi = 150, bg = "#0d1117"
  )
  message("Wrote ", outfile)

  metrics %>%
    mutate(scientificName = species) %>%
    select(
      scientificName, year, dataset,
      presence_hauls, absence_hauls, total_hauls,
      presence_pct, sightings
    ) %>%
    write_csv(sub("\\.png$", "_metrics.csv", outfile))
}

# Remove prior bar outputs in this folder (png + csv from older runs)
old_png <- list.files(OUT_DIR, pattern = "\\.png$", full.names = TRUE)
old_csv <- list.files(OUT_DIR, pattern = "\\.csv$", full.names = TRUE)
for (f in c(old_png, old_csv)) {
  file.remove(f)
}

message("Loading survey records...")
dat <- load_survey_data()
if (nrow(dat) == 0) stop("No survey records — check data paths.")

dat_summer <- dat %>% filter(month %in% SUMMER_MONTHS)

top10 <- dat_summer %>%
  count(sp_key, name = "n_summer") %>%
  arrange(desc(n_summer)) %>%
  slice_head(n = 10)

name_pick <- dat %>%
  count(sp_key, scientificName) %>%
  group_by(sp_key) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sp_key, scientificName)

species_top10 <- top10 %>%
  left_join(name_pick, by = "sp_key")

message("Top 10 species (summer record count):")
print(species_top10 %>% select(scientificName, n_summer))

write_csv(species_top10, file.path(OUT_DIR, "top10_summer_species.csv"))

all_metrics <- list()

for (i in seq_len(nrow(species_top10))) {
  sp <- species_top10$scientificName[i]
  key <- species_top10$sp_key[i]
  pa <- summer_year_pa(dat, key)
  slug <- slug_species(sp)
  plot_summer_year_bars(
    pa, sp,
    file.path(OUT_DIR, paste0("bar_summer_year_", slug, ".png"))
  )
  if (nrow(pa$metrics) > 0) {
    all_metrics[[slug]] <- pa$metrics %>% mutate(scientificName = sp, .before = 1)
  }
}

if (length(all_metrics) > 0) {
  write_csv(
    bind_rows(all_metrics),
    file.path(OUT_DIR, "bar_summer_year_all_species_metrics.csv")
  )
}

old_depth_png <- list.files(OUT_DEPTH_DIR, pattern = "\\.png$", full.names = TRUE)
old_depth_csv <- list.files(OUT_DEPTH_DIR, pattern = "\\.csv$", full.names = TRUE)
for (f in c(old_depth_png, old_depth_csv)) {
  file.remove(f)
}

depth_metrics <- list()

for (i in seq_len(nrow(species_top10))) {
  sp <- species_top10$scientificName[i]
  key <- species_top10$sp_key[i]
  slug <- slug_species(sp)
  dep <- summer_year_sightings_by_depth(dat, key)
  plot_summer_year_depth_bars(
    dep, sp,
    file.path(OUT_DEPTH_DIR, paste0("bar_summer_year_depth_", slug, ".png"))
  )
  if (nrow(dep$metrics) > 0) {
    depth_metrics[[slug]] <- dep$metrics %>% mutate(scientificName = sp, .before = 1)
  }
}

if (length(depth_metrics) > 0) {
  write_csv(
    bind_rows(depth_metrics),
    file.path(OUT_DEPTH_DIR, "bar_summer_year_depth_all_species_metrics.csv")
  )
}

message("\nDone. Summer-by-year bars in ", OUT_DIR, " (bar_summer_year_*.png)")
message("Depth-binned sightings in ", OUT_DEPTH_DIR, " (bar_summer_year_depth_*.png)")
