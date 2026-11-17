# =============================================================================
# Dataset inventory: sightings, columns, depth distribution, spatial extent
# Output: dataset_summary/ (markdown + CSV tables)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("gbif_paper.R")

MONTH_LAB <- month.abb[1:12]

OUT_DIR <- "dataset_summary"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BBOX_ISRAEL <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 34.5)
BBOX_MEDITS <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 31.0, lat_max = 35.5)
DEPTH_BIN_M <- 50L

file_cols <- function(path, delim = NULL) {
  if (!file.exists(path)) {
    return(tibble(file = path, column_name = character(), column_index = integer()))
  }
  if (is.null(delim)) {
    delim <- if (str_detect(path, "\\.tsv$|occurrence\\.txt$")) "\t" else ","
  }
  hdr <- if (delim == "\t") {
    names(read_tsv(path, n_max = 0, show_col_types = FALSE))
  } else {
    names(read_csv(path, n_max = 0, show_col_types = FALSE))
  }
  tibble(
    file = path,
    column_index = seq_along(hdr),
    column_name = hdr
  )
}

in_bbox <- function(df, box) {
  df %>%
    filter(
      decimalLongitude >= box$lon_min, decimalLongitude <= box$lon_max,
      decimalLatitude >= box$lat_min, decimalLatitude <= box$lat_max
    )
}

depth_band <- function(d) {
  d <- as.numeric(d)
  out <- rep(NA_character_, length(d))
  ok <- is.finite(d)
  if (!any(ok)) {
    return(out)
  }
  dd <- d[ok]
  lo <- pmax(0, floor(dd / DEPTH_BIN_M) * DEPTH_BIN_M)
  hi <- lo + DEPTH_BIN_M
  out[ok] <- ifelse(
    dd >= 500, "500+ m",
    ifelse(dd < 0, "<0 m", sprintf("%d–%d m", lo, hi))
  )
  out[!ok] <- "depth_unknown"
  out
}

order_depth_band <- function(x) {
  ord <- c(
    "depth_unknown", "<0 m",
    sprintf("%d–%d m", seq(0, 450, DEPTH_BIN_M), seq(DEPTH_BIN_M, 500, DEPTH_BIN_M)),
    "500+ m"
  )
  factor(x, levels = unique(c(ord, setdiff(unique(x), ord))))
}

loc_summary <- function(df, dataset, scope) {
  if (nrow(df) == 0) {
    return(tibble(
      dataset = dataset, scope = scope,
      n_with_coordinates = 0L, lon_min = NA_real_, lon_max = NA_real_,
      lat_min = NA_real_, lat_max = NA_real_, centroid_lon = NA_real_,
      centroid_lat = NA_real_, n_grid_cells_0p5deg = 0L
    ))
  }
  g <- df %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
  tibble(
    dataset = dataset,
    scope = scope,
    n_with_coordinates = nrow(g),
    lon_min = min(g$decimalLongitude),
    lon_max = max(g$decimalLongitude),
    lat_min = min(g$decimalLatitude),
    lat_max = max(g$decimalLatitude),
    centroid_lon = mean(g$decimalLongitude),
    centroid_lat = mean(g$decimalLatitude),
    n_grid_cells_0p5deg = n_distinct(
      paste0(
        floor(g$decimalLongitude * 2) / 2, "_",
        floor(g$decimalLatitude * 2) / 2
      )
    )
  )
}

depth_tables <- function(df, dataset, scope, level = c("sighting", "observation")) {
  level <- match.arg(level)
  if (nrow(df) == 0) {
    return(tibble(dataset = dataset, scope = scope, level = level,
                  depth_band = character(), n = integer()))
  }
  if (level == "sighting") {
    out <- df %>%
      mutate(depth_band = depth_band(depth)) %>%
      count(depth_band, name = "n")
  } else {
    if (!"obs_key" %in% names(df) || all(is.na(df$obs_key))) {
      out <- df %>%
        mutate(depth_band = depth_band(depth)) %>%
        count(depth_band, name = "n")
    } else {
      haul_d <- df %>%
        filter(!is.na(obs_key), obs_key != "") %>%
        group_by(obs_key) %>%
        summarise(haul_depth = median(depth, na.rm = TRUE), .groups = "drop")
      out <- haul_d %>%
        mutate(depth_band = depth_band(haul_depth)) %>%
        count(depth_band, name = "n")
    }
  }
  out %>%
    mutate(
      dataset = dataset,
      scope = scope,
      level = level,
      depth_band = order_depth_band(depth_band)
    ) %>%
    arrange(depth_band) %>%
    select(dataset, scope, level, depth_band, n)
}

medits_dd <- function(x) {
  x <- as.numeric(x)
  deg <- floor(x / 100)
  deg + (x - deg * 100) / 60
}

norm_species <- function(x) {
  x %>% as.character() %>% str_trim() %>%
    str_remove("\\s*\\([^)]*\\)$") %>% str_squish()
}

medits_haul_key <- function(country, area, vessel, year, haul_number) {
  paste(country, area, vessel, year, haul_number, sep = "|")
}

load_mit <- function(scope = c("all_rows", "israel_box")) {
  scope <- match.arg(scope)
  path <- "Combined_Menachem_Itai_corrected_May_2026_DBG.csv"
  if (!file.exists(path)) {
    return(tibble())
  }
  mit <- read_csv(path, show_col_types = FALSE) %>%
    transmute(
      scientificName = as.character(ID_SPECIES),
      decimalLatitude = as.numeric(LATITUDE),
      decimalLongitude = as.numeric(LONGITUDE),
      depth = as.numeric(DEPTH),
      dataset = "Israel_trawl_MIT",
      year = as.integer(YEAR),
      month = as.integer(MONTH),
      obs_key = as.character(SAMPLE_DESC)
    ) %>%
    filter(!is.na(scientificName))
  if (scope == "israel_box") {
    mit <- in_bbox(mit, BBOX_ISRAEL)
  }
  mit
}

load_edelist <- function(scope = c("all_rows", "israel_box", "actinopterygii_israel")) {
  scope <- match.arg(scope)
  path <- "data/raw/edelist_occurrence.txt"
  if (!file.exists(path)) {
    return(tibble())
  }
  ed <- read_tsv(path, show_col_types = FALSE, guess_max = 50000)
  id_col <- if ("id" %in% names(ed)) "id" else "occurrenceID"
  meta <- str_match(
    ed[[id_col]],
    "^([0-9.]+)_([0-9.]+)_([0-9.]+)_([0-9]+)_([0-9]+)_([0-9]{4})"
  )
  out <- ed %>%
    transmute(
      scientificName = as.character(scientificName),
      decimalLatitude = as.numeric(meta[, 2]),
      decimalLongitude = as.numeric(meta[, 3]),
      depth = as.numeric(meta[, 4]),
      dataset = "Edelist_OBIS",
      month = as.integer(meta[, 5]),
      year = as.integer(meta[, 6]),
      obs_key = str_replace(ed[[id_col]], "_Fish_.*$", ""),
      class = if ("class" %in% names(ed)) ed$class else NA_character_
    )
  if (scope == "actinopterygii_israel") {
    out <- out %>%
      filter(class == "Actinopterygii" | is.na(class)) %>%
      in_bbox(BBOX_ISRAEL) %>%
      filter(!is.na(scientificName), scientificName != "")
  } else if (scope == "israel_box") {
    out <- in_bbox(out, BBOX_ISRAEL)
  }
  out %>% select(-class)
}

load_medits_catch <- function(scope = c("raw_rows", "levant_gsa22_25")) {
  scope <- match.arg(scope)
  if (!file.exists("MEDITS/medits_raw_data.csv") || !file.exists("MEDITS/hauls_data.csv")) {
    return(tibble())
  }
  catch <- read_csv("MEDITS/medits_raw_data.csv", show_col_types = FALSE,
                    locale = locale(encoding = "Latin1"))
  hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE) %>%
    filter(validity == "V") %>%
    mutate(
      haul_month = as.integer(month),
      decimalLatitude = medits_dd(shooting_latitude),
      decimalLongitude = medits_dd(shooting_longitude),
      depth = (as.numeric(shooting_depth) + as.numeric(hauling_depth)) / 2,
      obs_key = medits_haul_key(country, area, vessel, year, haul_number)
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
      select(medits_code, scientificName, CATFAU)
  }
  catch <- catch %>% filter(nbtot > 0)
  catch <- catch %>% mutate(medits_code = paste0(genus, species))
  if (!is.null(taxa)) {
    catch <- catch %>% left_join(taxa, by = "medits_code")
    catch <- catch %>%
      filter(is.na(CATFAU) | CATFAU %in% c("Ao", "Ae", "C"))
  }
  out <- catch %>%
    inner_join(
      hauls %>% select(country, area, vessel, year, haul_number, haul_month,
                       decimalLatitude, decimalLongitude, depth, obs_key),
      by = c("country", "area", "vessel", "year", "haul_number")
    ) %>%
    transmute(
      scientificName = coalesce(scientificName, medits_code),
      decimalLatitude, decimalLongitude, depth,
      dataset = "MEDITS_raw",
      year = as.integer(year),
      month = as.integer(haul_month),
      gsa = area,
      obs_key
    )
  if (scope == "levant_gsa22_25") {
    out <- out %>% filter(gsa %in% c(22L, 25L)) %>% in_bbox(BBOX_MEDITS)
  }
  out %>% select(-gsa)
}

load_singer <- function(scope = c("global_depth_ok", "israel_box")) {
  scope <- match.arg(scope)
  path <- "data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv"
  if (!file.exists(path)) {
    return(tibble())
  }
  hdr <- names(read_csv(path, n_max = 0, show_col_types = FALSE))
  pick <- function(cands) cands[cands %in% hdr][1]
  sp_c <- pick(c("WoRMS_scientificName", "SOcurated_scientificName", "dwc_scientificName", "scientificName"))
  lat_c <- pick(c("iDigBio_dwc_decimalLatitude", "dwc_decimalLatitude", "decimalLatitude"))
  lon_c <- pick(c("iDigBio_dwc_decimalLongitude", "dwc_decimalLongitude", "decimalLongitude"))
  dep_c <- pick(c("SOcurated_depth", "depth"))
  yr_c <- pick(c("iDigBio_dwc_year", "dwc_year", "year"))
  mo_c <- pick(c("iDigBio_dwc_month", "dwc_month", "month"))
  if (any(is.na(c(sp_c, lat_c, lon_c, dep_c)))) {
    return(tibble())
  }
  use_cols <- unique(na.omit(c(sp_c, lat_c, lon_c, dep_c, yr_c, mo_c)))
  raw <- read_csv(path, show_col_types = FALSE, col_select = all_of(use_cols), progress = FALSE)
  out <- raw %>%
    transmute(
      scientificName = norm_species(as.character(.data[[sp_c]])),
      decimalLatitude = as.numeric(.data[[lat_c]]),
      decimalLongitude = as.numeric(.data[[lon_c]]),
      depth = as.numeric(.data[[dep_c]]),
      dataset = "Singer_Owens_iDigBio",
      year = if (!is.na(yr_c)) suppressWarnings(as.integer(.data[[yr_c]])) else NA_integer_,
      month = if (!is.na(mo_c)) suppressWarnings(as.integer(.data[[mo_c]])) else NA_integer_,
      obs_key = NA_character_
    ) %>%
    filter(
      !is.na(scientificName), scientificName != "",
      !scientificName %in% c("Unaccepted", "Equivocal"),
      is.finite(decimalLatitude), is.finite(decimalLongitude), is.finite(depth)
    )
  if (scope == "israel_box") {
    out <- in_bbox(out, BBOX_ISRAEL)
  }
  out
}

load_gbif_paper <- function() {
  gb <- load_gbif_paper_cleaned()
  if (nrow(gb) == 0) {
    return(tibble())
  }
  gb %>%
    transmute(
      scientificName = norm_species(scientificName),
      decimalLatitude = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      depth = as.numeric(depth),
      dataset = GBIF_PAPER_DATASET,
      year = suppressWarnings(as.integer(substr(eventDate, 1, 4))),
      month = suppressWarnings(as.integer(substr(eventDate, 6, 7))),
      obs_key = NA_character_
    )
}

# ---- Processed pipeline exports (data/processed/) ----------------------------

load_processed_occ_dwc <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  raw <- read_csv(path, show_col_types = FALSE)
  if (!all(c("scientificName", "decimalLatitude", "decimalLongitude") %in% names(raw))) {
    return(tibble())
  }
  raw %>%
    transmute(
      scientificName = norm_species(scientificName),
      decimalLatitude = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      depth = if ("depth" %in% names(raw)) as.numeric(depth) else NA_real_,
      year = if ("eventDate" %in% names(raw)) {
        suppressWarnings(as.integer(substr(eventDate, 1, 4)))
      } else if ("year" %in% names(raw)) {
        suppressWarnings(as.integer(year))
      } else {
        NA_integer_
      },
      month = if ("eventDate" %in% names(raw)) {
        suppressWarnings(as.integer(substr(eventDate, 6, 7)))
      } else if ("month" %in% names(raw)) {
        suppressWarnings(as.integer(month))
      } else {
        NA_integer_
      },
      obs_key = if ("station_id" %in% names(raw)) {
        as.character(station_id)
      } else if ("obs_key" %in% names(raw)) {
        as.character(obs_key)
      } else {
        NA_character_
      }
    )
}

load_processed_singer_subset <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      scientificName = norm_species(scientificName),
      decimalLatitude = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      depth = if ("depth" %in% names(.)) as.numeric(depth) else NA_real_,
      year = if ("year" %in% names(.)) suppressWarnings(as.integer(year)) else NA_integer_,
      month = NA_integer_,
      obs_key = NA_character_
    )
}

processed_layers_manifest <- function() {
  list(
    list(
      dataset = "proc_gbif_paper_cleaned",
      scope = "Paper GBIF only (`01` bbox + CoordinateCleaner)",
      file = GBIF_PAPER_CLEANED_CSV,
      notes = "Canonical GBIF for analysis · not `gbif_fish_levant.rds`",
      loader = function() load_processed_occ_dwc(GBIF_PAPER_CLEANED_CSV)
    ),
    list(
      dataset = "proc_med_fish_occurrences_combined",
      scope = "All rows in `01_med_fish_3d_pa_maps.R` combined export",
      file = "data/processed/med_fish_occurrences_combined.csv",
      notes = "GBIF + Edelist + optional Michmoret · columns include source, depth_bin",
      loader = function() load_processed_occ_dwc("data/processed/med_fish_occurrences_combined.csv")
    ),
    list(
      dataset = "proc_pa_focus_species",
      scope = "Presence / absence / pseudo-absence table for first 3D map species",
      file = "data/processed/pa_focus_species.csv",
      notes = "From `01_med_fish_3d_pa_maps.R` · pa, pa_type, station_id",
      loader = function() load_processed_occ_dwc("data/processed/pa_focus_species.csv")
    ),
    list(
      dataset = "proc_singer_owens_israel_box",
      scope = "Singer & Owens iDigBio — Israel coast box, depth required",
      file = "data/processed/singer_owens_israel_box.csv",
      notes = "Cached subset of Zenodo extract (script 03 / maps)",
      loader = function() load_processed_singer_subset("data/processed/singer_owens_israel_box.csv")
    ),
    list(
      dataset = "proc_singer_med_focus5",
      scope = "Singer & Owens — Mediterranean clip, FOCUS_FIVE species only",
      file = "data/processed/singer_med_focus5.csv",
      notes = "Cached for `04_med_2d_maps.R`",
      loader = function() load_processed_singer_subset("data/processed/singer_med_focus5.csv")
    )
  )
}

summarize_processed_by_source <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  raw <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"source" %in% names(raw)) {
    return(tibble())
  }
  raw %>%
    mutate(
      occurrenceStatus = if ("occurrenceStatus" %in% names(.)) {
        tolower(as.character(occurrenceStatus))
      } else {
        "present"
      }
    ) %>%
    count(source, occurrenceStatus, name = "rows", .drop = FALSE) %>%
    arrange(source, occurrenceStatus)
}

valid_year_month <- function(df) {
  df %>%
    filter(
      !is.na(year), year >= 1800L, year <= as.integer(format(Sys.Date(), "%Y")) + 1L,
      !is.na(month), month >= 1L, month <= 12L
    )
}

month_year_sightings <- function(df, dataset, scope) {
  if (nrow(df) == 0) {
    return(tibble(dataset = dataset, scope = scope, year = integer(),
                  month = integer(), sightings = integer()))
  }
  valid_year_month(df) %>%
    count(year, month, name = "sightings") %>%
    mutate(dataset = dataset, scope = scope, .before = 1)
}

month_year_observations <- function(df, dataset, scope) {
  if (nrow(df) == 0) {
    return(tibble(dataset = dataset, scope = scope, year = integer(),
                  month = integer(), observations = integer()))
  }
  base <- valid_year_month(df)
  if ("obs_key" %in% names(base) && any(!is.na(base$obs_key) & base$obs_key != "")) {
    base <- base %>%
      filter(!is.na(obs_key), obs_key != "") %>%
      group_by(obs_key) %>%
      summarise(year = first(year), month = first(month), .groups = "drop")
  }
  base %>%
    count(year, month, name = "observations") %>%
    mutate(dataset = dataset, scope = scope, .before = 1)
}

prep_joint_obs_id <- function(df, dataset) {
  if (nrow(df) == 0) {
    return(tibble())
  }
  df %>%
    mutate(
      dataset = dataset,
      obs_id = if_else(
        !is.na(.data$obs_key) & .data$obs_key != "",
        paste0(dataset, "::", .data$obs_key),
        NA_character_
      )
    ) %>%
    group_by(dataset) %>%
    mutate(
      obs_id = if_else(
        is.na(obs_id),
        paste0(dataset, "::occ::", row_number()),
        obs_id
      )
    ) %>%
    ungroup()
}

build_joint_observations <- function(layers) {
  all_rec <- bind_rows(lapply(layers, function(L) {
    prep_joint_obs_id(L$data, L$dataset)
  }))
  if (nrow(all_rec) == 0) {
    empty <- tibble()
    return(list(
      units = empty,
      by_year_month = empty,
      by_month = empty,
      by_year = empty,
      by_depth = empty,
      by_dataset_month = empty,
      by_dataset_year = empty,
      by_dataset_depth = empty
    ))
  }

  first_finite <- function(x) {
    y <- x[!is.na(x)]
    if (length(y)) y[[1]] else NA_integer_
  }

  units <- all_rec %>%
    group_by(obs_id, dataset) %>%
    summarise(
      year = first_finite(year),
      month = first_finite(month),
      haul_depth = median(depth, na.rm = TRUE),
      sightings = n(),
      .groups = "drop"
    )

  vy <- units %>%
    filter(
      !is.na(year), year >= 1800L, year <= as.integer(format(Sys.Date(), "%Y")) + 1L,
      !is.na(month), month >= 1L, month <= 12L
    )

  list(
    units = units,
    by_year_month = vy %>% count(year, month, name = "observations") %>% arrange(year, month),
    by_month = vy %>% count(month, name = "observations") %>%
      mutate(month_name = factor(MONTH_LAB[month], levels = MONTH_LAB)) %>%
      arrange(month),
    by_year = vy %>% count(year, name = "observations") %>% arrange(year),
    by_depth = units %>%
      mutate(depth_band = order_depth_band(depth_band(haul_depth))) %>%
      count(depth_band, name = "observations") %>%
      arrange(depth_band),
    by_dataset_year_month = vy %>%
      count(dataset, year, month, name = "observations") %>%
      arrange(dataset, year, month),
    by_dataset_month = vy %>%
      count(dataset, month, name = "observations") %>%
      mutate(month_name = factor(MONTH_LAB[month], levels = MONTH_LAB)),
    by_dataset_year = vy %>% count(dataset, year, name = "observations"),
    by_dataset_depth = units %>%
      mutate(depth_band = order_depth_band(depth_band(haul_depth))) %>%
      count(dataset, depth_band, name = "observations")
  )
}

plot_joint_stacked <- function(
    tab, x_var, title, subtitle, ylab = "Count", count_col = "observations"
) {
  if (nrow(tab) == 0 || !count_col %in% names(tab)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data", size = 4) +
        theme_void() +
        labs(title = title, subtitle = subtitle)
    )
  }
  ggplot(tab, aes(.data[[x_var]], .data[[count_col]], fill = .data$dataset)) +
    geom_col(position = "stack", width = 0.78) +
    scale_fill_brewer(palette = "Set2", name = "Database") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

build_joint_sightings <- function(layers) {
  all_rec <- bind_rows(lapply(layers, function(L) {
    L$data %>% mutate(dataset = L$dataset)
  }))
  if (nrow(all_rec) == 0) {
    empty <- tibble()
    return(list(
      by_year_month = empty,
      by_month = empty,
      by_year = empty,
      by_depth = empty,
      by_dataset_month = empty,
      by_dataset_year = empty,
      by_dataset_depth = empty,
      total = 0L
    ))
  }

  vy <- valid_year_month(all_rec)

  list(
    total = nrow(all_rec),
    by_year_month = vy %>% count(year, month, name = "sightings") %>% arrange(year, month),
    by_month = vy %>%
      count(month, name = "sightings") %>%
      mutate(month_name = factor(MONTH_LAB[month], levels = MONTH_LAB)) %>%
      arrange(month),
    by_year = vy %>% count(year, name = "sightings") %>% arrange(year),
    by_depth = all_rec %>%
      mutate(depth_band = order_depth_band(depth_band(depth))) %>%
      count(depth_band, name = "sightings") %>%
      arrange(depth_band),
    by_dataset_year_month = vy %>%
      count(dataset, year, month, name = "sightings") %>%
      arrange(dataset, year, month),
    by_dataset_month = vy %>%
      count(dataset, month, name = "sightings") %>%
      mutate(month_name = factor(MONTH_LAB[month], levels = MONTH_LAB)),
    by_dataset_year = vy %>% count(dataset, year, name = "sightings"),
    by_dataset_depth = all_rec %>%
      mutate(depth_band = order_depth_band(depth_band(depth))) %>%
      count(dataset, depth_band, name = "sightings")
  )
}

open_pdf_device <- function(out_path) {
  if (file.exists(out_path)) {
    removed <- tryCatch(
      { file.remove(out_path); TRUE },
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!removed) {
      out_path <- file.path(
        dirname(out_path),
        paste0(tools::file_path_sans_ext(basename(out_path)), "_latest.pdf")
      )
      message("Primary PDF locked (close in viewer); writing ", out_path)
    }
  }
  pdf(out_path, width = 11, height = 8.5, onefile = TRUE)
  out_path
}

print_combined_metric_block <- function(
    joint_part, value_col, label, subtitle_suffix
) {
  print(plot_month_year_heatmap(
    joint_part$by_year_month,
    value_col,
    paste("All databases —", label, "by year and month"),
    paste("Summed across surveys;", subtitle_suffix)
  ))
  cal <- joint_part$by_month %>%
    rename(n = .data[[value_col]]) %>%
    select(month, month_name, n)
  print(plot_month_calendar(
    cal,
    paste("All databases —", label, "by calendar month (all years)")
  ))
  print(plot_joint_stacked(
    joint_part$by_dataset_month,
    "month_name",
    paste(label, "by month — per database (stacked)"),
    "Heights sum to combined monthly total",
    label,
    count_col = value_col
  ))

  if (nrow(joint_part$by_year) > 0) {
    p_yr <- ggplot(joint_part$by_year, aes(.data$year, .data[[value_col]])) +
      geom_col(fill = "#d29922", width = 0.75) +
      scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
      labs(
        title = paste("All databases —", label, "by year"),
        x = "Year",
        y = label
      ) +
      theme_minimal(base_size = 10)
    print(p_yr)
    print(plot_joint_stacked(
      joint_part$by_dataset_year,
      "year",
      paste(label, "by year — per database (stacked)"),
      "Heights sum to combined yearly total",
      label,
      count_col = value_col
    ))
  }

  dep <- joint_part$by_depth %>% select(depth_band, n = all_of(value_col))
  print(plot_depth_distribution(
    dep,
    paste("All databases —", label, "by depth"),
    "50 m bins · record depth (sightings) or haul median (observations)",
    label,
    fill = if (value_col == "sightings") "#58a6ff" else "#3fb950"
  ))
  jd <- joint_part$by_dataset_depth %>%
    mutate(depth_band = order_depth_band(as.character(depth_band)))
  print(plot_joint_stacked(
    jd,
    "depth_band",
    paste(label, "by depth — per database (stacked)"),
    "Stacked contribution by depth band",
    label,
    count_col = value_col
  ))
}

write_combined_pdf <- function(sightings_summary, joint_obs, joint_sight, out_path) {
  out_path <- open_pdf_device(out_path)

  grid::grid.newpage()
  grid::grid.text(
    "All databases combined",
    y = 0.93,
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )
  grid::grid.text(
    "Master's project — pooled totals (MIT, Edelist, MEDITS, Singer, GBIF)",
    y = 0.88,
    gp = grid::gpar(fontsize = 11)
  )
  grid::grid.text(
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
    y = 0.84,
    gp = grid::gpar(fontsize = 10, col = "grey40")
  )

  y <- 0.76
  grid::grid.text("Per-database counts (analysis scope)", x = 0.05, y = y,
                  just = "left", gp = grid::gpar(fontsize = 11, fontface = "bold"))
  y <- y - 0.05
  for (i in seq_len(nrow(sightings_summary))) {
    r <- sightings_summary[i, ]
    grid::grid.text(
      sprintf(
        "%s — %s sightings, %s observations, %s species",
        r$dataset,
        format(r$sightings, big.mark = ","),
        format(r$observations, big.mark = ","),
        format(r$species, big.mark = ",")
      ),
      x = 0.05, y = y, just = "left", gp = grid::gpar(fontsize = 9.5)
    )
    y <- y - 0.035
  }

  y <- y - 0.04
  grid::grid.text(
    sprintf(
      "COMBINED: %s sightings · %s observation units",
      format(joint_sight$total, big.mark = ","),
      format(nrow(joint_obs$units), big.mark = ",")
    ),
    x = 0.05, y = y, just = "left",
    gp = grid::gpar(fontsize = 11, fontface = "bold", col = "#2f81f7")
  )

  grid::grid.newpage()
  grid::grid.text(
    "Part 1 — Sightings (fish records, all rows summed)",
    y = 0.95,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  print_combined_metric_block(
    joint_sight, "sightings", "Sightings",
    "each row is one record"
  )

  grid::grid.newpage()
  grid::grid.text(
    "Part 2 — Observations (hauls / stations, deduplicated per database)",
    y = 0.95,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  print_combined_metric_block(
    joint_obs, "observations", "Observations",
    "obs_id = database::haul or database::occ::n"
  )

  dev.off()
  message("Combined PDF saved: ", out_path)
  invisible(out_path)
}

month_calendar <- function(df, dataset, scope, level = c("sighting", "observation")) {
  level <- match.arg(level)
  if (level == "sighting") {
    tab <- month_year_sightings(df, dataset, scope)
    tab %>% rename(n = sightings)
  } else {
    tab <- month_year_observations(df, dataset, scope)
    tab %>% rename(n = observations)
  } %>%
    group_by(dataset, scope, month) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(month_name = factor(MONTH_LAB[month], levels = MONTH_LAB))
}

plot_month_year_heatmap <- function(tab, value_col, title, subtitle) {
  if (nrow(tab) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No year/month data", size = 5) +
        theme_void() +
        labs(title = title, subtitle = subtitle)
    )
  }
  tab <- tab %>%
    mutate(
      month_name = factor(MONTH_LAB[month], levels = MONTH_LAB),
      year_lab = factor(year)
    )
  show_lbl <- max(tab[[value_col]], na.rm = TRUE) <= 2500
  p <- ggplot(tab, aes(.data$month_name, .data$year_lab, fill = .data[[value_col]])) +
    geom_tile(color = "grey35", linewidth = 0.2)
  if (show_lbl) {
    p <- p + geom_text(aes(label = .data[[value_col]]), size = 2.8, color = "white")
  }
  p +
    scale_fill_viridis_c(option = "magma", name = NULL, trans = "sqrt") +
    labs(title = title, subtitle = subtitle, x = "Month", y = "Year") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_month_calendar <- function(cal, title) {
  if (nrow(cal) == 0 || !"month" %in% names(cal)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No month metadata", size = 4) +
        theme_void() +
        labs(title = title)
    )
  }
  cal <- cal %>%
    mutate(month_name = factor(MONTH_LAB[month], levels = MONTH_LAB))
  ggplot(cal, aes(.data$month_name, .data$n)) +
    geom_col(fill = "steelblue", width = 0.75) +
    labs(title = title, x = "Month (all years pooled)", y = "Count") +
    theme_minimal(base_size = 10)
}

plot_depth_distribution <- function(tab, title, subtitle, ylab, fill = "#2ea043") {
  if (nrow(tab) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No depth data", size = 4) +
        theme_void() +
        labs(title = title, subtitle = subtitle)
    )
  }
  tab <- tab %>% mutate(depth_band = order_depth_band(as.character(depth_band)))
  ggplot(tab, aes(.data$depth_band, .data$n)) +
    geom_col(fill = fill, width = 0.82, alpha = 0.92) +
    coord_flip() +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Depth band",
      y = ylab
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank())
}

write_datasets_pdf <- function(
    sightings_summary, loc_all, my_sight, my_obs,
    depth_sight, depth_obs, layers, joint, out_path,
    processed_summary = NULL
) {
  if (file.exists(out_path)) {
    removed <- tryCatch(
      { file.remove(out_path); TRUE },
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!removed) {
      out_path <- file.path(
        dirname(out_path),
        paste0(tools::file_path_sans_ext(basename(out_path)), "_latest.pdf")
      )
      message("Primary PDF locked (close in viewer); writing ", out_path)
    }
  }
  pdf(out_path, width = 11, height = 8.5, onefile = TRUE)

  # Cover: summary table as text
  grid::grid.newpage()
  grid::grid.text("Dataset summary — Master's project", y = 0.94, gp = grid::gpar(fontsize = 16, fontface = "bold"))
  grid::grid.text(paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")), y = 0.90, gp = grid::gpar(fontsize = 10))
  y <- 0.84
  for (i in seq_len(nrow(sightings_summary))) {
    r <- sightings_summary[i, ]
    tag <- if ("layer_type" %in% names(r) && !is.na(r$layer_type)) {
      if (grepl("Processed", r$layer_type, fixed = TRUE)) "processed" else "source"
    } else {
      "source"
    }
    line <- sprintf(
      "[%s] %s — %s sightings, %s observations, %s species",
      tag,
      r$dataset,
      format(r$sightings, big.mark = ","),
      format(r$observations, big.mark = ","),
      format(r$species, big.mark = ",")
    )
    grid::grid.text(line, y = y, just = "left", x = 0.05, gp = grid::gpar(fontsize = 9))
    y <- y - 0.038
    if (y < 0.14) {
      grid::grid.text("(continued in markdown)", y = 0.08, gp = grid::gpar(fontsize = 8, col = "grey40"))
      break
    }
  }
  grid::grid.text(
    paste(
      "Following pages: year x month heatmaps, calendar-month bars,",
      "and 50 m depth distributions (sightings & haul observations)."
    ),
    y = 0.12, gp = grid::gpar(fontsize = 9, col = "grey30")
  )

  if (!is.null(processed_summary) && nrow(processed_summary) > 0) {
    grid::grid.newpage()
    grid::grid.text(
      "Processed exports (data/processed/)",
      y = 0.94,
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    y <- 0.88
    for (i in seq_len(nrow(processed_summary))) {
      r <- processed_summary[i, ]
      line <- sprintf(
        "%s (%s rows): %s sightings · %s obs · %s spp · %s%% depth",
        r$dataset,
        format(r$file_rows, big.mark = ","),
        format(r$sightings, big.mark = ","),
        format(r$observations, big.mark = ","),
        format(r$species, big.mark = ","),
        ifelse(is.na(r$pct_depth_known), "—", as.character(r$pct_depth_known))
      )
      grid::grid.text(line, y = y, just = "left", x = 0.05, gp = grid::gpar(fontsize = 9))
      y <- y - 0.035
      if (nzchar(r$notes)) {
        grid::grid.text(
          paste0("  ", r$notes),
          y = y, just = "left", x = 0.06,
          gp = grid::gpar(fontsize = 8, col = "grey35")
        )
        y <- y - 0.03
      }
      if (y < 0.08) break
    }
  }

  for (L in layers) {
    ds <- L$dataset
    sc <- L$scope
    s_tab <- my_sight %>% filter(dataset == ds)
    o_tab <- my_obs %>% filter(dataset == ds)
    cal_s <- month_calendar(L$data, ds, sc, "sighting")
    cal_o <- month_calendar(L$data, ds, sc, "observation")

    print(plot_month_year_heatmap(
      s_tab, "sightings",
      paste(ds, "— sightings by year and month"),
      sc
    ))
    print(plot_month_calendar(
      cal_s,
      paste(ds, "— sightings pooled by calendar month")
    ))
    print(plot_month_year_heatmap(
      o_tab, "observations",
      paste(ds, "— observations (hauls) by year and month"),
      sc
    ))
    print(plot_month_calendar(
      cal_o,
      paste(ds, "— observations pooled by calendar month")
    ))

    dep_s <- depth_sight %>%
      filter(dataset == ds) %>%
      select(depth_band, n)
    dep_o <- depth_obs %>%
      filter(dataset == ds) %>%
      select(depth_band, n)

    print(plot_depth_distribution(
      dep_s,
      paste(ds, "— sightings by depth"),
      paste0(sc, " · 50 m bins · record depth"),
      "Sightings (records)",
      fill = "#58a6ff"
    ))
    print(plot_depth_distribution(
      dep_o,
      paste(ds, "— observations by depth"),
      paste0(
        sc,
        " · 50 m bins · median haul depth where haul key exists;",
        " otherwise record depth"
      ),
      "Observations (hauls or records)",
      fill = "#3fb950"
    ))
  }

  grid::grid.newpage()
  grid::grid.text(
    "All databases combined (summed)",
    y = 0.92,
    gp = grid::gpar(fontsize = 15, fontface = "bold")
  )
  tot_obs <- nrow(joint$units)
  tot_sight <- sum(joint$units$sightings, na.rm = TRUE)
  grid::grid.text(
    sprintf(
      "%s joint observations · %s records (sightings) · obs_id prefixed by database",
      format(tot_obs, big.mark = ","),
      format(tot_sight, big.mark = ",")
    ),
    y = 0.86,
    gp = grid::gpar(fontsize = 10)
  )

  print(plot_month_year_heatmap(
    joint$by_year_month,
    "observations",
    "All databases — observations by year and month (sum)",
    "Hauls/stations pooled; no cross-database deduplication"
  ))
  cal_joint <- joint$by_month %>%
    rename(n = observations) %>%
    select(month, month_name, n)
  print(plot_month_calendar(
    cal_joint,
    "All databases — observations by calendar month (all years)"
  ))
  print(plot_joint_stacked(
    joint$by_dataset_month,
    "month_name",
    "Observations by month — contribution per database",
    "Stacked counts sum to joint monthly total"
  ))

  if (nrow(joint$by_year) > 0) {
    p_yr <- ggplot(joint$by_year, aes(.data$year, .data$observations)) +
      geom_col(fill = "#d29922", width = 0.75) +
      scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
      labs(
        title = "All databases — observations by year (sum)",
        x = "Year",
        y = "Observations"
      ) +
      theme_minimal(base_size = 10)
    print(p_yr)
    print(plot_joint_stacked(
      joint$by_dataset_year,
      "year",
      "Observations by year — contribution per database",
      "Stacked counts sum to joint yearly total"
    ))
  }

  dep_j <- joint$by_depth %>% select(depth_band, n = observations)
  print(plot_depth_distribution(
    dep_j,
    "All databases — observations by depth (sum)",
    "Median depth per observation unit · 50 m bins",
    "Observations",
    fill = "#bc8cff"
  ))
  jd <- joint$by_dataset_depth %>%
    mutate(depth_band = order_depth_band(as.character(depth_band)))
  print(plot_joint_stacked(
    jd,
    "depth_band",
    "Observations by depth — contribution per database",
    "Stacked counts; haul median depth or record depth",
  ))

  dev.off()
  message("PDF saved: ", out_path)
  invisible(out_path)
}

count_summary <- function(df, dataset, scope, source_file, notes = "") {
  n_sightings <- nrow(df)
  n_obs <- if ("obs_key" %in% names(df) && any(!is.na(df$obs_key) & df$obs_key != "")) {
    df %>% filter(!is.na(obs_key), obs_key != "") %>% distinct(obs_key) %>% nrow()
  } else {
    n_sightings
  }
  n_species <- if ("scientificName" %in% names(df)) {
    n_distinct(df$scientificName)
  } else {
    NA_integer_
  }
  n_depth_known <- sum(is.finite(df$depth))
  tibble(
    dataset, scope, source_file, notes,
    sightings = n_sightings,
    observations = n_obs,
    species = n_species,
    sightings_with_depth = n_depth_known,
    pct_depth_known = if (n_sightings > 0) round(100 * n_depth_known / n_sightings, 1) else NA_real_
  )
}

message("Building column inventory...")
column_files <- bind_rows(
  file_cols("Combined_Menachem_Itai_corrected_May_2026_DBG.csv") %>% mutate(dataset = "Israel_trawl_MIT"),
  file_cols("data/raw/edelist_occurrence.txt") %>% mutate(dataset = "Edelist_OBIS"),
  file_cols("MEDITS/medits_raw_data.csv") %>% mutate(dataset = "MEDITS_catch"),
  file_cols("MEDITS/hauls_data.csv") %>% mutate(dataset = "MEDITS_hauls"),
  file_cols("MEDITS/taxa_medits.csv") %>% mutate(dataset = "MEDITS_taxa"),
  file_cols("data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv") %>%
    mutate(dataset = "Singer_Owens_iDigBio"),
  file_cols(GBIF_PAPER_CLEANED_CSV) %>%
    mutate(dataset = paste0(GBIF_PAPER_DATASET, " (processed)"))
)
proc_col <- bind_rows(lapply(processed_layers_manifest(), function(L) {
  file_cols(L$file) %>% mutate(dataset = L$dataset)
}))
column_files <- bind_rows(column_files, proc_col)
write_csv(column_files, file.path(OUT_DIR, "columns_by_file.csv"))

message("Loading standardized layers...")
layers <- list(
  list(
    dataset = "Israel_trawl_MIT",
    scope = "Israel coast box (lon 32–36.5°, lat 29.5–34.5°)",
    file = "Combined_Menachem_Itai_corrected_May_2026_DBG.csv",
    notes = "Trawl catch rows; observation = SAMPLE_DESC haul",
    data = load_mit("israel_box")
  ),
  list(
    dataset = "Edelist_OBIS",
    scope = "Actinopterygii + Israel box",
    file = "data/raw/edelist_occurrence.txt",
    notes = "Coords/depth/year parsed from id; observation = haul prefix before _Fish_",
    data = load_edelist("actinopterygii_israel")
  ),
  list(
    dataset = "MEDITS_raw",
    scope = "GSA 22 & 25 + Levant box (maps script 03)",
    file = "MEDITS/medits_raw_data.csv (+ hauls_data.csv)",
    notes = "Catch records with valid hauls; observation = country|area|vessel|year|haul",
    data = load_medits_catch("levant_gsa22_25")
  ),
  list(
    dataset = "Singer_Owens_iDigBio",
    scope = "Israel box, depth-resolved specimens",
    file = "data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv",
    notes = "Museum/iDigBio; each record treated as one observation",
    data = load_singer("israel_box")
  ),
  list(
    dataset = GBIF_PAPER_DATASET,
    scope = "Levant fish GBIF · paper pipeline (presence-only)",
    file = GBIF_PAPER_CLEANED_CSV,
    notes = "From `01_med_fish_3d_pa_maps.R` · CoordinateCleaner · sparse depth",
    data = load_gbif_paper()
  )
)

sightings_summary <- bind_rows(lapply(layers, function(L) {
  count_summary(L$data, L$dataset, L$scope, L$file, L$notes)
}))
write_csv(sightings_summary, file.path(OUT_DIR, "sightings_by_database_survey.csv"))

message("Processed folder exports...")
processed_manifest <- processed_layers_manifest()
processed_layers <- lapply(processed_manifest, function(L) {
  L$data <- L$loader()
  L
})
processed_discovered <- list.files(
  "data/processed", pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE
)
processed_known <- vapply(processed_manifest, function(L) L$file, character(1))
processed_orphans <- setdiff(processed_discovered, processed_known)
if (length(processed_orphans) > 0) {
  for (f in processed_orphans) {
    base <- tools::file_path_sans_ext(basename(f))
    processed_layers[[length(processed_layers) + 1]] <- list(
      dataset = paste0("proc_", base),
      scope = "Unlisted processed CSV (auto-included)",
      file = f,
      notes = "Add to processed_layers_manifest() for a clearer description",
      data = load_processed_occ_dwc(f)
    )
  }
  message("Also summarizing unlisted processed files: ", paste(basename(processed_orphans), collapse = ", "))
}

count_rows_file <- function(file) {
  if (!file.exists(file)) {
    return(NA_integer_)
  }
  as.integer(nrow(read_csv(file, show_col_types = FALSE, progress = FALSE)))
}

processed_summary <- bind_rows(lapply(processed_layers, function(L) {
  cs <- count_summary(L$data, L$dataset, L$scope, L$file, L$notes)
  cs %>% mutate(file_rows = count_rows_file(L$file), .before = 1)
}))
write_csv(processed_summary, file.path(OUT_DIR, "processed_files_summary.csv"))

sightings_summary_part1 <- bind_rows(
  sightings_summary %>% mutate(layer_type = "Survey / raw source"),
  processed_summary %>%
    select(-file_rows) %>%
    mutate(layer_type = "Processed (`data/processed/`)")
)
write_csv(sightings_summary_part1, file.path(OUT_DIR, "sightings_by_database.csv"))

processed_loc <- bind_rows(lapply(processed_layers, function(L) {
  loc_summary(L$data, L$dataset, L$scope)
}))
write_csv(processed_loc, file.path(OUT_DIR, "processed_location_extent.csv"))

processed_depth <- bind_rows(lapply(processed_layers, function(L) {
  depth_tables(L$data, L$dataset, L$scope, "sighting")
}))
write_csv(processed_depth, file.path(OUT_DIR, "processed_depth_distribution.csv"))

combined_path <- "data/processed/med_fish_occurrences_combined.csv"
if (file.exists(combined_path)) {
  write_csv(
    summarize_processed_by_source(combined_path),
    file.path(OUT_DIR, "processed_med_occurrences_by_source.csv")
  )
}

loc_all <- bind_rows(lapply(layers, function(L) {
  loc_summary(L$data, L$dataset, L$scope)
}))
write_csv(loc_all, file.path(OUT_DIR, "location_extent_by_database.csv"))

depth_sightings <- bind_rows(lapply(layers, function(L) {
  depth_tables(L$data, L$dataset, L$scope, "sighting")
}))
depth_obs <- bind_rows(lapply(layers, function(L) {
  depth_tables(L$data, L$dataset, L$scope, "observation")
}))
write_csv(depth_sightings, file.path(OUT_DIR, "depth_distribution_sightings.csv"))
write_csv(depth_obs, file.path(OUT_DIR, "depth_distribution_observations.csv"))

message("Year x month distributions...")
my_sightings <- bind_rows(lapply(layers, function(L) {
  month_year_sightings(L$data, L$dataset, L$scope)
}))
my_observations <- bind_rows(lapply(layers, function(L) {
  month_year_observations(L$data, L$dataset, L$scope)
}))
write_csv(my_sightings, file.path(OUT_DIR, "month_year_sightings.csv"))
write_csv(my_observations, file.path(OUT_DIR, "month_year_observations.csv"))

message("Joint (all databases) observation totals...")
joint <- build_joint_observations(layers)
write_csv(joint$by_year_month, file.path(OUT_DIR, "joint_observations_year_month.csv"))
write_csv(joint$by_month, file.path(OUT_DIR, "joint_observations_by_month.csv"))
write_csv(joint$by_year, file.path(OUT_DIR, "joint_observations_by_year.csv"))
write_csv(joint$by_depth, file.path(OUT_DIR, "joint_observations_by_depth.csv"))
write_csv(
  joint$by_dataset_year_month,
  file.path(OUT_DIR, "joint_observations_by_database_year_month.csv")
)
write_csv(
  joint$units,
  file.path(OUT_DIR, "joint_observation_units.csv")
)
joint_totals <- tibble(
  metric = c("observations", "sightings_on_units", "datasets"),
  value = c(
    nrow(joint$units),
    sum(joint$units$sightings, na.rm = TRUE),
    n_distinct(joint$units$dataset)
  )
)
write_csv(joint_totals, file.path(OUT_DIR, "joint_summary_totals.csv"))

joint_sight <- build_joint_sightings(layers)
write_csv(joint_sight$by_year_month, file.path(OUT_DIR, "joint_sightings_year_month.csv"))
write_csv(joint_sight$by_month, file.path(OUT_DIR, "joint_sightings_by_month.csv"))
write_csv(joint_sight$by_year, file.path(OUT_DIR, "joint_sightings_by_year.csv"))
write_csv(joint_sight$by_depth, file.path(OUT_DIR, "joint_sightings_by_depth.csv"))
write_csv(
  joint_sight$by_dataset_year_month,
  file.path(OUT_DIR, "joint_sightings_by_database_year_month.csv")
)

combined_pdf <- file.path(OUT_DIR, "DATASETS_COMBINED.pdf")
write_combined_pdf(sightings_summary, joint, joint_sight, combined_pdf)

pdf_path <- file.path(OUT_DIR, "DATASETS_SUMMARY.pdf")
write_datasets_pdf(
  sightings_summary_part1, loc_all, my_sightings, my_observations,
  depth_sightings, depth_obs, layers, joint, pdf_path,
  processed_summary = NULL
)
message("Wrote ", pdf_path)
message("Wrote ", combined_pdf)

# Ancillary: full-file row counts (unfiltered)
count_rows <- function(file) {
  if (!file.exists(file)) {
    return(NA_integer_)
  }
  if (str_detect(file, "occurrence\\.txt$")) {
    as.integer(nrow(read_tsv(file, show_col_types = FALSE, guess_max = 50000)))
  } else {
    as.integer(nrow(read_csv(file, show_col_types = FALSE, progress = FALSE)))
  }
}

raw_counts <- tibble(
  file = c(
    "Combined_Menachem_Itai_corrected_May_2026_DBG.csv",
    "data/raw/edelist_occurrence.txt",
    "MEDITS/medits_raw_data.csv",
    "MEDITS/hauls_data.csv",
    "data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv"
  ),
  description = c(
    "MIT trawl (all rows)",
    "Edelist OBIS export (all rows)",
    "MEDITS catch (all rows)",
    "MEDITS haul metadata (validity not filtered here)",
    "Singer & Owens Zenodo (all rows)"
  ),
  rows_in_file = vapply(file, count_rows, integer(1))
)
write_csv(raw_counts, file.path(OUT_DIR, "raw_file_row_counts.csv"))

processed_counts <- processed_summary %>%
  transmute(
    file = source_file,
    dataset = dataset,
    description = scope,
    rows_in_file = file_rows
  )
write_csv(processed_counts, file.path(OUT_DIR, "processed_file_row_counts.csv"))

md_path <- file.path(OUT_DIR, "DATASETS_SUMMARY.md")
md_buf <- character()

write_lines <- function(...) {
  lines <- unlist(lapply(list(...), function(x) if (length(x) > 1) x else as.character(x)))
  md_buf <<- c(md_buf, lines)
}

write_lines(
  "# Dataset summary (Master's project)",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "## Definitions",
  "",
  "- **Sightings**: individual catch/specimen **records** (fish rows).",
  "- **Observations**: unique **hauls/stations** where `obs_key` exists;",
  "  for Singer/GBIF each sighting is also one observation.",
  "- **Depth bands**: `50 m` bins (0–50, 50–100, …, 500+), plus `depth_unknown`.",
  "- **Observation depths**: median depth per haul (`obs_key`).",
  "",
  "## Sightings and observations (part 1 — sources + processed)",
  "",
  "Survey/raw layers (analysis scope) plus pipeline exports under `data/processed/`.",
  "Joint totals below use **survey / raw source** rows only (no double-counting).",
  "",
  "| Layer | Dataset | Scope | Sightings | Observations | Species | % depth | Source |",
  "|-------|---------|-------|----------:|-------------:|--------:|--------:|--------|"
)
for (i in seq_len(nrow(sightings_summary_part1))) {
  r <- sightings_summary_part1[i, ]
  write_lines(sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | `%s` |",
    r$layer_type, r$dataset, r$scope,
    format(r$sightings, big.mark = ","),
    format(r$observations, big.mark = ","),
    ifelse(is.na(r$species), "—", format(r$species, big.mark = ",")),
    ifelse(is.na(r$pct_depth_known), "—", paste0(r$pct_depth_known, "%")),
    r$source_file
  ))
}

write_lines(
  "",
  "## Raw file sizes (unfiltered)",
  "",
  "| File | Rows |",
  "|------|-----:|"
)
for (i in seq_len(nrow(raw_counts))) {
  write_lines(sprintf(
    "| `%s` | %s |",
    raw_counts$file[i],
    ifelse(is.na(raw_counts$rows_in_file[i]), "missing",
           format(raw_counts$rows_in_file[i], big.mark = ","))
  ))
}

write_lines(
  "",
  "## Processed files (`data/processed/`) — detail",
  "",
  "Same exports as **Processed** rows in [part 1](#sightings-and-observations-part-1--sources--processed).",
  "File row counts and breakdowns below.",
  ""
)

write_lines(
  "| File | Rows on disk |",
  "|------|-------------:|"
)
for (i in seq_len(nrow(processed_counts))) {
  write_lines(sprintf(
    "| `%s` | %s |",
    processed_counts$file[i],
    format(processed_counts$rows_in_file[i], big.mark = ",")
  ))
}

by_src_path <- file.path(OUT_DIR, "processed_med_occurrences_by_source.csv")
if (file.exists(by_src_path)) {
  by_src <- read_csv(by_src_path, show_col_types = FALSE)
  if (nrow(by_src) > 0) {
    write_lines(
      "",
      "### `med_fish_occurrences_combined.csv` — rows by `source`",
      "",
      "| source | occurrenceStatus | rows |",
      "|--------|------------------|-----:|"
    )
    for (j in seq_len(nrow(by_src))) {
      write_lines(sprintf(
        "| %s | %s | %s |",
        by_src$source[j], by_src$occurrenceStatus[j],
        format(by_src$rows[j], big.mark = ",")
      ))
    }
  }
}

if (nrow(processed_loc) > 0) {
  write_lines(
    "",
    "### Processed files — spatial extent",
    "",
    "| Dataset | N with coords | Lon min–max | Lat min–max |",
    "|---------|--------------:|-------------|-------------|"
  )
  for (i in seq_len(nrow(processed_loc))) {
    r <- processed_loc[i, ]
    if (r$n_with_coordinates == 0) next
    write_lines(sprintf(
      "| %s | %s | %.3f – %.3f | %.3f – %.3f |",
      r$dataset, format(r$n_with_coordinates, big.mark = ","),
      r$lon_min, r$lon_max, r$lat_min, r$lat_max
    ))
  }
}

write_lines(
  "",
  "CSVs: `processed_files_summary.csv`, `processed_file_row_counts.csv`,",
  "`processed_location_extent.csv`, `processed_depth_distribution.csv`.",
  ""
)

write_lines(
  "",
  "## Spatial extent (records with coordinates)",
  "",
  "| Dataset | N | Lon min–max | Lat min–max | Centroid | 0.5° grid cells |",
  "|---------|----:|-------------|-------------|----------|----------------:|"
)
for (i in seq_len(nrow(loc_all))) {
  r <- loc_all[i, ]
  if (r$n_with_coordinates == 0) next
  write_lines(sprintf(
    "| %s | %s | %.3f – %.3f | %.3f – %.3f | %.3f, %.3f | %s |",
    r$dataset, format(r$n_with_coordinates, big.mark = ","),
    r$lon_min, r$lon_max, r$lat_min, r$lat_max,
    r$centroid_lon, r$centroid_lat,
    format(r$n_grid_cells_0p5deg, big.mark = ",")
  ))
}

write_lines(
  "",
  "## Columns by file",
  "",
  "Full list: `columns_by_file.csv` (one row per column).",
  ""
)
for (ds in unique(column_files$dataset)) {
  cols <- column_files %>% filter(dataset == ds) %>% arrange(column_index)
  write_lines(paste0("### ", ds, " (`", cols$file[1], "`)"))
  write_lines("")
  write_lines(paste(cols$column_name, collapse = ", "))
  write_lines("")
}

write_lines(
  "## Temporal coverage (year x month)",
  "",
  "- Sightings: `month_year_sightings.csv`",
  "- Observations (hauls): `month_year_observations.csv`",
  "- **Joint (all DBs):** `joint_observations_year_month.csv`,",
  "  `joint_observations_by_month.csv`, `joint_observations_by_year.csv`,",
  "  `joint_observations_by_depth.csv`, `joint_summary_totals.csv`",
  "- Per-database report: `DATASETS_SUMMARY.pdf`",
  "- **All databases combined:** `DATASETS_COMBINED.pdf`",
  "- Joint sightings CSVs: `joint_sightings_*.csv`",
  ""
)

write_lines(
  "### All databases combined — observation counts",
  "",
  sprintf(
    "- **%s** joint observations · **%s** fish records on those units",
    format(nrow(joint$units), big.mark = ","),
    format(sum(joint$units$sightings, na.rm = TRUE), big.mark = ",")
  ),
  "",
  "| Year | Month | Observations |",
  "|-----:|------:|-------------:|"
)
for (j in seq_len(nrow(joint$by_year_month))) {
  r <- joint$by_year_month[j, ]
  write_lines(sprintf(
    "| %s | %s (%s) | %s |",
    r$year, r$month, MONTH_LAB[r$month],
    format(r$observations, big.mark = ",")
  ))
}
write_lines("")

for (ds in unique(my_sightings$dataset)) {
  sub <- my_sightings %>% filter(dataset == ds) %>% arrange(year, month)
  write_lines(paste0("### ", ds, " — sightings by year/month"))
  write_lines("")
  if (nrow(sub) == 0) {
    write_lines("_No month/year metadata in scope._")
  } else {
    write_lines("| Year | Month | Sightings |", "|-----:|------:|----------:|")
    for (j in seq_len(nrow(sub))) {
      write_lines(sprintf(
        "| %s | %s (%s) | %s |",
        sub$year[j], sub$month[j], MONTH_LAB[sub$month[j]],
        format(sub$sightings[j], big.mark = ",")
      ))
    }
  }
  write_lines("")
}

write_lines(
  "## Depth distributions",
  "",
  "- Sightings: `depth_distribution_sightings.csv`",
  "- Observations (haul median depth): `depth_distribution_observations.csv`",
  "",
  "### Observation counts by depth band (pooled datasets in analysis scope)",
  ""
)
for (ds in unique(depth_obs$dataset)) {
  sub <- depth_obs %>% filter(dataset == ds) %>% arrange(depth_band)
  write_lines(paste0("**", ds, "**"))
  write_lines("")
  write_lines("| Depth band | Observations |", "|------------|-------------:|")
  for (j in seq_len(nrow(sub))) {
    write_lines(sprintf("| %s | %s |", sub$depth_band[j], format(sub$n[j], big.mark = ",")))
  }
  write_lines("")
}

write_lines(
  "## Notes",
  "",
  "- **MEDITS** in maps uses GSA 22 & 25 inside the Levant box;",
  "  the raw catch file covers the full MEDITS programme.",
  "- **GBIF** for analysis is **`gbif_fish_levant_cleaned.csv`** (bbox + CoordinateCleaner);",
  "  mostly lacks depth (`pct_depth_known` low). Do not use `gbif_fish_levant.rds`.",
  "- **Part 1** lists processed exports alongside survey sources;",
  "  `proc_gbif_paper_cleaned` = same cleaned GBIF; `proc_med_fish_occurrences_combined` adds Edelist/Michmoret.",
  "- Regenerate: `source(\"06_dataset_summary.R\")` (markdown, CSV, PDF).",
  ""
)

write_md_ok <- tryCatch(
  {
    writeLines(md_buf, md_path, useBytes = TRUE)
    TRUE
  },
  error = function(e) {
    alt <- file.path(
      OUT_DIR,
      paste0(tools::file_path_sans_ext(basename(md_path)), "_latest.md")
    )
    message("Could not write ", md_path, " (", conditionMessage(e), "); trying ", alt)
    writeLines(md_buf, alt, useBytes = TRUE)
    FALSE
  }
)
if (!write_md_ok) {
  md_path <- file.path(OUT_DIR, "DATASETS_SUMMARY_latest.md")
}

message("Wrote ", md_path)
message("Outputs in ", OUT_DIR, " (incl. DATASETS_SUMMARY.pdf)")
