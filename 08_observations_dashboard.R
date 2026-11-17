# =============================================================================
# Observations dashboard (HTML + figures): big picture -> detail
# Output: dashboard/index.html, dashboard/figures/*.png, dashboard/*.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("gbif_paper.R")

OUT_DIR <- "dashboard"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

BBOX_ISRAEL <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 34.5)
BBOX_MEDITS <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 31.0, lat_max = 35.5)
BBOX_MED <- list(lon_min = -6.0, lon_max = 37.0, lat_min = 30.0, lat_max = 46.0)
DEPTH_BIN_M <- 50L
FOCUS_FIVE <- c(
  "Pagellus erythrinus",
  "Boops boops",
  "Spicara maena",
  "Upeneus moluccensis",
  "Trachurus mediterraneus"
)

# Israeli Hebrew common names (fisheries / field use)
FOCUS_FIVE_HEBREW <- c(
  "Pagellus erythrinus" = "\u05e4\u05d2\u05dc",
  "Boops boops" = "\u05d1\u05e7\u05dc\u05d4",
  "Spicara maena" = "\u05e1\u05e4\u05e8\u05d9\u05d3",
  "Upeneus moluccensis" = "\u05e1\u05e0\u05e4\u05d9\u05e8\u05d5\u05df",
  "Trachurus mediterraneus" = "\u05e1\u05e8\u05d2\u05d5\u05e1 \u05de\u05e6\u05d5\u05d9"
)

species_display <- function(latin) {
  h <- unname(FOCUS_FIVE_HEBREW[latin])
  if (length(h) == 0 || is.na(h) || !nzchar(h)) {
    return(latin)
  }
  paste0(latin, " (", h, ")")
}

SURVEY_MAP_DATASETS <- c("MEDITS_raw", "Israel_trawl_MIT", "Edelist_OBIS")
PA_EXCLUDED_DATASETS <- c("GBIF_paper_cleaned", "Singer_Owens_iDigBio")

filter_boops_surveys <- function(df) {
  df %>% filter(dataset %in% SURVEY_MAP_DATASETS)
}

theme_dash <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

in_bbox <- function(df, box) {
  df %>%
    filter(
      decimalLongitude >= box$lon_min, decimalLongitude <= box$lon_max,
      decimalLatitude >= box$lat_min, decimalLatitude <= box$lat_max
    )
}

medits_bottom_temp_c <- function(beg, end) {
  b <- suppressWarnings(as.numeric(beg))
  e <- suppressWarnings(as.numeric(end))
  bad <- function(x) !is.finite(x) | x <= 0
  b[bad(b)] <- NA_real_
  e[bad(e)] <- NA_real_
  out <- (b + e) / 2
  out[!is.finite(out) | out <= 0] <- NA_real_
  out
}

medits_dd <- function(x) {
  x <- as.numeric(x)
  deg <- floor(x / 100)
  deg + (x - deg * 100) / 60
}

depth_band50 <- function(d) {
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

order_depth_band50_factor <- function(x) {
  x <- as.character(x)
  incr <- sprintf("%d\u2013%d m", seq(0, 450, DEPTH_BIN_M), seq(DEPTH_BIN_M, 500, DEPTH_BIN_M))
  ord <- c(incr, "500+ m", "<0 m", "depth_unknown")
  lv <- ord[ord %in% unique(x)]
  lv <- c(lv, setdiff(unique(x), lv))
  factor(x, levels = lv)
}

gbif_is_present <- gbif_paper_is_present

gbif_match_key <- function(lon, lat, event_date) {
  paste0(
    round(as.numeric(lon), 6), "|",
    round(as.numeric(lat), 6), "|",
    substr(as.character(event_date), 1, 10)
  )
}

load_gbif_levant_records <- function(box = BBOX_ISRAEL) {
  raw <- load_gbif_paper_cleaned()
  if (nrow(raw) == 0) {
    message("GBIF survey layer: missing cleaned export (run 01_med_fish_3d_pa_maps.R)")
    return(tibble())
  }
  out <- raw %>%
    transmute(
      dataset = GBIF_PAPER_DATASET,
      obs_key = gbif_match_key(decimalLongitude, decimalLatitude, eventDate),
      year = suppressWarnings(as.integer(substr(eventDate, 1, 4))),
      month = suppressWarnings(as.integer(substr(eventDate, 6, 7))),
      depth = suppressWarnings(as.numeric(depth)),
      decimalLatitude = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      bottom_temp_c = NA_real_,
      sampling_protocol = NA_character_,
      basis_of_record = NA_character_,
      eventDate = as.character(eventDate)
    ) %>%
    distinct(obs_key, .keep_all = TRUE) %>%
    in_bbox(box) %>%
    select(-eventDate)
}

load_singer_survey_records <- function(box = BBOX_ISRAEL) {
  cache <- "data/processed/singer_owens_israel_box.csv"
  if (file.exists(cache)) {
    return(
      read_csv(cache, show_col_types = FALSE, progress = FALSE) %>%
        transmute(
          dataset = "Singer_Owens_iDigBio",
          obs_key = paste0(
            round(as.numeric(decimalLongitude), 6), "|",
            round(as.numeric(decimalLatitude), 6), "|",
            suppressWarnings(as.integer(year))
          ),
          year = suppressWarnings(as.integer(year)),
          month = NA_integer_,
          depth = suppressWarnings(as.numeric(depth)),
          decimalLatitude = as.numeric(decimalLatitude),
          decimalLongitude = as.numeric(decimalLongitude),
          bottom_temp_c = NA_real_,
          sampling_protocol = NA_character_,
          basis_of_record = NA_character_
        ) %>%
        filter(is.finite(decimalLatitude), is.finite(decimalLongitude)) %>%
        in_bbox(box)
    )
  }

  raw_path <- "data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv"
  if (!file.exists(raw_path)) {
    return(tibble())
  }
  hdr <- names(read_csv(raw_path, n_max = 0, show_col_types = FALSE))
  pick <- function(cands) cands[cands %in% hdr][1]
  lat_c <- pick(c("iDigBio_dwc_decimalLatitude", "dwc_decimalLatitude"))
  lon_c <- pick(c("iDigBio_dwc_decimalLongitude", "dwc_decimalLongitude"))
  dep_c <- pick(c("SOcurated_depth", "depth"))
  yr_c <- pick(c("iDigBio_dwc_year", "dwc_year", "year"))
  if (any(is.na(c(lat_c, lon_c, dep_c)))) {
    return(tibble())
  }
  read_csv(raw_path, show_col_types = FALSE,
           col_select = all_of(unique(na.omit(c(lat_c, lon_c, dep_c, yr_c)))),
           progress = FALSE) %>%
    transmute(
      dataset = "Singer_Owens_iDigBio",
      obs_key = paste0(
        round(as.numeric(.data[[lon_c]]), 6), "|",
        round(as.numeric(.data[[lat_c]]), 6), "|",
        if (!is.na(yr_c)) suppressWarnings(as.integer(.data[[yr_c]])) else NA_integer_
      ),
      year = if (!is.na(yr_c)) suppressWarnings(as.integer(.data[[yr_c]])) else NA_integer_,
      month = NA_integer_,
      depth = as.numeric(.data[[dep_c]]),
      decimalLatitude = as.numeric(.data[[lat_c]]),
      decimalLongitude = as.numeric(.data[[lon_c]]),
      bottom_temp_c = NA_real_,
      sampling_protocol = NA_character_,
      basis_of_record = NA_character_
    ) %>%
    filter(is.finite(decimalLatitude), is.finite(decimalLongitude), is.finite(depth)) %>%
    in_bbox(box)
}

load_observation_records <- function() {
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
        decimalLongitude = as.numeric(LONGITUDE),
        bottom_temp_c = NA_real_,
        sampling_protocol = NA_character_,
        basis_of_record = NA_character_
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
        decimalLongitude = as.numeric(meta[, 3]),
        bottom_temp_c = NA_real_,
        sampling_protocol = NA_character_,
        basis_of_record = NA_character_
      ) %>%
      filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
      in_bbox(BBOX_ISRAEL)
  }

  if (file.exists("MEDITS/medits_raw_data.csv") && file.exists("MEDITS/hauls_data.csv")) {
    hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE) %>%
      filter(validity == "V") %>%
      mutate(
        obs_key = paste(country, area, vessel, year, haul_number, sep = "|"),
        haul_month = as.integer(month),
        decimalLatitude = medits_dd(shooting_latitude),
        decimalLongitude = medits_dd(shooting_longitude),
        depth = (as.numeric(shooting_depth) + as.numeric(hauling_depth)) / 2,
        bottom_temp_c = medits_bottom_temp_c(
          bottom_temperature_beginning, bottom_temperature_end
        )
      )
    catch <- read_csv("MEDITS/medits_raw_data.csv", show_col_types = FALSE,
                      locale = locale(encoding = "Latin1")) %>%
      filter(nbtot > 0) %>%
      distinct(country, area, vessel, year, haul_number, .keep_all = FALSE)
    parts[[length(parts) + 1]] <- catch %>%
      inner_join(
        hauls %>% select(country, area, vessel, year, haul_number,
                         haul_month, decimalLatitude, decimalLongitude,
                         depth, obs_key, bottom_temp_c),
        by = c("country", "area", "vessel", "year", "haul_number")
      ) %>%
      transmute(
        dataset = "MEDITS_raw",
        obs_key, year = as.integer(year), month = haul_month,
        depth, decimalLatitude, decimalLongitude, bottom_temp_c,
        sampling_protocol = NA_character_,
        basis_of_record = NA_character_,
        gsa = as.integer(area)
      ) %>%
      filter(gsa %in% c(22L, 25L)) %>%
      in_bbox(BBOX_MEDITS) %>%
      select(-gsa)
  }

  singer <- load_singer_survey_records(BBOX_ISRAEL)
  if (nrow(singer) > 0) {
    parts[[length(parts) + 1]] <- singer
  }

  gbif <- load_gbif_levant_records(BBOX_ISRAEL)
  if (nrow(gbif) > 0) {
    parts[[length(parts) + 1]] <- gbif
  }

  bind_rows(parts)
}

norm_species <- function(x) {
  x %>% as.character() %>% str_trim() %>%
    str_remove("\\s*\\([^)]*\\)$") %>% str_squish()
}

slug_species <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_sub(1, 60)
}

in_bbox_med <- function(df) {
  df %>%
    filter(
      decimalLongitude >= BBOX_MED$lon_min, decimalLongitude <= BBOX_MED$lon_max,
      decimalLatitude >= BBOX_MED$lat_min, decimalLatitude <= BBOX_MED$lat_max
    )
}

load_boops_occurrences <- function() {
  target <- norm_species("Boops boops")
  parts <- list()

  mit_path <- "Combined_Menachem_Itai_corrected_May_2026_DBG.csv"
  if (file.exists(mit_path)) {
    parts[[length(parts) + 1]] <- read_csv(mit_path, show_col_types = FALSE) %>%
      transmute(
        scientificName = norm_species(ID_SPECIES),
        decimalLatitude = as.numeric(LATITUDE),
        decimalLongitude = as.numeric(LONGITUDE),
        year = as.integer(YEAR),
        dataset = "Israel_trawl_MIT"
      ) %>%
      filter(scientificName == target) %>%
      in_bbox_med()
  }

  if (file.exists("MEDITS/medits_raw_data.csv") && file.exists("MEDITS/hauls_data.csv")) {
    hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE) %>%
      filter(validity == "V") %>%
      mutate(
        country = as.character(country), area = as.integer(area),
        vessel = as.character(vessel), year = as.integer(year),
        haul_number = as.numeric(haul_number),
        decimalLatitude = medits_dd(shooting_latitude),
        decimalLongitude = medits_dd(shooting_longitude)
      ) %>%
      select(country, area, vessel, year, haul_number, decimalLatitude, decimalLongitude)

    catch <- read_csv("MEDITS/medits_raw_data.csv", show_col_types = FALSE,
                      locale = locale(encoding = "Latin1")) %>%
      filter(nbtot > 0) %>%
      mutate(
        country = as.character(country), area = as.integer(area),
        vessel = as.character(vessel), year = as.integer(year),
        haul_number = as.numeric(haul_number),
        scientificName = norm_species(paste0(genus, species))
      )

    taxa_path <- "MEDITS/taxa_medits.csv"
    if (file.exists(taxa_path)) {
      taxa <- read_csv(taxa_path, show_col_types = FALSE, locale = locale(encoding = "Latin1"))
      code_col <- names(taxa)[str_detect(names(taxa), "^MEDITS")][1]
      name_col <- names(taxa)[str_detect(names(taxa), "Scientific Name")][1]
      if (!is.na(code_col) && !is.na(name_col)) {
        taxa <- taxa %>%
          mutate(
            medits_code = .data[[code_col]],
            sci = norm_species(.data[[name_col]])
          ) %>%
          select(medits_code, sci)
        catch <- catch %>%
          mutate(medits_code = paste0(genus, species)) %>%
          left_join(taxa, by = "medits_code") %>%
          mutate(scientificName = coalesce(sci, scientificName))
      }
    }

    parts[[length(parts) + 1]] <- catch %>%
      filter(scientificName == target) %>%
      inner_join(hauls, by = c("country", "area", "vessel", "year", "haul_number")) %>%
      transmute(scientificName, decimalLatitude, decimalLongitude, year, dataset = "MEDITS_raw") %>%
      in_bbox_med()
  }

  if (file.exists("data/raw/edelist_occurrence.txt")) {
    ed <- read_tsv("data/raw/edelist_occurrence.txt", show_col_types = FALSE, guess_max = 50000)
    id_col <- if ("id" %in% names(ed)) "id" else "occurrenceID"
    parts[[length(parts) + 1]] <- ed %>%
      mutate(scientificName = norm_species(scientificName)) %>%
      filter(scientificName == target) %>%
      mutate(
        parsed = str_match(.data[[id_col]], "^([0-9.]+)_([0-9.]+)_([0-9.]+)")
      ) %>%
      transmute(
        scientificName,
        decimalLatitude = as.numeric(parsed[, 2]),
        decimalLongitude = as.numeric(parsed[, 3]),
        year = suppressWarnings(as.integer(str_match(.data[[id_col]], "_([0-9]{4})_")[, 2])),
        dataset = "Edelist_OBIS"
      ) %>%
      filter(is.finite(decimalLatitude), is.finite(decimalLongitude)) %>%
      in_bbox_med()
  }

  gbif_paper <- load_gbif_paper_cleaned()
  if (nrow(gbif_paper) > 0) {
    parts[[length(parts) + 1]] <- gbif_paper %>%
      mutate(scientificName = norm_species(scientificName)) %>%
      filter(scientificName == target) %>%
      transmute(
        scientificName,
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        year = suppressWarnings(as.integer(substr(eventDate, 1, 4))),
        dataset = GBIF_PAPER_DATASET
      ) %>%
      filter(is.finite(decimalLatitude), is.finite(decimalLongitude)) %>%
      in_bbox_med()
  }

  if (file.exists("data/processed/singer_med_focus5.csv")) {
    parts[[length(parts) + 1]] <- read_csv(
      "data/processed/singer_med_focus5.csv", show_col_types = FALSE, progress = FALSE
    ) %>%
      mutate(scientificName = norm_species(scientificName)) %>%
      filter(scientificName == target) %>%
      transmute(
        scientificName,
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        year = suppressWarnings(as.integer(year)),
        dataset = "Singer_Owens_iDigBio"
      ) %>%
      filter(is.finite(decimalLatitude), is.finite(decimalLongitude)) %>%
      in_bbox_med()
  }

  bind_rows(parts) %>%
    distinct(decimalLatitude, decimalLongitude, year, dataset, .keep_all = TRUE)
}

plot_boops_by_era <- function(df, bin_years, fig_name, surveys_only = FALSE) {
  df <- df %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude), !is.na(year), year >= 1900L)
  if (surveys_only) {
    df <- filter_boops_surveys(df)
  }
  if (nrow(df) == 0) {
    return(list(path = NULL, n = 0L))
  }

  df <- df %>%
    mutate(
      era_start = bin_years * floor(year / bin_years),
      era_label = paste0(era_start, "\u2013", era_start + bin_years - 1L)
    ) %>%
    mutate(era_label = factor(era_label, levels = sort(unique(era_label))))

  if (nrow(df) > 12000) {
    df <- df %>% group_by(era_label) %>% slice_sample(n = min(2500, n()), .groups = "drop")
  }

  n_eras <- n_distinct(df$era_label)
  scope_lab <- if (surveys_only) {
    "Surveys only (MIT + Edelist + MEDITS; no GBIF)"
  } else {
    "All sources"
  }
  p <- ggplot(df, aes(decimalLongitude, decimalLatitude, colour = dataset)) +
    geom_point(alpha = 0.55, size = 1.1) +
    facet_wrap(~era_label, ncol = min(3L, n_eras)) +
    coord_fixed(
      xlim = c(BBOX_MED$lon_min, BBOX_MED$lon_max),
      ylim = c(BBOX_MED$lat_min, BBOX_MED$lat_max),
      expand = FALSE
    ) +
    scale_colour_brewer(palette = "Set2", name = "Database") +
    labs(
      title = paste0("Boops boops \u2014 presence by ", bin_years, "-year period"),
      subtitle = paste0(scope_lab, " \u00b7 Mediterranean extent"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_dash() +
    theme(strip.text = element_text(face = "bold"))

  h <- max(6, 2.8 + 2.6 * ceiling(n_eras / 3))
  path <- save_fig(p, fig_name, w = 12, h = h)
  list(path = path, n = nrow(df), bin_years = bin_years, surveys_only = surveys_only)
}

sync_boops_pa_map <- function() {
  src <- file.path("primary_sdm", slug_species("Boops boops"), "map_2d_MED_PA_joint_surveys.png")
  if (!file.exists(src)) {
    return(NULL)
  }
  dest <- file.path(FIG_DIR, "07c_boops_pa_surveys.png")
  file.copy(src, dest, overwrite = TRUE)
  dest
}

load_processed_dashboard_records <- function() {
  parts <- list()

  load_dwc_processed <- function(path, ds_name) {
    if (!file.exists(path)) {
      return(tibble())
    }
    raw <- read_csv(path, show_col_types = FALSE, progress = FALSE)
    req <- c("decimalLatitude", "decimalLongitude")
    if (!all(req %in% names(raw))) {
      return(tibble())
    }
    raw %>%
      transmute(
        dataset = ds_name,
        obs_key = if ("station_id" %in% names(.)) {
          as.character(station_id)
        } else if ("obs_key" %in% names(.)) {
          as.character(obs_key)
        } else {
          NA_character_
        },
        year = if ("eventDate" %in% names(.)) {
          suppressWarnings(as.integer(substr(eventDate, 1, 4)))
        } else if ("year" %in% names(.)) {
          suppressWarnings(as.integer(year))
        } else {
          NA_integer_
        },
        month = if ("eventDate" %in% names(.)) {
          suppressWarnings(as.integer(substr(eventDate, 6, 7)))
        } else if ("month" %in% names(.)) {
          suppressWarnings(as.integer(month))
        } else {
          NA_integer_
        },
        depth = if ("depth" %in% names(.)) as.numeric(depth) else NA_real_,
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        bottom_temp_c = NA_real_,
        sampling_protocol = NA_character_,
        basis_of_record = NA_character_
      )
  }

  load_singer_processed <- function(path, ds_name) {
    if (!file.exists(path)) {
      return(tibble())
    }
    read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
      transmute(
        dataset = ds_name,
        obs_key = NA_character_,
        year = if ("year" %in% names(.)) suppressWarnings(as.integer(year)) else NA_integer_,
        month = NA_integer_,
        depth = if ("depth" %in% names(.)) as.numeric(depth) else NA_real_,
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        bottom_temp_c = NA_real_,
        sampling_protocol = NA_character_,
        basis_of_record = NA_character_
      )
  }

  parts[[length(parts) + 1]] <- load_dwc_processed(
    GBIF_PAPER_CLEANED_CSV,
    "proc_gbif_paper_cleaned"
  )
  parts[[length(parts) + 1]] <- load_dwc_processed(
    "data/processed/med_fish_occurrences_combined.csv",
    "proc_med_fish_occurrences_combined"
  )
  parts[[length(parts) + 1]] <- load_dwc_processed(
    "data/processed/pa_focus_species.csv",
    "proc_pa_focus_species"
  )
  parts[[length(parts) + 1]] <- load_singer_processed(
    "data/processed/singer_owens_israel_box.csv",
    "proc_singer_owens_israel_box"
  )
  parts[[length(parts) + 1]] <- load_singer_processed(
    "data/processed/singer_med_focus5.csv",
    "proc_singer_med_focus5"
  )

  bind_rows(parts)
}

sync_species_map_figures <- function() {
  dest <- file.path(FIG_DIR, "species_2d")
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  catalog <- list()
  for (sp in FOCUS_FIVE) {
    sp_dir <- file.path("primary_sdm", slug_species(sp))
    variants <- list(
      all = "map_2d_MED_all_datasets.png",
      panel = "map_2d_MED_panel.png",
      pa_surveys = "map_2d_MED_PA_joint_surveys.png"
    )
    entry <- list(label = sp, maps = list())
    for (nm in names(variants)) {
      src <- file.path(sp_dir, variants[[nm]])
      if (file.exists(src)) {
        dst_name <- paste0(slug_species(sp), "_", nm, ".png")
        dst <- file.path(dest, dst_name)
        file.copy(src, dst, overwrite = TRUE)
        entry$maps[[nm]] <- paste0("figures/species_2d/", dst_name)
      }
    }
    if (length(entry$maps) > 0) {
      catalog[[length(catalog) + 1]] <- entry
    }
  }
  catalog
}

collapse_observations <- function(rec) {
  rec %>%
    mutate(
      obs_id = if_else(
        !is.na(obs_key) & obs_key != "",
        paste0(dataset, "::", obs_key),
        NA_character_
      )
    ) %>%
    group_by(dataset) %>%
    mutate(
      obs_id = if_else(is.na(obs_id), paste0(dataset, "::occ::", row_number()), obs_id)
    ) %>%
    ungroup() %>%
    group_by(obs_id, dataset) %>%
    summarise(
      year = year[!is.na(year)][1],
      depth_m = median(depth, na.rm = TRUE),
      latitude = median(decimalLatitude, na.rm = TRUE),
      longitude = median(decimalLongitude, na.rm = TRUE),
      bottom_temp_c = median(bottom_temp_c, na.rm = TRUE),
      sampling_protocol = dplyr::first(na.omit(sampling_protocol)),
      basis_of_record = dplyr::first(na.omit(basis_of_record)),
      sightings = n(),
      .groups = "drop"
    ) %>%
    mutate(
      depth_band50 = depth_band50(depth_m),
      lon_bin = floor(longitude * 2) / 2,
      lat_bin = floor(latitude * 2) / 2,
      grid_id = paste0(lon_bin, "_", lat_bin)
    )
}

save_fig <- function(p, name, w = 10, h = 5.5) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 150, bg = "white")
  path
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

js_quote <- function(path_basename) {
  if (is.null(path_basename) || !nzchar(path_basename)) {
    return("null")
  }
  x <- paste0("figures/", path_basename)
  x <- gsub("\\\\", "\\\\\\\\", x, fixed = TRUE)
  x <- gsub("'", "\\'", x, fixed = TRUE)
  paste0("'", x, "'")
}

species_maps_to_js <- function(catalog) {
  if (length(catalog) == 0) {
    return("window.SPECIES_MAPS = {};")
  }
  js_str <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", x, fixed = TRUE)
    x <- gsub("'", "\\'", x, fixed = TRUE)
    paste0("'", x, "'")
  }
  lines <- vapply(
    catalog,
    function(e) {
      slug <- slug_species(e$label)
      map_bits <- vapply(
        names(e$maps),
        function(k) sprintf("%s: %s", js_str(k), js_str(e$maps[[k]])),
        character(1)
      )
      sprintf(
        "%s: { label: %s, maps: { %s } }",
        js_str(slug),
        js_str(species_display(e$label)),
        paste(map_bits, collapse = ", ")
      )
    },
    character(1)
  )
  paste0("window.SPECIES_MAPS = {", paste(lines, collapse = ",\n  "), "};")
}

write_html_dashboard <- function(obs, by_db, paths, stats, species_catalog, boops_era) {
  species_js <- species_maps_to_js(species_catalog)
  default_sp <- if (length(species_catalog)) slug_species(species_catalog[[1]]$label) else ""
  boops_pick <- function(x, y) {
    p <- boops_era[[x]][[y]]$path
    if (is.null(p)) NULL else basename(p)
  }
  boops_10_all <- boops_pick("all", "y10")
  boops_5_all <- boops_pick("all", "y5")
  boops_10_surv <- boops_pick("surveys", "y10")
  boops_5_surv <- boops_pick("surveys", "y5")
  boops_pa <- boops_era$pa_map

  html <- c(
    "<!DOCTYPE html>",
    "<html lang='en'><head>",
    "<meta charset='utf-8'/>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'/>",
    "<title>Observations dashboard</title>",
    "<style>",
    "body{font-family:Segoe UI,system-ui,sans-serif;margin:0;background:#f6f8fa;color:#1f2328}",
    "header{background:#0d1117;color:#f0f6fc;padding:24px 32px}",
    "nav{background:#fff;border-bottom:1px solid #d0d7de;padding:12px 32px;position:sticky;top:0}",
    "nav a{margin-right:16px;color:#0969da;text-decoration:none;font-size:14px}",
    "main{max-width:1100px;margin:0 auto;padding:24px 32px 64px}",
    "section{background:#fff;border:1px solid #d0d7de;border-radius:8px;padding:20px 24px;margin-bottom:24px}",
    "h2{margin:0 0 8px;font-size:1.35rem}",
    "p.sub{color:#656d76;margin:0 0 16px;font-size:14px}",
    ".kpis{display:flex;flex-wrap:wrap;gap:12px;margin:16px 0}",
    ".kpi{flex:1 1 160px;background:#f6f8fa;border:1px solid #d0d7de;border-radius:8px;padding:12px 16px}",
    ".kpi b{display:block;font-size:1.6rem}",
    "table{border-collapse:collapse;width:100%;font-size:14px}",
    "th,td{border:1px solid #d0d7de;padding:8px 10px;text-align:left}",
    "th{background:#f6f8fa}",
    "img.chart{max-width:100%;height:auto;border:1px solid #d0d7de;border-radius:6px}",
    ".pickers{display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin:12px 0 16px}",
    ".pickers label{font-size:14px;color:#656d76}",
    ".pickers select{padding:6px 10px;font-size:14px;border-radius:6px;border:1px solid #d0d7de}",
    "</style>",
    "<script>",
    species_js,
    "</script>",
    "</head><body>",
    "<header><h1>Observations dashboard</h1>",
    "<p>Survey &amp; cache data · Levant / Israel analysis scope · ",
    format(Sys.time(), "%Y-%m-%d %H:%M"), "</p></header>",
    "<nav>",
    "<a href='#overview'>Overview</a>",
    "<a href='#species-maps'>Species maps</a>",
    "<a href='#boops-era'>Boops eras</a>",
    "<a href='#temperature'>Temperature</a>",
    "<a href='#depth'>Depth</a>",
    "</nav><main>",

    "<section id='overview'>",
    "<h2>1. How many observations?</h2>",
    "<p class='sub'>Survey units: haul/station (<code>obs_key</code>), Singer &amp; Owens voucher, or",
    " GBIF present rows from <code>gbif_fish_levant_cleaned.csv</code> (01 pipeline: bbox + CoordinateCleaner).",
    " <strong>Processed</strong> exports listed separately; not in survey KPIs.</p>",
    "<div class='kpis'>",
    sprintf("<div class='kpi'><b>%s</b>Survey observation units</div>", format(stats$total, big.mark = ",")),
    sprintf("<div class='kpi'><b>%s</b>Sightings on survey units</div>", format(stats$sightings, big.mark = ",")),
    sprintf("<div class='kpi'><b>%s</b>Survey databases</div>", stats$datasets_survey),
    sprintf("<div class='kpi'><b>%s</b>Processed export rows</div>", format(stats$processed_rows, big.mark = ",")),
    "</div>",
    "<table><tr><th>Layer</th><th>Database</th><th>Observations</th><th>Sightings</th></tr>"
  )
  for (i in seq_len(nrow(by_db))) {
    html <- c(html, sprintf(
      "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
      by_db$layer[i],
      by_db$dataset[i],
      format(by_db$observations[i], big.mark = ","),
      format(by_db$sightings[i], big.mark = ",")
    ))
  }
  html <- c(html,
    "</table>",
    sprintf("<p><img class='chart' src='figures/%s' alt='by database'/></p>", basename(paths$by_db)),
    "</section>",

    "<section id='species-maps'>",
    "<h2>2. Focus species — Mediterranean 2D maps</h2>",
    "<p class='sub'>Maps from <code>04_med_2d_maps.R</code>: all sources, panels, or",
    " survey presence–absence (MIT + Edelist + MEDITS; no GBIF/Singer).</p>",
    if (length(species_catalog) == 0) {
      "<p><em>No species map PNGs found under <code>primary_sdm/</code>. Run <code>source('04_med_2d_maps.R')</code>.</em></p>"
    } else {
      c(
        "<div class='pickers'>",
        "<label>Species <select id='sp-species'>",
        paste(
          vapply(
            species_catalog,
            function(e) sprintf(
              "<option value='%s'>%s</option>",
              slug_species(e$label),
              html_escape(species_display(e$label))
            ),
            character(1)
          ),
          collapse = ""
        ),
        "</select></label>",
        "<label>View <select id='sp-variant'>",
        "<option value='all'>All datasets (joint)</option>",
        "<option value='panel'>Panel by dataset</option>",
        "<option value='pa_surveys'>Surveys P/A (excl. GBIF)</option>",
        "</select></label>",
        "</div>",
        sprintf("<p><img id='sp-map-img' class='chart' src='%s' alt='species map'/></p>",
                if (nchar(default_sp) && "all" %in% names(species_catalog[[1]]$maps)) {
                  species_catalog[[1]]$maps[["all"]]
                } else {
                  species_catalog[[1]]$maps[[1]]
                }),
        "<script>",
        "function refreshSpeciesMap(){",
        "const slug=document.getElementById('sp-species').value;",
        "const variant=document.getElementById('sp-variant').value;",
        "const entry=window.SPECIES_MAPS[slug];",
        "const img=document.getElementById('sp-map-img');",
        "if(!entry||!entry.maps){return;}",
        "let src=entry.maps[variant];",
        "if(!src){src=entry.maps.all||entry.maps.pa_surveys||Object.values(entry.maps)[0];}",
        "if(src){img.src=src;img.alt=entry.label+' Mediterranean map';}",
        "else if(variant==='pa_surveys'){img.alt='No survey P/A map for this species (run 04_med_2d_maps.R)';}",
        "}",
        "document.getElementById('sp-species').addEventListener('change',refreshSpeciesMap);",
        "document.getElementById('sp-variant').addEventListener('change',refreshSpeciesMap);",
        "</script>"
      )
    },
    "</section>",

    "<section id='boops-era'>",
    sprintf("<h2>3. %s — presence by time period</h2>", html_escape(species_display("Boops boops"))),
    "<p class='sub'>By era (all sources or surveys only) or joint survey presence–absence",
    " (MIT + Edelist + MEDITS; excludes GBIF &amp; Singer). Maps from <code>04_med_2d_maps.R</code>.</p>",
    if (is.null(boops_10_all) && is.null(boops_5_all) && is.null(boops_pa)) {
      "<p><em>No dated Boops boops records / PA map in scope.</em></p>"
    } else {
      boops_init <- if (!is.null(boops_10_all)) {
        boops_10_all
      } else if (!is.null(boops_5_all)) {
        boops_5_all
      } else {
        boops_pa
      }
      c(
        "<div class='pickers'>",
        "<label>Scope <select id='boops-scope'>",
        if (!is.null(boops_10_all) || !is.null(boops_5_all)) {
          "<option value='era_all' selected>All sources · by era</option>"
        } else "",
        if (!is.null(boops_10_surv) || !is.null(boops_5_surv)) {
          if (is.null(boops_10_all) && is.null(boops_5_all)) {
            "<option value='era_surveys' selected>Surveys only (no GBIF) · by era</option>"
          } else {
            "<option value='era_surveys'>Surveys only (no GBIF) · by era</option>"
          }
        } else "",
        if (!is.null(boops_pa)) {
          "<option value='pa_surveys'>Surveys P/A (excl. GBIF)</option>"
        } else "",
        "</select></label>",
        "<span id='boops-bin-wrap'>",
        "<label>Bin width <select id='boops-bin'>",
        if (!is.null(boops_10_all) || !is.null(boops_10_surv)) {
          "<option value='10' selected>10 years</option>"
        } else "",
        if (!is.null(boops_5_all) || !is.null(boops_5_surv)) {
          if (is.null(boops_10_all) && is.null(boops_10_surv)) {
            "<option value='5' selected>5 years</option>"
          } else {
            "<option value='5'>5 years</option>"
          }
        } else "",
        "</select></label>",
        "</span>",
        "</div>",
        sprintf(
          "<p><img id='boops-era-img' class='chart' src='figures/%s' alt='Boops map'/></p>",
          boops_init
        ),
        sprintf(
          paste0(
            "<script>",
            "var BOOPS_MAPS={",
            "era_all:{5:%s,10:%s},",
            "era_surveys:{5:%s,10:%s},",
            "pa_surveys:%s",
            "};",
            "function refreshBoopsMap(){",
            "var scope=document.getElementById('boops-scope').value,",
            "bin=document.getElementById('boops-bin').value,",
            "img=document.getElementById('boops-era-img'),",
            "wrap=document.getElementById('boops-bin-wrap');",
            "if(scope==='pa_surveys'){wrap.style.display='none';",
            "if(BOOPS_MAPS.pa_surveys){img.src=BOOPS_MAPS.pa_surveys;}",
            "return;}",
            "wrap.style.display='';",
            "var era=BOOPS_MAPS[scope];",
            "if(era&&era[bin]){img.src=era[bin];}",
            "}",
            "document.getElementById('boops-scope').addEventListener('change',refreshBoopsMap);",
            "document.getElementById('boops-bin').addEventListener('change',refreshBoopsMap);",
            "</script>"
          ),
          js_quote(boops_5_all),
          js_quote(boops_10_all),
          js_quote(boops_5_surv),
          js_quote(boops_10_surv),
          js_quote(boops_pa)
        )
      )
    },
    "</section>",

    "<section id='temperature'>",
    "<h2>4. Bottom temperature</h2>",
    "<p class='sub'>", stats$temp_note, "</p>",
    sprintf("<p><img class='chart' src='figures/%s' alt='temperature'/></p>", basename(paths$temp)),
    "</section>",

    "<section id='depth'>",
    "<h2>5. Distribution by depth</h2>",
    "<p class='sub'>Median depth per observation; 50 m bins.</p>",
    sprintf("<p><img class='chart' src='figures/%s' alt='depth'/></p>", basename(paths$depth)),
    "</section>",

    "</main></body></html>"
  )
  writeLines(html, file.path(OUT_DIR, "index.html"), useBytes = FALSE)
}

season_files <- c(
  file.path(OUT_DIR, "observations_by_month.csv"),
  file.path(OUT_DIR, "observations_summer_month_depth.csv"),
  file.path(FIG_DIR, "02_observations_by_month.png"),
  file.path(FIG_DIR, "06_summer_month_x_depth.png"),
  file.path(FIG_DIR, "06b_summer_depth_dodge.png")
)
for (f in season_files) {
  if (file.exists(f)) {
    file.remove(f)
  }
}

message("Loading records...")
rec <- load_observation_records()
obs <- collapse_observations(rec)
write_csv(obs, file.path(OUT_DIR, "observation_units.csv"))

rec_processed <- load_processed_dashboard_records()
obs_processed <- if (nrow(rec_processed) > 0) {
  collapse_observations(rec_processed)
} else {
  tibble()
}
if (nrow(obs_processed) > 0) {
  write_csv(obs_processed, file.path(OUT_DIR, "observation_units_processed.csv"))
}

by_db_survey <- obs %>%
  group_by(dataset) %>%
  summarise(
    observations = n(),
    sightings = sum(sightings),
    .groups = "drop"
  ) %>%
  mutate(layer = "Survey (Levant scope)")
by_db_processed <- if (nrow(obs_processed) > 0) {
  obs_processed %>%
    group_by(dataset) %>%
    summarise(
      observations = n(),
      sightings = sum(sightings),
      .groups = "drop"
    ) %>%
    mutate(layer = "Processed (`data/processed/`)")
} else {
  tibble(layer = character(), dataset = character(),
         observations = integer(), sightings = integer())
}
by_db <- bind_rows(by_db_survey, by_db_processed) %>%
  arrange(layer, desc(observations))
write_csv(by_db, file.path(OUT_DIR, "observations_by_database.csv"))

stats <- list(
  total = nrow(obs),
  sightings = sum(obs$sightings),
  datasets_survey = n_distinct(obs$dataset),
  processed_rows = nrow(rec_processed),
  temp_note = sprintf(
    "Bottom temperature (°C) available for %s of %s observations (MEDITS hauls only). Values ≤0 are missing (MEDITS uses −1 as placeholder).",
    sum(is.finite(obs$bottom_temp_c)),
    nrow(obs)
  )
)

paths <- list()

by_db$dataset <- factor(
  by_db$dataset,
  levels = by_db %>% arrange(layer, desc(observations)) %>% pull(dataset)
)
p_db <- ggplot(by_db, aes(dataset, observations, fill = layer)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = format(observations, big.mark = ",")), hjust = -0.12, size = 3.2) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(
    values = c("Survey (Levant scope)" = "#2ea043", "Processed (`data/processed/`)" = "#8250df"),
    name = NULL
  ) +
  labs(
    title = "Observation units by database",
    subtitle = "Survey layers (Levant bbox) vs processed pipeline exports",
    x = NULL, y = "Observation units"
  ) +
  theme_dash()
paths$by_db <- save_fig(p_db, "01_observations_by_database")

temp_df <- obs %>%
  filter(is.finite(bottom_temp_c), bottom_temp_c > 0, bottom_temp_c < 35)
if (nrow(temp_df) > 0) {
  p_temp <- ggplot(temp_df, aes(bottom_temp_c)) +
    geom_histogram(bins = 25, fill = "#d97706", colour = "white") +
    labs(
      title = "Bottom temperature (MEDITS observations)",
      x = "Temperature (°C)", y = "Observations"
    ) +
    theme_dash()
} else {
  p_temp <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No temperature data", size = 5) +
    theme_void()
}
paths$temp <- save_fig(p_temp, "04_temperature")

depth_ord <- obs %>%
  filter(!is.na(depth_band50)) %>%
  count(depth_band50, name = "observations") %>%
  mutate(depth_band50 = order_depth_band50_factor(depth_band50))

p_depth <- ggplot(depth_ord, aes(depth_band50, observations)) +
  geom_col(fill = "#2ea043") +
  coord_flip() +
  labs(title = "Observations by depth (50 m bins)", x = NULL, y = "Observations") +
  theme_dash()
paths$depth <- save_fig(p_depth, "05_depth", h = 7)

write_csv(
  obs %>% count(depth_band50, name = "observations"),
  file.path(OUT_DIR, "observations_by_depth.csv")
)

message("Species 2D maps for dashboard...")
species_catalog <- sync_species_map_figures()

message("Boops boops era maps...")
boops_df <- load_boops_occurrences()
write_csv(boops_df, file.path(OUT_DIR, "boops_boops_occurrences_med.csv"))
boops_era <- list(
  all = list(
    y10 = plot_boops_by_era(boops_df, 10L, "07_boops_boops_by_decade", surveys_only = FALSE),
    y5 = plot_boops_by_era(boops_df, 5L, "07b_boops_boops_by_5yr", surveys_only = FALSE)
  ),
  surveys = list(
    y10 = plot_boops_by_era(
      boops_df, 10L, "07d_boops_boops_by_decade_surveys", surveys_only = TRUE
    ),
    y5 = plot_boops_by_era(
      boops_df, 5L, "07e_boops_boops_by_5yr_surveys", surveys_only = TRUE
    )
  ),
  pa_map = {
    pa <- sync_boops_pa_map()
    if (is.null(pa)) NULL else basename(pa)
  }
)
if (!is.null(boops_era$all$y10$path)) {
  write_csv(
    boops_df %>%
      filter(!is.na(year)) %>%
      mutate(
        era_10 = paste0(10L * floor(year / 10L), "\u2013", 10L * floor(year / 10L) + 9L)
      ) %>%
      count(era_10, dataset, name = "records"),
    file.path(OUT_DIR, "boops_boops_records_by_decade.csv")
  )
}

write_html_dashboard(obs, by_db, paths, stats, species_catalog, boops_era)

message("Dashboard written: ", file.path(OUT_DIR, "index.html"))
message("Open in a browser: file://", normalizePath(file.path(OUT_DIR, "index.html"), winslash = "/"))
