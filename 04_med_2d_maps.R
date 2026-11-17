# =============================================================================
# 2D Mediterranean occurrence maps (lon/lat only — no depth colour)
# Fixed FOCUS_FIVE species · all datasets · full Med extent
# Outputs: primary_sdm/<species>/map_2d_MED_*.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

options(stringsAsFactors = FALSE)

source("gbif_paper.R")

FOCUS_FIVE <- c(
  "Pagellus erythrinus",
  "Boops boops",
  "Spicara maena",
  "Upeneus moluccensis",
  "Trachurus mediterraneus"
)

OUTPUT_DIR <- "primary_sdm"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Mediterranean Sea (incl. Alboran, Tyrrhenian, Adriatic, Levant)
BBOX_MED <- list(lon_min = -6.0, lon_max = 37.0, lat_min = 30.0, lat_max = 46.0)

DATASET_ORDER <- c(
  "MEDITS_raw",
  "Israel_trawl_MIT",
  "Edelist_OBIS",
  "Singer_Owens_iDigBio",
  GBIF_PAPER_DATASET
)

# Survey hauls for presence–absence (exclude GBIF + Singer/Owens iDigBio extract)
SURVEY_DATASETS <- c("MEDITS_raw", "Israel_trawl_MIT", "Edelist_OBIS")
PA_EXCLUDED <- c(GBIF_PAPER_DATASET, "Singer_Owens_iDigBio")

theme_dark_map <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.background   = element_rect(fill = "#0d1117", color = NA),
      panel.background  = element_rect(fill = "#161b22", color = NA),
      panel.border      = element_rect(fill = NA, color = "#30363d"),
      panel.grid.major  = element_line(color = "#21262d", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(color = "#f0f6fc", face = "bold"),
      plot.subtitle     = element_text(color = "#8b949e"),
      plot.caption      = element_text(color = "#8b949e", size = rel(0.85)),
      axis.text         = element_text(color = "#c9d1d9"),
      axis.title        = element_text(color = "#e6edf3"),
      legend.background = element_rect(fill = "#0d1117", color = NA),
      legend.key        = element_rect(fill = "#161b22", color = NA),
      legend.text       = element_text(color = "#c9d1d9"),
      legend.title      = element_text(color = "#e6edf3")
    )
}

MAP_LANDMARKS_MED <- tibble::tibble(
  label = c(
    "Gibraltar", "Barcelona", "Marseille", "Rome", "Athens", "Alexandria",
    "Tunis", "Haifa", "Tel Aviv", "Limassol", "Palermo", "Tripoli"
  ),
  decimalLongitude = c(
    -5.35, 2.17, 5.37, 12.50, 23.73, 29.92,
    10.18, 34.99, 34.78, 33.04, 13.36, 13.19
  ),
  decimalLatitude = c(
    36.14, 41.39, 43.30, 41.90, 37.98, 31.20,
    36.81, 32.79, 32.09, 34.68, 38.12, 32.89
  )
)

load_land_basemap <- function() {
  if (requireNamespace("maps", quietly = TRUE)) {
    return(list(kind = "df", data = ggplot2::map_data("world")))
  }
  if (requireNamespace("sf", quietly = TRUE) &&
      requireNamespace("rnaturalearth", quietly = TRUE)) {
    return(list(kind = "sf", data = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")))
  }
  message("Install package 'maps' for coastlines.")
  list(kind = "none", data = NULL)
}

MAP_BASE <- load_land_basemap()

sf_to_polygon_df <- function(land_sf) {
  if (is.null(land_sf) || nrow(land_sf) == 0) return(tibble())
  land_sf <- sf::st_cast(land_sf, "MULTIPOLYGON", warn = FALSE)
  purrr::map_dfr(seq_len(nrow(land_sf)), function(i) {
    coords <- sf::st_coordinates(land_sf$geometry[[i]])
    if (nrow(coords) == 0) return(tibble())
    tibble(
      long = coords[, "X"], lat = coords[, "Y"],
      group = i * 10000L + coords[, "L1"]
    )
  })
}

crop_land_for_limits <- function(lim, pad = 0.8) {
  xmin <- lim$xlim[1] - pad
  xmax <- lim$xlim[2] + pad
  ymin <- lim$ylim[1] - pad
  ymax <- lim$ylim[2] + pad
  if (MAP_BASE$kind == "sf") {
    bb <- sf::st_as_sfc(
      sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = 4326),
      crs = sf::st_crs(MAP_BASE$data)
    )
    land_df <- sf_to_polygon_df(suppressWarnings(sf::st_crop(MAP_BASE$data, bb)))
  } else if (MAP_BASE$kind == "df") {
    land_df <- MAP_BASE$data %>%
      filter(long >= xmin, long <= xmax, lat >= ymin, lat <= ymax)
  } else {
    return(NULL)
  }
  if (nrow(land_df) == 0) NULL else land_df
}

landmarks_for_limits <- function(lim, pad = 0.5) {
  MAP_LANDMARKS_MED %>%
    filter(
      decimalLongitude >= lim$xlim[1] - pad,
      decimalLongitude <= lim$xlim[2] + pad,
      decimalLatitude >= lim$ylim[1] - pad,
      decimalLatitude <= lim$ylim[2] + pad
    )
}

lim_from_map_box <- function(map_box) {
  list(xlim = c(map_box$lon_min, map_box$lon_max), ylim = c(map_box$lat_min, map_box$lat_max))
}

add_basemap_layers <- function(p, map_box = BBOX_MED) {
  lim <- lim_from_map_box(map_box)
  land <- crop_land_for_limits(lim)
  if (!is.null(land)) {
    p <- p + geom_polygon(
      data = land, aes(x = long, y = lat, group = group),
      fill = "#262b33", color = "#6e7681", linewidth = 0.2,
      inherit.aes = FALSE
    )
  }
  lm <- landmarks_for_limits(lim)
  if (nrow(lm) > 0) {
    p <- p +
      geom_point(data = lm, aes(x = decimalLongitude, y = decimalLatitude),
                 inherit.aes = FALSE, shape = 3, size = 1.4, color = "#8b949e") +
      geom_text(data = lm, aes(x = decimalLongitude, y = decimalLatitude, label = label),
                inherit.aes = FALSE, size = 2.3, color = "#e6edf3", vjust = -0.5,
                check_overlap = TRUE)
  }
  p
}

in_bbox <- function(df, box = BBOX_MED) {
  df %>%
    filter(
      decimalLongitude >= box$lon_min, decimalLongitude <= box$lon_max,
      decimalLatitude >= box$lat_min, decimalLatitude <= box$lat_max
    )
}

norm_species <- function(x) {
  x %>% as.character() %>% str_trim() %>%
    str_remove("\\s*\\([^)]*\\)$") %>% str_squish()
}

slug_species <- function(x) {
  x %>% str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>% str_sub(1, 60)
}

medits_dd <- function(x) {
  x <- as.numeric(x)
  deg <- floor(x / 100)
  deg + (x - deg * 100) / 60
}

thin_points <- function(df, max_n = 5000) {
  if (nrow(df) <= max_n) return(df)
  df %>% slice_sample(n = max_n)
}

pick_col <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit)) hit[[1]] else NA_character_
}

# ---- Loaders (Mediterranean clip) --------------------------------------------
load_medits_med <- function() {
  if (!file.exists("MEDITS/medits_raw_data.csv") || !file.exists("MEDITS/hauls_data.csv")) {
    return(tibble())
  }
  catch <- read_csv("MEDITS/medits_raw_data.csv", show_col_types = FALSE,
                    locale = locale(encoding = "Latin1"))
  hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE)
  taxa <- if (file.exists("MEDITS/taxa_medits.csv")) {
    read_csv("MEDITS/taxa_medits.csv", show_col_types = FALSE, locale = locale(encoding = "Latin1"))
  } else NULL

  catch <- catch %>%
    mutate(
      country = as.character(country), area = as.integer(area),
      vessel = as.character(vessel), year = as.integer(year),
      haul_number = as.numeric(haul_number)
    ) %>%
    filter(nbtot > 0) %>%
    mutate(medits_code = paste0(genus, species))

  if (!is.null(taxa)) {
    code_col <- names(taxa)[str_detect(names(taxa), "^MEDITS")]
    name_col <- names(taxa)[str_detect(names(taxa), "Scientific Name")]
    taxa <- taxa %>%
      mutate(
        medits_code = .data[[code_col[1]]],
        scientificName = .data[[name_col[1]]]
      ) %>%
      select(medits_code, scientificName, CATFAU)
    catch <- catch %>% left_join(taxa, by = "medits_code") %>%
      filter(is.na(CATFAU) | CATFAU %in% c("Ao", "Ae", "C"))
  }

  hauls <- hauls %>%
    mutate(
      country = as.character(country), area = as.integer(area),
      vessel = as.character(vessel), year = as.integer(year),
      haul_number = as.numeric(haul_number)
    ) %>%
    filter(validity == "V") %>%
    mutate(
      decimalLatitude = medits_dd(shooting_latitude),
      decimalLongitude = medits_dd(shooting_longitude),
      depth = (as.numeric(shooting_depth) + as.numeric(hauling_depth)) / 2
    )

  catch %>%
    inner_join(
      hauls %>% select(country, area, vessel, year, haul_number,
                       decimalLatitude, decimalLongitude, depth),
      by = c("country", "area", "vessel", "year", "haul_number")
    ) %>%
    transmute(
      scientificName = coalesce(scientificName, medits_code),
      decimalLatitude, decimalLongitude, depth,
      dataset = "MEDITS_raw", year = as.integer(year),
      obs_key = paste(country, area, vessel, year, haul_number, sep = "|")
    ) %>%
    in_bbox(BBOX_MED) %>%
    distinct()
}

load_singer_med <- function(focus_keys) {
  path <- "data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv"
  if (!file.exists(path)) return(tibble())
  cache <- "data/processed/singer_med_focus5.csv"
  if (file.exists(cache) && file.info(cache)$mtime > file.info(path)$mtime) {
    return(read_csv(cache, show_col_types = FALSE))
  }
  hdr <- names(read_csv(path, n_max = 0, show_col_types = FALSE))
  sp_c <- pick_col(hdr, c("WoRMS_scientificName", "SOcurated_scientificName", "dwc_scientificName"))
  lat_c <- pick_col(hdr, c("iDigBio_dwc_decimalLatitude", "dwc_decimalLatitude"))
  lon_c <- pick_col(hdr, c("iDigBio_dwc_decimalLongitude", "dwc_decimalLongitude"))
  if (any(is.na(c(sp_c, lat_c, lon_c)))) return(tibble())

  raw <- read_csv(path, show_col_types = FALSE,
                  col_select = all_of(c(sp_c, lat_c, lon_c)), progress = FALSE) %>%
    transmute(
      scientificName = norm_species(.data[[sp_c]]),
      decimalLatitude = as.numeric(.data[[lat_c]]),
      decimalLongitude = as.numeric(.data[[lon_c]]),
      depth = NA_real_,
      dataset = "Singer_Owens_iDigBio",
      year = NA_integer_
    ) %>%
    filter(
      scientificName %in% focus_keys,
      !scientificName %in% c("Unaccepted", "Equivocal"),
      is.finite(decimalLatitude), is.finite(decimalLongitude)
    ) %>%
    in_bbox(BBOX_MED)

  dir.create("data/processed", showWarnings = FALSE)
  write_csv(raw, cache)
  raw
}

load_all_med <- function(focus_keys) {
  out <- list()

  med <- load_medits_med()
  if (nrow(med) > 0) {
    out[["MEDITS_raw"]] <- med
    message("MEDITS (Med bbox): ", nrow(med), " records")
  }

  if (file.exists("Combined_Menachem_Itai_corrected_May_2026_DBG.csv")) {
    mit <- read_csv("Combined_Menachem_Itai_corrected_May_2026_DBG.csv", show_col_types = FALSE) %>%
      transmute(
        scientificName = as.character(ID_SPECIES),
        decimalLatitude = as.numeric(LATITUDE),
        decimalLongitude = as.numeric(LONGITUDE),
        depth = as.numeric(DEPTH),
        dataset = "Israel_trawl_MIT",
        year = as.integer(YEAR),
        obs_key = as.character(SAMPLE_DESC)
      ) %>%
      in_bbox(BBOX_MED)
    out[["Israel_trawl_MIT"]] <- mit
  }

  if (file.exists("data/raw/edelist_occurrence.txt")) {
    ed <- read_tsv("data/raw/edelist_occurrence.txt", show_col_types = FALSE, guess_max = 50000)
    id_col <- if ("id" %in% names(ed)) "id" else "occurrenceID"
    ed <- ed %>%
      filter(class == "Actinopterygii" | is.na(class)) %>%
      mutate(
        obs_key = str_replace(.data[[id_col]], "_Fish_.*$", ""),
        parsed = str_match(.data[[id_col]], "^([0-9.]+)_([0-9.]+)_([0-9.]+)")
      ) %>%
      transmute(
        scientificName = as.character(scientificName),
        decimalLatitude = as.numeric(parsed[, 2]),
        decimalLongitude = as.numeric(parsed[, 3]),
        depth = as.numeric(parsed[, 4]),
        dataset = "Edelist_OBIS",
        year = suppressWarnings(as.integer(str_match(.data[[id_col]], "_([0-9]{4})_")[, 2])),
        obs_key
      ) %>%
      in_bbox(BBOX_MED)
    out[["Edelist_OBIS"]] <- ed
  }

  singer <- load_singer_med(focus_keys)
  if (nrow(singer) > 0) out[["Singer_Owens_iDigBio"]] <- singer

  gb <- load_gbif_paper_cleaned()
  if (nrow(gb) > 0) {
    gb_std <- gb %>%
      transmute(
        scientificName = norm_species(scientificName),
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        depth = as.numeric(depth),
        dataset = GBIF_PAPER_DATASET,
        year = suppressWarnings(as.integer(substr(eventDate, 1, 4))),
        obs_key = NA_character_
      ) %>%
      in_bbox(BBOX_MED)
    out[[GBIF_PAPER_DATASET]] <- gb_std
    message(
      GBIF_PAPER_DATASET, ": ", nrow(gb_std),
      " paper-cleaned rows (Levant bbox; clipped to Med map extent)"
    )
  }

  bind_rows(out)
}

plot_2d_map <- function(df, title, subtitle, outfile, map_box = BBOX_MED,
                        color_by_dataset = TRUE) {
  df <- df %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
  if (nrow(df) == 0) {
    message("Skip empty 2D map: ", outfile)
    return(invisible(NULL))
  }
  lim <- lim_from_map_box(map_box)
  p <- add_basemap_layers(ggplot(), map_box)
  if (color_by_dataset) {
    p <- p +
      geom_point(
        data = df,
        aes(decimalLongitude, decimalLatitude, color = dataset),
        alpha = 0.65, size = 1.2
      ) +
      scale_color_brewer(palette = "Set2", name = "Dataset")
  } else {
    p <- p +
      geom_point(
        data = df,
        aes(decimalLongitude, decimalLatitude),
        color = "#58a6ff", alpha = 0.7, size = 1.3
      )
  }
  p <- p +
    coord_fixed(1.4, xlim = lim$xlim, ylim = lim$ylim) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude (°E)", y = "Latitude (°N)",
      caption = "2D presence map · Mediterranean extent · depth not shown"
    ) +
    theme_dark_map()
  ggsave(outfile, p, width = 11, height = 7.5, dpi = 150, bg = "#0d1117")
  message("Wrote ", outfile)
}

other_species_stations <- function(survey_df, focus_sp_key) {
  if (nrow(survey_df) == 0) return(tibble())

  station_flags <- survey_df %>%
    group_by(dataset, obs_key) %>%
    summarise(has_focus = any(sp_key == focus_sp_key), .groups = "drop")

  other_keys <- station_flags %>% filter(!has_focus)
  if (nrow(other_keys) == 0) return(tibble())

  survey_df %>%
    inner_join(other_keys, by = c("dataset", "obs_key")) %>%
    group_by(dataset, obs_key) %>%
    summarise(
      decimalLongitude = median(decimalLongitude, na.rm = TRUE),
      decimalLatitude = median(decimalLatitude, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
}

plot_2d_pa_map <- function(presence_df, absence_df, title, subtitle, outfile,
                           map_box = BBOX_MED) {
  presence_df <- presence_df %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
  if (nrow(presence_df) == 0) {
    message("Skip empty PA map: ", outfile)
    return(invisible(NULL))
  }

  lim <- lim_from_map_box(map_box)
  p <- add_basemap_layers(ggplot(), map_box)

  if (nrow(absence_df) > 0) {
    p <- p +
      geom_point(
        data = absence_df,
        aes(decimalLongitude, decimalLatitude),
        color = "#6e7681", shape = 4, size = 1.4, alpha = 0.45, stroke = 0.35
      )
  }

  p <- p +
    geom_point(
      data = presence_df,
      aes(decimalLongitude, decimalLatitude),
      color = "#3fb950", alpha = 0.75, size = 1.5
    ) +
    coord_fixed(1.4, xlim = lim$xlim, ylim = lim$ylim) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude (°E)", y = "Latitude (°N)",
      caption = "Green = presence · grey cross = survey haul with other fish, species absent (2D)"
    ) +
    theme_dark_map()

  ggsave(outfile, p, width = 11, height = 7.5, dpi = 150, bg = "#0d1117")
  message("Wrote ", outfile)
}

# ---- Main --------------------------------------------------------------------
message("Loading Mediterranean occurrences...")
focus_keys <- norm_species(FOCUS_FIVE)
all_med <- load_all_med(focus_keys) %>%
  mutate(sp_key = norm_species(scientificName))

datasets_present <- intersect(DATASET_ORDER, unique(all_med$dataset))
if (length(datasets_present) == 0) {
  stop("No datasets loaded for Mediterranean maps.")
}

survey_occ <- all_med %>%
  filter(
    dataset %in% SURVEY_DATASETS,
    !is.na(obs_key), obs_key != "",
    is.finite(decimalLongitude), is.finite(decimalLatitude)
  )

for (sp in FOCUS_FIVE) {
  sp_key <- norm_species(sp)
  sp_slug <- slug_species(sp)
  sp_dir <- file.path(OUTPUT_DIR, sp_slug)
  dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)

  sp_df <- all_med %>% filter(sp_key == .env$sp_key)
  if (nrow(sp_df) == 0) {
    message("No Mediterranean records for ", sp)
    next
  }

  for (ds in datasets_present) {
    ds_df <- sp_df %>% filter(dataset == ds) %>% thin_points(6000)
    if (nrow(ds_df) == 0) next
    plot_2d_map(
      ds_df,
      title = paste(sp, "— Mediterranean"),
      subtitle = paste0("Dataset: ", ds, " · 2D occurrence (no depth axis)"),
      outfile = file.path(sp_dir, paste0("map_2d_MED_", ds, ".png")),
      color_by_dataset = FALSE
    )
  }

  plot_2d_map(
    sp_df %>% thin_points(8000),
    title = paste(sp, "— Mediterranean (all datasets)"),
    subtitle = "Joint 2D map · colour = data source",
    outfile = file.path(sp_dir, "map_2d_MED_all_datasets.png"),
    color_by_dataset = TRUE
  )

  ds_with_data <- sp_df %>% count(dataset) %>% filter(n > 0) %>% pull(dataset)
  ds_with_data <- intersect(DATASET_ORDER, ds_with_data)
  if (length(ds_with_data) >= 2) {
    panel_df <- sp_df %>%
      filter(dataset %in% ds_with_data) %>%
      thin_points(4000)
    lim <- lim_from_map_box(BBOX_MED)
    land <- crop_land_for_limits(lim)
    lm <- landmarks_for_limits(lim)
    p <- ggplot()
    if (!is.null(land) && nrow(land) > 0) {
      p <- p + geom_polygon(
        data = land, aes(long, lat, group = group),
        fill = "#262b33", color = "#6e7681", linewidth = 0.2, inherit.aes = FALSE
      )
    }
    if (nrow(lm) > 0) {
      p <- p +
        geom_point(data = lm, aes(decimalLongitude, decimalLatitude),
                   inherit.aes = FALSE, shape = 3, size = 1.2, color = "#8b949e") +
        geom_text(data = lm, aes(decimalLongitude, decimalLatitude, label = label),
                  inherit.aes = FALSE, size = 2, color = "#e6edf3", vjust = -0.5,
                  check_overlap = TRUE)
    }
    p <- p +
      geom_point(
        data = panel_df,
        aes(decimalLongitude, decimalLatitude),
        color = "#58a6ff", alpha = 0.7, size = 0.9
      ) +
      facet_wrap(~dataset, ncol = 2, scales = "free") +
      coord_cartesian(expand = TRUE) +
      labs(
        title = paste(sp, "— Mediterranean by dataset"),
        subtitle = "2D panels · no depth",
        x = "Longitude (°E)", y = "Latitude (°N)"
      ) +
      theme_dark_map(10)
    ggsave(
      file.path(sp_dir, "map_2d_MED_panel.png"),
      p, width = 12, height = 3 + 2.2 * ceiling(length(ds_with_data) / 2),
      dpi = 150, bg = "#0d1117"
    )
    message("Wrote ", file.path(sp_dir, "map_2d_MED_panel.png"))
  }

  # Presence–absence (surveys only; no GBIF, no Singer/Owens iDigBio)
  pa_datasets <- intersect(SURVEY_DATASETS, unique(survey_occ$dataset))

  for (ds in pa_datasets) {
    survey_ds <- survey_occ %>% filter(dataset == ds)
    presence_ds <- sp_df %>%
      filter(dataset == ds) %>%
      thin_points(6000)
    if (nrow(presence_ds) == 0) next
    absence_ds <- other_species_stations(survey_ds, sp_key) %>% thin_points(5000)
    plot_2d_pa_map(
      presence_df = presence_ds,
      absence_df = absence_ds,
      title = paste(sp, "— presence / absence"),
      subtitle = paste0(
        ds, " · Mediterranean 2D · ",
        nrow(absence_ds), " hauls with other fish, ", sp, " absent"
      ),
      outfile = file.path(sp_dir, paste0("map_2d_MED_PA_", ds, ".png"))
    )
  }

  presence_joint <- sp_df %>%
    filter(dataset %in% pa_datasets) %>%
    thin_points(8000)
  if (nrow(presence_joint) > 0) {
    absence_joint <- other_species_stations(survey_occ, sp_key) %>% thin_points(6000)
    plot_2d_pa_map(
      presence_df = presence_joint,
      absence_df = absence_joint,
      title = paste(sp, "— joint survey presence / absence"),
      subtitle = paste0(
        "MIT + Edelist + MEDITS · ",
        nrow(absence_joint), " pseudo-absence hauls (other species only)"
      ),
      outfile = file.path(sp_dir, "map_2d_MED_PA_joint_surveys.png")
    )
  }
}

message("\nDone. 2D Mediterranean maps in ", OUTPUT_DIR, "/<species>/map_2d_MED_*.png")
message("Presence–absence: map_2d_MED_PA_*.png (surveys only; excludes GBIF & Singer_Owens_iDigBio)")
