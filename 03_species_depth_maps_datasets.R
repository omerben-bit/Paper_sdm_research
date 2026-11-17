# =============================================================================
# Depth-coloured lon/lat maps for 5 Israeli fish species
# One map per dataset + one joint map per species (dark mode)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

options(stringsAsFactors = FALSE)

source("gbif_paper.R")

# ---- Map extents ---------------------------------------------------------------
BBOX_ISRAEL <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 34.5)
# MEDITS GSA 22/25 (Cyprus / Levantine slope) sits north of Israel — separate extent
BBOX_MEDITS <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 33.2, lat_max = 35.3)
# Joint maps (Israel surveys + Levantine MEDITS)
BBOX_JOINT <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 35.3)

FOCUS_FIVE <- c(
  "Pagellus erythrinus",
  "Boops boops",
  "Spicara maena",
  "Upeneus moluccensis",
  "Trachurus mediterraneus"
)

OUTPUT_DIR <- "primary_sdm"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_dark_map <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.background   = element_rect(fill = "#0d1117", color = NA),
      panel.background  = element_rect(fill = "#161b22", color = NA),
      panel.border      = element_rect(fill = NA, color = "#30363d"),
      panel.grid.major  = element_line(color = "#21262d", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#21262d", color = "#30363d"),
      strip.text        = element_text(color = "#e6edf3", face = "bold"),
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

# ---- Basemap: land + coastal landmarks ---------------------------------------
MAP_LANDMARKS <- tibble::tibble(
  label = c(
    "Tel Aviv", "Haifa", "Ashdod", "Eilat", "Jerusalem",
    "Gaza", "Beirut", "Limassol", "Larnaca", "Nicosia"
  ),
  decimalLongitude = c(
    34.78, 34.99, 34.65, 34.95, 35.21,
    34.47, 35.50, 33.02, 33.63, 33.38
  ),
  decimalLatitude = c(
    32.09, 32.79, 31.80, 29.56, 31.77,
    31.50, 33.89, 34.68, 34.92, 35.19
  )
)

load_land_basemap <- function() {
  if (requireNamespace("maps", quietly = TRUE)) {
    return(list(kind = "df", data = ggplot2::map_data("world")))
  }
  if (requireNamespace("sf", quietly = TRUE) &&
      requireNamespace("rnaturalearth", quietly = TRUE)) {
    land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    return(list(kind = "sf", data = land))
  }
  message("Tip: install 'maps' or ('sf' + 'rnaturalearth') for coastlines on maps.")
  list(kind = "none", data = NULL)
}

MAP_BASE <- load_land_basemap()

sf_to_polygon_df <- function(land_sf) {
  if (is.null(land_sf) || nrow(land_sf) == 0) {
    return(tibble())
  }
  land_sf <- sf::st_cast(land_sf, "MULTIPOLYGON", warn = FALSE)
  idx <- seq_len(nrow(land_sf))
  purrr::map_dfr(idx, function(i) {
    geom_i <- land_sf$geometry[[i]]
    if (is.null(geom_i) || length(geom_i) == 0) {
      return(tibble())
    }
    coords <- sf::st_coordinates(geom_i)
    if (nrow(coords) == 0) {
      return(tibble())
    }
    tibble(
      long = coords[, "X"],
      lat = coords[, "Y"],
      group = i * 10000L + coords[, "L1"],
      order = seq_len(nrow(coords))
    )
  })
}

crop_land_for_limits <- function(lim, pad = 0.4) {
  xmin <- lim$xlim[1] - pad
  xmax <- lim$xlim[2] + pad
  ymin <- lim$ylim[1] - pad
  ymax <- lim$ylim[2] + pad

  if (MAP_BASE$kind == "sf") {
    bb <- sf::st_as_sfc(
      sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = 4326),
      crs = sf::st_crs(MAP_BASE$data)
    )
    cropped <- suppressWarnings(sf::st_crop(MAP_BASE$data, bb))
    land_df <- sf_to_polygon_df(cropped)
  } else if (MAP_BASE$kind == "df") {
    land_df <- MAP_BASE$data %>%
      filter(long >= xmin, long <= xmax, lat >= ymin, lat <= ymax)
  } else {
    return(NULL)
  }

  if (nrow(land_df) == 0) {
    return(NULL)
  }
  land_df
}

landmarks_for_limits <- function(lim, pad = 0.35) {
  MAP_LANDMARKS %>%
    filter(
      decimalLongitude >= lim$xlim[1] - pad,
      decimalLongitude <= lim$xlim[2] + pad,
      decimalLatitude >= lim$ylim[1] - pad,
      decimalLatitude <= lim$ylim[2] + pad
    )
}

lim_from_map_box <- function(map_box) {
  list(
    xlim = c(map_box$lon_min, map_box$lon_max),
    ylim = c(map_box$lat_min, map_box$lat_max)
  )
}

map_box_for_dataset <- function(ds) {
  if (ds == "MEDITS_raw") BBOX_MEDITS else BBOX_ISRAEL
}

build_faceted_land <- function(ds_list) {
  purrr::map_dfr(ds_list, function(ds) {
    lim <- lim_from_map_box(map_box_for_dataset(ds))
    land <- crop_land_for_limits(lim)
    if (is.null(land) || nrow(land) == 0) {
      return(tibble())
    }
    land %>%
      mutate(dataset = ds, poly_group = paste(ds, group, sep = ":"))
  })
}

build_faceted_landmarks <- function(ds_list) {
  purrr::map_dfr(ds_list, function(ds) {
    lim <- lim_from_map_box(map_box_for_dataset(ds))
    landmarks_for_limits(lim) %>% mutate(dataset = ds)
  })
}

add_basemap_layers <- function(p, lim, map_box = NULL) {
  if (!is.null(map_box)) {
    lim <- lim_from_map_box(map_box)
  }
  land <- crop_land_for_limits(lim)
  if (!is.null(land)) {
    p <- p + ggplot2::geom_polygon(
      data = land,
      aes(x = long, y = lat, group = group),
      fill = "#262b33",
      color = "#6e7681",
      linewidth = 0.25,
      inherit.aes = FALSE
    )
  }

  lm <- landmarks_for_limits(lim)
  if (nrow(lm) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = lm,
        aes(x = decimalLongitude, y = decimalLatitude),
        inherit.aes = FALSE,
        shape = 3,
        size = 1.6,
        color = "#8b949e"
      ) +
      ggplot2::geom_text(
        data = lm,
        aes(x = decimalLongitude, y = decimalLatitude, label = label),
        inherit.aes = FALSE,
        size = 2.6,
        color = "#e6edf3",
        vjust = -0.6,
        check_overlap = TRUE
      )
  }
  p
}

in_bbox <- function(df, box = BBOX_ISRAEL) {
  df %>%
    filter(
      decimalLongitude >= box$lon_min, decimalLongitude <= box$lon_max,
      decimalLatitude >= box$lat_min, decimalLatitude <= box$lat_max
    )
}

normalize_haul_keys <- function(df) {
  df %>%
    mutate(
      country = as.character(country),
      area = as.integer(area),
      vessel = as.character(vessel),
      year = as.integer(year),
      haul_number = as.numeric(haul_number)
    )
}

# Primary fish table: MEDITS/medits_raw_data.csv + haul coordinates
load_medits_raw <- function() {
  path_catch <- "MEDITS/medits_raw_data.csv"
  path_hauls <- "MEDITS/hauls_data.csv"
  if (!file.exists(path_catch) || !file.exists(path_hauls)) {
    return(tibble())
  }

  catch <- read_csv(path_catch, show_col_types = FALSE, locale = locale(encoding = "Latin1"))
  hauls <- read_csv(path_hauls, show_col_types = FALSE)
  taxa <- if (file.exists("MEDITS/taxa_medits.csv")) {
    read_csv("MEDITS/taxa_medits.csv", show_col_types = FALSE, locale = locale(encoding = "Latin1"))
  } else {
    NULL
  }

  catch <- normalize_haul_keys(catch) %>% filter(nbtot > 0)
  hauls <- normalize_haul_keys(hauls) %>%
    filter(validity == "V") %>%
    mutate(
      decimalLatitude = medits_dd(shooting_latitude),
      decimalLongitude = medits_dd(shooting_longitude),
      depth = (as.numeric(shooting_depth) + as.numeric(hauling_depth)) / 2
    )

  catch <- catch %>% mutate(medits_code = paste0(genus, species))

  if (!is.null(taxa)) {
    code_col <- names(taxa)[str_detect(names(taxa), "^MEDITS")]
    name_col <- names(taxa)[str_detect(names(taxa), "Scientific Name")]
    taxa <- taxa %>%
      mutate(
        medits_code = .data[[code_col[1]]],
        scientificName = .data[[name_col[1]]]
      ) %>%
      select(medits_code, scientificName, CATFAU)
    catch <- catch %>% left_join(taxa, by = "medits_code")
    catch <- catch %>%
      filter(is.na(CATFAU) | CATFAU %in% c("Ao", "Ae", "C"))
  }

  catch %>%
    inner_join(
      hauls %>%
        select(country, area, vessel, year, haul_number,
               decimalLatitude, decimalLongitude, depth),
      by = c("country", "area", "vessel", "year", "haul_number")
    ) %>%
    transmute(
      scientificName = coalesce(scientificName, medits_code),
      decimalLatitude,
      decimalLongitude,
      depth,
      dataset = "MEDITS_raw",
      year = as.integer(year),
      gsa = area,
      obs_key = paste(country, area, vessel, year, haul_number, sep = "|")
    ) %>%
    filter(gsa %in% c(22L, 25L)) %>%
    in_bbox(BBOX_MEDITS) %>%
    distinct()
}

slug_species <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_sub(1, 60)
}

# Match "Boops boops (Linnaeus, 1758)" with "Boops boops" across GBIF vs surveys
norm_species <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_remove("\\s*\\([^)]*\\)$") %>%
    str_squish()
}

pick_col <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit)) hit[[1]] else NA_character_
}

# Singer & Owens (2026) Scientific Data — depth-curated iDigBio specimens (Zenodo)
load_singer_owens <- function() {
  paths <- c(
    "data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv",
    "data/raw/processedFishDataWithDatesAndNames.csv"
  )
  path <- paths[file.exists(paths)][1]
  if (is.na(path)) {
    message(
      "Singer & Owens dataset not found. Download iDigBioProcessedFishOccurrences_SingerAndOwens.csv ",
      "from https://doi.org/10.5281/zenodo.22124920 into data/raw/"
    )
    return(tibble())
  }

  cache <- "data/processed/singer_owens_israel_box.csv"
  if (file.exists(cache) && file.info(cache)$mtime > file.info(path)$mtime) {
    return(read_csv(cache, show_col_types = FALSE))
  }

  hdr <- names(read_csv(path, n_max = 0, show_col_types = FALSE))
  sp_c <- pick_col(hdr, c(
    "WoRMS_scientificName", "SOcurated_scientificName", "dwc_scientificName", "scientificName"
  ))
  lat_c <- pick_col(hdr, c("iDigBio_dwc_decimalLatitude", "dwc_decimalLatitude", "decimalLatitude"))
  lon_c <- pick_col(hdr, c("iDigBio_dwc_decimalLongitude", "dwc_decimalLongitude", "decimalLongitude"))
  dep_c <- pick_col(hdr, c("SOcurated_depth", "depth"))
  yr_c <- pick_col(hdr, c("iDigBio_dwc_year", "dwc_year", "year"))

  if (any(is.na(c(sp_c, lat_c, lon_c, dep_c)))) {
    message("Singer & Owens CSV missing expected columns; check Zenodo file version.")
    return(tibble())
  }

  use_cols <- unique(c(sp_c, lat_c, lon_c, dep_c, yr_c))
  raw <- read_csv(path, show_col_types = FALSE, col_select = all_of(use_cols), progress = FALSE)

  out <- raw %>%
    transmute(
      scientificName = norm_species(as.character(.data[[sp_c]])),
      decimalLatitude = as.numeric(.data[[lat_c]]),
      decimalLongitude = as.numeric(.data[[lon_c]]),
      depth = as.numeric(.data[[dep_c]]),
      dataset = "Singer_Owens_iDigBio",
      year = if (!is.na(yr_c)) suppressWarnings(as.integer(.data[[yr_c]])) else NA_integer_,
      obs_key = NA_character_
    ) %>%
    filter(
      !is.na(scientificName), scientificName != "",
      !scientificName %in% c("Unaccepted", "Equivocal"),
      is.finite(decimalLatitude), is.finite(decimalLongitude), is.finite(depth)
    ) %>%
    in_bbox(BBOX_ISRAEL) %>%
    distinct()

  dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
  write_csv(out, cache)
  message("Singer_Owens_iDigBio records (Israel box, depth-resolved): ", nrow(out))
  out
}

# MEDITS coordinates: DDMM.mm (degrees + decimal minutes in one number)
medits_dd <- function(x) {
  x <- as.numeric(x)
  deg <- floor(x / 100)
  minutes <- x - deg * 100
  deg + minutes / 60
}

# ---- Load datasets -----------------------------------------------------------
load_all_datasets <- function() {
  out <- list()

  # 1) Moria folder (user uploads)
  moria_files <- list.files("moria", pattern = "\\.(csv|tsv|txt)$", full.names = TRUE, ignore.case = TRUE)
  moria_files <- moria_files[!str_detect(tolower(moria_files), "readme")]
  for (f in moria_files) {
    nm <- paste0("Moria_", tools::file_path_sans_ext(basename(f)))
    raw <- if (str_detect(f, "\\.tsv$")) read_tsv(f, show_col_types = FALSE) else read_csv(f, show_col_types = FALSE)
    out[[nm]] <- standardize_table(raw, nm)
  }

  # 2) Israel trawl (Menachem / Itai combined)
  mit_path <- "Combined_Menachem_Itai_corrected_May_2026_DBG.csv"
  if (file.exists(mit_path)) {
    mit <- read_csv(mit_path, show_col_types = FALSE)
    mit_std <- mit %>%
      transmute(
        scientificName = as.character(ID_SPECIES),
        decimalLatitude = as.numeric(LATITUDE),
        decimalLongitude = as.numeric(LONGITUDE),
        depth = as.numeric(DEPTH),
        dataset = "Israel_trawl_MIT",
        year = as.integer(YEAR),
        obs_key = as.character(SAMPLE_DESC)
      ) %>%
      filter(!is.na(scientificName), !is.na(decimalLatitude), !is.na(decimalLongitude))
    out[["Israel_trawl_MIT"]] <- mit_std
  }

  # 3) MEDITS — medits_raw_data.csv (+ hauls_data.csv for position/depth)
  medits_df <- load_medits_raw()
  if (nrow(medits_df) > 0) {
    out[["MEDITS_raw"]] <- medits_df
    message("MEDITS_raw records (GSA 22/25, Levantine box): ", nrow(medits_df))
  }

  # 4) Edelist Israeli Med trawls (coords + depth encoded in id field)
  if (file.exists("data/raw/edelist_occurrence.txt")) {
    ed <- read_tsv("data/raw/edelist_occurrence.txt", show_col_types = FALSE, guess_max = 50000)
    id_col <- if ("id" %in% names(ed)) "id" else "occurrenceID"
    parsed <- str_match(ed[[id_col]], "^([0-9.]+)_([0-9.]+)_([0-9.]+)")
    ed <- ed %>%
      mutate(
        decimalLatitude = as.numeric(parsed[, 2]),
        decimalLongitude = as.numeric(parsed[, 3]),
        depth = as.numeric(parsed[, 4]),
        scientificName = as.character(scientificName),
        dataset = "Edelist_OBIS",
        year = suppressWarnings(as.integer(str_match(ed[[id_col]], "_([0-9]{4})_")[, 2])),
        obs_key = str_replace(.data[[id_col]], "_Fish_.*$", "")
      ) %>%
      filter(
        !is.na(decimalLatitude), !is.na(decimalLongitude),
        !is.na(scientificName), scientificName != "",
        class == "Actinopterygii" | is.na(class)
      ) %>%
      select(scientificName, decimalLatitude, decimalLongitude, depth, dataset, year, obs_key)
    out[["Edelist_OBIS"]] <- ed
  }

  # 5) Singer & Owens (2026) — depth-curated iDigBio / museum specimens (Zenodo)
  singer_df <- load_singer_owens()
  if (nrow(singer_df) > 0) {
    out[["Singer_Owens_iDigBio"]] <- singer_df
  }

  # 6) Paper-cleaned GBIF (01 pipeline; presence-only, sparse depth)
  gb <- load_gbif_paper_cleaned()
  if (nrow(gb) > 0) {
    gb <- gb %>% mutate(dataset = GBIF_PAPER_DATASET)
    out[[GBIF_PAPER_DATASET]] <- gb %>%
      transmute(
        scientificName = norm_species(scientificName),
        decimalLatitude, decimalLongitude, depth,
        dataset, year = suppressWarnings(as.integer(substr(eventDate, 1, 4))),
        obs_key = NA_character_
      )
  }

  # Israel-coast datasets → Israel bbox; MEDITS + Singer pre-filtered
  skip_bbox <- c("MEDITS_raw", "Singer_Owens_iDigBio")
  israel_parts <- bind_rows(out[!names(out) %in% skip_bbox]) %>% in_bbox(BBOX_ISRAEL)
  singer_part <- out[["Singer_Owens_iDigBio"]]
  if (is.null(singer_part)) singer_part <- tibble()
  medits_part <- out[["MEDITS_raw"]]
  if (is.null(medits_part)) medits_part <- tibble()
  bind_rows(israel_parts, singer_part, medits_part)
}

standardize_table <- function(df, dataset_name) {
  nm <- names(df)
  pick <- function(cands) {
    hit <- cands[cands %in% nm]
    if (length(hit)) hit[[1]] else NA_character_
  }
  lon <- pick(c("decimalLongitude", "LONGITUDE", "longitude", "lon", "Lon"))
  lat <- pick(c("decimalLatitude", "LATITUDE", "latitude", "lat", "Lat"))
  sp  <- pick(c("scientificName", "ID_SPECIES", "species", "Species", "acceptedNameUsage"))
  dep <- pick(c("depth", "DEPTH", "depthInMeters", "minimumDepthInMeters"))
  yr  <- pick(c("year", "YEAR", "eventDate"))

  if (is.na(lon) || is.na(lat) || is.na(sp)) {
    return(tibble())
  }

  obs <- pick(c("obs_key", "SAMPLE_DESC", "sample_desc"))
  tibble(
    scientificName = as.character(df[[sp]]),
    decimalLatitude = as.numeric(df[[lat]]),
    decimalLongitude = as.numeric(df[[lon]]),
    depth = if (!is.na(dep)) as.numeric(df[[dep]]) else NA_real_,
    dataset = dataset_name,
    obs_key = if (!is.na(obs)) as.character(df[[obs]]) else NA_character_,
    year = if (!is.na(yr)) {
      if (yr == "eventDate") suppressWarnings(as.integer(substr(as.character(df[[yr]]), 1, 4)))
      else suppressWarnings(as.integer(df[[yr]]))
    } else {
      NA_integer_
    }
  ) %>%
    filter(!is.na(scientificName), !is.na(decimalLatitude), !is.na(decimalLongitude))
}

# ---- Fixed focus: 5 common Levant demersal species -----------------------------
message("Loading datasets...")
all_occ <- load_all_datasets() %>%
  mutate(sp_key = norm_species(scientificName))
write_csv(all_occ, file.path(OUTPUT_DIR, "all_occurrences_israel_box.csv"))

if (nrow(all_occ) == 0) stop("No occurrence records in Israel box. Check data paths / moria/ folder.")

SURVEY_DATASETS <- c("Israel_trawl_MIT", "Edelist_OBIS", "MEDITS_raw")
COLLECTION_DATASET <- "Singer_Owens_iDigBio" # Singer & Owens 2026 Scientific Data (Zenodo)

species_by_key <- all_occ %>%
  group_by(sp_key) %>%
  summarise(
    n_survey_datasets = n_distinct(
      dataset[dataset %in% SURVEY_DATASETS & is.finite(depth)]
    ),
    n_survey_depth = sum(dataset %in% SURVEY_DATASETS & is.finite(depth)),
    n_total = n(),
    .groups = "drop"
  )

pick5 <- tibble(
  scientificName = FOCUS_FIVE,
  sp_key = norm_species(FOCUS_FIVE)
) %>%
  left_join(species_by_key, by = "sp_key")

missing_sp <- pick5 %>% filter(is.na(n_total) | n_survey_depth == 0)
if (nrow(missing_sp) > 0) {
  message("Warning: no survey depth records for: ", paste(missing_sp$scientificName, collapse = ", "))
}

message("Focus species (fixed list of 5):")
print(pick5 %>% select(scientificName, n_survey_datasets, n_survey_depth, n_total))
write_csv(pick5, file.path(OUTPUT_DIR, "focus_5_species.csv"))

# Thin duplicate points for plotting (same lon/lat/depth/dataset)
thin_points <- function(df, max_n = 4000) {
  if (nrow(df) <= max_n) return(df)
  df %>% slice_sample(n = max_n)
}

plot_limits <- function(df, pad = 0.08) {
  x <- df$decimalLongitude
  y <- df$decimalLatitude
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  dx <- diff(xr) * pad + 0.05
  dy <- diff(yr) * pad + 0.05
  list(
    xlim = c(xr[1] - dx, xr[2] + dx),
    ylim = c(yr[1] - dy, yr[2] + dy)
  )
}

plot_depth_map <- function(df, title, subtitle, outfile, show_dataset = FALSE,
                           map_box = NULL, allow_na_depth = FALSE) {
  df <- df %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
  if (!allow_na_depth) {
    df <- df %>% filter(is.finite(depth))
  }
  if (nrow(df) == 0) {
    message("Skip empty map: ", outfile)
    return(invisible(NULL))
  }

  df_depth <- df %>% filter(is.finite(depth))
  df_nodepth <- df %>% filter(!is.finite(depth))

  lim_df <- if (nrow(df_depth) > 0) df_depth else df
  if (!is.null(map_box)) {
    lim <- lim_from_map_box(map_box)
  } else {
    lim <- plot_limits(lim_df)
  }

  depth_limits <- if (nrow(df_depth) > 0) range(df_depth$depth, na.rm = TRUE) else c(0, 1)

  p <- add_basemap_layers(ggplot(), lim, map_box = map_box)

  if (nrow(df_depth) > 0) {
    pt_aes <- if (show_dataset) {
      aes(decimalLongitude, decimalLatitude, color = depth, shape = dataset)
    } else {
      aes(decimalLongitude, decimalLatitude, color = depth)
    }
    p <- p + geom_point(data = df_depth, mapping = pt_aes, alpha = 0.75, size = 1.8)
  }
  if (nrow(df_nodepth) > 0) {
    nd_aes <- if (show_dataset) {
      aes(decimalLongitude, decimalLatitude, shape = dataset)
    } else {
      aes(decimalLongitude, decimalLatitude)
    }
    p <- p +
      geom_point(
        data = df_nodepth, mapping = nd_aes,
        color = "#8b949e", fill = NA, shape = 1, size = 1.6, alpha = 0.85, stroke = 0.4
      )
  }

  if (nrow(df_depth) > 0) {
    p <- p +
      scale_color_viridis_c(option = "C", direction = -1, name = "Depth (m)",
                            limits = depth_limits)
  }
  p <- p +
    coord_fixed(1.3, xlim = lim$xlim, ylim = lim$ylim) +
    labs(title = title, subtitle = subtitle, x = "Longitude (°E)", y = "Latitude (°N)") +
    theme_dark_map()

  if (show_dataset) {
    p <- p + scale_shape_discrete(name = "Dataset")
  }

  ggsave(outfile, p, width = 9, height = 6.5, dpi = 150, bg = "#0d1117")
  message("Wrote ", outfile)
}

# Survey stations where other fish were recorded but focus species was not
other_species_stations <- function(survey_df, focus_sp_key) {
  if (nrow(survey_df) == 0) {
    return(tibble())
  }

  station_flags <- survey_df %>%
    group_by(dataset, obs_key) %>%
    summarise(
      has_focus = any(sp_key == focus_sp_key),
      .groups = "drop"
    )

  other_keys <- station_flags %>% filter(!has_focus)
  if (nrow(other_keys) == 0) {
    return(tibble())
  }

  survey_df %>%
    inner_join(other_keys, by = c("dataset", "obs_key")) %>%
    group_by(dataset, obs_key) %>%
    summarise(
      decimalLongitude = median(decimalLongitude, na.rm = TRUE),
      decimalLatitude = median(decimalLatitude, na.rm = TRUE),
      depth = median(depth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
}

plot_survey_context_map <- function(focus_df, other_df, title, subtitle, outfile,
                                    map_box = BBOX_JOINT) {
  focus_df <- focus_df %>%
    filter(is.finite(decimalLongitude), is.finite(decimalLatitude), is.finite(depth))
  if (nrow(focus_df) == 0) {
    message("Skip empty context map: ", outfile)
    return(invisible(NULL))
  }

  lim <- lim_from_map_box(map_box)
  depth_limits <- range(focus_df$depth, na.rm = TRUE)

  p <- add_basemap_layers(ggplot(), lim, map_box = map_box)

  if (nrow(other_df) > 0) {
    p <- p +
      geom_point(
        data = other_df,
        aes(decimalLongitude, decimalLatitude),
        color = "#6e7681",
        shape = 4,
        size = 1.5,
        alpha = 0.45,
        stroke = 0.35
      )
  }

  p <- p +
    geom_point(
      data = focus_df,
      aes(decimalLongitude, decimalLatitude, color = depth),
      alpha = 0.8,
      size = 2
    ) +
    scale_color_viridis_c(option = "C", direction = -1, name = "Depth (m)",
                          limits = depth_limits) +
    coord_fixed(1.3, xlim = lim$xlim, ylim = lim$ylim) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Longitude (°E)", y = "Latitude (°N)",
      caption = "Crosses: survey haul/station with other fish, focus species absent"
    ) +
    theme_dark_map()

  ggsave(outfile, p, width = 9, height = 6.5, dpi = 150, bg = "#0d1117")
  message("Wrote ", outfile)
}

datasets_present <- all_occ %>% distinct(dataset) %>% pull(dataset)

survey_occ <- all_occ %>%
  filter(
    dataset %in% SURVEY_DATASETS,
    is.finite(depth),
    !is.na(obs_key), obs_key != ""
  )

for (i in seq_len(nrow(pick5))) {
  sp <- pick5$scientificName[i]
  sp_key <- pick5$sp_key[i]
  sp_slug <- slug_species(sp)
  sp_dir <- file.path(OUTPUT_DIR, sp_slug)
  dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)

  sp_all <- all_occ %>%
    filter(sp_key == .env$sp_key) %>%
    thin_points(8000)

  sp_depth <- sp_all %>% filter(is.finite(depth))
  sp_survey_depth <- sp_all %>%
    filter(dataset %in% SURVEY_DATASETS, is.finite(depth))
  if (nrow(sp_survey_depth) == 0) {
    message("No survey depth for ", sp, " — skipping maps")
    next
  }

  # Per-dataset maps
  for (ds in datasets_present) {
    ds_df <- sp_all %>% filter(dataset == ds) %>% thin_points(3000)
    if (ds == GBIF_PAPER_DATASET) {
      ds_df <- ds_df %>% filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
    } else {
      ds_df <- ds_df %>% filter(is.finite(depth))
    }
    if (nrow(ds_df) == 0) next
    map_box <- if (ds == "MEDITS_raw") BBOX_MEDITS else BBOX_ISRAEL
    plot_depth_map(
      ds_df,
      title = sp,
      subtitle = paste0(
        "Dataset: ", ds,
        if (ds == GBIF_PAPER_DATASET) {
          " · filled = depth · open circles = depth missing"
        } else if (ds == COLLECTION_DATASET) {
          " · Singer & Owens (2026) museum/iDigBio · colour = mean depth (m)"
        } else {
          " · colour = depth (m)"
        },
        if (ds == "MEDITS_raw") " · GSA 22/25, Levantine" else ""
      ),
      outfile = file.path(sp_dir, paste0("map_", ds, ".png")),
      show_dataset = FALSE,
      map_box = map_box,
      allow_na_depth = (ds == GBIF_PAPER_DATASET)
    )
  }

  # Joint map (all datasets, shared depth colour; shape = dataset)
  plot_depth_map(
    sp_all %>% thin_points(6000),
    title = paste(sp, "— all datasets"),
    subtitle = "Joint map · colour = depth · Singer = museum vouchers · open circles = paper GBIF without depth",
    outfile = file.path(sp_dir, "map_ALL_datasets_joint.png"),
    show_dataset = TRUE,
    map_box = BBOX_JOINT,
    allow_na_depth = TRUE
  )

  # Survey co-occurrence: focus records + stations with other fish (focus absent)
  focus_survey <- survey_occ %>%
    filter(sp_key == .env$sp_key) %>%
    thin_points(6000)
  other_stations <- other_species_stations(survey_occ, sp_key) %>% thin_points(4000)
  if (nrow(other_stations) == 0) {
    message("No other-species-only stations for ", sp, " — context map may show focus only")
  }
  plot_survey_context_map(
    focus_df = focus_survey,
    other_df = other_stations,
    title = paste(sp, "— survey context"),
    subtitle = paste0(
      "MIT + Edelist + MEDITS · colour = focus depth (m) · ",
      nrow(other_stations), " hauls with other fish, not ", sp
    ),
    outfile = file.path(sp_dir, "map_joint_survey_context.png"),
    map_box = BBOX_JOINT
  )

  # Panel: one row per dataset (same species, same depth scale)
  ds_list <- sp_all %>%
    filter(
      (dataset %in% c(GBIF_PAPER_DATASET, COLLECTION_DATASET)) |
        (dataset %in% SURVEY_DATASETS & is.finite(depth))
    ) %>%
    count(dataset, sort = TRUE) %>%
    filter(n > 0) %>%
    pull(dataset)
  if (length(ds_list) >= 2) {
    panel_df <- sp_all %>%
      filter(
        dataset %in% ds_list,
        is.finite(decimalLongitude), is.finite(decimalLatitude),
        dataset %in% c(GBIF_PAPER_DATASET, COLLECTION_DATASET) | is.finite(depth)
      ) %>%
      thin_points(5000)
    land_f <- build_faceted_land(ds_list)
    lm_f <- build_faceted_landmarks(ds_list)
    p_panel <- ggplot() +
      geom_polygon(
        data = land_f,
        aes(x = long, y = lat, group = poly_group),
        fill = "#262b33",
        color = "#6e7681",
        linewidth = 0.25,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = lm_f,
        aes(x = decimalLongitude, y = decimalLatitude),
        inherit.aes = FALSE,
        shape = 3, size = 1.4, color = "#8b949e"
      ) +
      geom_text(
        data = lm_f,
        aes(x = decimalLongitude, y = decimalLatitude, label = label),
        inherit.aes = FALSE,
        size = 2.4,
        color = "#e6edf3",
        vjust = -0.6,
        check_overlap = TRUE
      ) +
      geom_point(
        data = panel_df,
        aes(decimalLongitude, decimalLatitude, color = depth),
        alpha = 0.75, size = 1.4
      ) +
      scale_color_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
      facet_wrap(~dataset, ncol = 2, scales = "free") +
      coord_cartesian(expand = TRUE) +
      labs(
        title = paste(sp, "— compare datasets"),
        subtitle = "Same species · depth colour · MEDITS = medits_raw_data.csv (Levantine GSA 22/25)",
        x = "Longitude (°E)", y = "Latitude (°N)"
      ) +
      theme_dark_map(10)
    ggsave(
      file.path(sp_dir, "map_panel_all_datasets.png"),
      p_panel, width = 11, height = max(4, 2.5 * ceiling(length(ds_list) / 2)),
      dpi = 150, bg = "#0d1117"
    )
    message("Wrote ", file.path(sp_dir, "map_panel_all_datasets.png"))
  }
}

message("\nDone. Maps in ", OUTPUT_DIR, "/<species>/")
message("Datasets used: ", paste(datasets_present, collapse = ", "))
