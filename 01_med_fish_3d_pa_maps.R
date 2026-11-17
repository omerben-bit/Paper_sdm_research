# =============================================================================
# Mediterranean fish — start 3D presence / absence maps
# Sources:
#   1) GBIF (Actinopterygii + Chondrichthyes in Med / Levant box)
#   2) Israeli Eastern Med trawl dataset (Edelist et al. via MedOBIS) — has depth
#   3) Optional local Michmoret CSV in data/michmoret/
#
# IMPORTANT: GBIF is mostly PRESENCE-ONLY.
# True absences come from surveys (Edelist hauls / Michmoret if coded).
# Elsewhere we use depth-binned BACKGROUND points as pseudo-absence for mapping.
# =============================================================================

options(timeout = 600)
options(pkgType = "binary")
options(install.packages.check.source = "no")
if (is.null(getOption("repos")) || identical(getOption("repos")["CRAN"], "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

need <- c(
  "rgbif", "dplyr", "tidyr", "ggplot2", "sf", "terra",
  "CoordinateCleaner", "readr", "stringr", "jsonlite"
)
miss <- need[!need %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss, type = "binary", dependencies = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(readr)
library(stringr)

# ---- study area: Eastern Med / Israel focus (expand later) ------------------
# Full Med roughly: lon -6..37, lat 30..46
# Start Levantine for faster GBIF + clearer maps
BBOX <- list(xmin = 32.0, xmax = 35.5, ymin = 31.0, ymax = 34.0) # Israel EEZ-ish
DEPTH_BINS <- c(0, 10, 30, 50, 100, 200, 500, 1000) # metres

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs_med_3d_pa", showWarnings = FALSE, recursive = TRUE)

# ---- helpers ----------------------------------------------------------------
standardize_occ <- function(df, source_name) {
  nm <- names(df)
  pick <- function(...) {
    cands <- c(...)
    hit <- cands[cands %in% nm]
    if (length(hit)) hit[[1]] else NA_character_
  }
  lon_c <- pick("decimalLongitude", "decimallongitude", "longitude", "lon", "x")
  lat_c <- pick("decimalLatitude", "decimallatitude", "latitude", "lat", "y")
  sp_c  <- pick("scientificName", "species", "acceptedNameUsage", "scientificname")
  dep_c <- pick(
    "depth", "depthInMeters", "minimumDepthInMeters", "verbatimDepth",
    "coordinateUncertaintyInMeters" # never use this as depth — listed only as last-resort guard
  )
  # Prefer true depth fields only
  dep_c <- pick("depth", "depthInMeters", "minimumDepthInMeters", "verbatimDepth")
  st_c  <- pick("occurrenceStatus", "occurrencestatus")
  dt_c  <- pick("eventDate", "eventdate", "year")

  out <- tibble::tibble(
    scientificName = if (!is.na(sp_c)) as.character(df[[sp_c]]) else NA_character_,
    decimalLongitude = if (!is.na(lon_c)) as.numeric(df[[lon_c]]) else NA_real_,
    decimalLatitude = if (!is.na(lat_c)) as.numeric(df[[lat_c]]) else NA_real_,
    depth = if (!is.na(dep_c)) suppressWarnings(as.numeric(df[[dep_c]])) else NA_real_,
    occurrenceStatus = if (!is.na(st_c)) tolower(as.character(df[[st_c]])) else "present",
    eventDate = if (!is.na(dt_c)) as.character(df[[dt_c]]) else NA_character_,
    source = source_name
  )

  out <- out %>%
    mutate(
      occurrenceStatus = case_when(
        occurrenceStatus %in% c("present", "presence", "1", "true") ~ "present",
        occurrenceStatus %in% c("absent", "absence", "0", "false") ~ "absent",
        TRUE ~ "present"
      ),
      depth_bin = cut(
        depth,
        breaks = DEPTH_BINS,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    filter(
      !is.na(scientificName), scientificName != "",
      !is.na(decimalLongitude), !is.na(decimalLatitude),
      decimalLongitude >= BBOX$xmin, decimalLongitude <= BBOX$xmax,
      decimalLatitude >= BBOX$ymin, decimalLatitude <= BBOX$ymax
    )
  out
}

# ---- 1) GBIF fish (presence-only) -------------------------------------------
message("=== GBIF download (Actinopterygii in bbox) ===")
gbif_rds <- "data/raw/gbif_fish_levant.rds"

if (!file.exists(gbif_rds)) {
  # Class Actinopterygii taxonKey on GBIF is 204
  # Chondrichthyes = 121
  occ_bone <- rgbif::occ_search(
    taxonKey = 204,
    geometry = sprintf(
      "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
      BBOX$xmin, BBOX$ymin,
      BBOX$xmax, BBOX$ymin,
      BBOX$xmax, BBOX$ymax,
      BBOX$xmin, BBOX$ymax,
      BBOX$xmin, BBOX$ymin
    ),
    hasCoordinate = TRUE,
    limit = 5000,
    basisOfRecord = "HUMAN_OBSERVATION;MACHINE_OBSERVATION;OBSERVATION;OCCURRENCE;PRESERVED_SPECIMEN"
  )
  occ_cart <- rgbif::occ_search(
    taxonKey = 121,
    geometry = sprintf(
      "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
      BBOX$xmin, BBOX$ymin,
      BBOX$xmax, BBOX$ymin,
      BBOX$xmax, BBOX$ymax,
      BBOX$xmin, BBOX$ymax,
      BBOX$xmin, BBOX$ymin
    ),
    hasCoordinate = TRUE,
    limit = 2000
  )
  gbif_df <- bind_rows(occ_bone$data, occ_cart$data)
  saveRDS(gbif_df, gbif_rds)
  message("Saved ", nrow(gbif_df), " GBIF rows -> ", gbif_rds)
} else {
  gbif_df <- readRDS(gbif_rds)
  message("Loaded cached GBIF: ", nrow(gbif_df), " rows")
}

gbif_std <- standardize_occ(gbif_df, "GBIF")
# Clean obviously bad coords
if (requireNamespace("CoordinateCleaner", quietly = TRUE) && nrow(gbif_std) > 0) {
  gbif_std$species_tmp <- gbif_std$scientificName
  flags <- tryCatch(
    CoordinateCleaner::clean_coordinates(
      x = as.data.frame(gbif_std),
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "species_tmp",
      tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
      value = "flagged"
    ),
    error = function(e) rep(TRUE, nrow(gbif_std))
  )
  gbif_std <- gbif_std[flags, ] %>% dplyr::select(-species_tmp)
}

write_csv(gbif_std, "data/processed/gbif_fish_levant_cleaned.csv")
message("Wrote paper-cleaned GBIF -> data/processed/gbif_fish_levant_cleaned.csv (", nrow(gbif_std), " rows)")

message("GBIF kept: ", nrow(gbif_std),
        " | with depth: ", sum(!is.na(gbif_std$depth)))

# ---- 2) Israeli trawl / OBIS (Edelist) — better for depth + effort ----------
message("=== Edelist Israeli Med fish (MedOBIS DwC-A) ===")
edelist_csv <- "data/raw/edelist_occurrence.txt"
edelist_zip <- "data/raw/edelist_medobis.zip"
edelist_url <- "https://ipt.medobis.eu/archive.do?r=fish_invertebrates_israelimediterraneansea"

if (!file.exists(edelist_csv)) {
  message("Downloading MedOBIS archive...")
  tryCatch({
    download.file(edelist_url, edelist_zip, mode = "wb", quiet = FALSE)
    utils::unzip(edelist_zip, exdir = "data/raw/edelist_dwca")
    occ_file <- list.files(
      "data/raw/edelist_dwca",
      pattern = "occurrence",
      full.names = TRUE,
      ignore.case = TRUE
    )[1]
    if (is.na(occ_file)) stop("No occurrence file in DwC-A")
    file.copy(occ_file, edelist_csv, overwrite = TRUE)
  }, error = function(e) {
    message("Edelist download failed: ", conditionMessage(e))
    message("You can manually download from OBIS and place occurrence.txt at:\n  ", edelist_csv)
  })
}

edelist_std <- tibble()
if (file.exists(edelist_csv)) {
  edelist_raw <- readr::read_tsv(edelist_csv, show_col_types = FALSE, guess_max = 50000)
  # Keep fish-like class / phylum if columns exist
  if ("class" %in% names(edelist_raw)) {
    edelist_raw <- edelist_raw %>%
      filter(is.na(class) | str_detect(tolower(class), "actinopteryg|chondrichth|elasmobranch|osteichth"))
  }
  if ("kingdom" %in% names(edelist_raw)) {
    edelist_raw <- edelist_raw %>% filter(is.na(kingdom) | tolower(kingdom) == "animalia")
  }
  # Depth: Darwin Core often has minimumDepthInMeters / maximumDepthInMeters
  if (!"depth" %in% names(edelist_raw)) {
    if (all(c("minimumDepthInMeters", "maximumDepthInMeters") %in% names(edelist_raw))) {
      edelist_raw$depth <- rowMeans(
        cbind(
          as.numeric(edelist_raw$minimumDepthInMeters),
          as.numeric(edelist_raw$maximumDepthInMeters)
        ),
        na.rm = TRUE
      )
      edelist_raw$depth[!is.finite(edelist_raw$depth)] <- NA_real_
    } else if ("minimumDepthInMeters" %in% names(edelist_raw)) {
      edelist_raw$depth <- as.numeric(edelist_raw$minimumDepthInMeters)
    }
  }
  edelist_std <- standardize_occ(edelist_raw, "Edelist_OBIS_Israel")
  message("Edelist kept: ", nrow(edelist_std),
          " | with depth: ", sum(!is.na(edelist_std$depth)))
} else {
  message("Skipping Edelist (file missing).")
}

# ---- 3) Optional Michmoret local CSV ----------------------------------------
message("=== Michmoret local ===")
mich_files <- list.files("data/michmoret", pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
mich_std <- tibble()
if (length(mich_files)) {
  mf <- mich_files[[1]]
  message("Reading ", mf)
  mich_raw <- if (str_detect(mf, "\\.tsv$")) {
    readr::read_tsv(mf, show_col_types = FALSE)
  } else {
    readr::read_csv(mf, show_col_types = FALSE)
  }
  mich_std <- standardize_occ(mich_raw, "Michmoret_local")
  message("Michmoret kept: ", nrow(mich_std),
          " | absences: ", sum(mich_std$occurrenceStatus == "absent"))
} else {
  message("No CSV in data/michmoret/ yet — add your file when ready (see README).")
}

# ---- combine ----------------------------------------------------------------
occ_all <- bind_rows(gbif_std, edelist_std, mich_std)
write_csv(occ_all, "data/processed/med_fish_occurrences_combined.csv")
message("Combined rows: ", nrow(occ_all))

# Species ranking for first maps
top_species <- occ_all %>%
  filter(occurrenceStatus == "present") %>%
  count(scientificName, sort = TRUE) %>%
  slice_head(n = 12)
print(top_species)
write_csv(top_species, "outputs_med_3d_pa/top_species.csv")

# Pick first mapped species = most records with depth
sp_focus <- occ_all %>%
  filter(occurrenceStatus == "present", !is.na(depth)) %>%
  count(scientificName, sort = TRUE) %>%
  slice_head(n = 1) %>%
  pull(scientificName)

if (!length(sp_focus)) {
  sp_focus <- top_species$scientificName[[1]]
  message("No depth for top species; mapping without depth bins for: ", sp_focus)
} else {
  message("Focus species for first 3D PA maps: ", sp_focus)
}

# ---- build presence / absence table for focus species -----------------------
# Presence: records of the species
# Absence candidates:
#   A) explicit absences in Michmoret / Edelist if any
#   B) other survey stations (Edelist) in same depth bin where species not listed
#   C) random background in bbox (pseudo-absence) — labelled clearly

pres <- occ_all %>%
  filter(scientificName == sp_focus, occurrenceStatus == "present") %>%
  mutate(pa = 1L, pa_type = "presence")

# Survey stations from Edelist: unique haul locations (+ depth) without the species
edelist_stations <- edelist_std %>%
  filter(!is.na(depth)) %>%
  mutate(
    station_id = paste(
      round(decimalLongitude, 3),
      round(decimalLatitude, 3),
      depth_bin,
      sep = "_"
    )
  )

stations_with_sp <- edelist_stations %>%
  filter(scientificName == sp_focus) %>%
  distinct(station_id)

survey_abs <- edelist_stations %>%
  filter(!station_id %in% stations_with_sp$station_id) %>%
  distinct(decimalLongitude, decimalLatitude, depth, depth_bin, station_id) %>%
  mutate(
    scientificName = sp_focus,
    occurrenceStatus = "absent",
    eventDate = NA_character_,
    source = "Edelist_survey_absence",
    pa = 0L,
    pa_type = "survey_absence"
  )

# Pseudo-absence background (only where we lack enough survey absences)
set.seed(42)
n_pseudo <- max(200, 3 * nrow(pres))
pseudo <- tibble(
  scientificName = sp_focus,
  decimalLongitude = runif(n_pseudo, BBOX$xmin, BBOX$xmax),
  decimalLatitude = runif(n_pseudo, BBOX$ymin, BBOX$ymax),
  depth = if (all(is.na(pres$depth))) {
    sample(c(10, 25, 45, 80, 120), n_pseudo, replace = TRUE)
  } else {
    sample(pres$depth[!is.na(pres$depth)], n_pseudo, replace = TRUE)
  },
  occurrenceStatus = "absent",
  eventDate = NA_character_,
  source = "pseudo_absence_background",
  pa = 0L,
  pa_type = "pseudo_absence"
) %>%
  mutate(depth_bin = cut(depth, breaks = DEPTH_BINS, include.lowest = TRUE, right = FALSE))

pa_df <- bind_rows(
  pres %>% mutate(station_id = NA_character_),
  survey_abs,
  pseudo
)
write_csv(pa_df, "data/processed/pa_focus_species.csv")
write_csv(pa_df, file.path("outputs_med_3d_pa", paste0("pa_", make.names(sp_focus), ".csv")))

message(
  "PA table for ", sp_focus, ":\n",
  "  presence: ", sum(pa_df$pa == 1), "\n",
  "  survey absence: ", sum(pa_df$pa_type == "survey_absence"), "\n",
  "  pseudo-absence: ", sum(pa_df$pa_type == "pseudo_absence")
)

# ---- MAP 1: all sources overview --------------------------------------------
p_overview <- ggplot(
  occ_all %>% filter(occurrenceStatus == "present"),
  aes(decimalLongitude, decimalLatitude, color = source)
) +
  geom_point(alpha = 0.45, size = 1.2) +
  coord_fixed(1.2, xlim = c(BBOX$xmin, BBOX$xmax), ylim = c(BBOX$ymin, BBOX$ymax)) +
  labs(
    title = "Fish occurrence records — Levantine / Israel box",
    subtitle = "GBIF + Edelist OBIS (+ Michmoret if provided)",
    x = "Longitude", y = "Latitude", color = "Source"
  ) +
  theme_bw()
ggsave("outputs_med_3d_pa/01_overview_sources.png", p_overview, width = 8, height = 6, dpi = 140)

# ---- MAP 2: focus species presence vs survey absence (2D) -------------------
map_pts <- pa_df %>%
  filter(pa_type %in% c("presence", "survey_absence")) %>%
  mutate(pa_type = factor(pa_type, levels = c("presence", "survey_absence")))

if (nrow(map_pts) > 0) {
  p_pa2d <- ggplot(map_pts, aes(decimalLongitude, decimalLatitude, color = pa_type)) +
    geom_point(alpha = 0.7, size = 1.6) +
    scale_color_manual(values = c(presence = "#1b4332", survey_absence = "#bc4749")) +
    coord_fixed(1.2, xlim = c(BBOX$xmin, BBOX$xmax), ylim = c(BBOX$ymin, BBOX$ymax)) +
    labs(
      title = paste("Presence / survey-absence —", sp_focus),
      subtitle = "Absences = Edelist haul stations where species was not recorded",
      x = "Longitude", y = "Latitude", color = NULL
    ) +
    theme_bw()
  ggsave("outputs_med_3d_pa/02_pa_2d_focus.png", p_pa2d, width = 8, height = 6, dpi = 140)
}

# ---- MAP 3: 3D-style — facet by depth bin -----------------------------------
pa_depth <- pa_df %>%
  filter(pa_type %in% c("presence", "survey_absence"), !is.na(depth_bin))

if (nrow(pa_depth) > 0) {
  p_3d <- ggplot(pa_depth, aes(decimalLongitude, decimalLatitude, color = pa_type)) +
    geom_point(alpha = 0.75, size = 1.4) +
    scale_color_manual(values = c(presence = "#1b4332", survey_absence = "#bc4749")) +
    facet_wrap(~depth_bin) +
    coord_fixed(1.2) +
    labs(
      title = paste("3D presence / absence by depth bin —", sp_focus),
      subtitle = "Each panel is a depth layer (metres). This is the start of a 3D PA map.",
      x = "Longitude", y = "Latitude", color = NULL
    ) +
    theme_bw(base_size = 11)
  ggsave("outputs_med_3d_pa/03_pa_by_depth.png", p_3d, width = 11, height = 8, dpi = 140)
}

# ---- MAP 4: presence density by depth (all fish species) --------------------
with_depth <- occ_all %>%
  filter(occurrenceStatus == "present", !is.na(depth), !is.na(depth_bin))

if (nrow(with_depth) > 0) {
  p_hist <- ggplot(with_depth, aes(depth)) +
    geom_histogram(bins = 30, fill = "#2a6f97", color = "white") +
    labs(
      title = "Depth distribution of fish records (all species in box)",
      subtitle = "Records without depth are excluded — GBIF often lacks depth",
      x = "Depth (m)", y = "Records"
    ) +
    theme_bw()
  ggsave("outputs_med_3d_pa/04_depth_histogram_all.png", p_hist, width = 7.5, height = 4.2, dpi = 140)

  # Top 6 species × depth facets presence only
  top6 <- with_depth %>% count(scientificName, sort = TRUE) %>% slice_head(n = 6) %>% pull(scientificName)
  p_sp_depth <- with_depth %>%
    filter(scientificName %in% top6) %>%
    ggplot(aes(decimalLongitude, decimalLatitude)) +
    geom_point(aes(color = depth), alpha = 0.8, size = 1.5) +
    scale_color_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
    facet_wrap(~scientificName) +
    coord_fixed(1.2) +
    labs(
      title = "Top species in the box — coloured by depth",
      x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_size = 10)
  ggsave("outputs_med_3d_pa/05_top_species_depth.png", p_sp_depth, width = 10, height = 7, dpi = 140)
}

# ---- summary ----------------------------------------------------------------
summary_tbl <- tibble(
  n_gbif = nrow(gbif_std),
  n_edelist = nrow(edelist_std),
  n_michmoret = nrow(mich_std),
  n_combined = nrow(occ_all),
  n_with_depth = sum(!is.na(occ_all$depth)),
  focus_species = sp_focus,
  n_presence = sum(pa_df$pa == 1),
  n_survey_absence = sum(pa_df$pa_type == "survey_absence"),
  n_pseudo_absence = sum(pa_df$pa_type == "pseudo_absence")
)
print(summary_tbl)
write_csv(summary_tbl, "outputs_med_3d_pa/summary.csv")

message("\nDone. Open folder: outputs_med_3d_pa/")
message("Next: drop Michmoret CSV into data/michmoret/ and re-run.")
message("Then we can fit a depth-aware SDM (voluModel / mgcv) on this PA table.")
