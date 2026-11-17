# Canonical GBIF layer: `01_med_fish_3d_pa_maps.R` (bbox + CoordinateCleaner, presence-only).
# Do not use `data/raw/gbif_fish_levant.rds` for analysis (uncleaned download cache).

if (!isNamespaceLoaded("dplyr")) {
  suppressPackageStartupMessages(library(dplyr))
}

GBIF_PAPER_CLEANED_CSV <- "data/processed/gbif_fish_levant_cleaned.csv"
GBIF_PAPER_COMBINED_CSV <- "data/processed/med_fish_occurrences_combined.csv"
GBIF_PAPER_DATASET <- "GBIF_paper_cleaned"

gbif_paper_is_present <- function(status) {
  st <- tolower(as.character(status))
  st %in% c("present", "presence", "1", "true")
}

#' @param fallback_combined If cleaned CSV missing, read GBIF rows from combined export.
load_gbif_paper_cleaned <- function(fallback_combined = TRUE) {
  path <- GBIF_PAPER_CLEANED_CSV
  if (file.exists(path)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (!fallback_combined || !file.exists(GBIF_PAPER_COMBINED_CSV)) {
    warning("GBIF paper-cleaned file missing; run 01_med_fish_3d_pa_maps.R")
    return(tibble::tibble())
  }
  readr::read_csv(GBIF_PAPER_COMBINED_CSV, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(source == "GBIF", gbif_paper_is_present(occurrenceStatus))
}