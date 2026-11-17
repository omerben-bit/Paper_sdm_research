# =============================================================================
# Species co-occurrence (haul / observation units): positive vs negative pairs
# Output: species_selection/table_positive_correlations.csv
#         species_selection/table_negative_correlations.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

OUT_DIR <- "species_selection"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BBOX_ISRAEL <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 29.5, lat_max = 34.5)
BBOX_MEDITS <- list(lon_min = 32.0, lon_max = 36.5, lat_min = 31.0, lat_max = 35.5)

MIN_PRESENCE_HAULS <- 15L   # each species in >= this many observation units
MIN_POS_R <- 0.20            # keep positive pairs with r >= this
MIN_NEG_R <- 0.10            # keep negative pairs with r <= -this (rarer in trawl data)
TOP_N_EACH <- 500L           # cap rows per table (after filters)

norm_species <- function(x) {
  x %>% as.character() %>% str_trim() %>%
    str_remove("\\s*\\([^)]*\\)$") %>% str_squish()
}

in_bbox <- function(df, box) {
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

medits_haul_key <- function(country, area, vessel, year, haul_number) {
  paste(country, area, vessel, year, haul_number, sep = "|")
}

load_survey_records <- function() {
  parts <- list()

  if (file.exists("Combined_Menachem_Itai_corrected_May_2026_DBG.csv")) {
    parts[[length(parts) + 1]] <- read_csv(
      "Combined_Menachem_Itai_corrected_May_2026_DBG.csv", show_col_types = FALSE
    ) %>%
      transmute(
        scientificName = as.character(ID_SPECIES),
        decimalLatitude = as.numeric(LATITUDE),
        decimalLongitude = as.numeric(LONGITUDE),
        dataset = "Israel_trawl_MIT",
        obs_key = as.character(SAMPLE_DESC)
      ) %>%
      filter(!is.na(scientificName), scientificName != "") %>%
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
        scientificName = as.character(scientificName),
        decimalLatitude = as.numeric(meta[, 2]),
        decimalLongitude = as.numeric(meta[, 3]),
        dataset = "Edelist_OBIS",
        obs_key = str_replace(.data[[id_col]], "_Fish_.*$", "")
      ) %>%
      filter(!is.na(scientificName), scientificName != "") %>%
      in_bbox(BBOX_ISRAEL)
  }

  if (file.exists("MEDITS/medits_raw_data.csv") && file.exists("MEDITS/hauls_data.csv")) {
    catch <- read_csv("MEDITS/medits_raw_data.csv", show_col_types = FALSE,
                      locale = locale(encoding = "Latin1")) %>%
      filter(nbtot > 0)
    hauls <- read_csv("MEDITS/hauls_data.csv", show_col_types = FALSE) %>%
      filter(validity == "V") %>%
      mutate(
        decimalLatitude = medits_dd(shooting_latitude),
        decimalLongitude = medits_dd(shooting_longitude),
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
          scientificName = .data[[name_col[1]]],
          CATFAU = .data[[names(taxa)[str_detect(names(taxa), "CATFAU")][1]]]
        ) %>%
        select(medits_code, scientificName, CATFAU)
      catch <- catch %>%
        mutate(medits_code = paste0(genus, species)) %>%
        left_join(taxa, by = "medits_code") %>%
        filter(is.na(CATFAU) | CATFAU %in% c("Ao", "Ae", "C"))
    } else {
      catch <- catch %>% mutate(scientificName = paste0(genus, species))
    }
    parts[[length(parts) + 1]] <- catch %>%
      inner_join(
        hauls %>% select(country, area, vessel, year, haul_number,
                         decimalLatitude, decimalLongitude, obs_key),
        by = c("country", "area", "vessel", "year", "haul_number")
      ) %>%
      transmute(
        scientificName = coalesce(scientificName, medits_code),
        decimalLatitude, decimalLongitude,
        dataset = "MEDITS_raw",
        obs_key,
        gsa = area
      ) %>%
      filter(gsa %in% c(22L, 25L)) %>%
      in_bbox(BBOX_MEDITS) %>%
      select(-gsa)
  }

  if (!exists("load_gbif_paper_cleaned", mode = "function")) {
    source("gbif_paper.R")
  }
  gbif_paper <- load_gbif_paper_cleaned()
  if (nrow(gbif_paper) > 0) {
    parts[[length(parts) + 1]] <- gbif_paper %>%
      transmute(
        scientificName = norm_species(scientificName),
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        dataset = GBIF_PAPER_DATASET,
        obs_key = NA_character_
      ) %>%
      filter(!is.na(scientificName), scientificName != "")
  }

  bind_rows(parts) %>%
    mutate(sp_key = norm_species(scientificName)) %>%
    filter(
      !sp_key %in% c("", "Unaccepted", "Equivocal"),
      is.finite(decimalLatitude), is.finite(decimalLongitude)
    )
}

assign_obs_id <- function(df) {
  df %>%
    mutate(
      obs_id = if_else(
        !is.na(obs_key) & obs_key != "",
        paste0(dataset, "::", obs_key),
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

pairwise_cor_table <- function(pres_mat, sp_names) {
  p <- ncol(pres_mat)
  n_obs <- nrow(pres_mat)
  if (p < 2L) {
    return(tibble())
  }

  cors <- stats::cor(pres_mat, use = "pairwise.complete.obs")
  hauls_present <- colSums(pres_mat)

  rows <- list()
  k <- 0L
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      r <- cors[i, j]
      if (is.na(r)) {
        next
      }
      co <- sum(pres_mat[, i] == 1L & pres_mat[, j] == 1L)
      k <- k + 1L
      rows[[k]] <- tibble(
        species_a = sp_names[i],
        species_b = sp_names[j],
        r = round(r, 4),
        n_observations = n_obs,
        hauls_a = hauls_present[i],
        hauls_b = hauls_present[j],
        co_presence_hauls = co
      )
    }
  }
  if (k == 0L) {
    return(tibble())
  }
  bind_rows(rows)
}

message("Loading records...")
dat <- load_survey_records()
if (nrow(dat) == 0) {
  stop("No records loaded.")
}

dat <- assign_obs_id(dat)

message("Building presence matrix (observation units x species)...")
presence <- dat %>%
  distinct(obs_id, dataset, sp_key) %>%
  mutate(present = 1L)

sp_counts <- presence %>%
  count(sp_key, name = "n_hauls") %>%
  filter(n_hauls >= MIN_PRESENCE_HAULS) %>%
  arrange(desc(n_hauls))

name_pick <- dat %>%
  count(sp_key, scientificName) %>%
  group_by(sp_key) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

kept_sp <- sp_counts$sp_key
message("Species with >= ", MIN_PRESENCE_HAULS, " presence units: ", length(kept_sp))

obs_ids <- unique(presence$obs_id)
mat_long <- presence %>%
  filter(sp_key %in% kept_sp) %>%
  select(obs_id, sp_key, present)

mat_wide <- mat_long %>%
  pivot_wider(names_from = sp_key, values_from = present, values_fill = 0L)

obs_order <- unique(dat$obs_id)
mat_wide <- mat_wide %>%
  right_join(tibble(obs_id = obs_order), by = "obs_id")
mat_wide[is.na(mat_wide)] <- 0L

pres_mat <- as.matrix(mat_wide %>% select(-obs_id))
sp_cols <- colnames(pres_mat)
sp_labels <- name_pick %>%
  filter(sp_key %in% sp_cols) %>%
  arrange(match(sp_key, sp_cols)) %>%
  pull(scientificName)

if (length(sp_labels) != ncol(pres_mat)) {
  sp_labels <- sp_cols
}

message("Computing pairwise correlations (Pearson on presence/absence)...")
pairs <- pairwise_cor_table(pres_mat, sp_labels)

meta <- tibble(
  metric = c(
    "observation_units", "datasets", "species_in_matrix",
    "min_presence_hauls", "min_pos_r", "min_neg_r",
    "pairs_exported_positive", "pairs_exported_negative"
  ),
  value = c(
    nrow(pres_mat),
    n_distinct(dat$dataset),
    ncol(pres_mat),
    MIN_PRESENCE_HAULS,
    MIN_POS_R,
    MIN_NEG_R,
    NA_real_,
    NA_real_
  )
)

pos <- pairs %>%
  filter(r >= MIN_POS_R) %>%
  arrange(desc(r), desc(co_presence_hauls)) %>%
  slice_head(n = TOP_N_EACH) %>%
  mutate(rank = row_number())

neg <- pairs %>%
  filter(r <= -MIN_NEG_R) %>%
  arrange(r, desc(co_presence_hauls)) %>%
  slice_head(n = TOP_N_EACH) %>%
  mutate(rank = row_number())

meta$value[meta$metric == "pairs_exported_positive"] <- nrow(pos)
meta$value[meta$metric == "pairs_exported_negative"] <- nrow(neg)
write_csv(meta, file.path(OUT_DIR, "run_metadata.csv"))

write_csv(pos, file.path(OUT_DIR, "table_positive_correlations.csv"))
write_csv(neg, file.path(OUT_DIR, "table_negative_correlations.csv"))

message("Wrote ", file.path(OUT_DIR, "table_positive_correlations.csv"),
        " (", nrow(pos), " rows)")
message("Wrote ", file.path(OUT_DIR, "table_negative_correlations.csv"),
        " (", nrow(neg), " rows)")
message("Done.")
