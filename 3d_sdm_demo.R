# =============================================================================
# 3D Species Distribution Modelling — hands-on demo (R)
# =============================================================================
# Goal: play with the core stack for a reusable multi-species 3D SDM:
#   voluModel, terra, sf, mgcv / glm, blockCV, (+ optional sdmTMB, biomod2)
#
# Data: built-in voluModel example (Luminous Hake + WOA temperature/oxygen).
#       Swap later for Mediterranean MEDITS + CMEMS with the same pipeline.
#
# Run section-by-section in RStudio (Ctrl+Enter), or:
#   source("3d_sdm_demo.R", echo = TRUE)
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
needed <- c(
  "voluModel", "terra", "sf", "dplyr", "tidyr", "ggplot2", "tidyterra",
  "CoordinateCleaner", "blockCV", "mgcv"
)
optional <- c("sdmTMB", "biomod2", "rnaturalearth", "rnaturalearthdata")

install_if_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing)) {
    message("Installing: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
}

install_if_missing(needed)
# Optional packages — uncomment when you want them:
# install_if_missing(optional)

library(voluModel)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyterra)
library(mgcv)

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Output folder for plots / small tables
out_dir <- file.path("outputs_3d_sdm_demo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ---- 1. Occurrences (voluModel extdata) -------------------------------------
# Columns: decimalLongitude, decimalLatitude, depth  (depth in metres)
occs_raw <- read.csv(
  system.file("extdata/Steindachneria_argentea.csv", package = "voluModel")
)

occurrences <- occs_raw %>%
  dplyr::select(decimalLongitude, decimalLatitude, depth) %>%
  dplyr::filter(!is.na(depth), depth >= 1, depth <= 2000) %>%
  dplyr::distinct()

message("Occurrence rows: ", nrow(occurrences))
print(head(occurrences))

# Optional cleaning pattern (useful later for GBIF/OBIS Mediterranean pulls)
if (has_pkg("CoordinateCleaner")) {
  library(CoordinateCleaner)
  occ_cc <- occurrences %>%
    mutate(species = "Steindachneria_argentea")
  flags <- tryCatch(
    clean_coordinates(
      x = as.data.frame(occ_cc),
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "species",
      tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
      value = "flagged"
    ),
    error = function(e) {
      message("CoordinateCleaner skipped: ", conditionMessage(e))
      rep(TRUE, nrow(occ_cc))
    }
  )
  message("CoordinateCleaner kept: ", sum(flags), " / ", length(flags))
  occurrences <- occurrences[flags, , drop = FALSE]
}

# Land polygons for maps (fallback if rnaturalearth missing)
if (has_pkg("rnaturalearth")) {
  land <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
} else {
  # Minimal empty sf — maps still work without fancy coastlines
  land <- st_sf(geometry = st_sfc(st_polygon()), crs = 4326)
}

# ---- 2. Environmental rasters (temperature + AOU depth stacks) --------------
td <- tempdir()
unzip(
  system.file("extdata/woa18_decav_t00mn01_cropped.zip", package = "voluModel"),
  exdir = file.path(td, "temperature"),
  junkpaths = TRUE
)

temperature_pts <- vect(file.path(td, "temperature", "woa18_decav_t00mn01_cropped.shp"))
temperature_bottom <- bottomRaster(temperature_pts)
template <- centerPointRasterTemplate(temperature_pts)
temperature <- terra::rasterize(
  x = temperature_pts, y = template, field = names(temperature_pts)
)

envt_names <- gsub("[d,M]", "", names(temperature_pts))
envt_names[[1]] <- "0"
names(temperature) <- envt_names

oxygen <- rast(system.file("extdata/oxygenSmooth.tif", package = "voluModel"))
oxygen_bottom <- rast(system.file("extdata/oxygenBottom.tif", package = "voluModel"))
names(oxygen) <- names(temperature)

message("Temperature stack layers (depths m): ", paste(names(temperature), collapse = ", "))
message("Raster resolution (deg): ", paste(res(temperature), collapse = " x "))

# Quick surface vs bottom plot
png(file.path(out_dir, "01_temp_surface_vs_bottom.png"), width = 1100, height = 500)
par(mfrow = c(1, 2))
oneRasterPlot(temperature[[1]], land = land, landCol = "grey20",
              title = "Surface temperature (°C)")
oneRasterPlot(temperature_bottom, land = land, landCol = "grey20",
              title = "Bottom temperature (°C)")
dev.off()

# ---- 3. Match each occurrence to nearest depth layer + downsample -----------
layer_depths <- as.numeric(names(temperature))
occurrences$index <- vapply(
  occurrences$depth,
  FUN = function(z) which.min(abs(layer_depths - z)),
  FUN.VALUE = integer(1)
)

downsampled <- lapply(sort(unique(occurrences$index)), function(i) {
  pts <- occurrences[occurrences$index == i, ]
  pts <- downsample(pts, temperature[[1]], verbose = FALSE)
  pts$depth <- layer_depths[[i]]
  pts
})
occurrences <- bind_rows(downsampled) %>%
  dplyr::select(decimalLongitude, decimalLatitude, depth)

message("Downsampled occurrence voxels: ", nrow(occurrences))

# ---- 4. THE 3D idea: extract env at occurrence depth ------------------------
# Surface / bottom / true 3D extractions — compare distributions
vals_surface <- data.frame(
  temperature = terra::extract(
    temperature[[1]], occurrences[, c("decimalLongitude", "decimalLatitude")]
  )[, 2],
  aou = terra::extract(
    oxygen[[1]], occurrences[, c("decimalLongitude", "decimalLatitude")]
  )[, 2],
  method = "surface"
)

vals_bottom <- data.frame(
  temperature = terra::extract(
    temperature_bottom, occurrences[, c("decimalLongitude", "decimalLatitude")]
  )[, 2],
  aou = terra::extract(
    oxygen_bottom, occurrences[, c("decimalLongitude", "decimalLatitude")]
  )[, 2],
  method = "bottom"
)

vals_3d <- data.frame(
  temperature = xyzSample(occs = occurrences, temperature),
  aou = xyzSample(occs = occurrences, oxygen),
  method = "3D_at_depth"
)

compare_extract <- bind_rows(vals_surface, vals_bottom, vals_3d) %>%
  pivot_longer(c(temperature, aou), names_to = "variable", values_to = "value")

p_extract <- ggplot(compare_extract, aes(method, value, fill = method)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    title = "Why 3D extraction matters",
    subtitle = "Same lon/lat points; env taken at surface vs bottom vs true depth",
    x = NULL, y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "02_extract_surface_bottom_3d.png"),
       p_extract, width = 8, height = 4.5)

# Modelling table: presence points with 3D covariates
pres <- occurrences %>%
  mutate(
    temperature = xyzSample(occs = occurrences, temperature),
    aou = xyzSample(occs = occurrences, oxygen),
    presence = 1L
  )
pres <- pres[stats::complete.cases(pres), , drop = FALSE]

# ---- 5. Study region + 3D background / pseudo-absences ----------------------
study_region <- marineBackground(
  occurrences,
  buff = 1e6,
  clipToOcean = TRUE,
  alpha = 1,
  partCount = 1
)

# Background points: sample xy in study region, random depth from observed range
n_bg <- min(2000L, 5L * nrow(pres))
bg_xy <- spatSample(study_region, size = n_bg, method = "random")
# spatSample on SpatVector returns a data.frame of coordinates in recent terra
if (inherits(bg_xy, "SpatVector")) {
  bg_coords <- crds(bg_xy)
} else {
  bg_coords <- as.matrix(bg_xy[, 1:2])
}
colnames(bg_coords) <- c("decimalLongitude", "decimalLatitude")

bg <- as.data.frame(bg_coords) %>%
  mutate(
    depth = sample(pres$depth, size = n(), replace = TRUE),
    presence = 0L
  )

bg$temperature <- xyzSample(occs = bg, temperature)
bg$aou <- xyzSample(occs = bg, oxygen)
bg <- bg[stats::complete.cases(bg), , drop = FALSE]

model_df <- bind_rows(
  pres %>% dplyr::select(decimalLongitude, decimalLatitude, depth,
                         temperature, aou, presence),
  bg %>% dplyr::select(decimalLongitude, decimalLatitude, depth,
                       temperature, aou, presence)
)

message("Presence / background: ",
        sum(model_df$presence == 1), " / ", sum(model_df$presence == 0))

write.csv(model_df, file.path(out_dir, "model_df_presence_background.csv"),
          row.names = FALSE)

# Distance-to-coast pattern (demo): Euclidean distance to study-region boundary.
# For Mediterranean work, replace with a real coastline / EMODnet distance raster.
study_sf <- st_as_sf(study_region) %>% st_make_valid()
pts_sf <- st_as_sf(
  model_df,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE
)
boundary <- tryCatch(
  st_cast(st_boundary(study_sf), "MULTILINESTRING"),
  error = function(e) st_cast(study_sf, "MULTILINESTRING")
)
model_df$dist_boundary_m <- as.numeric(st_distance(pts_sf, boundary))

# ---- 6. Spatial block CV (defensible evaluation) ----------------------------
if (has_pkg("blockCV")) {
  library(blockCV)
  # Environmental / spatial blocks on presence+background points
  sb <- cv_spatial(
    x = st_as_sf(
      model_df,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326
    ),
    column = "presence",
    k = 5,
    size = 400000, # metres — tweak for Med basins later
    selection = "random",
    iteration = 50,
    seed = 42,
    plot = FALSE
  )
  model_df$fold <- sb$folds_ids
  message("blockCV folds assigned: ", paste(sort(unique(model_df$fold)), collapse = ", "))
} else {
  # Fallback: random folds
  model_df$fold <- sample(1:5, nrow(model_df), replace = TRUE)
  message("blockCV not installed — using random 5-fold (not ideal for space).")
}

# ---- 7. Fit a simple, reusable model: GAM (stats-friendly) ------------------
# Interpretable niche: smooths of temperature, AOU, depth
fit_gam_species <- function(dat) {
  mgcv::gam(
    presence ~ mgcv::s(temperature, k = 5) + mgcv::s(aou, k = 5) + mgcv::s(depth, k = 5),
    data = dat,
    family = binomial(link = "logit"),
    method = "REML"
  )
}

# Leave-one-fold-out AUC
auc_bin <- function(obs, pred) {
  # Mann–Whitney / Wilcoxon form of AUC (no pROC dependency)
  pos <- pred[obs == 1]
  neg <- pred[obs == 0]
  if (!length(pos) || !length(neg)) return(NA_real_)
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

cv_scores <- lapply(sort(unique(model_df$fold)), function(k) {
  train <- model_df %>% filter(fold != k)
  test  <- model_df %>% filter(fold == k)
  m <- fit_gam_species(train)
  p <- predict(m, newdata = test, type = "response")
  data.frame(fold = k, auc = auc_bin(test$presence, p), n_test = nrow(test))
}) %>% bind_rows()

message("GAM spatial-block CV AUC:")
print(cv_scores)
print(summarise(cv_scores, mean_auc = mean(auc, na.rm = TRUE),
                sd_auc = sd(auc, na.rm = TRUE)))

gam_full <- fit_gam_species(model_df)
print(summary(gam_full))

png(file.path(out_dir, "03_gam_smooths.png"), width = 900, height = 350)
par(mfrow = c(1, 3))
plot(gam_full, pages = 0, scheme = 1)
dev.off()

# ---- 8. Optional: sdmTMB (spatial GLMM — great with survey PA data) ---------
if (has_pkg("sdmTMB")) {
  library(sdmTMB)
  message("Fitting sdmTMB (spatial random field)...")
  mesh <- make_mesh(
    model_df,
    xy_cols = c("decimalLongitude", "decimalLatitude"),
    cutoff = 1.5
  )
  fit_tmb <- try(
    sdmTMB(
      presence ~ s(temperature) + s(aou) + s(depth),
      data = model_df,
      mesh = mesh,
      family = binomial(link = "logit"),
      spatial = "on"
    ),
    silent = TRUE
  )
  if (!inherits(fit_tmb, "try-error")) {
    print(fit_tmb)
    saveRDS(fit_tmb, file.path(out_dir, "sdmTMB_fit.rds"))
  } else {
    message("sdmTMB fit failed on this toy grid — OK for first play; retry with denser Med data.")
  }
} else {
  message("Skip sdmTMB (not installed). install.packages('sdmTMB') when ready.")
}

# ---- 9. Project niche in 3D (depth stack) -----------------------------------
# Envelope-style projection for intuition (voluModel vignette style),
# then also a GAM probability stack on a coarse subset of depths.

temp_lims <- quantile(pres$temperature, c(0, 1), na.rm = TRUE)
aou_lims  <- quantile(pres$aou, c(0, 1), na.rm = TRUE)

temp_suit <- classify(
  temperature,
  rcl = matrix(
    c(-Inf, temp_lims[[1]], 0,
      temp_lims[[1]], temp_lims[[2]], 1,
      temp_lims[[2]], Inf, 0),
    ncol = 3, byrow = TRUE
  )
)
aou_suit <- classify(
  oxygen,
  rcl = matrix(
    c(-Inf, aou_lims[[1]], 0,
      aou_lims[[1]], aou_lims[[2]], 1,
      aou_lims[[2]], Inf, 0),
    ncol = 3, byrow = TRUE
  )
)

envelope_3d <- temp_suit * aou_suit
envelope_3d <- mask(crop(envelope_3d, study_region), study_region)
names(envelope_3d) <- names(temperature)

# Depths actually occupied by the species in this demo
idx_occ <- sort(unique(vapply(
  occurrences$depth,
  function(z) which.min(abs(layer_depths - z)),
  integer(1)
)))

png(file.path(out_dir, "04_envelope_3d_layers.png"), width = 800, height = 700)
plotLayers(
  envelope_3d[[min(idx_occ):max(idx_occ)]],
  title = "3D envelope suitability (red=shallow, blue=deep)",
  land = land
)
dev.off()

# GAM probability at selected depths (multi-species template idea)
project_gam_depth <- function(model, temp_stack, aou_stack, depths_m, study) {
  out <- list()
  for (d in depths_m) {
    i <- which.min(abs(as.numeric(names(temp_stack)) - d))
    nd <- c(
      temperature = temp_stack[[i]],
      aou = aou_stack[[i]]
    )
    # Predict from SpatRaster via as.data.frame of values
    coords <- as.data.frame(temp_stack[[i]], xy = TRUE, na.rm = FALSE)
    names(coords)[1:3] <- c("decimalLongitude", "decimalLatitude", "temperature")
    coords$aou <- values(aou_stack[[i]], mat = FALSE)
    coords$depth <- d
    coords$pred <- as.numeric(
      predict(model, newdata = coords, type = "response")
    )
    r <- temp_stack[[i]]
    values(r) <- coords$pred
    r <- mask(crop(r, study), study)
    names(r) <- paste0("p_depth_", d)
    out[[as.character(d)]] <- r
  }
  rast(out)
}

# Use a few representative depths for a fast demo
demo_depths <- sort(unique(round(quantile(pres$depth, probs = c(0.2, 0.5, 0.8)))))
gam_stack <- project_gam_depth(
  gam_full, temperature, oxygen, demo_depths, study_region
)

png(file.path(out_dir, "05_gam_prob_by_depth.png"), width = 1100, height = 400)
plot(gam_stack, main = paste("GAM P(presence) @", paste(demo_depths, "m")))
dev.off()

writeRaster(
  gam_stack,
  file.path(out_dir, "gam_prob_depth_stack.tif"),
  overwrite = TRUE
)

# ---- 10. Multi-species wrapper (the pattern you want for your thesis) -------
# One function, many species: only the occurrence table changes.
fit_3d_sdm <- function(species_id,
                       occ_df,
                       env_list,
                       n_background = 1000,
                       k_folds = 5) {
  stopifnot(all(c("decimalLongitude", "decimalLatitude", "depth") %in% names(occ_df)))
  stopifnot(all(c("temperature", "aou") %in% names(env_list)))

  occ_df <- occ_df %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(depth)) %>%
    distinct()

  # 3D extract at presences
  occ_df$temperature <- xyzSample(occ_df, env_list$temperature)
  occ_df$aou <- xyzSample(occ_df, env_list$aou)
  occ_df$presence <- 1L
  occ_df <- occ_df %>% filter(!is.na(temperature), !is.na(aou))

  region <- marineBackground(occ_df, buff = 1e6, clipToOcean = TRUE,
                             alpha = 1, partCount = 1)

  bg_xy <- spatSample(region, size = n_background, method = "random")
  if (inherits(bg_xy, "SpatVector")) bg_xy <- crds(bg_xy)
  bg <- as.data.frame(bg_xy)
  names(bg)[1:2] <- c("decimalLongitude", "decimalLatitude")
  bg$depth <- sample(occ_df$depth, nrow(bg), replace = TRUE)
  bg$temperature <- xyzSample(bg, env_list$temperature)
  bg$aou <- xyzSample(bg, env_list$aou)
  bg$presence <- 0L
  bg <- bg %>% filter(!is.na(temperature), !is.na(aou))

  dat <- bind_rows(
    occ_df[, c("decimalLongitude", "decimalLatitude", "depth",
               "temperature", "aou", "presence")],
    bg[, c("decimalLongitude", "decimalLatitude", "depth",
           "temperature", "aou", "presence")]
  )
  dat$fold <- sample(seq_len(k_folds), nrow(dat), replace = TRUE)
  dat$species <- species_id

  model <- mgcv::gam(
    presence ~ mgcv::s(temperature, k = 5) + mgcv::s(aou, k = 5) + mgcv::s(depth, k = 5),
    data = dat,
    family = binomial(),
    method = "REML"
  )

  list(species = species_id, data = dat, model = model, region = region)
}

# Fake a second "species" by jittering depths — shows the multi-species loop
env_list <- list(temperature = temperature, aou = oxygen)
sp1 <- fit_3d_sdm("Steindachneria_argentea", occurrences, env_list, n_background = 800)

occ_sp2 <- occurrences %>%
  mutate(depth = pmax(50, depth * runif(n(), 0.6, 1.1)))
sp2 <- fit_3d_sdm("Demo_species_B", occ_sp2, env_list, n_background = 800)

multi_summary <- data.frame(
  species = c(sp1$species, sp2$species),
  n_pres = c(sum(sp1$data$presence == 1), sum(sp2$data$presence == 1)),
  n_bg = c(sum(sp1$data$presence == 0), sum(sp2$data$presence == 0)),
  deviance_explained = c(
    summary(sp1$model)$dev.expl,
    summary(sp2$model)$dev.expl
  )
)
print(multi_summary)
write.csv(multi_summary, file.path(out_dir, "multi_species_summary.csv"),
          row.names = FALSE)

# ---- 11. Optional biomod2 sketch (ensemble for many species) ----------------
# biomod2 wants a presence / PA table + env as SpatRaster (2D). For 3D you
# typically flatten to "points-with-data" (what we built) OR fit per depth.
# Below is a minimal 2D surface demo so you can see the API; keep voluModel
# for depth-aware covariates in the real Med workflow.
if (has_pkg("biomod2")) {
  library(biomod2)
  message("biomod2 is installed. Minimal 2D surface ensemble sketch...")

  biomod_occ <- occurrences[, c("decimalLongitude", "decimalLatitude")]
  names(biomod_occ) <- c("x", "y")
  # Presence vector (all rows are presences)
  my_resp <- rep(1, nrow(biomod_occ))

  my_biomod_data <- BIOMOD_FormatingData(
    resp.var = my_resp,
    expl.var = c(temperature[[1]], oxygen[[1]]),
    resp.xy = biomod_occ,
    resp.name = "LuminousHake",
    PA.nb.rep = 1,
    PA.nb.absences = min(500, n_bg),
    PA.strategy = "random"
  )

  my_biomod_model <- try(
    BIOMOD_Modeling(
      my_biomod_data,
      modeling.id = "demo",
      models = c("GLM", "GAM", "RF"),
      CV.strategy = "random",
      CV.nb.rep = 2,
      CV.perc = 0.8,
      metric.eval = c("TSS", "ROC"),
      var.import = 1
    ),
    silent = TRUE
  )

  if (!inherits(my_biomod_model, "try-error")) {
    print(my_biomod_model)
    saveRDS(my_biomod_model, file.path(out_dir, "biomod2_demo.rds"))
  } else {
    message("biomod2 modeling sketch failed — API versions vary; use help(BIOMOD_Modeling).")
  }
} else {
  message("Skip biomod2 (not installed). install.packages('biomod2') for ensembles.")
}

# ---- 12. Clean up + next steps ----------------------------------------------
unlink(td, recursive = TRUE)

message("\n=== Demo finished ===")
message("Artifacts written to: ", normalizePath(out_dir))
message("
Next for Mediterranean thesis work:
  1. Replace occurrences with MEDITS / OBIS (lon, lat, depth, species, year).
  2. Replace temperature/oxygen stacks with CMEMS 3D NetCDF via terra::rast().
  3. Add predictors: salinity, O2, chl, bathymetry, distance-to-coast.
  4. Loop fit_3d_sdm() / sdmTMB over species with the SAME env + CV design.
  5. Project historical vs future CMEMS temperature for climate sensitivity.
")

