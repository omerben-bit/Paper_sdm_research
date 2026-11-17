# Make the key teaching figures from the voluModel demo data.
# Run after packages are installed: source("plot_demo_figures.R")

library(voluModel)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)

out_dir <- file.path("outputs_3d_sdm_demo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- data ----
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", package = "voluModel")) %>%
  dplyr::select(decimalLongitude, decimalLatitude, depth) %>%
  dplyr::filter(!is.na(depth), depth >= 1, depth <= 2000) %>%
  dplyr::distinct()

td <- tempdir()
unzip(
  system.file("extdata/woa18_decav_t00mn01_cropped.zip", package = "voluModel"),
  exdir = file.path(td, "temperature"), junkpaths = TRUE
)
temperature_pts <- vect(file.path(td, "temperature", "woa18_decav_t00mn01_cropped.shp"))
temperature_bottom <- bottomRaster(temperature_pts)
template <- centerPointRasterTemplate(temperature_pts)
temperature <- terra::rasterize(temperature_pts, template, field = names(temperature_pts))
envt_names <- gsub("[d,M]", "", names(temperature_pts))
envt_names[[1]] <- "0"
names(temperature) <- envt_names

oxygen <- rast(system.file("extdata/oxygenSmooth.tif", package = "voluModel"))
names(oxygen) <- names(temperature)

# downsample to voxels
layer_depths <- as.numeric(names(temperature))
occs$index <- vapply(occs$depth, function(z) which.min(abs(layer_depths - z)), 1L)
occs <- bind_rows(lapply(sort(unique(occs$index)), function(i) {
  pts <- occs[occs$index == i, ]
  pts <- downsample(pts, temperature[[1]], verbose = FALSE)
  pts$depth <- layer_depths[[i]]
  pts
})) %>% dplyr::select(decimalLongitude, decimalLatitude, depth)

# ---- Fig 1: occurrence map ----
occ_df <- as.data.frame(occs)
p1 <- ggplot(occ_df, aes(decimalLongitude, decimalLatitude, color = depth)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_viridis_c(option = "C", direction = -1, name = "Depth (m)") +
  coord_fixed(1.1) +
  labs(
    title = "Where the fish were found (example species)",
    subtitle = "Luminous hake — each point has lon, lat, AND depth",
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "fig1_occurrence_map.png"), p1, width = 8, height = 5, dpi = 150)

# ---- Fig 2: surface vs bottom temperature maps ----
surf <- as.data.frame(temperature[[1]], xy = TRUE, na.rm = TRUE)
names(surf) <- c("x", "y", "temp")
bot <- as.data.frame(temperature_bottom, xy = TRUE, na.rm = TRUE)
names(bot) <- c("x", "y", "temp")
surf$which <- "Surface temperature"
bot$which <- "Bottom temperature"
tb <- bind_rows(surf, bot)

p2 <- ggplot(tb, aes(x, y, fill = temp)) +
  geom_raster() +
  scale_fill_viridis_c(option = "H", name = "°C") +
  facet_wrap(~which) +
  coord_fixed(1.1) +
  labs(
    title = "Same sea, very different temperatures",
    subtitle = "Surface vs seafloor — why depth matters for fish SDMs",
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(size = 8))
ggsave(file.path(out_dir, "fig2_surface_vs_bottom_temp.png"), p2, width = 10, height = 4.5, dpi = 150)

# ---- Fig 3: THE key teaching plot — extract method comparison ----
vals <- bind_rows(
  data.frame(
    temperature = terra::extract(temperature[[1]], occs[, 1:2])[, 2],
    method = "1. Surface only (2D)"
  ),
  data.frame(
    temperature = terra::extract(temperature_bottom, occs[, 1:2])[, 2],
    method = "2. Bottom only (2D)"
  ),
  data.frame(
    temperature = xyzSample(occs, temperature),
    method = "3. At true depth (3D)"
  )
) %>%
  dplyr::filter(is.finite(temperature))

p3 <- ggplot(vals, aes(method, temperature, fill = method)) +
  geom_boxplot(outlier.alpha = 0.35, width = 0.65) +
  labs(
    title = "The 3D SDM idea in one plot",
    subtitle = "Same fish locations — temperature depends on HOW you extract it",
    x = NULL, y = "Temperature at occurrence (°C)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10))
ggsave(file.path(out_dir, "fig3_why_3d_extraction.png"), p3, width = 8, height = 4.5, dpi = 150)

# ---- Fig 4: depth profile of occurrences ----
p4 <- ggplot(occ_df, aes(depth)) +
  geom_histogram(bins = 20, fill = "#2a6f97", color = "white") +
  labs(
    title = "Vertical habitat of the example species",
    subtitle = "3D SDMs treat depth as a real axis, not just a map colour",
    x = "Depth (m)", y = "Number of occurrence voxels"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "fig4_depth_histogram.png"), p4, width = 7, height = 4, dpi = 150)

# ---- Fig 5: simple GAM response curves ----
dat <- occs %>%
  mutate(
    temperature = xyzSample(occs, temperature),
    aou = xyzSample(occs, oxygen),
    presence = 1L
  )
dat <- dat[stats::complete.cases(dat), ]

set.seed(1)
region <- marineBackground(occs, buff = 1e6, clipToOcean = TRUE, alpha = 1, partCount = 1)
bg_xy <- spatSample(region, size = min(1500, 5 * nrow(dat)), method = "random")
if (inherits(bg_xy, "SpatVector")) bg_xy <- crds(bg_xy)
bg <- as.data.frame(bg_xy)
names(bg)[1:2] <- c("decimalLongitude", "decimalLatitude")
bg$depth <- sample(dat$depth, nrow(bg), replace = TRUE)
bg$temperature <- xyzSample(bg, temperature)
bg$aou <- xyzSample(bg, oxygen)
bg$presence <- 0L
bg <- bg[stats::complete.cases(bg), ]

model_df <- bind_rows(
  dat[, c("decimalLongitude", "decimalLatitude", "depth", "temperature", "aou", "presence")],
  bg[, c("decimalLongitude", "decimalLatitude", "depth", "temperature", "aou", "presence")]
)

# Use mgcv:: explicitly — another package's s()/gam() may be masking these
m <- mgcv::gam(
  presence ~ mgcv::s(temperature, k = 5) + mgcv::s(aou, k = 5) + mgcv::s(depth, k = 5),
  data = model_df, family = binomial(), method = "REML"
)

# predict smooths on a grid for ggplot
pred_smooth <- function(var, m, data) {
  rng <- range(data[[var]], na.rm = TRUE)
  new <- data.frame(
    temperature = mean(data$temperature),
    aou = mean(data$aou),
    depth = mean(data$depth)
  )
  new <- new[rep(1, 100), ]
  new[[var]] <- seq(rng[1], rng[2], length.out = 100)
  new$fit <- predict(m, newdata = new, type = "response")
  new$variable <- var
  new$x <- new[[var]]
  new[, c("variable", "x", "fit")]
}

curve_df <- bind_rows(
  pred_smooth("temperature", m, model_df),
  pred_smooth("aou", m, model_df),
  pred_smooth("depth", m, model_df)
)
curve_df$variable <- factor(
  curve_df$variable,
  levels = c("temperature", "aou", "depth"),
  labels = c("Temperature (°C)", "AOU (oxygen proxy)", "Depth (m)")
)

p5 <- ggplot(curve_df, aes(x, fit)) +
  geom_line(linewidth = 1.1, color = "#1b4332") +
  facet_wrap(~variable, scales = "free_x") +
  labs(
    title = "What the model learned (GAM response curves)",
    subtitle = "Probability of presence vs each predictor (others held at mean)",
    x = NULL, y = "P(presence)"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "fig5_gam_response_curves.png"), p5, width = 10, height = 3.8, dpi = 150)

# ---- Fig 6: suitability at a few depths (mini 3D stack) ----
# envelope from 3D-extracted limits
pres <- model_df %>% filter(presence == 1)
tlim <- quantile(pres$temperature, c(0, 1), na.rm = TRUE)
alim <- quantile(pres$aou, c(0, 1), na.rm = TRUE)

demo_depths <- as.numeric(names(temperature))[
  unique(vapply(quantile(pres$depth, c(0.25, 0.5, 0.75)), function(z) {
    which.min(abs(as.numeric(names(temperature)) - z))
  }, 1L))
]

stack_df <- lapply(demo_depths, function(d) {
  i <- which.min(abs(as.numeric(names(temperature)) - d))
  tlayer <- temperature[[i]]
  alayer <- oxygen[[i]]
  suit <- (tlayer >= tlim[1] & tlayer <= tlim[2]) *
    (alayer >= alim[1] & alayer <= alim[2])
  suit <- mask(crop(suit, region), region)
  df <- as.data.frame(suit, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "suitable"
  df$depth_m <- paste0(d, " m")
  df
}) %>% bind_rows()

p6 <- ggplot(stack_df, aes(x, y, fill = factor(suitable, labels = c("no", "yes")))) +
  geom_raster() +
  scale_fill_manual(values = c("no" = "#e9ecef", "yes" = "#1d3557"), name = "Suitable") +
  facet_wrap(~depth_m) +
  coord_fixed(1.1) +
  labs(
    title = "3D projection: suitable habitat by depth layer",
    subtitle = "Same niche rules, mapped at different depths (envelope model)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "fig6_suitability_by_depth.png"), p6, width = 10, height = 4.2, dpi = 150)

unlink(td, recursive = TRUE)

message("\nFigures written to: ", normalizePath(out_dir))
message(paste(" -", list.files(out_dir, pattern = "^fig.*\\.png$"), collapse = "\n"))
