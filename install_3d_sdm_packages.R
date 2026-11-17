# Lean Windows install for the 3D SDM demo
# Prefer CRAN binaries so Rtools does NOT have to compile everything.

options(timeout = 600)
options(pkgType = "binary")          # critical on Windows
options(install.packages.check.source = "no")

# Mirror that usually resolves; change if needed
options(repos = c(CRAN = "https://cloud.r-project.org"))

core <- c(
  "terra", "sf", "dplyr", "tidyr", "ggplot2", "tidyterra",
  "mgcv", "voluModel", "CoordinateCleaner", "blockCV"
)

# Install binaries only; dependencies = TRUE but binary-first
install.packages(core, dependencies = TRUE, type = "binary")

# Report what actually loaded
needed <- core
ok <- vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
print(data.frame(package = needed, installed = unname(ok)))

if (!all(ok)) {
  message("\nStill missing: ", paste(needed[!ok], collapse = ", "))
  message("Check internet / VPN, then re-run this script.")
} else {
  message("\nAll core packages OK. Next:\n  source('3d_sdm_demo.R', echo = TRUE)")
}
