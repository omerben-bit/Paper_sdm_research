# Paper SDM research — Mediterranean fish

R project for 3D/2D species distribution work: Israeli/Levantine surveys, paper-cleaned GBIF, maps, and dashboards.

## Layout

| Path | Role |
|------|------|
| `01_med_fish_3d_pa_maps.R` … `08_observations_dashboard.R` | Main pipelines (occurrence merge → maps → summary → dashboard) |
| `gbif_paper.R` | Load canonical cleaned GBIF (`data/processed/gbif_fish_levant_cleaned.csv`) |
| `data/processed/` | Exports (combined occurrences, PA table, Singer subsets) |
| `primary_sdm/` | Focus-species map outputs |
| `dashboard/` | HTML observations dashboard (`index.html`) |
| `dataset_summary/` | Inventory markdown/PDF from `06_dataset_summary.R` |

## Run

From this folder in RStudio or R:

```r
source("01_med_fish_3d_pa_maps.R")   # refresh GBIF + combined CSV
source("08_observations_dashboard.R")
```

## Large data (local only)

Git does not track some files (see `.gitignore`): Singer/Owens iDigBio CSV (~113 MB), `MEDITS/medits_raw_data.csv`, and `data/raw/*.rds` / `*.zip`. Keep these on disk and regenerate with the scripts where noted.

## Remote

- **GitHub:** [omerben-bit/Paper_sdm_research](https://github.com/omerben-bit/Paper_sdm_research)
- **Branch:** `main` — SDM code; `template` on the remote is the original phd-kb scaffold.

Set **default branch** to `main` on GitHub (Settings → General → Default branch).