# Moria — local fish datasets

Drop CSV / TSV files here (lon, lat, depth, species columns).

The mapping script `03_species_depth_maps_datasets.R` also loads:

- `Combined_Menachem_Itai_corrected_May_2026_DBG.csv` (Israel trawl)
- `MEDITS/` (hauls + catch + taxa)
- `data/raw/edelist_occurrence.txt`
- GBIF cache in `data/processed/`

Outputs: `primary_sdm/`
