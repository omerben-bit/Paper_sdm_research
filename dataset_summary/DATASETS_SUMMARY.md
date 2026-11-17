# Dataset summary (Master's project)

Generated: 2026-11-17 10:53 IST

## Definitions

- **Sightings**: individual catch/specimen **records** (fish rows).
- **Observations**: unique **hauls/stations** where `obs_key` exists;
  for Singer/GBIF each sighting is also one observation.
- **Depth bands**: `50 m` bins (0–50, 50–100, …, 500+), plus `depth_unknown`.
- **Observation depths**: median depth per haul (`obs_key`).

## Sightings and observations (part 1 — sources + processed)

Survey/raw layers (analysis scope) plus pipeline exports under `data/processed/`.
Joint totals below use **survey / raw source** rows only (no double-counting).

| Layer | Dataset | Scope | Sightings | Observations | Species | % depth | Source |
|-------|---------|-------|----------:|-------------:|--------:|--------:|--------|
| Survey / raw source | Israel_trawl_MIT | Israel coast box (lon 32–36.5°, lat 29.5–34.5°) | 88,662 | 213 | 176 | 100% | `Combined_Menachem_Itai_corrected_May_2026_DBG.csv` |
| Survey / raw source | Edelist_OBIS | Actinopterygii + Israel box | 0 | 0 | 0 | — | `data/raw/edelist_occurrence.txt` |
| Survey / raw source | MEDITS_raw | GSA 22 & 25 + Levant box (maps script 03) | 9,106 | 388 | 233 | 100% | `MEDITS/medits_raw_data.csv (+ hauls_data.csv)` |
| Survey / raw source | Singer_Owens_iDigBio | Israel box, depth-resolved specimens | 24 | 24 | 23 | 100% | `data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv` |
| Survey / raw source | GBIF_Levant_cache | Levant GBIF download (combined cache) | 802 | 802 | 100 | 46.3% | `data/processed/med_fish_occurrences_combined.csv` |
| Processed (`data/processed/`) | proc_med_fish_occurrences_combined | All rows in `01_med_fish_3d_pa_maps.R` combined export | 802 | 802 | 100 | 46.3% | `data/processed/med_fish_occurrences_combined.csv` |
| Processed (`data/processed/`) | proc_pa_focus_species | Presence / absence / pseudo-absence table for first 3D map species | 356 | 356 | 1 | 96.6% | `data/processed/pa_focus_species.csv` |
| Processed (`data/processed/`) | proc_singer_owens_israel_box | Singer & Owens iDigBio — Israel coast box, depth required | 23 | 23 | 23 | 100% | `data/processed/singer_owens_israel_box.csv` |
| Processed (`data/processed/`) | proc_singer_med_focus5 | Singer & Owens — Mediterranean clip, FOCUS_FIVE species only | 12 | 12 | 4 | 0% | `data/processed/singer_med_focus5.csv` |

## Raw file sizes (unfiltered)

| File | Rows |
|------|-----:|
| `Combined_Menachem_Itai_corrected_May_2026_DBG.csv` | 90,132 |
| `data/raw/edelist_occurrence.txt` | 13,941 |
| `MEDITS/medits_raw_data.csv` | 448,881 |
| `MEDITS/hauls_data.csv` | 13,528 |
| `data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv` | 354,622 |

## Processed files (`data/processed/`) — detail

Same exports as **Processed** rows in [part 1](#sightings-and-observations-part-1--sources--processed).
File row counts and breakdowns below.

| File | Rows on disk |
|------|-------------:|
| `data/processed/med_fish_occurrences_combined.csv` | 802 |
| `data/processed/pa_focus_species.csv` | 356 |
| `data/processed/singer_owens_israel_box.csv` | 23 |
| `data/processed/singer_med_focus5.csv` | 12 |

### `med_fish_occurrences_combined.csv` — rows by `source`

| source | occurrenceStatus | rows |
|--------|------------------|-----:|
| GBIF | present | 802 |

### Processed files — spatial extent

| Dataset | N with coords | Lon min–max | Lat min–max |
|---------|--------------:|-------------|-------------|
| proc_med_fish_occurrences_combined | 802 | 33.424 – 35.458 | 31.375 – 33.791 |
| proc_pa_focus_species | 356 | 32.001 – 35.488 | 31.004 – 33.972 |
| proc_singer_owens_israel_box | 23 | 34.780 – 35.489 | 29.500 – 33.901 |
| proc_singer_med_focus5 | 12 | 7.417 – 34.780 | 32.500 – 43.567 |

CSVs: `processed_files_summary.csv`, `processed_file_row_counts.csv`,
`processed_location_extent.csv`, `processed_depth_distribution.csv`.


## Spatial extent (records with coordinates)

| Dataset | N | Lon min–max | Lat min–max | Centroid | 0.5° grid cells |
|---------|----:|-------------|-------------|----------|----------------:|
| Israel_trawl_MIT | 88,662 | 34.224 – 34.481 | 31.402 – 31.781 | 34.480, 31.780 | 2 |
| MEDITS_raw | 9,106 | 32.339 – 34.088 | 34.533 – 35.108 | 33.181, 34.755 | 7 |
| Singer_Owens_iDigBio | 24 | 34.780 – 35.489 | 29.500 – 33.901 | 34.865, 31.569 | 3 |
| GBIF_Levant_cache | 802 | 33.424 – 35.458 | 31.375 – 33.791 | 34.817, 32.407 | 15 |

## Columns by file

Full list: `columns_by_file.csv` (one row per column).

### Israel_trawl_MIT (`Combined_Menachem_Itai_corrected_May_2026_DBG.csv`)

...1, SAMPLE_DESC, Date_obj, SESSION_ID, LATITUDE, LONGITUDE, DEPTH, Day.Night, ID_SPECIES, ABUNDANCE, BIOMASS, DAY, MONTH, YEAR, StartTime, EndTime, HaulDuration_min

### Edelist_OBIS (`data/raw/edelist_occurrence.txt`)

id, institutionCode, collectionCode, basisOfRecord, occurrenceID, catalogNumber, occurrenceStatus, eventID, fieldNumber, scientificNameID, scientificName, kingdom, phylum, class, order, family, genus, specificEpithet, scientificNameAuthorship, nomenclaturalCode, taxonRemarks

### MEDITS_catch (`MEDITS/medits_raw_data.csv`)

id, country, area, vessel, year, haul_number, codend_closing, partit, catfau, genus, species, lref, ptot, nbtot, nbfem, nbmal, nbind, tf, month, day, upload_id, name_of_survey

### MEDITS_hauls (`MEDITS/hauls_data.csv`)

id, country, area, vessel, gear, rigging, door, year, month, day, haul_number, codend_closing, shooting_time, shooting_quadrant, shooting_latitude, shooting_longitude, shooting_depth, hauling_time, hauling_quadrant, hauling_latitude, hauling_longitude, hauling_depth, hauling_duration, validity, course, recorded_species, distance, vertical_opening, wing_opening, geometrical_precision, bridles_length, warp_length, warp_diameter, hydrological_station, observations, type_of_file, bottom_temperature_beginning, bottom_temperature_end, measuring_system_temp, number_of_the_stratum, part_of_the_codend, upload_id, name_of_survey, bottom_salinity_beginning, bottom_salinity_end, measuring_system_sal

### MEDITS_taxa (`MEDITS/taxa_medits.csv`)

MEDITS    Code, Scientific Name  valid, Authorship, Source, Reference, Remarks, CATFAU, CODLON, GSAs, Year

### Singer_Owens_iDigBio (`data/raw/iDigBioProcessedFishOccurrences_SingerAndOwens.csv`)

SOcurated_scientificName, iDigBio_coreid, SOcurated_family, SOcurated_genus, iDigBio_dwc_verbatimDepth, SOcurated_maximumDepthInMeters, SOcurated_minimumDepthInMeters, iDigBio_dwc_decimalLongitude, iDigBio_dwc_decimalLatitude, iDigBio_dwc_year, iDigBio_dwc_month, iDigBio_dwc_day, iDigBio_dwc_basisOfRecord, iDigBio_dwc_samplingProtocol, iDigBio_dwc_coordinateUncertaintyInMeters, iDigBio_dwc_coordinatePrecision, iDigBio_dwc_institutionCode, iDigBio_dwc_collectionID, iDigBio_dwc_collectionCode, iDigBio_dwc_datasetName, iDigBio_dwc_catalogNumber, SOcurated_bibliographicCitation, SOcurated_license, SOcurated_depthUncertainty, SOcurated_bathymetry, SOcurated_depth, WoRMS_AphiaID, WoRMS_scientificName, WoRMS_family

### GBIF_Levant_cache (raw layer) (`data/processed/med_fish_occurrences_combined.csv`)

scientificName, decimalLongitude, decimalLatitude, depth, occurrenceStatus, eventDate, source, depth_bin

### proc_med_fish_occurrences_combined (`data/processed/med_fish_occurrences_combined.csv`)

scientificName, decimalLongitude, decimalLatitude, depth, occurrenceStatus, eventDate, source, depth_bin

### proc_pa_focus_species (`data/processed/pa_focus_species.csv`)

scientificName, decimalLongitude, decimalLatitude, depth, occurrenceStatus, eventDate, source, depth_bin, pa, pa_type, station_id

### proc_singer_owens_israel_box (`data/processed/singer_owens_israel_box.csv`)

scientificName, decimalLatitude, decimalLongitude, depth, dataset, year

### proc_singer_med_focus5 (`data/processed/singer_med_focus5.csv`)

scientificName, decimalLatitude, decimalLongitude, depth, dataset, year

## Temporal coverage (year x month)

- Sightings: `month_year_sightings.csv`
- Observations (hauls): `month_year_observations.csv`
- **Joint (all DBs):** `joint_observations_year_month.csv`,
  `joint_observations_by_month.csv`, `joint_observations_by_year.csv`,
  `joint_observations_by_depth.csv`, `joint_summary_totals.csv`
- Per-database report: `DATASETS_SUMMARY.pdf`
- **All databases combined:** `DATASETS_COMBINED.pdf`
- Joint sightings CSVs: `joint_sightings_*.csv`

### All databases combined — observation counts

- **1,427** joint observations · **98,594** fish records on those units

| Year | Month | Observations |
|-----:|------:|-------------:|
| 1926 | 1 (Jan) | 2 |
| 1934 | 3 (Mar) | 1 |
| 1936 | 1 (Jan) | 1 |
| 1936 | 3 (Mar) | 2 |
| 1953 | 5 (May) | 1 |
| 1954 | 3 (Mar) | 1 |
| 1954 | 5 (May) | 1 |
| 1954 | 9 (Sep) | 1 |
| 1955 | 5 (May) | 1 |
| 1958 | 2 (Feb) | 1 |
| 1958 | 11 (Nov) | 3 |
| 1959 | 2 (Feb) | 1 |
| 1959 | 3 (Mar) | 1 |
| 1959 | 5 (May) | 1 |
| 1959 | 8 (Aug) | 2 |
| 1961 | 11 (Nov) | 1 |
| 1962 | 7 (Jul) | 9 |
| 1963 | 5 (May) | 1 |
| 1963 | 8 (Aug) | 1 |
| 1964 | 2 (Feb) | 1 |
| 1964 | 7 (Jul) | 1 |
| 1965 | 1 (Jan) | 1 |
| 1965 | 4 (Apr) | 2 |
| 1965 | 9 (Sep) | 1 |
| 1966 | 7 (Jul) | 1 |
| 1968 | 2 (Feb) | 1 |
| 1969 | 8 (Aug) | 18 |
| 1971 | 5 (May) | 1 |
| 1972 | 6 (Jun) | 1 |
| 1972 | 12 (Dec) | 1 |
| 1977 | 3 (Mar) | 1 |
| 1977 | 4 (Apr) | 3 |
| 1977 | 5 (May) | 2 |
| 1979 | 5 (May) | 1 |
| 1979 | 7 (Jul) | 1 |
| 1979 | 9 (Sep) | 2 |
| 1980 | 5 (May) | 1 |
| 1980 | 7 (Jul) | 2 |
| 1982 | 6 (Jun) | 1 |
| 1982 | 9 (Sep) | 1 |
| 1984 | 3 (Mar) | 5 |
| 1984 | 9 (Sep) | 2 |
| 1985 | 6 (Jun) | 7 |
| 1985 | 9 (Sep) | 1 |
| 1986 | 11 (Nov) | 4 |
| 1987 | 1 (Jan) | 1 |
| 1987 | 2 (Feb) | 2 |
| 1987 | 5 (May) | 1 |
| 1987 | 6 (Jun) | 1 |
| 1988 | 3 (Mar) | 1 |
| 1988 | 12 (Dec) | 1 |
| 1989 | 3 (Mar) | 1 |
| 1989 | 5 (May) | 1 |
| 1990 | 7 (Jul) | 1 |
| 1990 | 12 (Dec) | 1 |
| 1991 | 1 (Jan) | 3 |
| 1991 | 2 (Feb) | 1 |
| 1991 | 3 (Mar) | 1 |
| 1991 | 4 (Apr) | 2 |
| 1991 | 5 (May) | 9 |
| 1991 | 6 (Jun) | 1 |
| 1991 | 7 (Jul) | 3 |
| 1991 | 10 (Oct) | 2 |
| 1992 | 1 (Jan) | 2 |
| 1992 | 2 (Feb) | 1 |
| 1992 | 4 (Apr) | 2 |
| 1992 | 5 (May) | 3 |
| 1992 | 6 (Jun) | 2 |
| 1992 | 7 (Jul) | 4 |
| 1992 | 8 (Aug) | 1 |
| 1992 | 9 (Sep) | 1 |
| 1992 | 11 (Nov) | 1 |
| 1992 | 12 (Dec) | 2 |
| 1993 | 1 (Jan) | 1 |
| 1993 | 2 (Feb) | 3 |
| 1993 | 5 (May) | 1 |
| 1993 | 6 (Jun) | 2 |
| 1993 | 8 (Aug) | 1 |
| 1993 | 10 (Oct) | 1 |
| 1993 | 12 (Dec) | 2 |
| 1994 | 2 (Feb) | 2 |
| 1994 | 3 (Mar) | 2 |
| 1994 | 4 (Apr) | 5 |
| 1994 | 6 (Jun) | 1 |
| 1994 | 8 (Aug) | 11 |
| 1994 | 10 (Oct) | 1 |
| 1995 | 3 (Mar) | 2 |
| 1996 | 12 (Dec) | 1 |
| 1998 | 1 (Jan) | 1 |
| 1999 | 1 (Jan) | 1 |
| 2000 | 2 (Feb) | 6 |
| 2000 | 6 (Jun) | 2 |
| 2000 | 7 (Jul) | 4 |
| 2000 | 8 (Aug) | 2 |
| 2000 | 11 (Nov) | 1 |
| 2001 | 6 (Jun) | 4 |
| 2004 | 5 (May) | 1 |
| 2005 | 8 (Aug) | 25 |
| 2005 | 12 (Dec) | 1 |
| 2006 | 3 (Mar) | 1 |
| 2006 | 6 (Jun) | 22 |
| 2006 | 7 (Jul) | 3 |
| 2007 | 6 (Jun) | 25 |
| 2008 | 4 (Apr) | 6 |
| 2008 | 5 (May) | 8 |
| 2008 | 6 (Jun) | 27 |
| 2008 | 7 (Jul) | 8 |
| 2008 | 8 (Aug) | 8 |
| 2008 | 11 (Nov) | 10 |
| 2008 | 12 (Dec) | 14 |
| 2009 | 1 (Jan) | 3 |
| 2009 | 2 (Feb) | 14 |
| 2009 | 3 (Mar) | 2 |
| 2009 | 4 (Apr) | 9 |
| 2009 | 5 (May) | 3 |
| 2009 | 6 (Jun) | 23 |
| 2009 | 7 (Jul) | 17 |
| 2009 | 8 (Aug) | 5 |
| 2009 | 9 (Sep) | 1 |
| 2009 | 11 (Nov) | 8 |
| 2009 | 12 (Dec) | 1 |
| 2010 | 1 (Jan) | 1 |
| 2010 | 2 (Feb) | 11 |
| 2010 | 3 (Mar) | 5 |
| 2010 | 4 (Apr) | 11 |
| 2010 | 5 (May) | 14 |
| 2010 | 6 (Jun) | 31 |
| 2010 | 7 (Jul) | 13 |
| 2010 | 8 (Aug) | 13 |
| 2010 | 9 (Sep) | 7 |
| 2010 | 10 (Oct) | 3 |
| 2010 | 11 (Nov) | 34 |
| 2010 | 12 (Dec) | 3 |
| 2011 | 1 (Jan) | 5 |
| 2011 | 2 (Feb) | 4 |
| 2011 | 3 (Mar) | 11 |
| 2011 | 4 (Apr) | 7 |
| 2011 | 5 (May) | 26 |
| 2011 | 6 (Jun) | 3 |
| 2011 | 7 (Jul) | 34 |
| 2011 | 8 (Aug) | 8 |
| 2011 | 9 (Sep) | 2 |
| 2011 | 10 (Oct) | 4 |
| 2011 | 11 (Nov) | 4 |
| 2011 | 12 (Dec) | 5 |
| 2012 | 1 (Jan) | 11 |
| 2012 | 2 (Feb) | 18 |
| 2012 | 3 (Mar) | 11 |
| 2012 | 4 (Apr) | 2 |
| 2012 | 5 (May) | 8 |
| 2012 | 6 (Jun) | 59 |
| 2012 | 7 (Jul) | 1 |
| 2012 | 8 (Aug) | 8 |
| 2012 | 10 (Oct) | 2 |
| 2012 | 11 (Nov) | 12 |
| 2013 | 4 (Apr) | 8 |
| 2013 | 6 (Jun) | 25 |
| 2013 | 7 (Jul) | 6 |
| 2013 | 8 (Aug) | 7 |
| 2014 | 5 (May) | 1 |
| 2014 | 6 (Jun) | 6 |
| 2014 | 8 (Aug) | 2 |
| 2014 | 12 (Dec) | 6 |
| 2015 | 6 (Jun) | 6 |
| 2015 | 8 (Aug) | 26 |
| 2015 | 12 (Dec) | 1 |
| 2016 | 6 (Jun) | 26 |
| 2016 | 8 (Aug) | 4 |
| 2016 | 11 (Nov) | 3 |
| 2016 | 12 (Dec) | 1 |
| 2017 | 4 (Apr) | 1 |
| 2017 | 6 (Jun) | 30 |
| 2017 | 10 (Oct) | 3 |
| 2017 | 11 (Nov) | 4 |
| 2017 | 12 (Dec) | 1 |
| 2018 | 1 (Jan) | 3 |
| 2018 | 3 (Mar) | 14 |
| 2018 | 4 (Apr) | 21 |
| 2018 | 5 (May) | 2 |
| 2018 | 7 (Jul) | 6 |
| 2018 | 8 (Aug) | 22 |
| 2018 | 10 (Oct) | 10 |
| 2018 | 11 (Nov) | 3 |
| 2018 | 12 (Dec) | 1 |
| 2019 | 1 (Jan) | 1 |
| 2019 | 2 (Feb) | 6 |
| 2019 | 5 (May) | 2 |
| 2019 | 6 (Jun) | 30 |
| 2019 | 8 (Aug) | 12 |
| 2019 | 10 (Oct) | 1 |
| 2019 | 11 (Nov) | 8 |
| 2019 | 12 (Dec) | 1 |
| 2020 | 2 (Feb) | 1 |
| 2020 | 4 (Apr) | 4 |
| 2020 | 5 (May) | 5 |
| 2020 | 6 (Jun) | 2 |
| 2020 | 7 (Jul) | 26 |
| 2020 | 8 (Aug) | 1 |
| 2020 | 9 (Sep) | 1 |
| 2020 | 10 (Oct) | 3 |
| 2020 | 11 (Nov) | 1 |
| 2020 | 12 (Dec) | 1 |
| 2021 | 1 (Jan) | 7 |
| 2021 | 4 (Apr) | 2 |
| 2021 | 5 (May) | 2 |
| 2021 | 6 (Jun) | 1 |
| 2021 | 9 (Sep) | 2 |
| 2021 | 10 (Oct) | 1 |
| 2021 | 11 (Nov) | 11 |
| 2022 | 1 (Jan) | 1 |
| 2022 | 2 (Feb) | 1 |
| 2022 | 4 (Apr) | 6 |
| 2022 | 5 (May) | 1 |
| 2022 | 6 (Jun) | 3 |
| 2022 | 7 (Jul) | 2 |
| 2022 | 9 (Sep) | 2 |
| 2022 | 10 (Oct) | 3 |
| 2022 | 11 (Nov) | 3 |
| 2022 | 12 (Dec) | 2 |
| 2023 | 3 (Mar) | 1 |
| 2023 | 9 (Sep) | 2 |
| 2023 | 10 (Oct) | 2 |
| 2023 | 11 (Nov) | 1 |
| 2023 | 12 (Dec) | 1 |
| 2024 | 1 (Jan) | 2 |
| 2024 | 5 (May) | 8 |
| 2024 | 6 (Jun) | 5 |
| 2024 | 7 (Jul) | 1 |
| 2024 | 9 (Sep) | 2 |
| 2024 | 10 (Oct) | 5 |
| 2024 | 12 (Dec) | 2 |
| 2025 | 1 (Jan) | 3 |
| 2025 | 2 (Feb) | 1 |
| 2025 | 3 (Mar) | 8 |
| 2025 | 4 (Apr) | 2 |
| 2025 | 5 (May) | 6 |
| 2025 | 6 (Jun) | 3 |
| 2025 | 9 (Sep) | 19 |
| 2025 | 10 (Oct) | 8 |
| 2025 | 11 (Nov) | 13 |
| 2025 | 12 (Dec) | 2 |
| 2026 | 2 (Feb) | 2 |
| 2026 | 3 (Mar) | 5 |
| 2026 | 5 (May) | 2 |
| 2026 | 6 (Jun) | 22 |
| 2026 | 7 (Jul) | 1 |
| 2026 | 8 (Aug) | 4 |
| 2026 | 9 (Sep) | 5 |
| 2026 | 10 (Oct) | 5 |

### Israel_trawl_MIT — sightings by year/month

| Year | Month | Sightings |
|-----:|------:|----------:|
| 2008 | 4 (Apr) | 1,961 |
| 2008 | 5 (May) | 3,861 |
| 2008 | 7 (Jul) | 5,097 |
| 2008 | 8 (Aug) | 5,196 |
| 2008 | 11 (Nov) | 2,731 |
| 2008 | 12 (Dec) | 5,835 |
| 2009 | 2 (Feb) | 3,993 |
| 2009 | 4 (Apr) | 3,795 |
| 2009 | 7 (Jul) | 3,562 |
| 2009 | 11 (Nov) | 4,135 |
| 2010 | 2 (Feb) | 3,970 |
| 2010 | 4 (Apr) | 3,597 |
| 2010 | 7 (Jul) | 3,724 |
| 2010 | 8 (Aug) | 3,972 |
| 2010 | 11 (Nov) | 10,547 |
| 2011 | 3 (Mar) | 4,212 |
| 2011 | 5 (May) | 7,505 |
| 2011 | 7 (Jul) | 843 |
| 2011 | 8 (Aug) | 474 |
| 2012 | 1 (Jan) | 3,378 |
| 2012 | 2 (Feb) | 2,955 |
| 2012 | 3 (Mar) | 3,032 |
| 2012 | 8 (Aug) | 131 |
| 2013 | 4 (Apr) | 156 |

### MEDITS_raw — sightings by year/month

| Year | Month | Sightings |
|-----:|------:|----------:|
| 2005 | 8 (Aug) | 574 |
| 2006 | 6 (Jun) | 507 |
| 2006 | 7 (Jul) | 74 |
| 2007 | 6 (Jun) | 594 |
| 2008 | 6 (Jun) | 528 |
| 2009 | 6 (Jun) | 454 |
| 2009 | 7 (Jul) | 185 |
| 2010 | 6 (Jun) | 567 |
| 2011 | 7 (Jul) | 618 |
| 2012 | 6 (Jun) | 625 |
| 2013 | 6 (Jun) | 598 |
| 2015 | 8 (Aug) | 650 |
| 2016 | 6 (Jun) | 626 |
| 2017 | 6 (Jun) | 660 |
| 2018 | 7 (Jul) | 127 |
| 2018 | 8 (Aug) | 497 |
| 2019 | 6 (Jun) | 591 |
| 2020 | 7 (Jul) | 631 |

### Singer_Owens_iDigBio — sightings by year/month

| Year | Month | Sightings |
|-----:|------:|----------:|
| 1958 | 11 (Nov) | 1 |
| 1969 | 8 (Aug) | 15 |
| 1972 | 6 (Jun) | 1 |
| 1977 | 4 (Apr) | 3 |
| 1977 | 5 (May) | 1 |
| 1986 | 11 (Nov) | 2 |
| 1987 | 5 (May) | 1 |

### GBIF_Levant_cache — sightings by year/month

| Year | Month | Sightings |
|-----:|------:|----------:|
| 1926 | 1 (Jan) | 2 |
| 1934 | 3 (Mar) | 1 |
| 1936 | 1 (Jan) | 1 |
| 1936 | 3 (Mar) | 2 |
| 1953 | 5 (May) | 1 |
| 1954 | 3 (Mar) | 1 |
| 1954 | 5 (May) | 1 |
| 1954 | 9 (Sep) | 1 |
| 1955 | 5 (May) | 1 |
| 1958 | 2 (Feb) | 1 |
| 1958 | 11 (Nov) | 2 |
| 1959 | 2 (Feb) | 1 |
| 1959 | 3 (Mar) | 1 |
| 1959 | 5 (May) | 1 |
| 1959 | 8 (Aug) | 2 |
| 1961 | 11 (Nov) | 1 |
| 1962 | 7 (Jul) | 9 |
| 1963 | 5 (May) | 1 |
| 1963 | 8 (Aug) | 1 |
| 1964 | 2 (Feb) | 1 |
| 1964 | 7 (Jul) | 1 |
| 1965 | 1 (Jan) | 1 |
| 1965 | 4 (Apr) | 2 |
| 1965 | 9 (Sep) | 1 |
| 1966 | 7 (Jul) | 1 |
| 1968 | 2 (Feb) | 1 |
| 1969 | 8 (Aug) | 3 |
| 1971 | 5 (May) | 1 |
| 1972 | 12 (Dec) | 1 |
| 1977 | 3 (Mar) | 1 |
| 1977 | 5 (May) | 1 |
| 1979 | 5 (May) | 1 |
| 1979 | 7 (Jul) | 1 |
| 1979 | 9 (Sep) | 2 |
| 1980 | 5 (May) | 1 |
| 1980 | 7 (Jul) | 2 |
| 1982 | 6 (Jun) | 1 |
| 1982 | 9 (Sep) | 1 |
| 1984 | 3 (Mar) | 5 |
| 1984 | 9 (Sep) | 2 |
| 1985 | 6 (Jun) | 7 |
| 1985 | 9 (Sep) | 1 |
| 1986 | 11 (Nov) | 2 |
| 1987 | 1 (Jan) | 1 |
| 1987 | 2 (Feb) | 2 |
| 1987 | 6 (Jun) | 1 |
| 1988 | 3 (Mar) | 1 |
| 1988 | 12 (Dec) | 1 |
| 1989 | 3 (Mar) | 1 |
| 1989 | 5 (May) | 1 |
| 1990 | 7 (Jul) | 1 |
| 1990 | 12 (Dec) | 1 |
| 1991 | 1 (Jan) | 3 |
| 1991 | 2 (Feb) | 1 |
| 1991 | 3 (Mar) | 1 |
| 1991 | 4 (Apr) | 2 |
| 1991 | 5 (May) | 9 |
| 1991 | 6 (Jun) | 1 |
| 1991 | 7 (Jul) | 3 |
| 1991 | 10 (Oct) | 2 |
| 1992 | 1 (Jan) | 2 |
| 1992 | 2 (Feb) | 1 |
| 1992 | 4 (Apr) | 2 |
| 1992 | 5 (May) | 3 |
| 1992 | 6 (Jun) | 2 |
| 1992 | 7 (Jul) | 4 |
| 1992 | 8 (Aug) | 1 |
| 1992 | 9 (Sep) | 1 |
| 1992 | 11 (Nov) | 1 |
| 1992 | 12 (Dec) | 2 |
| 1993 | 1 (Jan) | 1 |
| 1993 | 2 (Feb) | 3 |
| 1993 | 5 (May) | 1 |
| 1993 | 6 (Jun) | 2 |
| 1993 | 8 (Aug) | 1 |
| 1993 | 10 (Oct) | 1 |
| 1993 | 12 (Dec) | 2 |
| 1994 | 2 (Feb) | 2 |
| 1994 | 3 (Mar) | 2 |
| 1994 | 4 (Apr) | 5 |
| 1994 | 6 (Jun) | 1 |
| 1994 | 8 (Aug) | 11 |
| 1994 | 10 (Oct) | 1 |
| 1995 | 3 (Mar) | 2 |
| 1996 | 12 (Dec) | 1 |
| 1998 | 1 (Jan) | 1 |
| 1999 | 1 (Jan) | 1 |
| 2000 | 2 (Feb) | 6 |
| 2000 | 6 (Jun) | 2 |
| 2000 | 7 (Jul) | 4 |
| 2000 | 8 (Aug) | 2 |
| 2000 | 11 (Nov) | 1 |
| 2001 | 6 (Jun) | 4 |
| 2004 | 5 (May) | 1 |
| 2005 | 12 (Dec) | 1 |
| 2006 | 3 (Mar) | 1 |
| 2008 | 11 (Nov) | 2 |
| 2008 | 12 (Dec) | 6 |
| 2009 | 1 (Jan) | 3 |
| 2009 | 2 (Feb) | 6 |
| 2009 | 3 (Mar) | 2 |
| 2009 | 4 (Apr) | 1 |
| 2009 | 5 (May) | 3 |
| 2009 | 6 (Jun) | 3 |
| 2009 | 7 (Jul) | 2 |
| 2009 | 8 (Aug) | 5 |
| 2009 | 9 (Sep) | 1 |
| 2009 | 12 (Dec) | 1 |
| 2010 | 1 (Jan) | 1 |
| 2010 | 2 (Feb) | 3 |
| 2010 | 3 (Mar) | 5 |
| 2010 | 4 (Apr) | 3 |
| 2010 | 5 (May) | 14 |
| 2010 | 6 (Jun) | 4 |
| 2010 | 7 (Jul) | 5 |
| 2010 | 8 (Aug) | 5 |
| 2010 | 9 (Sep) | 7 |
| 2010 | 10 (Oct) | 3 |
| 2010 | 11 (Nov) | 11 |
| 2010 | 12 (Dec) | 3 |
| 2011 | 1 (Jan) | 5 |
| 2011 | 2 (Feb) | 4 |
| 2011 | 3 (Mar) | 3 |
| 2011 | 4 (Apr) | 7 |
| 2011 | 5 (May) | 2 |
| 2011 | 6 (Jun) | 3 |
| 2011 | 7 (Jul) | 4 |
| 2011 | 8 (Aug) | 3 |
| 2011 | 9 (Sep) | 2 |
| 2011 | 10 (Oct) | 4 |
| 2011 | 11 (Nov) | 4 |
| 2011 | 12 (Dec) | 5 |
| 2012 | 1 (Jan) | 3 |
| 2012 | 2 (Feb) | 10 |
| 2012 | 3 (Mar) | 3 |
| 2012 | 4 (Apr) | 2 |
| 2012 | 5 (May) | 8 |
| 2012 | 6 (Jun) | 33 |
| 2012 | 7 (Jul) | 1 |
| 2012 | 10 (Oct) | 2 |
| 2012 | 11 (Nov) | 12 |
| 2013 | 4 (Apr) | 1 |
| 2013 | 7 (Jul) | 6 |
| 2013 | 8 (Aug) | 7 |
| 2014 | 5 (May) | 1 |
| 2014 | 6 (Jun) | 6 |
| 2014 | 8 (Aug) | 2 |
| 2014 | 12 (Dec) | 6 |
| 2015 | 6 (Jun) | 6 |
| 2015 | 12 (Dec) | 1 |
| 2016 | 8 (Aug) | 4 |
| 2016 | 11 (Nov) | 3 |
| 2016 | 12 (Dec) | 1 |
| 2017 | 4 (Apr) | 1 |
| 2017 | 6 (Jun) | 4 |
| 2017 | 10 (Oct) | 3 |
| 2017 | 11 (Nov) | 4 |
| 2017 | 12 (Dec) | 1 |
| 2018 | 1 (Jan) | 3 |
| 2018 | 3 (Mar) | 14 |
| 2018 | 4 (Apr) | 21 |
| 2018 | 5 (May) | 2 |
| 2018 | 7 (Jul) | 1 |
| 2018 | 8 (Aug) | 1 |
| 2018 | 10 (Oct) | 10 |
| 2018 | 11 (Nov) | 3 |
| 2018 | 12 (Dec) | 1 |
| 2019 | 1 (Jan) | 1 |
| 2019 | 2 (Feb) | 6 |
| 2019 | 5 (May) | 2 |
| 2019 | 6 (Jun) | 5 |
| 2019 | 8 (Aug) | 12 |
| 2019 | 10 (Oct) | 1 |
| 2019 | 11 (Nov) | 8 |
| 2019 | 12 (Dec) | 1 |
| 2020 | 2 (Feb) | 1 |
| 2020 | 4 (Apr) | 4 |
| 2020 | 5 (May) | 5 |
| 2020 | 6 (Jun) | 2 |
| 2020 | 8 (Aug) | 1 |
| 2020 | 9 (Sep) | 1 |
| 2020 | 10 (Oct) | 3 |
| 2020 | 11 (Nov) | 1 |
| 2020 | 12 (Dec) | 1 |
| 2021 | 1 (Jan) | 7 |
| 2021 | 4 (Apr) | 2 |
| 2021 | 5 (May) | 2 |
| 2021 | 6 (Jun) | 1 |
| 2021 | 9 (Sep) | 2 |
| 2021 | 10 (Oct) | 1 |
| 2021 | 11 (Nov) | 11 |
| 2022 | 1 (Jan) | 1 |
| 2022 | 2 (Feb) | 1 |
| 2022 | 4 (Apr) | 6 |
| 2022 | 5 (May) | 1 |
| 2022 | 6 (Jun) | 3 |
| 2022 | 7 (Jul) | 2 |
| 2022 | 9 (Sep) | 2 |
| 2022 | 10 (Oct) | 3 |
| 2022 | 11 (Nov) | 3 |
| 2022 | 12 (Dec) | 2 |
| 2023 | 3 (Mar) | 1 |
| 2023 | 9 (Sep) | 2 |
| 2023 | 10 (Oct) | 2 |
| 2023 | 11 (Nov) | 1 |
| 2023 | 12 (Dec) | 1 |
| 2024 | 1 (Jan) | 2 |
| 2024 | 5 (May) | 8 |
| 2024 | 6 (Jun) | 5 |
| 2024 | 7 (Jul) | 1 |
| 2024 | 9 (Sep) | 2 |
| 2024 | 10 (Oct) | 5 |
| 2024 | 12 (Dec) | 2 |
| 2025 | 1 (Jan) | 3 |
| 2025 | 2 (Feb) | 1 |
| 2025 | 3 (Mar) | 8 |
| 2025 | 4 (Apr) | 2 |
| 2025 | 5 (May) | 6 |
| 2025 | 6 (Jun) | 3 |
| 2025 | 9 (Sep) | 19 |
| 2025 | 10 (Oct) | 8 |
| 2025 | 11 (Nov) | 13 |
| 2025 | 12 (Dec) | 2 |
| 2026 | 2 (Feb) | 2 |
| 2026 | 3 (Mar) | 5 |
| 2026 | 5 (May) | 2 |
| 2026 | 6 (Jun) | 22 |
| 2026 | 7 (Jul) | 1 |
| 2026 | 8 (Aug) | 4 |
| 2026 | 9 (Sep) | 5 |
| 2026 | 10 (Oct) | 5 |

## Depth distributions

- Sightings: `depth_distribution_sightings.csv`
- Observations (haul median depth): `depth_distribution_observations.csv`

### Observation counts by depth band (pooled datasets in analysis scope)

**Israel_trawl_MIT**

| Depth band | Observations |
|------------|-------------:|
| 0–50 m | 128 |
| 100–150 m | 34 |
| 250–300 m | 4 |
| 400–450 m | 5 |
| 50–100 m | 42 |

**MEDITS_raw**

| Depth band | Observations |
|------------|-------------:|
| 0–50 m | 69 |
| 100–150 m | 49 |
| 150–200 m | 25 |
| 200–250 m | 1 |
| 250–300 m | 20 |
| 300–350 m | 23 |
| 350–400 m | 1 |
| 500+ m | 63 |
| 50–100 m | 137 |

**Singer_Owens_iDigBio**

| Depth band | Observations |
|------------|-------------:|
| 0–50 m | 9 |
| 50–100 m | 15 |

**GBIF_Levant_cache**

| Depth band | Observations |
|------------|-------------:|
| 0–50 m | 218 |
| 100–150 m | 13 |
| 150–200 m | 3 |
| 200–250 m | 18 |
| 250–300 m | 9 |
| 300–350 m | 3 |
| 350–400 m | 5 |
| 400–450 m | 9 |
| 500+ m | 37 |
| 50–100 m | 56 |
| depth_unknown | 431 |

## Notes

- **MEDITS** in maps uses GSA 22 & 25 inside the Levant box;
  the raw catch file covers the full MEDITS programme.
- **GBIF** cache is Levant-focused and mostly lacks depth (`pct_depth_known` low).
- **Part 1** lists processed exports alongside survey sources;
  `GBIF_Levant_cache` = GBIF-present slice of `med_fish_occurrences_combined.csv`;
  `proc_med_fish_occurrences_combined` = all rows in that file.
- Regenerate: `source("06_dataset_summary.R")` (markdown, CSV, PDF).

