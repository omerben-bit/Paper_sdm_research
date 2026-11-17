# MEDITS (Mediterranean International Trawl Survey) - summary and README

Written by Claude on 2026-10-24 from the files in this folder, including the MEDITS Handbook (version 9, 2017). The collection method below is taken from that handbook.

## Files
| File | What it is |
|---|---|
| `hauls_data.csv` | One row per haul: 13,528 hauls x 46 columns. Date, position, depth, duration, gear geometry, bottom temperature and salinity. |
| `medits_raw_data.csv` | The catch table: 448,881 rows x 22 columns. One row per haul x taxon (and size/sex category). Encoding is Latin-1, so read it with `encoding="latin-1"`. |
| `taxa_medits.csv` | Species code key: 1,617 MEDITS codes with scientific name, authority, fauna category (`CATFAU`) and the GSAs where each is recorded. Also Latin-1. |
| `Medits_Handbook_2017_version_9_5-60417r.pdf` | The survey instruction manual (protocol, gear, strata). |

## What the dataset is
Bottom-trawl survey data covering the northern Mediterranean shelf and slope, 1994-2020: 13,528 hauls (13,526 MEDITS and 2 SOLEMON), 8 countries (ESP, FRA, ITA, SVN, HRV, GRC, MLT, CYP). Only 13 of the Mediterranean's GSAs (Geographical Sub-Areas) are present: 1, 2, 5, 6, 7, 8, 15, 17, 18, 20, 22, 23, 25. GSAs 16 and 19 are not in these files. Years per GSA range from 12 (GSA 2) to 27.

## Structure and how to join the tables
- **Haul key:** `country`, `area` (the GSA), `vessel`, `year`, `haul_number`. This is unique in `hauls_data.csv`. The `id` column differs between the two files and does not join them.
- **Join check:** 13,404 of the 13,435 haul keys in the catch file find a match in `hauls_data.csv`; 31 do not (probably an issue with one of the key fields, which I did not resolve).
- **Catch table:** `genus` (4 letters) + `species` (3 letters) make the MEDITS species code (compare with `MEDITS Code` in the taxa key). `nbtot` is the number of individuals, `ptot` the total weight in g, `nbfem`/`nbmal`/`nbind` the sex split. `catfau` is the fauna category (`Ao` is the largest group, 230,003 rows). `lref` shows which size-measurement reference was used.
- **Haul table:** `shooting_depth` / `hauling_depth` (m; median 119, range 0-877), `hauling_duration` (min), `distance` and `wing_opening` (for swept area), `bottom_temperature_beginning/end` (4% missing), `bottom_salinity_*` (value -1 means missing), `validity` (V = valid, 13,441; I = invalid, 87), `number_of_the_stratum`. Months are mostly May-July (11,391 hauls).

## Method used to collect the data (from the MEDITS Handbook)
- **Period:** centred on June (May-July), kept consistent across years.
- **Sampling design:** depth-stratified; positions are drawn at random within each stratum, with the number of hauls proportional to stratum area. Strata limits are fixed at 10-50, 51-100, 101-200, 201-500 and 501-800 m. In practice hauls are repeated in the same positions from year to year.
- **Gear:** a GOC73 bottom trawl with a cod-end mesh of 10 mm mesh side (about 20 mm opening). Speed is 3 knots.
- **Timing:** daylight only (from 30 min after sunrise to 30 min before sunset). Haul duration is 30 min at depths under 200 m and 60 min deeper. A haul is valid if at least 2/3 of the time or distance was achieved.
- **Catch processing:** the catch is identified to species, counted and weighed (`nbtot`, `ptot`). Sex, length, individual weight and maturity are collected for different species groups (MEDITS Group 1 and Group 2 lists; see section 3 of the handbook). The raw catch table keeps the counts by sex (`nbfem`, `nbmal`, `nbind`).

## Using it for the course
- **Abundance:** counts and weight per species per haul. The sampling effort differs between hauls, so your earlier subset kept hauls with swept area within +-10% of the median and did not divide counts by swept area, to keep integer counts.
- **Treatments:** GSA/region (west, central, east), depth stratum, year.
- **Continuous variables:** depth, bottom temperature, salinity.
- **Known issue:** temperature is almost completely confounded with area (from your earlier test runs: GSA explains about 36% of composition, temperature adds 0.7% after GSA and depth).
- **Cleaning needed:** drop `validity` = I hauls; `hydrological_station` is 55% missing and untidy; `observations` is a numeric code; use Latin-1 for the catch and taxa tables. Species names in `taxa_medits.csv` have damaged characters (such as "R?ppell") in the authority column.
