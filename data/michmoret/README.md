# Michmoret local database

Drop your Michmoret / Ruppin fish table here as CSV.

Required columns (names can vary; the script remaps them):

- `scientificName` (or `species`)
- `decimalLongitude` / `decimalLatitude` (or `lon` / `lat`)
- `depth` in metres (or `depthInMeters`)
- optional: `eventDate`, `occurrenceStatus` (`present` / `absent`), `individualCount`

If the file has true survey absences (`occurrenceStatus = absent`), those will be used.
Otherwise the pipeline treats rows as presences only.
