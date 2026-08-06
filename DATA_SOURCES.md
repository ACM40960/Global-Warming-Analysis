# Data Sources & Provenance
This file documents exactly where each raw file came from

| File | Dataset name / version | Source organization | URL | Accessed |
|---|---|---|---|---|
| `Land_Temperature.asc` | e.g. HadCRUT5 / GISTEMP land-only annual anomaly | e.g. Met Office Hadley Centre / NASA GISS | _fill in_ | _fill in_ |
| `Global_Temperature.asc` | e.g. HadCRUT5 or GISTEMP blended land+ocean | _fill in_ | _fill in_ | _fill in_ |
| `Sea_Temperature.nc` | ERSSTv5 gridded monthly SST | NOAA NCEI | _fill in_ | _fill in_ |
| `Solar_Irradiance.nc` | NOAA/LASP composite TSI reconstruction | NOAA / LASP | _fill in_ | _fill in_ |
| `Global_C02_Annual.csv` | NOAA GML global annual mean CO2 | NOAA Global Monitoring Laboratory | _fill in_ | _fill in_ |
| `Volcanic_Activity.tsv` | Global Volcanism Program eruption catalog | Smithsonian Institution | _fill in_ | _fill in_ |

## Notes on data curation choices

The Project raises two hard curation questions directly: how many
individual locations are "enough" to represent the globe (e.g. 100 points
gives only ~2 for the whole of Europe), and how to handle wildly
different record lengths per station (centuries of data for London or
Paris vs. a few decades for a station like Vostok). This project
sidesteps hand-picking individual weather stations by using two
already-gridded/pre-aggregated products instead:

- **Ocean coverage**: `Sea_Temperature.nc` is a full 89×180 lat/lon grid
  (up to ~16,000 cells), so ocean coverage is treated as complete and
  latitude-area-weighted rather than sampled from a handful of buoys
  this directly avoids the "only 2 points for Europe" undersampling
  problem for the ocean domain.
- **Land coverage**: `Land_Temperature.asc` is a pre-built global land
  average (not a set of individual station records), so the
  station-selection and record-length-mismatch problem (e.g. London vs.
  Vostok) was already solved upstream.
  This project reuses that curation rather than re-deriving it
  from raw station-level data.
