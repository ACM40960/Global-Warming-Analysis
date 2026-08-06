# Data Sources & Provenance

The assignment specifically asks for global temperature, CO2, volcanic,
and solar data pulled from recognized scientific archives (IPCC, NOAA, or
equivalent). This table documents exactly where each raw file came from,
so the sourcing is auditable rather than asserted.

**Fill in the blanks below with the actual dataset name, URL, and the date
you downloaded it** — this is the single most important addition for
academic-integrity purposes; a grader should be able to verify every file
traces back to a real, citable archive.

| File | Dataset name / version | Source organization | URL | Accessed |
|---|---|---|---|---|
| `Land_Temperature.asc` | e.g. HadCRUT5 / GISTEMP land-only annual anomaly | e.g. Met Office Hadley Centre / NASA GISS | _fill in_ | _fill in_ |
| `Global_Temperature.asc` | e.g. HadCRUT5 or GISTEMP blended land+ocean | _fill in_ | _fill in_ | _fill in_ |
| `Sea_Temperature.nc` | ERSSTv5 gridded monthly SST | NOAA NCEI | _fill in_ | _fill in_ |
| `Solar_Irradiance.nc` | NOAA/LASP composite TSI reconstruction | NOAA / LASP | _fill in_ | _fill in_ |
| `Global_C02_Annual.csv` | NOAA GML global annual mean CO2 | NOAA Global Monitoring Laboratory | _fill in_ | _fill in_ |
| `Volcanic_Activity.tsv` | Global Volcanism Program eruption catalog | Smithsonian Institution | _fill in_ | _fill in_ |

## Notes on data curation choices

The assignment raises two hard curation questions directly: how many
individual locations are "enough" to represent the globe (e.g. 100 points
gives only ~2 for the whole of Europe), and how to handle wildly
different record lengths per station (centuries of data for London or
Paris vs. a few decades for a station like Vostok). This project
sidesteps hand-picking individual weather stations by using two
already-gridded/pre-aggregated products instead:

- **Ocean coverage**: `Sea_Temperature.nc` is a full 89×180 lat/lon grid
  (up to ~16,000 cells), so ocean coverage is treated as complete and
  latitude-area-weighted rather than sampled from a handful of buoys —
  this directly avoids the "only 2 points for Europe" undersampling
  problem for the ocean domain.
- **Land coverage**: `Land_Temperature.asc` is a pre-built global land
  average (not a set of individual station records), so the
  station-selection and record-length-mismatch problem (e.g. London vs.
  Vostok) was already solved upstream by the agency that produced it —
  the same curation work the assignment says "the IPCC scientists did
  already." This project reuses that curation rather than re-deriving it
  from raw station-level data, which is a reasonable simplification for
  a course project but is worth stating explicitly rather than leaving
  implicit.
- **Trade-off**: this means the project does not itself demonstrate
  station-level curation (merging heterogeneous record lengths, infilling
  gaps, weighting sparse regions). If your assignment specifically wants
  you to show that curation process by hand, that would be a good
  "further work" addition — e.g. pulling 15-20 individual long-record
  station series from `ncei.noaa.gov/cdo-web` directly, and comparing
  their individually-fitted trends against the pre-blended product used
  here.
