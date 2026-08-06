# Data Sources & Provenance
This file documents exactly where each raw file came from

| File | Source organization | URL |
|---|---|---|
| `Land_Temperature.asc` | NOAA National Centers for Environmental Information | (https://www.ncei.noaa.gov/products/land-based-station/noaa-global-temp)|
| `Global_Temperature.asc` | NOAA National Centers for Environmental Information | (https://www.ncei.noaa.gov/products/land-based-station/noaa-global-temp) |
| `Sea_Temperature.nc` |  NOAA Physical Sciences Laboratory |(https://www.psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html) |
| `Solar_Irradiance.nc` |  NOAA National Centers for Environmental Information | (https://www.ncei.noaa.gov/products/climate-data-records/total-solar-irradiance) |
| `Global_C02_Annual.csv` |  NOAA Global Monitoring Laboratory |(https://gml.noaa.gov/ccgg/trends/gl_data.html) | 
| `Volcanic_Activity.tsv` |  NOAA National Centers for Environmental Information | (https://www.ngdc.noaa.gov/hazel/view/hazards/volcano/event-data?maxYear=2025&minYear=1850) |

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
