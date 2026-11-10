# Data dictionary

Regions are `benguela` and `gulf_stream`. "Ring" means one of the twenty 10 km distance-from-coast bands (0–10, 10–20, …, 190–200 km); in CSV columns a ring is written as `0-10`, `10-20`, etc. Months run from January 2003 to December 2020 (216 months, 18 years). Missing values are NaN. Longitudes are degrees east (west negative), latitudes degrees north. NPZ files open with `numpy.load(path, allow_pickle=False)`.

## `results/trends/` — trend estimates

`all_trends.csv` has one row per region and metric (260 rows); `key_results.csv` is the subset with the main results.

| Column | Meaning |
|---|---|
| `region`, `metric` | What was fitted (metric names below) |
| `estimate` | OLS trend, in `units` |
| `ci95_low`, `ci95_high` | 95% confidence interval (HAC standard errors) |
| `standard_error`, `p_value` | HAC standard error; two-sided P value, not adjusted for multiple comparisons |
| `n`, `hac_lags`, `model` | Number of observations, HAC lags, regression model |
| `descriptive_sen` | Theil–Sen slope of the STL trend (centroid rows only, km/year); shown as a dashed line in the eddy figure, not covered by the CI |
| `significant_unadjusted` | True when the CI excludes zero |

| Metric | Meaning | Units |
|---|---|---|
| `eddy_<type>\|centroid` | Cross-shore centroid of eddy occupancy; `<type>` is `all`, `cyclonic` or `anticyclonic` | km/year |
| `eddy_<type>\|share_0_60` | Share of eddy occupancy within 0–60 km | percentage points/decade |
| `eddy_<type>\|share_ring_<a-b>` | Share of eddy occupancy in one ring | percentage points/decade |
| `ekman_xshelf_15m\|ring_<a-b>`, `total_xshelf_15m\|ring_<a-b>` | Offshore Ekman / total (geostrophic + Ekman) velocity at 15 m in one ring | cm/s/decade |
| `ekman_xshelf_15m\|band_<a_b>`, `total_xshelf_15m\|band_<a_b>` | Same, area-weighted over the band 0–60, 60–200 or 30–60 km | cm/s/decade |
| `ekman_div_15m\|share_ring_<a-b>` | Share of positive Ekman divergence (upwelling proxy) in one ring, 30–200 km | percentage points/decade |
| `ekman_div_15m\|upwelling_share_30_60` | Same, summed over 30–60 km | percentage points/decade |

## `results/series/` — series the trends are fitted to

One file per region (`<region>_…`); one column per ring plus `month` or `year`.

| File | Content | Units |
|---|---|---|
| `<region>_eddy_<type>_occupancy_monthly.csv` | Fraction of each ring's ocean area inside eddies (monthly mean of daily values) | 0–1 |
| `<region>_eddy_<type>_centroid_monthly.csv` | Centroid (`centroid_km`) and its STL trend (`stl_trend_km`) | km |
| `<region>_eddy_<type>_share_annual.csv` | Each ring's share of the year's eddy occupancy (rings sum to 100); `nearshore_0_60` is the 0–60 km sum, not part of that total | % |
| `<region>_ekman_offshore_velocity_monthly.csv`, `<region>_total_offshore_velocity_monthly.csv` | Offshore velocity at 15 m, positive seaward | cm/s |
| `<region>_ekman_divergence_monthly.csv` | Ekman divergence at 15 m | s⁻¹ |
| `<region>_upwelling_share_annual.csv` | Each ring's share of positive divergence, 30–200 km (rings sum to 100) | % |

## `data/processed/` — figure inputs

| File: variables | Content | Units |
|---|---|---|
| `<region>_coastal_rings.npz`: `dist`, `lon`, `lat`, `res`, `area_km2` | 0.02° ocean grid with distance from the coast; ring areas | km, degrees, km² |
| `<region>_eddy_occupancy.npz`: `months`, `ring_area_km2`, `ring_limit_km`, `occupancy_all`, `occupancy_cyclonic`, `occupancy_anticyclonic` | Monthly eddy occupancy per ring (month × ring) | 0–1 |
| `<region>_eddy_occurrence_map.npz`: `lon`, `lat`, `freq_pct`, `dist_km`, `n_days` | Percentage of days inside an eddy, on a 0.04° grid | % |
| `<region>_currents.npz`: `ring_ekman_xshelf_15m`, `ring_total_xshelf_15m`, `ring_ekman_div_15m` | Monthly ring means of offshore velocity and divergence (month × ring) | m/s, s⁻¹ |
| `<region>_currents.npz`: `map_mean_*`, `cov_*`, `lon`, `lat` | Mean maps on the 0.25° grid; fraction of each ring covered by current cells | m/s, s⁻¹, 0–1 |
| `ne10m_coast_filtered_2km.npz`, `ne10m_land_coast_lines.npz` | Coastline vertices used for distances; land polygons for maps | degrees |

## `data/inputs/` — regional extracts of the source products

| File | Content |
|---|---|
| `meta4_<region>_<cyclonic\|anticyclonic>.npz` | META4.0 eddies in the regional box (±3°): `time` (days since 1950-01-01), centre `lon`/`lat`, 30-point effective contours |
| `globcurrent_benguela_chunk0-3.npz` | Monthly GlobCurrent offshore velocity and divergence grids at 15 m, in four time chunks |
| `globcurrent_gulf_stream_monthly.npz` | Monthly GlobCurrent velocity components at 15 m: `ue`/`ve` (Ekman), `uo`/`vo` (total), m/s |
