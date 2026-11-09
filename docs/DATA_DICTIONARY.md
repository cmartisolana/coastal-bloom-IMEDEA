# Data dictionary

All NPZ files are readable with `numpy.load(path, allow_pickle=False)`. Missing floating values are NaN. String times are ISO months unless stated otherwise. Geographic coordinates are degrees east/north (western longitudes negative).

| File/variable | Meaning and units |
|---|---|
| `data/processed/*_ring_monthly.npz`: `months` | 216 months, January 2003–December 2020 |
| `ring_area_km2`, `ring_limit_km` | Ocean area and outer distance edge for 20 coastal rings |
| `occupancy_all`, `occupancy_cyclonic`, `occupancy_anticyclonic` | Monthly mean daily covered-area fractions, shape month × ring, 0–1 |
| `occurrence_map_*.npz`: `lon`, `lat`, `freq_pct`, `dist_km`, `n_days` | Map axes, latitude × longitude eddy frequency (%), mean fine-grid distance (km), number of daily samples |
| `*_ekman_rings.npz`: `ring_ekman_xshelf_15m`, `ring_total_xshelf_15m` | Monthly ring offshore velocity, m/s (multiply by 100 for cm/s) |
| `ring_ekman_div_15m` | Monthly ring divergence proxy, s⁻¹ |
| `map_mean_*` | Mean monthly current/divergence map, m/s or s⁻¹ |
| `cov_*` | Fraction of ring ocean area represented by temporally supported current cells |
| `rings_*.npz`: `dist`, `lon`, `lat`, `res` | Frozen nominal 0.02° geometry, distance in km, NaN outside domain; `res` in degrees |
| `ne10m_*` | Prepared coast vertices, groups and land polygons; NaN separates polygon rings |
| `data/intermediate/meta4_*`: `time` | Days since 1950-01-01; centre and contour coordinates in degrees |
| `data/intermediate/mob_*` | Monthly current-grid inputs; velocity m/s, divergence s⁻¹; `ue/ve` Ekman and `uo/vo` total components |

`results/tables/all_trends.csv` is the canonical statistical table:

- `region`, `metric`: unique identifiers. `eddy_all`, `eddy_cyclonic`, `eddy_anticyclonic`; velocity component/depth identifiers; `ring_a-b` is distance a–b km.
- `estimate`, `standard_error`, `ci95_low`, `ci95_high`, `units`: OLS estimate and matching HAC inference, scaled as stated in units.
- `p_value`: two-sided asymptotic normal test of zero slope; unadjusted for multiple comparisons.
- `n`, `hac_lags`, `model`: sample size, monthly/annual lag count and model specification.
- `descriptive_sen`: Theil–Sen slope of STL trend for centroid rows only, km/year. Blank for other metrics. It does not share the OLS confidence interval.
- `significant_unadjusted`: CI excludes zero.

CSV names containing `occupancy_monthly` hold fractions; `centroid_monthly` hold km and the descriptive STL series; `monthly_cm_s` hold cm/s; `divergence_monthly_s-1` hold s⁻¹. Ring columns such as `0-10` denote distance in km. Annual share CSVs hold percentages, not fractions; their 20 (or 17) ring columns sum to 100 each year. `nearshore_0_60`, where present, is an additional summary column and must not be included again in that sum.
