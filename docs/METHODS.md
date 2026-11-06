# Datasets, processing and statistical methods

## Scope and spatial domain

Analyses cover 216 months (January 2003–December 2020), 18 annual observations and two rectangular regions. Benguela (`SW_Africa`): 9.83–19.77° E, 36.56–13.27° S. Gulf Stream (`NA_east`): 80.26–60.86° W, 31.03–46.09° N. Twenty distance bands span 0–200 km in 10 km increments. The outer edge (10, 20, …, 200 km) represents each band in centroid calculations. Geographic boxes are stored in `config/regions.json`.

The coast is a previously prepared Natural Earth 1:10 million coastline, filtered to long coastal segments (source workspace reports a 1428 km minimum) and resampled at about 2 km. The supplied coastline and land arrays are frozen inputs. A nominal 0.02° ocean grid is formed by rasterising Natural Earth land polygons. Great-circle distance to the filtered coastline determines the band; land and distances beyond 200 km are masked. The historical distance grid is stored as float16. Ocean-pixel area is approximated by `(R Δθ)^2 cos(latitude)`, with R = 6371.0088 km, when aggregating the archived grid.

These are approximate Table S1 box domains, not the manuscript authors' manually delineated coastal masks. In particular, Benguela includes the Agulhas Bank and the Gulf Stream box includes parts of the Gulf of Maine and Nova Scotia. Geographic agreement with the bloom analysis must not be described as exact. Nor is 10 km band spacing the effective observational resolution.

## Eddy occurrence

META4.0 delayed-time all-satellite effective eddy contours are used for both polarities. The original files span 1993-01-01–2023-09-08. Extraction selected eddy centres within each box padded by 3° and dates 2002-09-01–2020-12-31; only 2003–2020 is analysed here. Each contour has 30 vertices. No lifetime or trajectory filter is applied to these retained analyses.

Daily contour unions are rasterised at centres of 0.04° cells (2×2 of the nominal 0.02° ocean grid). Covered fine-grid ocean-pixel areas are summed by ring and divided by ring ocean area. The all-eddy union avoids double counting overlap; separate unions give cyclonic and anticyclonic coverage. Daily ring coverage fractions are averaged within each month.

The original occupancy preselection is preserved: distance is evaluated at contour vertices and their midpoints inside the regional box, with an initial centre-distance cutoff of 500 km; a contour is retained when its minimum such distance is ≤200 km. This sampling can miss an intersecting contour without sampled points in the domain. The occurrence map instead uses contour bounding-box intersection, so the map and ring series are related but do not use identical preselection. These historical choices are exposed in the upstream scripts, not silently harmonised.

The occurrence map is the fraction of all 6575 days for which a 0.04° location is inside any effective contour, multiplied by 100. It is not an eddy-centre count.

Monthly centroid: `C_m = Σ_k d_k f_mk / Σ_k f_mk`, where `f_mk` is ring coverage fraction and `d_k` its outer edge. Area is not multiplied again in this rate-weighted centroid. It is a centroid of monthly ring means, not a monthly median of daily centroids.

Annual share: `S_yk = 100 A_k Σ_(m in y) f_mk / Σ_j A_j Σ_(m in y) f_mj`. Each monthly mean receives equal weight, irrespective of month length. Thus this is a share of summed monthly occupied-area means, **not exact daily area-days**. Nearshore share is the sum for 0–60 km. Separate polarity shares are normalised within that polarity and are not additive across polarities.

## Offshore currents and positive-divergence proxy

GlobCurrent daily 0.25° product version 202411 provides Ekman (`ue`, `ve`) and total (`uo`, `vo`) velocity at 15 m (the primary depth). Total means geostrophic plus Ekman. Surface fields are present in upstream archives but are not part of the retained figures.

The offshore direction is a unit vector from the nearest point on the smoothed coastline to each current-grid centre, projected into local east and north. Coastline unit vectors are smoothed using 13 successive approximately 2 km vertices with edge padding, then normalised. Grid-centre offshore velocity is `u n_e + v n_n`, positive seaward. Contrary to one original configuration note, this calculation evaluates normals at current-grid centres, not ocean-pixel centroids.

Benguela uses daily projected velocities/divergence averaged to monthly grids in four source time chunks, combining split months by calendar-day weights. Gulf Stream uses monthly velocity components followed by projection and differencing. These linear operations commute with averaging only with consistent finite-value support; the products are mapped fields and missing-data handling is retained from the source implementation.

The divergence proxy is `∂u_E/∂x + ∂v_E/∂y`, using centred differences, Δy = R Δlatitude and Δx = R cos(latitude) Δlongitude. Outer grid edges are masked. This is the implemented local metric approximation; it omits the spherical `-v tan(latitude)/R` term and is not a measured vertical velocity or a vertically integrated Ekman pumping estimate.

Current cells finite in more than 95% of months define spatial support. Fine-grid ocean areas assign weights to the nearest 0.25° cell and ring. Remaining missing monthly values on that retained support are converted to zero in the archived implementation, with a fixed denominator; this can bias an average if residual gaps are present. The coverage fractions in the bundled arrays describe represented area, not resolution or an accuracy estimate.

Velocity means/trends retain all 20 rings (0–200 km). This corrects the proposal's blanket statement that all current metrics exclude 0–30 km: that exclusion is applied only to the divergence-share metric in the supplied figures. The inner 30 km are shaded in the publication plots as low reliability.

For upwelling shares, retain 30–200 km rings with finite monthly ring-mean divergence for every month. Clip the **ring-mean monthly** divergence to zero, then multiply by ring area, sum with equal monthly weights by year and normalise across retained rings to 100%. Clipping precedes annual averaging but follows spatial averaging; this is not an integral of pointwise positive divergence. The resulting metric is called a positive-divergence/upwelling proxy. Values are not vertical transport and do not establish nutrient fluxes.

## Trends and inference

All model fits are recalculated in `scripts/analyse.py`.

- Monthly centroid and velocity inference: OLS with intercept, linear time in years and eleven calendar-month indicators. Newey–West/Bartlett covariance uses 12 monthly lags, actual calendar separation and n/(n−13) scaling. No interpolation is used for these fits. There are 216 observations in complete series.
- Annual shares: OLS intercept plus centred year; Bartlett covariance with two annual lags and n/(n−2) scaling. There are 18 observations.
- Confidence intervals: estimate ± 1.96 HAC standard errors. Two-sided P value: `2 Φ(−|estimate / SE|)`. This is asymptotic normal inference, not a finite-sample t test.
- Descriptive centroid slope: Theil–Sen applied to the robust STL trend, with period=12, seasonal=999 and seasonal_deg=0, time in months/12. It is provided in `descriptive_sen` for centroid rows. The OLS interval and P value must not be called uncertainty for this different estimator.
- Figure significance colours use unadjusted 95% CIs excluding zero. No familywise or false-discovery-rate correction is applied. Neighbouring rings share mapped input cells, are correlated, and annual shares sum to 100%; ring tests are not independent.

The source workspace reported inflated rejection rates in simulations of its earlier STL/block-permutation approach, motivating separate OLS-HAC inference. Those exploratory simulations are not rerun or included as a validated result here. With only 18 years, borderline inference remains sensitive to covariance assumptions. HAC inference does not address sensor changes, uncertain coastal masks, effective resolution, product bias or physical causation.

Centroid slopes are km/year; current trends are cm/s/decade; share trends are percentage points/decade. Neither a relative percent trend nor a transport flux should be substituted for these units.
