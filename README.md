# Coastal bloom redistribution: regional physical diagnostics

Reproducible eddy-occurrence, offshore-current and Ekman-divergence analyses for Benguela and the Gulf Stream, January 2003–December 2020. This companion repository contains the three figures in the contribution by Cristina Martí-Solana and Simón Ruiz (IMEDEA, CSIC-UIB). It describes regional physical context, not a causal attribution of phytoplankton bloom changes.

## Handoff to manuscript authors

| Requested item | Location |
|---|---|
| Figures 1, 2 and S1 | `figures/`: PDF, editable SVG, 300 dpi PNG |
| Underlying figure data | `data/processed/`; readable monthly and annual CSVs in `results/tables/` |
| Trend estimates, 95% CIs, P values | `results/tables/all_trends.csv` (260 estimates) |
| Datasets, processing and statistics | [Methods](docs/METHODS.md), [datasets](docs/DATASETS.md), [data dictionary](docs/DATA_DICTIONARY.md) |
| Figure captions | [Captions](docs/CAPTIONS.md) |
| Reproducibility | Commands below, [validation](results/validation.txt), [source audit](docs/AUDIT.md) |

## Reproduce

Python 3.12 is the tested interpreter. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/reproduce.py
```

This recalculates all statistics, checks them against the original calculations, and regenerates the three figures without network access. PDF and SVG retain vector axes, text, curves and bars; the dense map fields are embedded rasters to keep them practical to edit. SVG text remains text. Plot inputs are supplied so maps can also be restyled.

To verify processing one level further upstream, using the included regional contour extracts and monthly current grids:

```bash
python scripts/rebuild_intermediates.py
```

This takes longer and writes to `results/rebuilt/`, leaving the bundled inputs intact. It rebuilds eddy occupancy, occurrence maps and current/divergence ring averages and compares them with the figure inputs. See [upstream processing](docs/UPSTREAM.md) for source-product provenance and reproduction boundaries. The large original global NetCDF products are not bundled.

## Scientific interpretation

The original point estimates and calculations are preserved. Centroid curves are robust STL trends; descriptive Theil–Sen slopes are reported separately from the OLS estimates that supply the confidence intervals and P values. Figure 1 labels now show **OLS estimates with matching OLS intervals**, correcting the original mixed-estimator presentation. Both sets of estimates are in the tables.

Annual contribution shares use **equal monthly weighting**, matching the submitted analysis. The approximate coastal domains are not the manuscript authors' exact masks. Velocity trends retain the original 0–200 km range, with 0–30 km shaded as low reliability; the positive-divergence share excludes 0–30 km. Read the methods and audit before quoting results.

This repository deliberately excludes EKE figures, other regions, particle/trajectory experiments, preliminary plotting variants, manuscript editorial correspondence and unrelated code. The source repository is unchanged.

## Authorship and use

Cristina Martí-Solana and Simón Ruiz, IMEDEA (CSIC-UIB), Esporles, Spain, as stated in the contribution document. No co-authorship acceptance or correspondence is implied by this repository. No software or derived-data licence has been assigned on the authors' behalf; see [rights and attribution](docs/RIGHTS.md). This handoff is intended for private collaboration under `cmartisolana`.
