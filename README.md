# Coastal bloom redistribution: eddies, currents and Ekman upwelling

Regional physical diagnostics for **Benguela** and the **Gulf Stream**, January 2003–December 2020, prepared by Cristina Martí-Solana and Simón Ruiz (IMEDEA, CSIC-UIB) for the manuscript *Global coastal phytoplankton bloom distributions are shifting offshore*. We ask whether eddies, offshore currents and wind-driven upwelling changed their distance from the coast over the same period as the blooms. These are regional physical diagnostics, not causal evidence for the bloom changes.

## Where to find what you asked for

| Requested | Where |
|---|---|
| Figures (editable vector) | [`figures/`](figures): PDF and SVG, plus 300 dpi PNG. Captions: [docs/CAPTIONS.md](docs/CAPTIONS.md) |
| Data underlying the figures | [`results/series/`](results/series) (CSV) and [`data/processed/`](data/processed) (maps and geometry). Described in [docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) |
| Trend estimates, 95% CIs and P values | [`results/trends/all_trends.csv`](results/trends/all_trends.csv) (all 260 estimates); main ones in [`key_results.csv`](results/trends/key_results.csv) |
| Datasets, processing and statistical methods | [docs/METHODS.md](docs/METHODS.md) |
| Code for reproducibility | [`scripts/`](scripts) and [`processing/`](processing); see "Reproduce" below |

## The three figures

| Figure | Shows |
|---|---|
| `eddies` | Where eddies occur (map), the cross-shore centroid of eddy occupancy over time (all, cyclonic and anticyclonic eddies), and the trend in each ring's share of eddy occupancy |
| `currents_Ekman` | Mean offshore Ekman velocity at 15 m (map), trends in offshore Ekman and total (geostrophic + Ekman) velocity per ring, and the trend in each ring's share of Ekman upwelling (30–200 km) |
| `eddy_polarity_share` | The ring-share trends of the eddy figure, separately for cyclonic and anticyclonic eddies |

## Key terms

- **Ring:** one of twenty 10 km bands of distance from the coast, 0–200 km.
- **Eddy occupancy:** fraction of a ring's ocean area inside eddy contours (META4.0), averaged over each month.
- **Centroid:** the occupancy-weighted mean distance from the coast; an increase means eddies sit further offshore.
- **Share:** a ring's percentage of the year's total over all rings; a rising share means that ring gains relative to the others.
- **Ekman upwelling proxy:** positive divergence of the 15 m Ekman current (wind-driven upwelling away from the coast).

## How the repository fits together

```
data/inputs/          regional extracts of META4.0 eddies and GlobCurrent currents
      │  processing/step1–4   (optional; scripts/run_processing.py)
      ▼
data/processed/       monthly values per ring, maps, coastal geometry
      │  scripts/compute_trends.py
      ▼
results/series/       monthly and annual series (CSV)
results/trends/       trend estimates, 95% CIs, P values
      │  scripts/plot_figures.py
      ▼
figures/              eddies, currents_Ekman, eddy_polarity_share
```

| Folder | Contents |
|---|---|
| `config/regions.json` | Region names and bounding boxes |
| `processing/` | Step 1 eddy contour distance to the coast; step 2 eddy occupancy per ring; step 3 eddy occurrence map; step 4 offshore currents and Ekman divergence per ring |
| `scripts/` | `reproduce.py` (trends + figures), `compute_trends.py`, `plot_figures.py`, `run_processing.py` (steps 1–4 with a check) |
| `src/coastal_bloom/trends.py` | Monthly trend regression with HAC standard errors |
| `docs/` | Methods, captions, data dictionary |

## Data sources

- [META4.0 DT allsat mesoscale eddy atlas](https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product/meta4-0-dt.html) (AVISO+)
- [Copernicus-GlobCurrent, MULTIOBS_GLO_PHY_MYNRT_015_003](https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_MYNRT_015_003/description) (Copernicus Marine)
- [Natural Earth 1:10m coastline and land](https://www.naturalearthdata.com/downloads/10m-physical-vectors/)

## Reproduce

Requires [uv](https://docs.astral.sh/uv/); Python and package versions are fixed in `.python-version`, `pyproject.toml` and `uv.lock`. From the repository root:

```bash
uv sync                                   # create the environment
uv run python scripts/reproduce.py        # trends (results/) and figures (figures/)
uv run python scripts/run_processing.py   # optional: rerun processing steps 1–4 (slower)
```

`reproduce.py` starts from `data/processed/` and needs no network access once the environment is installed. `run_processing.py` starts from `data/inputs/`, writes to `results/rebuilt/` and checks the result against `data/processed/`.

## Main caveats

- The coastal domains approximate the manuscript's regional boxes (Table S1); they are not the authors' coastal masks.
- Satellite products are least reliable within about 30 km of the coast; this band is shaded in the figures.
- P values are not adjusted for multiple comparisons, and neighbouring rings are not independent.
