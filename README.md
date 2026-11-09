# Coastal bloom redistribution: eddies, currents and Ekman upwelling

Regional physical diagnostics for Benguela and the Gulf Stream, January 2003–December 2020: eddy occurrence, offshore currents and an Ekman upwelling proxy on 10 km coastal distance bands (0–200 km). Prepared by Cristina Martí-Solana and Simón Ruiz (IMEDEA, CSIC-UIB) for the manuscript *Global coastal phytoplankton bloom distributions are shifting offshore*. These are regional physical diagnostics, not causal evidence for the bloom changes.

| Requested | Where |
|---|---|
| Figures (editable vector) | `figures/`: PDF and SVG, plus 300 dpi PNG; captions in [docs/CAPTIONS.md](docs/CAPTIONS.md) |
| Data underlying the figures | `data/processed/` (figure inputs) and `results/tables/` (monthly and annual series as CSV); see [docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) |
| Trend estimates, 95% CIs and P values | `results/tables/all_trends.csv`; main results in `results/tables/key_results.csv` |
| Datasets, processing and statistical methods | [docs/METHODS.md](docs/METHODS.md) |
| Code for reproducibility | `scripts/`, `processing/`, `src/`; commands below |

## Figures

- `eddies`: eddy occurrence map, cross-shore centroid of eddy occupancy (all, cyclonic, anticyclonic) and its trend, and trends in each distance band's share.
- `currents_Ekman`: mean offshore Ekman velocity at 15 m, trends in offshore Ekman and total (geostrophic + Ekman) velocity, and trends in the share of Ekman upwelling (30–200 km).
- `eddy_polarity_share`: distance-band share trends for cyclonic and anticyclonic eddies separately.

## Reproduce

Uses [uv](https://docs.astral.sh/uv/); the Python version and exact package versions are fixed in `.python-version`, `pyproject.toml` and `uv.lock`. From the repository root:

```bash
uv sync                                         # create .venv with the locked environment
uv run python scripts/reproduce.py              # statistics (results/tables/) and figures (figures/)
uv run python scripts/rebuild_intermediates.py  # optional: recompute figure inputs from the regional extracts
```

`reproduce.py` uses the bundled figure inputs and needs no network access once the environment is installed. `rebuild_intermediates.py` recomputes eddy occupancy, occurrence maps and current/upwelling band averages from `data/intermediate/` into `results/rebuilt/` and checks them against `data/processed/`.

## Main caveats

- The coastal domains approximate the manuscript's regional boxes (Table S1); they are not the authors' coastal masks.
- Satellite products are least reliable within about 30 km of the coast; this band is shaded in the figures.
- P values are not adjusted for multiple comparisons, and neighbouring distance bands are not independent.
