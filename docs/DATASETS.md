# Datasets

| Product | Used here |
|---|---|
| [META4.0 DT allsat mesoscale eddy atlas (AVISO+)](https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product/meta4-0-dt.html) | Effective contours of cyclonic and anticyclonic eddies, 2003–2020 (files `META4_DT_allsat_{anticyclonic,cyclonic}_19930101_20230908.nc`) |
| [Copernicus-GlobCurrent, MULTIOBS_GLO_PHY_MYNRT_015_003 (Copernicus Marine)](https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_MYNRT_015_003/description) | Dataset `cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m`, version 202411: total and Ekman currents at 15 m, 2003–2020 |
| [Natural Earth 1:10m physical vectors](https://www.naturalearthdata.com/downloads/10m-physical-vectors/) | Coastline and land, prepared as the frozen arrays in `data/processed/ne10m_*.npz` |

The global source files are not redistributed. The repository bundles regional extracts (`data/intermediate/`) and the derived figure inputs (`data/processed/`); see [UPSTREAM.md](UPSTREAM.md) for how they were produced. Bloom rasters and the manuscript authors' masks are not used.
