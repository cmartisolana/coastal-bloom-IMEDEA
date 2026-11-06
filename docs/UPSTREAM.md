# Upstream processing and reproducibility boundaries

## Offline reconstruction from included intermediates

Run `python scripts/rebuild_intermediates.py` after installing `requirements.txt`. For each region it runs only:

1. `upstream/contour_geometry.py`: the contour-distance preselection needed for eddy occupancy.
2. `upstream/eddy_occupancy.py`: daily contour unions and monthly ring occupancy for all, cyclonic and anticyclonic eddies.
3. `upstream/s19_occurrence_map.py`: the daily occurrence map.
4. `upstream/s28_ekman_metrics.py`: monthly grid-to-ring aggregation for 15 m Ekman velocity, total velocity and Ekman divergence.

Rebuilt arrays are compared with the frozen publication inputs using rtol=1e-7, atol=1e-10 and equal NaNs. A successful run writes `results/upstream_validation.txt`. Statistical calculations are centralised in `scripts/analyse.py`. No surface-depth sensitivity or unplotted centroid/divergence analysis is executed.

Outputs go to `results/rebuilt/`; the frozen publication inputs remain intact. No network access or extra dependencies are needed.

## Source-product provenance

Global products were processed in the original research workspace. The optional download and coast-grid reconstruction scripts have been removed from this publication repository. The reproducible workflow starts with the bundled regional contour extracts, monthly current grids and frozen coast geometry.

The original META4.0 source filenames are listed in `DATASETS.md`; extracts selected centres within each regional box padded by 3°, over 2002-09-01–2020-12-31, retaining effective contours. Only 2003–2020 enters the figures.

The original GlobCurrent versioned Zarr source was:

```
https://s3.waw3-1.cloudferro.com/mdl-arco-geo-037/arco/MULTIOBS_GLO_PHY_MYNRT_015_003/cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m_202411/geoChunked.zarr
```

The original readers decoded this snapshot's Blosc-compressed chunks. Benguela was projected and differentiated daily before monthly averaging in four time chunks; Gulf Stream monthly components were projected/differentiated during aggregation. The included intermediates preserve these acquisition routes. Global raw-data acquisition was not rerun for this handoff and is not part of the offline reproduction claim.

## Coast geometry

The supplied `ne10m_*` arrays and `rings_*` distance fields are fixed provenance inputs. The complete preparation of the Natural Earth source archive was not recoverable as a single versioned acquisition step. No end-to-end claim from an exact original Natural Earth release is made. Geometry and processing limitations are documented in `METHODS.md`.
