# Upstream processing and reproducibility boundaries

## Offline reconstruction from the included intermediates

Run `python scripts/rebuild_intermediates.py` after installing `requirements.txt`. For each region it runs:

1. `upstream/contour_geometry.py`: the contour-distance preselection from the original `s2_contact.py`, without its unrelated trajectory history calculations.
2. `upstream/eddy_occupancy.py`: daily contour unions and monthly ring occupancy, extracted from `s5_positional.py` without EKE, centre-density or contact-history metrics.
3. `upstream/s19_occurrence_map.py`: daily occurrence map.
4. `upstream/s28_ekman_metrics.py`: monthly grid-to-ring velocity/divergence aggregation.

Rebuilt arrays are compared with the frozen publication inputs using rtol=1e-7, atol=1e-10 and equal NaNs. A successful run writes `results/upstream_validation.txt`. The original upstream current aggregator still computes surface-depth and other diagnostic summaries internally; these are not used by `scripts/analyse.py` or the publication plots. The reference series are never overwritten by this check.

## Reacquiring source products

This is an optional network-dependent extension, not needed to rebuild the supplied figures. Install `requirements-upstream.txt`. Product credentials/access rules and snapshot availability depend on the provider.

Download the two META4.0 files named in `DATASETS.md` from AVISO+ into `data/raw/`. Then run:

```bash
python upstream/s2_extract_meta.py
```

This reader uses effective contours and the original padded boxes/time selection. It recreates regional extracts in `data/intermediate/` (with additional source metadata). Global files can be large; extraction was preserved from the source repository and full global acquisition was not rerun for this handoff.

The original GlobCurrent versioned Zarr base is:

```
https://s3.waw3-1.cloudferro.com/mdl-arco-geo-037/arco/MULTIOBS_GLO_PHY_MYNRT_015_003/cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m_202411/geoChunked.zarr
```

To reproduce the historical acquisition routes:

```bash
python upstream/s27_mob_ekman.py SW_Africa coords
python upstream/s27_mob_ekman.py SW_Africa 0
python upstream/s27_mob_ekman.py SW_Africa 1
python upstream/s27_mob_ekman.py SW_Africa 2
python upstream/s27_mob_ekman.py SW_Africa 3
python upstream/s27b_mob_components.py NA_east data/intermediate/mob_components_NA_east.npz
python scripts/rebuild_intermediates.py
```

The preserved readers decode the recorded Zarr chunk arrangement using Blosc. They are snapshot-specific, not general download clients. A newer endpoint or storage layout may require adaptation. HTTP 403/404 chunks are treated as unavailable/all fill in the historical readers; consequently, successful execution alone is not evidence of complete source coverage. Use the included numerical comparisons and inspect missing values when reacquiring.

## Coast geometry

The supplied `ne10m_*` arrays and `rings_*` distance fields are fixed provenance inputs. `upstream/eddy_ring_geometry.py --tag REGION` can reconstruct rings from the prepared coast/land arrays, but writes `data/processed/rings_REGION.npz`; use a disposable checkout for this optional geometry experiment. Matplotlib rasterisation and the float16 distance representation can affect exact boundary pixels. The complete preparation of the Natural Earth source archive was not recoverable as a single versioned acquisition step, so no end-to-end claim from an exact original Natural Earth release is made.
