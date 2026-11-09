# Datasets and provenance

## META4.0 eddy atlas

AVISO+ Mesoscale Eddy Trajectory Atlas META4.0, delayed-time all-satellite, DOI [10.24400/527896/a01-2026.001](https://doi.org/10.24400/527896/a01-2026.001). The provider's [catalogue](https://www.aviso.altimetry.fr/fr/donnees/produits/produits-a-valeur-ajoutee/trajectoires-de-tourbillons-mesoechelles.html) lists this product and DOI (checked 2026-11-06).

Original local filenames:

- `META4_DT_allsat_anticyclonic_19930101_20230908.nc`
- `META4_DT_allsat_cyclonic_19930101_20230908.nc`

These filenames identify the source snapshots; the temporal interval is broader than the analysis. Regional extracts under `data/intermediate/meta4_*.npz` retain time, centre coordinates and effective contours only. Global NetCDF files are not redistributed. Atlas contours arise from mapped altimetry; neither the contour vertex spacing nor the analytical band width establishes effective resolution.

## GlobCurrent

Copernicus Marine product `MULTIOBS_GLO_PHY_MYNRT_015_003`, dataset `cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m`, version `202411`. Daily 0.25° total and Ekman currents; primary analysis at 15 m. Exact product/version identifiers and acquisition endpoint are recorded in the original configuration and documented below.

- [Product user manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-MOB-PUM-015-003.pdf)
- [Quality information document](https://documentation.marine.copernicus.eu/QUID/CMEMS-MOB-QUID-015-003.pdf)
- [Copernicus Marine catalogue](https://data.marine.copernicus.eu/)

The provider documents were located on 2026-11-06. The private handoff preserves the **202411 snapshot**; no claim is made that it is the current recommended version. The original public ARCO endpoint is listed in `UPSTREAM.md`. It was not redownloaded during validation.

Included grids: four Benguela chunks of projected monthly current and divergence means; one Gulf Stream monthly component archive. Retaining these allows offline verification of ring aggregation and the manuscript figures. Full daily global currents are not bundled.

## Coastline and masks

Natural Earth 1:10 million coastline/land, supplied as the exact archived prepared NPZ arrays. The original workspace reports a coastline-length filter of 1428 km and resampling at 2 km. No upstream Natural Earth archive checksum or release version was recorded in the files inspected, so the frozen prepared coastline is the reproducible geometry input. Do not claim an exact Natural Earth release or an exact match to the manuscript masks.

## Excluded inputs

Bloom rasters and the manuscript authors' exact masks are not required to reproduce these three figures and are not bundled. No independent replication of the biological bloom analysis is claimed. The original draft's bloom-reference slopes have been removed from the new figure legend; this repository does not provide their source observations.
