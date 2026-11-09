# Data and code availability: draft statements

Drafts for the manuscript, in the format of its existing *Data availability* and *Code availability* sections. Replace the bracketed items before submission.

**Data availability.** Eddy contours are from the Mesoscale Eddy Trajectory Atlas META4.0 (delayed-time, all-satellite), distributed by AVISO+ (https://doi.org/10.24400/527896/a01-2026.001). Total, geostrophic and Ekman currents at 15 m are from the Copernicus Marine product MULTIOBS_GLO_PHY_MYNRT_015_003, dataset cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m, version 202411 (https://doi.org/10.48670/mds-00327). The coastline is the Natural Earth 1:10 m vector (https://www.naturalearthdata.com). The Copernicus Marine product was accessed in November 2026 and META4.0 in [month] 2026. Regional eddy-contour extracts, monthly current grids, coastal ring geometry and the monthly and annual series underlying Figs. [x], [y] and [S1] are archived with the code (see Code availability).

**Code availability.** Code to compute eddy occupancy, offshore velocity and Ekman-divergence diagnostics, their trends and the figures, together with the regional input data needed to reproduce them, is archived at [Zenodo DOI] and available at https://github.com/cmartisolana/coastal-bloom-publication.

## Required acknowledgements

- **META4.0 (provider's citation text):** "The altimetric Mesoscale Eddy Trajectories Atlas (META4.0exp DT) was produced by SSALTO/DUACS and distributed by AVISO+ (https://aviso.altimetry.fr) with support from CNES, in collaboration with IMEDEA (DOI: 10.24400/527896/a01-2026.001 for the META4.0exp DT allsat version)." Cite also Gamot et al. (2026), https://doi.org/10.5194/essd-18-7043-2026.
- **Copernicus Marine:** "This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/mds-00327". Cite also Rio et al. (2014), https://doi.org/10.1002/2014GL061773.

## Before submission

- Archive a tagged release (e.g. via the GitHub–Zenodo integration or a manual Zenodo upload) to obtain the DOI above. Zenodo reads `CITATION.cff` for the record's metadata.
- Assign a licence (see `RIGHTS.md`).
- Add ORCID iDs to `CITATION.cff`.
