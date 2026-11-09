# Data and code availability: draft statements

Drafts for the manuscript, in the format of its existing *Data availability* and *Code availability* sections. Replace the bracketed items before submission.

**Data availability.** Eddy contours are from the META4.0 mesoscale eddy atlas distributed by AVISO+ (https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product/meta4-0-dt.html). Total and Ekman currents at 15 m are from the Copernicus Marine product MULTIOBS_GLO_PHY_MYNRT_015_003 (https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_MYNRT_015_003/description). The coastline is the Natural Earth 1:10 m vector (https://www.naturalearthdata.com). Regional extracts, coastal ring geometry and the series underlying the figures are archived with the code (see Code availability).

**Code availability.** Code to compute eddy occupancy, offshore velocity and Ekman-divergence diagnostics, their trends and the figures, together with the regional input data needed to reproduce them, is archived at [Zenodo DOI] and available at https://github.com/cmartisolana/coastal-bloom-publication.

## Required acknowledgements

- **META4.0 (provider's citation text):** "The altimetric Mesoscale Eddy Trajectories Atlas (META4.0exp DT) was produced by SSALTO/DUACS and distributed by AVISO+ (https://aviso.altimetry.fr) with support from CNES, in collaboration with IMEDEA (DOI: 10.24400/527896/a01-2026.001 for the META4.0exp DT allsat version)."
- **Copernicus Marine:** "This study has been conducted using E.U. Copernicus Marine Service Information; https://doi.org/10.48670/mds-00327".

## Before submission

- Archive a tagged release (e.g. via the GitHub–Zenodo integration or a manual Zenodo upload) to obtain the DOI above. Zenodo reads `CITATION.cff` for the record's metadata.
- Assign a licence (see `RIGHTS.md`).
- Add ORCID iDs to `CITATION.cff`.
