# Selection and scientific audit

## Figure-to-code trace

| Contribution figure | Original code block | Retained implementation |
|---|---|---|
| Figure 1, eddies | `analysis_mesoscale_export/scripts/s30_definitive_figures.py`, Fig 1 | `scripts/plot_figures.py`, Fig 1 |
| Figure 2, currents | Same file, original Fig 3 | Renumbered Fig 2 |
| Figure S1, polarity shares | Same file, Fig S1 | Fig S1 |
| Statistics | `s17_equivalent_metrics.py`, `s26_polarity.py`, `s28_ekman_metrics.py`, `src/coastal_bloom/trends.py` | `scripts/analyse.py` and `src/coastal_bloom/trends.py` |
| Eddy occupancy input | Eddy section of `s5_positional.py`; contour-distance part of `s2_contact.py` | `upstream/eddy_occupancy.py`, `upstream/contour_geometry.py` |
| Occurrence map | `s19_occurrence_map.py` | Retained upstream |
| Current grids and rings | `s27_mob_ekman.py`, `s27b_mob_components.py`, `s28_ekman_metrics.py` | Frozen monthly grids plus the 15 m ring-aggregation portion of `s28_ekman_metrics.py` |

The original Fig 2 EKE block is excluded because the contribution explicitly excluded it. Other systems, kinetic-energy variants, jet-axis diagnostics, particle experiments, contact-history outcomes, synthetic examples, alternative figure scripts and manuscript editorial comments are excluded. No files in the original repository were deleted.

## Changes made for the handoff

1. Recomputed 260 relevant centroid, share and velocity estimates and added matching model-based two-sided P values.
2. Kept descriptive Theil–Sen slopes separately from OLS inference. Figure 1 legend now pairs OLS estimates with OLS CIs; the STL curves and descriptive Sen fit remain. This changes the displayed centroid estimate labels, not the original input data or calculation.
3. Corrected figure numbering to 1, 2 and S1. Added editable SVG and consistent low-reliability shading for near-coastal velocity trends.
4. Exported monthly/annual source series as CSV; bundled map arrays, coast geometry, eddy contour extracts and current-grid intermediates.
5. Removed unrelated arrays from eddy monthly/contour archives; converted trusted historical object-string arrays to Unicode to allow safe NumPy loading without pickle.
6. Added an offline rebuild, numerical regression checks against archived outputs, independent statsmodels covariance checks, dependency pins and CI.
7. Documented equal-month weighting, approximate domains, contour preselection, missing-value treatment, local divergence approximation and the actual 0–30 km exclusion.

## Limits of validation

Numerical reproduction establishes consistency with the submitted analysis. It does not establish that the coastal masks match the manuscript, that the proxy is vertical velocity, that mapped 10 km bins are independently resolved, or that the physical changes caused bloom redistribution. Source downloads are not automatically refreshed: a newer upstream product may produce different results. The exact global raw download path is provided for provenance, but full global-product acquisition is not part of the offline validation.

`source_manifest.json` records selected source paths and SHA-256 fingerprints; `array_schema.json` records every bundled NPZ array and shape. Adapted code is identified above; copied-source checksums describe its original source, not necessarily its edited final content. `checksums.sha256` records the actual final bundled data files.

## Further code reduction

Removed optional global-data acquisition scripts and the unused coast-grid reconstruction script. Removed surface-depth sensitivity calculations, unplotted divergence/centroid statistics, duplicate trend routines, unused eddy-centre allocations and unused plot imports/constants. All retained Python files are called by either the figure reproduction command or the offline upstream validation command. Frozen geometry and source/intermediate data remain included; full source-product downloading is outside this repository's execution workflow.
