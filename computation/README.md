# Numerical reproduction

From the repository root, install the pinned Python dependencies and run:

```sh
python -m pip install -r requirements.txt
python scripts/run_computations.py --output build/computations --scope all
```

The runner uses a single numerical thread and expects a new or empty output directory. The output includes JSON summaries, parent-level synthetic predictions, diagnostic plots, per-script logs and `reproduction_report.json`. The R renderers produce the manuscript figure panels.

`--scope main` is the default and runs the main-text calculations offline. Both `--scope all` and `--scope supplement` include those dependencies and run the supplementary calculations. Use `--skip-public` for a fully offline numerical run. This skips the published Amiri comparisons and their paired diagnostic plot.

| Result | Calculation directory |
| --- | --- |
| Figure 1 model equivalence and pulse-order challenge | `joint_fret_actin_mechanics/`, `frontier_conformation_chaos/run_frontier.py` |
| Figure 2 additional gate readout and optical references | `joint_fret_actin_mechanics/independent_gate_readout.py`, `frontier_conformation_chaos/run_frontier.py` |
| Figure 3 stationary availability, feedback stability and autonomous dynamics | `wasl_omics_gate/`, `frontier_conformation_chaos/feedback_boundary.py`, `autonomous_ensemble_chaos/` |
| Supplementary Figure 1 paired public U2OS traces and model dynamics | `mechanochemical_disorder/` |
| Supplementary Figure 2 autonomous dynamics and independent solver | `autonomous_ensemble_chaos/` |
| Supplementary Figure 3 molecular redistribution and Supplementary Figure 4 organization classifier | `molecular_organization_protrusion/` |
| Supplementary Figures 5 and 6 spatial and nonlinear amplitude challenges | `spatial_molecular_protrusion/` |
| Supplementary Figure 7 conformational-gate classifier | `conformational_ensemble/` |
| Supplementary lysed-sensor reference comparison | `wasl_paired_reconciliation/calibration_reference_screen.py` |

The runner compares selected outputs against the values reported by the original calculations. Counts and symbolic residuals are exact checks. Floating-point summaries have stated absolute tolerances in `reference_checks.json`. Chaotic trajectories can diverge pointwise across numerical environments, so the autonomous benchmark checks the finite-time exponents and the original numerical criteria. For Supplementary Figure 7 at m = 1.0, use `results_v2_fresh/results.json`, which contains the replacement held-out cohort. The first run still supplies the m = 0.3 comparison.

The default supplementary run fetches the public Amiri Source Data workbook from BioStudies and verifies SHA-256 `dcf21574fbcd5318d1581b3e1c80bd77893c9214dd871a71f174bda49f23b907`. To use a downloaded copy, pass `--amiri-workbook PATH`. Figure 2 of the source workbook supplies three representative U2OS time series, each from one cell. Its `SFig 8` sheet supplies the descriptive reference distributions. [Amiri et al., Nature Communications 14, 8011 (2023)](https://doi.org/10.1038/s41467-023-43612-5), [Source Data workbook](https://www.ebi.ac.uk/biostudies/files/S-EPMC10695988/41467_2023_43612_MOESM8_ESM.xlsx)

The abundance calculations use published expression context and specified hypothetical transfer and localization factors. The unpublished proteomic context is outside the numerical reproduction. [Human Protein Atlas v21, WASL](https://v21.proteinatlas.org/ENSG00000106299-WASL/single+cell+type/brain), [Mazzolini et al.](https://doi.org/10.1002/glia.24154)

Use `--include-exploratory` on the computation runner to replay abundance scenarios retained outside the article. The default run checks the stationary comparison used in Figure 3b.
