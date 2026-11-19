# Figure inputs

The frozen tables below contain the coordinates, counts or estimates used by the
R figure scripts. The synthetic parent identifiers in the held-out tables
preserve which candidate classes share a parameter draw.

| Input | Figure | Content |
| --- | --- | --- |
| `classification.json` | F01b, F01d | Exact latent gate discrepancy and synthetic classifier accuracies |
| `pulse_order_counts.csv` | F01e | Joint classification counts before and after reversing the pulses |
| `independent_gate_readout.json` | F02c–e | Synthetic decision fractions for the stated actin bounds, sampling intervals and reporter response times |
| `optical_reference.json` | F02f | Held-out median relative ratio errors before and after unmixing synthetic measurements |
| `chaos_exponents.json` | F03d, S02c | Largest finite-time Lyapunov exponent for each displayed model case |
| `S01_amiri_traces.json` | S01a,b | Three published representative U2OS traces, each from one cell |
| `S01_load.csv` | S01c | Synthetic dimensionless load by driving condition |
| `S02_state.csv` | F03d, S02a,b | Synthetic autonomous state trajectory |
| `S03_trajectory.csv` | S03 | Synthetic trajectories by class and assumed mean concentration |
| `S04_heldout_predictions.csv` | S04 | True and predicted class for each synthetic parent, channel and concentration |
| `S05_second_challenge_parents.csv` and `S05_third_challenge_parents.csv` | S05 | Spectral estimates, noisy force changes and fixed-rule calls in two synthetic challenges |
| `S06_fourth_challenge_parents.csv` | S06 | Nonlinear ratio and induced force change for synthetic parents |
| `S07_conformation.json` | S07 | Held-out synthetic accuracies and paired bootstrap gains, using the reserved replacement test at m = 1.0 |

## Published U2OS traces

Amiri S, Muresan C, Shang X, et al. Intracellular tension sensor reveals
mechanical anisotropy of the actin cytoskeleton. *Nature Communications* 14,
8011 (2023). [Article and Source Data](https://doi.org/10.1038/s41467-023-43612-5).

`S01_amiri_traces.json` contains values transcribed from the Source Data
workbook for Figure 2c, d and f, in the order time in minutes, mean traction
magnitude in pascals, and FRET efficiency on a unit scale. The renderer converts
efficiency to percent and leaves the 30-minute sampling gap unconnected. The
figure is redrawn from those values. The workbook checksum is retained in the
JSON.
