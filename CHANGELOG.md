# Changelog

## 1.8.0

- Added a periodic, linear reaction–diffusion ring with modal projections and the existing detector acquisition laws.
- Separated centroid displacement from radial deformation in public G01 membrane masks.
- Tested pathway recovery from modes 2–4 under shared and mode-specific gains.

## 1.7.0

- Added seeded finite-dimensional Itô dynamics with event-aligned Euler–Maruyama, shared or independent Wiener drivers and optional detector noise.
- Added a stochastic ensemble command and a direct-versus-mediated pulse challenge under white and temporally correlated noise.
- Retained the deterministic schemas and the original Biophysical Journal analysis reference.

## 1.6.0

- Added causal reporter response, exposure, delay and linear channel mixing through `ode-v2`.
- Enabled v2 models in simulation, comparison, calibration and scoring.
- Added a direct-versus-mediated pulse benchmark and replay of G01 and H1299 actin–motion recordings.

## 1.5.1

- Synchronized the title, text, significance statement and reading-copy layout with the final article prepared for Biophysical Journal.
- Linked the companion manuscript in the software citation and retained v1.1.0 as the reference for the original article analyses.

## 1.5.0

- Added a preparation-held-out edge forecast for signed velocity after a scheduled perturbation. Both fitted models receive two motion lags and the pulse history. One additionally receives local actin.
- Added paired synthetic perturbation controls and an extractor for a public two-channel vesicle recording. Biological forecast validation requires paired series with preparation identities.
- Preserved the paper and its original software citation.

## 1.4.0

- Added causal-prefix forecast fitting and held-out scoring on paired actin and motion records.
- Added the H1299 benchmark, with stationary, persistence, movement-history and actin-modulated forecasts.
- Preserved the paper and its original software citation.

## 1.3.0

- Added bounded calibration and separate scoring of frozen fits on held-out preparations.
- Added a checked external-data adapter and a preparation-level CK666 migration benchmark.
- Retained the paper's versioned reference results and original software citation.

## 1.2.0

- Added declarative ODE models with piecewise parameter modulations and point or interval-rate readouts.
- Added deterministic comparison of aligned observed trajectories and an external fibroblast benchmark with Rac1 photoactivation.
- Extended the software-use permission to research on laboratory-supplied models and data.

## 1.1.0

- Added the installable Protrusion Memory package and the simulate, compare and reproduce commands.
- Published the manuscript and supporting information, plotting materials and full reconstruction workflow as a versioned paper bundle.
- Shared the exact paired-model implementation with the paper's computation scripts. Frozen numerical reference values and figure data are preserved.

## 1.0.0

- Archived the original end-to-end reproduction workflow and figure inputs.
