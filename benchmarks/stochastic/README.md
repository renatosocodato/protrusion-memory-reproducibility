# Can a recorded pulse distinguish direct from actin-mediated displacement?

The direct response and the actin-mediated response are designed to have the same integrated displacement. The [acquisition benchmark](../acquisition/README.md) showed that their recorded shapes depend on when the camera captures the pulse. Here we add noise and allow an unknown nonnegative amplitude for each rival. The two rivals have the same actin input response. This contrast asks whether the sampled shape response can locate how the pulse drives motion.

The first six frames of the public G01 vesicle gave a median temporal radial SD of **0.225 µm** across its eight sectors. The frame interval was **65.589 s**. This is a scale of observed radial variation in one specimen, which includes biological fluctuations and image-extraction variation. It anchors an assigned SD for the intervention-minus-control samples. The [Dryad archive](https://doi.org/10.5061/dryad.tqjq2bw66) and [pinned G01 settings](https://github.com/basharif/GUV-symmetry-breaking-analysis/blob/56d7a4eb5cd3c2097c4e08ab138734dcf1e6edd3/target_settings.csv) document the acquisition. The H1299 dataset contributes a 300-s clock; its whole-cell area noise is incommensurate with radial displacement.

The synthetic study uses white detector noise or stationary, exposure-averaged Ornstein–Uhlenbeck displacement noise with a 20-s correlation time. They are matched on observed marginal variance. Detector optics, pulse phase and kinetic template shapes are assumed known, and a constant displacement offset and response amplitude are independently fitted under each candidate. Scores are evaluated on independent synthetic observations. See the [frozen design](DESIGN.md) and [confusion tables](results.json).

| Clock and noise scenario | Correct mediated assignment | Correct direct assignment |
| --- | ---: | ---: |
| 1 s, baseline scale, white or OU | 99.45–100% | 100% |
| 1 s, twice baseline scale, OU | 67.67–69.60% | 99.33–99.55% |
| 65.589 s, baseline scale, five regular phases and both noise forms | 0/4,000 in each setting | 0/4,000 in each setting |
| 65.589 s, baseline scale, near-frame alignment, both noise forms | 3.38–3.77% | 2.12–2.50% |
| 300 s, baseline scale, five phases and both noise forms | 0/4,000 in each setting | 0/4,000 in each setting |

At the G01 clock and baseline scale, some regular-phase settings detected a mechanical response in about 99% of traces before applying the mechanism-separation margin. The best-case classifier with known amplitudes reached 67% across the regular G01 phases. Detectability therefore says little about whether a direct or mediated pathway generated the response at that cadence. A separate known-amplitude result at 300 s reached 90.7% near the frame edge while the free-amplitude test still abstained. These are illustrative synthetic noise draws, not estimates of actual preparation-level performance.

## Reproduce

    python benchmarks/edge/guv_g01.py G01_mCh_memb.tif G01_YFP_actin.tif guv.csv
    PYTHONPATH=benchmarks/acquisition:benchmarks/edge python benchmarks/stochastic/resolve.py guv.csv

The benchmark uses the exact Gaussian covariance of the declared noise models and samples their two sufficient projections. The general-purpose SDE solver is checked independently with common-Wiener step refinement, OU moment comparisons and event-aligned acquisition tests.
