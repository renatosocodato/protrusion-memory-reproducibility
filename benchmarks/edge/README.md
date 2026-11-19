# Does edge-local actin improve a forecast of edge motion?

The primary comparison is between two forecasts of signed edge velocity, evaluated on the same origins in reserved preparations. Both receive two motion lags and the timing of the independently scheduled pulse. The augmented forecast receives one additional observation, the local actin signal at the same edge patch. See the [forecast specification](../../FORECAST.md#scheduled-edge-perturbations) for input and timing requirements.

## Synthetic perturbation challenge

Run from a checkout after installing Protrusion Memory 1.5.0.

    python benchmarks/edge/benchmark.py --output edge-synthetic-results.json

The generator defines two cases before scoring. In both, the pulse affects velocity and actin, so actin can mark the intervention even in the absence of a local contribution. Only the positive case also couples current local actin to next-frame velocity. Each of 40 realizations per case has six calibration and eight reserved preparations, each with two cells and two edge patches. Noise is shared across patches in each preparation. Reserved preparations have independently generated noise and later, longer pulses. The forecast horizon is 5 s. An RMSE reduction greater than 0.05 µm/min is the recovery criterion for these simulations.

| Synthetic generator | Median augmented minus motion-and-pulse RMSE | Range across 40 realizations | Reductions > 0.05 µm/min |
| --- | ---: | ---: | ---: |
| Pulse and motion only | +0.00007 µm/min | −0.00032 to +0.00091 µm/min | 0/40 |
| Added local actin term | −0.21609 µm/min | −0.24226 to −0.17007 µm/min | 40/40 |

The [results](synthetic-results.json) include every replicate contrast and the preparation-level differences for the first replicate. These rates characterize the two specified generators and noise regimes. Intermediate couplings, slow omitted pulse responses, optical confounds and biological transfer remain open.

## Independent experimental material

| Public source | Paired readouts and perturbation | What can be tested |
| --- | --- | --- |
| [PA-Rac1 mouse fibroblast edge data](https://github.com/DanuserLab/Protrusion_Retraction_Database) | Local edge velocity and optical stimulation. Public series lack paired actin and preparation IDs. | Motion reference only. |
| [RAW 264.7 opto-Cdc42, Movie S4](https://pmc.ncbi.nlm.nih.gov/articles/PMC4850032/) | Lifeact and activation in one representative cell. The red channel also reports GEF translocation. | A single visual example. |
| [Rapamycin GUVs](https://doi.org/10.5061/dryad.tqjq2bw66) | Paired actin and independent membrane marker. [The study](https://doi.org/10.1126/sciadv.adk9731) reports 18 globally stimulated vesicles from three independent experiments. The available metadata do not map those vesicles to their experiments. | Causal extraction demonstration on the accessible raw G01 movie. |

For the G01 extraction, use the pinned [original TIFF files](https://github.com/basharif/GUV-symmetry-breaking-analysis/tree/56d7a4eb5cd3c2097c4e08ab138734dcf1e6edd3/G01_input/G01_raw). The script verifies both hashes. The example requires Pillow in addition to the core package.

    python -m pip install Pillow
    python benchmarks/edge/guv_g01.py G01_mCh_memb.tif G01_YFP_actin.tif guv-g01.csv

In the source's 28-frame analysis window at 65.589 s per frame, the extractor returns 202 paired sector-frame readouts of actin and signed radial membrane velocity. For each frame it identifies the mCherry ridge along fixed rays around the first-frame center, with a search confined to the vicinity of the previous frame's ridge. It reads YFP-actin 1–2 pixels inward and differences the sector-average radius with the preceding frame. It leaves a value blank if any ray fails the search-window check. This is a demonstration of timestamped signal extraction in a reconstituted system. The apparent radial velocity can include vesicle translation. Sector measurements share one preparation, so no fit or performance estimate is reported for G01.

An empirical comparison needs dual-channel movies mapped to preparation and pulse timestamps. Freeze the edge definition, exposure and recovery windows, motion lags, reporter availability and forecast horizon on calibration preparations. Reserve whole independent preparations before scoring. With the two cohorts in place, run the same `forecast-fit` and `forecast-score` commands and report paired errors by preparation.
