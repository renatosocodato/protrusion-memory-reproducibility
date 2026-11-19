# Protrusion Memory

Companion software for “Conformational and actin memory can leave identical signatures in protrusion models.” The software simulates the exact paired models, compares the information supplied by a calibrated reporter and rebuilds the numerical tests, the three main figures, the seven supporting figures and the manuscript and supporting-information reading copies. The biological application motivates a prospective test. The quantitative checks compare synthetic or analytical results and previously published source data.

## Install and run

Download the wheel and paper bundle from [GitHub Release v1.2.0](https://github.com/renatosocodato/protrusion-memory-reproducibility/releases/tag/v1.2.0). Python 3.13 or 3.14 is supported by the pinned numerical dependencies. The two modeling commands need only the Python installation.

    python -m pip install protrusion_memory-1.2.0-py3-none-any.whl
    protrusion-memory simulate
    protrusion-memory compare

Commands write JSON to standard output or to a new file passed with `--output`. `simulate` returns the latent gate, actin, deformation and original optical precursor for both models. `compare` reports the smallest standardized mean separation over the chosen pulse orders and all integer sampling phases. It reports correct, wrong and inconclusive decision probabilities at the closest means when the intervals are separated. At touching means it reports the probability of either declaration under the common mean. Overlapping means or a zero sampled template leave no uniform family distinction.

To change the design, copy `examples/exact-pair.toml` or `examples/reporter-calibration.toml` from the paper bundle and pass it with `--config`. Pulse windows are half-open and must lie on the one-second grid. `simulate --order vu` swaps the signal windows and retains each input's amplitude. The ratio of the fast to slow time constants sets the reporter's touching bound. The calculations assume unit gain, exact subtraction of the common actin trace and independent Gaussian noise at a fixed variance per sample. These synthetic probabilities apply to the specified observation model.

## Use your own model

The [ode-v1 specification](PORTABLE_MODELS.md) describes state equations, parameter modulations and observation laws. Start from the [external Rac1 edge-velocity benchmark](benchmarks/rac1/README.md), with a configurable example, a stationary comparator and separate calibration and held-out cells.

## Rebuild the paper

Extract the paper bundle and then run:

    python -m pip install -r protrusion-memory-paper-1.2.0/requirements.txt
    Rscript --vanilla protrusion-memory-paper-1.2.0/scripts/install_r_dependencies.R
    protrusion-memory reproduce --target paper --bundle protrusion-memory-paper-1.2.0 --output build/paper

The paper bundle includes the reading copies, vector figures, panel exports, plotting inputs, all source, and a manifest checked against the listed files before reconstruction. Verify the downloaded archive against `SHA256SUMS` on the release page. The wheel supplies the modeling interface and entry point. The paper rebuild uses the R and document tools described below. A trusted Git checkout of the same version can serve as `--bundle`.

## Reproduce from a fresh clone

The reference run used Python 3.14.7 and R 4.5.2. For figure and document generation install R with Cairo support, Pandoc, XeLaTeX, Poppler and the Arial font or a compatible substitute. The PDF template also supports Liberation Sans, Liberation Serif and Latin Modern Math. On Debian or Ubuntu, the system packages commonly needed include fontconfig, fonts-liberation, texlive-xetex, texlive-latex-extra, poppler-utils, libcairo2-dev, libfreetype6-dev, libharfbuzz-dev, libfribidi-dev, libpng-dev, libjpeg-dev and libtiff-dev.

    git clone --branch v1.2.0 https://github.com/renatosocodato/protrusion-memory-reproducibility.git
    cd protrusion-memory-reproducibility
    python3 -m venv .venv
    . .venv/bin/activate
    python -m pip install . -r requirements.txt
    Rscript --vanilla scripts/install_r_dependencies.R
    protrusion-memory reproduce --output build/run-1

The full run downloads the published Amiri source workbook from EMBL-EBI and verifies its SHA-256 digest. Supply a locally downloaded original with `--amiri-workbook PATH` for a run without network access. The output directory for each build must be new or empty. The computation runner uses a single numerical thread per process and applies a 150-second limit to each script.

    python scripts/reproduce.py --artwork-only --output build/artwork
    python scripts/run_computations.py --scope main --output build/main-analysis

The first command uses the supplied figure input tables and produces the artwork and documents without running the models or downloading the workbook. The second runs just the main calculations and needs no external data. Use the wrapper for typesetting so the source files and reference figures resolve relative paths.

| Output | Contents |
| --- | --- |
| `computations/` | Numerical outputs and a per-script report with reference checks |
| `figure_inputs/` | Plotting quantities derived and compared with the supplied reference inputs |
| `article/submission_biophysj/figures/` | PDF, editable SVG and 600 dpi PNG for each figure and each panel |
| `article/submission_biophysj/` | Main manuscript and Supporting Information in PDF and Word |

See [Results and source map](RESULTS_MAP.md) for the computations behind the figures and the interpretation limits. The source-data credit for the paired U2OS traces is in [figure data](figures/data/README.md). The original workbook is fetched from its publisher deposit and retains its original terms.

## Numerical and visual tolerances

The conditional kinetic comparisons use fixed seeds and independent streams for test draws. Chaotic trajectories are sensitive to floating-point differences. The checks compare the early trace and late variability, the bounded state range and Lyapunov estimates. Long chaotic trajectories can separate pointwise on another machine. For a comparison with the exact artwork displayed in the manuscript, use the supplied figure inputs. The fonts and PDF metadata can also differ between systems. Reproduction was exercised on macOS. Linux-compatible fonts and tools are specified for a portable build.

The restricted mouse-brain proteomic dataset is cited as background. No restricted measurements or access credentials are needed for any computation or figure.

Use `--include-exploratory` on the computation runner to replay abundance scenarios retained outside the article. The default run checks the stationary comparison used in Figure 3b.


## Reuse

The [license](LICENSE) permits running the code to reproduce and verify the article and to research models and data from your own lab. The original Amiri data and transcribed values retain their [source attribution and terms](figures/data/README.md).
