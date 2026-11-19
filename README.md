# Conformational and actin memory can leave identical signatures in protrusion models

This repository rebuilds the numerical tests, the three main figures, the seven supporting figures and the manuscript and supporting-information reading copies. The biological application motivates a prospective test. The quantitative checks compare synthetic or analytical results and previously published source data.

## Reproduce from a fresh clone

The reference run used Python 3.14.7 and R 4.5.2. For figure and document generation install R with Cairo support, Pandoc, XeLaTeX, Poppler and the Arial font or a compatible substitute. The PDF template also supports Liberation Sans, Liberation Serif and Latin Modern Math. On Debian or Ubuntu, the system packages commonly needed include fontconfig, fonts-liberation, texlive-xetex, texlive-latex-extra, poppler-utils, libcairo2-dev, libfreetype6-dev, libharfbuzz-dev, libfribidi-dev, libpng-dev, libjpeg-dev and libtiff-dev.

    git clone --branch v1.0.0 https://github.com/renatosocodato/protrusion-memory-reproducibility.git
    cd protrusion-memory-reproducibility
    python3 -m venv .venv
    . .venv/bin/activate
    python -m pip install -r requirements.txt
    Rscript --vanilla scripts/install_r_dependencies.R
    python scripts/reproduce.py --output build/run-1

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

The [license](LICENSE) permits running the code to reproduce and verify the article. The original Amiri data and transcribed values retain their [source attribution and terms](figures/data/README.md).
