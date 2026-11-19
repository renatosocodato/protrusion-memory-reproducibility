# CK666 across independent ASPC1 acquisition series

The [cellmig study](https://doi.org/10.1371/journal.pcbi.1014472) provides an external test of preparation-level transfer. The software calibrates a condition-dependent log migration-speed shift for CK666 at 40 µM using the two Dataset 1 experiments, with eight plates in total. It then evaluates the frozen shift in the independent Dataset 2 screen. Dataset 2 contributes 18 plates, scored in 14 conservative preparation groups after pairing plates with matching g1 stems. Each plate's DMSO wells define a common baseline for both candidates.

The calibrated shift was −0.554 log units, giving a velocity ratio of 0.575 relative to same-plate DMSO. The RMSE in Dataset 2 was **0.117 log units with the CK666 shift** and **0.546 for the stationary reference**. The reduction was 78.6%, with lower mean squared error in all 14 preparation groups. Omitting each group in turn gave modulated-minus-stationary RMSE differences from −0.464 to −0.418 log units. This range describes sensitivity to one group at a time.

The paper already reported the compound's direction and magnitude in both datasets. This is retrospective transfer of a known positive control, conditional on a measured same-plate baseline. The endpoint is aggregate cell migration at one dose. The experiment does not resolve edge dynamics, temporal modulation or actin memory.

## Reproduce

Run from a checkout or extracted paper bundle with Protrusion Memory and base R installed. The adapter checks the original source hashes before loading each file. The two commands fetch only the source relevant to the requested split.

    python benchmarks/ck666/prepare.py --split calibration --output ck666-calibration.csv
    protrusion-memory calibrate --plan benchmarks/ck666/calibration.toml --data ck666-calibration.csv --output ck666-fit.json
    python benchmarks/ck666/prepare.py --split heldout --output ck666-heldout.csv
    protrusion-memory score --plan benchmarks/ck666/calibration.toml --fit ck666-fit.json --data ck666-heldout.csv --output ck666-score.json

Use `--source` to point the adapter to the pinned original RData file when offline. [Fixed design](DESIGN.md) · [Reference fit and preparation-level results](reference.json) · [Original data](https://github.com/snaketron/cellmig_paper_data/tree/6161c370d7a616dde35a12f2a83fff4ee6e713aa/data)
