# External benchmark: Rac1 and fibroblast edge velocity

This case uses mouse embryonic fibroblasts from Figure 7 of Sapoznik et al. (2020), [eLife 9:e57681](https://doi.org/10.7554/eLife.57681). The [Danuser Lab source](https://github.com/DanuserLab/Protrusion_Retraction_Database/tree/a78f2d8d46d8e03a85cf8d568b5ace10fbaa8e72/rawData) provides signed edge velocity in nm/s, frame and window labels, and a marker for illuminated windows. The example fetches two pinned source files or verifies local originals. It does not redistribute those datasets.

## Fixed design

Use layer 1, the authors' illuminated-window mask, and the matched group and cell index across the before and activation epochs. Reduce windows to one mean per frame with observed values only and exclude the first frame of each epoch, when the velocity has no preceding interval. Use 5 s intervals and treat the activation epoch as an effective sustained drive. The readout is the position difference over each interval divided by 5 s. The exact optical pulse sequence within the epoch and an intracellular Rac1 trajectory are unavailable.

The model approaches an assigned target velocity with a fixed 30 s time constant and integrates position. Subtract each cell's mean before illumination. Estimate a single effective drive amplitude from odd-numbered PA-Rac1 cells by averaging their individual least-squares amplitudes. Freeze that amplitude and score on even-numbered PA-Rac1 cells. The stationary comparator uses the same pre-light data and predicts zero deviation. Compute one mean squared error per held-out cell and average those errors with equal cell weights. Report the signed before-to-activation shift in even-numbered light-exposed controls as a separate descriptive check. A closed-form integrated response checks the portable numerical solver independently.

The benchmark tests portability to an external assay and withheld cells in that assay. The published HMM analysis uses a different endpoint. The assigned time constant and the relation between effective optical drive and Rac1 molecular dynamics are outside this test. Biological preparations and acquisition days are unavailable as grouping variables.

## Run and interpret

With the installed package and a source checkout or extracted paper bundle:

    python benchmarks/rac1/run.py --output rac1-report.json

To use verified files already downloaded from the source commit, pass their directory with `--data-dir`. The report records the source hashes and cell-level errors.

The numerical frame-average law agreed with its independent closed-form limit within 5.0 × 10⁻⁹ nm/s for the unit-amplitude example. The shared drive amplitude estimated from four calibration cells was 7.16 nm/s (individual amplitudes 4.05–11.60 nm/s). Across three reserved PA-Rac1 cells, pooled cell-weighted RMSE was 16.80 nm/s for the driven model and 16.96 nm/s for the stationary baseline. The driven model improved two of three cells. Their observed mean changes ranged from −0.31 to +7.37 nm/s. The three reserved light-only controls shifted by −1.28 to −0.16 nm/s. See [the reference report](reference.json) for individual values and the numeric check.

Portability here means that the released schema and solver can express and execute the external assay's forward model, map its signal to observed units and score against withheld cells. The small aggregate gain and the cell that worsened leave forecasting value unestablished outside this small cohort.
