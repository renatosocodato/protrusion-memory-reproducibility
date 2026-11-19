# Five-minute protrusion forecasts in H1299 cells

The source [publication](https://doi.org/10.1371/journal.pone.0135204) reports paired whole-cell RubyRed-LifeAct imaging and membrane motion at five-minute intervals. The [original workbook](https://doi.org/10.1371/journal.pone.0135204.s001) has 6,419 source observations from 122 cells across 19 independently acquired dates. Protrusion is the area of the segmented mask newly occupied over one interval, in the lab frame. This endpoint can also reflect translation or segmentation variation. LifeAct is a whole-cell intensity. It does not resolve local assembly at the edge.

The fixed split used the first 12 dates for calibration and the latest seven for scoring. All 122 cells contributed eligible origins, with 90 calibration cells and 32 held-out cells. The 1,131 scoring origins are repeated forecasts within those held-out cells. The primary contrasts use the same targets and inputs available up to the forecast origin.

| Predictor | Held-out RMSE (µm²) |
| --- | ---: |
| Stationary prefix mean | 47.759 |
| Delayed persistence | 54.665 |
| Fitted motion history | 45.954 |
| Fitted motion history and LifeAct | 46.066 |

The motion-plus-LifeAct forecast reduced RMSE by 1.693 µm² (3.5%) relative to the stationary mean and lowered squared error on five of seven dates. The RMSE difference ranged from −2.375 to −1.219 µm² after omitting each date in turn. Adding LifeAct increased RMSE by 0.112 µm² relative to fitted motion history, with lower squared error on three of seven dates. The leave-one-date-out difference ranged from −0.015 to +0.207 µm². These ranges are descriptive.

Motion history supplies useful short-horizon information in these cells. This whole-cell actin readout adds no consistent predictive gain under the fixed design. The source series are observational and were processed retrospectively. A localized actin measurement, an independent intervention or a different cadence would answer a separate question.

## Reproduce

Install the software and the workbook adapter dependency `pyxlsb==1.0.10`. Run from the checkout or extracted paper bundle.

    python benchmarks/h1299/prepare.py --split calibration --output h1299-calibration.csv
    protrusion-memory forecast-fit --plan benchmarks/h1299/forecast.toml --data h1299-calibration.csv --output h1299-fit.json
    python benchmarks/h1299/prepare.py --split heldout --output h1299-heldout.csv
    protrusion-memory forecast-score --plan benchmarks/h1299/forecast.toml --data h1299-heldout.csv --fit h1299-fit.json --output h1299-score.json

The adapter checks the source workbook hash and the expected date, cell and frame counts. Use `--source` with the original workbook for offline export. [Fixed design](DESIGN.md) · [Forecast plan](forecast.toml) · [Reference results](reference.json)
