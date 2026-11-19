# Forecasting protrusion from paired observations

Use `forecast-fit` to estimate short-horizon models on calibration preparations and `forecast-score` to evaluate the frozen coefficients on independent preparations. This workflow compares the information in a paired actin channel with both a stationary reference and a fitted movement-history reference.

## Inputs

Supply a `forecast-v1` TOML file. The [H1299 plan](benchmarks/h1299/forecast.toml) illustrates the fields. Set the frame interval in seconds, horizon in frames, number of consecutive baseline frames, availability delay of the motion label in frames, and motion unit. A label that requires one subsequent image needs at least a one-frame delay. Use disjoint lists of whole preparation IDs.

Each split has one CSV with header `preparation,unit,frame,actin,protrusion`. Units identify trajectories nested within preparations. Use original integer frame indices and the same cadence across a split. Positive actin values are required for the natural log. The protrusion measurement must be nonnegative. Leave missing values blank. No interpolation bridges gaps. A trajectory needs an initial usable sequence and at least one evaluable forecast. Every declared preparation must contribute an evaluable forecast.

## Models and scoring

For a forecast issued at frame t, each model uses the same mean protrusion from the first usable prefix. The prefix must have finished and all its motion labels must be available by t. The stationary prediction stays at this prefix mean. Delayed persistence uses the latest motion admitted by the configured delay. The motion-only model fits a coefficient multiplying the deviation of that delayed value from baseline. The actin-modulated model fits the same term plus a coefficient multiplying the current log actin signal minus its prefix mean. The last two models are fitted on calibration preparations only. Each horizon is fitted directly with fixed coefficients.

    protrusion-memory forecast-fit --plan forecast.toml --data calibration.csv --output fit.json
    protrusion-memory forecast-score --plan forecast.toml --data heldout.csv --fit fit.json --output score.json

Scoring checks the split, version and plan hash. It evaluates all models on identical origins. Each preparation gets equal weight, with equal weights across its trajectories and across eligible origins within each trajectory. The output includes counts at each level, RMSE in the user-declared motion unit, paired MSE differences by preparation and the range of the RMSE difference after leaving one preparation out. The range describes stability to omitting one preparation. Frames within a preparation remain dependent.

These are conditional forecasts refreshed with observed past motion and available actin at each origin. An improvement over the stationary model establishes temporal predictive value of the full predictor. The comparison against motion alone addresses the additional value of the measured actin signal. An observational predictor leaves the effect of a molecular perturbation unresolved.
