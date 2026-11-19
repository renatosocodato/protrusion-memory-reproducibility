# Calibration and held-out scoring

Supply two `ode-v1` models and a `calibration-v1` plan. The plan names the stationary reference and the modulated candidate. Their observation expressions may differ. Each scored channel has the same measurement unit and operator in both models. The software checks that every scored time is on the corresponding readout grid.

See [the complete CK666 plan](benchmarks/ck666/calibration.toml) for the file structure.

- `[split]` lists disjoint calibration and held-out preparation IDs. IDs must represent whole independent units for the scientific claim. Prefix IDs by experiment or source when labels are reused.
- `[scenarios.NAME]` declares known numeric overrides for model parameters. Scenario overrides cannot replace fitted values. Put piecewise interventions in the model TOMLs.
- `[channels.NAME]` declares the measurement `unit` and a fixed positive `scale` for scoring. Set the scale in the plan before scoring, using calibration information or a meaningful external scale.
- `[models.stationary]` and `[models.modulated]` give paths to the candidates. Add `[models.CANDIDATE.fit.PARAMETER]` tables with `initial`, `lower` and `upper` for selected model parameters. The stationary reference can also have fitted parameters. Fixed parameters stay at their declared values.

## Observation CSV

Use `preparation,unit,scenario,channel,time_s,value,measurement_unit` as the header. Each file contains exactly the preparation IDs declared for its split. Here `unit` identifies a sampling unit nested in its preparation, and `measurement_unit` is the physical or normalized readout unit. Each combination of preparation, sampling unit, scenario, channel and time must be unique. For deeper experimental hierarchies, aggregate lower levels at a scientifically defined sampling unit before export.

Use project-specific controls to define baseline inputs when appropriate. Apply the same baseline information to both candidates and keep these baseline records out of the scored target table unless they are an independently reserved outcome.

## Fit, then score

    protrusion-memory calibrate --plan plan.toml --data calibration.csv --output fit.json
    protrusion-memory score --plan plan.toml --fit fit.json --data heldout.csv --output score.json

The fit uses bounded least squares. Each preparation receives equal weight, with equal weights across sampling units within it and channels within each unit. Within each channel of a unit, observations receive equal weight. Residuals are divided by the declared channel scale before fitting or computing the pooled score. Channel-specific RMSE is in its native measurement unit and reweights the preparations and units supplying that channel.

The fit file records model and plan hashes, preparation lists, input hashes, parameter estimates and numerical diagnostics. Scoring checks the recorded software version, plan, model files and split. It uses the frozen estimates and outputs pooled and channel-specific RMSE, preparation-level results, and the leave-one-preparation-out range of the signed RMSE difference. A negative difference favours the modulated candidate on that endpoint. The range describes stability to omitting one preparation. Jacobian rank and active bounds are numerical diagnostics. Parameter identifiability requires an appropriate design and independent evidence.
