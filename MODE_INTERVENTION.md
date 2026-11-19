# Test a signed actin intervention in one deformation mode

This workflow analyzes a predeclared non-translational mode and time window under fixed upstream stimulation. Use a centroid-frame contour and fixed orientation. For mode m, rotating the added actin intervention by π/m reverses the projection on that mode. It does not reverse every mode. The same stimulation history and timing should be monitored across the four conditions.

## Inputs and commands

The `mode-test-v1` TOML plan specifies `mode`, a unit `axis = [cosine_weight, sine_weight]`, `window_s = [start, stop]`, and disjoint `[split].calibration` and `[split].reserved` preparation IDs. It also has `[criteria]` and `[resampling]` tables. The synthetic [generator](benchmarks/mode_intervention/challenge.py) writes an example plan and CSVs when run with `--datasets`.

The gain CSV columns are `preparation,mode,channel,known_cos,known_sin,read_cos,read_sin`. Each calibration preparation must independently span the two angular axes for each of `actin`, `shape` and `drive`. Known, baseline-corrected standards from separate calibrations define gains in both axes, including rotations or cross-sensitivity. Regression through zero and equal preparation weighting estimate the 2×2 measured-from-known gain. The frozen calibration records input and plan hashes and sufficient statistics by preparation.

The response CSV columns are `preparation,mode,arm,actin_cos,actin_sin,shape_cos,shape_sin,drive_cos,drive_sin`. Supply one precomputed window summary per preparation and arm after aligned, identical acquisition and baseline handling. Shape is expressed in µm after calibration. Record detector coordinates before applying the matrices. Each preparation contributes `active_plus`, `delivery_plus`, `active_minus` and `delivery_minus`. Use randomized order or allocation before forming summaries. Delivery is the inactive or sham intervention in the same orientation. Cells, sectors and frames stay within preparation.

    protrusion-memory mode-calibrate --plan plan.toml --data gains.csv --output frozen.json
    protrusion-memory mode-test --plan plan.toml --data reserved.csv --calibration frozen.json --output report.json

## Contrast and limits

For each sign, subtract matched delivery from active intervention. Define the signed response as half the difference of the two shape contrasts projected on the frozen axis. Require both actin signs to respond by a predeclared amount (`actin_min`). Outside the declared drive tolerance (`drive_max`, a norm for each arm), the fixed-drive comparison remains unresolved.

`drive_gain` is the independent upper bound on upstream-to-shape susceptibility for the window. `delivery_bound` is an independent upper bound on the orientation-odd active-specific mechanical effect of the intervention. With residual drive vectors ΔU+ and ΔU−, triangle inequality bounds the upstream contribution by drive_gain × (‖ΔU+‖+‖ΔU−‖)/2. Add the delivery bound to this allowance. Optional `drive_error` is an independent upper bound on the measurement error norm of each active-minus-delivery drive contrast. This amount is added to each drive norm in the guard and to its mean in the leakage allowance. With the default zero, the calculation conditions on the observed drive summaries and propagates angular gain uncertainty. This bound needs a measurable matched control that also constrains direct effects specific to the active intervention. The standard delivery arm removes shared effects but cannot by itself bound the active-specific effect.

`resampling.seed` and `resampling.draws` set independent preparation bootstrap draws in the two cohorts. The report gives the signed shape result, a bootstrap interval expanded by the leakage allowance, and preparation-level responses including the orthogonal shape component. A shape interval entirely beyond ±`shape_min` supports a contribution under the specified bounds. An interval entirely inside that band bounds the response as small over the tested window. Overlap leaves the response unresolved. This interval is conditional on the upstream and delivery bounds transferring from separate calibrations. Supply conservative margins from those calibrations. Small preparation counts can limit bootstrap coverage.

The software tests the declared mode and summary window. Verify the measured upstream waveform, other modes and possible nonlinear coupling in calibration before using the result to identify a pathway. Unknown optical drift, timing errors and unbounded mechanical effects can yield the same result under a different mechanism. See the [synthetic challenge and biological data qualification](benchmarks/mode_intervention/README.md).
