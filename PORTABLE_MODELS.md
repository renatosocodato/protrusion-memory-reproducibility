# Models from other laboratories

Use the `ode-v1` TOML schema for a finite system of deterministic ordinary differential equations. Each state has an initial value at time zero and exactly one derivative expression. Parameters are finite real constants unless modulated. Expressions can use the states, parameters, time `t` in seconds, numeric constants, `+ - * / **`, and the single-argument functions `abs exp log sqrt sin cos tanh`. The interpreter processes arithmetic expressions directly and does not load external code. Use unambiguous identifiers for states and parameters.

`[time]` declares `stop_s` and `sample_s`. The trajectory starts at zero and includes the final time even if the last interval is shorter. DOP853 uses fixed relative and absolute tolerances of 10⁻⁹ and 10⁻¹¹. The solver restarts at each modulation boundary with continuous states. A `[[modulations]]` table declares the parameter, half-open time window `start_s` and `stop_s`, and one of these operations on its declared baseline:

- `set` uses `value` during the window.
- `add` adds `value` to the baseline during the window.
- `scale` multiplies the baseline by `value` during the window.

Windows on the same parameter must be disjoint. Adjacent windows and overlapping windows on distinct parameters are supported.

Each `[observations.channel_name]` declares an `expression` and a `unit`. The default `point` operator evaluates the expression at the sample time with the parameter on the right side of any boundary. For an accumulated quantity such as position, `interval_rate` divides the difference of consecutive sample values by their time separation and labels it by the later time. Declare the final output unit (for example, nm/s) and an optional positive `noise_sd` in that unit. Units are labels. Users supply equations with consistent dimensions and choose the sampling interval.

    protrusion-memory simulate --config examples/rac1-edge.toml
    protrusion-memory compare --config examples/rac1-edge.toml --against examples/rac1-stationary.toml

Portable comparison requires matching channel names, units, operators, times and noise scales. It reports the RMSE and maximum absolute difference between expected observations and the RMSE divided by SD when declared. It uses the supplied clocks directly. It reports no probability of a mechanism. The original exact-pair commands and assumptions are preserved for configurations without a schema key.

The [independent Rac1 benchmark](benchmarks/rac1/README.md) provides a data alignment and scoring example. The `ode-v1` interface covers finite deterministic ODEs and the documented observation operators. Models with spatial state, stochastic forcing, or observation-history effects beyond interval rates need a further specification.

## Reporter and camera acquisition with `ode-v2`

The v2 schema retains the v1 time grid, equations, parameters and modulations. A named reporter defines an `expression` over model states, parameters and time. A positive `response_s` gives a first-order response towards the source expression, and requires an explicit `initial` value at time zero. Default `response_s = 0` is instantaneous.

Each `[observations.name]` supplies a `mix` table mapping reporter names to finite linear coefficients. Coefficients include gains or unit conversions where required, so every term must be expressible in the detector's declared output units. This is a forward detector model. Apply instrument backgrounds, clipping or nonlinear normalization only with a separate observation model. Optical mixing applies to signals, not to velocities or masks.

At reported time t, a detector with `delay_s = d` and `exposure_s = e` averages its mixed signal over [t−d−e, t−d]. With zero exposure it evaluates the point signal at t−d. Both settings default to zero. The time label is the availability or alignment time, so the shutter closes at t−d. Early samples without a complete window are omitted. The optional `interval_rate` differences two consecutive acquired readings and labels the result with the later time. Integration respects stimulation boundaries, and instantaneous terms use the right-hand parameter at a boundary.

    protrusion-memory simulate --config examples/pulse-acquisition.toml

`compare` checks matching channel names and acquisition clocks before comparing outputs. `calibrate` and `score` accept v2 models using fixed acquisition parameters in their configurations. Supply only measurement times supported by both rivals. A fitted optical or timing parameter needs a fixed common target set across values and an independent calibration design. `noise_sd` describes the final channel's scale. Filtered, averaged or differenced measurements can have correlated noise. DOP853 is explicit, so very fast reporter times or very stiff state kinetics may need a different solver.

See the [pulse and imaging benchmarks](benchmarks/acquisition/README.md).
