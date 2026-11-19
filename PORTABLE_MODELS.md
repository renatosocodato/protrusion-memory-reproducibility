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

## Seeded stochastic trajectories with `sde-v1`

Use the `ode-v2` reporters, detector mixing and causal shutter with an Itô SDE. The existing `[equations]` give the drift for each state. Each `[diffusion.driver_name]` maps affected states to expressions giving the coefficient of that independent Wiener driver. Listing several states under the same driver gives them a common increment. Unlisted coefficients are zero. State and parameter names and expression operations follow the deterministic schema.

In `[stochastic]`, supply a nonnegative integer `seed` and positive `max_step_s`. Euler–Maruyama scales each increment by the square root of its actual duration and aligns the mesh with stimulus changes, the output grid and the shutter endpoints. The shutter accumulates the left-endpoint detector signal on that mesh. Check convergence by reducing the step. The explicit method is intended for resolved, nonstiff trajectories. It does not enforce positivity or a reflecting boundary.

The optional `measurement_noise = true` draws independent Gaussian error of each output's declared `noise_sd` after acquisition and any interval-rate operation. Thus the declared noise describes the final output units. Process and detector noise use independent streams. Enabled detector noise does not change the latent trajectory for the same seed and mesh.

    protrusion-memory simulate --config examples/stochastic-pulse.toml
    protrusion-memory ensemble --config examples/stochastic-pulse.toml --replicates 64

The ensemble command returns independent path seeds, detector means, sample standard deviations and pointwise quantiles. The deterministic comparison and calibration commands retain their existing contracts. The [stochastic pulse challenge](benchmarks/stochastic/README.md) supplies a mechanism-specific discrimination test with atomic conditions and stated noise assumptions.

## Linear periodic fields with `ring-v1`

The ring has a fixed angular coordinate, evenly spaced nodes and periodic boundary conditions. Declare `ring.nodes` and the Fourier `ring.modes`. Each field has a nonnegative `decay_s_inv` and `diffusion_rad2_per_s` (both optional, default zero), optional initial cosine and sine coefficients, local linear coupling to other fields, and linear drives. An input multiplies a declared temporal parameter by separate cosine and sine patterns. Coupling and drive strengths may be real numbers or declared parameter names. Parameter modulation follows the deterministic interface.

The solver diagonalizes the second-order centered periodic difference operator. Mode n sees the nonnegative value 4 sin²(πn/N)/(2π/N)² as the negative Laplacian eigenvalue. Each field retains its declared real modes, with only a cosine coefficient for zero. Nyquist and aliased modes are excluded. Reconstruct a field by summing the cosine and sine terms. Units for angular diffusivity differ from diffusivity in physical distance. A fixed physical radius can convert between them. Enlarging the node count changes the discrete eigenvalues towards the continuous values n².

Reporters select a `field`, `mode` and `component` (cosine by default, spelled `cos`, or `sin`). The existing first-order response, detector mixing, exposure and delay are available after projection. This assumes the optical processing is spatially homogeneous and linear on the selected modes. The output contains the unfiltered mode coefficients and the detector readouts.

    protrusion-memory simulate --config examples/spatial-pulse.toml

`compare` accepts two ring configurations with the same grid and reporter definitions. `ring-v1` describes linear fields on a fixed ring. It has no moving boundary, nonlinear intermode coupling or spatial stochastic forcing. The [G01 and spatial pulse comparisons](benchmarks/spatial/README.md) give the assumptions and scientific use.
