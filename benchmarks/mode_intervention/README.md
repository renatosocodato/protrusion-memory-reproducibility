# Actin intervention at fixed drive

The v1.8 analysis showed that the spatial spectrum distinguishes pathways only when angular gains are constrained. This follow-up reverses the signed actin perturbation in the same mode at fixed upstream drive. Active intervention is paired with matched delivery at each orientation. The gain matrices are estimated from separate preparations. In the ring model mode 2 reads the same two orthogonal directions for actin and deformation. The frozen window is 10–30 s and the acquisition cadence is an assigned 1 s.

In the representative held-out synthetic cohort the mediated mechanism changes the signed shape coefficient by **+0.251 µm**. The 95% interval widened for the assigned upstream and delivery contributions is **+0.199 to +0.302 µm**. For the direct-only mechanism, with actin perturbed but no actin-to-shape coupling, the corresponding interval is **−0.045 to +0.058 µm**, within the ±0.08 µm band. A delivery effect of 0.12 µm remains unresolved when supplied with a valid independent upper bound. A mismatched drive and failed actin first stage are reported separately.

In 100 fresh cohorts per case, the mediated effect was identified in 100/100 (95% Wilson Monte Carlo interval 96.3–100%). The direct-only and covered-leakage cases gave 0/100 declarations beyond leakage (0–3.7%). A hidden active-specific effect of 0.24 µm with an underestimated 0.025 µm bound produced positive declarations in 100/100 even though the true actin-to-shape coupling was zero. These results characterize the simulated settings. They show the role of independently constraining the mechanical effects of delivery and actuator activity. [Design](DESIGN.md) · [Effects](results.json) · [Recovery](recovery-results.json)

## Biological eligibility

G01 supplies the geometry and a one-vesicle observation scale. In the same GUV study, latrunculin A prevented growth of F-actin patches and reduced the fluorescence of patches already present. The accompanying public archive lists those experiments. The available material does not establish the four randomized mode-specific arms, independent angular standards, a delivery bound, or identities that would establish independent calibration and reserved preparations. It supports a related actin perturbation observation and motivates the proposed assay. [Original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11168455/) · [Image archive](https://doi.org/10.5061/dryad.tqjq2bw66)

## Reproduce

    python benchmarks/mode_intervention/challenge.py --datasets synthetic-mode-test
    python benchmarks/mode_intervention/recovery.py
    protrusion-memory mode-calibrate --plan synthetic-mode-test/mediated/plan.toml --data synthetic-mode-test/mediated/gains.csv --output gains.json
    protrusion-memory mode-test --plan synthetic-mode-test/mediated/plan.toml --data synthetic-mode-test/mediated/reserved.csv --calibration gains.json --output effect.json

See the [input requirements and effect definitions](../../MODE_INTERVENTION.md).

Use R with jsonlite and Cairo to render the effect and allowance intervals:

    Rscript benchmarks/mode_intervention/render.R benchmarks/mode_intervention/results.json signed_mode_summary.pdf
