# Centroid motion, membrane deformation and pathway resolution

## Geometry of the public G01 vesicle

The [published G01 inclusion mask](https://github.com/basharif/GUV-symmetry-breaking-analysis/blob/56d7a4eb5cd3c2097c4e08ab138734dcf1e6edd3/G01_input/G01_acs_v0/G01_mask_inclusion.tif) is derived from membrane fluorescence. Its frame interval is 65.589 s, and 28 frames span 29.5 min. The image-plane area centroid gives a reference for orienting each radial profile. Radial coefficients use the same laboratory directions in all frames, with translation removed and rotation retained. The area-equivalent radius, the first harmonic and the cos/sin pairs of modes 2–4 are recorded separately. Actin is sampled in a narrow band on the inner side of the membrane contour. Association between actin intensity and shape can still depend on the imaging geometry and the optical channels.

All 28 masks yield a centroid. Of these, 27 yield complete current-centroid profiles and 25 yield complete profiles in both the laboratory and current-centroid coordinates. Profiles with multiple ray crossings are excluded from radial analysis. In the 25 paired frames, the median first-mode amplitude drops from 1.338 to 0.016 µm when the reference is moved to the current centroid. The largest change in the centered mode-2 cosine/sine vector relative to its early mean is 3.210 µm. The centroid excursion reaches 12.131 µm, and the projected area-equivalent radius changes from 18.562 to 9.969 µm. This is a change in the image-plane section and does not establish a change in vesicle volume.

A known subpixel shift of the first mask yields 0.155 µm centroid error and 0.0034 µm change in the centered mode-2 vector. An artificial local bulge on a circle with a fixed geometric reference shifts its area centroid by 0.628 µm. These controls illustrate why image-plane centroid displacement can include the geometric effect of asymmetric deformation. Subpixel-shift leakage is not a general segmentation uncertainty. [Geometry and controls](g01-geometry.json) · [Data archive](https://doi.org/10.5061/dryad.tqjq2bw66)

## Spatial pathway challenge

A single five-second drive excites modes 2–4 in a fixed periodic ring. The direct and actin-mediated alternatives share the shape operator and are matched for mode-2 integrated deformation. An extra actin filter gives the mediated/direct integrated-response ratio `(k_a + D_a λ_2)/(k_a + D_a λ_n)` in mode n, with the discrete Laplacian eigenvalue λ. With assigned actin diffusivity of 0.015 rad²/s the ratios for modes 2–4 are 1.00, 0.71 and 0.51. At zero actin diffusivity this spatial difference disappears. Separate gains per mode can absorb the ratios, leaving temporal information to distinguish the pathways.

In the recovery test, 4,000 synthetic draws per class and setting use white noise scaled to the variation of fitted radial-mode coefficients in five usable early G01 frames (0.098, 0.096 and 0.081 µm, respectively). We freeze detection using separate null realizations. A decision also requires a prespecified likelihood margin. Unresolved outcomes remain in the denominators. The rate below refers to correct assignments, not to accuracy conditional on having declared an assignment.

| G01 pulse phase within a frame | Mediated, modes 2–4 and common gain | Direct, modes 2–4 and common gain |
| --- | ---: | ---: |
| 0.00 | 0% | 0% |
| 0.25 | <0.1% | <0.1% |
| 0.50 | 7% | 4% |
| 0.75 | 53% | 56% |
| 0.99 | 0% | 0% |

The best regular phase detected all 4,000 responses per alternative before the mechanism margin. At that phase, 2.5% of mediated and 1.9% of direct cases received the wrong label. Mode 2 alone or separate gains per mode resulted in zero correct declarations at each of the five regular phases. The separate near-frame alignment yielded some recovery through timing. At the assigned 1-s cadence all cases were correctly declared in this simulation. That is a model-conditional result at known pulse timing, kinetic laws and optical settings, and under assigned noise. A single G01 recording at its actual clock provides no independent test of angular gains or pathway identity. [Fixed design](DESIGN.md) · [Confusion tables by phase](results.json)

## Reproduce

Requires Python, the package and Pillow for the TIFF extractor. Download the pinned inclusion mask and the matching [raw actin TIFF](https://github.com/basharif/GUV-symmetry-breaking-analysis/blob/56d7a4eb5cd3c2097c4e08ab138734dcf1e6edd3/G01_input/G01_raw/G01_YFP_actin.tif).

    python benchmarks/spatial/guv_geometry.py G01_mask_inclusion.tif G01_YFP_actin.tif > g01-geometry.json
    python benchmarks/spatial/discriminate.py g01-geometry.json

Use R with jsonlite and Cairo to render the three-panel vector summary:

    Rscript benchmarks/spatial/render.R g01-geometry.json benchmarks/spatial/results.json spatial_summary.pdf
