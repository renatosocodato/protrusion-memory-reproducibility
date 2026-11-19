# Imaging acquisition and actin-to-motion discrimination

The follow-up project asks whether local actin contains information about later edge motion and whether an independent perturbation can distinguish actin-mediated motion from a direct mechanical response. The observation laws are defined before comparing those mechanisms.

## Equal-area pulse challenge

[The model challenge](challenge.py) uses a five-second pulse and assigns actin and deformation relaxation times of 8 and 20 s. The actin channel has a 4-s response and a 0.1 actuator-signal contribution. Both channels use a 0.7-s exposure and a 0.6-s delay. A response mediated by actin and a direct response can produce the same integrated deformation. The test asks whether their time courses still separate after acquisition.

At 1 s sampling, the two acquired shape responses gave a ratio of the smaller to larger singular value of 0.255. The noiseless, correctly specified acquisition model recovered both components. The instantaneous fit assigned a spurious normalized direct coefficient of −0.117 to the mediated-only case. This compares receiver models with assigned, known optics; no empirical calibration is implied. The common actin traces constrain the input response but contribute no information to this particular split between the two mechanical branches.

At the 65.589-s G01 clock, the ratio across five assigned pulse phases ranged from 0.00016 to 0.01163. A separate alignment at onset 130 s landed near a recorded frame and gave 0.435. The unresolved quantity is the phase of a short pulse relative to the effective shutter and the noise of the resulting shape samples. Inverse recovery at one convenient phase gives no guarantee across exposures. [Full numerical results](synthetic-results.json)

## Replay of two existing biological series

- **G01 vesicle.** The [paired raw TIFFs](../edge/README.md) yield local actin, mCherry membrane intensity and signed radial membrane motion. We apply the observation law forward to recorded fluorescence, represented as constant between consecutive frames. Four initial frames set a separate scale in each sector and channel. We compare identical origins after frame 7 only. The 146 common sector-frame origins are repeated measurements from eight sectors of one vesicle. Median within-sector correlation of actin with next-frame radial velocity changes from −0.369 with the recorded point samples to −0.135 with an extra one-frame lag or exposure. Adding 0.1 of the normalized membrane intensity changes the median to −0.363. Applying all four operators together gives −0.268. These are assigned perturbations to the readout, with the original extracted geometry held fixed. They do not estimate the optics of the actual instrument or independent preparation-level prediction. [Results by sector](guv-results.json)
- **H1299 cells.** The existing [date-held-out dataset](../h1299/README.md) provides 32 cells, seven reserved dates and 1,067 common forecast origins after acquisition history is required. The readout is whole-cell LifeAct and protruded area, with five-minute frames. Against fitted movement history, adding LifeAct changes RMSE by +0.124 µm² for matched point samples and −0.045 µm² after a one-frame reporter response, trailing exposure and delay. The latter difference changes sign when an individual date is omitted (range −0.136 to +0.056 µm²). Thus a change in observation timing is enough to reverse the pooled direction of a small, unstable effect in this retrospective analysis. No second fluorescence channel was available for optical mixing in these exported trajectories. [Results by date](h1299-results.json)

## Reproduce

Install the package and run from the checkout or extracted paper bundle. The G01 extractor additionally uses Pillow and checks the hashes of both public TIFFs. Use the existing H1299 adapter to export calibration and held-out CSVs.

    python benchmarks/acquisition/challenge.py
    PYTHONPATH=benchmarks/edge python benchmarks/acquisition/guv.py G01_mCh_memb.tif G01_YFP_actin.tif
    python benchmarks/acquisition/h1299.py benchmarks/h1299/forecast.toml h1299-calibration.csv h1299-heldout.csv

The next discriminating biological test needs locally paired actin and independent contour channels, known pulse and shutter times, an intervention-only control, and whole independent preparations held for scoring. Use the simulator to select a cadence and pulse alignment from calibration measurements before opening those preparations.
