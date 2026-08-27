# ATLAS Run-2 top-antitop entanglement

[![tests](https://github.com/jrosasep/atlas-ttbar-entanglement-run2/actions/workflows/tests.yml/badge.svg)](https://github.com/jrosasep/atlas-ttbar-entanglement-run2/actions/workflows/tests.yml)

A transparent, reproducible study of the entanglement marker measured by the
ATLAS Collaboration in dileptonic top-antitop events at 13 TeV. The project
has three deliberately separate layers:

1. a **numerical reproduction** of the published Run-2 result and its quoted
   uncertainty budget;
2. a **toy Monte Carlo** implementation of the angular distribution and the
   estimator used by ATLAS, including closure tests;
3. an **ATLAS Open Data extension** for studying what can and cannot be
   reconstructed with the public 2015-2016 samples.

This repository does **not** claim to reproduce the full ATLAS measurement.
The publication uses 140 fb^-1 collected in 2015-2018, detector calibrations,
background estimates, systematic variations, and simulation-derived
calibration curves that are not all public. The paper states that derived data
and analysis configuration are available from ATLAS upon request.

## Physics target

ATLAS uses the angle between the charged leptons in the rest frames of their
parent top quarks:

```text
D = -3 <cos(phi)>
(1/sigma) d sigma / d cos(phi) = (1/2) [1 - D cos(phi)].
```

At parton level, `D < -1/3` is a sufficient entanglement condition. Because
the published result is fiducial and particle-level, ATLAS first maps this
parton-level boundary to particle level. For Powheg+Pythia the paper quotes
`D_limit = -0.322 +/- 0.009`; it is therefore incorrect to compare the measured
particle-level number directly with `-1/3` without this mapping.

The published signal-region result is

```text
340 < m_ttbar < 380 GeV
D_observed = -0.537 +/- 0.002 (stat.) +/- 0.019 (syst.)
D_expected = -0.470 +/- 0.002 (stat.) +/- 0.017 (syst.)
```

![Published ATLAS particle-level values](results/figures/published_regions.png)

## Reproduce the current results

On Windows PowerShell:

```powershell
py -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
.venv\Scripts\python scripts/reproduce_published_result.py
.venv\Scripts\python scripts/run_toy_closure.py
.venv\Scripts\python -m pytest
```

The scripts write machine-readable summaries and figures under `results/`.
All random studies use explicit seeds.

![Moment-estimator closure](results/figures/toy_closure.png)

## Open-data status

The most directly useful compact public sample identified so far is the 2025
ATLAS 2015+2016 `2J2LMET30` beta release: at least two jets, at least two tight
leptons, and missing transverse momentum above 30 GeV. It is enriched in
dileptonic top events and includes a nominal Powheg+Pythia `ttbar` simulation.

The raw ROOT files are intentionally excluded from Git. After downloading a
sample into `data/raw/`, inspect it with:

```powershell
.venv\Scripts\python scripts/inspect_open_data.py data/raw/sample.root
```

Whether the exact ATLAS observable can be formed depends on the available
truth history or on implementing a two-neutrino top reconstruction. The
published analysis primarily uses the analytic Ellipse method, falls back to
Neutrino Weighting, and then applies detector-to-particle and
parton-to-particle calibration curves. A detector-level proxy is not labelled
as the published measurement.

The inspected nominal `ttbar` file has 1,388,001 events and 118 branches. It
contains reconstructed lepton and jet four-vectors, missing transverse
momentum, event weights, and simplified truth-object collections. It does not
contain parent-top four-vectors, decay ancestry, or the two neutrinos
separately. The exact `cos(phi)` construction therefore requires a dileptonic
top reconstruction; it cannot be read directly from this flat ntuple.

The reproducible first open-data milestone is an `e mu` cutflow plus lab-frame
spin-correlation controls:

```powershell
.venv\Scripts\python scripts/analyse_open_data_baseline.py `
  data/raw/ttbar_nonallhad_2J2LMET30.root
```

The resulting dilepton azimuthal separation, opening angle, and invariant mass
are explicitly labelled as control observables, not as the ATLAS marker `D`.

![Selected open-data ttbar controls](results/figures/open_data_ttbar_baseline.png)

A development-only numerical solver for the two-neutrino constraints is also
included. It tries both lepton-jet assignments, enforces the two W and two top
mass constraints together with measured missing transverse momentum, and
chooses the real solution with the lowest reconstructed `m_ttbar`. Run a small
prototype sample with:

```powershell
.venv\Scripts\python scripts/reconstruct_open_data_proxy.py `
  data/raw/ttbar_nonallhad_2J2LMET30.root --max-selected 500
```

Its `D_detector_proxy` output remains a development quantity until the solver
is validated against generator history or the ATLAS Ellipse implementation and
the detector-to-particle calibration is available.

![Numerical dilepton reconstruction prototype](results/figures/reconstruction_proxy.png)

## Scope table

| Layer | Input | What is reproduced | Scientific status |
|---|---|---|---|
| Published numbers | ATLAS paper | Regions, `D`, uncertainty budget, folded boundary | Exact transcription and independent arithmetic checks |
| Toy Monte Carlo | Analytic angular PDF | Estimator, statistical scaling, reweighting, closure | Pedagogical validation |
| Open Data baseline | ATLAS 2015-2016 public ROOT files | `e mu` selection, cutflow, dilepton controls | Educational extension; not the 140 fb^-1 result |
| Open Data reconstruction | Same ntuples plus a two-neutrino solver | Detector-level `m_ttbar` and `cos(phi)` proxy | Planned; requires closure and calibration studies |

## Primary references

- ATLAS Collaboration, *Observation of quantum entanglement with top quarks at
  the ATLAS detector*, Nature **633** (2024) 542-547,
  [arXiv:2311.07288](https://arxiv.org/abs/2311.07288),
  [DOI:10.1038/s41586-024-07824-z](https://doi.org/10.1038/s41586-024-07824-z).
- [ATLAS public auxiliary page, TOPQ-2021-24](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2021-24/).
- [ATLAS 2015+2016 `2J2LMET30` Open Data](https://opendata.cern.ch/record/atlas-93913),
  DOI:10.7483/OPENDATA.ATLAS.NNF8.76IX.
- [ATLAS 2016 two-lepton educational release](https://opendata.cern.ch/record/atlas-15003),
  DOI:10.7483/OPENDATA.ATLAS.GQ1W.I9VI.

## Authorship and acknowledgements

Original analysis code and documentation: Jose Ignacio Rosas. The collision
data, simulated samples, published values, and experiment-specific methods are
the work of the ATLAS Collaboration. This independent educational project is
not reviewed or endorsed by ATLAS or CERN.
