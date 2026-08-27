# Data policy and manifest

Large ROOT files are never committed. Put local downloads under `data/raw/`;
the directory is ignored by Git.

## Published 2015-2018 measurement

The ATLAS article uses 140 fb^-1 of proton-proton data at 13 TeV, recorded in
2015-2018. The article's data-availability statement says that derived data are
available from the ATLAS Collaboration upon request. The full event-level
analysis input and all calibration variations are therefore not mirrored here.

## Recommended compact development sample

- Record: `atlas-93913`
- DOI: `10.7483/OPENDATA.ATLAS.NNF8.76IX`
- Years: 2015+2016
- Skim: `2J2LMET30`
- Nominal top sample used during development:
  `ODEO_FEB2025_v0_2J2LMET30_mc_410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.2J2LMET30.root`
- Size of that file: 699,911,723 bytes
- Adler-32: `8d5c9cdd`

This sample is designed for education and outreach. It is not a substitute for
the complete Run-2 data and analysis software used in the Nature article.

The ntuple was produced with the official `PhysLiteToOpenData` framework
(Zenodo DOI `10.5281/zenodo.15791091`). Its `jet_btag_quantile` stores the
continuous DL1dv01 working-point bin. The framework's own skim code uses
`jet_btag_quantile >= 4` for the 70% working point; this project uses
`>= 2` for the looser 85% bin required by the article's baseline selection.
The precise choice and any scale factors remain configuration-controlled and
are documented in every generated summary.

## Older 10 fb^-1 sample

Record `atlas-15003` contains 2016 data and simulation preselected with at
least two leptons. It is useful for teaching event selections, but the portal
explicitly marks it as educational and unsuitable for scientific publication.
