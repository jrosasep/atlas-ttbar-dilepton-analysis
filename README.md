# ATLAS dilepton angular distributions
## Distribuciones angulares dileptónicas con ATLAS Open Data

Independent educational analysis developed for **Laboratorio III, Universidad de Concepción**.
Análisis educativo independiente desarrollado para **Laboratorio III, Universidad de Concepción**.

### Abstract

This project compares the reconstructed electron-muon angular distributions in
ATLAS collision data with a Standard Model Monte Carlo prediction. The selected
sample is enriched in dileptonic top-antitop production. The observables are the
azimuthal separation, Δφ(e,μ), and the absolute pseudorapidity separation,
|Δη(e,μ)|.

### Resumen

Este proyecto compara las distribuciones angulares reconstruidas del electrón y
del muón en colisiones registradas por ATLAS con una predicción Monte Carlo del
Modelo Estándar. La muestra está enriquecida en producción dileptónica
top-antitop. Los observables son la separación azimutal Δφ(e,μ) y la separación
absoluta en pseudorrapidez |Δη(e,μ)|.

## Scientific question / Pregunta científica

> To what extent are the observed Δφ(e,μ) and |Δη(e,μ)| distributions described
> by the simulated top-antitop signal and the available background processes?

> ¿En qué medida las distribuciones observadas de Δφ(e,μ) y |Δη(e,μ)| son
> descritas por la señal top-antitop simulada y los procesos de fondo disponibles?

## Data and method / Datos y método

| Component / Componente | Definition / Definición |
|---|---|
| Collision data / Datos | ATLAS proton-proton collisions, 2015-2016, √s = 13 TeV |
| Simulation / Simulación | top-antitop, single-top, diboson, Z+jets, W+jets, ttV and H→WW |
| Input format / Formato | Educational ROOT ntuples, `2J2LMET30` skim |
| Main selection / Selección | One accepted electron and one accepted muon with opposite charge and pT > 25 GeV; at least two accepted jets and at least two b-tags |
| Observables | Δφ(e,μ) and |Δη(e,μ)| |
| Comparison / Comparación | Event yields, normalized shapes and data/MC ratios |

Monte Carlo events are weighted using the integrated luminosity, cross sections,
generator weights and the available correction factors. The analysis is performed
with reconstructed objects and reads the ROOT files in chunks.

## Preliminary result / Resultado preliminar

For the main selection, the current analysis finds:

| Collision data / Datos | Total MC / MC total | MC statistical uncertainty / Incertidumbre estadística MC | Predicted tt̄ fraction / Fracción tt̄ predicha |
|---:|---:|---:|---:|
| 36,945 | 36,691.3 | 131.9 | 96.57% |

These numbers describe the selected sample as modelled by the available
simulations. They do not identify the physical origin of individual collision
events. The principal angular plots are available in
[results/](results/), together with the [method](results/METHOD.md) and
[numerical summary](results/RESULTS.md).

Estos valores describen la muestra según las simulaciones disponibles, pero no
identifican el origen físico de cada colisión real.

## Reproduce / Reproducir

Python 3.12 is recommended. From the repository root:

~~~text
python -m pip install -r requirements-lock.txt
python run.py check
python run.py plot
~~~

These commands run tests that do not require local ROOT files and regenerate the
figures from the stored numerical results.

To download and process the complete input used here (63 files, approximately
14.97 GB):

~~~text
python run.py download
python run.py check --with-data
python run.py analyze
~~~

## Repository structure / Estructura

~~~text
analysis/   event selection, weights, observables, tests and provenance
results/    figures, tables and explanation of the method
licenses/   software notices and data attribution
run.py      single command-line entry point
~~~

The processing sequence is:

~~~text
ATLAS ROOT files → event selection → weights and observables
                 → histograms → data/MC comparison
~~~

## Current scope / Alcance actual

The current uncertainty bands are statistical. A complete treatment of
systematic uncertainties and non-prompt or misidentified-lepton backgrounds is
still pending. The distributions have not been corrected for detector acceptance
and resolution. Therefore, this is a reconstructed-level educational comparison,
not a precision measurement.

Las bandas actuales representan incertidumbres estadísticas. Todavía falta un
tratamiento completo de las incertidumbres sistemáticas y de los fondos con
leptones no prompt o mal identificados. Las distribuciones no han sido corregidas
por aceptación y resolución del detector.

## Data, references and credit / Datos, referencias y créditos

- [ATLAS collision data, DOI 10.7483/OPENDATA.ATLAS.0CJR.N7ZT](https://doi.org/10.7483/OPENDATA.ATLAS.0CJR.N7ZT)
- [ATLAS Monte Carlo simulation, DOI 10.7483/OPENDATA.ATLAS.NNF8.76IX](https://doi.org/10.7483/OPENDATA.ATLAS.NNF8.76IX)
- [ATLAS Open Data documentation](https://opendata.atlas.cern/docs/data/for_education/13TeV25_details)
- [Educational TTbarDilepAnalysis reference](https://github.com/atlas-outreach-data-tools/atlas-outreach-cpp-framework-13tev/tree/ff71d6ba82f2afd5a45d2e9b7f80bf915ceb1c84/Analysis/TTbarDilepAnalysis)

The original ROOT files are hosted by CERN and are not duplicated in this
repository. Their exact paths and checksums are recorded in
[analysis/manifest.json](analysis/manifest.json). ATLAS Open Data are released
under CC0; citation and acknowledgement of the ATLAS Collaboration are requested.
Neither ATLAS nor CERN endorses this analysis.

The preliminary implementation was developed with AI assistance and requires
continued scientific review by the author. Software and data notices are listed
in [licenses/](licenses/).
