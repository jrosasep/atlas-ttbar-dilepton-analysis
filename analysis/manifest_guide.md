# Cómo leer manifest.json

Este archivo es la lista de entradas del análisis: no contiene eventos ni
histogramas. Se conserva sin comentarios porque JSON no admite comentarios.

| Campo | Qué significa |
|---|---|
| `records` | Registros de CERN Open Data de los que proceden los archivos. |
| `luminosity_pb` | Luminosidad usada para normalizar el Monte Carlo, en pb⁻¹. |
| `luminosity_source`, `luminosity_note` | Fuente y precisión del valor adoptado. |
| `nominal_groups` | Identificadores de las simulaciones incluidas en cada proceso. |
| `files` | Lista de archivos de colisiones reales y simulaciones. |
| `key` | Nombre del archivo ROOT; se busca en `analysis/data/`. |
| `group` | `data` para colisiones; otros nombres identifican procesos simulados. |
| `dsid` | Identificador de una muestra de simulación, no de un evento individual. |
| `url`, `size`, `checksum` | Enlace, tamaño en bytes y control Adler-32 del archivo. |
| `metadata` | Parámetros de la muestra, no medidas de cada partícula. |

Dentro de `metadata`, `crossSection_pb` es la sección eficaz en pb,
`genFiltEff` la eficiencia del filtro de generación, `kFactor` una corrección
de normalización y `sumOfWeights` la suma original de pesos de generación.
Se usan juntos para normalizar MC a la luminosidad adoptada. No hay que
normalizar cada archivo únicamente por los eventos que pasan nuestra selección.

`fetch_inputs.py` verifica tamaño y checksum. Después registra en
`provenance/verified_inputs.json` la huella SHA256 y la ubicación local de cada
ROOT. `analysis.py` guarda qué manifiesto y programas utilizó. Un manifiesto
verificado describe una descarga; no es prueba de que el archivo siga presente
en cualquier computador donde se clone el repositorio.

Fuentes: [colisiones reales](https://opendata.cern.ch/record/93934),
[simulación](https://opendata.cern.ch/record/93913) y
[documentación educativa](https://opendata.atlas.cern/docs/data/for_education/13TeV25_details).

Los textos en `provenance/source_records/` conservan las fuentes exactas de la
corrida guardada. Son registros de procedencia, no programas alternativos para
ejecutar. Permiten redibujar sus histogramas sin fingir que se calcularon con
el código modificado hoy. Las figuras nuevas registran también el código del
postprocesamiento en `results/plot_provenance.json`.
