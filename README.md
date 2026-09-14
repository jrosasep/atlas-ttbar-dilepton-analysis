# Lab III — Distribuciones angulares eμ con ATLAS Open Data

Proyecto de José Ignacio Rosas Sepúlveda. Análisis educativo independiente:
no es un resultado oficial de la colaboración ATLAS.

## 1. Qué estudiamos

¿Cómo se comparan las distribuciones de separación del electrón y del muón en
colisiones reales, seleccionadas para enriquecer la producción top–antitop,
con la simulación del Modelo Estándar y los fondos disponibles?

Calculamos Δφ (separación alrededor del haz) y |Δη| (separación absoluta en
pseudorrapidez) a partir de los leptones reconstruidos.

## 2. Por dónde empezar

1. Abre [el gráfico principal](results/angulos_conteos_SR2b.png): puntos negros =
   colisiones; colores = simulaciones; panel inferior = cociente datos/MC.
2. Lee [los resultados y sus límites](results/RESULTS.md).
3. En [analysis.py](analysis/analysis.py), empieza por la función
   **select_events**: decide qué eventos entran. Sigue por **event_weights** y **observables**.
4. Consulta [las ecuaciones y cómo leer las figuras](results/METHOD.md).

No necesitas estudiar cinco programas a la vez. **run.py es la única entrada
que debes ejecutar.** Los otros archivos separan responsabilidades, no versiones.

## 3. Qué contiene cada carpeta

| Carpeta | Contenido | Para qué abrirla |
|---|---|---|
| analysis/ | Código y archivos técnicos del cálculo. | Aprender o modificar el análisis. |
| results/ | Figuras, tablas CSV, método y explicación de resultados. | Estudiar, presentar y comprobar números. |
| licenses/ | Avisos del proyecto previo y de la referencia ATLAS. | Revisar atribuciones antes de publicar. |

Dentro de **analysis/**:

- **manifest.json**: lista de los 63 archivos, grupos físicos, direcciones, tamaños
  y metadatos. No contiene eventos.
- **provenance/**: metadatos oficiales, comprobaciones de descarga y cálculo, y copia
  de la referencia C++. No son análisis alternativos.
- **results/angular_v2/**: registros de **la única corrida incluida**. **samples/**
  tiene un JSON por ROOT; **run.json** registra la ejecución y **summary.json** agrupa
  las muestras. Permiten dibujar sin ROOT. El nombre angular_v2 es un identificador
  interno; no significa que se conserven aquí versiones antiguas.
- **data/**: aparecerá si descargas o enlazas los ROOT. No se publica en GitHub.

Los JSON internos son el registro detallado del cálculo; las tablas de results/
son su presentación legible, no otra versión del análisis.

En la raíz, **requirements-lock.txt** lista dependencias. Los archivos ocultos de Git
excluyen datos/cachés y preservan los bytes utilizados para verificar el código.

## 4. Qué hace cada programa y en qué orden

| Etapa | Archivo dentro de analysis/ | Responsabilidad |
|---|---|---|
| 1. Descargar | fetch_inputs.py | Obtiene archivos públicos y verifica tamaño y checksum. |
| 2. Analizar | analysis.py | Lee bloques, comprueba entradas, selecciona eμ, calcula pesos y ángulos, llena histogramas. |
| 3. Representar | plots_and_summary.py | Suma procesos, dibuja datos frente a MC y exporta tablas. |
| Control auxiliar | audit_met_inputs.py | Comprueba consistencia del momento faltante redondeado. El análisis importa su función; no debes ejecutarlo aparte. |
| Pruebas | test_analysis.py | Comprueba selección, fórmulas y entradas con ejemplos conocidos. |

**Archivos → selección → pesos y ángulos → histogramas → comparación.**

run.py coordina estos programas y coloca las figuras y tablas en results/.
pytest.ini configura las pruebas.

## 5. Cómo ejecutar

Con Python 3.12 y un entorno virtual, desde la raíz:

~~~text
python -m pip install -r requirements-lock.txt
python run.py check
python run.py plot
~~~

Esto prueba funciones y regenera gráficos desde los resultados guardados:
**no descarga datos ni vuelve a analizar los eventos**. Sin argumentos,
run.py muestra ayuda. Las dependencias son las del entorno probado, no una
garantía de compatibilidad con cualquier sistema.

Solo para repetir desde los eventos:

~~~text
python run.py download
python run.py check --with-data
python run.py analyze
~~~

La descarga ocupa **14,97 GB (13,94 GiB)**. Si ya tienes los ROOT, no los
descargues otra vez: se puede enlazar analysis/data a su carpeta existente.
Este paquete no incluye enlaces específicos al computador del autor.
analyze también genera las figuras. Se leen bloques de 100 000 eventos;
se necesita espacio adicional y tiempo. No usar Python con -O, que desactiva
controles. Las huellas del código detectan incluso cambios de finales de línea.

## 6. De dónde vienen los datos

- [Colisiones reales: CERN 93934](https://opendata.cern.ch/record/93934).
- [Simulación: CERN 93913](https://opendata.cern.ch/record/93913).
- [Formato y unidades](https://opendata.atlas.cern/docs/data/for_education/13TeV25_details).
- [Metadatos](https://opendata.atlas.cern/docs/data/for_education/13TeV25_metadata).

Archivos educativos 2J2LMET30 de la entrega 2025; colisiones de 2015–2016 a
13 TeV, **no todo el Run 2**. Ya tienen filtros previos.
La luminosidad es el valor educativo redondeado de 36 fb⁻¹.
MC incluye ttbar, single_top, diboson, Zjets, Wjets, ttV y HWW.

## 7. Qué usamos realmente de ATLAS

La referencia implementada es
[TTbarDilepAnalysis.C, commit ff71d6b](https://github.com/atlas-outreach-data-tools/atlas-outreach-cpp-framework-13tev/blob/ff71d6ba82f2afd5a45d2e9b7f80bf915ceb1c84/Analysis/TTbarDilepAnalysis/TTbarDilepAnalysis.C).
Su copia está en analysis/provenance/official/.

| Criterio de referencia | Dónde se aplica |
|---|---|
| Trigger electrónico o muónico; leptones tight, aislados y asociados al trigger. | select_events. |
| pT > 25 GeV y aceptación en η; exclusión de la transición del calorímetro para electrones. | Selección de objetos. |
| Dos leptones de distinta especie y cargas opuestas. | Un electrón y un muón aceptados de cargas opuestas. |
| Jets con JVT y jet_btag_quantile >= 2; al menos dos jets b. | Región principal SR2b. SR1b es una extensión auxiliar anidada. |
| Pesos MC y factores de eficiencia. | event_weights, normalizado explícitamente con metadatos oficiales. |

**No ejecutamos el framework C++ ni comprobamos equivalencia evento a evento
con él.** Adaptamos criterios a Python sobre la estructura preliminar del proyecto.
El propio C++ advierte que es educativo, no una reproducción de resultados
experimentales publicables de ATLAS.

El [notebook stop recomendado en el correo de ATLAS](https://github.com/IoPapadopoulos/notebooks-collection-opendata/blob/fd47966c6efa51cf9b0bea7104c41e85a8c97468/13-TeV-examples/uproot_python/stop_analysis.ipynb)
se revisó como referencia relacionada y usa el mismo skim.
**No lo incorporamos como motor ejecutado**: sus cortes y objetivo difieren.
Su recomendación no significa aprobación de nuestro código.

Aquí se desarrollaron los histogramas angulares, la región auxiliar, controles
de calidad y duplicados, pruebas, tablas y comparación de formas con propagación
de covarianza al normalizar. No deben atribuirse como implementación oficial.

## 8. ¿Falta código?

La cadena **descarga → selección → histogramas → comparación** está implementada
y ejecutada. No falta un módulo para producir las figuras actuales.

Sí falta desarrollo para conclusiones más exigentes: fondos no prompt/falsos,
variaciones sistemáticas y validación más amplia de la selección. Los errores
actuales son estadísticos. La comparación se realiza a nivel reconstruido;
no se han corregido las distribuciones por los efectos de aceptación y resolución
del detector.

La región principal contiene 36 945 colisiones frente a 36 691,3 eventos MC
esperados. La proximidad de los conteos totales no basta para establecer el
acuerdo de las distribuciones angulares ni validar todos los fondos.
SR2b está contenida en SR1b; no son muestras independientes.

## Créditos y estado del proyecto

Resultados preliminares para Laboratorio III. El código fue desarrollado con
asistencia de IA; las pruebas automatizadas no sustituyen la revisión física
del autor. Las limitaciones del análisis se indican arriba.

Se agradece a ATLAS Collaboration la publicación de datos y ejemplos educativos.
Los [créditos, DOI y condiciones de uso](licenses/DATA_SOURCES.md) identifican las
fuentes exactas. El código adaptado se distribuye bajo EUPL-1.1; se conserva
el aviso MIT del código previo en licenses/. Los datos tienen su licencia propia.
No se incluyen ROOT, correos ni documentos personales.
