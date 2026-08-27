# Lista de revisión antes de usar académicamente el proyecto

Este documento no certifica calidad científica. Es una lista de preguntas que
José debe poder responder con sus propias palabras y comprobar en el código.

## 1. Pregunta física

- [ ] Explicar qué significa entrelazamiento de espines en un par top-antitop.
- [ ] Explicar por qué los leptones cargados actúan como analizadores de espín.
- [ ] Derivar o justificar `D = -3 <cos(phi)>` desde la distribución angular
      utilizada.
- [ ] Distinguir nivel partónico, nivel de partículas y nivel de detector.
- [ ] Explicar por qué la frontera partónica `-1/3` no se compara directamente
      con el valor fiducial publicado.
- [ ] Leer el artículo completo y el material auxiliar TOPQ-2021-24.

## 2. Cifras publicadas y estadística

- [ ] Verificar cada número de `config/published_atlas.yaml` contra su fuente.
- [ ] Revisar la combinación por cuadratura de incertidumbres sistemáticas.
- [ ] Entender por qué `z_simple` no es la significancia oficial de ATLAS.
- [ ] Identificar correlaciones, parámetros de nuisance y supuestos omitidos.
- [ ] Decidir si el control gaussiano aporta valor o debe eliminarse.

## 3. Monte Carlo sintético

- [ ] Revisar la normalización y el dominio físico de la densidad angular.
- [ ] Entender el muestreo aleatorio y la fórmula de incertidumbre del estimador.
- [ ] Reproducir las pruebas de cierre con distintas semillas y tamaños.
- [ ] Explicar qué demuestra la reponderación y, sobre todo, qué no demuestra.

## 4. Muestra pública y selección

- [ ] Confirmar la página, DOI, campaña y generador de la muestra.
- [ ] Poder explicar que la ejecución incluida usa simulación, no datos reales.
- [ ] Revisar las 118 ramas y justificar cada variable empleada.
- [ ] Comprobar unidades, cortes, cargas, multiplicidades y definición de b-tag.
- [ ] Comparar el flujo de selección con una referencia física apropiada.
- [ ] Revisar pesos y normalización; el flujo actual usa conteos sin pesar.

## 5. Reconstrucción

- [ ] Derivar las seis restricciones cinemáticas del sistema dileptónico.
- [ ] Revisar ambos emparejamientos leptón-jet y el criterio de elección.
- [ ] Diagnosticar por qué solo el 43.4% de la prueba pequeña obtiene solución.
- [ ] Validar tops y neutrinos reconstruidos contra verdad de generador.
- [ ] Comparar con Ellipse y Neutrino Weighting descritos por ATLAS.
- [ ] Estudiar resolución, sesgo, backgrounds y sistemáticas.
- [ ] No interpretar `D_detector_proxy` hasta completar estas validaciones.

## 6. Reproducibilidad y autoría

- [ ] Crear el entorno desde cero y ejecutar todos los scripts.
- [ ] Confirmar que las pruebas automáticas pasan.
- [ ] Revisar manualmente cada figura y tabla.
- [ ] Reescribir o comentar las partes que José decida adoptar.
- [ ] Mantener la declaración de asistencia de IA completa y exacta.
- [ ] Separar claramente las presentaciones históricas de este trabajo posterior.

## 7. Revisión académica

- [ ] Preparar una explicación oral de cinco minutos sin apoyarse en el README.
- [ ] Pedir a Francisca Garay una revisión del objetivo y del sentido físico.
- [ ] Pedir a Guillermo Rubilar una revisión metodológica si corresponde.
- [ ] Corregir todas las observaciones antes de hacerlo público.
- [ ] Solo después decidir si el proyecto merece una línea prudente en el CV.
