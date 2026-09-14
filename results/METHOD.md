# Método y lectura de los resultados

## 1. De los eventos al histograma

Un evento es un registro de una colisión o de su simulación. Las ramas de leptones
y jets contienen objetos reconstruidos, no quarks libres observados directamente.
Una etiqueta b identifica un jet compatible con contener productos de hadrones
con quark bottom; no es una identificación infalible.

`select_events` construye máscaras booleanas: cada requisito indica qué eventos
pasan. `cutflow.csv` registra el conteo tras cada filtro. No se comparan eventos
de simulación individualmente con colisiones individuales: se comparan distribuciones.

## 2. Qué ángulos calculamos

El azimut φ describe la dirección alrededor del haz. Para dos leptones:

$$\Delta\phi=\arccos[\cos(\phi_e-\phi_\mu)]\in[0,\pi].$$

Esta expresión toma la separación menor y respeta la periodicidad. La
pseudorrapidez se define por $\eta=-\ln\tan(\theta/2)$, donde θ es el ángulo con
el haz. Usamos las ramas η reconstruidas y calculamos:

$$|\Delta\eta|=|\eta_e-\eta_\mu|.$$

Ambas operaciones están en `observables`. No son el ángulo entre leptones en los
reposos de sus tops: no pueden sustituirse sin más en una fórmula del testigo D.
Las pruebas incluyen periodicidad, valores conocidos y permutación de leptones.

## 3. Por qué se pondera la simulación

Los datos reales cuentan con peso uno. El número de eventos simulados depende
del tamaño de la producción informática, no de cuántas colisiones observó ATLAS.
`event_weights` convierte cada evento MC en una contribución esperada:

$$w_i=\frac{\mathcal L\,\sigma\,k\,\epsilon_{\rm filtro}}
{\sum_{\rm producción}w_{\rm generador}}
w_{{\rm generador},i}\prod_a SF_{a,i}.$$

Aquí $\mathcal L=36000\,\mathrm{pb}^{-1}$ es la luminosidad educativa redondeada,
σ la sección eficaz en pb, k el factor de corrección de normalización,
ε la eficiencia del filtro del generador, y SF los factores de corrección
disponibles para reconstrucción/selección. La suma del denominador corresponde
a la producción original, no a los eventos que sobreviven nuestros cortes.
Se conservan los pesos negativos. No se ajusta la normalización para forzar acuerdo.

En cada intervalo B, `histogram` acumula:

$$N_B=\sum_{i\in B}1,\qquad Y_B=\sum_{i\in B}w_i,
\qquad V_B=\sum_{i\in B}w_i^2.$$

Para contribuciones estadísticas independientes, el error estimado es
$\sqrt{V_B}$. Con peso uno recuperamos $\sqrt{N_B}$, una aproximación de conteo
que requiere cautela en intervalos poco poblados. `histograms.csv` permite
comprobar cada barra sin abrir el código de dibujo.

## 4. Cómo leer las figuras

- **Conteos:** eje horizontal = ángulo; vertical = eventos por intervalo.
  Puntos negros = colisiones; componentes de color = predicción MC por proceso.
  La banda muestra incertidumbre estadística de MC, no todas las incertidumbres.
- **Cociente inferior:** datos/MC. Uno significa igualdad de valores centrales,
  no una prueba de que el modelo sea correcto. Las barras de datos y la banda de
  MC se presentan por separado.
- **Formas:** cada distribución se divide por su propio total. Se compara la
  forma perdiendo información sobre el rendimiento total; no se restan fondos.

Para fracciones $p_i=Y_i/S$, $S=\sum_iY_i$, `normalized_shape` propaga la
incertidumbre del total mediante:

$$J_{ij}=\frac{\delta_{ij}}S-\frac{Y_i}{S^2},\qquad
C=J\,\operatorname{diag}(V)\,J^{\mathsf T}.$$

Así se incluyen las correlaciones creadas al normalizar. El diagnóstico de forma
usa $\chi^2=\delta^{\mathsf T}(C_{\rm datos}+C_{\rm MC})^+\delta$, donde δ es la
diferencia de fracciones y + indica pseudoinversa. La covarianza es singular
porque las fracciones suman uno. El valor p asintótico es **solo estadístico**;
no es una significancia de descubrimiento ni una probabilidad de que el modelo
sea verdadero. No se combinan como independientes las dos regiones anidadas.

## 5. Controles y lo que sigue faltando

Se verifican tamaño y checksum de entradas, ramas, longitudes, finitud,
normalización y consistencia cinemática teniendo en cuenta el redondeo publicado.
Se detectan duplicados con identificadores de eventos durante la corrida local;
esos identificadores auxiliares no se incluyen en el paquete de publicación.
Los resultados guardan configuración y huellas de código e inputs.

Esto no sustituye una validación experimental: quedan incertidumbres de
eficiencias, etiquetado b, escalas/resoluciones, luminosidad, modelado y fondos no
prompt/falsos. No hay unfolding ni extracción de una correlación de espines.
Las selecciones ya aplicadas por ATLAS limitan la población estudiada.

Para estudiar el procedimiento, leer `select_events`, `event_weights`,
`observables`, `histogram`, `normalized_shape` y `process_file`, en ese orden.
Contrastar con la [documentación del formato](https://opendata.atlas.cern/docs/data/for_education/13TeV25_details)
y los [metadatos](https://opendata.atlas.cern/docs/data/for_education/13TeV25_metadata).
