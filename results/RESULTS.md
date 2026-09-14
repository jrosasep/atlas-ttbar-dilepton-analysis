# Resultado de la preparación local

La corrida angular procesó los 63 archivos completos: 6 242 521 eventos reales
y 27 100 357 simulados. Son eventos del skim educativo, no todas las colisiones
registradas por ATLAS. Se verificó que los eventos reales no están duplicados.

| Selección | Colisiones seleccionadas | Predicción MC | Error estadístico MC |
|---|---:|---:|---:|
| Al menos 1 etiqueta b (auxiliar) | 75 218 | 76 247,3 | 189,8 |
| Al menos 2 etiquetas b (principal) | 36 945 | 36 691,3 | 131,9 |

La incertidumbre estadística aproximada del conteo real es la raíz cuadrada del
conteo; no está incluida en la columna de error MC. Estos errores no incluyen
las incertidumbres sistemáticas. La región principal está contenida en la auxiliar.

La contribución ttbar representa el 96,57 % de la predicción de la región
principal. Esta fracción pertenece al modelo disponible: no demuestra que ese
porcentaje exacto de las colisiones reales sea ttbar.

## Verificaciones efectuadas

- 18 pruebas completas pasan en el entorno con los ROOT.
- 13 pruebas pasan sin ROOT; las 5 de integración se excluyen explícitamente.
- Los 280 histogramas anteriores, conteos por corte y rendimientos permanecen
  idénticos al añadir |Δη|. Los resultados anteriores se conservaron aparte.
- Los nuevos histogramas |Δη| recuperan los totales de todas las selecciones,
  sin eventos fuera del intervalo angular definido.
- Las figuras se regeneraron también desde la carpeta preparada, sin ROOT.
- Se revisaron visualmente las figuras angulares principales de conteos y formas.
- El ZIP se comprobó íntegro. Incluye huellas SHA256 y omite ROOT, identificadores
  auxiliares de eventos, entorno Python y documentos personales.

Esto demuestra que el cálculo descriptivo está implementado y es reproducible
en el entorno probado. No constituye una validación completa del modelo físico,
del tratamiento de fondos o de todas las incertidumbres.

## Trabajo científico pendiente

Profundizar la revisión bibliográfica, validar la selección frente a las
referencias, evaluar los fondos incompletos y las incertidumbres sistemáticas.
Las pruebas automatizadas no sustituyen estas tareas. La procedencia y licencia
de los datos se documentan en ../licenses/DATA_SOURCES.md.
