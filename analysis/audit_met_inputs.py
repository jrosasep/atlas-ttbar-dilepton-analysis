"""Comprobar el momento transversal faltante (MET) y su redondeo.

El ROOT guarda MET y sus componentes x e y como float32. Su redondeo puede
producir una diferencia pequeña entre el módulo y hypot(x, y).
analysis.py importa interval_gap para comprobar esta coherencia.
Ejecutar este archivo directamente audita los ROOT, sin cambiar sus eventos.
"""
import json                        # Guardar los diagnósticos.
from pathlib import Path           # Construir rutas.
import numpy as np                 # Operaciones numéricas.
import uproot                      # Lectura de ROOT.

ruta_del_codigo = Path(__file__)
ruta_completa_del_codigo = ruta_del_codigo.resolve()
BASE = ruta_completa_del_codigo.parent


def rounding_interval(values):
    """Delimitar el intervalo de un redondeo a float32.

    nextafter encuentra el número representable vecino. Los puntos medios
    hacia los dos vecinos se calculan con precisión float64.
    """
    valores = np.asarray(values, dtype=np.float32)
    vecino_inferior = np.nextafter(valores, np.float32(-np.inf))
    vecino_superior = np.nextafter(valores, np.float32(np.inf))
    valores_double = valores.astype(float)
    inferior_double = vecino_inferior.astype(float)
    superior_double = vecino_superior.astype(float)
    inferior = (valores_double + inferior_double) / 2
    superior = (valores_double + superior_double) / 2
    return inferior, superior


def interval_gap(met, mx, my):
    """Separación entre los intervalos de MET y del módulo de (mx, my).

    Devuelve la separación y una cota aritmética. Cero significa que los
    intervalos se tocan. analysis.py permite además 0.01 GeV, como antes.
    """
    intervalos = []
    for componente in (mx, my):
        inferior, superior = rounding_interval(componente)
        contiene_cero = (inferior <= 0) & (superior >= 0)
        modulo_inferior = np.abs(inferior)
        modulo_superior = np.abs(superior)
        menor_extremo = np.minimum(modulo_inferior, modulo_superior)
        menor_modulo = np.where(contiene_cero, 0., menor_extremo)
        mayor_modulo = np.maximum(modulo_inferior, modulo_superior)
        intervalos.append((menor_modulo, mayor_modulo))
    min_x, max_x = intervalos[0]
    min_y, max_y = intervalos[1]
    radio_inferior = np.hypot(min_x, min_y)
    radio_superior = np.hypot(max_x, max_y)
    met_inferior, met_superior = rounding_interval(met)
    separacion_superior = met_inferior - radio_superior
    separacion_inferior = radio_inferior - met_superior
    mayor_separacion = np.maximum(separacion_superior, separacion_inferior)
    separacion = np.maximum(0., mayor_separacion)
    precision_double = np.finfo(float).eps
    escala = np.maximum(radio_superior, 1.)
    cota_aritmetica = 8 * precision_double * escala
    return separacion, cota_aritmetica


def main():
    """Leer las tres ramas MET de cada archivo y guardar sus controles."""
    ruta_manifiesto = BASE / 'manifest.json'
    texto = ruta_manifiesto.read_text()
    manifest = json.loads(texto)
    auditoria = []
    for sample in manifest['files']:
        resultado = {
            'file': sample['key'], 'dsid': sample['dsid'], 'entries': 0,
            'max_residual_float64_gev': 0., 'max_interval_gap_gev': 0.,
            'nonoverlapping_single_rounding_intervals': 0,
            'violations_with_0p01_gev_allowance': 0,
            'met_over_13tev': 0, 'bad_examples': [],
        }
        ruta_root = BASE / 'data' / sample['key']
        # with cierra el archivo al terminar este bloque.
        with uproot.open(ruta_root) as archivo:
            arbol = archivo['analysis']
            bloques = arbol.iterate(
                ['met', 'met_mpx', 'met_mpy'], step_size=500000, library='np'
            )
            for bloque in bloques:
                modulo = bloque['met']
                x = bloque['met_mpx']
                y = bloque['met_mpy']
                x_double = x.astype(float)
                y_double = y.astype(float)
                modulo_double = modulo.astype(float)
                calculado = np.hypot(x_double, y_double)
                residuo = np.abs(calculado - modulo_double)
                separacion, cota = interval_gap(modulo, x, y)
                separados = separacion > cota
                posiciones = np.flatnonzero(separados)
                for posicion in posiciones[:3]:
                    if len(resultado['bad_examples']) < 3:
                        ejemplo = {
                            'entry': int(resultado['entries'] + posicion),
                            'met': float(modulo[posicion]),
                            'mx': float(x[posicion]), 'my': float(y[posicion]),
                            'residual_gev': float(residuo[posicion]),
                            'interval_gap_gev': float(separacion[posicion]),
                        }
                        resultado['bad_examples'].append(ejemplo)
                resultado['entries'] += len(modulo)
                resultado['max_residual_float64_gev'] = max(
                    resultado['max_residual_float64_gev'], float(np.max(residuo))
                )
                resultado['max_interval_gap_gev'] = max(
                    resultado['max_interval_gap_gev'], float(np.max(separacion))
                )
                resultado['nonoverlapping_single_rounding_intervals'] += int(separados.sum())
                viola_contrato = separacion > .01 + cota
                resultado['violations_with_0p01_gev_allowance'] += int(viola_contrato.sum())
                resultado['met_over_13tev'] += int((modulo > 13000).sum())
        auditoria.append(resultado)
        viola = resultado['violations_with_0p01_gev_allowance'] > 0
        visible = resultado['max_residual_float64_gev'] > .005
        if viola or visible:
            print(json.dumps(resultado), flush=True)
    carpeta = BASE / 'provenance'
    carpeta.mkdir(parents=True, exist_ok=True)
    salida = carpeta / 'met_intervals_audit.json'
    salida.write_text(json.dumps(auditoria, indent=2))
    total_violaciones = 0
    for resultado in auditoria:
        total_violaciones += resultado['violations_with_0p01_gev_allowance']
    print('VIOLACIONES DEL CONTRATO', total_violaciones, flush=True)


if __name__ == '__main__':
    main()
