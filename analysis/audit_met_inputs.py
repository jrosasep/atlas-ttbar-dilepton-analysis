"""Audita TODOS los MET con intervalos de redondeo de sus columnas float32.

Si x es un número redondeado, el valor anterior al redondeo se encuentra entre
los puntos medios hacia sus dos vecinos representables. Se propagan esos
intervalos al módulo sqrt(x²+y²), sin ajustar tolerancias a un histograma.
El contrato admite una separación residual de 0.01 GeV después de considerar
esa cuantización. El solapamiento exacto no se exige: la producción del ntuple
puede implicar más operaciones que un único redondeo final.
"""
import json
from pathlib import Path
import numpy as np
import uproot

BASE = Path(__file__).resolve().parent


def rounding_interval(values):
    """Intervalo de números reales que pueden redondearse a cada float32."""
    x = np.asarray(values, dtype=np.float32)
    lower = (x.astype(float) + np.nextafter(x, np.float32(-np.inf)).astype(float))/2
    upper = (x.astype(float) + np.nextafter(x, np.float32(np.inf)).astype(float))/2
    return lower, upper


def interval_gap(met, mx, my):
    """Distancia entre el intervalo de MET y el de su módulo; cero si se tocan."""
    ranges = []
    for values in (mx, my):
        lower, upper = rounding_interval(values)
        smallest = np.where((lower <= 0) & (upper >= 0), 0., np.minimum(abs(lower), abs(upper)))
        largest = np.maximum(abs(lower), abs(upper))
        ranges.append((smallest, largest))
    radius_lower = np.hypot(ranges[0][0], ranges[1][0])
    radius_upper = np.hypot(ranges[0][1], ranges[1][1])
    met_lower, met_upper = rounding_interval(met)
    gap = np.maximum(0., np.maximum(met_lower-radius_upper, radius_lower-met_upper))
    # Cota conservadora del cálculo auxiliar, que sí se hace en float64.
    arithmetic = 8*np.finfo(float).eps*np.maximum(radius_upper, 1.)
    return gap, arithmetic


def main():
    manifest = json.loads((BASE / 'manifest.json').read_text())
    audit = []
    for sample in manifest['files']:
        result = {'file': sample['key'], 'dsid': sample['dsid'], 'entries': 0,
                  'max_residual_float64_gev': 0., 'max_interval_gap_gev': 0.,
                  'nonoverlapping_single_rounding_intervals': 0,
                  'violations_with_0p01_gev_allowance': 0,
                  'met_over_13tev': 0, 'bad_examples': []}
        with uproot.open(BASE / 'data' / sample['key']) as root:
            for a in root['analysis'].iterate(['met', 'met_mpx', 'met_mpy'], step_size=500000, library='np'):
                m, x, y = a['met'], a['met_mpx'], a['met_mpy']
                residual = abs(np.hypot(x.astype(float), y.astype(float))-m.astype(float))
                gap, arithmetic = interval_gap(m, x, y)
                bad = gap > arithmetic
                for i in np.flatnonzero(bad)[:3]:
                    if len(result['bad_examples']) < 3:
                        result['bad_examples'].append({'entry': int(result['entries']+i),
                            'met': float(m[i]), 'mx': float(x[i]), 'my': float(y[i]),
                            'residual_gev': float(residual[i]), 'interval_gap_gev': float(gap[i])})
                result['entries'] += len(m)
                result['max_residual_float64_gev'] = max(result['max_residual_float64_gev'], float(np.max(residual)))
                result['max_interval_gap_gev'] = max(result['max_interval_gap_gev'], float(np.max(gap)))
                result['nonoverlapping_single_rounding_intervals'] += int(bad.sum())
                result['violations_with_0p01_gev_allowance'] += int((gap > .01+arithmetic).sum())
                result['met_over_13tev'] += int((m > 13000).sum())
        audit.append(result)
        if result['violations_with_0p01_gev_allowance'] or result['max_residual_float64_gev'] > .005:
            print(json.dumps(result), flush=True)
    (BASE / 'provenance/met_intervals_audit.json').write_text(json.dumps(audit, indent=2))
    print('VIOLACIONES DEL CONTRATO', sum(r['violations_with_0p01_gev_allowance'] for r in audit), flush=True)


if __name__ == '__main__':
    main()
