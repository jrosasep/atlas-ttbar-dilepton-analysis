"""Entrada única del proyecto. Ejecutar sin argumentos muestra los pasos.

Este archivo organiza la ejecución; no cambia cortes, pesos ni resultados físicos.
Los cinco módulos originales se conservan byte por byte para mantener trazabilidad.
"""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CODE = ROOT / 'analysis'


def execute(*args):
    """Usa el mismo Python, detiene el proceso ante errores y evita cachés locales."""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([sys.executable, '-B', *args], cwd=CODE, env=env, check=True)


def export_results():
    """Tras dibujar, mueve las salidas para lectura a results/, sin duplicarlas."""
    destination = ROOT / 'results'
    destination.mkdir(exist_ok=True)
    for path in (CODE / 'figures/angular_v2').iterdir():
        if path.suffix in ('.png', '.pdf'):
            shutil.move(str(path), str(destination / path.name))
    for name in ('cutflow.csv', 'histograms.csv'):
        shutil.move(str(CODE / 'results/angular_v2' / name), str(destination / name))
    (CODE / 'figures/angular_v2').rmdir()
    (CODE / 'figures').rmdir()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step', nargs='?', choices=('check', 'plot', 'download', 'analyze'))
    parser.add_argument('--with-data', action='store_true', help='Incluye pruebas que requieren ROOT locales')
    args = parser.parse_args()
    if args.step is None:
        parser.print_help()
        print('\nEmpieza: python run.py check\nDespués: python run.py plot\n'
              'Solo para reprocesar eventos: download (14,97 GB) y luego analyze.')
    elif args.step == 'check':
        selection = [] if args.with_data else ['-m', 'not integration']
        execute('-m', 'pytest', '-p', 'no:cacheprovider', 'test_analysis.py', '-q', *selection)
    elif args.step == 'plot':
        execute('plots_and_summary.py')
        export_results()
        print('Figuras y tablas actualizadas en results/. Sin descargar datos.')
    elif args.step == 'download':
        print('Descarga explícita: 63 archivos, aproximadamente 14,97 GB.', flush=True)
        execute('fetch_inputs.py', 'download')
    elif args.step == 'analyze':
        if not (CODE / 'data').is_dir():
            parser.error('Faltan los ROOT en analysis/data/. Consulta README.md antes de descargar.')
        execute('analysis.py')
        execute('plots_and_summary.py')
        export_results()


if __name__ == '__main__':
    main()
