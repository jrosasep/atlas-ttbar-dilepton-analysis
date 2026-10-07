"""Elegir qué paso del proyecto ejecutar desde la carpeta principal.

check: comprobar el código mediante ejemplos pequeños, sin descargar datos.
plot: dibujar los histogramas ya guardados, verificando su procedencia.
download: obtener los ROOT públicos; es una descarga grande y explícita.
analyze: leer los ROOT locales, seleccionar eventos y crear nuevos resultados.

Este archivo organiza los pasos; no define cortes ni ecuaciones físicas.
"""
import argparse          # Leer la acción escrita después de run.py.
import os                # Conservar la configuración del entorno de Python.
from pathlib import Path # Construir rutas sin depender de dónde abrimos la terminal.
import subprocess        # Ejecutar otro programa y esperar su terminación.
import sys               # Usar el mismo intérprete de Python en todos los pasos.

# __file__ es el nombre de este programa. resolve obtiene la ruta completa;
# parent obtiene su carpeta. El símbolo / une una carpeta con un nombre.
ruta_del_programa = Path(__file__)
ruta_completa_del_programa = ruta_del_programa.resolve()
ROOT = ruta_completa_del_programa.parent
CODE = ROOT / 'analysis'


def execute(*args):
    """Ejecutar un paso; un error detiene la secuencia, no se oculta.

    *args reúne los argumentos recibidos. No significa multiplicación aquí.
    cwd fija la carpeta de ejecución para encontrar los otros programas.
    check=True exige que el programa termine sin errores.
    """
    entorno = dict(os.environ)
    entorno['PYTHONDONTWRITEBYTECODE'] = '1'
    comando = [sys.executable, '-B']
    comando.extend(args)
    subprocess.run(comando, cwd=CODE, env=entorno, check=True)


def main():
    """Relacionar la acción solicitada con el programa correspondiente."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'step', nargs='?', choices=('check', 'plot', 'download', 'analyze')
    )
    parser.add_argument(
        '--with-data', action='store_true',
        help='Añadir las pruebas de integración que requieren ROOT locales'
    )
    args = parser.parse_args()

    if args.step is None:
        parser.print_help()
        print('\nPrimero: python run.py check')
        print('Para ver los resultados guardados: python run.py plot')
        print('Para una nueva corrida: ROOT en analysis/data/ y python run.py analyze')
        print('download sólo es necesario si no tiene los ROOT; no se inicia solo.')

    elif args.step == 'check':
        # Las pruebas normales usan ejemplos construidos para comprobar cálculos.
        # Las marcadas integration necesitan archivos reales en analysis/data/.
        seleccion_de_pruebas = ['-m', 'not integration']
        if args.with_data:
            seleccion_de_pruebas = []
        execute(
            '-m', 'pytest', '-p', 'no:cacheprovider', 'test_analysis.py', '-q',
            *seleccion_de_pruebas
        )

    elif args.step == 'plot':
        # Los JSON guardados conservan la huella del código que los produjo.
        # Este modo permite redibujarlos tras comprobar las fuentes registradas;
        # no vuelve a procesar ROOT ni cambia su selección original.
        execute('plots_and_summary.py', '--stored-run')
        print('Figuras, tablas y procedencia del redibujo en results/.')

    elif args.step == 'download':
        print('Descarga explícita: 63 archivos, aproximadamente 15 GB.', flush=True)
        execute('fetch_inputs.py', 'download')

    elif args.step == 'analyze':
        carpeta_datos = CODE / 'data'
        if not carpeta_datos.is_dir():
            parser.error('Faltan los ROOT en analysis/data/. Consulte README.md.')
        execute('analysis.py')
        # Para resultados nuevos se exige coincidencia con el código actual.
        # Si el análisis falla, execute impide dibujar como si hubiese terminado.
        execute('plots_and_summary.py')


if __name__ == '__main__':
    main()
