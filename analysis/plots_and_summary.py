"""Combinar resultados numéricos y dibujar datos frente a MC.

Entrada: run.json y los JSON de cada muestra, producidos por analysis.py.
Salida: gráficos, tablas y resumen en la carpeta results/ del repositorio.
No vuelve a seleccionar eventos ROOT ni ajusta la predicción a los datos.

Modo normal: exige el código de análisis que produjo la corrida.
--stored-run: redibuja una corrida anterior tras verificar su fuente registrada.
Este último modo conserva las huellas originales y registra el código utilizado
ahora para el postprocesamiento. No convierte una corrida antigua en una nueva.

El contraste incluye estadística; aún faltan sistemáticas y fondos instrumentales.
"""
from __future__ import annotations
import argparse                    # Opciones de ejecución.
import csv                         # Tablas legibles de números.
import hashlib                     # Huellas de archivos y programas.
import json                        # Leer los resultados guardados.
from pathlib import Path           # Rutas y nombres de archivos.
import matplotlib
matplotlib.use('Agg')               # Guardar gráficos sin abrir ventanas.
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2        # Cuantiles de conteo y valores p aproximados.
from analysis import BASE, RESULTS, BINS, add_hist, normalized_shape, pipeline_hashes

LABELS = {
    'data': 'Colisiones reales', 'ttbar': r'$t\bar{t}$',
    'single_top': 'Top individual', 'diboson': 'Dos bosones',
    'Zjets': 'Z + jets', 'Wjets': 'W + jets',
    'ttV': r'$t\bar{t}$ + bosones', 'HWW': r'Higgs $\to WW$',
}
COLORS = {
    'ttbar': '#4b83b5', 'single_top': '#e6a657', 'diboson': '#70a78a',
    'Zjets': '#bd83b8', 'Wjets': '#d17474', 'ttV': '#bdc766', 'HWW': '#929292',
}
ORDER = ['HWW', 'ttV', 'Wjets', 'Zjets', 'diboson', 'single_top', 'ttbar']
RUN_LABEL = ''
LUMINOSITY_FB = 36.0
XLABELS = {
    'electron_pt': r'$p_T(e)$ [GeV]', 'muon_pt': r'$p_T(\mu)$ [GeV]',
    'n_jets': 'Número de jets aceptados',
    'delta_phi': r'$|\Delta\phi(e,\mu)|$ [rad]',
    'delta_eta': r'$|\Delta\eta(e,\mu)|$ (sin unidades)',
    'm_emu': r'$m(e,\mu)$ [GeV]', 'met': 'Momento transversal faltante [GeV]',
    'n_bjets': 'Número de jets con etiqueta b',
}


def file_hash(path):
    """Leer los bytes y devolver su huella SHA256, sin modificar el archivo."""
    contenido = path.read_bytes()
    calculo = hashlib.sha256(contenido)
    huella = calculo.hexdigest()
    return huella


def recorded_source_directory(run):
    """Ubicar las fuentes de procedencia de una corrida registrada."""
    huella = run['code_sha256']
    caracteres_hexadecimales = '0123456789abcdef'
    if len(huella) != 64 or any(c not in caracteres_hexadecimales for c in huella):
        raise ValueError('Huella de código registrada inválida')
    carpeta_registros = BASE / 'provenance' / 'source_records'
    return carpeta_registros / huella


def verify_recorded_sources(run):
    """Comprobar las fuentes archivadas; nunca ejecutar ese código archivado.

    Los archivos .py.txt son evidencia de procedencia de los histogramas.
    Cada huella debe coincidir exactamente con la que registró analysis.py.
    """
    huellas = run.get('pipeline_hashes', {})
    nombres_esperados = {'analysis.py', 'audit_met_inputs.py', 'fetch_inputs.py'}
    if set(huellas) != nombres_esperados:
        raise ValueError('La corrida no identifica todas sus fuentes de cálculo')
    if run['code_sha256'] != huellas['analysis.py']:
        raise ValueError('Las huellas de analysis.py de la corrida son contradictorias')
    carpeta = recorded_source_directory(run)
    for nombre, huella_esperada in huellas.items():
        ruta_fuente = carpeta / (nombre + '.txt')
        if not ruta_fuente.exists():
            raise ValueError('Falta la fuente registrada: ' + str(ruta_fuente))
        if file_hash(ruta_fuente) != huella_esperada:
            raise ValueError('La fuente registrada fue modificada: ' + nombre)


def verified_run_manifest(run):
    """Recuperar la lista de entradas que coincide con la huella de la corrida."""
    candidatos = [
        BASE / 'provenance' / 'verified_inputs.json',
        BASE / 'manifest.json',
        recorded_source_directory(run) / 'input_manifest.json.txt',
    ]
    for ruta in candidatos:
        if ruta.exists() and file_hash(ruta) == run['manifest_sha256']:
            texto = ruta.read_text(encoding='utf-8')
            return json.loads(texto)
    raise ValueError('No se encuentra el manifiesto que corresponde a la corrida')


def validate_run(run, stored_run=False, partial=False):
    """Distinguir una corrida actual de un redibujo documentado de la anterior."""
    if not partial:
        if not run['complete']:
            raise ValueError('No se pueden publicar figuras completas con una corrida parcial')
        if not run.get('global_checks', {}).get('all_checks_passed'):
            raise ValueError('La corrida no completó los controles globales')
    huellas_actuales = pipeline_hashes()
    mismo_codigo = run['code_sha256'] == file_hash(BASE / 'analysis.py')
    mismos_controles = run.get('pipeline_hashes') == huellas_actuales
    if not mismo_codigo or not mismos_controles:
        if not stored_run:
            raise ValueError(
                'El código actual difiere de la corrida. Reprocese los ROOT o '
                'use --stored-run para un redibujo con procedencia verificada.'
            )
        verify_recorded_sources(run)
    manifiesto = verified_run_manifest(run)
    return manifiesto


def poisson_errors(counts):
    """Barras inferior y superior de Garwood al 68.27 % para conteos Poisson.

    Son asimétricas y conservan incertidumbre superior cuando se observan cero
    eventos. ppf obtiene un cuantil de la distribución chi-cuadrado.
    """
    n = np.asarray(counts, float)
    alpha = 1 - .682689492137
    limite_inferior = np.zeros_like(n)
    positivos = n > 0
    probabilidad_inferior = alpha / 2
    grados_inferiores = 2 * n[positivos]
    limite_inferior[positivos] = .5 * chi2.ppf(probabilidad_inferior, grados_inferiores)
    probabilidad_superior = 1 - alpha / 2
    grados_superiores = 2 * (n + 1)
    limite_superior = .5 * chi2.ppf(probabilidad_superior, grados_superiores)
    error_inferior = n - limite_inferior
    error_superior = limite_superior - n
    return np.array([error_inferior, error_superior])


def aggregate(run):
    """Sumar archivos por proceso físico sin mezclar códigos ni selecciones."""
    groups = {}
    for ruta_relativa in run['files']:
        # Los registros antiguos pueden guardar separadores de Windows.
        ruta_portable = Path(ruta_relativa.replace('\\', '/'))
        if ruta_portable.is_absolute() or '..' in ruta_portable.parts:
            raise ValueError('Ruta de resultado fuera de la carpeta del análisis')
        ruta_muestra = BASE / ruta_portable
        muestra = json.loads(ruta_muestra.read_text())
        if muestra['code_sha256'] != run['code_sha256']:
            raise ValueError('Se intentó mezclar versiones de código')
        if muestra['config'] != run['config']:
            raise ValueError('Se intentó mezclar selecciones')
        if muestra.get('pipeline_hashes') != run.get('pipeline_hashes'):
            raise ValueError('La muestra no corresponde a los controles de la corrida')

        nombre_grupo = muestra['group']
        if nombre_grupo not in groups:
            groups[nombre_grupo] = {
                'entries': 0, 'files': 0, 'cutflow': {},
                'regions': {}, 'histograms': {}, 'negative_weights': 0,
            }
        grupo = groups[nombre_grupo]
        grupo['entries'] += muestra['entries']
        grupo['files'] += 1
        grupo['negative_weights'] += muestra['negative_final_weights']
        for categoria in ('cutflow', 'regions'):
            for etiqueta, contador in muestra[categoria].items():
                if etiqueta not in grupo[categoria]:
                    grupo[categoria][etiqueta] = {'raw': 0, 'sumw': 0., 'sumw2': 0.}
                destino = grupo[categoria][etiqueta]
                for campo in destino:
                    destino[campo] += contador[campo]
        for region, histogramas in muestra['histograms'].items():
            if region not in grupo['histograms']:
                grupo['histograms'][region] = {}
            destino = grupo['histograms'][region]
            for observable, histograma in histogramas.items():
                if observable in destino:
                    add_hist(destino[observable], histograma)
                else:
                    destino[observable] = histograma
    return groups


def check_saved_results(run, manifest, groups):
    """Comprobar la conservación de totales y la coherencia de los registros."""
    if len(set(run['files'])) != len(run['files']):
        raise ValueError('Una muestra se contó dos veces en la corrida')
    catalogo = {}
    for entrada in manifest['files']:
        catalogo[entrada['key']] = entrada
    entradas_vistas = set()
    for ruta in run['files']:
        ruta = Path(ruta.replace('\\', '/'))
        muestra = json.loads((BASE / ruta).read_text())
        if muestra['file'] not in catalogo:
            raise ValueError('Muestra ausente del manifiesto registrado')
        if muestra['file'] in entradas_vistas:
            raise ValueError('Dos resultados contienen la misma entrada ROOT')
        entradas_vistas.add(muestra['file'])
        entrada = catalogo[muestra['file']]
        if muestra['group'] != entrada['group'] or muestra['metadata'] != entrada['metadata']:
            raise ValueError('Grupo o metadatos distintos del manifiesto')
        if entrada.get('sha256') and muestra['input_sha256'] != entrada['sha256']:
            raise ValueError('Huella de entrada distinta del manifiesto')
        if muestra['processed'] != muestra['entries']:
            raise ValueError('La muestra no procesó todas sus entradas')
    if run.get('complete') and entradas_vistas != set(catalogo):
        raise ValueError('La corrida completa no contiene todas las entradas del manifiesto')
    for grupo in groups.values():
        for region, histogramas in grupo['histograms'].items():
            contador = grupo['regions'][region]
            for observable, h in histogramas.items():
                cantidad = len(BINS[observable]) - 1
                for campo in ('raw', 'sumw', 'sumw2'):
                    valores = np.asarray(h[campo])
                    if len(valores) != cantidad or not np.all(np.isfinite(valores)):
                        raise ValueError('Dimensiones o valores inválidos en el histograma')
                if sum(h['raw']) != contador['raw']:
                    raise ValueError('El histograma perdió conteos')
                if not np.isclose(sum(h['sumw']), contador['sumw']):
                    raise ValueError('El histograma perdió pesos')
                if not np.isclose(sum(h['sumw2']), contador['sumw2']):
                    raise ValueError('El histograma perdió varianza')


def total_mc(groups, region, variable):
    """Sumar la predicción nominal y sus varianzas, conservando los signos."""
    numero_intervalos = len(BINS[variable]) - 1
    suma_pesos = np.zeros(numero_intervalos)
    suma_varianzas = np.zeros(numero_intervalos)
    for nombre_grupo in ORDER:
        if nombre_grupo in groups:
            histograma = groups[nombre_grupo]['histograms'][region][variable]
            suma_pesos += histograma['sumw']
            suma_varianzas += histograma['sumw2']
    return suma_pesos, suma_varianzas


def comparison_panel(ax, ratio, groups, region, variable):
    """Dibujar conteos y cociente: ax y ratio son dos paneles del gráfico."""
    bordes = BINS[variable]
    centros = (bordes[:-1] + bordes[1:]) / 2
    acumulado = np.zeros(len(centros))
    for nombre_grupo in ORDER:
        if nombre_grupo not in groups:
            continue
        histograma = groups[nombre_grupo]['histograms'][region][variable]
        contribucion = np.array(histograma['sumw'])
        techo = acumulado + contribucion
        ax.stairs(techo, bordes, baseline=acumulado, fill=True,
                  color=COLORS[nombre_grupo], label=LABELS[nombre_grupo], linewidth=.4)
        acumulado += contribucion

    mc, varianza = total_mc(groups, region, variable)
    error_mc = np.sqrt(varianza)
    ax.stairs(mc + error_mc, bordes, baseline=mc - error_mc, fill=True,
              facecolor='none', edgecolor='#444444', hatch='////',
              linewidth=0, label='Estadística MC')
    datos = np.array(groups['data']['histograms'][region][variable]['raw'])
    errores_datos = poisson_errors(datos)
    ax.errorbar(centros, datos, yerr=errores_datos, fmt='o', color='#161616',
                markersize=3, linewidth=.8, label='Datos reales')
    ax.set_ylabel('Eventos / intervalo')
    max_datos = float(np.max(datos + errores_datos[1]))
    max_mc = float(np.max(mc + error_mc))
    techo_eje = max(max_datos, max_mc, 1) * 1.28
    ax.set_ylim(bottom=0, top=techo_eje)
    ax.grid(axis='y', alpha=.14)

    # Se divide sólo donde la predicción es positiva.
    validos = mc > 0
    error_relativo = np.divide(error_mc, mc, out=np.zeros_like(error_mc), where=validos)
    ratio.stairs(1 + error_relativo, bordes, baseline=1 - error_relativo,
                 fill=True, facecolor='#cccccc', linewidth=0)
    cociente = datos[validos] / mc[validos]
    errores_cociente = errores_datos[:, validos] / mc[validos]
    ratio.errorbar(centros[validos], cociente, yerr=errores_cociente,
                   fmt='o', color='black', markersize=3, linewidth=.8)
    ratio.axhline(1, color='#4b5563', lw=.8)
    ratio.set_ylabel('Datos / MC')
    ratio.set_xlabel(XLABELS[variable])
    if validos.any():
        minimo_datos = float(np.min((datos[validos] - errores_datos[0, validos]) / mc[validos])) - .04
        minimo_mc = float(np.min(1 - error_relativo[validos])) - .04
        maximo_datos = float(np.max((datos[validos] + errores_datos[1, validos]) / mc[validos])) + .04
        maximo_mc = float(np.max(1 + error_relativo[validos])) + .04
        inferior = min(.8, minimo_datos, minimo_mc)
        superior = max(1.2, maximo_datos, maximo_mc)
        ratio.set_ylim(max(0, inferior), superior)
    else:
        ratio.set_ylim(0, 2)
    ratio.grid(axis='y', alpha=.2)
    if variable in ('n_jets', 'n_bjets'):
        posiciones = centros.astype(int)
        etiquetas = []
        for posicion in posiciones[:-1]:
            etiquetas.append(str(posicion))
        etiquetas.append(f'>={posiciones[-1]}')
        ratio.set_xticks(posiciones, etiquetas)
    elif variable not in ('delta_phi', 'delta_eta'):
        texto = f'Último intervalo incluye valores >= {bordes[-1]:g}'
        ratio.text(.98, -.62, texto, transform=ratio.transAxes, ha='right', fontsize=7)


def draw_comparisons(groups, region, out):
    """Dibujar cuatro observables de control, cada uno con su cociente."""
    fig = plt.figure(figsize=(12, 9.5), layout='constrained')
    grid = fig.add_gridspec(4, 2, height_ratios=[3, 1, 3, 1])
    axes = []
    observables = ('electron_pt', 'muon_pt', 'n_jets', 'delta_phi')
    for posicion, variable in enumerate(observables):
        fila = 2 * (posicion // 2)
        columna = posicion % 2
        ax = fig.add_subplot(grid[fila, columna])
        ratio = fig.add_subplot(grid[fila + 1, columna], sharex=ax)
        plt.setp(ax.get_xticklabels(), visible=False)
        comparison_panel(ax, ratio, groups, region, variable)
        axes.append(ax)
    axes[0].legend(fontsize=9, ncols=3, loc='upper right')
    if region == 'SR1b':
        etiqueta = 'al menos una'
    else:
        etiqueta = 'al menos dos'
    titulo = RUN_LABEL + f'ATLAS Open Data: e-mu, cargas opuestas, {etiqueta} etiqueta(s) b\n'
    titulo += f'2015-2016 | 13 TeV | L = {LUMINOSITY_FB:g} fb$^{{-1}}$ (redondeada) | predicción nominal sin ajuste'
    fig.suptitle(titulo, fontsize=12)
    fig.savefig(out / f'datos_mc_{region}.png', dpi=175)
    plt.close(fig)


def draw_angles(groups, region, out, normalized=False):
    """Dibujar Delta-phi y |Delta-eta| en conteos o fracciones por intervalo.

    Cada histograma se normaliza por su propio total. Los fondos siguen
    presentes. Las covarianzas se calculan con normalized_shape de analysis.py.
    """
    fig = plt.figure(figsize=(11, 5.8), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1])
    for columna, observable in enumerate(('delta_phi', 'delta_eta')):
        ax = fig.add_subplot(grid[0, columna])
        inferior = fig.add_subplot(grid[1, columna], sharex=ax)
        if not normalized:
            comparison_panel(ax, inferior, groups, region, observable)
        else:
            bordes = BINS[observable]
            centros = (bordes[:-1] + bordes[1:]) / 2
            datos = np.asarray(groups['data']['histograms'][region][observable]['raw'], float)
            mc, varianza = total_mc(groups, region, observable)
            fraccion_datos, cov_datos = normalized_shape(datos, datos)
            fraccion_mc, cov_mc = normalized_shape(mc, varianza)
            var_datos = np.maximum(np.diag(cov_datos), 0)
            var_mc = np.maximum(np.diag(cov_mc), 0)
            error_datos = np.sqrt(var_datos)
            error_mc = np.sqrt(var_mc)
            ax.stairs(fraccion_mc, bordes, color='#3476a8', label='MC nominal + fondos')
            ax.stairs(fraccion_mc + error_mc, bordes, baseline=fraccion_mc - error_mc,
                      fill=True, color='#3476a8', alpha=.2, label='Incertidumbre estadística MC')
            ax.errorbar(centros, fraccion_datos, yerr=error_datos,
                        fmt='o', color='black', markersize=4, label='Colisiones reales')
            ax.set_ylabel('Fracción de eventos / intervalo')
            validos = fraccion_mc > 0
            cociente = fraccion_datos[validos] / fraccion_mc[validos]
            error_cociente = error_datos[validos] / fraccion_mc[validos]
            inferior.errorbar(centros[validos], cociente, yerr=error_cociente,
                              fmt='o', color='black', markersize=3)
            banda = np.divide(error_mc, fraccion_mc, out=np.zeros_like(fraccion_mc), where=validos)
            inferior.stairs(1 + banda, bordes, baseline=1 - banda,
                            fill=True, color='#3476a8', alpha=.2)
            inferior.axhline(1, color='gray', lw=.8)
            inferior.set_ylabel('Datos / MC')
            inferior.set_xlabel(XLABELS[observable])
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.grid(axis='y', alpha=.15)
        inferior.grid(axis='y', alpha=.15)
        if columna == 0:
            handles, labels = ax.get_legend_handles_labels()
    etiqueta = 'al menos 1' if region == 'SR1b' else 'al menos 2'
    modo = 'formas normalizadas' if normalized else 'conteos absolutos'
    titulo = RUN_LABEL + f'ATLAS Open Data | e-mu, {etiqueta} b-tags | {modo}\n'
    titulo += f'2015-2016 | 13 TeV | {LUMINOSITY_FB:g} fb$^{{-1}}$ (redondeada) | nivel reconstruido'
    fig.suptitle(titulo, fontsize=11)
    fig.legend(handles, labels, loc='outside lower center', ncols=3, fontsize=8)
    nombre_modo = 'formas' if normalized else 'conteos'
    nombre = f'angulos_{nombre_modo}_{region}'
    fig.savefig(out / (nombre + '.png'), dpi=175)
    fig.savefig(out / (nombre + '.pdf'))
    plt.close(fig)


def diagnostics(groups):
    """Calcular rendimientos, fracción ttbar y contraste estadístico de formas.

    La fracción ttbar pertenece al modelo MC. El valor p aproximado no mide
    la probabilidad de que el modelo sea verdadero ni una nueva física.
    """
    results = {}
    for region in ('PRETAG', 'ZEROb', 'SR1b', 'SR2b', 'SS1b'):
        numero_datos = groups['data']['regions'][region]['raw']
        contribuciones = []
        varianzas = []
        for nombre, grupo in groups.items():
            if nombre != 'data':
                contribuciones.append(grupo['regions'][region]['sumw'])
                varianzas.append(grupo['regions'][region]['sumw2'])
        # sum conserva la misma operación de la implementación anterior.
        # En Python moderno su suma de floats reduce el error de redondeo.
        prediccion = sum(contribuciones)
        varianza = sum(varianzas)
        señal = groups['ttbar']['regions'][region]['sumw']
        if prediccion > 0:
            fraccion_ttbar = señal / prediccion
            cociente_global = numero_datos / prediccion
        else:
            fraccion_ttbar = None
            cociente_global = None
        results[region] = {
            'data': numero_datos, 'mc': prediccion,
            'mc_stat': float(np.sqrt(varianza)),
            'ttbar_fraction_in_model': fraccion_ttbar,
            'data_over_mc': cociente_global, 'shapes': {},
        }
        for observable in BINS:
            datos = np.array(groups['data']['histograms'][region][observable]['raw'], float)
            mc, var_mc = total_mc(groups, region, observable)
            if mc.sum() <= 0 or datos.sum() <= 0:
                continue
            p_datos, c_datos = normalized_shape(datos, datos)
            p_mc, c_mc = normalized_shape(mc, var_mc)
            covarianza = c_datos + c_mc
            rango = int(np.linalg.matrix_rank(covarianza))
            diferencia = p_datos - p_mc
            pseudoinversa = np.linalg.pinv(covarianza)
            producto = diferencia @ pseudoinversa
            valor = float(producto @ diferencia)
            valor_p = float(chi2.sf(valor, rango)) if rango else None
            results[region]['shapes'][observable] = {
                'chi2_stat_shape': valor, 'rank': rango,
                'p_approx_stat_only': valor_p,
                'covariance_data': c_datos.tolist(),
                'covariance_mc': c_mc.tolist(),
                'warning': 'Diagnóstico asintótico solo estadístico; sin sistemáticas, no es significancia de descubrimiento.',
            }
    return results


def draw_control(groups, out):
    """Dibujar controles auxiliares de la selección y del modelado."""
    fig = plt.figure(figsize=(11, 5), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1])
    controles = (('PRETAG', 'n_bjets'), ('SS1b', 'delta_phi'))
    for columna, (region, observable) in enumerate(controles):
        ax = fig.add_subplot(grid[0, columna])
        ratio = fig.add_subplot(grid[1, columna])
        comparison_panel(ax, ratio, groups, region, observable)
        titulo = 'Antes de exigir etiqueta b' if region == 'PRETAG' else 'Control: cargas iguales'
        ax.set_title(titulo, fontsize=11)
        if columna == 0:
            handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8.5, ncols=5, loc='outside lower center')
    fig.suptitle(RUN_LABEL + 'Controles de selección y de modelado | incertidumbres estadísticas', fontsize=12)
    fig.savefig(out / 'controles.png', dpi=175)
    plt.close(fig)


def draw_extra(groups, out):
    """Dibujar masa e-mu y MET como controles, sin retocar la selección."""
    fig = plt.figure(figsize=(11, 5), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1])
    for columna, observable in enumerate(('m_emu', 'met')):
        ax = fig.add_subplot(grid[0, columna])
        ratio = fig.add_subplot(grid[1, columna])
        comparison_panel(ax, ratio, groups, 'SR1b', observable)
        if columna == 0:
            handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8.5, ncols=5, loc='outside lower center')
    fig.suptitle(RUN_LABEL + 'SR1b: masa electrón-muón y momento faltante | errores estadísticos', fontsize=12)
    fig.savefig(out / 'variables_adicionales.png', dpi=175)
    plt.close(fig)


def write_tables(groups, out, partial=False):
    """Exportar cada conteo y barra del gráfico a tablas CSV."""
    with (out / 'cutflow.csv').open('w', newline='', encoding='utf-8') as archivo:
        writer = csv.writer(archivo)
        writer.writerow(['group', 'stage', 'raw', 'sumw', 'sumw2', 'stat_error'])
        for nombre, grupo in groups.items():
            for etapa, contador in grupo['cutflow'].items():
                error = np.sqrt(contador['sumw2'])
                writer.writerow([nombre, etapa, contador['raw'], contador['sumw'], contador['sumw2'], error])
    nombre_tabla = 'histograms_partial.csv' if partial else 'histograms.csv'
    with (out / nombre_tabla).open('w', newline='', encoding='utf-8') as archivo:
        writer = csv.writer(archivo)
        writer.writerow(['group', 'region', 'observable', 'bin_low', 'bin_high',
                         'last_bin_includes_overflow', 'raw', 'sumw', 'sumw2'])
        for nombre, grupo in groups.items():
            for region, histogramas in grupo['histograms'].items():
                for observable, histograma in histogramas.items():
                    bordes = BINS[observable]
                    for posicion in range(len(bordes) - 1):
                        es_ultimo = posicion == len(bordes) - 2
                        incluye_desborde = es_ultimo and observable not in ('delta_phi', 'delta_eta')
                        writer.writerow([nombre, region, observable, bordes[posicion], bordes[posicion + 1],
                                         incluye_desborde, histograma['raw'][posicion],
                                         histograma['sumw'][posicion], histograma['sumw2'][posicion]])


def main():
    """Leer una corrida, verificarla y guardar sólo resultados derivados."""
    global RUN_LABEL, BINS, LUMINOSITY_FB
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partial', action='store_true')
    parser.add_argument('--stored-run', action='store_true',
                        help='Redibuja una corrida anterior con fuentes de procedencia verificadas')
    args = parser.parse_args()
    RUN_LABEL = 'CORRIDA PARCIAL - faltan muestras\n' if args.partial else ''
    nombre_corrida = 'run_partial.json' if args.partial else 'run.json'
    ruta_corrida = RESULTS / nombre_corrida
    run = json.loads(ruta_corrida.read_text())
    manifiesto = validate_run(run, stored_run=args.stored_run, partial=args.partial)

    # Los bordes y la luminosidad pertenecen a la corrida leída, no a una
    # configuración que el estudiante pueda haber cambiado después.
    bordes_registrados = {}
    for observable, bordes in run['bins'].items():
        bordes_registrados[observable] = np.asarray(bordes, float)
    BINS = bordes_registrados
    LUMINOSITY_FB = run['config']['luminosity_pb'] / 1000
    groups = aggregate(run)
    check_saved_results(run, manifiesto, groups)
    diagnosticos = diagnostics(groups)

    out = BASE.parent / 'results'
    if args.partial:
        out = out / 'partial'
    out.mkdir(parents=True, exist_ok=True)
    for region in ('SR1b', 'SR2b'):
        draw_comparisons(groups, region, out)
        draw_angles(groups, region, out)
        draw_angles(groups, region, out, normalized=True)
    draw_control(groups, out)
    draw_extra(groups, out)
    write_tables(groups, out, partial=args.partial)

    # Mantener separada la fuente del procesamiento y la del dibujo actual.
    ruta_programa = Path(__file__)
    procedencia = {
        'stored_run_mode': args.stored_run,
        'run_file_sha256': file_hash(ruta_corrida),
        'data_processing_code_sha256': run['code_sha256'],
        'data_processing_pipeline_hashes': run['pipeline_hashes'],
        'input_manifest_sha256': run['manifest_sha256'],
        'plotting_code_sha256': file_hash(ruta_programa),
        'current_analysis_helper_sha256': file_hash(BASE / 'analysis.py'),
        'config': run['config'], 'bins': run['bins'],
    }
    summary = {
        'complete': run['complete'], 'config': run['config'],
        'groups': groups, 'diagnostics': diagnosticos,
        'n_files': len(run['files']), 'provenance': procedencia,
    }
    # Se conservan intactos run.json, los JSON de muestras y su resumen anterior.
    (out / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    (out / 'plot_provenance.json').write_text(json.dumps(procedencia, indent=2), encoding='utf-8')
    resumen_terminal = {}
    for region, resultado in diagnosticos.items():
        resumen_terminal[region] = {}
        for campo, valor in resultado.items():
            if campo != 'shapes':
                resumen_terminal[region][campo] = valor
    print('Gráficos a partir de la corrida registrada:', run['code_sha256'])
    print(json.dumps(resumen_terminal, indent=2))
    print('Salida:', out)


if __name__ == '__main__':
    main()
