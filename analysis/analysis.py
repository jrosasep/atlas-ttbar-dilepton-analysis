"""Análisis de eventos electrón-muón con ATLAS Open Data.

PREGUNTA FÍSICA
    Comparar las distribuciones Delta-phi y |Delta-eta| de leptones
    reconstruidos en colisiones reales con top-antitop y fondos simulados.

RECORRIDO DEL PROGRAMA
    main: elige los archivos y la configuración.
    process_file: abre un ROOT y recorre sus eventos por bloques.
    validate_chunk: comprueba que los datos del bloque sean coherentes.
    event_weights: calcula cuánto aporta cada evento a la predicción.
    select_events: decide qué objetos y eventos cumplen los requisitos.
    observables: calcula las cantidades que llevaremos a los histogramas.
    histogram: acumula conteos, pesos y pesos al cuadrado.

Los gráficos se construyen por separado, en plots_and_summary.py.
"""

# Las anotaciones, como «float» en Config, describen los tipos esperados.
# Esta instrucción permite conservarlas sin evaluarlas inmediatamente.
from __future__ import annotations

# Herramientas incluidas en Python.
import argparse                    # Opciones al ejecutar el programa.
from dataclasses import asdict, dataclass  # Agrupar y registrar parámetros.
import hashlib                     # Huellas para identificar archivos.
import importlib.metadata          # Versiones de las bibliotecas instaladas.
import json                        # Leer manifiestos y guardar resultados.
from pathlib import Path           # Construir rutas de archivos y carpetas.
import time                        # Medir cuánto tarda cada procesamiento.

# Herramientas para los datos y los cálculos numéricos.
import awkward as ak               # Listas de objetos de diferente longitud.
import numpy as np                 # Operaciones numéricas e histogramas.
import uproot                      # Lectura de archivos ROOT.

# Funciones de otros archivos de este mismo proyecto.
from fetch_inputs import verify    # Comprobar la integridad de un ROOT.
from audit_met_inputs import interval_gap  # Comprobar el redondeo del MET.


# ---------------------------------------------------------------------------
# 1. UBICACIONES DE LOS ARCHIVOS
# ---------------------------------------------------------------------------

# __file__ contiene la ubicación de este archivo, analysis.py.
ruta_del_codigo = Path(__file__)

# resolve() obtiene su ubicación absoluta.
ruta_absoluta_del_codigo = ruta_del_codigo.resolve()

# parent obtiene la carpeta que contiene analysis.py.
BASE = ruta_absoluta_del_codigo.parent

# El operador / une las partes de una ruta representada mediante Path.
carpeta_general_de_resultados = BASE / 'results'

# Se conserva la ruta que utilizan los otros programas y los registros actuales.
# Su nombre es una etiqueta de carpeta: no define una selección física distinta.
RESULTS = carpeta_general_de_resultados / 'angular_v2'


# ---------------------------------------------------------------------------
# 2. NOMBRES DE LAS RAMAS QUE LEEREMOS DEL ROOT
# ---------------------------------------------------------------------------

# Estas tuplas contienen NOMBRES, no los valores de los eventos.
# Cada ScaleFactor corrige una parte de la predicción MC. El factor de pileup
# repondera las interacciones simultáneas; otros corrigen eficiencias.
SCALE_FACTORS = (
    'ScaleFactor_PILEUP',
    'ScaleFactor_ELE',
    'ScaleFactor_MUON',
    'ScaleFactor_LepTRIGGER',
    'ScaleFactor_FTAG',
    'ScaleFactor_JVT',
)

# Por evento, cada una de estas ramas contiene una lista de leptones.
# Las posiciones correspondientes describen al mismo leptón.
LEPTON_FIELDS = (
    'lep_pt',              # Momento transversal [GeV, con c = 1].
    'lep_eta',             # Pseudorrapidez, adimensional.
    'lep_phi',             # Ángulo azimutal [rad].
    'lep_e',               # Energía [GeV].
    'lep_type',            # Especie: el código usa 11 para e y 13 para mu.
    'lep_charge',          # Carga: +1 o -1.
    'lep_isTightID',        # Cumple identificación exigente.
    'lep_isTightIso',       # Cumple el requisito de aislamiento.
    'lep_isTrigMatched',    # Está asociado a un objeto del trigger.
)

# Lo mismo para los jets reconstruidos.
JET_FIELDS = (
    'jet_pt',
    'jet_eta',
    'jet_phi',
    'jet_e',
    'jet_jvt',             # Criterio para reducir jets de otras interacciones.
    'jet_btag_quantile',   # Categoría del identificador de jets b.
)

# Ramas con un valor por evento y metadatos de la muestra.
EVENT_FIELDS = (
    'lep_n',
    'jet_n',
    'met',
    'met_phi',
    'met_mpx',
    'met_mpy',
    'mcWeight',
    'trigE',
    'trigM',
    'runNumber',
    'eventNumber',
    'channelNumber',
    'xsec',
    'kfac',
    'filteff',
    'sum_of_weights',
    'num_events',
)

# Construimos una lista en el mismo orden que en el análisis anterior.
BRANCHES = []
for nombre_rama in LEPTON_FIELDS:
    BRANCHES.append(nombre_rama)
for nombre_rama in JET_FIELDS:
    BRANCHES.append(nombre_rama)
for nombre_rama in SCALE_FACTORS:
    BRANCHES.append(nombre_rama)
for nombre_rama in EVENT_FIELDS:
    BRANCHES.append(nombre_rama)


# ---------------------------------------------------------------------------
# 3. BORDES DE LOS HISTOGRAMAS
# ---------------------------------------------------------------------------

# Un histograma con K intervalos necesita K + 1 bordes.
# Se mantienen los bordes del análisis existente.
bordes_pt_electron = np.array([25, 35, 45, 60, 80, 110, 150, 200, 300.])
bordes_pt_muon = np.array([25, 35, 45, 60, 80, 110, 150, 200, 300.])

# arange(inicio, final, paso) no incluye el extremo final.
# Los bordes semienteros permiten contar cantidades enteras de jets.
bordes_numero_jets = np.arange(1.5, 8.5, 1)
bordes_numero_jets_b = np.arange(-.5, 5.5, 1)

# linspace(inicio, final, cantidad) sí incluye ambos extremos.
# Trece bordes producen doce intervalos de Delta-phi, entre 0 y pi.
bordes_delta_phi = np.linspace(0, np.pi, 13)

# Once bordes producen diez intervalos de |Delta-eta|, entre 0 y 5.
bordes_delta_eta = np.linspace(0, 5, 11)

# Variables auxiliares para controles del análisis.
bordes_masa_emu = np.array([0, 25, 50, 75, 100, 125, 150, 200, 300, 450.])
bordes_met = np.array([30, 45, 60, 80, 110, 150, 200, 300.])

BINS = {
    'electron_pt': bordes_pt_electron,
    'muon_pt': bordes_pt_muon,
    'n_jets': bordes_numero_jets,
    'delta_phi': bordes_delta_phi,
    'delta_eta': bordes_delta_eta,
    'm_emu': bordes_masa_emu,
    'met': bordes_met,
    'n_bjets': bordes_numero_jets_b,
}


# ---------------------------------------------------------------------------
# 4. PARÁMETROS DEL ANÁLISIS
# ---------------------------------------------------------------------------

# dataclass agrupa valores y permite crear Config(...).
# frozen=True impide reasignar sus atributos después de crear la configuración.
@dataclass(frozen=True)
class Config:
    """Parámetros de selección y de normalización.

    Una configuración tiene atributos como cfg.lepton_pt.
    Las anotaciones float/int describen el tipo esperado; no son cortes físicos.
    Cambiar un umbral modifica la selección. Cambiar la luminosidad modifica
    la normalización de la simulación.
    """

    lepton_pt: float = 25.0       # Umbral estricto pT > 25 GeV.
    jet_pt: float = 25.0         # Umbral estricto pT > 25 GeV.
    btag_quantile: int = 2       # Se acepta categoría >= 2.
    luminosity_pb: float = 36000.0  # 36 fb^-1 = 36 000 pb^-1, redondeados.


def np1(array):
    """Convertir una columna escalar por evento de Awkward a NumPy.

    Entrada típica: una cantidad por evento, por ejemplo lep_n.
    No se utiliza para forzar listas irregulares de leptones a una matriz.
    La conversión no cambia la cantidad física que contiene la columna.
    """
    columna_numpy = ak.to_numpy(array)
    return columna_numpy


def pipeline_hashes():
    """Identificar el contenido de los tres archivos que afectan al cálculo.

    SHA256 convierte el contenido en una huella de 64 caracteres.
    También cambia si modificamos comentarios; por eso los resultados
    guardados antes de esta reescritura conservan otra huella de código.
    """
    nombres_archivos = ('analysis.py', 'audit_met_inputs.py', 'fetch_inputs.py')
    huellas = {}
    for nombre in nombres_archivos:
        ruta_archivo = BASE / nombre
        contenido = ruta_archivo.read_bytes()
        calculo_sha256 = hashlib.sha256(contenido)
        huella = calculo_sha256.hexdigest()
        huellas[nombre] = huella
    return huellas


# ---------------------------------------------------------------------------
# 5. COMPROBAR LOS DATOS LEÍDOS DEL ROOT
# ---------------------------------------------------------------------------

def validate_chunk(a, sample):
    """Comprobar un bloque de eventos antes de hacer el análisis físico.

    a: bloque leído mediante Uproot, representado como Awkward Array.
    sample: descripción de su archivo, procedencia y metadatos.
    Devuelve el mayor residuo entre MET y el módulo de sus componentes.

    Estas comprobaciones son controles de entrada, no filtros escogidos
    para mejorar el acuerdo de las figuras. Una incoherencia detiene la corrida.
    """

    # Cada lista de propiedades debe tener el número de objetos anunciado.
    grupos_objetos = ((LEPTON_FIELDS, 'lep_n'), (JET_FIELDS, 'jet_n'))
    for campos, nombre_conteo in grupos_objetos:
        for campo in campos:
            lista_objetos = a[campo]
            longitudes = ak.num(lista_objetos, axis=1)
            conteos_declarados = a[nombre_conteo]
            longitudes_correctas = longitudes == conteos_declarados
            todas_correctas = ak.all(longitudes_correctas)
            if not todas_correctas:
                raise ValueError(f'Longitud incompatible en {campo}')

            valores_finitos = np.isfinite(lista_objetos)
            todos_finitos = ak.all(valores_finitos)
            if not todos_finitos:
                raise ValueError(f'Valor no finito en {campo}')

    # Para las demás ramas hay un valor por evento.
    campos_objetos = LEPTON_FIELDS + JET_FIELDS
    nombres_campos_objetos = set(campos_objetos)
    nombres_todas_ramas = set(BRANCHES)
    campos_escalares = nombres_todas_ramas - nombres_campos_objetos
    for campo in campos_escalares:
        valores = np1(a[campo])
        valores_finitos = np.isfinite(valores)
        if not np.all(valores_finitos):
            raise ValueError(f'Valor no finito en {campo}')

    # Rangos básicos de momentos, ángulos, categorías y cargas.
    pt_lepton_valido = a.lep_pt >= 0
    phi_lepton_valido = abs(a.lep_phi) <= np.pi + 1e-5
    leptones_validos = pt_lepton_valido & phi_lepton_valido
    if not ak.all(leptones_validos):
        raise ValueError('Momento o rango angular inválido')

    pt_jet_valido = a.jet_pt >= 0
    categoria_minima = a.jet_btag_quantile >= 1
    categoria_maxima = a.jet_btag_quantile <= 5
    jets_validos = pt_jet_valido & categoria_minima & categoria_maxima
    if not ak.all(jets_validos):
        raise ValueError('pT de jet o categoría b-tag inválida')

    carga_positiva = a.lep_charge == 1
    carga_negativa = a.lep_charge == -1
    cargas_validas = carga_positiva | carga_negativa
    if not ak.all(cargas_validas):
        raise ValueError('Carga leptónica inesperada')

    # Comprobar el filtro previo 2J2LMET30 de ATLAS.
    # 1e-5 reconoce el redondeo cerca de los umbrales; no añade otro corte.
    jets_del_skim = a.jet_pt >= 20 - 1e-5
    numero_jets_skim = ak.sum(jets_del_skim, axis=1)
    dos_jets_skim = numero_jets_skim >= 2

    pt_leptones_skim = a.lep_pt >= 7 - 1e-5
    leptones_del_skim = pt_leptones_skim & a.lep_isTightID
    numero_leptones_skim = ak.sum(leptones_del_skim, axis=1)
    dos_leptones_skim = numero_leptones_skim >= 2
    met_del_skim = a.met >= 30 - 1e-5

    cumple_skim = dos_jets_skim & dos_leptones_skim & met_del_skim
    if not ak.all(cumple_skim):
        raise ValueError('Eventos incompatibles con el skim público 2J2LMET30')

    # MET es el módulo de un vector transversal con componentes mpx y mpy.
    met_original = np1(a.met)
    met_px_original = np1(a.met_mpx)
    met_py_original = np1(a.met_mpy)
    if len(a) > 0:
        met_px_double = met_px_original.astype(float)
        met_py_double = met_py_original.astype(float)
        met_double = met_original.astype(float)
        modulo_calculado = np.hypot(met_px_double, met_py_double)
        diferencias = np.abs(modulo_calculado - met_double)
        mayor_residuo = np.max(diferencias)
    else:
        mayor_residuo = 0.0

    # El control interval_gap reconoce que las columnas originales son float32.
    # Se conserva la tolerancia de 0.01 GeV del análisis anterior.
    separacion_intervalos, tolerancia_aritmetica = interval_gap(
        met_original, met_px_original, met_py_original
    )
    limite_permitido = .01 + tolerancia_aritmetica
    met_incoherente = separacion_intervalos > limite_permitido
    if np.any(met_incoherente):
        raise ValueError('MET incompatible con sus componentes')

    # Los metadatos MC del ROOT deben coincidir con el catálogo utilizado.
    if sample['group'] != 'data':
        correspondencia = {
            'xsec': 'crossSection_pb',
            'filteff': 'genFiltEff',
            'kfac': 'kFactor',
            'sum_of_weights': 'sumOfWeights',
            'num_events': 'nEvents',
        }
        for rama, columna in correspondencia.items():
            valor_catalogo = float(sample['metadata'][columna])
            valores_root = np1(a[rama])
            coinciden = np.allclose(
                valores_root, valor_catalogo, rtol=5e-5, atol=1e-8
            )
            if not coinciden:
                raise ValueError(
                    f'Metadato {rama} contradice CSV oficial para {sample["dsid"]}'
                )

        identificadores_mc = np1(a.channelNumber)
        identificador_esperado = sample['dsid']
        identificadores_correctos = identificadores_mc == identificador_esperado
        if not np.all(identificadores_correctos):
            raise ValueError('Identificador MC distinto al del manifiesto')

    return float(mayor_residuo)


# ---------------------------------------------------------------------------
# 6. PESOS: CUÁNTO APORTA UN EVENTO A LA PREDICCIÓN
# ---------------------------------------------------------------------------

def event_weights(a, sample, cfg=Config()):
    """Calcular un peso por evento, siguiendo la normalización de la Pauta 3.

    Datos reales: cada evento tiene peso 1.
    MC: L * sigma * eficiencia_filtro * k / suma_pesos_original,
    multiplicado por mcWeight y los factores de corrección disponibles.

    La suma del denominador pertenece a la producción anterior al skim.
    No se utiliza la suma posterior a nuestros cortes ni se ajusta a los datos.
    Se conservan los pesos negativos y los pesos cero.
    """
    if sample['group'] == 'data':
        pesos_datos = np.ones(len(a), dtype=np.float64)
        return pesos_datos

    metadatos = sample['metadata']
    suma_pesos_original = float(metadatos['sumOfWeights'])
    if suma_pesos_original <= 0:
        raise ValueError('Normalización original no positiva')

    luminosidad = cfg.luminosity_pb
    seccion_eficaz = float(metadatos['crossSection_pb'])
    eficiencia_filtro = float(metadatos['genFiltEff'])
    factor_k = float(metadatos['kFactor'])

    # Se conserva el orden de las multiplicaciones del código anterior.
    normalizacion = luminosidad * seccion_eficaz
    normalizacion = normalizacion * eficiencia_filtro
    normalizacion = normalizacion * factor_k
    normalizacion = normalizacion / suma_pesos_original

    pesos_generador = np1(a.mcWeight)
    pesos_generador = pesos_generador.astype(np.float64)
    pesos = pesos_generador * normalizacion

    for nombre_factor in SCALE_FACTORS:
        factor = np1(a[nombre_factor])
        factor = factor.astype(np.float64)
        pesos = pesos * factor

    pesos_finitos = np.isfinite(pesos)
    if not np.all(pesos_finitos):
        raise ValueError('Peso final no finito: el análisis debe detenerse')
    return pesos


# ---------------------------------------------------------------------------
# 7. SELECCIÓN DE OBJETOS Y EVENTOS
# ---------------------------------------------------------------------------

def select_events(a, cfg=Config()):
    """Construir máscaras de objetos y luego de eventos completos.

    Una máscara es un conjunto de valores True/False que marca qué se acepta.
    Las máscaras de objetos tienen una lista por evento; las de eventos tienen
    un solo True/False por evento.

    Se exige exactamente un electrón y un muón entre los objetos aceptados.
    Otros objetos que fallen los requisitos no se vetan por su sola presencia.
    SR2b es la región principal; SR1b es auxiliar y contiene a SR2b.
    SS1b es un control de cargas iguales, no un fondo restado automáticamente.

    Salida, en el mismo orden usado por el resto del proyecto:
    etapas, regiones, máscara electrones, máscara muones, número jets, número b.
    """
    eta_absoluta = abs(a.lep_eta)

    # Calidad común de los leptones.
    cumple_pt_lepton = a.lep_pt > cfg.lepton_pt
    cumple_identificacion = a.lep_isTightID
    cumple_aislamiento = a.lep_isTightIso
    calidad = cumple_pt_lepton & cumple_identificacion
    calidad = calidad & cumple_aislamiento

    # Electrón: tipo 11, aceptación y exclusión de la transición del calorímetro.
    es_electron = abs(a.lep_type) == 11
    electron_en_barril = eta_absoluta < 1.37
    fuera_transicion = eta_absoluta > 1.52
    dentro_extremo_electron = eta_absoluta < 2.47
    electron_en_extremo = fuera_transicion & dentro_extremo_electron
    aceptacion_electron = electron_en_barril | electron_en_extremo
    electron_base = calidad & es_electron
    electron_base = electron_base & aceptacion_electron
    electron = electron_base & a.lep_isTrigMatched

    # Muón: tipo 13 y aceptación |eta| < 2.5.
    es_muon = abs(a.lep_type) == 13
    aceptacion_muon = eta_absoluta < 2.5
    muon_base = calidad & es_muon
    muon_base = muon_base & aceptacion_muon
    muon = muon_base & a.lep_isTrigMatched

    # Jets aceptados y subconjunto con etiqueta b.
    cumple_pt_jet = a.jet_pt > cfg.jet_pt
    cumple_eta_jet = abs(a.jet_eta) < 2.5
    jet = cumple_pt_jet & cumple_eta_jet
    jet = jet & a.jet_jvt
    cumple_etiqueta_b = a.jet_btag_quantile >= cfg.btag_quantile
    bjet = jet & cumple_etiqueta_b

    # axis=1 suma dentro de la lista de objetos de cada evento.
    # True aporta 1 y False aporta 0: obtenemos cuántos objetos se aceptaron.
    conteo_jets = ak.sum(jet, axis=1)
    conteo_jets_b = ak.sum(bjet, axis=1)
    nj = np1(conteo_jets)
    nb = np1(conteo_jets_b)

    trigger_aceptado = a.trigE | a.trigM
    trigger = np1(trigger_aceptado)
    numero_electrones = ak.sum(electron, axis=1)
    numero_muones = ak.sum(muon, axis=1)
    exactamente_un_electron = numero_electrones == 1
    exactamente_un_muon = numero_muones == 1
    pareja_aceptada = exactamente_un_electron & exactamente_un_muon
    pair = np1(pareja_aceptada)

    # Las propiedades enmascaradas conservan la carga del leptón aceptado.
    cargas_electrones = a.lep_charge * electron
    cargas_muones = a.lep_charge * muon
    suma_carga_electron = ak.sum(cargas_electrones, axis=1)
    suma_carga_muon = ak.sum(cargas_muones, axis=1)
    carga_electron = np1(suma_carga_electron)
    carga_muon = np1(suma_carga_muon)
    producto_cargas = carga_electron * carga_muon

    pair_trigger = trigger & pair
    cargas_opuestas = producto_cargas < 0
    opposite = pair_trigger & cargas_opuestas
    pretag = opposite & (nj >= 2)
    una_etiqueta_b = pretag & (nb >= 1)
    dos_etiquetas_b = pretag & (nb >= 2)

    # Estas etapas son acumulativas. Permiten construir la tabla cutflow.
    stages = {
        'input': np.ones(len(a), dtype=bool),
        'trigger': trigger,
        'one_e_one_mu_matched': pair_trigger,
        'opposite_sign': opposite,
        'two_jets_jvt': pretag,
        'one_bjet': una_etiqueta_b,
        'two_bjets': dos_etiquetas_b,
    }

    cargas_iguales = producto_cargas > 0
    control_cargas_iguales = pair_trigger & cargas_iguales
    control_cargas_iguales = control_cargas_iguales & (nj >= 2)
    control_cargas_iguales = control_cargas_iguales & (nb >= 1)
    regions = {
        'SR1b': stages['one_bjet'],
        'SR2b': stages['two_bjets'],
        'PRETAG': pretag,
        'ZEROb': pretag & (nb == 0),
        'SS1b': control_cargas_iguales,
    }
    return stages, regions, electron, muon, nj, nb


# ---------------------------------------------------------------------------
# 8. OBSERVABLES DE LOS OBJETOS ACEPTADOS
# ---------------------------------------------------------------------------

def observables(a, mask, electron, muon, nj, nb):
    """Calcular cantidades del laboratorio para una región seleccionada.

    mask marca los eventos de esa región; electron/muon marcan sus objetos.
    Salida: diccionario de valores y diagnóstico de redondeo de la masa auxiliar.
    Delta-phi está entre 0 y pi; |Delta-eta| es adimensional.
    m_emu es la masa del par de leptones, no la masa del sistema top-antitop.
    """
    values = {}
    objetos = (('electron', electron), ('muon', muon))
    for etiqueta, mascara_objetos in objetos:
        for propiedad in ('pt', 'eta', 'phi', 'e'):
            nombre_rama = 'lep_' + propiedad
            listas_del_bloque = a[nombre_rama]
            listas_objetos_aceptados = listas_del_bloque[mascara_objetos]
            listas_eventos_aceptados = listas_objetos_aceptados[mask]

            # Hay exactamente un objeto de cada especie en estos eventos.
            # firsts toma ese valor único de cada lista y conserva el orden.
            valor_por_evento = ak.firsts(listas_eventos_aceptados)
            columna = np1(valor_por_evento)
            columna = columna.astype(float)
            nombre_salida = etiqueta + '_' + propiedad
            values[nombre_salida] = columna

    phi_electron = values['electron_phi']
    phi_muon = values['muon_phi']
    diferencia_phi = phi_electron - phi_muon
    seno_diferencia = np.sin(diferencia_phi)
    coseno_diferencia = np.cos(diferencia_phi)
    diferencia_periodica = np.arctan2(seno_diferencia, coseno_diferencia)
    delta_phi = np.abs(diferencia_periodica)

    eta_electron = values['electron_eta']
    eta_muon = values['muon_eta']
    diferencia_eta = eta_electron - eta_muon
    delta_eta = np.abs(diferencia_eta)

    # Control auxiliar: sumar los momentos tridimensionales de ambos leptones.
    px = 0
    py = 0
    pz = 0
    for etiqueta in ('electron', 'muon'):
        momento_transversal = values[etiqueta + '_pt']
        azimut = values[etiqueta + '_phi']
        pseudorrapidez = values[etiqueta + '_eta']
        componente_x = momento_transversal * np.cos(azimut)
        componente_y = momento_transversal * np.sin(azimut)
        componente_z = momento_transversal * np.sinh(pseudorrapidez)
        px = px + componente_x
        py = py + componente_y
        pz = pz + componente_z

    energia_total = values['electron_e'] + values['muon_e']
    energia_cuadrada = energia_total**2
    momento_cuadrado = px**2 + py**2 + pz**2
    masa_cuadrada = energia_cuadrada - momento_cuadrado
    escala_redondeo = energia_cuadrada + momento_cuadrado
    masa, redondeo = checked_mass(masa_cuadrada, escala_redondeo)

    met_del_bloque = np1(a.met)
    resultado = {
        'electron_pt': values['electron_pt'],
        'muon_pt': values['muon_pt'],
        'delta_phi': delta_phi,
        'delta_eta': delta_eta,
        'm_emu': masa,
        'met': met_del_bloque[mask],
        'n_jets': nj[mask],
        'n_bjets': nb[mask],
    }
    return resultado, redondeo


def checked_mass(m2, scale):
    """Controlar el redondeo de m² = E² - |p|² antes de tomar la raíz.

    Las entradas originales son float32. Restar números grandes y cercanos
    puede dar una negatividad pequeña por precisión numérica. Sólo se admite
    la tolerancia ya fijada: 32 * epsilon_float32 * max(escala, 1).
    Esos casos se registran y reciben masa cero; otros detienen el análisis.
    """
    precision_float32 = np.finfo(np.float32).eps
    factor_tolerancia = 32 * precision_float32
    escala_controlada = np.maximum(scale, 1.)
    tolerancia = factor_tolerancia * escala_controlada
    demasiado_negativa = m2 < -tolerancia
    if np.any(demasiado_negativa):
        raise ValueError('Masa invariante inconsistente más allá de la precisión float32')

    negativos = m2 < 0
    cantidad_negativos = int(negativos.sum())
    if len(m2) > 0:
        menor_masa_cuadrada = float(min(0., np.min(m2)))
    else:
        menor_masa_cuadrada = 0.0
    if negativos.any():
        negatividades_relativas = -m2[negativos] / scale[negativos]
        mayor_negatividad_relativa = float(np.max(negatividades_relativas))
    else:
        mayor_negatividad_relativa = 0.0

    diagnostic = {
        'count': cantidad_negativos,
        'minimum_m2_gev2': menor_masa_cuadrada,
        'max_relative_negative': mayor_negatividad_relativa,
    }
    masa_cuadrada_controlada = np.maximum(m2, 0)
    masa = np.sqrt(masa_cuadrada_controlada)
    return masa, diagnostic


# ---------------------------------------------------------------------------
# 9. HISTOGRAMAS: CONTEOS, PREDICCIONES Y VARIANZAS
# ---------------------------------------------------------------------------

def histogram(values, weights, edges):
    """Acumular un observable en los intervalos definidos por edges.

    values: un valor del observable por evento seleccionado.
    weights: su contribución; 1 en datos y peso calculado en MC.
    edges: bordes, con un elemento más que el número de intervalos.

    sumw es suma(w); sumw2 es suma(w²); raw es el número sin ponderar.
    Los desbordes se registran y se incluyen en los intervalos extremos.
    """
    limite_inferior = edges[0]
    limite_superior = edges[-1]

    # nextafter toma el número representable inmediatamente anterior al borde
    # final. Así los valores de overflow se acumulan en el último intervalo.
    ultimo_valor_interior = np.nextafter(limite_superior, limite_inferior)
    valores_para_histograma = np.clip(
        values, limite_inferior, ultimo_valor_interior
    )

    histograma_pesos, _ = np.histogram(
        valores_para_histograma, edges, weights=weights
    )
    pesos_cuadrados = weights**2
    histograma_varianza, _ = np.histogram(
        valores_para_histograma, edges, weights=pesos_cuadrados
    )
    histograma_conteos, _ = np.histogram(valores_para_histograma, edges)

    debajo_rango = values < limite_inferior
    encima_rango = values >= limite_superior
    resultado = {
        'sumw': histograma_pesos.tolist(),
        'sumw2': histograma_varianza.tolist(),
        'raw': histograma_conteos.tolist(),
        'underflow_raw': int(np.sum(debajo_rango)),
        'overflow_raw': int(np.sum(encima_rango)),
    }
    return resultado


def add_hist(target, update):
    """Añadir al histograma acumulado la contribución del bloque actual.

    target se modifica; update no se modifica. Se suman los intervalos
    correspondientes sin normalizar cada bloque por separado.
    """
    for campo in ('sumw', 'sumw2', 'raw'):
        valores_acumulados = np.array(target[campo])
        valores_del_bloque = np.array(update[campo])
        nueva_suma = valores_acumulados + valores_del_bloque
        target[campo] = nueva_suma.tolist()
    for campo in ('underflow_raw', 'overflow_raw'):
        target[campo] = target[campo] + update[campo]


def normalized_shape(s, v):
    """Calcular fracciones normalizadas y su matriz de covarianza.

    s: conteo o suma de pesos de cada intervalo.
    v: varianza estadística de cada intervalo.
    p_i = s_i / S, con S = suma(s).
    J_ij = delta_ij / S - s_i / S².
    C = J diag(v) J^T: también propaga la fluctuación del total.
    """
    sumas = np.asarray(s, float)
    varianzas = np.asarray(v, float)
    total = sumas.sum()
    if total <= 0:
        raise ValueError('Área no positiva en histograma normalizado')

    numero_intervalos = len(sumas)
    matriz_identidad = np.eye(numero_intervalos)
    primer_termino = matriz_identidad / total

    # [:, None] representa las sumas como una columna. NumPy repite esta
    # columna al operar con la matriz: es el término común de cada fila.
    columna_sumas = sumas[:, None]
    total_cuadrado = total**2
    segundo_termino = columna_sumas / total_cuadrado
    jacobiano = primer_termino - segundo_termino

    # [None, :] representa las varianzas como una fila.
    fila_varianzas = varianzas[None, :]
    jacobiano_ponderado = jacobiano * fila_varianzas
    jacobiano_transpuesto = jacobiano.T
    covarianza = jacobiano_ponderado @ jacobiano_transpuesto
    fracciones = sumas / total
    return fracciones, covarianza


# ---------------------------------------------------------------------------
# 10. PROCESAR UN ARCHIVO COMPLETO POR BLOQUES
# ---------------------------------------------------------------------------

def process_file(sample, cfg, chunk_size=100000):
    """Leer un ROOT, aplicar el análisis y devolver su registro numérico.

    sample describe el archivo; cfg contiene sus parámetros.
    chunk_size es el máximo de eventos leído por bloque.
    El diccionario devuelto conserva conteos, histogramas y diagnósticos.
    """
    inicio = time.monotonic()

    # A. Ubicar el ROOT y comprobar su integridad.
    carpeta_datos = BASE / 'data'
    nombre_archivo = sample['key']
    path = carpeta_datos / nombre_archivo
    huella_entrada = verify(path, sample)
    archivo_verificado = huella_entrada is not None
    if 'sha256' in sample:
        misma_huella = huella_entrada == sample['sha256']
    else:
        misma_huella = True
    if not archivo_verificado or not misma_huella:
        raise ValueError('El archivo cambió después de la verificación de descarga')

    # B. Abrir el archivo y recuperar el árbol llamado analysis.
    archivo_root = uproot.open(path)
    tree = archivo_root['analysis']
    nombres_disponibles = set(tree.keys())
    nombres_necesarios = set(BRANCHES)
    nombres_ausentes = nombres_necesarios - nombres_disponibles
    missing = sorted(nombres_ausentes)
    if missing:
        raise ValueError(f'Ramas requeridas ausentes: {missing}')

    # C. Preparar contadores e histogramas vacíos.
    result = {
        'file': nombre_archivo,
        'dsid': sample['dsid'],
        'group': sample['group'],
        'entries': tree.num_entries,
        'branches': tree.typenames(),
        'cutflow': {},
        'histograms': {},
        'regions': {},
        'processed': 0,
        'negative_generator_weights': 0,
        'negative_final_weights': 0,
        'zero_final_weights': 0,
        'max_met_residual_gev': 0.,
        'max_abs_weight': 0.,
        'metadata': sample['metadata'],
        'mass_rounding': {},
        'extreme_met': {
            'input_over_13tev': 0,
            'max_input_gev': 0.,
            'regions_over_13tev': {},
        },
    }
    result['input_sha256'] = huella_entrada
    nombres_etapas = (
        'input', 'trigger', 'one_e_one_mu_matched', 'opposite_sign',
        'two_jets_jvt', 'one_bjet', 'two_bjets'
    )
    for etapa in nombres_etapas:
        result['cutflow'][etapa] = {'raw': 0, 'sumw': 0., 'sumw2': 0.}

    nombres_regiones = ('SR1b', 'SR2b', 'PRETAG', 'ZEROb', 'SS1b')
    for region in nombres_regiones:
        result['regions'][region] = {'raw': 0, 'sumw': 0., 'sumw2': 0.}
        result['histograms'][region] = {}
        for observable, bordes in BINS.items():
            valores_vacios = np.array([])
            pesos_vacios = np.array([])
            histograma_vacio = histogram(valores_vacios, pesos_vacios, bordes)
            result['histograms'][region][observable] = histograma_vacio
        result['mass_rounding'][region] = {
            'count': 0,
            'minimum_m2_gev2': 0.,
            'max_relative_negative': 0.,
        }
        result['extreme_met']['regions_over_13tev'][region] = 0

    identificadores_por_bloque = []

    # D. Leer las ramas solicitadas en bloques representados por Awkward.
    bloques = tree.iterate(
        expressions=BRANCHES,
        step_size=chunk_size,
        library='ak',
    )
    for a in bloques:
        # El recorrido físico de cada bloque empieza aquí.
        residuo_met = validate_chunk(a, sample)
        pesos = event_weights(a, sample, cfg)
        stages, regions, electron, muon, nj, nb = select_events(a, cfg)

        # Diagnóstico de cantidades reconstruidas anómalas: no añade un corte.
        met_bloque = np1(a.met)
        met_extremo = met_bloque > 13000.
        cantidad_extremos = int(met_extremo.sum())
        result['extreme_met']['input_over_13tev'] += cantidad_extremos
        mayor_met_bloque = float(np.max(met_bloque))
        mayor_met_anterior = result['extreme_met']['max_input_gev']
        result['extreme_met']['max_input_gev'] = max(
            mayor_met_anterior, mayor_met_bloque
        )

        # E. Acumular el cutflow. Los conjuntos deben estar anidados.
        # Los rendimientos ponderados pueden no disminuir con pesos negativos.
        mascara_anterior = np.ones(len(a), dtype=bool)
        for etiqueta, mascara in stages.items():
            eventos_fuera_etapa_anterior = mascara & ~mascara_anterior
            if np.any(eventos_fuera_etapa_anterior):
                raise AssertionError('Cutflow no acumulativo')
            mascara_anterior = mascara
            valores_iniciales = {'raw': 0, 'sumw': 0., 'sumw2': 0.}
            contador = result['cutflow'].setdefault(etiqueta, valores_iniciales)
            pesos_aceptados = pesos[mascara]
            contador['raw'] += int(mascara.sum())
            contador['sumw'] += float(pesos_aceptados.sum())
            pesos_cuadrados = pesos_aceptados**2
            contador['sumw2'] += float(np.sum(pesos_cuadrados))

        # F. Calcular observables y llenar histogramas para cada región.
        for region, mascara in regions.items():
            extremos_aceptados = met_extremo & mascara
            result['extreme_met']['regions_over_13tev'][region] += int(
                np.sum(extremos_aceptados)
            )
            valores_iniciales = {'raw': 0, 'sumw': 0., 'sumw2': 0.}
            contador = result['regions'].setdefault(region, valores_iniciales)
            pesos_aceptados = pesos[mascara]
            contador['raw'] += int(mascara.sum())
            contador['sumw'] += float(pesos_aceptados.sum())
            pesos_cuadrados = pesos_aceptados**2
            contador['sumw2'] += float(np.sum(pesos_cuadrados))

            histogramas_region = result['histograms'].setdefault(region, {})
            cantidades, redondeo = observables(a, mascara, electron, muon, nj, nb)
            registro_redondeo = result['mass_rounding'][region]
            registro_redondeo['count'] += redondeo['count']
            registro_redondeo['minimum_m2_gev2'] = min(
                registro_redondeo['minimum_m2_gev2'],
                redondeo['minimum_m2_gev2'],
            )
            registro_redondeo['max_relative_negative'] = max(
                registro_redondeo['max_relative_negative'],
                redondeo['max_relative_negative'],
            )

            for nombre, bordes in BINS.items():
                valores = cantidades[nombre]
                contribucion_bloque = histogram(valores, pesos_aceptados, bordes)
                if nombre in histogramas_region:
                    add_hist(histogramas_region[nombre], contribucion_bloque)
                else:
                    histogramas_region[nombre] = contribucion_bloque

        # G. Registrar controles del bloque y sus identificadores de datos.
        result['processed'] += len(a)
        pesos_generador = np1(a.mcWeight)
        result['negative_generator_weights'] += int(np.sum(pesos_generador < 0))
        result['negative_final_weights'] += int(np.sum(pesos < 0))
        result['zero_final_weights'] += int(np.sum(pesos == 0))
        mayor_peso_bloque = float(np.max(np.abs(pesos)))
        result['max_abs_weight'] = max(result['max_abs_weight'], mayor_peso_bloque)
        result['max_met_residual_gev'] = max(
            result['max_met_residual_gev'], residuo_met
        )
        if sample['group'] == 'data':
            # Un identificador consta de runNumber y eventNumber.
            tipo_identificador = [('run', '<u4'), ('event', '<u8')]
            identificadores = np.empty(len(a), dtype=tipo_identificador)
            identificadores['run'] = np1(a.runNumber)
            identificadores['event'] = np1(a.eventNumber)
            identificadores_por_bloque.append(identificadores)

    # H. Comprobar que el resultado conserva los totales del archivo.
    if result['processed'] != tree.num_entries:
        raise AssertionError('No se leyó el archivo completo')
    for region, histogramas_region in result['histograms'].items():
        totales_region = result['regions'][region]
        for histograma in histogramas_region.values():
            conteo_histograma = sum(histograma['raw'])
            peso_histograma = sum(histograma['sumw'])
            varianza_histograma = sum(histograma['sumw2'])
            assert conteo_histograma == totales_region['raw']
            assert np.isclose(peso_histograma, totales_region['sumw'])
            assert np.isclose(varianza_histograma, totales_region['sumw2'])

    if identificadores_por_bloque:
        identificadores_archivo = np.concatenate(identificadores_por_bloque)
        identificadores_unicos = np.unique(identificadores_archivo)
        if len(identificadores_unicos) != len(identificadores_archivo):
            raise ValueError('Identificadores repetidos en archivo de colisiones')
        carpeta_identificadores = RESULTS / 'event_ids'
        nombre_identificadores = path.stem + '.npy'
        ruta_identificadores = carpeta_identificadores / nombre_identificadores
        np.save(ruta_identificadores, identificadores_archivo)

    duracion = time.monotonic() - inicio
    result['seconds'] = round(duracion, 2)
    return result


# ---------------------------------------------------------------------------
# 11. COORDINAR LA EJECUCIÓN Y GUARDAR LOS RESULTADOS
# ---------------------------------------------------------------------------

def main():
    """Elegir las muestras, procesarlas y guardar la descripción de la corrida.

    main se ejecuta desde el bloque final del archivo.
    Los datos ROOT se buscan en BASE/data; no se descargan automáticamente.
    --resume sólo admite resultados de este mismo código y configuración.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--only-dsid', type=int)
    parser.add_argument('--chunk-size', type=int, default=100000)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument(
        '--available-only',
        action='store_true',
        help='Procesa archivos ya descargados; marca salida como parcial',
    )
    args = parser.parse_args()

    # Elegir el manifiesto de archivos verificados.
    carpeta_procedencia = BASE / 'provenance'
    manifest_path = carpeta_procedencia / 'verified_inputs.json'
    if not manifest_path.exists():
        if not args.available_only:
            raise ValueError('La descarga completa aún no está verificada')
        manifest_path = BASE / 'manifest.json'
    texto_manifiesto = manifest_path.read_text(encoding='utf-8')
    manifest = json.loads(texto_manifiesto)

    carpeta_muestras = RESULTS / 'samples'
    carpeta_muestras.mkdir(parents=True, exist_ok=True)
    carpeta_identificadores = RESULTS / 'event_ids'
    carpeta_identificadores.mkdir(exist_ok=True)

    cfg = Config(luminosity_pb=manifest['luminosity_pb'])
    configuracion_registrada = asdict(cfg)
    contenido_codigo = ruta_del_codigo.read_bytes()
    huella_codigo = hashlib.sha256(contenido_codigo)
    code_sha = huella_codigo.hexdigest()
    hashes = pipeline_hashes()
    contenido_manifiesto = manifest_path.read_bytes()
    huella_manifiesto = hashlib.sha256(contenido_manifiesto)

    versiones = {}
    for paquete in ('numpy', 'uproot', 'awkward', 'matplotlib', 'scipy'):
        versiones[paquete] = importlib.metadata.version(paquete)
    bordes_registrados = {}
    for nombre, bordes in BINS.items():
        bordes_registrados[nombre] = bordes.tolist()

    run = {
        'config': configuracion_registrada,
        'code_sha256': code_sha,
        'pipeline_hashes': hashes,
        'manifest_sha256': huella_manifiesto.hexdigest(),
        'versions': versiones,
        'bins': bordes_registrados,
        'files': [],
    }

    for sample in manifest['files']:
        ruta_root = BASE / 'data' / sample['key']
        if args.available_only and not ruta_root.exists():
            continue
        if args.only_dsid is not None and sample['dsid'] != args.only_dsid:
            continue

        nombre_resultado = sample['key'] + '.json'
        dest = carpeta_muestras / nombre_resultado
        if args.resume and dest.exists():
            # No se cambia la procedencia de una corrida guardada para hacerla
            # pasar por una ejecución de un código modificado.
            texto_resultado = dest.read_text()
            saved = json.loads(texto_resultado)
            mismo_codigo = saved.get('code_sha256') == code_sha
            misma_configuracion = saved.get('config') == configuracion_registrada
            mismos_controles = saved.get('pipeline_hashes') == hashes
            if not mismo_codigo or not misma_configuracion or not mismos_controles:
                raise ValueError('Caché de una configuración/código diferente; volver a ejecutar')
            if saved.get('metadata') != sample['metadata']:
                raise ValueError('Los metadatos cambiaron respecto a la corrida guardada')
            if sample.get('sha256'):
                misma_entrada = saved.get('input_sha256') == sample['sha256']
                if not misma_entrada:
                    raise ValueError('El archivo verificado no es el usado para el resultado guardado')
        else:
            saved = process_file(sample, cfg, args.chunk_size)
            saved['code_sha256'] = code_sha
            saved['config'] = configuracion_registrada
            saved['pipeline_hashes'] = hashes
            texto_resultado = json.dumps(saved, indent=2)
            dest.write_text(texto_resultado, encoding='utf-8')

        ruta_relativa = dest.relative_to(BASE)
        run['files'].append(str(ruta_relativa))
        etiqueta_archivo = sample['dsid'] or sample['key']
        # Se mantiene el mensaje de la ejecución anterior: muestra SR1b auxiliar.
        print(etiqueta_archivo, saved['processed'], saved['regions']['SR1b'], flush=True)

    run['complete'] = len(run['files']) == len(manifest['files'])
    if run['complete']:
        # Cada colisión real debe aparecer una sola vez también ENTRE archivos.
        rutas_identificadores = sorted(carpeta_identificadores.glob('*.npy'))
        listas_identificadores = []
        for ruta in rutas_identificadores:
            identificadores_archivo = np.load(ruta)
            listas_identificadores.append(identificadores_archivo)
        todos_identificadores = np.concatenate(listas_identificadores)
        identificadores_unicos = np.unique(todos_identificadores)
        if len(identificadores_unicos) != len(todos_identificadores):
            raise ValueError('Eventos de datos duplicados ENTRE archivos')
        cantidad_datos = len(todos_identificadores)
        if cantidad_datos != 6242521:
            raise ValueError('Conteo global de datos distinto al catálogo')
        run['global_checks'] = {
            'unique_data_events': cantidad_datos,
            'expected_data_events': 6242521,
            'all_checks_passed': True,
        }
        print('CONTROL GLOBAL: 6242521 eventos reales únicos', flush=True)

    if run['complete']:
        nombre_corrida = 'run.json'
    else:
        nombre_corrida = 'run_partial.json'
    destination = RESULTS / nombre_corrida
    temporary = destination.with_suffix('.json.part')
    texto_corrida = json.dumps(run, indent=2)
    temporary.write_text(texto_corrida, encoding='utf-8')
    temporary.replace(destination)


# Python define primero las funciones de arriba. Cuando ejecutamos este .py
# directamente, entra en main(). Importarlo desde otro programa no llama main.
if __name__ == '__main__':
    main()
