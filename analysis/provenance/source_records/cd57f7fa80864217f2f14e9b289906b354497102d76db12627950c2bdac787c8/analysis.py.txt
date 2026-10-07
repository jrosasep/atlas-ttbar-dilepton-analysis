"""Análisis reproducible e-mu de colisiones REALES y Monte Carlo de ATLAS.

Orden de lectura recomendado: select_events -> event_weights -> process_file.
El diseño y la selección básica parten de open_data.py del repositorio privado;
se agregan los triggers, el matching y JVT del ejemplo oficial TTbarDilepAnalysis.
Los números se calculan con los archivos ROOT descargados, nunca con las figuras
o el resumen guardado del trabajo anterior.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time

import awkward as ak
import numpy as np
import uproot
from fetch_inputs import verify
from audit_met_inputs import interval_gap

BASE = Path(__file__).resolve().parent
# La versión angular escribe aparte para conservar la corrida anterior intacta.
RESULTS = BASE / 'results' / 'angular_v2'
SCALE_FACTORS = ('ScaleFactor_PILEUP', 'ScaleFactor_ELE', 'ScaleFactor_MUON',
                 'ScaleFactor_LepTRIGGER', 'ScaleFactor_FTAG', 'ScaleFactor_JVT')
LEPTON_FIELDS = ('lep_pt', 'lep_eta', 'lep_phi', 'lep_e', 'lep_type', 'lep_charge',
                'lep_isTightID', 'lep_isTightIso', 'lep_isTrigMatched')
JET_FIELDS = ('jet_pt', 'jet_eta', 'jet_phi', 'jet_e', 'jet_jvt', 'jet_btag_quantile')
BRANCHES = list(LEPTON_FIELDS + JET_FIELDS + SCALE_FACTORS) + [
    'lep_n', 'jet_n', 'met', 'met_phi', 'met_mpx', 'met_mpy', 'mcWeight',
    'trigE', 'trigM', 'runNumber', 'eventNumber', 'channelNumber', 'xsec',
    'kfac', 'filteff', 'sum_of_weights', 'num_events']

# Los bordes se fijan antes de obtener resultados. El último bin recoge overflow
# y se rotula con >= cuando corresponda, para no perder eventos seleccionados.
BINS = {
    'electron_pt': np.array([25, 35, 45, 60, 80, 110, 150, 200, 300.]),
    'muon_pt': np.array([25, 35, 45, 60, 80, 110, 150, 200, 300.]),
    'n_jets': np.arange(1.5, 8.5, 1),
    'delta_phi': np.linspace(0, np.pi, 13),
    # |eta_e - eta_mu| no tiene unidades. La aceptación impone un máximo <5.
    # Fijamos intervalos de 0.5 antes de inspeccionar su distribución completa.
    'delta_eta': np.linspace(0, 5, 11),
    'm_emu': np.array([0, 25, 50, 75, 100, 125, 150, 200, 300, 450.]),
    'met': np.array([30, 45, 60, 80, 110, 150, 200, 300.]),
    'n_bjets': np.arange(-.5, 5.5, 1),
}


@dataclass(frozen=True)
class Config:
    """Parámetros visibles: una modificación define otra selección física."""
    lepton_pt: float = 25.0
    jet_pt: float = 25.0
    btag_quantile: int = 2
    luminosity_pb: float = 36000.0


def np1(array):
    """Convierte una columna escalar por evento en NumPy, sin valores ausentes."""
    return ak.to_numpy(array)


def pipeline_hashes():
    """Identifica también los controles importados: cambiarlos invalida la caché."""
    return {name: hashlib.sha256((BASE / name).read_bytes()).hexdigest()
            for name in ('analysis.py', 'audit_met_inputs.py', 'fetch_inputs.py')}


def select_events(a, cfg=Config()):
    """Forma máscaras por objeto y luego máscaras acumulativas por colisión.

    Un electrón/muón aceptado tiene pT>25 GeV, identificación y aislamiento
    tight, aceptación geométrica y matching de trigger. Se exige exactamente
    uno de cada especie entre los objetos aceptados; no se vetan otros objetos
    que no cumplan estos criterios. JVT se entrega como booleano en estos ROOT.
    SR1b implementa la propuesta; SR2b es la selección final del ejemplo oficial
    aplicada al skim 2J2LMET30. SS1b es un control de cargas iguales, no un fondo
    que se reste automáticamente.
    """
    eta = abs(a.lep_eta)
    quality = (a.lep_pt > cfg.lepton_pt) & a.lep_isTightID & a.lep_isTightIso
    ele_base = quality & (abs(a.lep_type) == 11) & ((eta < 1.37) | ((eta > 1.52) & (eta < 2.47)))
    mu_base = quality & (abs(a.lep_type) == 13) & (eta < 2.5)
    electron = ele_base & a.lep_isTrigMatched
    muon = mu_base & a.lep_isTrigMatched
    jet = (a.jet_pt > cfg.jet_pt) & (abs(a.jet_eta) < 2.5) & a.jet_jvt
    bjet = jet & (a.jet_btag_quantile >= cfg.btag_quantile)
    nj, nb = np1(ak.sum(jet, axis=1)), np1(ak.sum(bjet, axis=1))
    trigger = np1(a.trigE | a.trigM)
    pair = np1((ak.sum(electron, axis=1) == 1) & (ak.sum(muon, axis=1) == 1))
    q_e = np1(ak.sum(a.lep_charge * electron, axis=1))
    q_mu = np1(ak.sum(a.lep_charge * muon, axis=1))
    charge_product = q_e * q_mu
    pair_trigger = trigger & pair
    opposite = pair_trigger & (charge_product < 0)
    pretag = opposite & (nj >= 2)
    stages = {
        'input': np.ones(len(a), dtype=bool),
        'trigger': trigger,
        'one_e_one_mu_matched': pair_trigger,
        'opposite_sign': opposite,
        'two_jets_jvt': pretag,
        'one_bjet': pretag & (nb >= 1),
        'two_bjets': pretag & (nb >= 2),
    }
    regions = {'SR1b': stages['one_bjet'], 'SR2b': stages['two_bjets'],
               'PRETAG': pretag, 'ZEROb': pretag & (nb == 0),
               'SS1b': pair_trigger & (charge_product > 0) & (nj >= 2) & (nb >= 1)}
    return stages, regions, electron, muon, nj, nb


def event_weights(a, sample, cfg=Config()):
    """Calcula eventos esperados para la luminosidad, SIN ajustar a los datos.

    Datos: peso 1. MC: L[pb^-1] * sigma[pb] * filtro * k / sumW_original
    multiplicado por mcWeight y los factores de eficiencia disponibles.
    sumW_original es el total ANTES del skim, obtenido del CSV oficial y
    contrastado con cada bloque ROOT. Nunca se normaliza por la suma post-corte.
    Se conservan los pesos negativos de generadores NLO y los pesos cero.
    """
    if sample['group'] == 'data':
        return np.ones(len(a), dtype=np.float64)
    m = sample['metadata']
    denominator = float(m['sumOfWeights'])
    if denominator <= 0:
        raise ValueError('Normalización original no positiva')
    norm = cfg.luminosity_pb * float(m['crossSection_pb']) * float(m['genFiltEff']) * float(m['kFactor']) / denominator
    w = np1(a.mcWeight).astype(np.float64) * norm
    for field in SCALE_FACTORS:
        w *= np1(a[field]).astype(np.float64)
    if not np.all(np.isfinite(w)):
        raise ValueError('Peso final no finito: el análisis debe detenerse')
    return w


def validate_chunk(a, sample):
    """Contrato de entrada: longitudes, finitud, unidades y metadatos.

    Una rama ausente, un NaN o una normalización discordante abortan la corrida.
    No se reparan silenciosamente ni se eliminan eventos para mejorar el acuerdo.
    """
    for fields, count in ((LEPTON_FIELDS, 'lep_n'), (JET_FIELDS, 'jet_n')):
        for field in fields:
            if not ak.all(ak.num(a[field], axis=1) == a[count]):
                raise ValueError(f'Longitud incompatible en {field}')
            if not ak.all(np.isfinite(a[field])):
                raise ValueError(f'Valor no finito en {field}')
    for field in set(BRANCHES) - set(LEPTON_FIELDS + JET_FIELDS):
        if not np.all(np.isfinite(np1(a[field]))):
            raise ValueError(f'Valor no finito en {field}')
    if not ak.all((a.lep_pt >= 0) & (abs(a.lep_phi) <= np.pi + 1e-5)):
        raise ValueError('Momento o rango angular inválido')
    if not ak.all((a.jet_pt >= 0) & (a.jet_btag_quantile >= 1) & (a.jet_btag_quantile <= 5)):
        raise ValueError('pT de jet o categoría b-tag inválida')
    if not ak.all((a.lep_charge == 1) | (a.lep_charge == -1)):
        raise ValueError('Carga leptónica inesperada')
    # Los umbrales oficiales están en GeV: se comprueba el skim declarado.
    skim = (ak.sum(a.jet_pt >= 20 - 1e-5, axis=1) >= 2) & (ak.sum((a.lep_pt >= 7 - 1e-5) & a.lep_isTightID, axis=1) >= 2) & (a.met >= 30 - 1e-5)
    if not ak.all(skim):
        raise ValueError('Eventos incompatibles con el skim público 2J2LMET30')
    # Calcular en float64 y reconocer la resolución de las columnas float32.
    # El límite de 0.01 GeV se aplica a la separación mínima de sus intervalos
    # de redondeo, no a números puntuales con distinta resolución absoluta.
    met_residual = np.max(np.abs(np.hypot(np1(a.met_mpx).astype(float),
                                        np1(a.met_mpy).astype(float)) - np1(a.met).astype(float))) if len(a) else 0.
    met_gap, arithmetic = interval_gap(np1(a.met), np1(a.met_mpx), np1(a.met_mpy))
    if np.any(met_gap > .01 + arithmetic):
        raise ValueError('MET incompatible con sus componentes')
    if sample['group'] != 'data':
        mapping = {'xsec': 'crossSection_pb', 'filteff': 'genFiltEff', 'kfac': 'kFactor',
                   'sum_of_weights': 'sumOfWeights', 'num_events': 'nEvents'}
        for branch, column in mapping.items():
            reference = float(sample['metadata'][column])
            if not np.allclose(np1(a[branch]), reference, rtol=5e-5, atol=1e-8):
                raise ValueError(f'Metadato {branch} contradice CSV oficial para {sample["dsid"]}')
        if not np.all(np1(a.channelNumber) == sample['dsid']):
            raise ValueError('Identificador MC distinto al del manifiesto')
    return float(met_residual)


def observables(a, mask, electron, muon, nj, nb):
    """Deriva cantidades medibles en el laboratorio, sin reconstruir tops.

    Delta-phi usa atan2(sin(delta),cos(delta)) para respetar periodicidad.
    |Delta-eta| es abs(eta_e-eta_mu), no abs(eta_e)-abs(eta_mu).
    Ambos observables son del laboratorio; no son el ángulo usado para medir D.
    m(e,mu)^2=(Ee+Emu)^2-|pe+pmu|^2. No es la masa del par top-antitop.
    """
    values = {}
    for label, obj in (('electron', electron), ('muon', muon)):
        for field in ('pt', 'eta', 'phi', 'e'):
            values[f'{label}_{field}'] = np1(ak.firsts(a[f'lep_{field}'][obj][mask])).astype(float)
    delta = values['electron_phi'] - values['muon_phi']
    px = sum(values[f'{l}_pt'] * np.cos(values[f'{l}_phi']) for l in ('electron', 'muon'))
    py = sum(values[f'{l}_pt'] * np.sin(values[f'{l}_phi']) for l in ('electron', 'muon'))
    pz = sum(values[f'{l}_pt'] * np.sinh(values[f'{l}_eta']) for l in ('electron', 'muon'))
    e2 = (values['electron_e'] + values['muon_e'])**2
    p2 = px**2 + py**2 + pz**2
    mass, rounding = checked_mass(e2-p2, e2+p2)
    return {'electron_pt': values['electron_pt'], 'muon_pt': values['muon_pt'],
            'delta_phi': np.abs(np.arctan2(np.sin(delta), np.cos(delta))),
            'delta_eta': np.abs(values['electron_eta'] - values['muon_eta']),
            'm_emu': mass, 'met': np1(a.met)[mask],
            'n_jets': nj[mask], 'n_bjets': nb[mask]}, rounding


def checked_mass(m2, scale):
    """Distingue cancelación numérica de una masa claramente inconsistente.

    Los insumos son float32 aunque calculemos en float64. En pares casi
    colineales E²-|p|² resta números grandes y próximos. Se acepta solo una
    negatividad <=32*epsilon_float32*(E²+|p|²), fijada por precisión numérica,
    no por el ajuste a datos. Se asigna masa cero y se registra el caso; no se
    elimina el evento. Inconsistencias mayores abortan el análisis.
    """
    tolerance = 32 * np.finfo(np.float32).eps * np.maximum(scale, 1.)
    if np.any(m2 < -tolerance):
        raise ValueError('Masa invariante inconsistente más allá de la precisión float32')
    negative = m2 < 0
    diagnostic = {'count': int(negative.sum()),
                  'minimum_m2_gev2': float(min(0., np.min(m2))) if len(m2) else 0.,
                  'max_relative_negative': float(np.max(-m2[negative]/scale[negative])) if negative.any() else 0.}
    return np.sqrt(np.maximum(m2, 0)), diagnostic


def histogram(values, weights, edges):
    """Almacena sum(w), sum(w²) y conteo bruto, incluidos los desbordes.

    Varianza MC = sum(w²), no sum(w), ni sqrt del número generado sin pesos.
    Los flujos se conservan separados y se pliegan también en los bins extremos.
    """
    clipped = np.clip(values, edges[0], np.nextafter(edges[-1], edges[0]))
    s = np.histogram(clipped, edges, weights=weights)[0]
    v = np.histogram(clipped, edges, weights=weights**2)[0]
    n = np.histogram(clipped, edges)[0]
    return {'sumw': s.tolist(), 'sumw2': v.tolist(), 'raw': n.tolist(),
            'underflow_raw': int(np.sum(values < edges[0])),
            'overflow_raw': int(np.sum(values >= edges[-1]))}


def add_hist(target, update):
    """Suma bloques independientes sin volver a normalizar cada bloque."""
    for key in ('sumw', 'sumw2', 'raw'):
        target[key] = (np.array(target[key]) + np.array(update[key])).tolist()
    for key in ('underflow_raw', 'overflow_raw'):
        target[key] += update[key]


def normalized_shape(s, v):
    """Fracciones y covarianza: propaga también la fluctuación del total.

    p_i=s_i/S; J_ij=delta_ij/S-s_i/S²; Cov(p)=J diag(v) J^T.
    Esto evita ignorar las correlaciones inducidas al normalizar un histograma.
    """
    s, v = np.asarray(s, float), np.asarray(v, float)
    total = s.sum()
    if total <= 0:
        raise ValueError('Área no positiva en histograma normalizado')
    jac = np.eye(len(s)) / total - s[:, None] / total**2
    return s / total, (jac * v[None, :]) @ jac.T


def process_file(sample, cfg, chunk_size=100000):
    """Lee TODO el archivo por bloques y escribe un resultado auditable."""
    start = time.monotonic()
    path = BASE / 'data' / sample['key']
    digest = verify(path, sample)
    if digest is None or ('sha256' in sample and digest != sample['sha256']):
        raise ValueError('El archivo cambió después de la verificación de descarga')
    tree = uproot.open(path)['analysis']
    missing = sorted(set(BRANCHES) - set(tree.keys()))
    if missing:
        raise ValueError(f'Ramas requeridas ausentes: {missing}')
    result = {'file': sample['key'], 'dsid': sample['dsid'], 'group': sample['group'],
              'entries': tree.num_entries, 'branches': tree.typenames(), 'cutflow': {},
              'histograms': {}, 'regions': {}, 'processed': 0,
              'negative_generator_weights': 0, 'negative_final_weights': 0,
              'zero_final_weights': 0, 'max_met_residual_gev': 0., 'max_abs_weight': 0.,
              'metadata': sample['metadata'], 'mass_rounding': {},
              'extreme_met': {'input_over_13tev': 0, 'max_input_gev': 0., 'regions_over_13tev': {}}}
    result['input_sha256'] = digest
    for stage in ('input', 'trigger', 'one_e_one_mu_matched', 'opposite_sign', 'two_jets_jvt', 'one_bjet', 'two_bjets'):
        result['cutflow'][stage] = {'raw': 0, 'sumw': 0., 'sumw2': 0.}
    for region in ('SR1b', 'SR2b', 'PRETAG', 'ZEROb', 'SS1b'):
        result['regions'][region] = {'raw': 0, 'sumw': 0., 'sumw2': 0.}
        result['histograms'][region] = {k: histogram(np.array([]), np.array([]), b) for k, b in BINS.items()}
        result['mass_rounding'][region] = {'count': 0, 'minimum_m2_gev2': 0., 'max_relative_negative': 0.}
        result['extreme_met']['regions_over_13tev'][region] = 0
    ids = []
    for a in tree.iterate(expressions=BRANCHES, step_size=chunk_size, library='ak'):
        residual = validate_chunk(a, sample)
        w = event_weights(a, sample, cfg)
        stages, regions, electron, muon, nj, nb = select_events(a, cfg)
        # Diagnóstico, NO corte añadido: una cantidad reconstruida patológica
        # no se interpreta como energía física disponible en la colisión.
        extreme_met = np1(a.met) > 13000.
        result['extreme_met']['input_over_13tev'] += int(extreme_met.sum())
        result['extreme_met']['max_input_gev'] = max(result['extreme_met']['max_input_gev'], float(np.max(np1(a.met))))
        # Con pesos negativos el rendimiento ponderado no tiene que decrecer;
        # sí deben ser anidados los conjuntos y sus conteos brutos.
        previous = np.ones(len(a), dtype=bool)
        for label, mask in stages.items():
            if np.any(mask & ~previous):
                raise AssertionError('Cutflow no acumulativo')
            previous = mask
            stat = result['cutflow'].setdefault(label, {'raw': 0, 'sumw': 0., 'sumw2': 0.})
            stat['raw'] += int(mask.sum())
            stat['sumw'] += float(w[mask].sum())
            stat['sumw2'] += float(np.sum(w[mask]**2))
        for region, mask in regions.items():
            result['extreme_met']['regions_over_13tev'][region] += int(np.sum(extreme_met & mask))
            stat = result['regions'].setdefault(region, {'raw': 0, 'sumw': 0., 'sumw2': 0.})
            stat['raw'] += int(mask.sum())
            stat['sumw'] += float(w[mask].sum())
            stat['sumw2'] += float(np.sum(w[mask]**2))
            dest = result['histograms'].setdefault(region, {})
            obs, rounding = observables(a, mask, electron, muon, nj, nb)
            audit = result['mass_rounding'][region]
            audit['count'] += rounding['count']
            audit['minimum_m2_gev2'] = min(audit['minimum_m2_gev2'], rounding['minimum_m2_gev2'])
            audit['max_relative_negative'] = max(audit['max_relative_negative'], rounding['max_relative_negative'])
            for name, edges in BINS.items():
                h = histogram(obs[name], w[mask], edges)
                if name in dest:
                    add_hist(dest[name], h)
                else:
                    dest[name] = h
        result['processed'] += len(a)
        result['negative_generator_weights'] += int(np.sum(np1(a.mcWeight) < 0))
        result['negative_final_weights'] += int(np.sum(w < 0))
        result['zero_final_weights'] += int(np.sum(w == 0))
        result['max_abs_weight'] = max(result['max_abs_weight'], float(np.max(np.abs(w))))
        result['max_met_residual_gev'] = max(result['max_met_residual_gev'], residual)
        if sample['group'] == 'data':
            pair = np.empty(len(a), dtype=[('run', '<u4'), ('event', '<u8')])
            pair['run'], pair['event'] = np1(a.runNumber), np1(a.eventNumber)
            ids.append(pair)
    if result['processed'] != tree.num_entries:
        raise AssertionError('No se leyó el archivo completo')
    for region, hists in result['histograms'].items():
        for h in hists.values():
            assert sum(h['raw']) == result['regions'][region]['raw']
            assert np.isclose(sum(h['sumw']), result['regions'][region]['sumw'])
            assert np.isclose(sum(h['sumw2']), result['regions'][region]['sumw2'])
    if ids:
        all_ids = np.concatenate(ids)
        if len(np.unique(all_ids)) != len(all_ids):
            raise ValueError('Identificadores repetidos en archivo de colisiones')
        np.save(RESULTS / 'event_ids' / (path.stem + '.npy'), all_ids)
    result['seconds'] = round(time.monotonic() - start, 2)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--only-dsid', type=int)
    parser.add_argument('--chunk-size', type=int, default=100000)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--available-only', action='store_true', help='Procesa archivos ya descargados; marca salida como parcial')
    args = parser.parse_args()
    manifest_path = BASE / 'provenance' / 'verified_inputs.json'
    if not manifest_path.exists():
        if not args.available_only:
            raise ValueError('La descarga completa aún no está verificada')
        manifest_path = BASE / 'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    (RESULTS / 'samples').mkdir(parents=True, exist_ok=True)
    (RESULTS / 'event_ids').mkdir(exist_ok=True)
    cfg = Config(luminosity_pb=manifest['luminosity_pb'])
    code_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    hashes = pipeline_hashes()
    run = {'config': asdict(cfg), 'code_sha256': code_sha, 'pipeline_hashes': hashes,
           'manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
           'versions': {p: importlib.metadata.version(p) for p in ('numpy', 'uproot', 'awkward', 'matplotlib', 'scipy')},
           'bins': {k: v.tolist() for k, v in BINS.items()}, 'files': []}
    for sample in manifest['files']:
        if args.available_only and not (BASE / 'data' / sample['key']).exists():
            continue
        if args.only_dsid is not None and sample['dsid'] != args.only_dsid:
            continue
        dest = RESULTS / 'samples' / (sample['key'] + '.json')
        if args.resume and dest.exists():
            saved = json.loads(dest.read_text())
            if saved.get('code_sha256') != code_sha or saved.get('config') != asdict(cfg) or saved.get('pipeline_hashes') != hashes:
                raise ValueError('Caché de una configuración/código diferente; volver a ejecutar')
            if saved.get('metadata') != sample['metadata']:
                raise ValueError('Los metadatos cambiaron respecto a la corrida guardada')
            if sample.get('sha256') and saved.get('input_sha256') != sample['sha256']:
                raise ValueError('El archivo verificado no es el usado para el resultado guardado')
        else:
            saved = dict(process_file(sample, cfg, args.chunk_size), code_sha256=code_sha, config=asdict(cfg), pipeline_hashes=hashes)
            dest.write_text(json.dumps(saved, indent=2), encoding='utf-8')
        run['files'].append(str(dest.relative_to(BASE)))
        print(sample['dsid'] or sample['key'], saved['processed'], saved['regions']['SR1b'], flush=True)
    run['complete'] = len(run['files']) == len(manifest['files'])
    if run['complete']:
        ids = np.concatenate([np.load(p) for p in sorted((RESULTS / 'event_ids').glob('*.npy'))])
        if len(np.unique(ids)) != len(ids):
            raise ValueError('Eventos de datos duplicados ENTRE archivos')
        if len(ids) != 6242521:
            raise ValueError('Conteo global de datos distinto al catálogo')
        run['global_checks'] = {'unique_data_events': len(ids), 'expected_data_events': 6242521,
                                'all_checks_passed': True}
        print('CONTROL GLOBAL: 6242521 eventos reales únicos', flush=True)
    destination = RESULTS / ('run.json' if run['complete'] else 'run_partial.json')
    temporary = destination.with_suffix('.json.part')
    temporary.write_text(json.dumps(run, indent=2), encoding='utf-8')
    temporary.replace(destination)


if __name__ == '__main__':
    main()
