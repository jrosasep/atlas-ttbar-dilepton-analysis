"""Ejemplos pequeños para comprobar el análisis y sus controles.

Una prueba compara lo calculado con una respuesta conocida independientemente.
Las pruebas normales no descargan datos. Las de integración abren ROOT locales.
Estos eventos inventados comprueban código; no son evidencia física de ATLAS.
"""
import copy                  # Cambiar un ejemplo sin alterar el original.
import hashlib               # Comprobar las huellas de procedencia.
import json                  # Leer manifiestos y registros.
import zlib                  # Calcular el checksum Adler-32 de una entrada.

import awkward as ak         # Representar listas de objetos por evento.
import numpy as np           # Cálculos y comparaciones numéricas.
import pytest                # Organizar pruebas y comprobar errores esperados.
import uproot                # Abrir archivos ROOT en las pruebas de integración.

import fetch_inputs
import plots_and_summary as plots
from analysis import (
    BASE, BRANCHES, SCALE_FACTORS, Config, checked_mass, event_weights,
    histogram, normalized_shape, observables, pipeline_hashes, select_events,
    validate_chunk,
)
from audit_met_inputs import interval_gap


def example():
    """Un evento inventado con e+, mu-, dos jets y dos etiquetas b.

    Cada lista contiene las propiedades de los objetos del mismo evento.
    Las posiciones deben corresponder: lep_pt[0] y lep_eta[0] describen
    el mismo leptón. Los campos de jets siguen su propia lista de objetos.
    """
    return {
        'lep_pt': [50., 40.],
        'lep_eta': [.3, -.6],
        'lep_phi': [3.1, -3.1],
        'lep_type': [11, 13],
        'lep_charge': [1, -1],
        'lep_isTightID': [True, True],
        'lep_isTightIso': [True, True],
        'lep_isTrigMatched': [True, True],
        'jet_pt': [60., 30.],
        'jet_eta': [.1, -.9],
        'jet_jvt': [True, True],
        'jet_btag_quantile': [2, 5],
        'trigE': True,
        'trigM': False,
    }


def test_selection_hand_cases():
    """Cambiar una condición por evento y comprobar qué selección supera."""
    rows = []
    for _ in range(7):
        rows.append(example())
    rows[1]['trigE'] = False
    rows[2]['lep_pt'][0] = 25.0  # Se exige pT > 25, no pT >= 25.
    rows[3]['lep_charge'][1] = 1 # Cargas iguales: región de control SS.
    rows[4]['jet_jvt'][1] = False
    rows[5]['jet_btag_quantile'][1] = 1
    rows[6]['lep_isTrigMatched'][0] = False

    eventos = ak.Array(rows)
    stages, regions, electron, muon, nj, nb = select_events(eventos)
    indices_1b = np.flatnonzero(regions['SR1b']).tolist()
    indices_2b = np.flatnonzero(regions['SR2b']).tolist()
    indices_ss = np.flatnonzero(regions['SS1b']).tolist()
    assert indices_1b == [0, 5]
    assert indices_2b == [0]
    assert indices_ss == [3]


def test_weights_have_physical_normalization_and_keep_negative_sign():
    """w = L*sigma*eficiencia*k/suma_original * peso_MC * factores.

    En este ejemplo: 1000*10*0.5*2/100 = 100. Los pesos finales son
    100*2 y 100*(-1). Un peso negativo no se reemplaza por su absoluto.
    """
    columnas = {'mcWeight': [2., -1.]}
    for nombre in SCALE_FACTORS:
        columnas[nombre] = [1., 1.]
    eventos = ak.Array(columnas)
    muestra = {
        'group': 'ttbar',
        'metadata': {
            'crossSection_pb': '10', 'genFiltEff': '.5',
            'kFactor': '2', 'sumOfWeights': '100',
        },
    }
    configuracion = Config(luminosity_pb=1000)
    pesos_mc = event_weights(eventos, muestra, configuracion)
    pesos_datos = event_weights(eventos, {'group': 'data'})
    assert np.allclose(pesos_mc, [200, -100])
    assert np.all(pesos_datos == 1)


def test_histogram_flow_and_signed_variance():
    """Los valores fuera del rango se guardan en los bins de los extremos.

    Primer bin: pesos 1+2; varianza 1²+2². Último: -3+4; varianza
    (-3)²+4². La varianza nunca se calcula como el cuadrado de la suma.
    """
    valores = np.array([-1, .2, 1.5, 2.])
    pesos = np.array([1., 2., -3., 4.])
    bordes = np.array([0., 1., 2.])
    h = histogram(valores, pesos, bordes)
    assert h['sumw'] == [3., 1.]
    assert h['sumw2'] == [5., 25.]
    assert h['raw'] == [2, 2]
    assert h['underflow_raw'] == h['overflow_raw'] == 1


def test_normalized_covariance_matches_multinomial():
    """Para conteos sin pesos, la covarianza tiene una expresión conocida.

    Normalizar introduce correlaciones: si una fracción aumenta, otras
    deben disminuir porque la suma permanece en uno. Sólo hay K-1 grados
    de libertad para K bins. Comprobamos esas dos propiedades.
    """
    counts = np.array([20., 30., 50.])
    p, cov = normalized_shape(counts, counts)
    expected = (np.diag(p) - np.outer(p, p)) / counts.sum()
    assert np.allclose(cov, expected)
    assert np.allclose(cov.sum(axis=0), 0)
    assert np.linalg.matrix_rank(cov) == 2


def scalar_reference(row):
    """Comprobar SR2b con bucles, independientemente de las máscaras Awkward.

    Tener dos implementaciones con la misma respuesta ayuda a detectar
    errores de índices. No demuestra que el modelo físico sea completo.
    """
    if not (row['trigE'] or row['trigM']):
        return False

    leptons = []
    for i, pt in enumerate(row['lep_pt']):
        pasa_pt = pt > 25
        pasa_identificacion = row['lep_isTightID'][i]
        pasa_aislamiento = row['lep_isTightIso'][i]
        pasa_trigger = row['lep_isTrigMatched'][i]
        if not (pasa_pt and pasa_identificacion and pasa_aislamiento and pasa_trigger):
            continue
        eta = abs(row['lep_eta'][i])
        kind = abs(row['lep_type'][i])
        electron_aceptado = kind == 11 and eta < 2.47 and (eta < 1.37 or eta > 1.52)
        muon_aceptado = kind == 13 and eta < 2.5
        if electron_aceptado or muon_aceptado:
            leptons.append(i)

    if len(leptons) != 2:
        return False
    i, j = leptons
    tipos = {abs(row['lep_type'][i]), abs(row['lep_type'][j])}
    if tipos != {11, 13}:
        return False
    if row['lep_charge'][i] * row['lep_charge'][j] >= 0:
        return False

    nb = 0
    for k, pt in enumerate(row['jet_pt']):
        pasa_pt = pt > 25
        pasa_eta = abs(row['jet_eta'][k]) < 2.5
        pasa_jvt = row['jet_jvt'][k]
        tiene_etiqueta_b = row['jet_btag_quantile'][k] >= 2
        if pasa_pt and pasa_eta and pasa_jvt and tiene_etiqueta_b:
            nb += 1
    return nb >= 2


def real_data_path():
    """Localizar el ROOT de colisiones usado por las pruebas de integración."""
    files = list((BASE / 'data').glob('*data15_periodD*.root'))
    assert files, 'Se requiere el ROOT real verificado en analysis/data/'
    return files[0]


@pytest.mark.integration
def test_real_file_vectorization_matches_scalar_reference():
    """La selección por columnas debe coincidir con la selección evento a evento."""
    with uproot.open(real_data_path()) as archivo:
        eventos = archivo['analysis'].arrays(BRANCHES, entry_stop=10000, library='ak')
    expected = []
    for evento in ak.to_list(eventos):
        expected.append(scalar_reference(evento))
    actual = select_events(eventos)[1]['SR2b']
    assert np.array_equal(np.array(expected), actual)

    selecciones_por_bloque = []
    for inicio in range(0, len(eventos), 791):
        bloque = eventos[inicio:inicio + 791]
        selecciones_por_bloque.append(select_events(bloque)[1]['SR2b'])
    chunked = np.concatenate(selecciones_por_bloque)
    assert np.array_equal(actual, chunked)


@pytest.mark.integration
def test_corrupted_met_is_rejected():
    """Modificar sólo MET rompe su consistencia con las componentes x e y."""
    with uproot.open(real_data_path()) as archivo:
        eventos = archivo['analysis'].arrays(BRANCHES, entry_stop=20, library='ak')
    broken = ak.with_field(eventos, eventos.met + 50, 'met')
    with pytest.raises(ValueError, match='MET incompatible'):
        validate_chunk(broken, {'group': 'data'})


@pytest.mark.integration
def test_data_event_count_small_period_is_known():
    """La muestra de referencia contiene 10085 entradas, antes de los cortes."""
    with uproot.open(real_data_path()) as archivo:
        assert archivo['analysis'].num_entries == 10085


@pytest.mark.integration
def test_wrong_original_weight_sum_is_rejected():
    """No aceptar un metadato de normalización incompatible con el ROOT."""
    manifest = json.loads((BASE / 'manifest.json').read_text(encoding='utf-8'))
    sample = next(s for s in manifest['files'] if s['dsid'] == 410645)
    path = BASE / 'data' / sample['key']
    with uproot.open(path) as archivo:
        eventos = archivo['analysis'].arrays(BRANCHES, entry_stop=20, library='ak')
    validate_chunk(eventos, sample)
    broken = copy.deepcopy(sample)
    suma_original = float(broken['metadata']['sumOfWeights'])
    broken['metadata']['sumOfWeights'] = str(2 * suma_original)
    with pytest.raises(ValueError, match='sum_of_weights'):
        validate_chunk(eventos, broken)


def test_checksum_rejects_same_size_modified_input(tmp_path):
    """El tamaño correcto no basta: cambiar un byte debe detectarse."""
    path = tmp_path / 'test.root'
    content = b'input verificable de prueba'
    path.write_bytes(content)
    info = {'size': len(content), 'checksum': f'adler32:{zlib.adler32(content):08x}'}
    assert fetch_inputs.verify(path, info) == hashlib.sha256(content).hexdigest()
    path.write_bytes(b'X' + content[1:])
    with pytest.raises(ValueError, match='Checksum'):
        fetch_inputs.verify(path, info)


def test_missing_or_incomplete_input_returns_none(tmp_path):
    """Distinguir un archivo ausente/incompleto de uno íntegro."""
    path = tmp_path / 'test.root'
    info = {'size': 20, 'checksum': 'adler32:00000000'}
    assert fetch_inputs.verify(path, info) is None
    path.write_bytes(b'incompleto')
    assert fetch_inputs.verify(path, info) is None


def test_mass_rounding_is_logged_but_large_inconsistency_is_rejected():
    """Una masa al cuadrado apenas negativa puede resultar del redondeo.

    Sólo se permite la corrección dentro de la tolerancia documentada;
    una inconsistencia grande debe detener el análisis.
    """
    mass, audit = checked_mass(np.array([-.01514, 100.]), np.array([1712482., 1000.]))
    assert mass.tolist() == [0., 10.]
    assert audit['count'] == 1
    with pytest.raises(ValueError, match='precisión'):
        checked_mass(np.array([-20.]), np.array([1712482.]))


def test_observables_match_independent_massless_formula_and_periodic_angle():
    """m² = 2*pTe*pTmu*(cosh(delta_eta)-cos(delta_phi)) si las masas son cero.

    Esta expresión permite revisar el cálculo por cuadrivectores con otro
    cálculo. Los ángulos 3.1 y -3.1 están cerca, no separados por 6.2 rad.
    """
    row = example()
    momentos = np.array(row['lep_pt'])
    pseudorrapideces = np.array(row['lep_eta'])
    row['lep_e'] = (momentos * np.cosh(pseudorrapideces)).tolist()
    row['met'] = 80.
    eventos = ak.Array([row])
    _, regions, electron, muon, nj, nb = select_events(eventos)
    obs, rounding = observables(eventos, regions['SR1b'], electron, muon, nj, nb)
    delta = 2 * np.pi - 6.2
    mass_squared = 2 * 50 * 40 * (np.cosh(.9) - np.cos(delta))
    assert np.allclose(obs['delta_phi'], [delta], atol=1e-12)
    assert np.allclose(obs['delta_eta'], [.9], atol=1e-12)
    assert np.allclose(obs['m_emu'] ** 2, [mass_squared], atol=1e-9)
    assert obs['n_jets'].tolist() == [2]
    assert obs['n_bjets'].tolist() == [2]
    assert rounding['count'] == 0


def test_poisson_interval_keeps_uncertainty_for_zero_observations():
    """Cero sucesos observados no significa incertidumbre cero."""
    errors = plots.poisson_errors(np.array([0, 100]))
    assert errors[0, 0] == 0
    assert np.isclose(errors[1, 0], 1.841021645, atol=1e-6)
    assert 9 < errors[0, 1] < 11
    assert 10 < errors[1, 1] < 12


def test_large_float32_met_rounding_is_not_a_physical_failure():
    """Calcular en float64 no debe confundir redondeo previo con incoherencia."""
    x = np.array([1e6], dtype=np.float32)
    radius = np.hypot(x, x)
    recalculado = np.hypot(x.astype(float), x.astype(float))
    assert abs(recalculado - radius)[0] > .01
    gap, arithmetic = interval_gap(radius, x, x)
    assert gap[0] <= arithmetic[0]


@pytest.mark.integration
def test_extreme_met_input_is_checked_in_double_precision_and_fails_selection():
    """Registro real anómalo: no relajar tolerancias ni borrar el evento."""
    manifest = json.loads((BASE / 'manifest.json').read_text(encoding='utf-8'))
    sample = next(s for s in manifest['files'] if s['dsid'] == 700325)
    path = BASE / 'data' / sample['key']
    with uproot.open(path) as archivo:
        eventos = archivo['analysis'].arrays(
            BRANCHES, entry_start=6869582, entry_stop=6869583, library='ak'
        )
    assert .008 < validate_chunk(eventos, sample) < .01
    for mask in select_events(eventos)[1].values():
        assert not mask.any()


@pytest.mark.parametrize('etas, expected', [
    ([.4, .4], 0.), ([.4, -.4], .8),
    ([-2.4, 2.49], 4.89), ([2.4, -2.49], 4.89),
])
def test_abs_delta_eta_geometry_and_lepton_order(etas, expected):
    """Reflejar el haz o intercambiar los leptones no cambia |delta_eta|."""
    for order in ([0, 1], [1, 0]):
        row = example()
        row['lep_eta'] = etas
        row['lep_e'] = (np.array(row['lep_pt']) * np.cosh(etas)).tolist()
        for name, values in list(row.items()):
            if name.startswith('lep_'):
                row[name] = [values[i] for i in order]
        row['met'] = 80.
        eventos = ak.Array([row])
        _, regions, electron, muon, nj, nb = select_events(eventos)
        obs, _ = observables(eventos, regions['SR2b'], electron, muon, nj, nb)
        assert obs['delta_eta'].shape == (1,)
        assert obs['delta_eta'][0] == pytest.approx(expected)


def test_pipeline_fingerprint_covers_imported_numerical_checks():
    """La procedencia debe incluir los módulos que controlan las entradas."""
    hashes = pipeline_hashes()
    assert set(hashes) == {'analysis.py', 'audit_met_inputs.py', 'fetch_inputs.py'}
    for value in hashes.values():
        assert len(value) == 64


class FakeResponse:
    """Respuesta de red inventada: permite probar sin acceder a Internet."""
    def __init__(self, content, status=206, headers=None):
        self.content = content
        self.status = status
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self):
        return self.content


def test_ranged_download_keeps_the_requested_bytes(tmp_path, monkeypatch):
    """Una respuesta parcial válida debe escribirse sin cambiar sus bytes."""
    content = b'contenido de prueba para un rango'
    headers = {'Content-Range': f'bytes 0-{len(content) - 1}/{len(content)}'}

    def simulated_response(request, timeout):
        assert request.get_header('Range') == f'bytes=0-{len(content) - 1}'
        return FakeResponse(content, headers=headers)

    monkeypatch.setattr(fetch_inputs.urllib.request, 'urlopen', simulated_response)
    path = tmp_path / 'download.part'
    fetch_inputs.ranged_download('https://example.invalid/test.root', path, len(content))
    assert path.read_bytes() == content


@pytest.mark.parametrize('status, content_range, payload', [
    (200, 'bytes 0-2/3', b'abc'),  # El servidor ignoró Range.
    (206, 'bytes 1-3/3', b'abc'),  # Devolvió otro rango.
    (206, 'bytes 0-2/3', b'ab'),   # Faltan bytes aunque el encabezado coincida.
])
def test_invalid_range_response_is_rejected(tmp_path, monkeypatch, status, content_range, payload):
    """No aceptar una descarga incompleta o de otra posición del archivo."""
    def simulated_response(request, timeout):
        return FakeResponse(payload, status, {'Content-Range': content_range})

    monkeypatch.setattr(fetch_inputs.urllib.request, 'urlopen', simulated_response)
    monkeypatch.setattr(fetch_inputs.time, 'sleep', lambda seconds: None)
    with pytest.raises(ValueError):
        fetch_inputs.ranged_download('https://example.invalid/test.root', tmp_path / 'bad.part', 3)


def test_references_use_commit_for_content_and_tree_for_listing(tmp_path, monkeypatch):
    """Un commit y un árbol de Git no son la misma referencia."""
    revision = fetch_inputs.OFFICIAL_REVISION
    tree_sha = '1' * 40
    calls = []

    def simulated_fetch(url):
        calls.append(url)
        if '/commits/' in url:
            return json.dumps({'sha': revision, 'commit': {'tree': {'sha': tree_sha}}}).encode()
        if '/git/trees/' in url:
            tree = {'tree': [{'type': 'blob', 'path': 'LICENSE'}]}
            return json.dumps(tree).encode()
        return b'contenido de referencia'

    monkeypatch.setattr(fetch_inputs, 'BASE', tmp_path)
    monkeypatch.setattr(fetch_inputs, 'fetch', simulated_fetch)
    fetch_inputs.references()
    assert any('/git/trees/' + tree_sha in url for url in calls)
    raw_urls = [url for url in calls if 'raw.githubusercontent.com' in url]
    assert len(raw_urls) == 1
    assert '/' + revision + '/LICENSE' in raw_urls[0]


def source_record_example(tmp_path, monkeypatch):
    """Preparar una corrida ficticia con sus tres fuentes y su manifiesto."""
    hashes = {}
    for name in ('analysis.py', 'audit_met_inputs.py', 'fetch_inputs.py'):
        hashes[name] = hashlib.sha256(name.encode()).hexdigest()
    code_hash = hashes['analysis.py']
    directory = tmp_path / 'provenance' / 'source_records' / code_hash
    directory.mkdir(parents=True)
    for name in hashes:
        (directory / (name + '.txt')).write_bytes(name.encode())
    manifest_bytes = b'{"files": []}'
    (directory / 'input_manifest.json.txt').write_bytes(manifest_bytes)
    (tmp_path / 'analysis.py').write_bytes(b'codigo actual distinto')
    run = {
        'code_sha256': code_hash, 'pipeline_hashes': hashes,
        'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
        'complete': True, 'global_checks': {'all_checks_passed': True},
    }
    monkeypatch.setattr(plots, 'BASE', tmp_path)
    monkeypatch.setattr(plots, 'pipeline_hashes', lambda: {'actual': 'otra huella'})
    return run, directory


def test_stored_run_requires_explicit_mode_and_original_sources(tmp_path, monkeypatch):
    """Redibujar una corrida antigua requiere reconocerlo y verificar su fuente."""
    run, directory = source_record_example(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match='difiere'):
        plots.validate_run(run)
    assert plots.validate_run(run, stored_run=True) == {'files': []}
    (directory / 'audit_met_inputs.py.txt').write_bytes(b'modificado')
    with pytest.raises(ValueError, match='modificada'):
        plots.validate_run(run, stored_run=True)


def test_duplicate_saved_sample_is_rejected():
    """La misma muestra no puede contarse dos veces."""
    run = {'files': ['sample.json', 'sample.json']}
    with pytest.raises(ValueError, match='dos veces'):
        plots.check_saved_results(run, {'files': []}, {})
