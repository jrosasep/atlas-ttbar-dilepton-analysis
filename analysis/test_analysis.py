"""Pruebas con resultados calculables a mano e implementación escalar independiente."""
import copy
import json
from pathlib import Path

import awkward as ak
import numpy as np
import pytest
import uproot

from analysis import (BASE, Config, event_weights, histogram, normalized_shape,
                      select_events, validate_chunk, BRANCHES, checked_mass, observables)


def example():
    """Colisión de prueba con e+, mu-, dos jets y dos etiquetas b."""
    return dict(lep_pt=[50., 40.], lep_eta=[.3, -.6], lep_phi=[3.1, -3.1],
                lep_type=[11, 13], lep_charge=[1, -1], lep_isTightID=[True, True],
                lep_isTightIso=[True, True], lep_isTrigMatched=[True, True],
                jet_pt=[60., 30.], jet_eta=[.1, -.9], jet_jvt=[True, True],
                jet_btag_quantile=[2, 5], trigE=True, trigM=False)


def test_selection_hand_cases():
    rows = [example() for _ in range(7)]
    rows[1]['trigE'] = False
    rows[2]['lep_pt'][0] = 25.0  # La desigualdad oficial es estricta.
    rows[3]['lep_charge'][1] = 1  # Debe ir al control SS.
    rows[4]['jet_jvt'][1] = False
    rows[5]['jet_btag_quantile'][1] = 1
    rows[6]['lep_isTrigMatched'][0] = False
    stages, regions, *_ = select_events(ak.Array(rows))
    assert np.flatnonzero(regions['SR1b']).tolist() == [0, 5]
    assert np.flatnonzero(regions['SR2b']).tolist() == [0]
    assert np.flatnonzero(regions['SS1b']).tolist() == [3]


def test_weights_have_physical_normalization_and_keep_negative_sign():
    a = ak.Array({'mcWeight': [2., -1.], **{f: [1., 1.] for f in (
        'ScaleFactor_PILEUP', 'ScaleFactor_ELE', 'ScaleFactor_MUON',
        'ScaleFactor_LepTRIGGER', 'ScaleFactor_FTAG', 'ScaleFactor_JVT')}})
    m = {'group': 'ttbar', 'metadata': {'crossSection_pb': '10', 'genFiltEff': '.5',
                                      'kFactor': '2', 'sumOfWeights': '100'}}
    assert np.allclose(event_weights(a, m, Config(luminosity_pb=1000)), [200, -100])
    assert np.all(event_weights(a, {'group': 'data'}) == 1)


def test_histogram_flow_and_signed_variance():
    h = histogram(np.array([-1, .2, 1.5, 2.]), np.array([1., 2., -3., 4.]), np.array([0., 1., 2.]))
    assert h['sumw'] == [3., 1.]
    assert h['sumw2'] == [5., 25.]
    assert h['raw'] == [2, 2]
    assert h['underflow_raw'] == h['overflow_raw'] == 1


def test_normalized_covariance_matches_multinomial():
    counts = np.array([20., 30., 50.])
    p, cov = normalized_shape(counts, counts)
    expected = (np.diag(p) - np.outer(p, p)) / counts.sum()
    assert np.allclose(cov, expected)
    assert np.allclose(cov.sum(axis=0), 0)
    assert np.linalg.matrix_rank(cov) == 2


def scalar_reference(row):
    """Traducción literal en bucles del criterio oficial, sin máscaras Awkward."""
    if not (row['trigE'] or row['trigM']):
        return False
    leptons = []
    for i, pt in enumerate(row['lep_pt']):
        if not (pt > 25 and row['lep_isTightID'][i] and row['lep_isTightIso'][i] and row['lep_isTrigMatched'][i]):
            continue
        eta, kind = abs(row['lep_eta'][i]), abs(row['lep_type'][i])
        if (kind == 11 and eta < 2.47 and (eta < 1.37 or eta > 1.52)) or (kind == 13 and eta < 2.5):
            leptons.append(i)
    if len(leptons) != 2:
        return False
    i, j = leptons
    if row['lep_type'][i] == row['lep_type'][j] or row['lep_charge'][i] * row['lep_charge'][j] >= 0:
        return False
    nb = sum(pt > 25 and abs(row['jet_eta'][k]) < 2.5 and row['jet_jvt'][k] and row['jet_btag_quantile'][k] >= 2
             for k, pt in enumerate(row['jet_pt']))
    return nb >= 2


@pytest.mark.integration
def test_real_file_vectorization_matches_scalar_reference():
    files = list((BASE / 'data').glob('*data15_periodD*.root'))
    assert files, 'Se requiere el archivo real verificado para esta prueba de integración'
    a = uproot.open(files[0])['analysis'].arrays(BRANCHES, entry_stop=10000, library='ak')
    expected = np.array([scalar_reference(r) for r in ak.to_list(a)])
    actual = select_events(a)[1]['SR2b']
    assert np.array_equal(expected, actual)
    # También debe dar lo mismo analizar todo de una vez o en bloques.
    chunked = np.concatenate([select_events(a[i:i+791])[1]['SR2b'] for i in range(0, len(a), 791)])
    assert np.array_equal(actual, chunked)


@pytest.mark.integration
def test_corrupted_met_is_rejected():
    path = next((BASE / 'data').glob('*data15_periodD*.root'))
    a = uproot.open(path)['analysis'].arrays(BRANCHES, entry_stop=20, library='ak')
    broken = ak.with_field(a, a.met + 50, 'met')
    with pytest.raises(ValueError, match='MET incompatible'):
        validate_chunk(broken, {'group': 'data'})


@pytest.mark.integration
def test_data_event_count_small_period_is_known():
    path = next((BASE / 'data').glob('*data15_periodD*.root'))
    assert uproot.open(path)['analysis'].num_entries == 10085


@pytest.mark.integration
def test_wrong_original_weight_sum_is_rejected():
    manifest = json.loads((BASE / 'manifest.json').read_text(encoding='utf-8'))
    sample = next(s for s in manifest['files'] if s['dsid'] == 410645)
    path = BASE / 'data' / sample['key']
    a = uproot.open(path)['analysis'].arrays(BRANCHES, entry_stop=20, library='ak')
    validate_chunk(a, sample)
    broken = copy.deepcopy(sample)
    broken['metadata']['sumOfWeights'] = str(2 * float(broken['metadata']['sumOfWeights']))
    with pytest.raises(ValueError, match='sum_of_weights'):
        validate_chunk(a, broken)


def test_checksum_rejects_same_size_modified_input(tmp_path):
    from fetch_inputs import verify
    import zlib
    path = tmp_path / 'test.root'
    content = b'input verificable de prueba'
    path.write_bytes(content)
    info = {'size': len(content), 'checksum': f'adler32:{zlib.adler32(content):08x}'}
    assert verify(path, info)
    path.write_bytes(b'X' + content[1:])
    with pytest.raises(ValueError, match='Checksum'):
        verify(path, info)


def test_mass_rounding_is_logged_but_large_inconsistency_is_rejected():
    mass, audit = checked_mass(np.array([-.01514, 100.]), np.array([1712482., 1000.]))
    assert mass.tolist() == [0., 10.]
    assert audit['count'] == 1
    with pytest.raises(ValueError, match='precisión'):
        checked_mass(np.array([-20.]), np.array([1712482.]))


def test_observables_match_independent_massless_formula_and_periodic_angle():
    """Una pareja sin masa tiene una fórmula cerrada que evita sumar cuadrivectores."""
    row = example()
    row['lep_e'] = (np.array(row['lep_pt']) * np.cosh(row['lep_eta'])).tolist()
    row['met'] = 80.
    a = ak.Array([row])
    _, regions, electron, muon, nj, nb = select_events(a)
    obs, rounding = observables(a, regions['SR1b'], electron, muon, nj, nb)
    delta = 2*np.pi - 6.2
    mass_squared = 2*50*40*(np.cosh(.9)-np.cos(delta))
    assert np.allclose(obs['delta_phi'], [delta], atol=1e-12)
    assert np.allclose(obs['delta_eta'], [.9], atol=1e-12)
    assert np.allclose(obs['m_emu']**2, [mass_squared], atol=1e-9)
    assert obs['n_jets'].tolist() == [2]
    assert obs['n_bjets'].tolist() == [2]
    assert rounding['count'] == 0


def test_poisson_interval_keeps_uncertainty_for_zero_observations():
    """Cero sucesos observados no significa incertidumbre cero."""
    from plots_and_summary import poisson_errors
    errors = poisson_errors(np.array([0, 100]))
    assert errors[0, 0] == 0
    assert np.isclose(errors[1, 0], 1.841021645, atol=1e-6)
    assert 9 < errors[0, 1] < 11
    assert 10 < errors[1, 1] < 12


@pytest.mark.integration
def test_extreme_met_input_is_checked_in_double_precision_and_fails_selection():
    """Regresión del registro anómalo: no alterar tolerancia ni borrar el evento."""
    manifest = json.loads((BASE / 'manifest.json').read_text(encoding='utf-8'))
    sample = next(s for s in manifest['files'] if s['dsid'] == 700325)
    a = uproot.open(BASE / 'data' / sample['key'])['analysis'].arrays(
        BRANCHES, entry_start=6869582, entry_stop=6869583, library='ak')
    assert .008 < validate_chunk(a, sample) < .01
    assert all(not mask.any() for mask in select_events(a)[1].values())
    from audit_met_inputs import interval_gap
    x = np.array([1e6], dtype=np.float32)
    radius = np.hypot(x, x)
    assert abs(np.hypot(x.astype(float), x.astype(float))-radius)[0] > .01
    gap, arithmetic = interval_gap(radius, x, x)
    assert gap[0] <= arithmetic[0]  # Un redondeo válido no es un fallo físico.


@pytest.mark.parametrize('etas, expected', [([.4, .4], 0.), ([.4, -.4], .8),
                                         ([-2.4, 2.49], 4.89), ([2.4, -2.49], 4.89)])
def test_abs_delta_eta_geometry_and_lepton_order(etas, expected):
    """Diferencia absoluta, reflexión del haz y orden de los leptones."""
    for order in ([0, 1], [1, 0]):
        row = example()
        row['lep_eta'] = etas
        row['lep_e'] = (np.array(row['lep_pt']) * np.cosh(etas)).tolist()
        for name, values in list(row.items()):
            if name.startswith('lep_'):
                row[name] = [values[i] for i in order]
        row['met'] = 80.
        a = ak.Array([row])
        _, regions, ele, mu, nj, nb = select_events(a)
        obs, _ = observables(a, regions['SR2b'], ele, mu, nj, nb)
        assert obs['delta_eta'].shape == (1,)
        assert obs['delta_eta'][0] == pytest.approx(expected)


def test_pipeline_fingerprint_covers_imported_numerical_checks():
    from analysis import pipeline_hashes
    hashes = pipeline_hashes()
    assert set(hashes) == {'analysis.py', 'audit_met_inputs.py', 'fetch_inputs.py'}
    assert all(len(value) == 64 for value in hashes.values())
