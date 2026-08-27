import numpy as np
import pytest

from ttbar_entanglement.reconstruction import (
    invariant_mass,
    massless_four_vector,
    reconstruct_dilepton_event,
)


def boost_from_rest(vector: np.ndarray, beta: np.ndarray) -> np.ndarray:
    beta2 = float(np.dot(beta, beta))
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    energy = gamma * (vector[0] + np.dot(beta, vector[1:]))
    if beta2 == 0:
        return vector.copy()
    space = vector[1:] + (
        (gamma - 1.0) * np.dot(beta, vector[1:]) / beta2 + gamma * vector[0]
    ) * beta
    return np.concatenate(([energy], space))


def synthetic_ttbar_at_rest() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mt = 172.5
    mw = 80.4
    top_decay_p = (mt**2 - mw**2) / (2.0 * mt)
    w_energy = (mt**2 + mw**2) / (2.0 * mt)
    lepton_w = np.array([mw / 2.0, mw / 2.0, 0.0, 0.0])
    neutrino_w = np.array([mw / 2.0, -mw / 2.0, 0.0, 0.0])

    beta_plus = np.array([0.0, 0.0, top_decay_p / w_energy])
    beta_minus = -beta_plus
    lepton_plus = boost_from_rest(lepton_w, beta_plus)
    neutrino = boost_from_rest(neutrino_w, beta_plus)
    lepton_minus = boost_from_rest(lepton_w, beta_minus)
    antineutrino = boost_from_rest(neutrino_w, beta_minus)
    b = np.array([top_decay_p, 0.0, 0.0, -top_decay_p])
    antib = np.array([top_decay_p, 0.0, 0.0, top_decay_p])
    leptons = np.stack([lepton_plus, lepton_minus])
    neutrinos = np.stack([neutrino, antineutrino])
    jets = np.stack([b, antib])
    charges = np.array([1, -1])
    return leptons, charges, jets, neutrinos


def test_numeric_reconstruction_closes_on_exact_event() -> None:
    leptons, charges, jets, neutrinos = synthetic_ttbar_at_rest()
    met = np.sum(neutrinos[:, 1:3], axis=0)
    result = reconstruct_dilepton_event(leptons, charges, jets, met[0], met[1])
    assert result is not None
    assert invariant_mass(result.top) == pytest.approx(172.5, abs=0.05)
    assert invariant_mass(result.antitop) == pytest.approx(172.5, abs=0.05)
    assert result.maximum_scaled_residual < 2.0e-3
    assert -1.0 <= result.cos_phi <= 1.0


def test_massless_four_vector() -> None:
    vector = massless_four_vector(3.0, 4.0, 12.0)
    assert invariant_mass(vector) == pytest.approx(0.0, abs=1e-7)
