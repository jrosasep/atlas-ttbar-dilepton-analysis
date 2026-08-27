import numpy as np
import pytest

from ttbar_entanglement.observable import (
    analytic_d_uncertainty,
    boost_to_rest,
    combine_in_quadrature,
    cosphi_from_four_vectors,
    estimate_d,
)
from ttbar_entanglement.toy import angular_pdf, reweight_d, sample_cos_phi


def test_quadrature() -> None:
    assert combine_in_quadrature(3.0, 4.0) == pytest.approx(5.0)


def test_moment_estimator() -> None:
    assert estimate_d([-1.0, 0.0, 1.0]) == pytest.approx(0.0)
    assert estimate_d([0.2, 0.4], weights=[1.0, 3.0]) == pytest.approx(-1.05)
    assert estimate_d([0.0, 1.0], weights=[2.0, -1.0]) == pytest.approx(3.0)


def test_pdf_is_normalized_and_has_correct_mean() -> None:
    grid = np.linspace(-1.0, 1.0, 100_001)
    d_value = -0.537
    density = angular_pdf(grid, d_value)
    assert np.trapezoid(density, grid) == pytest.approx(1.0, abs=1e-9)
    mean = np.trapezoid(grid * density, grid)
    assert -3.0 * mean == pytest.approx(d_value, abs=1e-9)


def test_sampling_and_reweighting_close() -> None:
    rng = np.random.default_rng(12345)
    sample = sample_cos_phi(-0.47, 400_000, rng)
    assert estimate_d(sample) == pytest.approx(-0.47, abs=0.01)
    weights = reweight_d(sample, -0.47, -0.537)
    assert estimate_d(sample, weights) == pytest.approx(-0.537, abs=0.01)


def test_analytic_uncertainty() -> None:
    assert analytic_d_uncertainty(0.0, 300) == pytest.approx(0.1)


def test_boost_parent_to_rest() -> None:
    mass = 2.0
    momentum = 3.0
    energy = np.sqrt(mass**2 + momentum**2)
    parent = np.array([energy, momentum, 0.0, 0.0])
    boosted = boost_to_rest(parent, parent)
    assert boosted == pytest.approx([mass, 0.0, 0.0, 0.0], abs=1e-12)


def test_cosphi_simple_ttbar_rest_frame() -> None:
    top = np.array([172.5, 0.0, 0.0, 0.0])
    antitop = np.array([172.5, 0.0, 0.0, 0.0])
    lepton_a = np.array([10.0, 10.0, 0.0, 0.0])
    lepton_b = np.array([10.0, 0.0, 10.0, 0.0])
    assert cosphi_from_four_vectors(top, antitop, lepton_a, lepton_b) == pytest.approx(0.0)
