import numpy as np
import pytest

from scripts.analyse_open_data_baseline import absolute_delta_phi, three_momentum


def test_absolute_delta_phi_wraps_at_pi() -> None:
    result = absolute_delta_phi(np.array([3.0]), np.array([-3.0]))
    assert result[0] == pytest.approx(2.0 * np.pi - 6.0)


def test_three_momentum_at_zero_eta_phi() -> None:
    px, py, pz = three_momentum(
        np.array([25.0]), np.array([0.0]), np.array([0.0])
    )
    assert px == pytest.approx([25.0])
    assert py == pytest.approx([0.0])
    assert pz == pytest.approx([0.0])
