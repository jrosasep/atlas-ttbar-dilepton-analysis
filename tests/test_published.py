from pathlib import Path

import pytest

from ttbar_entanglement.published import (
    gaussian_distance_from_boundary,
    load_config,
    measurement_from_dict,
    systematic_quadrature,
)


PROJECT = Path(__file__).resolve().parents[1]


def test_published_uncertainty_budget_rounds_to_quoted_value() -> None:
    config = load_config(PROJECT / "config" / "published_atlas.yaml")
    calculated = systematic_quadrature(config)
    assert calculated == pytest.approx(0.018520259, rel=1e-6)
    assert round(calculated, 3) == 0.019


def test_signal_result_is_well_beyond_five_sigma_in_simple_cross_check() -> None:
    config = load_config(PROJECT / "config" / "published_atlas.yaml")
    observed = measurement_from_dict(config["regions"]["signal"]["observed"])
    boundary = config["particle_level_boundaries"]["powheg_pythia"]
    distance = gaussian_distance_from_boundary(
        observed, boundary["value"], boundary["uncertainty"]
    )
    assert distance > 5.0
