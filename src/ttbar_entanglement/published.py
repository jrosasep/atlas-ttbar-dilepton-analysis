"""Load and independently check the published ATLAS values."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .observable import combine_in_quadrature


@dataclass(frozen=True)
class Measurement:
    value: float
    stat: float
    syst: float

    @property
    def total(self) -> float:
        return combine_in_quadrature(self.stat, self.syst)


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError("published configuration must be a mapping")
    return config


def measurement_from_dict(values: dict[str, float]) -> Measurement:
    return Measurement(
        value=float(values["value"]),
        stat=float(values["stat"]),
        syst=float(values["syst"]),
    )


def gaussian_distance_from_boundary(
    measurement: Measurement,
    boundary_value: float,
    boundary_uncertainty: float = 0.0,
) -> float:
    """Simple Gaussian distance below a no-entanglement boundary.

    This is an arithmetic cross-check, not the official ATLAS likelihood-based
    significance. Positive values indicate that the measurement lies below
    the boundary.
    """

    denominator = combine_in_quadrature(measurement.total, boundary_uncertainty)
    return (boundary_value - measurement.value) / denominator


def systematic_quadrature(config: dict[str, Any]) -> float:
    values = [float(value) for value in config["signal_systematics"].values()]
    return combine_in_quadrature(*values)
