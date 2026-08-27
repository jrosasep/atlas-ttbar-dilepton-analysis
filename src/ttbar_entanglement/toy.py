"""Toy Monte Carlo for the normalized cos(phi) distribution."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .observable import estimate_d


def angular_pdf(cos_phi: ArrayLike, d_value: float) -> NDArray[np.float64]:
    """Evaluate p(x|D) = 0.5*(1-D*x) for x in [-1, 1]."""

    if abs(d_value) > 1:
        raise ValueError("the linear angular PDF requires |D| <= 1")
    x = np.asarray(cos_phi, dtype=float)
    if np.any(np.abs(x) > 1):
        raise ValueError("cos_phi must lie in [-1, 1]")
    return 0.5 * (1.0 - d_value * x)


def sample_cos_phi(
    d_value: float, n_events: int, rng: np.random.Generator
) -> NDArray[np.float64]:
    """Draw exact samples by rejection sampling from the linear angular PDF."""

    if abs(d_value) > 1:
        raise ValueError("the linear angular PDF requires |D| <= 1")
    if n_events <= 0:
        raise ValueError("n_events must be positive")

    accepted: list[NDArray[np.float64]] = []
    remaining = n_events
    envelope = 1.0 + abs(d_value)
    while remaining > 0:
        batch_size = max(1024, 2 * remaining)
        candidates = rng.uniform(-1.0, 1.0, batch_size)
        keep = rng.random(batch_size) < (1.0 - d_value * candidates) / envelope
        selected = candidates[keep][:remaining]
        accepted.append(selected)
        remaining -= selected.size
    return np.concatenate(accepted)


def reweight_d(
    cos_phi: ArrayLike, nominal_d: float, target_d: float
) -> NDArray[np.float64]:
    """Return event weights that change the angular slope from nominal to target."""

    if abs(nominal_d) > 1 or abs(target_d) > 1:
        raise ValueError("both D values must satisfy |D| <= 1")
    x = np.asarray(cos_phi, dtype=float)
    denominator = 1.0 - nominal_d * x
    if np.any(denominator <= 0):
        raise ValueError("nominal density is zero for at least one event")
    return (1.0 - target_d * x) / denominator


def closure_trials(
    d_value: float,
    n_events: int,
    n_trials: int,
    seed: int,
) -> NDArray[np.float64]:
    """Run repeated pseudo-experiments and return their D estimates."""

    if n_trials <= 0:
        raise ValueError("n_trials must be positive")
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_trials, dtype=float)
    for index in range(n_trials):
        estimates[index] = estimate_d(sample_cos_phi(d_value, n_events, rng))
    return estimates
