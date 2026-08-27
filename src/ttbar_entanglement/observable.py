"""Definition and estimators for the ATLAS entanglement marker D."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray


def combine_in_quadrature(*uncertainties: float) -> float:
    """Return the quadrature sum of independent non-negative uncertainties."""

    values = np.asarray(uncertainties, dtype=float)
    if np.any(values < 0) or np.any(~np.isfinite(values)):
        raise ValueError("uncertainties must be finite and non-negative")
    return float(np.sqrt(np.sum(values**2)))


def estimate_d(cos_phi: ArrayLike, weights: ArrayLike | None = None) -> float:
    """Estimate D using the moment estimator D = -3 <cos(phi)>.

    Parameters
    ----------
    cos_phi:
        One-dimensional values in [-1, 1].
    weights:
        Optional finite event weights. Signed weights from NLO Monte Carlo are
        supported as long as their total is non-zero.
    """

    x = np.asarray(cos_phi, dtype=float)
    if x.ndim != 1 or x.size == 0:
        raise ValueError("cos_phi must be a non-empty one-dimensional array")
    if np.any(~np.isfinite(x)) or np.any(np.abs(x) > 1.0 + 1e-12):
        raise ValueError("cos_phi values must be finite and lie in [-1, 1]")

    if weights is None:
        mean = float(np.mean(x))
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != x.shape:
            raise ValueError("weights and cos_phi must have the same shape")
        if np.any(~np.isfinite(w)) or abs(float(np.sum(w))) <= 1e-15:
            raise ValueError("weights must be finite and have a non-zero sum")
        mean = float(np.average(x, weights=w))
    return -3.0 * mean


def estimate_d_uncertainty(
    cos_phi: ArrayLike, weights: ArrayLike | None = None
) -> float:
    """Estimate the statistical standard error of the moment estimator.

    For unweighted samples this is `3 * s / sqrt(N)`. For weighted samples,
    the unbiased weighted variance and the Kish effective sample size are used.
    """

    x = np.asarray(cos_phi, dtype=float)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("at least two one-dimensional observations are required")

    if weights is None:
        return float(3.0 * np.std(x, ddof=1) / math.sqrt(x.size))

    w = np.asarray(weights, dtype=float)
    if w.shape != x.shape or np.any(~np.isfinite(w)):
        raise ValueError("valid weights with the same shape as cos_phi are required")
    sum_w = float(np.sum(w))
    sum_w2 = float(np.sum(w**2))
    if abs(sum_w) <= 1e-15:
        raise ValueError("weights must have a non-zero sum")
    mean = float(np.average(x, weights=w))

    # NLO Monte Carlo can contain negative event weights. A Kish effective
    # sample size is not defined for that case, so use the standard sum(w^2)
    # propagation for the weighted mean.
    if np.any(w < 0):
        variance_of_mean = float(np.sum(w**2 * (x - mean) ** 2) / sum_w**2)
        return 3.0 * math.sqrt(variance_of_mean)

    if sum_w**2 <= sum_w2:
        raise ValueError("weights do not define at least two effective observations")
    variance_unbiased = float(
        np.sum(w * (x - mean) ** 2) / (sum_w - sum_w2 / sum_w)
    )
    n_eff = sum_w**2 / sum_w2
    return 3.0 * math.sqrt(variance_unbiased / n_eff)


def analytic_d_uncertainty(d_value: float, n_events: int) -> float:
    """Analytic standard error for samples from 0.5*(1-D*x)."""

    if abs(d_value) > 1:
        raise ValueError("the linear angular PDF requires |D| <= 1")
    if n_events <= 0:
        raise ValueError("n_events must be positive")
    return math.sqrt((3.0 - d_value**2) / n_events)


def boost_to_rest(vectors: ArrayLike, parent: ArrayLike) -> NDArray[np.float64]:
    """Lorentz-boost four-vectors into the rest frame of `parent`.

    Four-vectors use the order `(E, px, py, pz)` and may be a single vector or
    arrays with matching leading dimensions. Natural units with c=1 are used.
    """

    v = np.asarray(vectors, dtype=float)
    p = np.asarray(parent, dtype=float)
    if v.shape[-1] != 4 or p.shape[-1] != 4:
        raise ValueError("four-vectors must end in four components (E, px, py, pz)")

    energy = p[..., :1]
    if np.any(energy <= 0):
        raise ValueError("parent energies must be positive")
    beta = p[..., 1:] / energy
    beta2 = np.sum(beta**2, axis=-1, keepdims=True)
    if np.any(beta2 >= 1):
        raise ValueError("parent four-vector must be timelike")

    gamma = 1.0 / np.sqrt(1.0 - beta2)
    v_energy = v[..., :1]
    v_space = v[..., 1:]
    beta_dot_v = np.sum(beta * v_space, axis=-1, keepdims=True)

    safe_beta2 = np.where(beta2 > 0, beta2, 1.0)
    factor = ((gamma - 1.0) * beta_dot_v / safe_beta2) - gamma * v_energy
    boosted_space = v_space + factor * beta
    boosted_energy = gamma * (v_energy - beta_dot_v)
    return np.concatenate((boosted_energy, boosted_space), axis=-1)


def cosphi_from_four_vectors(
    top: ArrayLike,
    antitop: ArrayLike,
    lepton_from_top: ArrayLike,
    lepton_from_antitop: ArrayLike,
) -> NDArray[np.float64]:
    """Compute cos(phi) following the two-stage rest-frame construction.

    All input four-vectors are provided in one common frame. They are first
    boosted into the ttbar rest frame. Each charged lepton is then boosted into
    the rest frame of its parent top or antitop. The returned value is the dot
    product of the two resulting unit directions.
    """

    t = np.asarray(top, dtype=float)
    tb = np.asarray(antitop, dtype=float)
    lp = np.asarray(lepton_from_top, dtype=float)
    lm = np.asarray(lepton_from_antitop, dtype=float)
    ttbar = t + tb

    t_tt = boost_to_rest(t, ttbar)
    tb_tt = boost_to_rest(tb, ttbar)
    lp_tt = boost_to_rest(lp, ttbar)
    lm_tt = boost_to_rest(lm, ttbar)

    lp_top = boost_to_rest(lp_tt, t_tt)
    lm_antitop = boost_to_rest(lm_tt, tb_tt)
    a = lp_top[..., 1:]
    b = lm_antitop[..., 1:]
    norm_a = np.linalg.norm(a, axis=-1)
    norm_b = np.linalg.norm(b, axis=-1)
    if np.any(norm_a == 0) or np.any(norm_b == 0):
        raise ValueError("lepton directions must have non-zero spatial momentum")
    result = np.sum(a * b, axis=-1) / (norm_a * norm_b)
    return np.clip(result, -1.0, 1.0)
