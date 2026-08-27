"""Prototype kinematic reconstruction for dileptonic ttbar events.

The published ATLAS analysis uses an analytic Ellipse method. This module uses
the same on-shell mass and missing-transverse-momentum constraints, solved
numerically, as an independently testable development implementation. It is
not presented as a bit-for-bit reproduction of the ATLAS reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .observable import cosphi_from_four_vectors


@dataclass(frozen=True)
class DileptonSolution:
    top: NDArray[np.float64]
    antitop: NDArray[np.float64]
    neutrino_from_top: NDArray[np.float64]
    neutrino_from_antitop: NDArray[np.float64]
    m_ttbar: float
    cos_phi: float
    maximum_scaled_residual: float
    pairing: tuple[int, int]
    function_evaluations: int


def four_vector_from_pt_eta_phi_e(
    pt: float, eta: float, phi: float, energy: float
) -> NDArray[np.float64]:
    return np.array(
        [energy, pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)],
        dtype=float,
    )


def invariant_mass_squared(vector: ArrayLike) -> float:
    p = np.asarray(vector, dtype=float)
    return float(p[0] ** 2 - np.dot(p[1:], p[1:]))


def invariant_mass(vector: ArrayLike) -> float:
    return float(np.sqrt(max(invariant_mass_squared(vector), 0.0)))


def massless_four_vector(px: float, py: float, pz: float) -> NDArray[np.float64]:
    return np.array([np.sqrt(px**2 + py**2 + pz**2), px, py, pz], dtype=float)


def _scaled_residuals(
    variables: NDArray[np.float64],
    lepton_a: NDArray[np.float64],
    lepton_b: NDArray[np.float64],
    bjet_a: NDArray[np.float64],
    bjet_b: NDArray[np.float64],
    met_x: float,
    met_y: float,
    w_mass: float,
    top_mass: float,
) -> NDArray[np.float64]:
    nu_a = massless_four_vector(variables[0], variables[1], variables[2])
    nu_b = massless_four_vector(
        met_x - variables[0], met_y - variables[1], variables[3]
    )
    return np.array(
        [
            (invariant_mass_squared(lepton_a + nu_a) - w_mass**2) / w_mass**2,
            (invariant_mass_squared(lepton_b + nu_b) - w_mass**2) / w_mass**2,
            (
                invariant_mass_squared(lepton_a + nu_a + bjet_a) - top_mass**2
            )
            / top_mass**2,
            (
                invariant_mass_squared(lepton_b + nu_b + bjet_b) - top_mass**2
            )
            / top_mass**2,
        ],
        dtype=float,
    )


def default_starting_points(met_x: float, met_y: float) -> list[NDArray[np.float64]]:
    """Deterministic starting points spanning common neutrino pz solutions."""

    starts = []
    for fraction, pz_a, pz_b in [
        (0.50, 0.0, 0.0),
        (0.50, 100.0, -100.0),
        (0.50, -100.0, 100.0),
        (0.25, 50.0, 50.0),
        (0.75, -50.0, -50.0),
        (0.50, 250.0, -250.0),
    ]:
        starts.append(
            np.array([fraction * met_x, fraction * met_y, pz_a, pz_b])
        )
    return starts


def reconstruct_dilepton_event(
    leptons: ArrayLike,
    charges: ArrayLike,
    bjets: ArrayLike,
    met_x: float,
    met_y: float,
    *,
    w_mass: float = 80.4,
    top_mass: float = 172.5,
    maximum_scaled_residual: float = 2.0e-3,
    max_function_evaluations: int = 300,
) -> DileptonSolution | None:
    """Solve the dilepton mass constraints and choose the lowest-mtt solution.

    Both assignments of the two selected jets to the two leptons are tried.
    A positive charged lepton is associated with the top quark, and a negative
    charged lepton with the antitop quark.
    """

    lepton_vectors = np.asarray(leptons, dtype=float)
    lepton_charges = np.asarray(charges, dtype=int)
    jet_vectors = np.asarray(bjets, dtype=float)
    if lepton_vectors.shape != (2, 4) or jet_vectors.shape != (2, 4):
        raise ValueError("exactly two lepton and two jet four-vectors are required")
    if lepton_charges.shape != (2,) or set(lepton_charges.tolist()) != {-1, 1}:
        raise ValueError("lepton charges must be one +1 and one -1")

    candidates: list[DileptonSolution] = []
    starts = default_starting_points(met_x, met_y)
    for pairing in [(0, 1), (1, 0)]:
        bjet_a = jet_vectors[pairing[0]]
        bjet_b = jet_vectors[pairing[1]]
        for starting_point in starts:
            fit = least_squares(
                _scaled_residuals,
                starting_point,
                args=(
                    lepton_vectors[0],
                    lepton_vectors[1],
                    bjet_a,
                    bjet_b,
                    met_x,
                    met_y,
                    w_mass,
                    top_mass,
                ),
                method="lm",
                xtol=1e-10,
                ftol=1e-10,
                gtol=1e-10,
                max_nfev=max_function_evaluations,
            )
            residuals = _scaled_residuals(
                fit.x,
                lepton_vectors[0],
                lepton_vectors[1],
                bjet_a,
                bjet_b,
                met_x,
                met_y,
                w_mass,
                top_mass,
            )
            quality = float(np.max(np.abs(residuals)))
            if quality > maximum_scaled_residual:
                continue

            nu_a = massless_four_vector(fit.x[0], fit.x[1], fit.x[2])
            nu_b = massless_four_vector(met_x - fit.x[0], met_y - fit.x[1], fit.x[3])
            parent_a = lepton_vectors[0] + nu_a + bjet_a
            parent_b = lepton_vectors[1] + nu_b + bjet_b
            ttbar = parent_a + parent_b

            if lepton_charges[0] > 0:
                top, antitop = parent_a, parent_b
                lepton_top, lepton_antitop = lepton_vectors[0], lepton_vectors[1]
                nu_top, nu_antitop = nu_a, nu_b
            else:
                top, antitop = parent_b, parent_a
                lepton_top, lepton_antitop = lepton_vectors[1], lepton_vectors[0]
                nu_top, nu_antitop = nu_b, nu_a

            cos_phi = float(
                cosphi_from_four_vectors(top, antitop, lepton_top, lepton_antitop)
            )
            candidates.append(
                DileptonSolution(
                    top=top,
                    antitop=antitop,
                    neutrino_from_top=nu_top,
                    neutrino_from_antitop=nu_antitop,
                    m_ttbar=invariant_mass(ttbar),
                    cos_phi=cos_phi,
                    maximum_scaled_residual=quality,
                    pairing=pairing,
                    function_evaluations=int(fit.nfev),
                )
            )

    if not candidates:
        return None
    return min(candidates, key=lambda candidate: candidate.m_ttbar)
