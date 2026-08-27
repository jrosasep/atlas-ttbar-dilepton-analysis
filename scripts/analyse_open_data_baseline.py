"""Run a transparent e-mu selection on the public ATLAS ttbar MC sample.

This script intentionally stops before reconstructing the two neutrinos. Its
outputs are lab-frame control distributions, not the published ATLAS D value.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


BRANCHES = [
    "mcWeight",
    "lep_type",
    "lep_pt",
    "lep_eta",
    "lep_phi",
    "lep_e",
    "lep_charge",
    "lep_isTightID",
    "lep_isTightIso",
    "jet_pt",
    "jet_eta",
    "jet_btag_quantile",
    "met",
]


def absolute_delta_phi(phi_a: np.ndarray, phi_b: np.ndarray) -> np.ndarray:
    difference = (phi_a - phi_b + np.pi) % (2.0 * np.pi) - np.pi
    return np.abs(difference)


def three_momentum(
    pt: np.ndarray, eta: np.ndarray, phi: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)


def process_chunk(arrays: ak.Array, btag_quantile_min: int) -> dict[str, object]:
    lepton_pt = arrays["lep_pt"]
    lepton_eta = arrays["lep_eta"]
    lepton_type = abs(arrays["lep_type"])
    selected_lepton = (
        (lepton_pt >= 25.0)
        & arrays["lep_isTightID"]
        & arrays["lep_isTightIso"]
    )
    electron_acceptance = (abs(lepton_eta) < 1.37) | (
        (abs(lepton_eta) > 1.52) & (abs(lepton_eta) < 2.47)
    )
    muon_acceptance = abs(lepton_eta) < 2.5
    electron = selected_lepton & (lepton_type == 11) & electron_acceptance
    muon = selected_lepton & (lepton_type == 13) & muon_acceptance

    emu = (ak.sum(electron, axis=1) == 1) & (ak.sum(muon, axis=1) == 1)
    selected_charge_sum = ak.sum(arrays["lep_charge"] * (electron | muon), axis=1)
    opposite_sign = emu & (selected_charge_sum == 0)

    selected_jet = (arrays["jet_pt"] >= 25.0) & (abs(arrays["jet_eta"]) < 2.5)
    at_least_two_jets = opposite_sign & (ak.sum(selected_jet, axis=1) >= 2)
    selected_bjet = selected_jet & (
        arrays["jet_btag_quantile"] >= btag_quantile_min
    )
    final = at_least_two_jets & (ak.sum(selected_bjet, axis=1) >= 1)

    counts = {
        "input": len(arrays),
        "exactly_one_selected_e_and_mu": int(ak.sum(emu)),
        "opposite_sign": int(ak.sum(opposite_sign)),
        "at_least_two_jets": int(ak.sum(at_least_two_jets)),
        "at_least_one_btag_85": int(ak.sum(final)),
    }

    electron_pt = ak.to_numpy(ak.firsts(arrays["lep_pt"][electron][final]))
    electron_eta = ak.to_numpy(ak.firsts(arrays["lep_eta"][electron][final]))
    electron_phi = ak.to_numpy(ak.firsts(arrays["lep_phi"][electron][final]))
    electron_energy = ak.to_numpy(ak.firsts(arrays["lep_e"][electron][final]))
    muon_pt = ak.to_numpy(ak.firsts(arrays["lep_pt"][muon][final]))
    muon_eta = ak.to_numpy(ak.firsts(arrays["lep_eta"][muon][final]))
    muon_phi = ak.to_numpy(ak.firsts(arrays["lep_phi"][muon][final]))
    muon_energy = ak.to_numpy(ak.firsts(arrays["lep_e"][muon][final]))

    e_px, e_py, e_pz = three_momentum(electron_pt, electron_eta, electron_phi)
    m_px, m_py, m_pz = three_momentum(muon_pt, muon_eta, muon_phi)
    e_norm = np.sqrt(e_px**2 + e_py**2 + e_pz**2)
    m_norm = np.sqrt(m_px**2 + m_py**2 + m_pz**2)
    cos_opening = (e_px * m_px + e_py * m_py + e_pz * m_pz) / (e_norm * m_norm)
    mass_squared = (
        (electron_energy + muon_energy) ** 2
        - (e_px + m_px) ** 2
        - (e_py + m_py) ** 2
        - (e_pz + m_pz) ** 2
    )

    return {
        "counts": counts,
        "delta_phi": absolute_delta_phi(electron_phi, muon_phi),
        "cos_opening_lab": np.clip(cos_opening, -1.0, 1.0),
        "m_emu": np.sqrt(np.clip(mass_squared, 0.0, None)),
        "met": ak.to_numpy(arrays["met"][final]),
        "weights": ak.to_numpy(arrays["mcWeight"][final]),
    }


def add_counts(total: dict[str, int], update: dict[str, int]) -> None:
    for key, value in update.items():
        total[key] = total.get(key, 0) + value


def normalized_histogram(
    values: np.ndarray, weights: np.ndarray, bins: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    histogram, _ = np.histogram(values, bins=bins, weights=weights)
    sum_weights_squared, _ = np.histogram(values, bins=bins, weights=weights**2)
    widths = np.diff(bins)
    normalization = float(np.sum(histogram))
    if normalization <= 0:
        raise ValueError("selected generator weights do not have a positive sum")
    return histogram / (normalization * widths), np.sqrt(sum_weights_squared) / (
        normalization * widths
    )


def draw_distribution(
    ax: plt.Axes,
    values: np.ndarray,
    weights: np.ndarray,
    bins: np.ndarray,
    xlabel: str,
) -> None:
    density, uncertainty = normalized_histogram(values, weights, bins)
    centres = 0.5 * (bins[:-1] + bins[1:])
    ax.stairs(density, bins, color="#1d4ed8", linewidth=1.8)
    ax.errorbar(
        centres,
        density,
        yerr=uncertainty,
        fmt="none",
        ecolor="#1e3a8a",
        elinewidth=0.8,
        alpha=0.75,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized events")
    ax.grid(alpha=0.22)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--tree", default="analysis")
    parser.add_argument("--step-size", default="100 MB")
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--btag-quantile-min", type=int, default=2)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    project = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir or project / "results"
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    cutflow: dict[str, int] = {}
    collected: dict[str, list[np.ndarray]] = {
        "delta_phi": [],
        "cos_opening_lab": [],
        "m_emu": [],
        "met": [],
        "weights": [],
    }
    processed = 0
    source = f"{args.root_file}:{args.tree}"
    for arrays in uproot.iterate(
        source,
        expressions=BRANCHES,
        step_size=args.step_size,
        library="ak",
    ):
        if args.max_events is not None:
            remaining = args.max_events - processed
            if remaining <= 0:
                break
            arrays = arrays[:remaining]
        chunk = process_chunk(arrays, args.btag_quantile_min)
        add_counts(cutflow, chunk["counts"])
        for key in collected:
            collected[key].append(chunk[key])
        processed += len(arrays)

    values = {key: np.concatenate(parts) for key, parts in collected.items()}
    weights = values["weights"].astype(float)
    if values["delta_phi"].size == 0:
        raise RuntimeError("no events passed the selection")

    summary = {
        "input_file": str(args.root_file),
        "tree": args.tree,
        "events_processed": processed,
        "selection": {
            "leptons": (
                "exactly one tight+isolated electron and one tight+isolated muon, "
                "pT >= 25 GeV, ATLAS eta acceptance, opposite charge"
            ),
            "jets": "at least two jets with pT >= 25 GeV and |eta| < 2.5",
            "b_tag": (
                "at least one jet with DL1dv01 continuous quantile >= "
                f"{args.btag_quantile_min} (85% working-point bin)"
            ),
            "trigger": "not emulated in this baseline",
        },
        "cutflow_unweighted": cutflow,
        "selected_fraction": cutflow["at_least_one_btag_85"] / cutflow["input"],
        "control_observables": {
            "mean_abs_delta_phi": float(np.average(values["delta_phi"], weights=weights)),
            "mean_cos_opening_lab": float(
                np.average(values["cos_opening_lab"], weights=weights)
            ),
            "mean_m_emu_gev": float(np.average(values["m_emu"], weights=weights)),
            "mean_met_gev": float(np.average(values["met"], weights=weights)),
        },
        "physics_status": (
            "These are reconstructed lab-frame controls in simulated ttbar events. "
            "They are not D and are not corrected to particle level."
        ),
    }
    (output_dir / "open_data_baseline_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    draw_distribution(
        axes[0], values["delta_phi"], weights, np.linspace(0, np.pi, 17),
        "|Delta phi(e, mu)| [rad]",
    )
    draw_distribution(
        axes[1], values["cos_opening_lab"], weights, np.linspace(-1, 1, 17),
        "cos(opening angle) in lab",
    )
    draw_distribution(
        axes[2], values["m_emu"], weights, np.linspace(0, 300, 21),
        "m(e, mu) [GeV]",
    )
    fig.suptitle(
        "ATLAS Open Data ttbar MC: selected e-mu control distributions",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.005,
        "2015-2016 educational MC; generator-weighted shapes; not the particle-level D measurement",
        ha="center",
        fontsize=8,
        color="#4b5563",
    )
    fig.savefig(figures_dir / "open_data_ttbar_baseline.png", dpi=180)
    fig.savefig(figures_dir / "open_data_ttbar_baseline.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
