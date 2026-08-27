"""Prototype two-neutrino reconstruction on selected public ttbar MC events."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src"))

from ttbar_entanglement.observable import estimate_d, estimate_d_uncertainty  # noqa: E402
from ttbar_entanglement.reconstruction import (  # noqa: E402
    four_vector_from_pt_eta_phi_e,
    reconstruct_dilepton_event,
)


BRANCHES = [
    "eventNumber",
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
    "jet_phi",
    "jet_e",
    "jet_btag_quantile",
    "met_mpx",
    "met_mpy",
]


def event_selection(arrays: ak.Array, btag_quantile_min: int) -> ak.Array:
    lepton_pt = arrays["lep_pt"]
    lepton_eta = arrays["lep_eta"]
    lepton_type = abs(arrays["lep_type"])
    tight = (
        (lepton_pt >= 25.0)
        & arrays["lep_isTightID"]
        & arrays["lep_isTightIso"]
    )
    electron = (
        tight
        & (lepton_type == 11)
        & ((abs(lepton_eta) < 1.37) | ((abs(lepton_eta) > 1.52) & (abs(lepton_eta) < 2.47)))
    )
    muon = tight & (lepton_type == 13) & (abs(lepton_eta) < 2.5)
    emu = (ak.sum(electron, axis=1) == 1) & (ak.sum(muon, axis=1) == 1)
    opposite_sign = emu & (ak.sum(arrays["lep_charge"] * (electron | muon), axis=1) == 0)
    jet = (arrays["jet_pt"] >= 25.0) & (abs(arrays["jet_eta"]) < 2.5)
    bjet = jet & (arrays["jet_btag_quantile"] >= btag_quantile_min)
    return opposite_sign & (ak.sum(jet, axis=1) >= 2) & (ak.sum(bjet, axis=1) >= 1)


def select_two_jets(event: ak.Record, btag_quantile_min: int) -> list[int]:
    pt = np.asarray(event["jet_pt"])
    eta = np.asarray(event["jet_eta"])
    quantile = np.asarray(event["jet_btag_quantile"])
    accepted = np.flatnonzero((pt >= 25.0) & (np.abs(eta) < 2.5))
    tagged = accepted[quantile[accepted] >= btag_quantile_min]
    tagged = tagged[np.argsort(pt[tagged])[::-1]]
    if tagged.size >= 2:
        return [int(tagged[0]), int(tagged[1])]
    untagged = accepted[accepted != tagged[0]]
    untagged = untagged[np.argsort(pt[untagged])[::-1]]
    return [int(tagged[0]), int(untagged[0])]


def build_four_vectors(event: ak.Record, indices: list[int], prefix: str) -> np.ndarray:
    vectors = []
    for index in indices:
        vectors.append(
            four_vector_from_pt_eta_phi_e(
                float(event[f"{prefix}_pt"][index]),
                float(event[f"{prefix}_eta"][index]),
                float(event[f"{prefix}_phi"][index]),
                float(event[f"{prefix}_e"][index]),
            )
        )
    return np.asarray(vectors)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--tree", default="analysis")
    parser.add_argument("--max-selected", type=int, default=500)
    parser.add_argument("--btag-quantile-min", type=int, default=2)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    output_dir = args.output_dir or PROJECT / "results"
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    chosen_events: list[ak.Record] = []
    for arrays in uproot.iterate(
        f"{args.root_file}:{args.tree}",
        expressions=BRANCHES,
        step_size="80 MB",
        library="ak",
    ):
        selected = arrays[event_selection(arrays, args.btag_quantile_min)]
        needed = args.max_selected - len(chosen_events)
        chosen_events.extend(ak.to_list(selected[:needed]))
        if len(chosen_events) >= args.max_selected:
            break

    reconstructed_mtt: list[float] = []
    reconstructed_cosphi: list[float] = []
    reconstructed_weights: list[float] = []
    qualities: list[float] = []
    failures = 0
    for raw_event in chosen_events:
        event = ak.Record(raw_event)
        lepton_indices = list(range(len(event["lep_pt"])))
        accepted_leptons = [
            index
            for index in lepton_indices
            if float(event["lep_pt"][index]) >= 25.0
            and bool(event["lep_isTightID"][index])
            and bool(event["lep_isTightIso"][index])
            and (
                (
                    abs(int(event["lep_type"][index])) == 11
                    and (
                        abs(float(event["lep_eta"][index])) < 1.37
                        or 1.52 < abs(float(event["lep_eta"][index])) < 2.47
                    )
                )
                or (
                    abs(int(event["lep_type"][index])) == 13
                    and abs(float(event["lep_eta"][index])) < 2.5
                )
            )
        ]
        # event_selection guarantees one accepted electron and one accepted muon.
        lepton_vectors = build_four_vectors(event, accepted_leptons, "lep")
        charges = np.asarray([int(event["lep_charge"][i]) for i in accepted_leptons])
        jet_indices = select_two_jets(event, args.btag_quantile_min)
        jet_vectors = build_four_vectors(event, jet_indices, "jet")

        solution = reconstruct_dilepton_event(
            lepton_vectors,
            charges,
            jet_vectors,
            float(event["met_mpx"]),
            float(event["met_mpy"]),
        )
        if solution is None:
            failures += 1
            continue
        reconstructed_mtt.append(solution.m_ttbar)
        reconstructed_cosphi.append(solution.cos_phi)
        reconstructed_weights.append(float(event["mcWeight"]))
        qualities.append(solution.maximum_scaled_residual)

    mtt = np.asarray(reconstructed_mtt)
    cos_phi = np.asarray(reconstructed_cosphi)
    weights = np.asarray(reconstructed_weights)
    if mtt.size == 0:
        raise RuntimeError("the numerical reconstruction found no accepted solutions")

    regions = {
        "signal": (mtt > 340.0) & (mtt < 380.0),
        "validation_mid": (mtt >= 380.0) & (mtt < 500.0),
        "validation_high": mtt >= 500.0,
    }
    region_results: dict[str, object] = {}
    for name, mask in regions.items():
        if np.sum(mask) < 2:
            region_results[name] = {"events": int(np.sum(mask)), "d_proxy": None}
            continue
        region_results[name] = {
            "events": int(np.sum(mask)),
            "d_detector_proxy": estimate_d(cos_phi[mask], weights[mask]),
            "statistical_uncertainty_proxy": estimate_d_uncertainty(
                cos_phi[mask], weights[mask]
            ),
        }

    summary = {
        "input_file": str(args.root_file),
        "selected_events_attempted": len(chosen_events),
        "solutions": int(mtt.size),
        "failures": failures,
        "solution_efficiency": float(mtt.size / len(chosen_events)),
        "median_maximum_scaled_residual": float(np.median(qualities)),
        "regions": region_results,
        "status": (
            "Development proxy from a numerical on-shell constraint solver. "
            "It is not the ATLAS Ellipse implementation, has no background subtraction, "
            "and has no detector-to-particle calibration."
        ),
    }
    (output_dir / "reconstruction_proxy_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), constrained_layout=True)
    axes[0].hist(mtt, bins=np.linspace(300, 700, 41), weights=weights, histtype="step", linewidth=1.8)
    axes[0].axvspan(340, 380, color="#f59e0b", alpha=0.2, label="ATLAS signal window")
    axes[0].set_xlabel("Reconstructed m(ttbar) [GeV]")
    axes[0].set_ylabel("Generator-weighted events")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.22)
    signal = regions["signal"]
    axes[1].hist(
        cos_phi[signal],
        bins=np.linspace(-1, 1, 17),
        weights=weights[signal],
        histtype="step",
        linewidth=1.8,
        color="#1d4ed8",
    )
    axes[1].set_xlabel("Reconstructed cos(phi), signal window")
    axes[1].set_ylabel("Generator-weighted events")
    axes[1].grid(alpha=0.22)
    fig.suptitle("Dileptonic ttbar numerical reconstruction prototype")
    fig.savefig(figures_dir / "reconstruction_proxy.png", dpi=180)
    fig.savefig(figures_dir / "reconstruction_proxy.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
