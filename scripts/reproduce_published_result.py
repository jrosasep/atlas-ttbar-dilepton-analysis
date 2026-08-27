"""Reproduce the published-number arithmetic and make summary figures."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src"))

from ttbar_entanglement.published import (  # noqa: E402
    gaussian_distance_from_boundary,
    load_config,
    measurement_from_dict,
    systematic_quadrature,
)
from ttbar_entanglement.toy import angular_pdf  # noqa: E402


def main() -> None:
    config = load_config(PROJECT / "config" / "published_atlas.yaml")
    results_dir = PROJECT / "results"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    labels: list[str] = []
    observed = []
    observed_error = []
    expected = []
    expected_error = []
    output_regions: dict[str, object] = {}

    for name, region in config["regions"].items():
        obs = measurement_from_dict(region["observed"])
        exp = measurement_from_dict(region["expected"])
        labels.append(region["label"])
        observed.append(obs.value)
        observed_error.append(obs.total)
        expected.append(exp.value)
        expected_error.append(exp.total)
        output_regions[name] = {
            "label": region["label"],
            "observed": {"value": obs.value, "total_uncertainty": obs.total},
            "expected": {"value": exp.value, "total_uncertainty": exp.total},
        }

    signal_observed = measurement_from_dict(
        config["regions"]["signal"]["observed"]
    )
    signal_expected = measurement_from_dict(
        config["regions"]["signal"]["expected"]
    )
    boundary = config["particle_level_boundaries"]["powheg_pythia"]
    boundary_value = float(boundary["value"])
    boundary_uncertainty = float(boundary["uncertainty"])

    summary = {
        "source": config["publication"],
        "regions": output_regions,
        "signal_systematic_quadrature_unrounded": systematic_quadrature(config),
        "quoted_signal_systematic": signal_observed.syst,
        "powheg_pythia_particle_boundary": boundary,
        "gaussian_cross_check_not_official_significance": {
            "observed_sigma": gaussian_distance_from_boundary(
                signal_observed, boundary_value, boundary_uncertainty
            ),
            "expected_sigma": gaussian_distance_from_boundary(
                signal_expected, boundary_value, boundary_uncertainty
            ),
            "note": (
                "Independent Gaussian quadrature only; the ATLAS paper reports "
                "observed and expected significances well beyond five sigma."
            ),
        },
    }
    (results_dir / "published_result_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    x_positions = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9.0, 5.3), constrained_layout=True)
    ax.errorbar(
        x_positions - 0.08,
        observed,
        yerr=observed_error,
        fmt="o",
        color="#111827",
        capsize=4,
        label="ATLAS observed",
    )
    ax.errorbar(
        x_positions + 0.08,
        expected,
        yerr=expected_error,
        fmt="s",
        color="#c2410c",
        capsize=4,
        label="Powheg+Pythia expected",
    )
    ax.fill_between(
        [-0.42, 0.42],
        boundary_value - boundary_uncertainty,
        boundary_value + boundary_uncertainty,
        color="#2563eb",
        alpha=0.18,
        label="Folded no-entanglement boundary (signal region)",
    )
    ax.hlines(boundary_value, -0.42, 0.42, color="#2563eb", linestyles="--")
    ax.axvline(0.5, color="#d1d5db", linewidth=1)
    plot_labels = [
        r"$340 < m_{t\bar{t}} < 380$ GeV",
        r"$380 < m_{t\bar{t}} < 500$ GeV",
        r"$m_{t\bar{t}} > 500$ GeV",
    ]
    ax.set_xticks(x_positions, plot_labels)
    ax.set_ylabel("Entanglement marker D")
    ax.set_title(r"ATLAS Run-2 $t\bar{t}$ entanglement: published particle-level values")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(figures_dir / "published_regions.png", dpi=180)
    fig.savefig(figures_dir / "published_regions.pdf")
    plt.close(fig)

    grid = np.linspace(-1.0, 1.0, 500)
    fig, ax = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    for d_value, label, color in [
        (-0.537, "Observed signal-region D", "#111827"),
        (-0.470, "SM expected D", "#c2410c"),
        (-1.0 / 3.0, "Parton separability boundary", "#2563eb"),
    ]:
        ax.plot(grid, angular_pdf(grid, d_value), label=label, color=color, linewidth=2)
    ax.set_xlabel("cos(phi)")
    ax.set_ylabel("Normalized density")
    ax.set_title("Angular model underlying the D moment estimator")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(figures_dir / "angular_pdf.png", dpi=180)
    fig.savefig(figures_dir / "angular_pdf.pdf")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
