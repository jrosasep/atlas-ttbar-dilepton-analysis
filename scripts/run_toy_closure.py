"""Run deterministic pseudo-experiment closure and reweighting checks."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src"))

from ttbar_entanglement.observable import (  # noqa: E402
    analytic_d_uncertainty,
    estimate_d,
)
from ttbar_entanglement.toy import (  # noqa: E402
    closure_trials,
    reweight_d,
    sample_cos_phi,
)


def main() -> None:
    seed = 20260827
    n_events = 20_000
    n_trials = 600
    truth_values = [-0.60, -0.537, -0.470, -1.0 / 3.0, -0.10, 0.0]
    results_dir = PROJECT / "results"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, float | int]] = []
    signal_estimates: np.ndarray | None = None
    for index, truth in enumerate(truth_values):
        estimates = closure_trials(truth, n_events, n_trials, seed + index)
        if truth == -0.537:
            signal_estimates = estimates
        rows.append(
            {
                "truth_d": truth,
                "mean_estimated_d": float(np.mean(estimates)),
                "bias": float(np.mean(estimates) - truth),
                "empirical_std": float(np.std(estimates, ddof=1)),
                "analytic_std": analytic_d_uncertainty(truth, n_events),
                "n_events": n_events,
                "n_trials": n_trials,
            }
        )

    with (results_dir / "toy_closure.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    rng = np.random.default_rng(seed + 100)
    nominal_sample = sample_cos_phi(-0.470, 2_000_000, rng)
    weights = reweight_d(nominal_sample, -0.470, -0.537)
    reweighted_estimate = estimate_d(nominal_sample, weights)
    reweighting_summary = {
        "nominal_truth": -0.470,
        "target_truth": -0.537,
        "unweighted_estimate": estimate_d(nominal_sample),
        "reweighted_estimate": reweighted_estimate,
        "absolute_closure_error": abs(reweighted_estimate - (-0.537)),
        "seed": seed + 100,
        "n_events": nominal_sample.size,
    }
    (results_dir / "reweighting_closure.json").write_text(
        json.dumps(reweighting_summary, indent=2) + "\n", encoding="utf-8"
    )

    truth = np.array([float(row["truth_d"]) for row in rows])
    means = np.array([float(row["mean_estimated_d"]) for row in rows])
    errors = np.array([float(row["empirical_std"]) for row in rows])
    fig, ax = plt.subplots(figsize=(6.2, 5.6), constrained_layout=True)
    ax.errorbar(truth, means, yerr=errors, fmt="o", capsize=4, color="#111827")
    span = np.linspace(min(truth) - 0.03, max(truth) + 0.03, 100)
    ax.plot(span, span, "--", color="#c2410c", label="Ideal closure")
    ax.set_xlabel("Generated D")
    ax.set_ylabel("Mean estimated D")
    ax.set_title(f"Moment-estimator closure ({n_trials} trials, N={n_events:,})")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(figures_dir / "toy_closure.png", dpi=180)
    fig.savefig(figures_dir / "toy_closure.pdf")
    plt.close(fig)

    assert signal_estimates is not None
    expected_std = analytic_d_uncertainty(-0.537, n_events)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.hist(signal_estimates, bins=30, density=True, color="#93c5fd", edgecolor="#1e3a8a")
    ax.axvline(-0.537, color="#c2410c", linewidth=2, label="Generated D = -0.537")
    ax.set_xlabel("Estimated D")
    ax.set_ylabel("Pseudo-experiment density")
    ax.set_title(f"Signal-point toy distribution; analytic sigma={expected_std:.4f}")
    ax.legend(frameon=False)
    fig.savefig(figures_dir / "toy_signal_distribution.png", dpi=180)
    fig.savefig(figures_dir / "toy_signal_distribution.pdf")
    plt.close(fig)

    print(json.dumps({"closure": rows, "reweighting": reweighting_summary}, indent=2))


if __name__ == "__main__":
    main()
