"""Inspect an ATLAS Open Data ROOT file without assuming a fixed schema."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import uproot


def categorize(branches: list[str]) -> dict[str, list[str]]:
    terms = {
        "top_or_truth": ("truth", "top", "parton", "parent", "pdgid"),
        "leptons": ("lep", "electron", "muon", "el_", "mu_"),
        "jets_and_btag": ("jet", "btag", "b_tag"),
        "missing_momentum": ("met", "missing"),
        "weights": ("weight", "xsec", "scale_factor"),
    }
    lowered = {branch: branch.lower() for branch in branches}
    return {
        category: sorted(
            branch
            for branch, lower in lowered.items()
            if any(term in lower for term in needles)
        )
        for category, needles in terms.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    with uproot.open(args.root_file) as root_file:
        objects = {key: root_file[key].classname for key in root_file.keys()}
        tree_keys = [key for key in root_file.keys() if hasattr(root_file[key], "keys")]
        trees = {}
        for key in tree_keys:
            obj = root_file[key]
            branches = list(obj.keys())
            trees[key] = {
                "class": obj.classname,
                "entries": int(obj.num_entries) if hasattr(obj, "num_entries") else None,
                "number_of_branches": len(branches),
                "categories": categorize(branches),
                "all_branches": branches,
            }

    report = {
        "file": str(args.root_file),
        "objects": objects,
        "trees": trees,
    }
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
