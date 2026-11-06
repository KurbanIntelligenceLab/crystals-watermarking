"""Estimate symmetry-class survival from historical paired structures, without audit keys."""

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from basinmark.io import JsonStore
from basinmark.symmetry_audit import eligible_recovery


def fit(root):
    prior = root / "v0.5"
    inputs = [
        (root / "v0.4/CrystalFormer-second.json.gz", "CF-dev"),
        (prior / "CrystalFormer-confirm.json.gz", "CF-confirm"),
        (prior / "CF-A.json.gz", "CF-A"),
        (prior / "CF-B.json.gz", "CF-B"),
    ]
    counts = defaultdict(Counter)
    for source, name in inputs:
        data = JsonStore(source).read_object()
        views = JsonStore(prior / (name + "-scores.json.gz")).read_object()
        if views["input_sha256"] != hashlib.sha256(source.read_bytes()).hexdigest():
            raise ValueError("historical measurement digest differs")
        col = views["tolerances"].index(0.30)
        rows = [r for r in data["rows"] if "producer" in r]
        for row, view in zip(rows, views["rows"], strict=True):
            if row["id"] != view["id"]:
                raise ValueError("historical source identities differ")
            stage = "strain_rerelaxed"
            if stage not in view["stages"]:
                continue
            recovered = eligible_recovery(view["stages"][stage][col])
            if recovered:
                sg = recovered.partition("|")[0]
                agree = int(row["stages"][stage]["no_cell"] == recovered)
                correct = int(row["producer"]["no_cell"] == recovered)
                counts[sg, agree][correct] += 1
    return {f"{sg}/{agree}": [c[0], c[1]] for (sg, agree), c in sorted(counts.items())}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("results"))
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    values = fit(args.root)
    if args.verify:
        from basinmark.reliability import survival_counts

        if values != survival_counts():
            raise ValueError("historical counts differ")
        print("Verified historical survival counts")
    else:
        print(json.dumps(values, indent=2))
