"""Measure bounded shear removal using up to eight binary queries to the audit."""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np

from basinmark.collections import prepared, rate
from basinmark.io import JsonStore
from basinmark.reliability import reliable_audit
from basinmark.symmetry_audit import eligible_recovery


def read_inputs(root):
    cohort = root / "CrystalFormer-confirm.json.gz"
    _, rows = prepared(cohort)
    measurements = JsonStore(root / "CF-confirm-scores.json.gz").read_object()
    sheared = JsonStore(root / "CF-shear-scored.json.gz").read_object()["rows"]
    trials = JsonStore(root / "CF-confirm-ranked.json.gz").read_object()
    if hashlib.sha256(cohort.read_bytes()).hexdigest() != measurements["input_sha256"]:
        raise ValueError("ordinary observations do not match the input cohort")
    physical = []
    for r, view, shear in zip(rows, measurements["rows"], sheared, strict=True):
        if r["id"] != view["id"] or r["id"] != shear["id"] or r["producer"] != shear["producer"]:
            raise ValueError("unpaired physical records")
        ordinary = None
        if "rerelaxed" in r["stages"]:
            ordinary = (
                r["stages"]["rerelaxed"]["no_cell"],
                eligible_recovery(
                    view["stages"]["rerelaxed"][measurements["tolerances"].index(0.30)]
                ),
                view["scores"]["rerelaxed"],
            )
        alternate = None
        if (
            ordinary
            and r.get("valid_after")
            and shear["accepted"]
            and shear.get("valid_after")
            and abs(shear["delta_eV_per_atom"]) <= 0.010
            and "shear_rerelaxed" in shear["stages"]
        ):
            alternate = (
                shear["stages"]["shear_rerelaxed"]["no_cell"],
                eligible_recovery(shear["views"]["shear_rerelaxed"]),
                shear["scores"]["shear_rerelaxed"],
            )
            if abs(alternate[2] - ordinary[2]) > 0.010:
                alternate = None
        physical.append((r["id"], ordinary, alternate))
    return physical, trials


def proposed_sets(physical, indices):
    eligible = [i for i in indices if physical[i][2]]
    # Proposal construction sees only measurable structures and computed energies.
    proposals = [[], eligible]
    for field in (0, 1):
        proposals.append([i for i in eligible if physical[i][1][field] != physical[i][2][field]])
    for salt in range(4):
        proposals.append(
            [
                i
                for i in eligible
                if hashlib.sha256(f"{salt}|{physical[i][1][0]}".encode()).digest()[0] < 128
            ]
        )
    return proposals


def compare(root, phonon_file=None, rule=None):
    physical, trials = read_inputs(root)
    if bool(phonon_file) != bool(rule) or rule not in (None, "both-pass", "preserve"):
        raise ValueError("supply a phonon file and a valid screening rule together")
    if phonon_file:
        screen = JsonStore(phonon_file).read_object()
        cutoff = screen["phonon_pass_minimum_THz"]
        expected = {
            i for t in trials["records"] for i in t["selected"] if physical[i][2] is not None
        }
        if expected != {r["index"] for r in screen["records"]} or len(expected) != len(
            screen["records"]
        ):
            raise ValueError("phonon coverage differs from eligible emitted structures")
        allowed = set()
        for item in screen["records"]:
            i = item["index"]
            if item["id"] != physical[i][0]:
                raise ValueError("phonon identity differs")
            old = item["ordinary"]["minimum_THz"] >= cutoff
            new = item["edited"]["minimum_THz"] >= cutoff
            permissible = old and new if rule == "both-pass" else (not old or new)
            if permissible:
                allowed.add(i)
        physical = [
            (rid, obs, alt if i in allowed else None) for i, (rid, obs, alt) in enumerate(physical)
        ]
    counts = Counter()
    records = []
    for trial in trials["records"]:
        key = np.random.default_rng(trial["seed"]).bytes(16)
        indices = trial["selected"]
        if len(indices) != 32 or len(set(indices)) != 32:
            raise ValueError("invalid emission set")
        proposals = proposed_sets(physical, indices)
        outcomes = []
        for proposal in proposals:
            edits = set(proposal)
            obs = [physical[i][2 if i in edits else 1] for i in indices]
            result = reliable_audit([row for row in obs if row is not None], key)
            outcomes.append(int(result.fires))
        if trial["stages"]["rerelaxed"][-1] and not outcomes[0]:
            raise ValueError("ordinary baseline lost a protected detection")
        chosen = outcomes.index(0) if 0 in outcomes else 0
        counts.update(
            ordinary=outcomes[0],
            fixed_bounded_shear=outcomes[1],
            adaptive=min(outcomes),
            escaped=outcomes[0] - min(outcomes),
            gained=max(0, min(outcomes) - outcomes[0]),
        )
        records.append(
            {
                "seed": trial["seed"],
                "selected": indices,
                "decisions": outcomes,
                "chosen": chosen,
                "edited": proposals[chosen],
                "queries": chosen + 1 if 0 in outcomes else len(proposals),
            }
        )
    return {
        "source": "CrystalFormer first confirmation",
        "protocol": "bounded-shear-eight-queries",
        **(
            {"phonon_sha256": hashlib.sha256(phonon_file.read_bytes()).hexdigest(), "rule": rule}
            if phonon_file
            else {}
        ),
        "input_sha256": {
            f: hashlib.sha256((root / f).read_bytes()).hexdigest()
            for f in (
                "CrystalFormer-confirm.json.gz",
                "CF-confirm-scores.json.gz",
                "CF-shear-scored.json.gz",
                "CF-confirm-ranked.json.gz",
            )
        },
        "emissions": 32,
        "candidates": 256,
        "trials": len(records),
        "eligible_structures": sum(alt is not None for _, _, alt in physical),
        "sum_queries": sum(r["queries"] for r in records),
        "conditional_removal": rate(counts["escaped"], counts["ordinary"]),
        "summary": dict(counts),
        "records": records,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("results/v0.5"))
    ap.add_argument("--phonons", type=Path)
    ap.add_argument("--rule", choices=["both-pass", "preserve"])
    dest = ap.add_mutually_exclusive_group(required=True)
    dest.add_argument("--out", type=Path)
    dest.add_argument("--verify", type=Path)
    a = ap.parse_args()
    if a.out and a.out.exists():
        ap.error("supply a new output")
    result = compare(a.root, a.phonons, a.rule)
    if a.verify:
        if JsonStore(a.verify).read_object() != result:
            raise ValueError("removal replay differs")
    else:
        JsonStore(a.out).write(result)
    print(result["summary"])
