"""Trace the published strain loss and compare recovery settings on the design cohort."""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np

from basinmark.collections import prepared
from basinmark.io import JsonStore
from basinmark.symmetry_audit import eligible_recovery, paired_audit


def diagnose(corpus, audits, views):
    _, rows = prepared(corpus)
    data = JsonStore(audits).read_object()
    extra = JsonStore(views).read_object()
    if extra["input_sha256"] != hashlib.sha256(corpus.read_bytes()).hexdigest() or [
        r["id"] for r in extra["rows"]
    ] != [r["id"] for r in rows]:
        raise ValueError("recovery views do not match the source")
    transitions, changes, paired = Counter(), Counter(), Counter()
    scans = {str(t): {"all": 0, "omit_triclinic": 0} for t in extra["tolerances"]}
    stage_sums = {}
    for trial, seed in zip(data["records"], data["trial_seeds"], strict=True):
        key = np.random.default_rng(seed).bytes(16)
        baseline = trial["budget625"]
        for stage in ["rerelaxed", "strain_rerelaxed"]:
            v = baseline["stages"][stage]
            stage_sums.setdefault(stage, Counter()).update(
                deposited=v[0], classes=v[1], green=v[2], detections=v[3]
            )
        paired[str([baseline["stages"][s][3] for s in ["rerelaxed", "strain_rerelaxed"]])] += 1
        primary, alternatives = [], []
        for i in baseline["selected"]:
            stages = rows[i]["stages"]
            a, b = stages.get("rerelaxed"), stages.get("strain_rerelaxed")
            if a is None or b is None:
                changes["missing_either"] += 1
            else:
                left, right = a["no_cell"].split("|"), b["no_cell"].split("|")
                fields = [
                    f
                    for f, x, y in zip(("group", "formula", "sites"), left, right, strict=True)
                    if x != y
                ]
                changes["+".join(fields) if fields else "identical"] += 1
                if fields:
                    transitions[left[0] + " to " + right[0]] += 1
            if b is not None:
                primary.append(b["no_cell"])
                alternatives.append(extra["rows"][i]["stages"]["strain_rerelaxed"])
        for j, t in enumerate(extra["tolerances"]):
            alt = [v[j] for v in alternatives if v[j] is not None]
            for name, secondary in [
                ("all", alt),
                ("omit_triclinic", [fp for fp in alt if eligible_recovery(fp)]),
            ]:
                result = paired_audit(primary, secondary, key)
                if result.primary_fires != bool(baseline["stages"]["strain_rerelaxed"][3]):
                    raise ValueError("original decision differs")
                scans[str(t)][name] += result.fires
    return {
        "trials": data["trials"],
        "stage_sums": stage_sums,
        "paired": dict(paired),
        "selected_changes": dict(changes),
        "group_transitions": dict(transitions),
        "recovery_sweep": scans,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("corpus", "audits", "views", "out"):
        ap.add_argument("--" + name, required=True, type=Path)
    a = ap.parse_args()
    if a.out.exists():
        ap.error("choose a new output")
    JsonStore(a.out).write(diagnose(a.corpus, a.audits, a.views))


if __name__ == "__main__":
    main()
