"""Replay fixed-budget paired audits on fresh or historical measured structures."""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
from benchmark_symmetry_audit import measure
from fit_survival import fit
from split_recovery_cohort import split

from basinmark.collections import prepared, rate
from basinmark.core import is_green
from basinmark.io import JsonStore
from basinmark.reliability import reliable_audit, survival_counts
from basinmark.symmetry_audit import eligible_recovery

NEW = Path("results/v0.6")
PRIOR = Path("results/v0.5")
CASES = {
    "A": (PRIOR, "CF-A", "CF-A-scores", "CF-A-audits"),
    "B": (PRIOR, "CF-B", "CF-B-scores", "CF-B-audits"),
    "first-confirm": (PRIOR, "CrystalFormer-confirm", "CF-confirm-scores", "CF-confirm-ranked"),
    "dev": (Path("results/v0.4"), "CrystalFormer-second", "CF-dev-scores", "CF-dev-ranked"),
    "shear": (PRIOR, "CF-shear-scored", None, "CF-shear-ranked"),
    "mattergen": (PRIOR, "MatterGen-confirm", "MG-confirm-scores", "MG-confirm-ranked"),
}
FRESH = {"C": (271573, 413833), "D": (862317, 739021)}


def path(base, stem):
    return base / (stem + ".json.gz")


def check_sources(directory):
    parent = path(directory, "CF-reserved")
    prior_sources = [
        path(Path("results/v0.3/measurements"), "CrystalFormer-corpus"),
        path(Path("results/v0.4"), "CrystalFormer-external"),
        path(Path("results/v0.4"), "CrystalFormer-second"),
        path(PRIOR, "CrystalFormer-confirm"),
        path(PRIOR, "CF-reserved"),
    ]
    ids, classes = set(), set()
    for index, source in enumerate(
        prior_sources + [path(directory, "CF-C"), path(directory, "CF-D")]
    ):
        rows = JsonStore(source).read_object()["rows"]
        found_ids = {r["id"] for r in rows}
        found_classes = {r["producer"]["no_cell"] for r in rows if "producer" in r}
        if index >= len(prior_sources) and (found_ids & ids or found_classes & classes):
            raise ValueError("fresh source overlaps a prior ID or class")
        ids.update(found_ids)
        classes.update(found_classes)
    for cohort, part in zip(FRESH, split(parent), strict=True):
        if JsonStore(path(directory, "CF-" + cohort)).read_object() != part:
            raise ValueError("reserved split differs")
    if fit(Path("results")) != survival_counts():
        raise ValueError("frozen training counts differ")


def benchmark(case, directory):
    if case in FRESH:
        source = path(directory, "CF-" + case)
        views = path(directory, "CF-" + case + "-scores")
        seed, null_seed = FRESH[case]
        reference = None
    else:
        base, corpus, scores, old = CASES[case]
        source = path(base, corpus)
        views = path(PRIOR, scores) if scores else None
        reference = JsonStore(path(PRIOR, old)).read_object()
        seed, null_seed = reference["seed"], reference["control_seed"]
    base = measure(source, views, 0.30, seed, 500, null_seed, 0)
    _, rows = prepared(source)
    if views:
        raw = JsonStore(views).read_object()
        col = raw["tolerances"].index(0.30)
        stages = [({s: f[col] for s, f in v["stages"].items()}, v["scores"]) for v in raw["rows"]]
    else:
        stages = [(r["views"], r["scores"]) for r in rows]
    observed = []
    for r, (secondary, scores) in zip(rows, stages, strict=True):
        deposited = {"producer": r["producer"], **r["stages"]}
        observed.append(
            {
                s: (fp["no_cell"], eligible_recovery(secondary[s]), scores[s])
                for s, fp in deposited.items()
            }
        )
    counts = {s: Counter() for s in base["summary"]}
    diagnosis = {outcome: Counter() for outcome in ("miss", "hit")}
    records = []
    for i, trial in enumerate(base["records"]):
        if reference is not None and trial != reference["records"][i]:
            raise ValueError("legacy replay differs in selection, material outcomes or audit")
        selected = trial["selected"]
        key = np.random.default_rng(trial["seed"]).bytes(16)
        rng = np.random.default_rng(trial["control_seed"])
        wrong = [rng.bytes(16) for _ in range(20)]
        outcome = {**trial, "stages": {}}
        for stage, previous in trial["stages"].items():
            indices = [j for j in selected if stage in observed[j]]
            obs = [observed[j][stage] for j in indices]
            decision = reliable_audit(obs, key)
            if decision.ranked_fires != bool(previous[-1]) or (previous[-1] and not decision.fires):
                raise ValueError("legacy decision changed")
            if decision.null_rate > 0.001 + 1e-14:
                raise ValueError("conditional null probability exceeds the audit level")
            guards = [reliable_audit(obs, k) for k in wrong]
            false = sum(a.fires for a in guards)
            old_false = sum(a.ranked_fires for a in guards)
            counts[stage].update(
                old=decision.ranked_fires,
                new=decision.fires,
                gained=decision.fires and not decision.ranked_fires,
                lost=decision.ranked_fires and not decision.fires,
                null=false,
                old_null=old_false,
                size=decision.null_rate,
            )
            outcome["stages"][stage] = [
                len(indices),
                int(decision.ranked_fires),
                decision.reliable_classes,
                decision.reliable_green,
                int(decision.fires),
                old_false,
                false,
                decision.null_rate,
            ]
            if stage == "strain_rerelaxed" and case in ("A", "B"):
                group = diagnosis["hit" if decision.ranked_fires else "miss"]
                group.update(
                    trials=1,
                    producer_detected=trial["stages"]["producer"][-1],
                    ordinary_detected=trial["stages"]["rerelaxed"][-1],
                    missing=32 - len(indices),
                    producer_non_green=sum(
                        not is_green(rows[j]["producer"]["no_cell"], key, 0.625) for j in selected
                    ),
                )
                ordered = sorted(
                    indices,
                    key=lambda j: (
                        observed[j][stage][2],
                        observed[j][stage][0],
                        observed[j][stage][1] or "",
                        j,
                    ),
                )
                for rank, j in enumerate(ordered):
                    p, s, _ = observed[j][stage]
                    original = rows[j]["producer"]["no_cell"]
                    group["matched_primary"] += p == original
                    group["matched_secondary"] += s == original
                    group["matched_both"] += p == s == original
                    group["neither"] += original not in (p, s)
                    quartile = rank * 4 // len(ordered)
                    group[f"q{quartile}_matches"] += original in (p, s)
                    group[f"q{quartile}_total"] += 1
        records.append(outcome)
    return {
        "cohort": case,
        "source_sha256": base["corpus_sha256"],
        "views_sha256": base["views_sha256"],
        "seed": seed,
        "control_seed": null_seed,
        "outputs": 32,
        "candidates": 256,
        "level": 0.001,
        "controls_per_trial": 20,
        "columns": [
            "deposited",
            "ranked",
            "reliable_classes",
            "reliable_green",
            "new",
            "wrong_ranked",
            "wrong_new",
            "null_probability",
        ],
        "summary": {
            s: {
                "ranked": rate(c["old"], 500),
                "reliable": rate(c["new"], 500),
                "gained": c["gained"],
                "lost": c["lost"],
                "wrong_ranked": c["old_null"],
                "wrong_reliable": c["null"],
                "mean_null_probability": c["size"] / 500,
            }
            for s, c in counts.items()
        },
        **({"diagnosis": diagnosis} if case in ("A", "B") else {}),
        "records": records,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", choices=[*CASES, *FRESH], required=True)
    ap.add_argument("--directory", type=Path, default=NEW)
    ap.add_argument("--write", action="store_true", help="write one new archive, otherwise verify")
    a = ap.parse_args()
    target = path(a.directory, a.case + "-reliability")
    if a.write and target.exists():
        ap.error("choose a directory without an existing output")
    if a.case in FRESH:
        check_sources(a.directory)
    result = benchmark(a.case, a.directory)
    if a.write:
        JsonStore(target).write(result)
    elif result != JsonStore(target).read_object():
        raise ValueError("archived results differ")
    print("Verified", a.case, hashlib.sha256(target.read_bytes()).hexdigest(), flush=True)
