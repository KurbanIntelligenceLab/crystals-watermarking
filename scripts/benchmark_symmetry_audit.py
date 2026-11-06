"""Pair original and symmetry-recovery decisions on exactly the same emitted structures."""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np

from basinmark.allocation import select_collection
from basinmark.collections import prepared, rate
from basinmark.io import JsonStore
from basinmark.symmetry_audit import RECOVERY_SYMPREC, eligible_recovery, paired_audit, ranked_audit

COLUMNS = ["K_primary", "K_secondary", "K_shared", "G_primary", "G_secondary", "old", "new"]


def measure(corpus, view_file, tolerance, seed, trials, control_seed, controls):
    meta, rows = prepared(corpus)
    digest = hashlib.sha256(corpus.read_bytes()).hexdigest()
    if view_file is None:
        views = {
            "tolerances": [meta["embedded_recovery_tolerance"]],
            "input_sha256": digest,
            "rows": [
                {
                    "id": r["id"],
                    "stages": {s: [fp] for s, fp in r["views"].items()},
                    **({"scores": r["scores"]} if "scores" in r else {}),
                }
                for r in rows
            ],
        }
    else:
        views = JsonStore(view_file).read_object()
    if views["input_sha256"] != digest or [r["id"] for r in rows] != [
        v["id"] for v in views["rows"]
    ]:
        raise ValueError("view file does not correspond to these producer records")
    col = views["tolerances"].index(tolerance)
    stages = sorted({s for r in views["rows"] for s in r["stages"]})
    pairs = []
    for r, v in zip(rows, views["rows"], strict=True):
        original = {"producer": r["producer"], **r["stages"]}
        if original.keys() != v["stages"].keys():
            raise ValueError("deposit availability differs between views")
        pairs.append(
            {
                stage: (fp["no_cell"], eligible_recovery(v["stages"][stage][col]))
                for stage, fp in original.items()
            }
        )
    scored = "scores" in views["rows"][0]
    records = []
    total = {st: Counter() for st in stages}
    seeds = np.random.default_rng(seed).integers(0, 2**53, size=trials).tolist()
    null = np.random.default_rng(control_seed)
    for s in seeds:
        rng = np.random.default_rng(s)
        key = rng.bytes(16)
        rng.bytes(16)  # Match the secondary-key draw in the original paired benchmark.
        batch = rng.permutation(len(rows))[:256]
        energy = [rows[i]["energy"] for i in batch]
        chosen = select_collection([rows[i]["producer"]["no_cell"] for i in batch], energy, key)
        indices = batch[chosen.indices].tolist()
        outputs = [rows[i] for i in indices]
        row = {
            "seed": s,
            "selected": indices,
            "mean_energy": float(np.mean([r["energy"] for r in outputs])),
            "valid": sum(r["valid"] for r in outputs),
            "distinct_formulas": len({r["formula"] for r in outputs}),
            "stages": {},
        }
        for metric in ["dft_hull", "phonon_stable"]:
            known = [r[metric] for r in outputs if r.get(metric) is not None]
            if known:
                row[metric] = [
                    len(known),
                    sum(x <= 0.1 if metric == "dft_hull" else x for x in known),
                ]
        replicate_seed = int(null.integers(0, 2**53))
        control_rng = np.random.default_rng(replicate_seed)
        keys = [control_rng.bytes(16) for _ in range(controls)]
        row["control_seed"] = replicate_seed
        for st in stages:
            observations = [pairs[i][st] for i in indices if st in pairs[i]]
            primary = {x[0] for x in observations}
            secondary = {x[1] for x in observations if x[1] is not None}
            a = paired_audit(primary, secondary, key)
            wrong = [paired_audit(primary, secondary, k).fires for k in keys]
            row["stages"][st] = [
                a.primary_classes,
                a.secondary_classes,
                a.shared_classes,
                a.primary_green,
                a.secondary_green,
                int(a.primary_fires),
                int(a.fires),
            ]
            if scored:
                observed = [
                    (*pairs[i][st], views["rows"][i]["scores"][st])
                    for i in indices
                    if st in pairs[i]
                ]
                extended = ranked_audit(observed, key)
                if extended.geometric != a:
                    raise ValueError("rank calibration changed the geometric decision")
                wrong_ranked = sum(ranked_audit(observed, k).fires for k in keys)
                row["stages"][st].extend(
                    [extended.ranked_classes, extended.ranked_green, int(extended.fires)]
                )
                total[st].update(
                    ranked=extended.fires,
                    wrong_ranked=wrong_ranked,
                    ranked_mass=extended.null_rate,
                    lost_geometric=a.fires and not extended.fires,
                )
            total[st].update(
                old=a.primary_fires,
                new=a.fires,
                gained=a.fires and not a.primary_fires,
                lost=a.primary_fires and not a.fires,
                wrong=sum(wrong),
                null_mass=a.null_rate,
                deposited=len(observations),
            )
        records.append(row)
    summary = {}
    for st, c in total.items():
        summary[st] = {
            "old": rate(c["old"], trials),
            "new": rate(c["new"], trials),
            "gained": c["gained"],
            "lost": c["lost"],
            "wrong_key": rate(c["wrong"], controls * trials) if controls else None,
            "attainable_null": c["null_mass"] / trials,
            "coverage": c["deposited"] / (32 * trials),
        }
        if scored:
            summary[st].update(
                ranked=rate(c["ranked"], trials),
                lost_geometric=c["lost_geometric"],
                ranked_attainable=c["ranked_mass"] / trials,
                ranked_wrong_key=rate(c["wrong_ranked"], trials * controls) if controls else None,
            )
    return {
        "dataset": meta["dataset"],
        "corpus_sha256": digest,
        "views_sha256": hashlib.sha256(view_file.read_bytes()).hexdigest() if view_file else None,
        "tolerance": tolerance,
        "seed": seed,
        "trials": trials,
        "control_seed": control_seed,
        "controls": controls,
        "columns": COLUMNS + (["K_ranked", "G_ranked", "ranked"] if scored else []),
        "summary": summary,
        "records": records,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--views", type=Path, help="separate view file for a standard cohort")
    ap.add_argument("--tolerance", type=float, choices=[RECOVERY_SYMPREC], default=RECOVERY_SYMPREC)
    ap.add_argument("--seed", type=int, default=734201)
    ap.add_argument("--trials", type=int, default=500)
    ap.add_argument("--control-seed", type=int, default=486209)
    ap.add_argument("--controls", type=int, default=20)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--out", type=Path)
    group.add_argument("--verify", type=Path)
    a = ap.parse_args()
    if (a.out and a.out.exists()) or a.trials < 1 or a.controls < 0:
        ap.error("supply a new output and valid trial counts")
    result = measure(a.corpus, a.views, a.tolerance, a.seed, a.trials, a.control_seed, a.controls)
    if a.verify:
        if result != JsonStore(a.verify).read_object():
            raise ValueError("replayed results differ")
        print("Verified", a.verify.name)
    else:
        JsonStore(a.out).write(result)


if __name__ == "__main__":
    main()
