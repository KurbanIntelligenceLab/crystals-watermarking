"""Replay paired production with and without a separate deformation preview."""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
from benchmark_reliability import CASES, FRESH, NEW, PRIOR, check_sources, path
from split_recovery_cohort import split

from basinmark.allocation import select_collection
from basinmark.collections import prepared, rate
from basinmark.io import JsonStore
from basinmark.reliability import reliable_audit
from basinmark.symmetry_audit import eligible_recovery

CONFIRM = {"E": (368219, 537619), "F": (653291, 938261)}


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def verify_confirmation(directory):
    check_sources(directory)
    seen_ids, seen_classes = set(), set()
    for file in (
        path(Path("results/v0.3/measurements"), "CrystalFormer-corpus"),
        path(Path("results/v0.4"), "CrystalFormer-external"),
        path(Path("results/v0.4"), "CrystalFormer-second"),
        path(PRIOR, "CrystalFormer-confirm"),
        path(PRIOR, "CF-reserved"),
        path(directory, "CF-reserved"),
    ):
        rows = JsonStore(file).read_object()["rows"]
        seen_ids.update(r["id"] for r in rows)
        seen_classes.update(r["producer"]["no_cell"] for r in rows if "producer" in r)
    for name, part in zip(CONFIRM, split(path(directory, "EF-reserved")), strict=True):
        if part != JsonStore(path(directory, "CF-" + name)).read_object():
            raise ValueError("prospective split differs")
        ids = {r["id"] for r in part["rows"]}
        classes = {r["producer"]["no_cell"] for r in part["rows"] if "producer" in r}
        if ids & seen_ids or classes & seen_classes:
            raise ValueError("prospective IDs or producer classes overlap")
        seen_ids.update(ids)
        seen_classes.update(classes)


def benchmark(case, directory):
    reference = None
    if case in CONFIRM or case in FRESH:
        source = path(directory, "CF-" + case)
        views = path(directory, "CF-" + case + "-scores")
        if case in CONFIRM:
            verify_confirmation(directory)
            seed, null_seed = CONFIRM[case]
        else:
            seed, null_seed = FRESH[case]
    else:
        base, stem, score, old = CASES[case]
        source = path(base, stem)
        views = path(PRIOR, score) if score else None
        reference = JsonStore(path(PRIOR, old)).read_object()
        seed, null_seed = reference["seed"], reference["control_seed"]
    _, rows = prepared(source)
    preview_file = path(directory, ("first-confirm" if case == "shear" else case) + "-preview")
    measured = JsonStore(preview_file).read_object()
    producer_source = path(PRIOR, "CrystalFormer-confirm") if case == "shear" else source
    if measured["source_sha256"] != sha(producer_source):
        raise ValueError("preview was measured on different inputs")
    raw = JsonStore(views).read_object() if views else None
    if raw and raw["input_sha256"] != sha(source):
        raise ValueError("downstream views were measured on different inputs")
    producer = []
    previews = []
    observations = []
    for i, (r, pv) in enumerate(zip(rows, measured["rows"], strict=True)):
        if r["id"] != pv["id"] or (raw and raw["rows"][i]["id"] != r["id"]):
            raise ValueError("physical records are out of order")
        producer.append(r["producer"]["no_cell"])
        previews.append((pv.get("primary"), pv.get("secondary")))
        if raw:
            item = raw["rows"][i]
            secondary = {s: v[raw["tolerances"].index(0.30)] for s, v in item["stages"].items()}
            scores = item["scores"]
        else:
            secondary, scores = r["views"], r["scores"]
        observations.append(
            {
                s: (val["no_cell"], eligible_recovery(secondary[s]), scores[s])
                for s, val in {"producer": r["producer"], **r["stages"]}.items()
            }
        )
    stages = sorted(set.union(*(set(r) for r in observations)))
    counts = {s: Counter() for s in stages}
    energies = np.array([r["energy"] for r in rows])
    seeds = np.random.default_rng(seed).integers(0, 2**53, 500).tolist()
    wrong_rng = np.random.default_rng(null_seed)
    records = []
    for t, s in enumerate(seeds):
        rng = np.random.default_rng(s)
        key = rng.bytes(16)
        rng.bytes(16)
        batch = rng.permutation(len(rows))[:256]
        fps = [producer[i] for i in batch]
        base = select_collection(fps, energies[batch], key)
        new = select_collection(
            fps, energies[batch], key, policy="preview", previews=[previews[i] for i in batch]
        )
        selections = [batch[x.indices].tolist() for x in (base, new)]
        if reference and selections[0] != reference["records"][t]["selected"]:
            raise ValueError("control release differs from archived reference")
        ctrl_seed = int(wrong_rng.integers(0, 2**53))
        wrng = np.random.default_rng(ctrl_seed)
        controls = [wrng.bytes(16) for _ in range(20)] if case in CONFIRM else []
        measures = []
        oracle_free_seconds = sum(measured["rows"][i]["seconds"] for i in batch)
        for ids in selections:
            known = [rows[i]["dft_hull"] for i in ids if rows[i].get("dft_hull") is not None]
            label = (
                sum(rows[i]["phonon_stable"] for i in ids)
                if "phonon_stable" in rows[ids[0]]
                else sum(x <= 0.1 for x in known)
            )
            measures.append(
                {
                    "mean_energy": float(energies[ids].mean()),
                    "valid": sum(rows[i]["valid"] for i in ids),
                    "formulas": len({rows[i]["formula"] for i in ids}),
                    "stable": label,
                }
            )
        row = {
            "seed": s,
            "control_seed": ctrl_seed,
            "selected": selections,
            "quality": measures,
            "preview_worker_seconds": oracle_free_seconds,
            "stages": {},
        }
        for stage in stages:
            decisions = []
            false = []
            for ids in selections:
                obs = [observations[i][stage] for i in ids if stage in observations[i]]
                decision = reliable_audit(obs, key)
                decisions.append(int(decision.fires))
                false.append(sum(reliable_audit(obs, k).fires for k in controls))
            a, b = decisions
            counts[stage].update(
                previous=a,
                preview=b,
                gained=(not a and b),
                lost=(a and not b),
                old_wrong=false[0],
                new_wrong=false[1],
            )
            row["stages"][stage] = [*decisions, *false]
        records.append(row)
    return {
        "cohort": case,
        "source_sha256": sha(source),
        "preview_sha256": sha(preview_file),
        "views_sha256": sha(views) if views else None,
        "seed": seed,
        "null_seed": null_seed,
        "controls": 20 if case in CONFIRM else 0,
        "trials": 500,
        "candidates": 256,
        "emissions": 32,
        "stage_columns": ["budget", "preview", "budget_wrong", "preview_wrong"],
        "summary": {
            stage: {**c, "preview_rate": rate(c["preview"], 500)} for stage, c in counts.items()
        },
        "records": records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=NEW)
    parser.add_argument("--case", choices=[*CASES, *FRESH, *CONFIRM], required=True)
    parser.add_argument("--write", action="store_true")
    a = parser.parse_args()
    target = path(a.directory, a.case + "-producer")
    if a.write and target.exists():
        parser.error("output already exists")
    result = benchmark(a.case, a.directory)
    if a.write:
        JsonStore(target).write(result)
    elif result != JsonStore(target).read_object():
        raise ValueError("producer replay differs")
    print(a.case, {s: r["preview"] for s, r in result["summary"].items()}, flush=True)
