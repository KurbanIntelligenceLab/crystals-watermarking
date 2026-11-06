"""Replay collection-level selection, material quality and downstream audits."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from basinmark.allocation import _select
from basinmark.collections import prepared, rate, tally
from basinmark.core import is_green
from basinmark.io import JsonStore

POLICIES = ("single", "minimum", "pooled_minimum", "unbounded", "fixed", "reserve", "pool")


def trial(rows, indices, key, audit_key, stages, followup):
    candidates = [rows[i] for i in indices]
    energy = np.array([r["energy"] for r in candidates])
    green = np.array([is_green(r["producer"]["no_cell"], key) for r in candidates])
    choices = _select(energy, green, 8, 0.1, POLICIES)
    green625 = None
    if followup:
        choices.update(_select(energy, green, 8, 0.1, ["budget"]))
        green625 = np.array([is_green(r["producer"]["no_cell"], key, 0.625) for r in candidates])
        choices.update(
            {
                p + "625": s
                for p, s in _select(energy, green625, 8, 0.1, ["fixed", "pool", "budget"]).items()
            }
        )
    results = {}
    for name, choice in choices.items():
        gamma = 0.625 if name.endswith("625") else 0.5
        labels = green625 if gamma == 0.625 else green
        ix = choice.indices
        selected = [candidates[i] for i in ix]
        diff = energy[ix] - choice.reference_energies
        result = {
            "selected": indices[ix].tolist(),
            "mean_energy": float(energy[ix].mean()),
            "mean_excess": float(diff.mean()),
            "max_excess": float(diff.max()),
            "unmarked_emissions": int((~labels[ix]).sum()),
            "valid": sum(r["valid"] for r in selected),
            "distinct_formulas": len({r["formula"] for r in selected}),
            "stages": {},
        }
        if choice.batch_sizes is not None:
            result["batch_sizes"] = choice.batch_sizes.tolist()
        for stage in stages:
            outputs = [
                r["producer"] if stage == "producer" else r["stages"].get(stage) for r in selected
            ]
            fps = [o["no_cell"] for o in outputs if o is not None]
            k, g, fired = tally(fps, key, gamma, 0.001)
            _, other_g, other_fired = tally(fps, audit_key, gamma, 0.001)
            result["stages"][stage] = [len(fps), k, g, fired, other_g, other_fired]
        for metric in ("dft_hull", "phonon_stable"):
            values = [r[metric] for r in selected if r.get(metric) is not None]
            if values:
                result[metric] = [
                    len(values),
                    sum(v <= 0.1 if metric == "dft_hull" else int(v) for v in values),
                ]
        results[name] = result
    return results


def summarize(records, stages):
    n = len(records)
    summary = {}
    for policy in records[0]:
        samples = [r[policy] for r in records]
        item = {
            f"mean_{k}": float(np.mean([r[k] for r in samples]))
            for k in ("unmarked_emissions", "valid", "distinct_formulas")
        }
        for k in ("mean_energy", "mean_excess"):
            item[k + "_meV_per_atom"] = 1000 * float(np.mean([r[k] for r in samples]))
        item["max_excess_meV_per_atom"] = 1000 * max(r["max_excess"] for r in samples)
        item["stages"] = {}
        for stage in stages:
            deposited, distinct, _, hits, _, wrong = np.array(
                [r["stages"][stage] for r in samples]
            ).T
            baseline = np.array([r["fixed"]["stages"][stage][3] for r in records])
            item["stages"][stage] = {
                **rate(int(sum(hits)), n),
                "wrong_key": rate(int(sum(wrong)), n),
                "mean_K": float(distinct.mean()),
                "coverage": float(deposited.mean() / 32),
                "gained_vs_fixed": int(np.sum((hits == 1) & (baseline == 0))),
                "lost_vs_fixed": int(np.sum((hits == 0) & (baseline == 1))),
            }
        for metric in ("dft_hull", "phonon_stable"):
            measured, stable = np.sum([r.get(metric, [0, 0]) for r in samples], axis=0)
            if measured:
                item[metric] = {"coverage": measured / (n * 32), "stable_yield": stable / (n * 32)}
        summary[policy] = item
    return summary


def benchmark(corpus, seed, trials, followup=False, drop_lowest=False):
    meta, rows = prepared(corpus)
    omitted = None
    if drop_lowest:
        lowest = min(rows, key=lambda row: row["energy"])
        omitted = lowest["id"]
        rows = [r for r in rows if r["id"] != omitted]
    if len(rows) < 256 or trials < 1:
        raise ValueError("expected at least 256 eligible records and a positive trial count")
    stages = ["producer", "cif", *meta.get("additional_stages", [])]
    if meta["engine"]:
        stages.append("rerelaxed")
    if meta["paired_dft_source"]:
        stages.append("paired_dft")
    seeds = np.random.default_rng(seed).integers(0, 2**53, size=trials).tolist()
    records = []
    for s in seeds:
        rng = np.random.default_rng(s)
        key, other = rng.bytes(16), rng.bytes(16)
        indices = rng.permutation(len(rows))[:256]
        records.append(trial(rows, indices, key, other, stages, followup))
    return {
        **({"gamma625": 0.625, "followup": True} if followup else {}),
        **({"excluded_lowest_id": omitted} if omitted is not None else {}),
        "dataset": meta["dataset"],
        "seed": seed,
        "trials": trials,
        "m": 8,
        "emissions": 32,
        "gamma": 0.5,
        "alpha": 0.001,
        "energy_window_eV_per_atom": 0.1,
        "eligible_records": len(rows),
        "eligible_ids_sha256": hashlib.sha256(
            json.dumps([r["id"] for r in rows]).encode()
        ).hexdigest(),
        "stage_columns": ["deposited", "K", "G", "reject", "wrong_G", "wrong_reject"],
        "quality_columns": ["assessed", "stable"],
        "trial_seeds": seeds,
        "summary": summarize(records, stages),
        "records": records,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--corpus", type=Path, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--trials", type=int, required=True)
    p.add_argument("--followup", action="store_true")
    p.add_argument("--drop-lowest", action="store_true")
    outputs = p.add_mutually_exclusive_group(required=True)
    outputs.add_argument("--out", type=Path)
    outputs.add_argument("--verify", type=Path)
    a = p.parse_args()
    if a.out and a.out.exists():
        p.error("choose a new output")
    result = benchmark(a.corpus, a.seed, a.trials, a.followup, a.drop_lowest)
    if a.verify:
        if result != JsonStore(a.verify).read_object():
            raise ValueError("recorded allocation results differ")
        print("Verified", a.verify.name)
    else:
        JsonStore(a.out).write(result)


if __name__ == "__main__":
    main()
