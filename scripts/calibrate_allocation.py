"""Calibrate the deposited pool-selected collections with independent audit keys."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from basinmark.collections import cutoff, prepared, rate, tally
from basinmark.core import binom_sf
from basinmark.io import JsonStore


def evaluate(corpus, audits, seed, per_collection):
    _, rows = prepared(corpus)
    source = JsonStore(audits).read_object()

    digest = hashlib.sha256(json.dumps([r["id"] for r in rows]).encode()).hexdigest()
    if digest != source["eligible_ids_sha256"]:
        raise ValueError("audit indices refer to a different corpus")
    rng = np.random.default_rng(seed)
    records, hits, counts = [], Counter(), Counter()
    policies = ["fixed", "pool", "pooled_minimum"]
    if source.get("followup"):
        policies.extend(["fixed625", "budget625"])
    for trial in source["records"]:
        fps = []
        for policy in policies:
            selected = [rows[i]["stages"].get("rerelaxed") for i in trial[policy]["selected"]]
            fps.append(sorted({s["no_cell"] for s in selected if s is not None}))
        keys_seed = int(rng.integers(0, 2**53))
        keys = np.random.default_rng(keys_seed)
        ks = list(map(len, fps))
        observations = []
        for _ in range(per_collection):
            key = keys.bytes(16)
            g = []
            for policy, classes in zip(policies, fps, strict=True):
                gamma = 0.625 if policy.endswith("625") else 0.5
                k, green, rejected = tally(classes, key, gamma, 0.001)
                g.append(green)
                hits[policy] += rejected
                counts[policy, k] += 1
            observations.append(g)
        records.append([keys_seed, ks, observations])
    n = len(records) * per_collection
    summary = {}
    for policy in policies:
        gamma = 0.625 if policy.endswith("625") else 0.5
        attainable = (
            sum(
                count * binom_sf(cutoff(k, gamma, 0.001) - 1, k, gamma)
                for (p, k), count in counts.items()
                if p == policy
            )
            / n
        )
        summary[policy] = {**rate(hits[policy], n), "mean_exact_null_rate": attainable}
    return {
        "seed": seed,
        "keys_per_collection": per_collection,
        "policies": policies,
        "records": records,
        "columns": ["keys_seed", "K_by_policy", "G_by_key_and_policy"],
        "summary": summary,
        "eligible_ids_sha256": digest,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ("corpus", "audits", "out"):
        ap.add_argument("--" + arg, type=Path, required=True)
    ap.add_argument("--seed", type=int, default=54827)
    ap.add_argument("--keys-per-collection", type=int, default=20)
    a = ap.parse_args()
    if a.out.exists() or a.keys_per_collection <= 0:
        ap.error("choose a new output and positive keys per collection")
    JsonStore(a.out).write(evaluate(a.corpus, a.audits, a.seed, a.keys_per_collection))


if __name__ == "__main__":
    main()
