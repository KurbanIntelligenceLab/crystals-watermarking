"""Measure producer selection CPU time on exactly the same prepared candidates."""

import argparse
import platform
import time
from pathlib import Path

import numpy as np

from basinmark import select_collection
from basinmark.collections import prepared
from basinmark.io import JsonStore


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        ap.error("choose a new output")
    meta, rows = prepared(a.corpus)
    if len(rows) < 256:
        ap.error("expected 256 or more eligible candidates")
    rng = np.random.default_rng(3119)
    results = []
    for i in range(201):
        ix = rng.permutation(len(rows))[:256]
        key = rng.bytes(16)
        fps = [rows[j]["producer"]["no_cell"] for j in ix]
        energy = [rows[j]["energy"] for j in ix]
        times = {}
        for policy in rng.permutation(["fixed", "reserve", "pool", "budget625"]):
            start = time.process_time_ns()
            if policy == "budget625":
                select_collection(fps, energy, key, policy="budget", gamma=0.625)
            else:
                select_collection(fps, energy, key, policy=policy)
            times[str(policy)] = time.process_time_ns() - start
        if i:
            results.append(times)
    JsonStore(a.out).write(
        {
            "clock": "process_time_ns",
            "unit": "nanoseconds per 32 emissions",
            "dataset": meta["dataset"],
            "processor": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "seed": 3119,
            "median_ms": {p: float(np.median([r[p] for r in results]) / 1e6) for p in results[0]},
            "trials": results,
        }
    )


if __name__ == "__main__":
    main()
