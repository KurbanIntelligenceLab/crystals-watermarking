"""Profile symmetry recovery and optional deposited-energy ranking separately."""

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
import torch
from measure_symmetry_views import structures
from pymatgen.core import Structure

from basinmark.collections import prepared
from basinmark.core import fingerprint_structure
from basinmark.deposit_energy import DepositedEnergy
from basinmark.io import JsonStore
from basinmark.symmetry_audit import recovery_class


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ("corpus", "archive", "weights", "hull", "out"):
        ap.add_argument("--" + arg, required=True, type=Path)
    a = ap.parse_args()
    if a.out.exists():
        ap.error("choose a new output")
    torch.set_num_threads(1)
    meta, rows = prepared(a.corpus)
    sample = [rows[i] for i in np.random.default_rng(18761).permutation(len(rows))[:65]]
    inputs, _ = structures(None, a.archive, meta, sample)
    start = time.perf_counter()
    ranker = DepositedEnergy(a.weights, a.hull)
    timings = []
    for n, r in enumerate(sample):
        s = Structure.from_dict(json.loads(inputs[r["id"]]))
        t = time.perf_counter()
        p = fingerprint_structure(s, use_cell_shape=False)
        primary = time.perf_counter() - t
        if p != r["producer"]["no_cell"]:
            raise ValueError("wrong source fingerprint")
        t = time.perf_counter()
        recovery_class(s)
        geometric = time.perf_counter() - t
        t = time.perf_counter()
        score = ranker(s)
        rank = time.perf_counter() - t
        if not np.isclose(score, r["energy"], rtol=1e-9, atol=1e-8):
            raise ValueError("independent public ranking differs from source")
        if n:
            timings.append([r["id"], primary, geometric, rank])
        else:
            warmup = time.perf_counter() - start
    JsonStore(a.out).write(
        {
            "seed": 18761,
            "cpu": platform.machine(),
            "threads": 1,
            "warmup_seconds": warmup,
            "columns": ["id", "primary_seconds", "secondary_seconds", "rank_seconds"],
            "records": timings,
            "median_seconds": np.median([r[1:] for r in timings], axis=0).tolist(),
        }
    )


if __name__ == "__main__":
    main()
