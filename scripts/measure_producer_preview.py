"""Measure an independent producer-side stress test on clean input geometries."""

import argparse
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from measure_symmetry_views import structures
from prepare_collection import REGISTRY, deposit
from prepare_phononbench import _r2scan_weights
from pymatgen.core import Lattice, Structure

from basinmark.collections import prepared
from basinmark.core import fingerprint_structure
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine
from basinmark.symmetry_audit import eligible_recovery


def preview(job):
    import torch

    rid, encoded = job
    torch.set_num_threads(1)
    structure = Structure.from_dict(json.loads(encoded))
    # Preview changes a copy. Each selected output is still the original candidate.
    matrix = np.diag([1.008, 1 / 1.008, 1.0]) @ structure.lattice.matrix
    start = Structure(Lattice(matrix), structure.species, structure.frac_coords)
    tic = time.perf_counter()
    try:
        result = RelaxationEngine(REGISTRY, "chgnet-r2scan", fmax=0.05, steps=150).run(start)
        if not result.converged:
            return {"id": rid, "status": "unconverged", "seconds": time.perf_counter() - tic}
        endpoint = deposit(result.structure, 5)
        primary = fingerprint_structure(endpoint, use_cell_shape=False)
        recovery = eligible_recovery(
            fingerprint_structure(endpoint, symprec=0.3, use_cell_shape=False)
        )
        return {
            "id": rid,
            "primary": primary,
            "secondary": recovery,
            "status": "ok",
            "seconds": time.perf_counter() - tic,
        }
    except (RuntimeError, ArithmeticError, ValueError, TypeError) as exc:
        return {"id": rid, "status": type(exc).__name__, "seconds": time.perf_counter() - tic}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("corpus", "directory", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    a = parser.parse_args()
    if a.out.exists() or a.workers < 1:
        parser.error("supply a new destination and positive worker count")
    meta, rows = prepared(a.corpus)
    inputs, _ = structures(a.directory, a.archive, meta, rows)
    jobs = [(r["id"], inputs[r["id"]]) for r in rows]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        observations = list(pool.map(preview, jobs, chunksize=2))
    JsonStore(a.out).write(
        {
            "source_sha256": hashlib.sha256(a.corpus.read_bytes()).hexdigest(),
            "weights_sha256": hashlib.sha256(_r2scan_weights().read_bytes()).hexdigest(),
            "strain": 0.008,
            "fmax": 0.05,
            "steps": 150,
            "rows": observations,
        }
    )
