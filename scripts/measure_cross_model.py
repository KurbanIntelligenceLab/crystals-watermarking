"""Re-relax the fixed E or F emissions with MACE-MPA-0 and measure deposited views."""

import argparse
import warnings
import zipfile
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from importlib.metadata import version
from pathlib import Path

import numpy as np
from benchmark_cross_model import MPA
from benchmark_producer_preview import NEW, path, sha
from prepare_collection import deposit
from prepare_phononbench import ARCHIVE, PREFIX
from pymatgen.core import Lattice, Structure

from basinmark.collections import prepared
from basinmark.core import fingerprint_structure
from basinmark.deposit_energy import WEIGHTS, energy_per_atom, verify
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine
from basinmark.symmetry_audit import recovery_class


class ModelFileRegistry:
    """Load the supplied relaxation checkpoint once in each CPU worker."""

    @staticmethod
    @lru_cache(maxsize=1)
    def get(filename):
        from mace.calculators import MACECalculator

        return MACECalculator(model_paths=filename, device="cpu", default_dtype="float64")


def measure(job):
    import torch

    index, row, cif, relaxation_weights, ranking_weights = job
    torch.set_num_threads(1)
    source = Structure.from_str(cif, fmt="cif")
    if fingerprint_structure(source, use_cell_shape=False) != row["producer"]["no_cell"]:
        raise ValueError(f"producer differs for {row['id']}")
    engine = RelaxationEngine(ModelFileRegistry(), str(relaxation_weights), fmax=0.05, steps=150)
    strained = Structure(
        Lattice(np.diag([1.01, 1 / 1.01, 1.0]) @ source.lattice.matrix),
        source.species,
        source.frac_coords,
    )
    starts = {
        "rerelaxed": source,
        "strain_rerelaxed": strained if row["strain_accepted"] else source,
    }
    stages = {}
    for stage, start in starts.items():
        result = engine.run(start)
        state = {"steps": result.steps, "converged": result.converged}
        if result.converged:
            endpoint = deposit(result.structure, 5)
            if (
                endpoint.composition.fractional_composition
                != source.composition.fractional_composition
            ):
                raise ValueError(f"composition differs for {row['id']}")
            score = energy_per_atom(endpoint, ranking_weights) - row["hull_eV_per_atom"]
            if not np.isfinite(score):
                raise ValueError(f"nonfinite deposited energy for {row['id']}")
            state.update(
                primary=fingerprint_structure(endpoint, use_cell_shape=False),
                secondary=recovery_class(endpoint),
                energy=score,
                valid=bool(endpoint.is_valid()),
            )
        stages[stage] = state
    return {
        "index": index,
        "id": row["id"],
        "strain_accepted": row["strain_accepted"],
        "stages": stages,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["E", "F"], required=True)
    parser.add_argument("--directory", type=Path, default=NEW)
    for flag in ("archive", "relax-weights", "rank-weights", "out"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("supply a positive number of workers")
    verify(args.relax_weights, MPA)
    verify(args.rank_weights, WEIGHTS)
    verify(args.archive, ARCHIVE)
    source = path(args.directory, "CF-" + args.case)
    producer = path(args.directory, args.case + "-producer")
    _, rows = prepared(source)
    trials = JsonStore(producer).read_object()
    if trials["source_sha256"] != sha(source):
        raise ValueError("release file and source differ")
    ids = sorted({i for t in trials["records"] for arm in t["selected"] for i in arm})
    metadata = {
        "case": args.case,
        "archive_sha256": ARCHIVE,
        "relax_weights_sha256": MPA,
        "rank_weights_sha256": WEIGHTS,
        "source_sha256": sha(source),
        "producer_sha256": sha(producer),
        "versions": {
            p: version(p) for p in ("torch", "ase", "mace-torch", "pymatgen", "spglib", "numpy")
        },
        "fmax": 0.05,
        "steps": 150,
        "indices": ids,
    }
    store = JsonStore(args.out)
    data = store.read_object() if args.out.exists() else {"metadata": metadata, "records": []}
    if data["metadata"] != metadata:
        raise ValueError("checkpoint inputs or environment differ")
    done = [r["index"] for r in data["records"]]
    if done != ids[: len(done)] or any(
        r["id"] != rows[i]["id"] for i, r in zip(done, data["records"], strict=True)
    ):
        raise ValueError("checkpoint identities differ")
    with zipfile.ZipFile(args.archive) as z:
        jobs = [
            (
                i,
                rows[i],
                z.read(PREFIX + rows[i]["id"]).decode(),
                args.relax_weights,
                args.rank_weights,
            )
            for i in ids[len(done) :]
        ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for record in pool.map(measure, jobs):
            data["records"].append(record)
            store.write(data)
    print(f"Measured {len(data['records'])} unique structures on {args.case}")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("once")
        main()
