"""Profile producer relaxation and energy-bounded strain on public MatterGen records."""

import argparse
import hashlib
import json
import math
import platform
import time
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import distribution, version
from pathlib import Path

import numpy as np
from prepare_collection import REGISTRY, SOURCES, deposit, process, read_rows, variants
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor

from basinmark.core import select_green
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine

RAW = (65498772, "c0e033305802d56c1b4e90a826cefa5b1c169e2f05e7fad68fea9b6fe8961b0c")


def process_record(item):
    import torch

    torch.set_num_threads(1)
    rid, energy, encoded, raw_data = item
    record = process((rid, energy, encoded, None, 5, None))
    source = Structure.from_dict(json.loads(encoded))
    calc = REGISTRY.get("chgnet")
    if raw_data is not None:
        raw = Structure.from_dict(json.loads(raw_data))
        if source.composition.reduced_formula != raw.composition.reduced_formula:
            raise ValueError(f"Composition mismatch for {rid}")
        engine = RelaxationEngine(REGISTRY, "chgnet", fmax=0.05, steps=150)
        start = time.perf_counter()
        try:
            run = engine.run(raw)
            record["profile_steps"] = run.steps
            record["profile_converged"] = run.converged
        except Exception as exc:
            record["errors"]["raw_profile"] = type(exc).__name__
        record["profile_seconds"] = time.perf_counter() - start
    try:
        lattice = np.diag([1.01, 1 / 1.01, 1.0]) @ source.lattice.matrix
        strained = Structure(Lattice(lattice), source.species, source.frac_coords)

        def energy_per_atom(s):
            atoms = AseAtomsAdaptor.get_atoms(s)
            atoms.calc = calc
            return float(atoms.get_potential_energy() / len(atoms))

        delta = energy_per_atom(strained) - energy_per_atom(source)
        accepted = abs(delta) <= 0.010 and bool(strained.is_valid())
        record["strain_delta_eV_per_atom"] = delta
        record["strain_accepted"] = accepted
        record["stages"]["bounded_strain"] = (
            variants(deposit(strained, 5)) if accepted else record["stages"]["cif"]
        )
    except Exception as exc:
        record["errors"]["bounded_strain"] = type(exc).__name__
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--seed", type=int, default=4291)
    args = parser.parse_args()
    if args.out.exists() or args.workers < 1:
        parser.error("choose a new output and positive worker count")
    key, file_id, digest = SOURCES["MatterGen"]
    source = []
    for row in read_rows(args.directory, file_id, digest):
        score = float(row["chgnet_e_above_hull_corrected"] or "nan")
        if math.isfinite(score):
            source.append((row[key], score, row["structure"]))
    indices = np.random.default_rng(args.seed).choice(len(source), 256, replace=False)
    sample = [source[i] for i in indices]
    raw_ids = {f"gen_{rid}.cif" for rid, _, _ in sample[:128]}
    raw = {
        row["material_id"]: row["structure"]
        for row in read_rows(args.directory, *RAW)
        if row["material_id"] in raw_ids
    }
    if raw.keys() != raw_ids:
        raise ValueError("Missing raw timing inputs")
    items = [(rid, e, s, raw.get(f"gen_{rid}.cif")) for rid, e, s in sample]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(process_record, items, chunksize=2))
    timed = [r for r in rows if "profile_seconds" in r]
    times = np.array([r["profile_seconds"] for r in timed]).reshape(-1, 8)
    # Time one marked decision per batch including symmetry analysis on all eight inputs.
    # Structure parsing and model initialization are outside both timing comparisons.
    geometries = [Structure.from_dict(json.loads(encoded)) for _, _, encoded in sample]
    audit_key = np.random.default_rng(args.seed).bytes(16)
    costs = []
    for i in range(0, len(sample), 8):
        batch_energies = [r[1] for r in sample[i : i + 8]]
        start = time.perf_counter()
        fps = [variants(s)["no_cell"] for s in geometries[i : i + 8]]
        select_green(fps, batch_energies, audit_key, energy_window=0.1)
        costs.append(time.perf_counter() - start)
    folder = distribution("chgnet").locate_file("chgnet/pretrained/0.3.0")
    (checkpoint,) = folder.glob("*.pth.tar")
    meta = {
        "dataset": "MatterGen",
        "raw_source": RAW,
        "producer_source": [file_id, digest],
        "seed": args.seed,
        "engine": None,
        "paired_dft_source": None,
        "additional_stages": ["bounded_strain"],
        "strain_factors": [1.01, 1 / 1.01, 1],
        "max_absolute_energy_change_eV_per_atom": 0.010,
        "digits": 5,
        "weights_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        "workers": args.workers,
        "cpu_architecture": platform.machine(),
        "force_threshold": 0.05,
        "step_limit": 150,
        "energy_field": "chgnet_e_above_hull_corrected",
        "versions": {n: version(n) for n in ["chgnet", "torch", "ase", "spglib", "pymatgen"]},
    }
    JsonStore(args.out).write(
        {
            "metadata": meta,
            "rows": rows,
            "profile": {
                "timed_raw_structures": len(timed),
                "converged": sum(r.get("profile_converged", False) for r in timed),
                "mean_fingerprint_and_mark_seconds_per_batch": float(np.mean(costs)),
                "mean_worker_seconds_single": float(times[:, 0].mean()),
                "mean_worker_seconds_best_of_eight": float(times.sum(axis=1).mean()),
                "mean_worker_seconds_marked_of_eight_relaxation": float(times.sum(axis=1).mean()),
                "strain_accepted": sum(r.get("strain_accepted", False) for r in rows),
            },
        }
    )


if __name__ == "__main__":
    main()
