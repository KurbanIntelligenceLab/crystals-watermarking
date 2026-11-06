"""Time both pipelines on the first registered E and F draws on one CPU thread."""

import argparse
import json
import platform
import time
from importlib.metadata import version
from pathlib import Path
from statistics import median

import numpy as np
import torch
from benchmark_producer_preview import NEW, path, sha
from compare_method_versions import ARMS, candidate_draw
from measure_symmetry_views import structures
from prepare_collection import REGISTRY, deposit
from prepare_phononbench import _r2scan_weights
from pymatgen.core import Lattice, Structure

from basinmark.allocation import select_collection
from basinmark.collections import prepared, tally
from basinmark.core import fingerprint_structure
from basinmark.deposit_energy import DepositedEnergy
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine
from basinmark.symmetry_audit import audit_deposited


def profile_case(case, args, engine, ranker):
    source = path(args.directory, "CF-" + case)
    meta, rows = prepared(source)
    if sha(_r2scan_weights()) != meta["r2scan_weights_sha256"]:
        raise ValueError("deposition checkpoint differs")
    comparison = JsonStore(path(args.comparisons, case + args.suffix)).read_object()
    if comparison["input_sha256"]["source"] != sha(source):
        raise ValueError("comparison does not use this cohort")
    trial = comparison["records"][0]
    key, batch = candidate_draw(trial["seed"], len(rows))
    inputs, _ = structures(None, args.archive, meta, [rows[i] for i in batch])
    geometry = {i: Structure.from_dict(json.loads(inputs[rows[i]["id"]])) for i in batch}
    energies = [rows[i]["energy"] for i in batch]
    preview = JsonStore(path(args.directory, case + "-preview")).read_object()
    views = [
        (preview["rows"][i].get("primary"), preview["rows"][i].get("secondary")) for i in batch
    ]
    # Generate exactly the strained deposits the frozen decisions used; stage costs are separate.
    unique = sorted({i for a in trial["arms"].values() for i in a["ids"]})
    endpoints, deposit_times = {}, {}
    for i in unique:
        start = geometry[i]
        if rows[i]["strain_accepted"]:
            deform = np.diag([1.01, 1 / 1.01, 1]) @ start.lattice.matrix
            start = Structure(Lattice(deform), start.species, start.frac_coords)
            start = Structure.from_dict(json.loads(start.to_json()))
        begin = time.process_time()
        result = engine.run(start)
        endpoint = deposit(result.structure, 5) if result.converged else None
        deposit_times[i] = time.process_time() - begin
        baseline = rows[i]["stages"].get("strain_rerelaxed")
        if (endpoint is None) != (baseline is None):
            raise ValueError(f"deposition convergence differs for {rows[i]['id']}")
        if endpoint is not None:
            if fingerprint_structure(endpoint) != baseline["full"]:
                raise ValueError(f"deposited fingerprint differs for {rows[i]['id']}")
            endpoints[i] = endpoint
    # Reread exactly the original no-cell outputs with the current auditor at their original gamma.
    arms = {arm: trial["arms"][arm] for arm in ARMS}
    arms["reader_only"] = arms["original_no_cell"]
    raw = {arm: [] for arm in arms}
    rng = np.random.default_rng(trial["seed"])
    for _ in range(3):
        for arm in rng.permutation(list(arms)).tolist():
            original = arm.startswith("original")
            full = arm == "original_full"
            begin = time.process_time()
            fps = [fingerprint_structure(geometry[i], use_cell_shape=full) for i in batch]
            choice = select_collection(
                fps,
                energies,
                key,
                policy="unbounded" if original or arm == "reader_only" else arm,
                previews=views if arm == "preview" else None,
            )
            selection_seconds = time.process_time() - begin
            ids = batch[choice.indices].tolist()
            if ids != arms[arm]["ids"]:
                raise ValueError("timed selection differs from paired experiment")
            deposits = [endpoints[i] for i in ids if i in endpoints]
            # Each audit starts with a preloaded model and a fresh composition cache.
            ranker.hulls.clear()
            begin = time.process_time()
            if original:
                fps = [fingerprint_structure(s, use_cell_shape=full) for s in deposits]
                fires = tally(fps, key, 0.5, 0.001)[2]
            else:
                gamma = 0.5 if arm == "reader_only" else 0.625
                fires = audit_deposited(deposits, key, ranker=ranker, gamma=gamma).fires
            audit_seconds = time.process_time() - begin
            column = 5 if arm == "reader_only" else 0
            if int(fires) != arms[arm]["stages"]["CHGNet_strain_rerelaxed"][column]:
                raise ValueError("timed decision differs from paired experiment")
            raw[arm].append([selection_seconds, audit_seconds])
    return {
        "seed": trial["seed"],
        "source_sha256": sha(source),
        "unique_attempted": len(unique),
        "unique_converged": len(endpoints),
        "timing_columns": ["selection_with_fingerprints", "deposited_audit"],
        "arms": {
            arm: {
                "trials": raw[arm],
                "median_selection": median(s for s, _ in raw[arm]),
                "median_audit": median(s for _, s in raw[arm]),
                "relax_and_deposit": sum(deposit_times[i] for i in a["ids"]),
                "available_outputs": sum(i in endpoints for i in a["ids"]),
            }
            for arm, a in arms.items()
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=NEW)
    parser.add_argument("--comparisons", type=Path, default=NEW)
    parser.add_argument("--suffix", default="-method-comparison")
    for flag in ("archive", "rank-weights", "hull", "out"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("choose a new output")
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    start = time.process_time()
    ranker = DepositedEnergy(args.rank_weights, args.hull)
    engine = RelaxationEngine(REGISTRY, "chgnet-r2scan", fmax=0.05, steps=150)
    meta, rows = prepared(path(args.directory, "CF-E"))
    warm, _ = structures(None, args.archive, meta, rows[:1])
    warm_structure = Structure.from_dict(json.loads(warm[rows[0]["id"]]))
    ranker(warm_structure)
    engine.run(warm_structure)
    profile = {
        "clock": "process CPU seconds",
        "machine": platform.machine(),
        "threads": 1,
        "versions": {
            p: version(p) for p in ("torch", "chgnet", "mace-torch", "numpy", "pymatgen", "ase")
        },
        "warmup_seconds": time.process_time() - start,
        "scope": "preloaded weights and hull; starting candidate structures and energies supplied",
        "cohorts": {},
    }
    for case in ("E", "F"):
        profile["cohorts"][case] = profile_case(case, args, engine, ranker)
        print(case, profile["cohorts"][case]["arms"], flush=True)
    JsonStore(args.out).write(profile)
