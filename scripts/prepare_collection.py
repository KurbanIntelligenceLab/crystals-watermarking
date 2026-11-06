"""Prepare keyed-selection inputs from the WyFormer benchmark's public CSV archives."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import re
import time
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import distribution, version
from pathlib import Path

import numpy as np
from pymatgen.core import Composition, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from basinmark.core import fingerprint_structure
from basinmark.io import JsonStore
from basinmark.relaxation import CalculatorRegistry, RelaxationEngine

# File identifiers and SHA-256 digests for Figshare 29145101, version 6.
SOURCES = {
    "MatterGen": (
        "id",
        65498775,
        "a0a4ebfa441425eb72d4c39d6f1c950ef574f2751b26e5b37683aca3562bb69e",
    ),
    "CrystalFormer": (
        "material_id",
        65498610,
        "56bbe1b7737e5ab1454f923e42fe42db587419cdd99694df8df10838fd66dea4",
    ),
    "WyCryst": (
        "material_id",
        65498871,
        "6cfb544307a77f360d0e4a5e06ea00cc55b873b78b5a2783c086cea7f4729a12",
    ),
}
DFT = (65498808, "a25ef2c678cb2b24eff79e3b7ee1a814c20b0d05eaee22110d26f2d87718a25f")
REGISTRY = CalculatorRegistry()


def read_rows(directory, file_id, digest):
    path = directory / f"{file_id}.csv.gz"
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Unexpected SHA-256 for {path.name}")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def joined_candidates(dataset, rows):
    """Join the energy table to the geometry ID and verify each retained composition."""
    geometry_id = SOURCES[dataset][0]
    energy_id = {"MatterGen": "id", "CrystalFormer": "str_ind", "WyCryst": "folder_ind"}[dataset]
    energies = {}
    stats = Counter()
    for row in rows:
        rid = row[energy_id]
        if rid in energies:
            raise ValueError(f"Duplicate energy ID {rid}")
        value = float(row["chgnet_e_above_hull_corrected"] or "nan")
        if dataset == "MatterGen":
            formula = row["filename"].split("/")[0].split("_")[1]
        else:
            formula = row["formula"]
            match = re.search(r"\(e_h:\s*([^()\s]+)\)", row["error_str"])
            if match:
                token = match.group(1)
                precise = float(token)
                if not math.isclose(precise, value, abs_tol=1e-10):
                    mantissa, sep, _ = token.lower().partition("e")
                    if not sep or not math.isclose(float(mantissa), value, abs_tol=1e-10):
                        raise ValueError(f"Unresolved energy inconsistency at {rid}")
                    value = precise
                    stats["restored_scientific_exponents"] += 1
        energies[rid] = (value, Composition(formula).fractional_composition)
    candidates = []
    ids = set()
    for row in rows:
        rid = row[geometry_id]
        if rid in ids:
            raise ValueError(f"Duplicate geometry ID {rid}")
        ids.add(rid)
        if rid not in energies:
            stats["missing_energy_id"] += 1
            continue
        energy, expected = energies[rid]
        structure = Structure.from_dict(json.loads(row["structure"]))
        if not structure.composition.fractional_composition.almost_equals(expected):
            raise ValueError(f"Composition mismatch after energy ID join at {rid}")
        if math.isfinite(energy):
            candidates.append((rid, energy, row["structure"]))
        else:
            stats["missing_finite_energy"] += 1
    return candidates, ids, {"energy_id": energy_id, "geometry_id": geometry_id, **stats}


def variants(structure):
    full = fingerprint_structure(structure)
    # The final five fields are the lattice ratios and angles in this encoding.
    return {"full": full, "no_cell": full.rsplit("|", 5)[0]}


def deposit(structure, digits):
    primitive = SpacegroupAnalyzer(structure, symprec=0.1).get_primitive_standard_structure()
    cif = str(CifWriter(primitive, significant_figures=digits))
    return Structure.from_str(cif, fmt="cif", primitive=False)


def process(item):
    rid, energy, encoded, endpoint, digits, engine = item
    out = {"id": rid, "energy": energy, "stages": {}, "errors": {}}
    try:
        producer = Structure.from_dict(json.loads(encoded))
        out["producer"] = variants(producer)
        out["formula"] = producer.composition.reduced_formula
        out["valid"] = bool(producer.is_valid())
    except Exception as exc:
        out["errors"]["producer"] = type(exc).__name__
        return out
    try:
        out["stages"]["cif"] = variants(deposit(producer, digits))
    except Exception as exc:
        out["errors"]["cif"] = type(exc).__name__
    if endpoint:
        try:
            ref = Structure.from_dict(json.loads(endpoint[0]))
            if ref.composition.reduced_formula != out["formula"]:
                raise ValueError("Reference composition differs")
            out["stages"]["paired_dft"] = variants(deposit(ref, digits))
            out["dft_hull"] = endpoint[1]
        except Exception as exc:
            out["errors"]["paired_dft"] = type(exc).__name__
    if engine:
        import torch

        torch.set_num_threads(1)
        start = time.perf_counter()
        try:
            run = RelaxationEngine(REGISTRY, engine, fmax=0.05, steps=150).run(producer)
            out["steps"] = run.steps
            out["converged"] = run.converged
            if run.converged:
                out["stages"]["rerelaxed"] = variants(deposit(run.structure, digits))
                out["valid_after"] = bool(run.structure.is_valid())
            else:
                out["errors"]["rerelaxed"] = "step_limit"
        except Exception as exc:
            out["errors"]["rerelaxed"] = type(exc).__name__
        out["relax_seconds"] = time.perf_counter() - start
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=SOURCES, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--exclude-cohort", type=Path, action="append", default=[],
                        help="reserve previously sampled IDs; repeat for multiple cohorts")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--digits", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--paired-dft", action="store_true")
    parser.add_argument("--engine", choices=["chgnet", "chgnet-r2scan"])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.workers < 1 or args.digits < 1:
        parser.error("choose a new output and positive worker/precision settings")
    if args.paired_dft and args.dataset != "MatterGen":
        parser.error("the paired DFT mapping is defined for MatterGen")
    _, file_id, digest = SOURCES[args.dataset]
    refs = {}
    if args.paired_dft:
        for row in read_rows(args.directory, *DFT):
            rid = re.fullmatch(r"gen_(\d+)\.cif", row["material_id"]).group(1)
            if rid in refs:
                raise ValueError(f"Duplicate DFT ID {rid}")
            score = float(row["dft_e_above_hull_corrected"] or "nan")
            refs[rid] = (row["structure"], score if math.isfinite(score) else None)
    candidates, ids, join = joined_candidates(
        args.dataset, list(read_rows(args.directory, file_id, digest))
    )
    items = [
        (rid, energy, encoded, refs.get(rid), args.digits, args.engine)
        for rid, energy, encoded in candidates
    ]
    finite_energies = len(items)
    excluded = set()
    for excluded_file in args.exclude_cohort:
        previous = JsonStore(excluded_file).read_object()
        if previous["metadata"]["sha256"] != digest:
            parser.error("excluded cohort belongs to another source")
        excluded.update(r["id"] for r in previous["rows"])
        if not excluded <= ids:
            parser.error("excluded cohort contains unknown IDs")
    items = [row for row in items if row[0] not in excluded]
    usable = len(items)
    if args.limit is not None:
        if not 0 < args.limit <= usable:
            parser.error("--limit must be between 1 and the finite-energy population size")
        indices = np.random.default_rng(args.seed).choice(usable, args.limit, replace=False)
        items = [items[i] for i in indices]
    modules = ["numpy", "pymatgen", "spglib", "ase"]
    weights = None
    if args.engine:
        modules += ["torch", "chgnet"]
        model = "r2scan" if args.engine == "chgnet-r2scan" else "0.3.0"
        folder = distribution("chgnet").locate_file(f"chgnet/pretrained/{model}")
        (checkpoint,) = folder.glob("*.pth.tar")
        weights = {
            "name": checkpoint.name,
            "sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        }
    metadata = {
        "dataset": args.dataset,
        "file_id": file_id,
        "sha256": digest,
        "source_records": len(ids),
        "energy_join": join,
        "finite_energy_records": finite_energies,
        "excluded_cohort_size": len(excluded),
        "excluded_ids_sha256": hashlib.sha256("\n".join(sorted(excluded)).encode()).hexdigest(),
        "seed": args.seed,
        "sample_size": len(items),
        "digits": args.digits,
        "paired_dft_source": DFT if args.paired_dft else None,
        "engine": args.engine,
        "weights": weights,
        "fmax": 0.05 if args.engine else None,
        "step_limit": 150 if args.engine else None,
        "workers": args.workers,
        "versions": {name: version(name) for name in modules},
        "python": platform.python_version(),
        "machine": platform.machine(),
    }
    warnings.filterwarnings("ignore", category=UserWarning)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(process, items, chunksize=2))
    failures = Counter(k for row in rows for k in row["errors"])
    JsonStore(args.out).write({"metadata": metadata, "errors_by_stage": failures, "rows": rows})


if __name__ == "__main__":
    main()
