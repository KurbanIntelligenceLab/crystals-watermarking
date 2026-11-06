"""Recompute independent symmetry views of the same archived deposited geometries."""

import argparse
import hashlib
import json
import re
import zipfile
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import numpy as np
from prepare_collection import DFT, REGISTRY, SOURCES, deposit, read_rows
from pymatgen.core import Lattice, Structure

from basinmark.collections import prepared
from basinmark.core import fingerprint_structure
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine

PHONON_PREFIX = "PhononBench/CrystalFormer-gen/mp20/relaxed/"


def structures(directory, archive, meta, rows):
    ids = {r["id"] for r in rows}
    if meta["dataset"].startswith("PhononBench"):
        if hashlib.sha256(archive.read_bytes()).hexdigest() != meta["source_sha256"]:
            raise ValueError("unexpected PhononBench archive")
        with zipfile.ZipFile(archive) as z:
            return {
                i: Structure.from_str(z.read(PHONON_PREFIX + i).decode(), fmt="cif").to_json()
                for i in sorted(ids)
            }, {}
    field, fid, digest = SOURCES[meta["dataset"]]
    inputs = {
        r[field]: r["structure"] for r in read_rows(directory, fid, digest) if r[field] in ids
    }
    refs = {}
    if meta["paired_dft_source"]:
        for r in read_rows(directory, *DFT):
            i = re.fullmatch(r"gen_(\d+)\.cif", r["material_id"]).group(1)
            if i in ids:
                refs[i] = r["structure"]
    return inputs, refs


def measure(item):
    import torch

    row, encoded, dft, tolerances, weights = item
    torch.set_num_threads(1)
    source = Structure.from_dict(json.loads(encoded))
    endpoints = {"producer": source}
    baseline = {"producer": row["producer"], **row["stages"]}
    if "cif" in baseline:
        endpoints["cif"] = deposit(source, 5)
    if "paired_dft" in baseline:
        endpoints["paired_dft"] = deposit(Structure.from_dict(json.loads(dft)), 5)

    def relaxed(start):
        result = RelaxationEngine(REGISTRY, "chgnet-r2scan", fmax=0.05, steps=150).run(start)
        if not result.converged:
            raise ValueError("previously converged endpoint did not converge")
        return deposit(result.structure, 5)

    if "rerelaxed" in baseline:
        endpoints["rerelaxed"] = relaxed(source)
    if "strain_rerelaxed" in baseline or "bounded_strain" in baseline:
        strained = Structure(
            Lattice(np.diag([1.01, 1 / 1.01, 1.0]) @ source.lattice.matrix),
            source.species,
            source.frac_coords,
        )
        if "strain_rerelaxed" in baseline:
            endpoints["strain_rerelaxed"] = (
                relaxed(Structure.from_dict(json.loads(strained.to_json())))
                if row["strain_accepted"]
                else endpoints["rerelaxed"]
            )
        if "bounded_strain" in baseline:
            endpoints["bounded_strain"] = deposit(strained if row["strain_accepted"] else source, 5)
    if set(endpoints) != set(baseline):
        raise ValueError("unknown downstream stage")
    stages, scores = {}, {}
    for name, structure in endpoints.items():
        original = fingerprint_structure(structure, use_cell_shape=False)
        if original != baseline[name]["no_cell"]:
            raise ValueError(f"deposited class differs for {row['id']} at {name}")
        values = []
        for tolerance in tolerances:
            try:
                fp = fingerprint_structure(structure, symprec=tolerance, use_cell_shape=False)
            except (ValueError, RuntimeError, TypeError):
                fp = None
            values.append(fp)
        stages[name] = values
        if weights:
            from prepare_phononbench import calculator
            from pymatgen.io.ase import AseAtomsAdaptor

            if (
                structure.composition.fractional_composition
                != source.composition.fractional_composition
            ):
                raise ValueError("deposition changed the composition used for the hull reference")
            if row.get("hull_eV_per_atom") is None:
                scores[name] = None
                continue
            atoms = AseAtomsAdaptor.get_atoms(structure)
            atoms.calc = calculator(weights)
            score = float(atoms.get_potential_energy() / len(atoms) - row["hull_eV_per_atom"])
            if not np.isfinite(score):
                raise ValueError("nonfinite deposited energy")
            scores[name] = score
    return {"id": row["id"], "stages": stages, **({"scores": scores} if weights else {})}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("corpus", "directory", "out"):
        ap.add_argument("--" + name, type=Path, required=True)
    ap.add_argument("--archive", type=Path)
    ap.add_argument("--weights", type=Path, help="also measure deposited MACE-OMAT hull energies")
    ap.add_argument("--hull", type=Path, help="pinned hull to score a different producer source")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--tolerances", nargs="+", type=float, required=True)
    a = ap.parse_args()
    if a.out.exists() or a.workers < 1 or any(t <= 0 for t in a.tolerances):
        ap.error("choose a new output and positive parameters")
    meta, rows = prepared(a.corpus)
    if a.weights:
        from prepare_phononbench import WEIGHTS, verify

        verify(a.weights, WEIGHTS)
        from basinmark.deposit_energy import HULL

        if a.hull:
            verify(a.hull, HULL)
        elif meta.get("hull_sha256") != HULL:
            raise ValueError("supply the pinned hull for this producer source")
    inputs, refs = structures(a.directory, a.archive, meta, rows)
    if a.weights and a.hull:
        from basinmark.deposit_energy import hull_at, references

        hull = references(a.hull)
        for r in rows:
            comp = Structure.from_dict(json.loads(inputs[r["id"]])).composition
            try:
                r["hull_eV_per_atom"] = hull_at(comp, hull, "linear")
            except ValueError:
                r["hull_eV_per_atom"] = None
    jobs = [(r, inputs[r["id"]], refs.get(r["id"]), a.tolerances, a.weights) for r in rows]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        views = list(pool.map(measure, jobs, chunksize=2))
    JsonStore(a.out).write(
        {
            "input_sha256": hashlib.sha256(a.corpus.read_bytes()).hexdigest(),
            "dataset": meta["dataset"],
            "tolerances": a.tolerances,
            "rows": views,
            **(
                {"weights_sha256": WEIGHTS, "hull_sha256": HULL} if a.weights else {}
            ),
            "versions": {n: version(n) for n in ("pymatgen", "spglib", "chgnet", "ase", "torch")},
        }
    )


if __name__ == "__main__":
    main()
