"""Paired finite-displacement phonon screen on bounded-removal endpoints."""

import argparse
import hashlib
import zipfile
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import numpy as np
from ase import Atoms
from benchmark_bounded_removal import read_inputs
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from prepare_collection import REGISTRY, deposit
from prepare_phononbench import ARCHIVE, PREFIX, _r2scan_weights
from pymatgen.core import Lattice, Structure

from basinmark.core import fingerprint_structure
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine
from basinmark.symmetry_audit import recovery_class


def phonons(structure):
    cell = PhonopyAtoms(
        symbols=[str(el) for el in structure.species],
        cell=structure.lattice.matrix,
        scaled_positions=structure.frac_coords,
    )
    phonon = Phonopy(cell, np.diag([2, 2, 2]), primitive_matrix=np.eye(3), symprec=1e-5)
    phonon.generate_displacements(distance=0.01, is_plusminus=True)
    forces = []
    calculator = REGISTRY.get("chgnet-r2scan")
    for displaced in phonon.supercells_with_displacements:
        atoms = Atoms(
            symbols=displaced.symbols,
            cell=displaced.cell,
            scaled_positions=displaced.scaled_positions,
            pbc=True,
        )
        atoms.calc = calculator
        force = atoms.get_forces()
        forces.append(force - force.mean(axis=0))
    phonon.forces = forces
    phonon.produce_force_constants(fc_calculator="traditional")
    phonon.symmetrize_force_constants()
    grid = np.indices((4, 4, 4)).reshape(3, -1).T / 4
    phonon.run_qpoints(grid)
    frequencies = phonon.qpoints.frequencies
    if not np.isfinite(frequencies).all():
        raise ValueError("nonfinite phonon frequency")
    return {
        "minimum_THz": float(frequencies.min()),
        "gamma_minimum_THz": float(frequencies[0].min()),
        "supercell_atoms": len(phonon.supercell),
        "displacements": len(forces),
        "wavevectors": len(grid),
    }


def calculate(job):
    import torch

    torch.set_num_threads(1)
    index, rid, cif, ordinary, edited = job
    s = Structure.from_str(cif, fmt="cif")
    shear = np.eye(3)
    shear[0, 1] = 0.015
    alt = Structure(Lattice(shear @ s.lattice.matrix), s.species, s.frac_coords)
    row = {"index": index, "id": rid}
    for name, initial, expected in [("ordinary", s, ordinary), ("edited", alt, edited)]:
        result = RelaxationEngine(REGISTRY, "chgnet-r2scan", fmax=0.05, steps=150).run(initial)
        if not result.converged or not result.structure.is_valid():
            raise ValueError(f"invalid reconstruction for {rid}")
        endpoint = deposit(result.structure, 5)
        if (
            fingerprint_structure(endpoint, use_cell_shape=False),
            recovery_class(endpoint),
        ) != expected[:2]:
            raise ValueError(f"reconstructed audit views differ for {rid}")
        row[name] = phonons(endpoint)
    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("results/v0.5"))
    for arg in ["removal", "archive", "out"]:
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    a = parser.parse_args()
    if a.out.exists() or a.workers < 1:
        parser.error("choose a new output and positive worker count")
    if hashlib.sha256(a.archive.read_bytes()).hexdigest() != ARCHIVE:
        raise ValueError("unexpected source archive")
    physical, _ = read_inputs(a.root)
    trials = JsonStore(a.removal).read_object()
    indices = sorted(
        {i for row in trials["records"] for i in row["selected"] if physical[i][2] is not None}
    )
    with zipfile.ZipFile(a.archive) as archive:
        jobs = [
            (i, physical[i][0], archive.read(PREFIX + physical[i][0]).decode(), *physical[i][1:])
            for i in indices
        ]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        records = list(pool.map(calculate, jobs))
    JsonStore(a.out).write(
        {
            "eligible_unique_candidates": len(indices),
            "indices": indices,
            "removal_sha256": hashlib.sha256(a.removal.read_bytes()).hexdigest(),
            "source_sha256": ARCHIVE,
            "model_sha256": hashlib.sha256(_r2scan_weights().read_bytes()).hexdigest(),
            "versions": {
                k: version(k) for k in ("chgnet", "phonopy", "spglib", "torch", "ase", "pymatgen")
            },
            "force_model": "CHGNet-r2SCAN",
            "supercell": [2, 2, 2],
            "mesh": [4, 4, 4],
            "phonon_pass_minimum_THz": -0.3,
            "displacement_angstrom": 0.01,
            "records": records,
        }
    )
