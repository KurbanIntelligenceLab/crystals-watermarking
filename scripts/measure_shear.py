"""Apply a bounded simple shear and independent relaxation to a PhononBench cohort."""

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from measure_symmetry_views import structures
from prepare_collection import REGISTRY, deposit, variants
from prepare_phononbench import WEIGHTS, calculator, verify
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor

from basinmark.collections import prepared
from basinmark.deposit_energy import energy_per_atom
from basinmark.io import JsonStore
from basinmark.relaxation import RelaxationEngine
from basinmark.symmetry_audit import RECOVERY_SYMPREC, recovery_class


def measure(job):
    import torch

    row, encoded, weights, scored = job
    torch.set_num_threads(1)
    s = Structure.from_dict(json.loads(encoded))
    out = {k: row[k] for k in ("id", "energy", "producer", "formula", "valid", "phonon_stable")}
    if variants(s)["no_cell"] != row["producer"]["no_cell"]:
        raise ValueError("producer fingerprint mismatch")
    out.update(stages={}, errors={}, views={"producer": recovery_class(s)})
    if scored:
        out["scores"] = {"producer": energy_per_atom(s, weights) - row["hull_eV_per_atom"]}
    shear = np.eye(3)
    shear[0, 1] = 0.015
    candidate = Structure(Lattice(shear @ s.lattice.matrix), s.species, s.frac_coords)
    atoms = AseAtomsAdaptor.get_atoms(candidate)
    atoms.calc = calculator(weights)
    try:
        delta = float(atoms.get_potential_energy() / len(atoms) - row["mace_eV_per_atom"])
        accepted = bool(candidate.is_valid() and np.isfinite(delta) and abs(delta) <= 0.01)
        out.update(delta_eV_per_atom=delta, accepted=accepted)
        run = RelaxationEngine(REGISTRY, "chgnet-r2scan", fmax=0.05, steps=150).run(
            candidate if accepted else s
        )
        out.update(steps=run.steps, converged=run.converged)
        if run.converged:
            endpoint = deposit(run.structure, 5)
            out["stages"]["shear_rerelaxed"] = variants(endpoint)
            out["views"]["shear_rerelaxed"] = recovery_class(endpoint)
            out["valid_after"] = bool(run.structure.is_valid())
            if scored:
                out["scores"]["shear_rerelaxed"] = (
                    energy_per_atom(endpoint, weights) - row["hull_eV_per_atom"]
                )
        else:
            out["errors"]["shear_rerelaxed"] = "step_limit"
    except (RuntimeError, ValueError, TypeError) as exc:
        out["errors"]["shear_rerelaxed"] = type(exc).__name__
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for n in ["corpus", "archive", "weights", "out"]:
        ap.add_argument("--" + n, type=Path, required=True)
    ap.add_argument("--score", action="store_true", help="also rank the deposited endpoints")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.out.exists() or a.workers < 1:
        ap.error("choose a new output and positive worker count")
    verify(a.weights, WEIGHTS)
    meta, rows = prepared(a.corpus)
    inputs, _ = structures(None, a.archive, meta, rows)
    jobs = [(r, inputs[r["id"]], a.weights, a.score) for r in rows]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        result = list(pool.map(measure, jobs, chunksize=2))
    JsonStore(a.out).write(
        {
            "metadata": {
                **meta,
                "parent_sha256": hashlib.sha256(a.corpus.read_bytes()).hexdigest(),
                "shear_xy": 0.015,
                "max_delta_eV_per_atom": 0.01,
                "embedded_recovery_tolerance": RECOVERY_SYMPREC,
                "additional_stages": ["shear_rerelaxed"],
            },
            "rows": result,
        }
    )


if __name__ == "__main__":
    main()
