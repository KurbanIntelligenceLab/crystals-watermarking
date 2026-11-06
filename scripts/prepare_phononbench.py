"""Prepare class-disjoint CrystalFormer structures with independent phonon labels."""

import argparse
import hashlib
import json
import math
import warnings
import zipfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import numpy as np
from prepare_collection import process, variants
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor

from basinmark.deposit_energy import HULL, WEIGHTS, calculator, hull_at, references, verify
from basinmark.io import JsonStore

ARCHIVE = "d786c0618af0eac9192f276adabaddd22068f2f4cca3811ba7035e18709ad620"
PREFIX = "PhononBench/CrystalFormer-gen/mp20/relaxed/"


def measure(item):
    import torch

    rid, encoded, stable, reference, weights = item
    torch.set_num_threads(1)
    s = Structure.from_dict(json.loads(encoded))

    def energy_per_atom(structure):
        atoms = AseAtomsAdaptor.get_atoms(structure)
        atoms.calc = calculator(weights)
        return float(atoms.get_potential_energy() / len(atoms))

    if reference is None:
        return {"id": rid, "errors": {"producer": "missing_hull_reference"}}
    try:
        absolute = energy_per_atom(s)
        energy = absolute - reference
        if not math.isfinite(energy):
            raise ValueError("invalid energy")
    except (ValueError, RuntimeError, KeyError) as exc:
        return {"id": rid, "errors": {"producer": type(exc).__name__}}
    row = process((rid, energy, encoded, None, 5, "chgnet-r2scan"))
    row.update(phonon_stable=stable, mace_eV_per_atom=absolute, hull_eV_per_atom=reference)
    if "producer" not in row:
        return row
    try:
        strained = Structure(
            Lattice(np.diag([1.01, 1 / 1.01, 1.0]) @ s.lattice.matrix),
            s.species,
            s.frac_coords,
        )
        delta = energy_per_atom(strained) - absolute
        accepted = abs(delta) <= 0.010 and bool(strained.is_valid())
        row.update(strain_delta_eV_per_atom=delta, strain_accepted=accepted)
        stress = (
            process((rid, energy, strained.to_json(), None, 5, "chgnet-r2scan"))
            if (accepted)
            else row
        )
        if "rerelaxed" in stress["stages"]:
            row["stages"]["strain_rerelaxed"] = stress["stages"]["rerelaxed"]
        else:
            row["errors"]["strain_rerelaxed"] = stress["errors"].get("rerelaxed", "unavailable")
    except Exception as exc:
        row["errors"]["strain_rerelaxed"] = type(exc).__name__
    return row


def source_labels(archive):
    labels = {}
    for name in archive.namelist():
        if name.startswith(PREFIX) and name.endswith("/Lable.txt"):
            for line in archive.read(name).decode().splitlines():
                if not line.strip():
                    continue
                stem, status = line.split()
                path = name.rsplit("/", 1)[0] + "/" + stem + "_relaxed.cif"
                if path in labels or status not in {"Stable", "unStable"}:
                    raise ValueError("duplicate identity or unexpected stability label")
                labels[path] = status == "Stable"
    return labels


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ["archive", "hull", "weights", "out"]:
        ap.add_argument("--" + arg, type=Path, required=True)
    ap.add_argument("--exclude", type=Path, action="append", required=True)
    ap.add_argument("--hull-solver", choices=["diagram", "linear"], default="diagram")
    ap.add_argument("--limit", type=int, default=384)
    ap.add_argument("--seed", type=int, default=47017)
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    if a.out.exists() or min(a.workers, a.limit) < 1:
        ap.error("supply a new output and positive limits")
    for path, digest in [(a.archive, ARCHIVE), (a.hull, HULL), (a.weights, WEIGHTS)]:
        verify(path, digest)
    excluded = set()
    previous_ids = set()
    for path in a.exclude:
        old = JsonStore(path).read_object()
        excluded.update(r["producer"]["no_cell"] for r in old["rows"] if "producer" in r)
        if old["metadata"].get("source_sha256") == ARCHIVE:
            previous_ids.update(r["id"] for r in old["rows"])
    groups = references(a.hull)
    counts = Counter()
    items = []
    warnings.filterwarnings("ignore", category=UserWarning)
    with zipfile.ZipFile(a.archive) as archive:
        names = sorted(n for n in archive.namelist() if n.startswith(PREFIX) and n.endswith(".cif"))
        labels = source_labels(archive)
        if set(names) - labels.keys() or len(set(names)) != len(names):
            raise ValueError("unmatched CIFs or duplicate ZIP member")
        for i in np.random.default_rng(a.seed).permutation(len(names)):
            name = names[i]
            if name.removeprefix(PREFIX) in previous_ids:
                counts["previous_id"] += 1
                continue
            try:
                structure = Structure.from_str(archive.read(name).decode(), fmt="cif")
                fp = variants(structure)
            except Exception:
                counts["unresolved_geometry"] += 1
                continue
            if fp["no_cell"] in excluded:
                counts["previous_class"] += 1
                continue
            try:
                reference = hull_at(structure.composition, groups, a.hull_solver)
            except ValueError:
                reference = None
            items.append(
                (name.removeprefix(PREFIX), structure.to_json(), labels[name], reference, a.weights)
            )
            if len(items) == a.limit:
                break
        if len(items) != a.limit:
            raise ValueError("insufficient class-disjoint structures")
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(measure, items, chunksize=2))
    metadata = {
        "dataset": "PhononBench CrystalFormer-MP20",
        "sample_size": len(rows),
        "seed": a.seed,
        "source_records": len(names),
        "source_sha256": ARCHIVE,
        "hull_sha256": HULL,
        **({"hull_solver": a.hull_solver} if a.hull_solver != "diagram" else {}),
        "ranking": "MACE-OMAT-0-medium energy minus LeMat MACE-OMAT hull energy per atom",
        "weights_sha256": WEIGHTS,
        "engine": "chgnet-r2scan",
        "paired_dft_source": None,
        "r2scan_weights_sha256": hashlib.sha256(_r2scan_weights().read_bytes()).hexdigest(),
        "additional_stages": ["strain_rerelaxed"],
        "force_threshold": 0.05,
        "step_limit": 150,
        "cif_digits": 5,
        "excluded_classes": len(excluded),
        "screening_exclusions": dict(counts),
        "excluded_class_sha256": hashlib.sha256("\n".join(sorted(excluded)).encode()).hexdigest(),
        "versions": {
            n: version(n)
            for n in [
                "ase",
                "pymatgen",
                "chgnet",
                "mace-torch",
                "torch",
                "numpy",
                "pyarrow",
                "spglib",
            ]
        },
    }
    JsonStore(a.out).write(
        {
            "metadata": metadata,
            "rows": rows,
            "errors_by_stage": dict(Counter(k for r in rows for k in r["errors"])),
        }
    )


def _r2scan_weights():
    from importlib.metadata import distribution

    (path,) = distribution("chgnet").locate_file("chgnet/pretrained/r2scan").glob("*.pth.tar")
    return path


if __name__ == "__main__":
    main()
