"""Test an optional reference-geometry check on all repeated-formula MP-20 pairs."""

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np
from catalogue_overlap import MATCH, MP
from prepare_collection import read_rows, variants
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure

from basinmark.io import JsonStore
from basinmark.reference import neighbor_profile, profile_distance

SPLITS = {
    "test": MP,
    "validation": (65498670, "301af7e18d32138612fb832c848f8e1d9237a7dd91dfb8ce1b9c7becbbe076b2"),
}
THRESHOLD = 0.03


def compare(group):
    matcher = StructureMatcher(**MATCH)
    prepared = []
    strain = []
    for rid, encoded in group:
        s = Structure.from_dict(json.loads(encoded))
        fp = variants(s)
        profile = neighbor_profile(s)
        changed = Structure(
            Lattice(np.diag([1.01, 1 / 1.01, 1]) @ s.lattice.matrix), s.species, s.frac_coords
        )
        strain.append(profile_distance(profile, neighbor_profile(changed)))
        prepared.append((rid, s, fp, profile))
    pairs = []
    for (i, a, fa, pa), (j, b, fb, pb) in combinations(prepared, 2):
        pairs.append(
            [
                i,
                j,
                bool(matcher.fit(a, b)),
                fa["full"] == fb["full"],
                fa["no_cell"] == fb["no_cell"],
                profile_distance(pa, pb),
            ]
        )
    return pairs, strain


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--split", choices=SPLITS, required=True)
    ap.add_argument("--directory", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    if a.out.exists() or a.workers <= 0:
        ap.error("choose a new output and a positive worker count")
    groups = defaultdict(list)
    ids = set()
    for r in read_rows(a.directory, *SPLITS[a.split]):
        if not np.isfinite(float(r["e_above_hull"] or "nan")):
            continue
        rid = r["material_id"]
        if rid in ids:
            raise ValueError("duplicate MP-20 ID")
        ids.add(rid)
        s = Structure.from_dict(json.loads(r["structure"]))
        groups[s.composition.reduced_formula].append((rid, r["structure"]))
    pairs, strain = [], []
    repeated = [g for g in groups.values() if len(g) > 1]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for compared, deformed in pool.map(compare, repeated, chunksize=4):
            pairs.extend(compared)
            strain.extend(deformed)
    counts = defaultdict(Counter)
    for _, _, equivalent, full, coarse, distance in pairs:
        group = counts["equivalent" if equivalent else "distinct"]
        group.update(
            total=1,
            full=int(full),
            no_cell=int(coarse),
            no_cell_plus_profile=int(coarse and distance <= THRESHOLD),
        )
    JsonStore(a.out).write(
        {
            "source": SPLITS[a.split],
            "split": a.split,
            "threshold": THRESHOLD,
            "matcher": MATCH,
            "source_ids_sha256": hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest(),
            "eligible_structures": len(ids),
            "formula_groups_with_pairs": len(repeated),
            "summary": counts,
            "pairs": pairs,
            "pair_columns": [
                "id_1",
                "id_2",
                "matcher_equivalent",
                "full_equal",
                "no_cell_equal",
                "profile_rms",
            ],
            "profile_strain": {
                "tested": len(strain),
                "passed": sum(d <= THRESHOLD for d in strain),
                "max_distance": max(strain),
            },
        }
    )


if __name__ == "__main__":
    main()
