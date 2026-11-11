"""Measure within-formula collisions for the primary and recovery views on labeled MP-20 pairs."""

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from benchmark_reference import SPLITS
from prepare_collection import read_rows
from pymatgen.core import Structure

from basinmark.core import fingerprint_structure
from basinmark.io import JsonStore
from basinmark.symmetry_audit import recovery_class


def compare(directory, split, reference):
    original = JsonStore(reference).read_object()
    pairs = original["pairs"]
    ids = {i for a, b, *_ in pairs for i in (a, b)}
    data = {}
    for r in read_rows(directory, *SPLITS[split]):
        rid = r["material_id"]
        if rid in ids:
            s = Structure.from_dict(json.loads(r["structure"]))
            primary = fingerprint_structure(s, use_cell_shape=False)
            data[rid] = primary, recovery_class(s)
    if set(data) != ids:
        raise ValueError("missing source structure")
    counts = defaultdict(Counter)
    records = []
    for a, b, equiv, _, original_match, _ in pairs:
        p, q = data[a], data[b]
        same = p[0] == q[0]
        if same != original_match:
            raise ValueError("original same-formula comparison differs")
        left = p[0] == q[1]
        right = q[0] == p[1]
        records.append([a, b, equiv, same, left, right])
        counts["equivalent" if equiv else "distinct"].update(
            pairs=1,
            primary=same,
            either_direction=same or left or right,
            added=(not same) and (left or right),
        )
    return {
        "split": split,
        "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
        "columns": [
            "id_a",
            "id_b",
            "equivalent",
            "primary_match",
            "a_to_secondary_b",
            "b_to_secondary_a",
        ],
        "summary": dict(counts),
        "records": records,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--split", choices=SPLITS, required=True)
    for n in ("directory", "reference", "out"):
        ap.add_argument("--" + n, type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        ap.error("choose a new output")
    JsonStore(a.out).write(compare(a.directory, a.split, a.reference))


if __name__ == "__main__":
    main()
