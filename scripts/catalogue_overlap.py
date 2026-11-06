"""Measure within-formula fingerprint collisions and blind transfer across catalogues."""

import argparse
import hashlib
import heapq
import json
import math
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from itertools import combinations
from pathlib import Path

import ijson
import numpy as np
from prepare_collection import read_rows, variants
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from basinmark.collections import rate, tally
from basinmark.core import select_green
from basinmark.io import JsonStore

MP = (65498655, "bc4507c7e8569e9c684e031e2ef3a9effd4e84fe1a7d25ebec579fbe98d3972d")
OQMD = (26790182, "3b1f3112ab9c4a9388937b60a802411072574605b0bc48e32876404c7539b85d")
MATCH = dict(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=True)


def oq_structure(atoms):
    return Structure(
        atoms["lattice_mat"],
        [e.strip() for e in atoms["elements"]],
        atoms["coords"],
        coords_are_cartesian=atoms["cartesian"],
    )


def describe(item):
    rid, data, energy = item
    try:
        structure = Structure.from_dict(data) if "sites" in data else oq_structure(data)
        return (rid, structure.composition.reduced_formula, variants(structure), energy)
    except Exception as exc:
        return (rid, None, type(exc).__name__, energy)


def compare_formula(group):
    matcher = StructureMatcher(**MATCH)
    rows = [(rid, Structure.from_dict(data), fp) for rid, data, fp in group]
    counts = Counter()
    for (_, a, fa), (_, b, fb) in combinations(rows, 2):
        match = matcher.fit(a, b)
        label = "equivalent" if match else "distinct"
        counts[label] += 1
        for mode in fa:
            if fa[mode] == fb[mode]:
                counts[f"{label}_{mode}_collision"] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--external-size", type=int, default=27000)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=31415)
    args = parser.parse_args()
    if args.out.exists() or min(args.external_size, args.workers, args.trials) <= 0:
        parser.error("choose a new output and positive sample, worker and trial counts")
    source = []
    for row in read_rows(args.directory, *MP):
        energy = float(row["e_above_hull"] or "nan")
        if math.isfinite(energy):
            source.append((row["material_id"], json.loads(row["structure"]), energy))
    path = args.directory / f"{OQMD[0]}.zip"
    if hashlib.sha256(path.read_bytes()).hexdigest() != OQMD[1]:
        raise ValueError("OQMD archive SHA-256 differs")
    with zipfile.ZipFile(path) as archive:
        with archive.open("all_oqmd.json") as handle:
            id_counts = Counter(ijson.items(handle, "item._oqmd_entry_id"))
        with archive.open("all_oqmd.json") as handle:
            raw = ijson.items(handle, "item", use_float=True)
            unique = (r for r in raw if id_counts[r["_oqmd_entry_id"]] == 1)
            sample = heapq.nsmallest(
                args.external_size,
                unique,
                key=lambda r: hashlib.sha256(
                    f"{args.seed}:{r['_oqmd_entry_id']}".encode()
                ).digest(),
            )
    external = [(str(r["_oqmd_entry_id"]), r["atoms"], None) for r in sample]
    if len(external) != args.external_size:
        raise ValueError("requested more records than the external archive contains")
    if any(len({x[0] for x in records}) != len(records) for records in (source, external)):
        raise ValueError("duplicate catalogue record ID")
    mp_geometry = {rid: data for rid, data, _ in source}
    oq_geometry = {rid: data for rid, data, _ in external}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        mp = [r for r in pool.map(describe, source, chunksize=32) if r[1] is not None]
        oq = [r for r in pool.map(describe, external, chunksize=32) if r[1] is not None]
        formulas = defaultdict(list)
        for rid, formula, fp, _ in mp:
            formulas[formula].append((rid, mp_geometry[rid], fp))
        pooled = Counter()
        groups = [v for v in formulas.values() if len(v) > 1]
        for counts in pool.map(compare_formula, groups, chunksize=1):
            pooled.update(counts)

    matcher = StructureMatcher(**MATCH)
    coarse = {r[2]["no_cell"] for r in mp}
    formula_indices = defaultdict(list)
    for i, (_, formula, _, _) in enumerate(mp):
        formula_indices[formula].append(i)
    overlap = []
    for rid, formula, fp, _ in oq:
        if fp["no_cell"] not in coarse:
            continue
        other = oq_structure(oq_geometry[rid])
        equivalent = {
            i
            for i in formula_indices[formula]
            if matcher.fit(other, Structure.from_dict(mp_geometry[mp[i][0]]))
        }
        overlap.append((rid, fp, equivalent))

    summary = {}
    trials = {}
    for mode in ("full", "no_cell"):
        source_classes = {r[2][mode] for r in mp}
        shared = [r for r in overlap if r[1][mode] in source_classes]
        rng = np.random.default_rng(args.seed)
        records = []
        for seed in rng.integers(0, 2**53, size=args.trials):
            trial = np.random.default_rng(seed)
            key = trial.bytes(16)
            batches = trial.permutation(len(mp))[: len(mp) // 8 * 8].reshape(-1, 8)
            selected = set()
            for batch in batches:
                idx, _ = select_green(
                    [mp[i][2][mode] for i in batch], [mp[i][3] for i in batch], key
                )
                selected.add(batch[idx])
            observed = {mp[i][2][mode] for i in selected}
            matched = [(fp[mode], eq) for _, fp, eq in shared if fp[mode] in observed]
            entry = [int(seed)]
            for pool_fps in (
                [f for f, _ in matched],
                [f for f, eq in matched if selected.isdisjoint(eq)],
            ):
                external_audit = trial.choice(pool_fps, 200, replace=True) if pool_fps else []
                entry.extend((len(pool_fps), *tally(external_audit, key, 0.5, 0.001)))
            records.append(entry)
        trials[mode] = records
        summary[mode] = {
            "shared_classes": len({fp[mode] for _, fp, _ in shared}),
            "matching_external_records": len(shared),
            "cross_catalogue_reuse": rate(sum(r[4] for r in records), args.trials),
            "matcher_distinct_reuse": rate(sum(r[8] for r in records), args.trials),
        }
    JsonStore(args.out).write(
        {
            "sources": {"mp20_test": MP, "oqmd": OQMD},
            "seed": args.seed,
            "mp_finite_energy": len(source),
            "mp_fingerprinted": len(mp),
            "external_archive_rows": sum(id_counts.values()),
            "external_repeated_id_rows_excluded": sum(n for n in id_counts.values() if n > 1),
            "external_sampled": len(external),
            "external_fingerprinted": len(oq),
            "same_formula_groups_with_pairs": len(groups),
            "within_formula_pairs": dict(pooled),
            "external_sample_rule": "lowest SHA256(seed:entry_id)",
            "matcher": MATCH,
            "versions": {name: version(name) for name in ("pymatgen", "spglib", "numpy", "ijson")},
            "m": 8,
            "gamma": 0.5,
            "alpha": 0.001,
            "collection_size": 200,
            "observed_outputs": len(mp) // 8,
            "summary": summary,
            "trial_columns": [
                "seed",
                "matched_pool",
                "K",
                "G",
                "reject",
                "distinct_pool",
                "distinct_K",
                "distinct_G",
                "distinct_reject",
            ],
            "trial_records": trials,
        }
    )


if __name__ == "__main__":
    main()
