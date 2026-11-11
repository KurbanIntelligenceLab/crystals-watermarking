"""Cross-catalogue transfer success as a function of observed marked outputs.

Extends the fixed 256-observation catalogue test to a sweep over observation
counts for the full and no-cell fingerprints, and compares each empirical
curve with the exposure law P{Bin(N, s) >= 10}, where s is the share of
source records whose class also occurs in the external catalogue.
"""

import argparse
import hashlib
import heapq
import json
import math
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import ijson
import numpy as np
from catalogue_overlap import MATCH, MP, OQMD, describe, oq_structure
from prepare_collection import read_rows
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from scipy.stats import binom

from basinmark.collections import rate, tally
from basinmark.core import critical_value, select_green
from basinmark.io import JsonStore

DEFAULT_OBSERVATIONS = (16, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1130)


def exposure_law(n_observed, share, gamma=0.5, alpha=0.001, collection=200):
    """Probability that N observed outputs expose at least c distinct matched classes.

    c is the smallest class count whose all-green audit rejects at level alpha.
    Each observed output is matched with probability `share`; distinct matched
    classes are approximated by the matched count, which makes this an upper
    bound when observed classes repeat.
    """
    needed = next(k for k in range(1, collection + 1) if critical_value(k, gamma, alpha) <= k)
    return float(binom.sf(needed - 1, n_observed, share)), needed


def load_catalogues(directory, external_size, seed, workers, cache):
    if cache.exists():
        data = JsonStore(cache).read_object()
        return data["mp"], data["oq"], data["meta"]
    source = []
    for row in read_rows(directory, *MP):
        energy = float(row["e_above_hull"] or "nan")
        if math.isfinite(energy):
            source.append((row["material_id"], json.loads(row["structure"]), energy))
    path = directory / f"{OQMD[0]}.zip"
    if hashlib.sha256(path.read_bytes()).hexdigest() != OQMD[1]:
        raise ValueError("OQMD archive SHA-256 differs")
    with zipfile.ZipFile(path) as archive:
        with archive.open("all_oqmd.json") as handle:
            id_counts = Counter(ijson.items(handle, "item._oqmd_entry_id"))
        with archive.open("all_oqmd.json") as handle:
            raw = ijson.items(handle, "item", use_float=True)
            unique = (r for r in raw if id_counts[r["_oqmd_entry_id"]] == 1)
            sample = heapq.nsmallest(
                external_size,
                unique,
                key=lambda r: hashlib.sha256(f"{seed}:{r['_oqmd_entry_id']}".encode()).digest(),
            )
    external = [(str(r["_oqmd_entry_id"]), r["atoms"], None) for r in sample]
    if len(external) != external_size:
        raise ValueError("requested more records than the external archive contains")
    mp_geometry = {rid: data for rid, data, _ in source}
    oq_geometry = {rid: data for rid, data, _ in external}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        mp = [r for r in pool.map(describe, source, chunksize=32) if r[1] is not None]
        oq = [r for r in pool.map(describe, external, chunksize=32) if r[1] is not None]
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
        equivalent = sorted(
            i for i in formula_indices[formula]
            if matcher.fit(other, Structure.from_dict(mp_geometry[mp[i][0]]))
        )
        overlap.append([rid, fp, equivalent])
    meta = {
        "sources": {"mp20_test": MP, "oqmd": OQMD},
        "seed": seed,
        "mp_fingerprinted": len(mp),
        "external_sampled": len(external),
        "external_fingerprinted": len(oq),
        "external_sample_rule": "lowest SHA256(seed:entry_id)",
        "matcher": MATCH,
        "versions": {name: version(name) for name in ("pymatgen", "spglib", "numpy", "ijson")},
    }
    mp_rows = [[rid, formula, fp, energy] for rid, formula, fp, energy in mp]
    JsonStore(cache).write({"mp": mp_rows, "oq": overlap, "meta": meta})
    return mp_rows, overlap, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True,
                        help="fingerprint and overlap cache; reused when present")
    parser.add_argument("--external-size", type=int, default=27000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=31415)
    parser.add_argument("--observations", type=int, nargs="+", default=list(DEFAULT_OBSERVATIONS))
    args = parser.parse_args()
    if args.out.exists():
        parser.error("choose a new output path")
    mp, overlap, meta = load_catalogues(
        args.directory, args.external_size, args.seed, args.workers, args.cache
    )
    max_observations = len(mp) // 8
    if any(not 1 <= n <= max_observations for n in args.observations):
        parser.error(f"observations must lie in 1..{max_observations}")
    gamma, alpha, collection = 0.5, 0.001, 200
    results = {}
    for mode in ("full", "no_cell"):
        source_classes = {r[2][mode] for r in mp}
        shared = [r for r in overlap if r[1][mode] in source_classes]
        shared_classes = {r[1][mode] for r in shared}
        matched_records = sum(1 for r in mp if r[2][mode] in shared_classes)
        share = matched_records / len(mp)
        # Trial-free predictor of the emitted matched share: the lowest-energy record of each
        # key-free batch of eight. Energy preference enriches catalogued classes.
        energies = np.array([r[3] for r in mp], dtype=float)
        matched_flag = np.array([r[2][mode] in shared_classes for r in mp])
        draws = np.random.default_rng(args.seed)
        picks = []
        for _ in range(300):
            batch = draws.permutation(len(mp))[: max_observations * 8].reshape(-1, 8)
            picks.append(matched_flag[batch[np.arange(len(batch)), np.argmin(energies[batch], axis=1)]].mean())
        min_energy_share = float(np.mean(picks))
        curve = []
        for n_obs in args.observations:
            rng = np.random.default_rng(args.seed + n_obs)
            records = []
            for seed in rng.integers(0, 2**53, size=args.trials):
                trial = np.random.default_rng(seed)
                key = trial.bytes(16)
                batches = trial.permutation(len(mp))[: n_obs * 8].reshape(-1, 8)
                selected = set()
                for batch in batches:
                    idx, _ = select_green(
                        [mp[i][2][mode] for i in batch], [mp[i][3] for i in batch], key, gamma
                    )
                    selected.add(int(batch[idx]))
                observed = {mp[i][2][mode] for i in selected}
                matched = [(fp[mode], eq) for _, fp, eq in shared if fp[mode] in observed]
                entry = [int(seed), len(observed)]
                for pool_fps in (
                    [f for f, _ in matched],
                    [f for f, eq in matched if selected.isdisjoint(eq)],
                ):
                    audit = trial.choice(pool_fps, collection, replace=True) if pool_fps else []
                    entry.extend((len(pool_fps), *tally(audit, key, gamma, alpha)))
                records.append(entry)
            law, needed = exposure_law(n_obs, share, gamma, alpha, collection)
            law_emitted, _ = exposure_law(n_obs, min_energy_share, gamma, alpha, collection)
            curve.append({
                "observed_outputs": n_obs,
                "replay_allowed": rate(sum(r[5] for r in records), args.trials),
                "matcher_distinct": rate(sum(r[9] for r in records), args.trials),
                "mean_matched_classes": float(np.mean([r[2] for r in records])),
                "exposure_law_population_share": law,
                "exposure_law_min_energy_share": law_emitted,
                "trial_records": records,
            })
        results[mode] = {
            "shared_classes": len(shared_classes),
            "matching_external_records": len(shared),
            "matched_source_records": matched_records,
            "matched_share": share,
            "min_energy_of_8_matched_share": min_energy_share,
            "classes_needed": needed,
            "curve": curve,
        }
        print(f"{mode:8s} share={share:.4f} min-energy share={min_energy_share:.4f} needed={needed}")
        for c in curve:
            print(f"  N={c['observed_outputs']:5d} replay={c['replay_allowed']['rate']:.3f} "
                  f"law(pop)={c['exposure_law_population_share']:.3f} "
                  f"law(emitted)={c['exposure_law_min_energy_share']:.3f} "
                  f"distinct={c['matcher_distinct']['rate']:.3f}")
    JsonStore(args.out).write({
        **meta, "m": 8, "gamma": gamma, "alpha": alpha, "collection_size": collection,
        "trials": args.trials, "observations": args.observations,
        "trial_columns": ["seed", "observed_classes", "matched_pool", "K", "G", "reject",
                          "distinct_pool", "distinct_K", "distinct_G", "distinct_reject"],
        "results": results,
    })


if __name__ == "__main__":
    main()
