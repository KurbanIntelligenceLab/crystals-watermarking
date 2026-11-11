"""Collection experiments using measured structural fingerprints and producer energies."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

from .core import (
    binom_sf,
    clopper_pearson,
    critical_value,
    is_green,
    n_forge,
    select_green,
)
from .io import FingerprintCorpus, JsonStore


def prepared(path):
    data = JsonStore(Path(path)).read_object()
    rows = data.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("expected a nonempty prepared collection")
    if not isinstance(data.get("metadata"), dict) or any(
        not isinstance(r, dict) or not isinstance(r.get("id"), str) for r in rows
    ):
        raise ValueError("expected metadata and one string ID per record")
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("duplicate record IDs")
    if any("producer" not in r and "producer" not in r.get("errors", {}) for r in rows):
        raise ValueError("each record needs producer fingerprints or a recorded failure")
    rows = [r for r in rows if "producer" in r]
    if not rows or any(not math.isfinite(r["energy"]) for r in rows):
        raise ValueError("expected usable producer fingerprints and finite energies")
    return data["metadata"], rows


def rate(hits, n):
    return {
        "rejections": hits,
        "trials": n,
        "rate": hits / n,
        "ci95": list(clopper_pearson(hits, n)),
    }


@lru_cache(maxsize=4096)
def cutoff(k, gamma, alpha):
    return critical_value(k, gamma, alpha)


def tally(fps, key, gamma, alpha):
    units = set(fps)
    green = sum(is_green(f, key, gamma) for f in units)
    return len(units), green, int(green >= cutoff(len(units), gamma, alpha))


def matching_indices(corpus, selected, disjoint):
    """Construct the record pool using only released identities and fingerprints."""
    exposed = {corpus[i] for i in selected}
    excluded = set(selected) if disjoint else set()
    return [i for i, f in enumerate(corpus) if f in exposed and i not in excluded]


def reuse(a):
    source = Path(a.corpus or a.fingerprints)
    if a.corpus:
        _, rows = prepared(source)
        fps = [r["producer"][a.representation] for r in rows]
        energies = np.array([r["energy"] for r in rows])
    else:
        fps = FingerprintCorpus(source).load()
        energies = np.zeros(len(fps))
    if max(a.observations) * a.m > len(fps):
        raise ValueError("each observation requires a distinct batch of m candidates")
    rng = np.random.default_rng(a.seed)
    out = {
        "schema": 2,
        "ranked_input_sha256": hashlib.sha256(
            json.dumps(
                list(zip(fps, energies.tolist())),
                ensure_ascii=False,
                separators=(",", ":"),
                allow_nan=False,
            ).encode()
        ).hexdigest(),
        "ranking": "producer_energy" if a.corpus else "input_order_with_equal_energies",
        "representation": a.representation if a.corpus else "supplied",
        "seed": a.seed,
        "m": a.m,
        "gamma": a.gamma,
        "alpha": a.alpha,
        "collection_size": a.collection,
        "N_forge_predicted": n_forge(a.gamma, a.alpha),
        "trial_columns": [
            "trial_seed",
            "selected_row_indices",
            "fallbacks",
            "replay_pool",
            "replay_K",
            "replay_green",
            "replay_fired",
            "disjoint_pool",
            "disjoint_K",
            "disjoint_green",
            "disjoint_fired",
        ],
        "rows": {},
    }
    for n in a.observations:
        trials = []
        for seed in rng.integers(0, 2**53, size=a.trials):
            trial = np.random.default_rng(seed)
            key = trial.bytes(16)
            batches = trial.choice(len(fps), n * a.m, replace=False).reshape(n, a.m)
            selected = []
            fallbacks = 0
            for batch in batches:
                idx, green = select_green([fps[i] for i in batch], energies[batch], key, a.gamma)
                selected.append(int(batch[idx]))
                fallbacks += int(not green)
            record = [int(seed), selected, fallbacks]
            for disjoint in (False, True):
                matches = matching_indices(fps, selected, disjoint)
                sample = trial.choice(matches, a.collection, replace=True) if matches else []
                record.extend(
                    (len(matches), *tally([fps[i] for i in sample], key, a.gamma, a.alpha))
                )
            trials.append(record)
        replay = rate(sum(r[6] for r in trials), a.trials)
        separate = rate(sum(r[10] for r in trials), a.trials)
        out["rows"][n] = {
            "detected": replay["rate"],
            "ci95": replay["ci95"],
            "trials": a.trials,
            "replay_allowed": replay,
            "disjoint_records": separate,
            "trial_records": trials,
            "infeasible_disjoint": sum(r[7] == 0 for r in trials),
        }
    JsonStore(Path(a.out)).write(out)


def key_controls(a):
    """Independent audit keys for fixed, key-independent batches of measured records."""
    meta, rows = prepared(a.corpus)
    rng = np.random.default_rng(a.seed)
    fps = np.array([r["producer"][a.representation] for r in rows])
    energies = np.array([r["energy"] for r in rows])
    perm = rng.permutation(len(fps))[: len(fps) // a.m * a.m]
    if not len(perm):
        raise ValueError("fewer than m usable producer records")
    batches = perm.reshape(-1, a.m)
    ranked = np.take_along_axis(batches, np.argsort(energies[batches], kind="stable"), axis=1)
    classes, ids = np.unique(fps, return_inverse=True)
    ranked = ids[ranked]
    unmarked = np.unique(ranked[:, 0])
    # Fixed batches are shared across trials. Both keys are freshly sampled per trial.
    trials = []
    for seed in rng.integers(0, 2**53, size=a.keys):
        trial = np.random.default_rng(seed)
        mark_key, other_key = trial.bytes(16), trial.bytes(16)
        labels = np.fromiter((is_green(f, mark_key, a.gamma) for f in classes), dtype=bool)
        batch_labels = labels[ranked]
        # First green in energy order or the lowest energy when all labels are false.
        selected = np.unique(ranked[np.arange(len(ranked)), batch_labels.argmax(axis=1)])
        null_g = int(labels[unmarked].sum())
        other_g = sum(is_green(classes[i], other_key, a.gamma) for i in selected)
        trials.append([int(seed), len(unmarked), null_g, len(selected), other_g])
    summary = {}
    for name, ki, gi in (("unmarked", 1, 2), ("independent_key", 3, 4)):
        counts = Counter(r[ki] for r in trials)
        attainable = (
            sum(
                n * binom_sf(cutoff(k, a.gamma, a.alpha) - 1, k, a.gamma) for k, n in counts.items()
            )
            / a.keys
        )
        summary[name] = {
            **rate(sum(r[gi] >= cutoff(r[ki], a.gamma, a.alpha) for r in trials), a.keys),
            "mean_exact_null_rate": attainable,
        }
    JsonStore(Path(a.out)).write(
        {
            "schema": 1,
            "source": meta,
            "representation": a.representation,
            "m": a.m,
            "seed": a.seed,
            "gamma": a.gamma,
            "alpha": a.alpha,
            "batches": len(ranked),
            "trial_columns": ["trial_seed", "unmarked_K", "unmarked_G", "marked_K", "other_key_G"],
            "summary": summary,
            "trial_records": trials,
        }
    )


def _arm_trial(rows, batches, key, a, policy, representation, stages):
    selected = []
    penalties = []
    fallbacks = 0
    for batch in batches:
        candidates = [rows[i] for i in batch]
        energies = [r["energy"] for r in candidates]
        best = int(np.argmin(energies))
        if policy == "single":
            index = 0
        elif policy == "minimum":
            index = best
        else:
            window = a.energy_window if policy == "bounded" else None
            index, green = select_green(
                [r["producer"][representation] for r in candidates],
                energies,
                key,
                a.gamma,
                energy_window=window,
            )
            fallbacks += int(not green)
        selected.append(candidates[index])
        penalties.append(energies[index] - energies[best])
    audits = {}
    for stage in stages:
        observed = [
            r["producer"] if stage == "producer" else r["stages"].get(stage) for r in selected
        ]
        fps = [value[representation] for value in observed if value is not None]
        audits[stage] = [len(fps), *tally(fps, key, a.gamma, a.alpha)]
    hull = [r["dft_hull"] for r in selected if r.get("dft_hull") is not None]
    return {
        "fallbacks": fallbacks if policy in ("marked", "bounded") else None,
        "energy_penalty_eV_per_atom": penalties,
        "producer_valid": sum(r["valid"] for r in selected),
        "rerelaxed_valid": sum(r.get("valid_after", False) for r in selected),
        "dft_evaluated": len(hull),
        "dft_stable": sum(e <= a.stability_threshold for e in hull),
        "audits": audits,
    }


def collection_benchmark(a):
    meta, rows = prepared(a.corpus)
    sizes = a.sizes
    if min(sizes) < 1 or len(set(sizes)) != len(sizes) or max(sizes) * a.m > len(rows):
        raise ValueError("audit sizes must be unique positive counts supported by the corpus")
    if not math.isfinite(a.energy_window) or a.energy_window < 0:
        raise ValueError("energy window must be finite and nonnegative")
    if not math.isfinite(a.stability_threshold):
        raise ValueError("stability threshold must be finite")
    stages = ["producer", "cif", *meta.get("additional_stages", [])]
    if meta["engine"]:
        stages.append("rerelaxed")
    if meta["paired_dft_source"]:
        stages.append("paired_dft")
    out = {
        "schema": 1,
        "source": meta,
        "eligible_records": len(rows),
        "m": a.m,
        "seed": a.seed,
        "gamma": a.gamma,
        "alpha": a.alpha,
        "energy_window_eV_per_atom": a.energy_window,
        "stability_threshold_eV_per_atom": a.stability_threshold,
        "audit_columns": ["deposited_records", "K", "green_classes", "rejected"],
        "summary": {},
        "trial_seeds": [],
        "trials": {},
    }
    seeds = np.random.default_rng(a.seed).integers(0, 2**53, size=a.trials)
    out["trial_seeds"] = [int(s) for s in seeds]
    for size in sizes:
        cells = {}
        for seed in seeds:
            trial = np.random.default_rng(seed)
            key = trial.bytes(16)
            # Nested collection sizes and paired methods within each replicate.
            batches = trial.permutation(len(rows))[: size * a.m].reshape(size, a.m)
            for representation in ("full", "no_cell"):
                policies = ["marked", "minimum", "single"]
                if representation == "no_cell":
                    policies.append("bounded")
                for policy in policies:
                    name = f"{representation}/{policy}"
                    cell = cells.setdefault(name, [])
                    cell.append(_arm_trial(rows, batches, key, a, policy, representation, stages))
        out["trials"][size] = cells
        summary = {}
        for name, records in cells.items():
            penalty = np.array([r["energy_penalty_eV_per_atom"] for r in records])
            n = len(records) * size
            result = {
                "producer_candidates_per_emission": 1 if name.endswith("/single") else a.m,
                "mean_penalty_meV_per_atom": float(penalty.mean() * 1000),
                "p95_penalty_meV_per_atom": float(np.percentile(penalty, 95) * 1000),
                "max_penalty_meV_per_atom": float(penalty.max() * 1000),
                "valid_producer_fraction": sum(r["producer_valid"] for r in records) / n,
                "stages": {},
            }
            if records[0]["fallbacks"] is not None:
                result["fallback_fraction"] = sum(r["fallbacks"] for r in records) / n
            if meta["paired_dft_source"]:
                result["dft_quality_coverage"] = sum(r["dft_evaluated"] for r in records) / n
                result["confirmed_stable_yield"] = sum(r["dft_stable"] for r in records) / n
            if meta["engine"]:
                result["valid_rerelaxed_yield"] = sum(r["rerelaxed_valid"] for r in records) / n
            for stage in stages:
                obs, k, g, hits = np.array([r["audits"][stage] for r in records]).T
                result["stages"][stage] = {
                    **rate(int(hits.sum()), len(records)),
                    "record_coverage": float(obs.mean() / size),
                    "mean_K": float(k.mean()),
                    "class_green_fraction": float(g.sum() / k.sum()) if k.sum() else None,
                }
            summary[name] = result
        out["summary"][size] = summary
    JsonStore(Path(a.out)).write(out)
