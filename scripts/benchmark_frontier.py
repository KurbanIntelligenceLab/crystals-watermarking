"""Reconstruct the exploratory energy and audit comparisons on the first external cohort."""

import argparse
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.special import gammaincc
from scipy.stats import irwinhall

from basinmark.allocation import _fixed
from basinmark.collections import cutoff, prepared
from basinmark.core import prf_uniform
from basinmark.io import JsonStore


@lru_cache(maxsize=64)
def uniform_sum_cutoff(k):
    return float(irwinhall.ppf(0.001, k))


def evaluate(corpus):
    _, rows = prepared(corpus)
    samples = {}
    for seed in range(88211, 88411):
        rng = np.random.default_rng(seed)
        key = rng.bytes(16)
        batch = [rows[i] for i in rng.permutation(len(rows))[:256]]
        energy = np.array([r["energy"] for r in batch])
        scores = np.array([prf_uniform(r["producer"]["no_cell"], key) for r in batch])
        order = np.argsort(energy, kind="stable")
        blocks = np.arange(256).reshape(32, 8)
        reference = energy[order[:32]]
        minima = np.sort(energy[blocks].min(axis=1))

        def observed(indices):
            fps = {
                batch[i]["stages"]["rerelaxed"]["no_cell"]
                for i in indices
                if "rerelaxed" in batch[i]["stages"]
            }
            return np.array([prf_uniform(f, key) for f in sorted(fps)])

        def record(name, indices, hit):
            samples.setdefault(name, []).append([int(hit), float(energy[indices].mean())])

        for gamma in (0.125, 0.25, 0.375, 0.5, 0.625, 0.75):
            labels = scores < gamma
            fixed = _fixed(energy, labels, blocks, 0.1)
            fx = energy[fixed.indices]
            strict = np.minimum(reference + 0.1, np.sort(fx))
            greens, reds = order[labels[order]], order[~labels[order]]
            for mode, ceilings, budget in (
                ("strict", strict, np.inf),
                ("slots", minima + 0.1, min(fx.mean(), reference.mean() + 0.1)),
                ("original", minima + 0.1, fx.mean()),
            ):
                for n in range(min(32, len(greens)), max(0, 32 - len(reds)) - 1, -1):
                    ids = np.concatenate((greens[:n], reds[: 32 - n]))
                    if np.all(np.sort(energy[ids]) <= ceilings) and (
                        energy[ids].sum() <= 32 * budget + 1e-9
                    ):
                        break
                u = observed(ids)
                record(f"{mode}-{gamma}", ids, sum(u < gamma) >= cutoff(len(u), gamma, 0.001))
        fixed = _fixed(energy, scores < 0.5, blocks, 0.1)
        for mode, ceilings in (
            ("strict", np.minimum(reference + 0.1, np.sort(energy[fixed.indices]))),
            ("slots", np.minimum(minima + 0.1, np.sort(energy[fixed.indices]))),
            ("pooled_window", reference + 0.1),
        ):
            chosen = []
            available = np.ones(256, dtype=bool)
            for ceiling in ceilings:
                feasible = np.flatnonzero(available & (energy <= ceiling))
                pick = feasible[np.argmin(scores[feasible])]
                available[pick] = False
                chosen.append(pick)
            u = observed(chosen)
            record(mode + "-fisher", chosen, gammaincc(len(u), -np.log(u).sum()) <= 0.001)
            record(mode + "-uniform_sum", chosen, u.sum() <= uniform_sum_cutoff(len(u)))
    return {
        "first_seed": 88211,
        "trials": 200,
        "columns": ["reject", "mean_energy_eV_per_atom"],
        "records": samples,
        "rejections": {k: sum(v[0] for v in t) for k, t in samples.items()},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--corpus", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    if a.out.exists():
        p.error("choose a new output")
    JsonStore(a.out).write(evaluate(a.corpus))


if __name__ == "__main__":
    main()
