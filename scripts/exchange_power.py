"""Check an affordable-class power bound for fixed and pooled selection on archived trials.

For each archived trial the script rebuilds the candidate draw and key, counts
the distinct classes within the energy window of each block minimum (fixed
rule) and of the global minimum (pooled rule), and compares the green counts
implied by those counts with the green counts the selectors actually emitted.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.stats import binom

from basinmark.allocation import select_collection
from basinmark.collections import cutoff
from basinmark.core import is_green
from basinmark.io import JsonStore

COHORTS = {
    "E": ("results/v0.6/CF-E.json.gz", 368219),
    "F": ("results/v0.6/CF-F.json.gz", 653291),
    "CF2": ("results/v0.4/CrystalFormer-second.json.gz", 91571),
    "MG2": ("results/v0.4/MatterGen-second.json.gz", 91571),
}


def affordable_counts(energy, fps, batches, window):
    """Distinct classes within the window of each block minimum and of the global minimum."""
    per_block = []
    for batch in batches:
        floor = energy[batch].min()
        per_block.append(len({fps[i] for i in batch if energy[i] <= floor + window}))
    floor = energy.min()
    pooled = len({fps[i] for i in range(len(fps)) if energy[i] <= floor + window})
    return per_block, pooled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=COHORTS, required=True)
    parser.add_argument("--trials", type=int, default=500)
    parser.add_argument("--window", type=float, default=0.1)
    parser.add_argument("--gamma", type=float, default=0.625)
    parser.add_argument("--alpha", type=float, default=0.001)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    source, seed = COHORTS[args.cohort]
    rows = JsonStore(Path(source)).read_object()["rows"]
    n, m = 32, 8
    seeds = np.random.default_rng(seed).integers(0, 2**53, args.trials).tolist()
    records = []
    for s in seeds:
        rng = np.random.default_rng(s)
        key = rng.bytes(16)
        rng.bytes(16)
        draw = rng.permutation(len(rows))[: n * m]
        fps = [rows[i]["producer"]["no_cell"] for i in draw]
        energy = np.array([rows[i]["energy"] for i in draw], dtype=float)
        batches = np.arange(n * m).reshape(n, m)
        per_block, pooled = affordable_counts(energy, fps, batches, args.window)
        green = {fp: is_green(fp, key, args.gamma) for fp in set(fps)}
        floor = energy.min()
        affordable_green = len({fps[i] for i in range(n * m)
                                if energy[i] <= floor + args.window and green[fps[i]]})
        bound = min(n, affordable_green)
        # Nested model: slot i (sorted block minima b_i) may draw any class below b_i + window.
        minima = np.sort(energy[batches].min(axis=1))
        nested_green = [len({fps[i] for i in range(n * m)
                             if energy[i] <= b + args.window and green[fps[i]]}) for b in minima]
        nested_model = max([g for g in range(1, n + 1) if nested_green[g - 1] >= g], default=0)
        # Which cap binds when the budget rule emits fewer than n greens: test the n cheapest
        # distinct green classes against the per-slot ceilings and the mean-energy ceiling.
        order = np.argsort(energy, kind="stable")
        seen, cheapest_green = set(), []
        for i in order:
            if green[fps[i]] and fps[i] not in seen:
                seen.add(fps[i]); cheapest_green.append(i)
            if len(cheapest_green) == n:
                break
        if len(cheapest_green) == n:
            e = np.sort(energy[cheapest_green])
            fixed_sel = select_collection(fps, energy, key, policy="fixed", m=m,
                                          energy_window=args.window, gamma=args.gamma)
            fixed_mean = float(np.sort(energy[fixed_sel.indices]).mean())
            pooled_mean = float(np.sort(energy)[:n].mean())
            slot_ok = bool(np.all(e <= minima + args.window))
            mean_ok = bool(e.mean() <= min(fixed_mean, pooled_mean + args.window))
        else:
            slot_ok = mean_ok = False
        caps = {"n_green_available": len(cheapest_green) == n, "slot_ok": slot_ok, "mean_ok": mean_ok}
        emitted = {}
        for policy in ("fixed", "budget", "pool"):
            sel = select_collection(fps, energy, key, policy=policy, m=m,
                                    energy_window=args.window, gamma=args.gamma)
            classes = {fps[i] for i in sel.indices}
            emitted[policy] = [len(classes), sum(green[c] for c in classes)]
        records.append({
            "seed": int(s), "per_block_affordable": per_block, "pooled_affordable": pooled,
            "affordable_green_classes": affordable_green, "bound": bound,
            "nested_green": nested_green, "nested_model": nested_model, "caps": caps,
            "emitted": emitted,
        })
    per_block = np.array([r["per_block_affordable"] for r in records])
    pooled = np.array([r["pooled_affordable"] for r in records])
    fixed_g = np.array([r["emitted"]["fixed"][1] for r in records])
    budget_g = np.array([r["emitted"]["budget"][1] for r in records])
    pool_g = np.array([r["emitted"]["pool"][1] for r in records])
    bound = np.array([r["bound"] for r in records])
    nested = np.array([r["nested_model"] for r in records])
    short = [r for r in records if r["emitted"]["budget"][1] < n]
    cap_attribution = {
        "trials_below_32_green": len(short),
        "slot_cap_violated": sum(not r["caps"]["slot_ok"] for r in short),
        "mean_cap_violated": sum(not r["caps"]["mean_ok"] for r in short),
        "fewer_than_32_green_classes": sum(not r["caps"]["n_green_available"] for r in short),
    }
    c = cutoff(n, args.gamma, args.alpha)
    # Model-implied producer power: fixed = P{sum_i Bernoulli(1-(1-gamma)^{a_i}) >= c};
    # pooled = P{min(n, Bin(a, gamma)) >= c}, with a the pooled affordable class count.
    pooled_power_model = float(np.mean([
        binom.sf(c - 1, a, args.gamma) if a >= c else 0.0 for a in pooled
    ]))
    fixed_power_model = float(np.mean([
        _poisson_binomial_tail(1 - (1 - args.gamma) ** np.array(a, dtype=float), c)
        for a in per_block
    ]))
    summary = {
        "cohort": args.cohort, "trials": args.trials, "window": args.window, "gamma": args.gamma,
        "cutoff_at_32": int(c),
        "per_block_affordable_mean": float(per_block.mean()),
        "per_block_affordable_is_one_fraction": float((per_block == 1).mean()),
        "pooled_affordable_mean": float(pooled.mean()),
        "pooled_affordable_min": int(pooled.min()),
        "bound_holds_budget": float((budget_g >= bound).mean()),
        "bound_holds_pool": float((pool_g >= bound).mean()),
        "bound_slack_budget_mean": float((budget_g - bound).mean()),
        "nested_model_mean": float(nested.mean()),
        "nested_holds_budget": float((budget_g >= nested).mean()),
        "nested_holds_pool": float((pool_g >= nested).mean()),
        "nested_slack_budget_mean": float((budget_g - nested).mean()),
        "nested_slack_pool_mean": float((pool_g - nested).mean()),
        "nested_power_model": float((nested >= c).mean()),
        "cap_attribution": cap_attribution,
        "fixed_green_mean": float(fixed_g.mean()),
        "budget_green_mean": float(budget_g.mean()),
        "pool_green_mean": float(pool_g.mean()),
        "fixed_producer_power_observed": float((fixed_g >= c).mean()),
        "budget_producer_power_observed": float((budget_g >= c).mean()),
        "fixed_power_model": fixed_power_model,
        "pooled_power_model": pooled_power_model,
    }
    for k, v in summary.items():
        print(f"{k:40s} {v}")
    if args.out:
        JsonStore(args.out).write({"summary": summary, "records": records})


def _poisson_binomial_tail(p, c):
    """P{sum of independent Bernoulli(p_i) >= c} by direct convolution."""
    dist = np.zeros(len(p) + 1)
    dist[0] = 1.0
    for q in p:
        dist[1:] = dist[1:] * (1 - q) + dist[:-1] * q
        dist[0] *= 1 - q
    return float(dist[c:].sum())


if __name__ == "__main__":
    main()
