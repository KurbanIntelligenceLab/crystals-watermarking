"""Draw cross-catalogue transfer success against observed outputs for two fingerprints."""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
from scipy.stats import binom

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from basinmark.io import JsonStore  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = JsonStore(args.results).read_object()
    results = data["results"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    styles = {"full": ("#2379A8", "o", "With cell shape"), "no_cell": ("#C66736", "s", "No cell shape")}
    grid = np.arange(1, 1131)
    for mode, (color, marker, label) in styles.items():
        r = results[mode]
        ns = np.array([c["observed_outputs"] for c in r["curve"]])
        rate = np.array([c["replay_allowed"]["rate"] for c in r["curve"]]) * 100
        lo = np.array([c["replay_allowed"]["ci95"][0] for c in r["curve"]]) * 100
        hi = np.array([c["replay_allowed"]["ci95"][1] for c in r["curve"]]) * 100
        ax.errorbar(ns, rate, yerr=[rate - lo, hi - rate], fmt=marker, color=color, ms=5,
                    capsize=2, lw=1, label=f"{label}, observed")
        share = r["min_energy_of_8_matched_share"]
        law = binom.sf(r["classes_needed"] - 1, grid, share) * 100
        ax.plot(grid, law, color=color, lw=1.2, ls="--",
                label=f"{label}, law (s = {share:.3f})")
        distinct = np.array([c["matcher_distinct"]["rate"] for c in r["curve"]]) * 100
        ax.plot(ns, distinct, marker=marker, color=color, lw=0, ms=4, mfc="white",
                label=f"{label}, exposed geometry excluded")
    ax.set_xscale("log")
    ax.set_xlim(14, 1300)
    ax.set_ylim(-3, 103)
    ax.set_xlabel("Observed marked outputs, N")
    ax.set_ylabel("Transfer success (%)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
