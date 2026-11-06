"""Rebuild the observation/reuse curves directly from genuine emitted-record trials."""

import argparse
from pathlib import Path

from basinmark.collections import rate
from basinmark.io import JsonStore


def series(root, dataset, permitted):
    source = JsonStore(root / f"{dataset}-reuse.json.gz").read_object()
    column = source["trial_columns"].index("replay_fired" if permitted else "disjoint_fired")
    points = []
    for n, entry in sorted(source["rows"].items(), key=lambda item: int(item[0])):
        records = entry["trial_records"]
        observed = rate(sum(r[column] for r in records), len(records))
        expected = entry["replay_allowed" if permitted else "disjoint_records"]
        if observed != expected:
            raise ValueError("trial outcomes disagree with the stored interval or fraction")
        points.append((int(n), observed))
    return points


def plot(root, output):
    import matplotlib.pyplot as plt

    plt.rcParams["svg.hashsalt"] = "basinmark-reuse-v1"

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.3), sharey=True, layout="constrained")
    for ax, permitted in zip(axes, (True, False), strict=True):
        for name, marker, col in [("MatterGen", "o", "#17598b"), ("CrystalFormer", "s", "#ad4120")]:
            points = series(root, name, permitted)
            x = [k for k, _ in points]
            y = [v["rate"] for _, v in points]
            errors = [
                [v["rate"] - v["ci95"][0] for _, v in points],
                [v["ci95"][1] - v["rate"] for _, v in points],
            ]
            ax.errorbar(x, y, yerr=errors, label=name, marker=marker, color=col, capsize=2)
        ax.set(
            xlabel="Outputs the observer saw",
            xscale="log",
            ylim=(-0.05, 1.06),
            title="Exposed records allowed" if permitted else "Exposed records excluded",
        )
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Fraction of collections detected")
    axes[0].legend(frameon=False)
    fig.savefig(output, metadata={"Date": None})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("results/v0.3/trials"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix.lower() != ".svg":
        parser.error("supply a new SVG output")
    plot(args.root, args.out)
