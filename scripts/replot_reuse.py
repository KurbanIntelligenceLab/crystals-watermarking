"""Rebuild the reuse panel and, optionally, export the actual observer inputs."""

import argparse
import hashlib
from pathlib import Path

import numpy as np

from basinmark.collections import prepared, rate
from basinmark.core import select_green
from basinmark.io import JsonStore

GENERATORS = (("MatterGen", "o", "#17598b", -0.22), ("CrystalFormer", "s", "#ad4120", 0.22))


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


def observer_records(trials_path, measurements_path):
    trials = JsonStore(trials_path).read_object()
    _, rows = prepared(measurements_path)
    fps = [r["producer"][trials["representation"]] for r in rows]
    energies = np.array([r["energy"] for r in rows])
    columns = trials["trial_columns"]
    by_n = {}
    for n, entry in trials["rows"].items():
        records = []
        for raw in entry["trial_records"]:
            row = dict(zip(columns, raw, strict=True))
            rng = np.random.default_rng(row["trial_seed"])
            key = rng.bytes(16)
            m = trials["m"]
            batches = rng.choice(len(rows), int(n) * m, replace=False).reshape(-1, m)
            indices, fallbacks = [], []
            for position, batch in enumerate(batches, start=1):
                selected, green = select_green(
                    [fps[i] for i in batch], energies[batch], key, trials["gamma"]
                )
                indices.append(int(batch[selected]))
                if not green:
                    fallbacks.append(position)
            if indices != row["selected_row_indices"] or len(fallbacks) != row["fallbacks"]:
                raise ValueError(f"selection differs from archived trial {row['trial_seed']}")
            records.append(
                {
                    "trial_seed": row["trial_seed"],
                    "observed_ids": [rows[i]["id"] for i in indices],
                    "fallback_positions_1based": fallbacks,
                }
            )
        by_n[n] = records
    return {
        "trials_sha256": hashlib.sha256(trials_path.read_bytes()).hexdigest(),
        "measurements_sha256": hashlib.sha256(measurements_path.read_bytes()).hexdigest(),
        "by_observations": by_n,
    }


def plot(root, output):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    plt.rcParams["svg.hashsalt"] = "basinmark-reuse-v1"
    plt.rcParams["svg.fonttype"] = "none"
    fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.45), sharey=True, layout="constrained")
    for ax, permitted in zip(axes, (True, False), strict=True):
        ticks = set()
        for name, marker, color, offset in GENERATORS:
            points = series(root, name, permitted)
            ticks.update(k for k, _ in points)
            x = [k + offset for k, _ in points]
            y = [v["rate"] for _, v in points]
            errors = [
                [v["rate"] - v["ci95"][0] for _, v in points],
                [v["ci95"][1] - v["rate"] for _, v in points],
            ]
            ax.errorbar(
                x,
                y,
                yerr=errors,
                label=name,
                fmt=marker,
                color=color,
                markersize=3.5,
                capsize=2,
                linewidth=1,
            )
        ax.set(xticks=sorted(ticks), xlabel="Outputs observed", ylim=(-0.06, 1.06))
        ax.set_title(
            "Exposed records allowed" if permitted else "Exposed records excluded", fontsize=10
        )
        ax.set_yticks([0, 0.25, 0.50, 0.75, 1.0])
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
        ax.grid(axis="y", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("Audits detecting reuse", fontsize=9)
    axes[1].legend(loc="center", frameon=False, fontsize=8)
    fig.savefig(output, metadata={"Date": None})
    plt.close(fig)
    svg = output.read_text(encoding="utf-8")
    output.write_text(
        "\n".join(line.rstrip() for line in svg.splitlines()) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("results/v0.3/trials"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--observations-out",
        type=Path,
        help="write IDs and retrospective fallback positions as JSON",
    )
    parser.add_argument("--measurements", type=Path, default=Path("results/v0.3/measurements"))
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix.lower() != ".svg":
        parser.error("supply a new SVG output")
    if args.observations_out:
        if args.observations_out.exists() or not args.observations_out.name.endswith(
            (".json", ".json.gz")
        ):
            parser.error("supply a new .json or .json.gz for the observation records")
        observations = {
            name: observer_records(
                args.root / f"{name}-reuse.json.gz", args.measurements / f"{name}-corpus.json.gz"
            )
            for name, *_ in GENERATORS
        }
    plot(args.root, args.out)
    if args.observations_out:
        JsonStore(args.observations_out).write(observations)
