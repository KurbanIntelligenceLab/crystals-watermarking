"""Recompute statistical trial records from the bundled immutable measurements."""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

from basinmark.analyses import main as analyse
from basinmark.io import JsonStore


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="include development windows and stress")
    args = parser.parse_args()
    archive = Path(__file__).resolve().parents[1] / "results" / "v0.3"
    jobs = []
    for name in ("MatterGen", "CrystalFormer"):
        jobs.append((f"{name}-confirm", f"{name}-384", "collection-benchmark"))
        jobs.append((f"{name}-reuse", f"{name}-corpus", "class-reuse"))
        if args.all:
            jobs.append((f"{name}-development", f"{name}-512", "collection-benchmark"))
    jobs.append(("WyCryst-controls", "WyCryst-corpus", "key-controls"))
    if args.all:
        for window in ("0.10", "0.15"):
            jobs.append((f"MatterGen-window-{window}", "MatterGen-512", "collection-benchmark"))
        jobs.append(("strain", "processing", "collection-benchmark"))
    for label, inputs, command in jobs:
        expected = JsonStore(archive / "trials" / f"{label}.json.gz").read_object()
        corpus = archive / "measurements" / f"{inputs}.json.gz"
        base = [
            command,
            "--corpus",
            str(corpus),
            "--seed",
            str(expected["seed"]),
            "--m",
            str(expected["m"]),
            "--gamma",
            str(expected["gamma"]),
            "--alpha",
            str(expected["alpha"]),
        ]
        if command == "collection-benchmark":
            base += [
                "--sizes",
                *expected["summary"],
                "--energy-window",
                str(expected["energy_window_eV_per_atom"]),
                "--stability-threshold",
                str(expected["stability_threshold_eV_per_atom"]),
                "--trials",
                str(len(expected["trial_seeds"])),
            ]
        elif command == "key-controls":
            base += [
                "--keys",
                str(len(expected["trial_records"])),
                "--representation",
                expected["representation"],
            ]
        else:
            first = next(iter(expected["rows"].values()))
            base += [
                "--observations",
                *expected["rows"],
                "--trials",
                str(first["trials"]),
                "--collection",
                str(expected["collection_size"]),
                "--representation",
                expected["representation"],
            ]
        with TemporaryDirectory() as temp:
            output = Path(temp) / "result.json"
            analyse([*base, "--out", str(output)])
            if JsonStore(output).read_object() != expected:
                raise ValueError(f"Replay differs for {label}")
        print(f"Verified {label}")


if __name__ == "__main__":
    main()
