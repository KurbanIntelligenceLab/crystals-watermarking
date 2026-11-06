"""Replay the symmetry-recovery decisions and strain diagnosis from public measurements."""

from pathlib import Path

from benchmark_symmetry_audit import measure
from diagnose_strain import diagnose

from basinmark.io import JsonStore


def main():
    root = Path(__file__).resolve().parents[1] / "results"
    new, old, past = root / "v0.5", root / "v0.4", root / "v0.3" / "measurements"

    def file(directory, stem):
        return directory / (stem + ".json.gz")

    cases = {
        "CF-dev": (old, "CrystalFormer-second"),
        "CF-first": (old, "CrystalFormer-external"),
        "CF-confirm": (new, "CrystalFormer-confirm"),
        "MG-confirm": (new, "MatterGen-confirm"),
        "MG-prior": (old, "MatterGen-second"),
        "MG-strain": (past, "processing"),
        "WyCryst": (old, "WyCryst-transfer"),
        "CF-A": (new, "CF-A"),
        "CF-B": (new, "CF-B"),
    }
    for result in sorted(list(new.glob("*-audits.json.gz")) + list(new.glob("*-ranked.json.gz"))):
        stem, kind = result.name.removesuffix(".json.gz").rsplit("-", 1)
        if stem == "CF-shear":
            corpus = file(new, "CF-shear-scored" if kind == "ranked" else "CF-shear")
            views = None
        else:
            corpus = file(*cases[stem])
            suffix = "-scores" if kind == "ranked" or stem in {"CF-A", "CF-B"} else "-views"
            views = file(new, stem + suffix)
        expected = JsonStore(result).read_object()
        actual = measure(
            corpus,
            views,
            expected["tolerance"],
            expected["seed"],
            expected["trials"],
            expected["control_seed"],
            expected["controls"],
        )
        if actual != expected:
            raise ValueError(f"recorded results differ in {result.name}")
        print("Verified", result.name, flush=True)
    diagnosis = diagnose(
        file(old, "CrystalFormer-second"),
        file(old, "CrystalFormer-second-audits"),
        file(new, "CF-dev-views"),
    )
    if diagnosis != JsonStore(file(new, "strain-diagnosis")).read_object():
        raise ValueError("recorded strain diagnosis differs")
    print("Verified design diagnosis")


if __name__ == "__main__":
    main()
