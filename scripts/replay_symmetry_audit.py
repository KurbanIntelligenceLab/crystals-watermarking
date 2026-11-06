"""Replay the symmetry-recovery decisions and strain diagnosis from public measurements."""

from pathlib import Path

from benchmark_symmetry_audit import measure
from diagnose_strain import diagnose
from split_recovery_cohort import split

from basinmark.io import JsonStore


def main():
    root = Path(__file__).resolve().parents[1] / "results"
    new, old, past = root / "v0.5", root / "v0.4", root / "v0.3" / "measurements"

    def file(directory, stem):
        return directory / (stem + ".json.gz")

    halves = split(file(new, "CF-reserved"))
    for name, expected in zip(("CF-A", "CF-B"), halves, strict=True):
        if expected != JsonStore(file(new, name)).read_object():
            raise ValueError(f"reserved split differs for {name}")

    def identities(path):
        rows = JsonStore(path).read_object()["rows"]
        return {r["id"] for r in rows}, {r["producer"]["no_cell"] for r in rows if "producer" in r}

    old_ids, old_classes = set(), set()
    for path in [
        file(old, "CrystalFormer-external"),
        file(old, "CrystalFormer-second"),
        file(new, "CrystalFormer-confirm"),
    ]:
        ids, classes = identities(path)
        old_ids.update(ids)
        old_classes.update(classes)
    old_classes.update(identities(file(past, "CrystalFormer-corpus"))[1])
    for name in ["CF-A", "CF-B"]:
        ids, classes = identities(file(new, name))
        if ids & old_ids or classes & old_classes:
            raise ValueError("reserved source IDs or producer classes intersect earlier data")
        old_ids.update(ids)
        old_classes.update(classes)
    print("Verified reserved source identities and class exclusions", flush=True)

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
    for cohort, stem in [("confirm", "CrystalFormer-confirm"), ("A", "CF-A"), ("B", "CF-B")]:
        for m in [6, 8]:
            path = file(new, f"size-{cohort}-{m}")
            expected = JsonStore(path).read_object()
            actual = measure(
                file(new, stem),
                file(new, f"CF-{cohort}-scores"),
                expected["tolerance"],
                expected["seed"],
                expected["trials"],
                expected["control_seed"],
                expected["controls"],
                40,
                m,
                True,
            )
            if actual != expected:
                raise ValueError(f"release-size comparison differs in {path.name}")
            print("Verified", path.name, flush=True)
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
