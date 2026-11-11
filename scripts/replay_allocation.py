"""Verify the published allocation audits and independent-key calibration records."""

from pathlib import Path

from benchmark_allocation import benchmark
from benchmark_frontier import evaluate as frontier
from calibrate_allocation import evaluate

from basinmark.io import JsonStore


def main():
    root = Path(__file__).resolve().parents[1] / "results"
    new, old = root / "v0.4", root / "v0.3" / "measurements"
    inputs = {
        "MatterGen-design": old / "MatterGen-512.json.gz",
        "CrystalFormer-design": old / "CrystalFormer-512.json.gz",
        "CrystalFormer-retrospective-audits": old / "CrystalFormer-384.json.gz",
        "MatterGen-strain-audits": old / "processing.json.gz",
        "CrystalFormer-first-sensitivity": new / "CrystalFormer-external.json.gz",
        "CrystalFormer-second-sensitivity": new / "CrystalFormer-second.json.gz",
        "MatterGen-second-sensitivity": new / "MatterGen-second.json.gz",
        **{
            s + "-audits": new / (s + ".json.gz")
            for s in (
                "MatterGen-fresh",
                "CrystalFormer-external",
                "WyCryst-transfer",
                "CrystalFormer-second",
                "MatterGen-second",
            )
        },
    }
    for name, corpus in inputs.items():
        expected = JsonStore(new / (name + ".json.gz")).read_object()
        actual = benchmark(
            corpus,
            expected["seed"],
            expected["trials"],
            expected.get("followup", False),
            "excluded_lowest_id" in expected,
        )
        if actual != expected:
            raise ValueError(f"mismatch in {name}")
        print("Verified", name, flush=True)
    for name, stem in [
        ("MatterGen", "MatterGen-fresh"),
        ("CrystalFormer", "CrystalFormer-external"),
        ("MatterGen-second", "MatterGen-second"),
        ("CrystalFormer-second", "CrystalFormer-second"),
    ]:
        expected = JsonStore(new / (name + "-null.json.gz")).read_object()
        actual = evaluate(
            inputs[stem + "-audits"],
            new / (stem + "-audits.json.gz"),
            expected["seed"],
            expected["keys_per_collection"],
        )
        if actual != expected:
            raise ValueError(f"mismatch in {name} null")
        print("Verified", name, "null calibration", flush=True)

    if (
        frontier(new / "CrystalFormer-external.json.gz")
        != JsonStore(new / "CrystalFormer-frontier.json.gz").read_object()
    ):
        raise ValueError("exploratory frontier differs")
    print("Verified exploratory frontier")


if __name__ == "__main__":
    main()
