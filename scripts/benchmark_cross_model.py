"""Replay paired collection audits after changing the downstream relaxation model."""

import argparse
from pathlib import Path

import numpy as np
from benchmark_producer_preview import NEW, path, sha

from basinmark.collections import prepared, rate
from basinmark.io import JsonStore
from basinmark.reliability import reliable_audit

MPA = "75428afe3a1d7d8062e19bcaabd5c433623cabf308242ec9fb493e38604fb638"


def benchmark(case, directory, physics):
    source = path(directory, "CF-" + case)
    producer = path(directory, case + "-producer")
    scores = path(directory, "CF-" + case + "-scores")
    source_meta, rows = prepared(source)
    releases = JsonStore(producer).read_object()
    views = JsonStore(scores).read_object()
    measured = JsonStore(physics).read_object()
    meta = measured["metadata"]
    if not (
        meta["case"] == case
        and meta["archive_sha256"] == source_meta["source_sha256"]
        and meta["rank_weights_sha256"] == views["weights_sha256"]
        and meta["relax_weights_sha256"] == MPA
        and (meta["fmax"], meta["steps"]) == (0.05, 150)
        and sha(source) == meta["source_sha256"] == releases["source_sha256"]
        and sha(producer) == meta["producer_sha256"]
        and sha(scores) == releases["views_sha256"]
        and views["input_sha256"] == meta["source_sha256"]
    ):
        raise ValueError("incompatible measurements and archived releases")
    indices = sorted({i for t in releases["records"] for arm in t["selected"] for i in arm})
    if [r["index"] for r in measured["records"]] != indices or indices != meta["indices"]:
        raise ValueError("coverage differs from the union of emitted records")
    endpoints = {}
    for r in measured["records"]:
        i = r["index"]
        if r["id"] != rows[i]["id"] or r["id"] != views["rows"][i]["id"]:
            raise ValueError("structure identity differs")
        endpoints[i] = r["stages"]

    records, summary = [], {}
    for trial in releases["records"]:
        key = np.random.default_rng(trial["seed"]).bytes(16)
        control_rng = np.random.default_rng(trial["control_seed"])
        controls = [control_rng.bytes(16) for _ in range(releases["controls"])]
        record = {"seed": trial["seed"], "control_seed": trial["control_seed"], "results": {}}
        for stage in ("rerelaxed", "strain_rerelaxed"):
            record["results"][stage] = {}
            for a, (policy, selected) in enumerate(
                zip(("budget", "preview"), trial["selected"], strict=True)
            ):
                obs = []
                delta_sum = 0.0
                paired = valid = 0
                for i in selected:
                    m = endpoints[i][stage]
                    if not m["converged"]:
                        continue
                    obs.append((m["primary"], m["secondary"], m["energy"]))
                    valid += m["valid"]
                    old_energy = views["rows"][i]["scores"].get(stage)
                    if old_energy is not None:
                        delta_sum += m["energy"] - old_energy
                        paired += 1
                current = reliable_audit(obs, key)
                if current.null_rate > 0.001 + 1e-14:
                    raise ValueError("conditional size exceeds the audit level")
                old = trial["stages"][stage][a]
                passed = int(current.fires)
                wrong = sum(reliable_audit(obs, k).fires for k in controls)
                record["results"][stage][policy] = [
                    old,
                    passed,
                    trial["stages"][stage][a + 2],
                    wrong,
                    len(obs),
                    valid,
                    paired,
                    delta_sum,
                ]
        records.append(record)
    for stage in ("rerelaxed", "strain_rerelaxed"):
        summary[stage] = {}
        for policy in ("budget", "preview"):
            arms = [t["results"][stage][policy] for t in records]
            prev, now, old_wrong, wrong, converged, valid, paired, delta = map(
                sum, zip(*arms, strict=True)
            )
            summary[stage][policy] = {
                "previous": prev,
                "alternate": rate(now, len(arms)),
                "gained": sum(o[0] < o[1] for o in arms),
                "lost": sum(o[0] > o[1] for o in arms),
                "previous_wrong": old_wrong,
                "alternate_wrong": wrong,
                "wrong_tests": releases["controls"] * len(arms),
                "attempted_emissions": len(arms) * releases["emissions"],
                "converged_emissions": converged,
                "valid_converged_emissions": valid,
                "jointly_scored_emissions": paired,
                "mean_paired_energy_change_meV_per_atom": 1000 * delta / paired if paired else None,
            }
    return {
        "case": case,
        "physics_sha256": sha(physics),
        "producer_sha256": sha(producer),
        "scores_sha256": sha(scores),
        "controls_per_trial": releases["controls"],
        "unique_emitted": len(indices),
        "summary": summary,
        "stage_columns": [
            "original_fires",
            "alternate_fires",
            "original_wrong",
            "alternate_wrong",
            "converged",
            "valid",
            "jointly_scored",
            "sum_energy_change_eV_per_atom",
        ],
        "records": records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["E", "F"], required=True)
    parser.add_argument("--directory", type=Path, default=NEW)
    parser.add_argument(
        "--physics", type=Path, help="defaults to the archived measurement for the selected cohort"
    )
    parser.add_argument(
        "--out", type=Path, help="write a new result instead of verifying the archive"
    )
    args = parser.parse_args()
    if args.out and args.out.exists():
        parser.error("choose a new output")
    result = benchmark(
        args.case,
        args.directory,
        args.physics or path(args.directory, args.case + "-cross-model-physics"),
    )
    if args.out:
        JsonStore(args.out).write(result)
    elif result != JsonStore(path(args.directory, args.case + "-cross-model")).read_object():
        raise ValueError("cross-model replay differs")
    print(
        args.case,
        {
            s: {p: v["alternate"]["rejections"] for p, v in arms.items()}
            for s, arms in result["summary"].items()
        },
    )
