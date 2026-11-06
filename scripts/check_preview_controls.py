"""Repeat the wrong-key audit on every prospective strain-treated preview release."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from benchmark_producer_preview import NEW, path, sha

from basinmark.collections import prepared
from basinmark.io import JsonStore
from basinmark.reliability import reliable_audit
from basinmark.symmetry_audit import eligible_recovery

SEEDS = {"E": 774530019, "F": 306417957}
KEYS_PER_RELEASE = 200


def check_cohort(case, directory, cross_model=False):
    stage = "strain_rerelaxed"
    source = path(directory, "CF-" + case)
    views = path(directory, "CF-" + case + "-scores")
    releases = path(directory, case + "-producer")
    trial_data = JsonStore(releases).read_object()
    measurements = JsonStore(views).read_object()
    if (
        trial_data["source_sha256"] != sha(source)
        or trial_data["views_sha256"] != sha(views)
        or measurements["input_sha256"] != sha(source)
    ):
        raise ValueError("prospective input hashes differ")
    _, rows = prepared(source)
    secondary = measurements["tolerances"].index(0.30)
    observations = []
    for r, m in zip(rows, measurements["rows"], strict=True):
        if r["id"] != m["id"]:
            raise ValueError("different physical records")
        obs = None
        if stage in r["stages"]:
            obs = (
                r["stages"][stage]["no_cell"],
                eligible_recovery(m["stages"][stage][secondary]),
                m["scores"][stage],
            )
        observations.append(obs)

    expected = [trial["stages"][stage][1] for trial in trial_data["records"]]
    extra = {}
    if cross_model:
        physics = path(directory, case + "-cross-model-physics")
        audit_file = path(directory, case + "-cross-model")
        alternate = JsonStore(audit_file).read_object()
        if alternate["physics_sha256"] != sha(physics) or alternate["producer_sha256"] != sha(
            releases
        ):
            raise ValueError("alternate physical records or releases differ")
        observations = [None] * len(rows)
        for record in JsonStore(physics).read_object()["records"]:
            i = record["index"]
            if rows[i]["id"] != record["id"]:
                raise ValueError("alternate structure identity differs")
            ep = record["stages"][stage]
            if ep["converged"]:
                observations[i] = (ep["primary"], ep["secondary"], ep["energy"])
        expected = []
        for alt, trial in zip(alternate["records"], trial_data["records"], strict=True):
            if alt["seed"] != trial["seed"]:
                raise ValueError("alternate audit trial order differs")
            expected.append(alt["results"][stage]["preview"][1])
        extra = {"cross_model_audit_sha256": sha(audit_file)}

    seeds = np.random.default_rng(SEEDS[case]).integers(0, 2**53, len(trial_data["records"]))
    records = []
    for release, control_seed, decision in zip(trial_data["records"], seeds, expected, strict=True):
        deposited = [observations[i] for i in release["selected"][1] if observations[i] is not None]
        marking_key = np.random.default_rng(release["seed"]).bytes(16)
        marked = reliable_audit(deposited, marking_key)
        if int(marked.fires) != decision or marked.null_rate > 0.001 + 1e-14:
            raise ValueError("true decision or conditional audit size differs")
        rng = np.random.default_rng(control_seed)
        hits = sum(reliable_audit(deposited, rng.bytes(16)).fires for _ in range(KEYS_PER_RELEASE))
        records.append([release["seed"], int(control_seed), hits, marked.null_rate])
    return {
        **extra,
        "cohort": case,
        "releases_sha256": sha(releases),
        "seed": SEEDS[case],
        "trial_columns": ["release_seed", "wrong_key_seed", "wrong_decisions", "conditional_size"],
        "records": records,
        "wrong_decisions": sum(row[2] for row in records),
        "wrong_keys": KEYS_PER_RELEASE * len(records),
        "expected_wrong": KEYS_PER_RELEASE * sum(row[3] for row in records),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=NEW)
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--cross-model", action="store_true", help="use matched MACE-MPA-0 endpoints"
    )
    args = parser.parse_args()
    target = path(
        args.directory, "cross-model-controls" if args.cross_model else "preview-controls"
    )
    if args.write and target.exists():
        parser.error("output already exists")
    with ProcessPoolExecutor(max_workers=len(SEEDS)) as pool:
        runs = list(
            pool.map(
                check_cohort, SEEDS, [args.directory] * len(SEEDS), [args.cross_model] * len(SEEDS)
            )
        )
    result = {"keys_per_release": KEYS_PER_RELEASE, "cohorts": runs}
    if args.write:
        JsonStore(target).write(result)
    elif result != JsonStore(target).read_object():
        raise ValueError("wrong-key replay differs")
    for run in runs:
        print(run["cohort"], run["wrong_decisions"], "/", run["wrong_keys"], flush=True)
