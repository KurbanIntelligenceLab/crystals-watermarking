"""Compare fixed-batch and pooled watermarking on identical candidate draws and keys."""

import argparse
from pathlib import Path
from statistics import mean

import numpy as np
from benchmark_cross_model import MPA
from benchmark_producer_preview import NEW, path, sha

from basinmark.allocation import select_collection
from basinmark.collections import prepared, tally
from basinmark.core import clopper_pearson
from basinmark.io import JsonStore
from basinmark.reliability import reliable_audit
from basinmark.symmetry_audit import eligible_recovery

ARMS = ("original_full", "original_no_cell", "budget", "preview")


def candidate_draw(seed, size):
    rng = np.random.default_rng(seed)
    key = rng.bytes(16)
    rng.bytes(16)  # The original release protocol reserves a second key before drawing.
    return key, rng.permutation(size)[:256]


def endpoint_tables(directory, case, rows, scores, alternate):
    tables = {}
    tolerance = scores["tolerances"].index(0.30)
    for stage in ("rerelaxed", "strain_rerelaxed"):
        table = []
        for row, view in zip(rows, scores["rows"], strict=True):
            if row["id"] != view["id"]:
                raise ValueError("score and candidate identities differ")
            point = row["stages"].get(stage)
            table.append(
                None
                if point is None
                else (
                    point["full"],
                    point["no_cell"],
                    eligible_recovery(view["stages"][stage][tolerance]),
                    view["scores"][stage],
                )
            )
        tables["CHGNet_" + stage] = table
    if alternate is None:
        return tables
    pool = JsonStore(alternate).read_object()
    meta = pool["metadata"]
    if not (
        meta["case"] == case
        and meta["relax_weights_sha256"] == MPA
        and meta["source_sha256"] == scores["input_sha256"]
        and meta["rank_weights_sha256"] == scores["weights_sha256"]
        and (meta["fmax"], meta["steps"]) == (0.05, 150)
        and meta["indices"] == list(range(len(rows)))
    ):
        raise ValueError("alternate measurements differ from the registered pool or treatment")
    previous = JsonStore(path(directory, case + "-cross-model-physics")).read_object()
    old = {r["index"]: r["stages"]["strain_rerelaxed"] for r in previous["records"]}
    table = []
    for i, (r, p) in enumerate(zip(rows, pool["records"], strict=True)):
        if (p["index"], p["id"], p["strain_accepted"]) != (i, r["id"], r["strain_accepted"]):
            raise ValueError("alternate endpoint identity or strain guard differs")
        m = p["stages"]["strain_rerelaxed"]
        if i in old and any(m[k] != v for k, v in old[i].items()):
            raise ValueError("regenerated alternate endpoint differs from the published release")
        table.append(
            (m["full"], m["primary"], m["secondary"], m["energy"]) if m["converged"] else None
        )
    tables["MACE_strain_rerelaxed"] = table
    return tables


def compare(case, directory, alternate):
    source = path(directory, "CF-" + case)
    trials_file = path(directory, case + "-producer")
    scores_file = path(directory, "CF-" + case + "-scores")
    preview_file = path(directory, case + "-preview")
    files = {
        "source": source,
        "releases": trials_file,
        "scores": scores_file,
        "preview": preview_file,
    }
    if alternate:
        files["alternate"] = alternate
    digests = {k: sha(v) for k, v in files.items()}
    _, rows = prepared(source)
    releases = JsonStore(trials_file).read_object()
    scores = JsonStore(scores_file).read_object()
    previews = JsonStore(preview_file).read_object()
    if not (
        releases["source_sha256"]
        == previews["source_sha256"]
        == scores["input_sha256"]
        == digests["source"]
        and releases["views_sha256"] == digests["scores"]
        and releases["preview_sha256"] == digests["preview"]
        and releases["controls"] == 20
        and (releases["candidates"], releases["emissions"], len(releases["records"]))
        == (256, 32, 500)
        and [r["id"] for r in rows] == [r["id"] for r in previews["rows"]]
    ):
        raise ValueError("incompatible source, preview, or release records")
    tables = endpoint_tables(directory, case, rows, scores, alternate)
    energies = np.array([r["energy"] for r in rows])
    preview_seconds = np.array([r["seconds"] for r in previews["rows"]])
    records = []
    for trial in releases["records"]:
        key, batch = candidate_draw(trial["seed"], len(rows))
        control_rng = np.random.default_rng(trial["control_seed"])
        control_keys = [control_rng.bytes(16) for _ in range(releases["controls"])]
        trial_energy = energies[batch]
        cheapest = float(np.sort(trial_energy)[:32].mean())
        measured_seconds = float(sum(preview_seconds[batch]))
        if measured_seconds != trial["preview_worker_seconds"]:
            raise ValueError("preview timing differs from the original candidate draw")
        record = {
            "seed": trial["seed"],
            "control_seed": trial["control_seed"],
            "cheapest_mean": cheapest,
            "preview_worker_seconds": measured_seconds,
            "arms": {},
        }
        for arm in ARMS:
            original = arm.startswith("original")
            variant = "full" if arm == "original_full" else "no_cell"
            gamma = 0.5 if original else 0.625
            fps = [rows[i]["producer"][variant] for i in batch]
            views = (
                [
                    (previews["rows"][i].get("primary"), previews["rows"][i].get("secondary"))
                    for i in batch
                ]
                if arm == "preview"
                else None
            )
            chosen = select_collection(
                fps,
                trial_energy,
                key,
                gamma=gamma,
                policy="unbounded" if original else arm,
                previews=views,
            )
            ids = batch[chosen.indices].tolist()
            release_mean = float(energies[ids].mean())
            if not original:
                j = ("budget", "preview").index(arm)
                if (
                    ids != trial["selected"][j]
                    or release_mean != trial["quality"][j]["mean_energy"]
                ):
                    raise ValueError(
                        "revised selections or energies differ from the recorded trial"
                    )
            item = {
                "ids": ids,
                "release_mean": release_mean,
                "stable": sum(rows[i]["phonon_stable"] for i in ids),
                "strain_accepted": sum(rows[i]["strain_accepted"] for i in ids),
                "stages": {},
            }
            for label, table in tables.items():
                points = [table[i] for i in ids if table[i] is not None]
                primary = [p[0 if variant == "full" else 1] for p in points]
                simple = tally(primary, key, gamma, 0.001)[2]
                detection = simple
                upgraded = None
                if variant == "no_cell":
                    result = reliable_audit([p[1:] for p in points], key, gamma=gamma)
                    upgraded = int(result.fires)
                    if simple > upgraded or result.null_rate > 0.001 + 1e-14:
                        raise ValueError("nesting or audit size violated")
                if not original:
                    detection = upgraded
                    if (
                        label.startswith("CHGNet_")
                        and detection != trial["stages"][label.removeprefix("CHGNet_")][j]
                    ):
                        raise ValueError("frozen detector disagrees with archived decision")
                if arm == "original_no_cell":
                    # Both readers see the original release and each independent control key.
                    audit_inputs = [p[1:] for p in points]
                    item.setdefault("wrong_keys", {})[label] = [
                        sum(tally(primary, wrong, gamma, 0.001)[2] for wrong in control_keys),
                        sum(
                            reliable_audit(audit_inputs, wrong, gamma=gamma).fires
                            for wrong in control_keys
                        ),
                    ]
                observed_scores = [p[3] for p in points if p[3] is not None]
                item["stages"][label] = [
                    detection,
                    simple,
                    len(points),
                    len(observed_scores),
                    sum(observed_scores),
                    upgraded,
                ]
            record["arms"][arm] = item
        records.append(record)
    summary = {}
    for arm in ARMS:
        samples = [r["arms"][arm] for r in records]
        summary[arm] = {
            "release_excess_meV": 1000
            * mean(
                x["release_mean"] - r["cheapest_mean"]
                for x, r in zip(samples, records, strict=True)
            ),
            "strain_accepted_fraction": sum(x["strain_accepted"] for x in samples)
            / (32 * len(samples)),
            "stable_fraction": sum(x["stable"] for x in samples) / (32 * len(samples)),
            "preview_worker_seconds_per_emission": mean(
                r["preview_worker_seconds"] for r in records
            )
            / 32
            if arm == "preview"
            else 0,
            "stages": {},
        }
        for label in tables:
            decisions, simple, covered, scored, total, upgraded = zip(
                *(x["stages"][label] for x in samples)
            )
            reference = [r["arms"]["original_no_cell"]["stages"][label][0] for r in records]
            gains = sum(a > b for a, b in zip(decisions, reference, strict=True))
            losses = sum(a < b for a, b in zip(decisions, reference, strict=True))
            gain_ci = clopper_pearson(gains, len(records), alpha=0.025)
            loss_ci = clopper_pearson(losses, len(records), alpha=0.025)
            summary[arm]["stages"][label] = {
                "detected": sum(decisions),
                "primary_only": sum(simple),
                "coverage": sum(covered),
                "scored": sum(scored),
                "mean_deposited_meV": 1000 * sum(total) / sum(scored),
                "gained_vs_no_cell": gains,
                "lost_vs_no_cell": losses,
                "paired_delta_ci95": [gain_ci[0] - loss_ci[1], gain_ci[1] - loss_ci[0]],
                "upgraded_on_same_outputs": sum(upgraded) if arm != "original_full" else None,
            }
            if arm == "original_no_cell":
                paired_controls = [x["wrong_keys"][label] for x in samples]
                summary[arm]["stages"][label]["wrong_keys"] = {
                    "primary": sum(x[0] for x in paired_controls),
                    "upgraded": sum(x[1] for x in paired_controls),
                    "trials": releases["controls"] * len(records),
                }
    return {
        "case": case,
        "input_sha256": digests,
        "candidates": 256,
        "outputs": 32,
        "trials": len(records),
        "alpha": 0.001,
        "summary": summary,
        "records": records,
        "stage_columns": [
            "detected",
            "primary_only_at_producer_gamma",
            "covered",
            "scored",
            "sum_deposited_eV_per_atom",
            "upgraded_at_producer_gamma",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=NEW)
    parser.add_argument("--case", choices=["E", "F"], required=True)
    parser.add_argument("--alternate", type=Path, help="all-candidate strained MACE measurements")
    parser.add_argument(
        "--out", type=Path, help="write a new result instead of checking the archive"
    )
    args = parser.parse_args()
    if args.out and args.out.exists():
        parser.error("choose a new output")
    measurements = args.alternate
    if not args.out and measurements is None:
        measurements = path(args.directory, args.case + "-strain-pool")
    comparison = compare(args.case, args.directory, measurements)
    if args.out:
        JsonStore(args.out).write(comparison)
    elif (
        comparison
        != JsonStore(path(args.directory, args.case + "-method-comparison")).read_object()
    ):
        raise ValueError("method comparison differs from archived trials")
    for policy, result in comparison["summary"].items():
        print(args.case, policy, result, flush=True)
