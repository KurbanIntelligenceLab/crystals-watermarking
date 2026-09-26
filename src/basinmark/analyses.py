"""Statistical and structure-based analysis commands."""

from __future__ import annotations

import argparse
import hashlib
import math
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from . import __version__
from .classification import ClassificationConfig, FlipClassifier
from .core import (
    audit,
    audit_size,
    clopper_pearson,
    design_effect,
    fingerprint_structure,
    holm_bonferroni,
    is_green,
    n_eff_debiased,
    n_forge,
    q_green,
    select_green,
    self_test,
)
from .io import FingerprintCorpus, JsonStore, PairManifest
from .relaxation import CalculatorRegistry, RelaxationEngine

GAMMA, ALPHA, M = 0.5, 1e-3, 8


def save(obj, path):
    JsonStore(Path(path)).write(obj)


def load_pairs(path, fields):
    return PairManifest(Path(path)).load(fields)


def load_structures(paths):
    """Read CIF/POSCAR files into pymatgen Structures."""
    from pymatgen.core import Structure

    out = []
    for p in paths:
        try:
            out.append((p, Structure.from_file(p)))
        except Exception as exc:  # noqa: BLE001
            print(f"  skipped {p}: {exc}")
    return out


_CALCULATORS = CalculatorRegistry()


def relax(structure, engine: str, fmax: float = 0.02, steps: int = 500):
    result = RelaxationEngine(_CALCULATORS, engine, fmax, steps).run(structure)
    return result.structure, result.energy_per_atom, result.steps


def rmsd_to_reference(a, b):
    """Geometry mismatch between two relaxations of one structure, the quantity
    Matbench Discovery reports for geometry optimisation."""
    from pymatgen.analysis.structure_matcher import StructureMatcher

    sm = StructureMatcher(primitive_cell=True, attempt_supercell=False)
    try:
        val = sm.get_rms_dist(a, b)
        return None if val is None else float(val[0])
    except Exception:  # noqa: BLE001
        return None


def fp(structure, **kw):
    return fingerprint_structure(structure, **kw)


def generator_controls(a):
    rng = np.random.default_rng(a.seed)
    out = {}
    for gen_dir in a.gen_dirs:
        directory = Path(gen_dir)
        name = directory.name
        files = sorted(
            path for path in directory.iterdir() if path.suffix.lower() in {".cif", ".vasp"}
        )
        structs = load_structures(files)
        rows, energies = [], []
        for _, s in structs:
            rs, e, _ = relax(s, a.engine, fmax=a.fmax)
            rows.append(fp(rs, symprec=a.symprec, cell_bin=a.cell_bin))
            energies.append(e)
        slots = [list(range(i, min(i + a.m, len(rows)))) for i in range(0, len(rows), a.m)]
        slots = [s for s in slots if len(s) == a.m]

        if not slots:
            raise ValueError(f"{gen_dir}: fewer than {a.m} usable structures")
        marked, unmarked, shuffled = [], [], []
        green_fractions, unmarked_pvalues = [], []
        for _ in range(a.keys):
            key = rng.bytes(16)
            selected = [
                select_green([rows[i] for i in slot], [energies[i] for i in slot], key, a.gamma)
                for slot in slots
            ]
            emitted = [rows[slot[index]] for slot, (index, _) in zip(slots, selected, strict=True)]
            green_fractions.append(float(np.mean([green for _, green in selected])))
            marked.append(audit(emitted, key, a.gamma, a.alpha).fires)
            control = audit([rows[slot[0]] for slot in slots], key, a.gamma, a.alpha)
            unmarked.append(control.fires)
            unmarked_pvalues.append(control.p_value)
            shuffled.append(audit(emitted, rng.bytes(16), a.gamma, a.alpha).fires)

        gf = float(np.mean(green_fractions))
        out[name] = dict(
            slots=len(slots),
            green_fraction_marked=gf,
            detect_marked=float(np.mean(marked)),
            detect_unmarked=float(np.mean(unmarked)),
            detect_shuffled=float(np.mean(shuffled)),
            q_predicted=q_green(a.gamma, a.m),
            holm_survives_unmarked=bool(holm_bonferroni(unmarked_pvalues, 0.05).any()),
        )
        print(
            f"  {name}: marked {out[name]['detect_marked']:.3f}  "
            f"unmarked {out[name]['detect_unmarked']:.4f}"
        )
    save(out, a.out or "results/generator_controls.json")


def relaxation_stability(a):
    pairs = load_pairs(a.pairs, ("mlip", "dft"))
    rng = np.random.default_rng(a.seed)
    by_system, flips = defaultdict(list), []
    for rec in pairs:
        from pymatgen.core import Structure

        sm, sd = Structure.from_file(rec["mlip"]), Structure.from_file(rec["dft"])
        if sm.composition.reduced_formula != sd.composition.reduced_formula:
            continue  # pairing check, see Appendix
        f1 = fp(sm, symprec=a.symprec, cell_bin=a.cell_bin)
        f2 = fp(sd, symprec=a.symprec, cell_bin=a.cell_bin)
        flipped = int(f1 != f2)
        flips.append(flipped)
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        by_system[SpacegroupAnalyzer(sd, symprec=a.symprec).get_crystal_system()].append(flipped)

    def block(v):
        v = np.asarray(v)
        p = float(v.mean())
        boot = [float(np.mean(rng.choice(v, v.size, replace=True))) for _ in range(5000)]
        return dict(
            pairs=int(v.size),
            p=p,
            ci95=[float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
            N_needed=audit_size(a.gamma, a.m, p, alpha=a.alpha),
        )

    if not flips:
        raise ValueError("No composition-matched pairs are available")
    save(
        dict(overall=block(flips), by_system={k: block(v) for k, v in by_system.items()}),
        a.out or "results/relaxation_stability.json",
    )


def potential_comparison(a):
    """Measure per-potential fingerprint flips and matched geometry distances."""
    recs = load_pairs(a.pairs, ("raw", "dft"))
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    # Include contents, not just manifest paths, in the resume identity.
    structure_digest = hashlib.sha256()
    for rec in recs:
        for field in ("raw", "dft"):
            structure_digest.update(hashlib.sha256(Path(rec[field]).read_bytes()).digest())
    out_path = Path(a.out or "results/potential_comparison.json")
    out = {}
    if out_path.exists():
        out = JsonStore(Path(out_path)).read_object()
    for engine in a.engines:
        signature = dict(
            engine=engine,
            pairs_sha256=hashlib.sha256(Path(a.pairs).read_bytes()).hexdigest(),
            structures_sha256=structure_digest.hexdigest(),
            steps=500,
            gamma=a.gamma,
            alpha=a.alpha,
            m=a.m,
            symprec=a.symprec,
            cell_bin=a.cell_bin,
            fmax=a.fmax,
        )
        if engine in out:
            if not isinstance(out[engine], dict) or "per_pair" not in out[engine]:
                raise ValueError("Cannot resume an aggregate summary; choose a new output path")
            if out[engine].get("_configuration") != signature:
                raise ValueError(
                    f"{engine}: checkpoint inputs or settings differ; choose a new output"
                )
            print(f"  {engine}: matching completed run in {out_path}, skipping")
            continue
        # Resume completed pair records from the last saved checkpoint.
        partial_key = engine + "__partial"
        flips, rmsds, sysflip, skipped_mismatch = [], [], defaultdict(list), 0
        per_pair = []
        resume_from = 0
        if partial_key in out:
            prev = out[partial_key]
            if not isinstance(prev, dict):
                raise ValueError(f"{engine}: malformed checkpoint")
            if prev.get("_configuration") != signature:
                raise ValueError(
                    f"{engine}: checkpoint inputs or settings differ; choose a new output"
                )
            per_pair = prev.get("per_pair", [])
            resume_from = prev.get("done", len(per_pair))
            if not isinstance(per_pair, list) or not isinstance(resume_from, int):
                raise ValueError(f"{engine}: malformed checkpoint")
            if resume_from < 0 or resume_from > len(recs):
                raise ValueError(f"{engine}: checkpoint position is outside the pair manifest")
            skipped_mismatch = prev.get("skipped_mismatch", 0)
            for row in per_pair:
                required = {"flipped", "rmsd", "crystal_system"}
                if not isinstance(row, dict) or not required.issubset(row):
                    raise ValueError(f"{engine}: checkpoint contains an incomplete pair record")
                flips.append(row["flipped"])
                if row["rmsd"] is not None:
                    rmsds.append(row["rmsd"])
                sysflip[row["crystal_system"]].append(row["flipped"])
            print(
                f"  {engine}: resuming from checkpoint at {resume_from}/{len(recs)} "
                f"({len(per_pair)} pairs already scored)"
            )
        t_engine0 = time.time()
        for i, rec in enumerate(recs):
            if i < resume_from:
                continue
            raw, dft = Structure.from_file(rec["raw"]), Structure.from_file(rec["dft"])
            # Reject composition-mismatched pairs before relaxation.
            if raw.composition.reduced_formula != dft.composition.reduced_formula:
                skipped_mismatch += 1
                continue
            try:
                rs, _, _ = relax(raw, engine, fmax=a.fmax)
            except Exception as exc:  # noqa: BLE001
                print(f"  {engine}: relaxation failed on {rec.get('pair_id')}, {exc}")
                continue
            f_mlip = fp(rs, symprec=a.symprec, cell_bin=a.cell_bin)
            f_dft = fp(dft, symprec=a.symprec, cell_bin=a.cell_bin)
            flipped = int(f_mlip != f_dft)
            flips.append(flipped)
            r = rmsd_to_reference(rs, dft)
            if r is not None:
                rmsds.append(r)
            csys = SpacegroupAnalyzer(dft, symprec=a.symprec).get_crystal_system()
            sysflip[csys].append(flipped)
            per_pair.append(
                dict(
                    pair_id=rec.get("pair_id"),
                    generator=rec.get("generator"),
                    crystal_system=csys,
                    flipped=flipped,
                    rmsd=r,
                )
            )
            if (i + 1) % 100 == 0:
                elapsed = time.time() - t_engine0
                print(
                    f"  {engine}: {i + 1}/{len(recs)}  p_so_far={np.mean(flips):.3f}  "
                    f"elapsed={elapsed:.0f}s  mean_s/pair={elapsed / (i + 1):.2f}"
                )
                out[partial_key] = dict(
                    n=len(flips),
                    p=float(np.mean(flips)) if flips else None,
                    done=i + 1,
                    total=len(recs),
                    per_pair=per_pair,
                    skipped_mismatch=skipped_mismatch,
                    _configuration=signature,
                )
                save(out, out_path)
        if not flips:
            raise ValueError(f"{engine}: no composition-matched pairs completed successfully")
        p = float(np.mean(flips))
        out.pop(engine + "__partial", None)
        out[engine] = dict(
            n=len(flips),
            p=p,
            skipped_formula_mismatch=skipped_mismatch,
            median_rmsd=float(np.median(rmsds)) if rmsds else None,
            n_rmsd=len(rmsds),
            N_needed=audit_size(a.gamma, a.m, p, alpha=a.alpha) if flips else None,
            N_needed_basis="observed_pair_pool",
            by_system={k: dict(n=len(v), p=float(np.mean(v))) for k, v in sysflip.items()},
            per_pair=per_pair,
            wall_s=time.time() - t_engine0,
            _configuration=signature,
        )
        print(
            f"  {engine:12s} p={p:.3f}  median RMSD={out[engine]['median_rmsd']}  "
            f"wall={out[engine]['wall_s']:.0f}s"
        )
        save(out, out_path)  # checkpoint after each engine
    save(out, out_path)


def class_reuse(a):
    """Simulate reuse of observed green fingerprint classes."""
    pool = FingerprintCorpus(Path(a.fingerprints)).load()
    rng = np.random.default_rng(a.seed)
    out = {"N_forge_predicted": n_forge(a.gamma, a.alpha), "rows": {}}
    for N in a.observations:
        hits = []
        for _ in range(a.trials):
            key = rng.bytes(16)
            observed = [
                f for f in rng.choice(pool, N * 4, replace=False) if is_green(f, key, a.gamma)
            ][:N]
            if len(observed) < N:
                hits.append(False)
                continue
            forged = list(rng.choice(observed, a.collection, replace=True))
            hits.append(audit(forged, key, a.gamma, a.alpha).fires)
        k = int(np.sum(hits))
        lo, hi = clopper_pearson(k, len(hits))
        out["rows"][N] = dict(detected=k / len(hits), ci95=[lo, hi], trials=len(hits))
        print(f"  N={N:6d}  forged detected {k / len(hits):.3f}  [{lo:.3f}, {hi:.3f}]")
    save(out, a.out or "results/class_reuse.json")


def generator_calibration(a):
    rng = np.random.default_rng(a.seed)
    out = {}
    for path in a.fingerprint_files:
        name = Path(path).stem
        fps = FingerprintCorpus(Path(path)).load()
        counts = list(Counter(fps).values())
        naive = dedup = 0
        for _ in range(a.keys):
            key = rng.bytes(16)
            naive += audit(fps, key, a.gamma, a.alpha, deduplicate=False).fires
            dedup += audit(fps, key, a.gamma, a.alpha, deduplicate=True).fires
        out[name] = dict(
            n=len(fps),
            distinct=len(counts),
            n_eff=n_eff_debiased(counts),
            design_effect=design_effect(counts),
            fpr_uncorrected=naive / a.keys,
            fpr_deduplicated=dedup / a.keys,
            ci_uncorrected=list(clopper_pearson(naive, a.keys)),
            ci_deduplicated=list(clopper_pearson(dedup, a.keys)),
        )
        print(
            f"  {name:28s} n_eff {out[name]['n_eff']:9.1f}  "
            f"uncorrected FPR {out[name]['fpr_uncorrected']:.4f}"
        )
    save(out, a.out or "results/generator_calibration.json")


def fingerprint_ablation(a):
    from pymatgen.core import Structure

    recs = load_pairs(a.pairs, ("raw", "relaxed"))
    out = {}
    settings = [("no cell shape", dict(use_cell_shape=False))] + [
        (f"cell bin {b}", dict(cell_bin=b)) for b in a.cell_bins
    ]
    for label, kw in settings:
        flips, fps = [], []
        for rec in recs:
            raw, rel = Structure.from_file(rec["raw"]), Structure.from_file(rec["relaxed"])
            f1 = fp(raw, symprec=a.symprec, **kw)
            f2 = fp(rel, symprec=a.symprec, **kw)
            flips.append(int(f1 != f2))
            fps.append(f2)
        counts = list(Counter(fps).values())
        p = float(np.mean(flips))
        out[label] = dict(
            p=p,
            n_eff=n_eff_debiased(counts),
            distinct=len(counts),
            N_needed=audit_size(a.gamma, a.m, p, alpha=a.alpha),
        )
        print(
            f"  {label:18s} p={p:.3f}  n_eff={out[label]['n_eff']:.0f}  N={out[label]['N_needed']}"
        )
    save(out, a.out or "results/fingerprint_ablation.json")


def cross_database_matching(a):
    ext = FingerprintCorpus(Path(a.external)).load()
    marked = FingerprintCorpus(Path(a.marked)).load()
    rng = np.random.default_rng(a.seed)
    fired = sum(audit(ext, rng.bytes(16), a.gamma, a.alpha).fires for _ in range(a.keys))
    key = a.key.encode()
    r_ext = audit(ext, key, a.gamma, a.alpha)
    r_mark = audit(marked, key, a.gamma, a.alpha)
    marked_set = set(marked)
    shared = sorted(set(ext) & marked_set)
    colliding = [f for f in ext if f in marked_set]
    r_coll = audit(colliding, key, a.gamma, a.alpha) if colliding else None
    save(
        dict(
            n_external=len(ext),
            keys=a.keys,
            fired_random_keys=fired,
            external_under_marking_key=dict(
                green=r_ext.green_fraction, p=r_ext.p_value, fires=r_ext.fires
            ),
            marked=dict(green=r_mark.green_fraction, p=r_mark.p_value, fires=r_mark.fires),
            shared_fingerprints=len(shared),
            shared_fraction=len(shared) / max(1, len(set(ext))),
            colliding_subset=None
            if r_coll is None
            else dict(
                n=len(colliding), green=r_coll.green_fraction, p=r_coll.p_value, fires=r_coll.fires
            ),
            _warning=(
                "The colliding subset is green by construction and fires "
                "against a provider that generated none of it.  An audit "
                "must be run on a sample the auditor did not select on "
                "fingerprint."
            ),
        ),
        a.out or "results/cross_database_matching.json",
    )


def classify_flips(a):
    config = ClassificationConfig(
        engine=a.engine,
        potential_full=a.potential_full,
        pairs_meta=a.pairs_meta,
        cifs_dir=a.cifs_dir,
        output=a.out,
        symprec=a.symprec,
        cell_bin=a.cell_bin,
        fmax=a.fmax,
    )
    FlipClassifier(config).run()


def run_self_test(_):
    if self_test():
        raise RuntimeError("self-test failed")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="basinmark", description="Analyses supporting the manuscript and supplement."
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="analysis", required=True)

    def common(parser):
        parser.add_argument(
            "--out",
            required=True,
            help="output JSON path; choose a new file to preserve reported results",
        )
        parser.add_argument("--gamma", type=float, default=GAMMA)
        parser.add_argument("--alpha", type=float, default=ALPHA)

    p = sub.add_parser("generator-controls", help="run marked and unmarked generator controls")
    common(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--m", type=int, default=M)
    p.add_argument("--symprec", type=float, default=0.1)
    p.add_argument("--cell-bin", type=float, default=0.05)
    p.add_argument("--engine", default="chgnet")
    p.add_argument("--fmax", type=float, default=0.02)
    p.add_argument("--gen-dirs", nargs="+", required=True)
    p.add_argument("--keys", type=int, default=200)
    p.set_defaults(func=generator_controls)

    p = sub.add_parser("relaxation-stability", help="compare relaxed structures with references")
    common(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--m", type=int, default=M)
    p.add_argument("--symprec", type=float, default=0.1)
    p.add_argument("--cell-bin", type=float, default=0.05)
    p.add_argument("--pairs", required=True)
    p.set_defaults(func=relaxation_stability)

    p = sub.add_parser("potential-comparison", help="compare relaxation potentials")
    common(p)
    p.add_argument("--m", type=int, default=M)
    p.add_argument("--symprec", type=float, default=0.1)
    p.add_argument("--cell-bin", type=float, default=0.05)
    p.add_argument("--pairs", required=True)
    p.add_argument("--fmax", type=float, default=0.02)
    p.add_argument(
        "--engines", nargs="+", default=["chgnet", "mace-mpa", "orb-v3", "mattersim", "sevennet"]
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="resume the same inputs and settings from an existing checkpoint",
    )
    p.set_defaults(func=potential_comparison)

    p = sub.add_parser("class-reuse", help="simulate observed-class reuse")
    common(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--fingerprints", required=True)
    p.add_argument("--observations", type=int, nargs="+", default=[5, 7, 9, 10, 11, 15, 100])
    p.add_argument("--trials", type=int, default=200)
    p.add_argument("--collection", type=int, default=200)
    p.set_defaults(func=class_reuse)

    p = sub.add_parser("generator-calibration", help="measure naive and deduplicated calibration")
    common(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--fingerprint-files", nargs="+", required=True)
    p.add_argument("--keys", type=int, default=1000)
    p.set_defaults(func=generator_calibration)

    p = sub.add_parser("fingerprint-ablation", help="compare fingerprint representations")
    common(p)
    p.add_argument("--m", type=int, default=M)
    p.add_argument("--symprec", type=float, default=0.1)
    p.add_argument("--pairs", required=True)
    p.add_argument("--cell-bins", type=float, nargs="+", default=[0.30, 0.10, 0.05, 0.02])
    p.set_defaults(func=fingerprint_ablation)

    p = sub.add_parser("cross-database-matching", help="evaluate shared fingerprint classes")
    common(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--key", default="basinmark-demo-key")
    p.add_argument("--external", required=True)
    p.add_argument("--marked", required=True)
    p.add_argument("--keys", type=int, default=200)
    p.set_defaults(func=cross_database_matching)

    p = sub.add_parser("classify-flips", help="classify full-fingerprint flips by mechanism")
    p.add_argument("--engine", required=True)
    p.add_argument("--potential-full", type=Path, required=True)
    p.add_argument("--cifs-dir", type=Path, required=True)
    p.add_argument("--pairs-meta", type=Path, required=True)
    p.add_argument("--symprec", type=float, default=0.1)
    p.add_argument("--cell-bin", type=float, default=0.05)
    p.add_argument("--fmax", type=float, default=0.02)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--resume",
        action="store_true",
        help="resume a checkpoint with matching inputs and settings",
    )
    p.set_defaults(func=classify_flips)

    p = sub.add_parser("self-test", help="run fast statistical and selection checks")
    p.set_defaults(func=run_self_test)

    args = ap.parse_args(argv)
    if hasattr(args, "gamma") and not 0 < args.gamma < 1:
        ap.error("--gamma must be strictly between 0 and 1")
    if hasattr(args, "alpha") and not 0 < args.alpha < 1:
        ap.error("--alpha must be strictly between 0 and 1")
    for name in ("m", "keys", "trials", "collection"):
        if hasattr(args, name) and getattr(args, name) <= 0:
            ap.error(f"--{name} must be positive")
    for name in ("symprec", "cell_bin", "fmax"):
        if hasattr(args, name):
            value = getattr(args, name)
            if not math.isfinite(value) or value <= 0:
                ap.error(f"--{name.replace('_', '-')} must be finite and positive")
    if hasattr(args, "cell_bins") and any(
        not math.isfinite(value) or value <= 0 for value in args.cell_bins
    ):
        ap.error("--cell-bins must be finite and positive")
    if hasattr(args, "observations") and any(value <= 0 for value in args.observations):
        ap.error("--observations must be positive")
    if hasattr(args, "out") and Path(args.out).exists():
        resumable = args.analysis in {"potential-comparison", "classify-flips"}
        if not (resumable and args.resume):
            ap.error(
                "output already exists; choose a new --out path "
                "(potential-comparison supports --resume)"
            )
    try:
        args.func(args)
    except (ImportError, OSError, RuntimeError, ValueError) as exc:
        ap.error(str(exc))


if __name__ == "__main__":
    main()
