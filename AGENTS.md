# Repository guide for coding agents

## Purpose and scope

BasinMark provides structural fingerprints, keyed candidate selection, statistical audits, relaxation comparisons, and Lean proofs for crystal watermarking. Read [README.md](README.md) for the paper and command overview and [formal/README.md](formal/README.md) for the precise proof coverage. This guide applies throughout the repository.

Use the checked-out code, configuration, and command help as the source of truth. Preserve existing user changes. Keep public documentation focused on the published methods, supported workflows, and reproducibility; do not describe private workspace organization or internal maintenance history.

## Setup and quick validation

Run commands from the repository root. The package requires Python 3.10 or newer. Select a compatible interpreter before creating an environment; do not replace an existing environment unnecessarily.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]' ruff
python -m basinmark --help
python -m basinmark self-test
python -m pytest
python -m pip check
ruff check src tests scripts
```

For structure handling and the real-structure invariance test:

```bash
python -m pip install -e '.[test,structures]'
python -m pytest
```

The base test suite does not require experimental datasets or model weights. The real-structure test skips when pymatgen is absent; report skips explicitly. Potential-comparison regression tests use a substituted relaxer and do not establish correctness of an installed model. Do not claim empirical reproduction from unit tests or the statistical self-test.

Potential packages and their weights are separate dependencies. Install only the engines needed for the task, following their upstream compatibility requirements. Keep model imports lazy: importing `basinmark`, statistical commands, and CLI help should not require torch or load model weights.

## Code map

| Path | Responsibility |
| --- | --- |
| `src/basinmark/core.py` | Fingerprints, keyed selection, statistical functions, and self-test |
| `src/basinmark/io.py` | Manifest and fingerprint readers; atomic JSON output |
| `src/basinmark/analyses.py` | CLI arguments, analysis workflows, and potential checkpoints |
| `src/basinmark/classification.py` | Flip classification and validated resume state |
| `src/basinmark/relaxation.py` | Lazy calculator construction and ASE relaxation |
| `src/basinmark/__init__.py` | Public API and package version |
| `tests/` | API, I/O, checkpoint, structure-invariance, and archive regressions |
| `results/` | Reported numerical summaries and tables |
| `formal/` | Lean sources, pinned dependencies, and compiled-declaration audit |
| `scripts/build_release.py` | Distribution build and source-archive metadata cleanup |

Keep statistical logic independent of model loading and CLI parsing. Reuse the existing I/O and relaxation helpers. Follow the Ruff configuration in `pyproject.toml` and preserve Python 3.10 compatibility for the base package. Add focused regression tests for changed behavior; avoid unrelated formatting or API changes.

## Inputs and reproducibility

Input structures and manifests are not shipped. Do not fabricate missing inputs, silently substitute datasets, or present synthetic examples as paper reproduction. Use `python -m basinmark COMMAND --help` before preparing a run.

- `potential-comparison` needs a nonempty JSON list with `raw` and `dft` structure paths. Supply unique string `pair_id` values for records that will be classified later.
- `relaxation-stability` needs `mlip` and `dft` paths, with potential outputs already relaxed.
- `fingerprint-ablation` needs `raw` and `relaxed` paths; `dft` is accepted as the relaxed reference when `relaxed` is absent.
- Relative structure paths resolve from the working directory first, then relative to the manifest. Use an unambiguous layout and check path resolution before expensive runs.
- Fingerprint corpora are nonempty text files with one fingerprint per line and no blank lines. Preserve repeated records when measuring duplicate effects.
- `class-reuse` samples without replacement from the input corpus before filtering green classes. It needs at least four times the largest requested observation count in input records.
- `cross-database-matching` needs the actual marking key to evaluate an independently marked collection. Its default key is a demonstration value.
- `classify-flips` needs per-pair potential output, pair metadata, and the matching CIF files. Aggregate summaries alone are insufficient. Classification reruns relaxation for flipped pairs.

For a reproducibility run, record the code revision, Python and dependency versions, input hashes, potential version and weight identity, random seed, CLI arguments, and output location. Run a small representative input first. Distinguish a smoke test from a full experiment and state exactly which engines and datasets were exercised.

## Results and checkpoints

Treat tracked `results/` files as reference evidence. Do not overwrite or regenerate them unless that is the requested task. Put exploratory outputs in an ignored location, such as `data/runs/`, and provide an explicit new `--out` path.

Use one JSON output file per run for progress and final results. Preserve atomic writes and configuration checks. Resume only with the same inputs, settings, and model weights, using `--resume` on supported commands. Potential checkpoints validate the manifest, structure contents, and settings; classification checkpoints also validate their input metadata and records. Model weights are externally managed. Never edit stored hashes or configuration fields to bypass a mismatch; start a new output instead.

Preserve fingerprint semantics, deduplication defaults, statistical thresholds, and reported values unless the task explicitly changes them. Check statistical changes against the mathematical assumptions and add tests that expose the intended difference. A structural-class audit is not a unique attribution of generator provenance.

## Formal verification

Use the versions pinned in `formal/lean-toolchain`, `formal/lakefile.toml`, and `formal/lake-manifest.json`. Install elan when needed; initial dependency setup can require network access.

```bash
python formal/verify.py --clean
```

The clean option removes this project's generated build directory while retaining dependency caches. The verifier builds the imported modules and audits compiled declarations. Keep new modules imported by `formal/BasinMark.lean`. Do not add `sorry`, new axioms, or weaken statements merely to make a build pass. Preserve upstream licenses and adaptation notes.

Update the coverage map when changing formal claims. Lean verification does not verify the Python implementation, model behavior, or empirical results; describe only the propositions and assumptions actually checked.

## Release and Git hygiene

For distributable archives, install `build` and use the supplied wrapper:

```bash
python -m pip install build
python scripts/build_release.py
```

Inspect both archives in `dist/`, including file lists, package metadata, tar ownership fields, and gzip headers. The wrapper strips machine-specific source-archive metadata; it does not anonymize README contents or Git history. Confirm that required formal sources, the toolchain pin, and licenses are included, and that datasets, environments, caches, credentials, and personal filesystem paths are absent. Check an installed wheel outside the checkout when changing packaging.

Respect `.gitignore` and `MANIFEST.in`; do not force-add ignored files. Preserve required third-party attribution. Keep this branch anonymous: exclude author names, affiliations, contact details, ORCID links, and identifying repository URLs from its contents. Use `Anonymous <anonymous@example.invalid>` for both Git author and committer on this branch, without assistant co-author trailers. Commit and push only within the user's requested scope; do not rewrite history or force-push as routine cleanup.

Before committing, inspect `git status`, the staged file list and diff, and run `git diff --cached --check`. For code changes, run the relevant tests, self-test, lint, and dependency check. For proof changes, run the formal verifier. For documentation-only changes, verify commands, links, and claims against the implementation without running unrelated expensive experiments. Repeat release verification and artifact scans against the final files after fixes.

Report the changes, checks actually run, skipped or unavailable checks, and any remaining reproduction requirements. Write documentation in Markdown and ground scientific claims in repository evidence and primary references.
