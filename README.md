# Physical survival and provenance limits on watermarking AI-generated crystals

Anonymous authors.

## BasinMark 0.6

A survival-ranked audit and an optional producer-side stress test are evaluated on independently processed releases at a fixed candidate and output budget. The comparisons include a bounded removal attack with audit feedback and paired phonon checks on edited structures. [Results and replay](results/v0.6/README.md)

## Statistical follow-ups

Two statistics-only analyses on the archived measurements: a reuse-exposure law for cross-catalogue transfer as a function of observed outputs, and a selection-power check for the fixed capped and pooled budget selectors. [Results and scripts](results/stats/README.md)

## BasinMark 0.5

Version 0.5 adds a calibrated recovery audit for strained, independently relaxed crystals. It reads the deposited geometry at two symmetry tolerances. An optional rank calculated from the deposited coordinates adds evidence from higher-energy classes. Every stage protects the earlier audit decisions under the same nominal false-positive level. The producer uses the 0.4 energy-budgeted selector. [Paired treatment results, material resolution and reproducibility](results/v0.5/README.md)

## BasinMark 0.4

Version 0.4 adds fixed-budget selection across interchangeable output slots. It includes paired baselines, downstream treatment of selected structures, independent material-quality labels, wrong-key controls and an optional reference-geometry comparison. [Results, energy guarantees and replay instructions](results/v0.4/README.md)

## BasinMark 0.3

Version 0.3 adds an optional energy ceiling to keyed selection, collection audits after independent re-relaxation, and trial-level calibration and class-reuse benchmarks. The public measurement files support offline replay of the statistical results. See the [evaluation, results and reproduction instructions](results/v0.3/README.md).

## Abstract

Generated crystals are relaxed and standardised before scientific use, so a provenance watermark must survive that processing. Exact relaxation to a common basin minimum erases coordinate marks, and selection audits fail when relaxation flips fine-grained fingerprint bins or when duplicate structures inflate false positives. **BasinMark** is a structure-only protocol that selects generated candidates by keyed structural fingerprints and audits deposited files at the level of distinct fingerprint classes. Across the evaluated generators and interatomic potentials, cell-shape bins account for most fingerprint instability; removing them substantially improves label survival while preserving structural diversity. Energy-budgeted selection that exchanges candidates across interchangeable output slots, combined with class-deduplicated binomial calibration, keeps collections detectable under exact false-positive control. With 32 structures released from 256 relaxed candidates, the audit detects 499/500 MatterGen and 476/500 CrystalFormer collections after independent re-relaxation, standardization and CIF deposition, at 11.9 and 37.3 meV/atom above the cheapest unmarked release, and 496/500 and 497/500 held-out collections both after bounded strain and after a change of relaxation model. Because detection operates on fingerprint classes, collections with identical fingerprints are indistinguishable, which bounds the provenance claim the audit supports.

## Code and results

Source code and reported results for the manuscript and supplementary information. The package contains fingerprint construction, statistical audits, relaxation comparisons, and flip classification.

Analysis results are available in `results/`: JSON summaries and CSV tables at the top level, and gzipped measurement archives with checksums and per-version guides under `results/v0.3` to `results/v0.6` and `results/stats`. The supporting code is in `src/`, with benchmark, measurement and replay drivers in `scripts/`.

## Package layout

```text
src/basinmark/   Installable library and command-line interface
results/         Reported summaries, archived measurements, checksums and per-version guides
scripts/         Benchmark, measurement and replay drivers for the archived results
formal/          Pinned Lean proofs, coverage map, and axiom-checking runner
tests/           Fast tests for the public API and validated I/O
README.md        Installation, usage, inputs, and result map
pyproject.toml   Python package metadata
```

## Installation

Requires Python 3.10 or newer. Run the following from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
basinmark --help
basinmark self-test
```

The base installation provides NumPy and SciPy for statistical calculations. Structure-based commands additionally require:

```bash
python -m pip install -e '.[structures]'
```

The `benchmark` extra adds the catalogue and collection drivers in `scripts/`, and the `phononbench` extra pins the MACE and CHGNet versions used for the archived phonon measurements.

Install the selected potential package and obtain its checkpoint in the same environment: CHGNet, MACE, ORB, SevenNet, or MatterSim. These are optional, separately managed dependencies.

MACE-MPA-0 requires MACE 0.3.10 or newer. The current MatterSim release requires Python 3.12 or newer, so use a separate compatible environment when selecting that engine; the BasinMark statistical package itself supports Python 3.10 or newer.

The same CLI is available as `python -m basinmark`. Paths in the examples below assume the repository root is the working directory.

## Running analyses

Every analysis requires an explicit `--out` path. Choose a new filename so the reported summaries are preserved. Existing output files are rejected, except when explicitly resuming a potential comparison or flip classification.

```bash
basinmark potential-comparison --help
```

With a raw/DFT pair manifest and a configured CHGNet environment, run a potential comparison:

```bash
basinmark potential-comparison \
  --pairs data/potential_pairs.json \
  --engines chgnet \
  --out results/potential_comparison_run.json
```

This command performs structure relaxation and can be computationally expensive. It writes incremental checkpoints. To continue an interrupted run, repeat the command with `--resume`, using the same input files, potential checkpoint, and analysis settings. Resume checks validate the manifest, structure contents, and analysis settings, including already completed engines. Older checkpoints without these hashes require a new output file. Potential weights are managed externally; use the same installed potential version and model weights when resuming.

Compare fingerprint representations on already paired structures without running a potential:

```bash
basinmark fingerprint-ablation \
  --pairs data/potential_pairs.json \
  --out results/fingerprint_ablation_run.json
```

| Command | Required input | Purpose |
| --- | --- | --- |
| `select-collection` | `--candidates`: prepared candidate records, `--key-file`: raw key bytes | Select interchangeable outputs from candidates under a fixed, reserve, pool, budget or preview policy |
| `audit-deposit` | `--structures`: deposited structure files, `--key-file`: raw marking key | Audit deposited structures with symmetry recovery and optional deposited-energy ranking |
| `generator-controls` | `--gen-dirs`: directories of CIF/VASP structures | Marked, unmarked, and shuffled-key controls after relaxation |
| `relaxation-stability` | `--pairs`: JSON records with `mlip` and `dft` paths | Compare already-relaxed structures with DFT references |
| `potential-comparison` | `--pairs`: JSON records with `raw` and `dft` paths | Relax raw structures with selected potentials and measure flips and distances |
| `fingerprint-ablation` | `--pairs`: JSON records with `raw` and `relaxed` or `dft` paths | Compare fingerprints with and without cell shape |
| `generator-calibration` | `--fingerprint-files`: text files, one fingerprint per line | Compare naive and deduplicated rejection rates |
| `class-reuse` | `--corpus`: prepared records with energies, or `--fingerprints`: equal-energy classes | Compare replay and disjoint-record reuse from unfiltered producer observations |
| `key-controls` | `--corpus`: prepared records | Compare unmarked and independent-key audit rates with attainable null rates |
| `collection-benchmark` | `--corpus`: prepared producer and deposition records | Compare full and no-cell fingerprints, bounded selection, and matched-budget controls |
| `cross-database-matching` | `--external` and `--marked`: fingerprint text files | Evaluate shared classes under the supplied marking key |
| `classify-flips` | potential result, pair metadata, and CIF directory | Classify full-fingerprint flips by mechanism |
| `self-test` | none | Run fast statistical and selection checks |

Use `basinmark COMMAND --help` for additional settings. Fingerprint text inputs must be nonempty and contain no blank lines. For cross-database matching, supply the actual marking key with `--key`; the default demonstration key cannot reproduce an independently marked collection.

Statistical commands expose `--gamma` and `--alpha`. Commands that model candidate selection also expose `--m`. Options not used by a command are omitted from its help output.

## Input data

The repository includes the archived measurements for versions 0.3 to 0.6 and the statistical follow-ups; the per-version guides in `results/` give their replay commands. Original structure inputs and manifests are obtained separately. To run structure-based analyses, supply these inputs in `data/` using the layout below. The statistical self-test and unit tests do not require the experimental datasets.

```text
data/
  potential_pairs.json      Runner-ready raw/DFT pair manifest
  structures/
    pairs.json              Pair metadata with paths relative to this directory
    cifs/                   Raw and DFT-reference structures
```

`potential_pairs.json` includes `pair_id`, `raw`, `dft`, and population metadata. Paths are resolved from the working directory first, then relative to the manifest if necessary. The raw/DFT manifest is not an `mlip`/DFT manifest: `relaxation-stability` requires already-relaxed potential outputs.

## Flip classification

After producing per-pair potential results, classify fingerprint changes into discrete changes and changes confined to cell-shape bins:

```bash
basinmark classify-flips \
  --engine chgnet \
  --potential-full results/potential_comparison_run.json \
  --cifs-dir data/structures/cifs \
  --pairs-meta data/structures/pairs.json \
  --out results/flip_classification_run.json
```

The input must contain a `chgnet` object with `per_pair` records. Classification reruns relaxation for flipped pairs and requires the same potential and fingerprint settings. It saves progress atomically. Add `--resume` to continue a checkpoint; checkpoints with different input hashes or settings are rejected.

## Code organization

- `core.py`: reusable fingerprint, selection, and statistical functions.
- `io.py`: validated inputs and atomic JSON output.
- `allocation.py`: budget-matched selection of interchangeable outputs from one candidate collection.
- `collections.py`: collection experiments using measured fingerprints and producer energies.
- `deposit_energy.py`: source-pinned single-point ranking for deposited structures.
- `reference.py`: scale-free local geometry summaries for comparing a deposited record to a reference.
- `reliability.py`: calibrated reliability-ranked extension of the deposited-energy audit.
- `symmetry_audit.py`: paired, deduplicated class audit for two views of a deposited collection.
- `relaxation.py`: lazy potential construction and relaxation.
- `classification.py`: resumable flip-mechanism classification.
- `analyses.py`: analysis commands and CLI wiring.

The statistical functions remain independent of model loading. Importing the package or requesting CLI help does not load potential checkpoints.

## Reported results

`results/` contains summaries and tables for analyses reported in the manuscript or supplementary information: representation invariance, relaxation stability, potential comparison, fingerprint granularity, continuous descriptors, flip mechanisms, calibration, controls, class reuse, and cross-database matching. The archived collection trials, deposition audits, removal and cross-model checks, and the reuse-exposure and selection-power statistics are under `results/v0.3` to `results/v0.6` and `results/stats`, each with a README, checksums and replay instructions.

## Formal mathematical verification

Run `python3 formal/verify.py --clean` to rebuild and audit the Lean proofs. See [the formal coverage map](formal/README.md) for the proved results and the remaining analytic inputs. The formalization checks mathematical statements separately from the empirical analysis scripts.

## Building release archives

Install the build tool and use the release script from the repository root:

```bash
python -m pip install build
python scripts/build_release.py
```

The script creates the wheel and source distribution in `dist/` with machine-specific ownership and timestamp metadata removed from the source archive.
