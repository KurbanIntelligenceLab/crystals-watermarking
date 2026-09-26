# Physical survival and provenance limits on watermarking AI-generated crystals

[Can Polat](https://orcid.org/0000-0002-1458-302X)<sup>1</sup>, [Mustafa Kurban](https://orcid.org/0000-0002-7263-0234)<sup>2,3,*</sup>, [Erchin Serpedin](https://orcid.org/0000-0001-9069-770X)<sup>1</sup>, [Hasan Kurban](https://orcid.org/0000-0003-3142-2866)<sup>4,*</sup>

1. Department of Electrical and Computer Engineering, Texas A&M University, College Station, Texas, USA.
2. Department of Electrical and Computer Engineering, Texas A&M University at Qatar, Doha, Qatar.
3. Department of Prosthetics and Orthotics, Ankara University, Ankara, Turkey.
4. College of Science and Engineering, Hamad Bin Khalifa University, Doha, Qatar.

**Corresponding authors:** Mustafa Kurban ([kurbanm@ankara.edu.tr](mailto:kurbanm@ankara.edu.tr)) and Hasan Kurban ([hkurban@hbku.edu.qa](mailto:hkurban@hbku.edu.qa)). Names link to ORCID profiles.

## Abstract

When crystals are generated, they are relaxed and standardised before scientific use. However, this can erase or alter structural evidence of their production. We ask what a structure-only watermark can preserve through this processing, how an auditor should calibrate the surviving evidence, and which provenance claims the resulting decision supports. **BasinMark** connects the deposited observation to a collection audit by selecting generated candidates according to keyed structural fingerprints. Exact relaxation to the same basin minimum erases coordinate marks, so the protocol instead selects among alternatives whose fingerprint distinctions may survive deposition. Across the evaluated generators and interatomic potentials, cell-shape bins account for much of the observed fingerprint instability; removing them in two generator populations substantially improves fingerprint-label survival while preserving structural diversity. To ensure rigorous statistical verification, under a key-independent null, counting each unique fingerprint class once provides exact false-positive control. However, because detection operates at the structural class level, collections with identical fingerprints remain indistinguishable—allowing watermark evidence to be copied or reused across datasets. Ultimately, physical processing, fingerprint choices, and structural duplicate handling determine how much watermark evidence actually survives.

## Code and results

Source code and reported results for the manuscript and supplementary information. The package contains fingerprint construction, statistical audits, relaxation comparisons, and flip classification.

Analysis results are available in `results/` as JSON summaries and CSV tables. The supporting code is in `src/`.

## Package layout

```text
src/basinmark/   Installable library and command-line interface
results/         Reported JSON summaries and CSV tables
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
| `generator-controls` | `--gen-dirs`: directories of CIF/VASP structures | Marked, unmarked, and shuffled-key controls after relaxation |
| `relaxation-stability` | `--pairs`: JSON records with `mlip` and `dft` paths | Compare already-relaxed structures with DFT references |
| `potential-comparison` | `--pairs`: JSON records with `raw` and `dft` paths | Relax raw structures with selected potentials and measure flips and distances |
| `fingerprint-ablation` | `--pairs`: JSON records with `raw` and `relaxed` or `dft` paths | Compare fingerprints with and without cell shape |
| `generator-calibration` | `--fingerprint-files`: text files, one fingerprint per line | Compare naive and deduplicated rejection rates |
| `class-reuse` | `--fingerprints`: one fingerprint per line | Simulate reuse of observed green fingerprint classes |
| `cross-database-matching` | `--external` and `--marked`: fingerprint text files | Evaluate shared classes under the supplied marking key |
| `classify-flips` | potential result, pair metadata, and CIF directory | Classify full-fingerprint flips by mechanism |
| `self-test` | none | Run fast statistical and selection checks |

Use `basinmark COMMAND --help` for additional settings. Fingerprint text inputs must be nonempty and contain no blank lines. For cross-database matching, supply the actual marking key with `--key`; the default demonstration key cannot reproduce an independently marked collection.

Statistical commands expose `--gamma` and `--alpha`. Commands that model candidate selection also expose `--m`. Options not used by a command are omitted from its help output.

## Input data

Input structures and manifests are not included in the repository. To run structure-based analyses, supply these inputs in `data/` using the layout below. The statistical self-test and unit tests do not require the experimental datasets.

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
- `relaxation.py`: lazy potential construction and relaxation.
- `classification.py`: resumable flip-mechanism classification.
- `analyses.py`: analysis commands and CLI wiring.

The statistical functions remain independent of model loading. Importing the package or requesting CLI help does not load potential checkpoints.

## Reported results

`results/` contains summaries and tables for analyses reported in the manuscript or supplementary information: representation invariance, relaxation stability, potential comparison, fingerprint granularity, continuous descriptors, flip mechanisms, calibration, controls, class reuse, and cross-database matching.

## Formal mathematical verification

Run `python3 formal/verify.py --clean` to rebuild and audit the Lean proofs. See [the formal coverage map](formal/README.md) for the proved results and the remaining analytic inputs. The formalization checks mathematical statements separately from the empirical analysis scripts.

## Building release archives

Install the build tool and use the release script from the repository root:

```bash
python -m pip install build
python scripts/build_release.py
```

The script creates the wheel and source distribution in `dist/` with machine-specific ownership and timestamp metadata removed from the source archive.
