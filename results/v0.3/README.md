# BasinMark 0.3 collection evaluation

BasinMark selects one crystal per candidate batch using a keyed fingerprint class and an energy ranking. Version 0.3 adds an energy window. A green candidate is eligible when its ranking energy is within the window above the batch minimum. The producer takes the lowest-energy eligible green candidate, or the batch minimum if none qualifies. The public API preserves the original unconstrained selection by default.

## Post-deposition detection

Each confirmation cohort contains 384 sampled producer records with no source ID in its generator’s 512-record development cohort. Each trial draws a fresh key and permutes the eligible cohort into batches of eight. The same batches and key are used for every method. The primary audit contains 32 emissions. The nominal level is 0.001 with green probability 0.5. The table reports rejecting collections out of 500 trials, conditional on each archived cohort. The intervals quantify variation over keys and batches on these fixed records. Generalization to fresh generator runs requires new candidate populations.

| Generator | Producer policy | Independent re-relaxation | Paired DFT archive | Mean producer penalty (meV/atom) | Confirmed stable yield at DFT |
|---|---|---:|---:|---:|---:|
| MatterGen | Full fingerprint | 61/500 | 10/500 | 29.1 | 71.7% |
| MatterGen | No cell | 500/500 | 500/500 | 29.4 | 71.9% |
| MatterGen | No cell, 100 meV/atom window | 459/500 | 454/500 | 14.0 | 75.2% |
| MatterGen | Unmarked lowest of eight | 0/500 | 1/500 | 0.0 | 80.5% |
| MatterGen | Unmarked single candidate | 0/500 | 1/500 | 99.2 | 42.0% |
| CrystalFormer | Full fingerprint | 428/500 | — | 105.6 | — |
| CrystalFormer | No cell | 500/500 | — | 103.5 | — |
| CrystalFormer | No cell, 100 meV/atom window | 344/500 | — | 9.3 | — |
| CrystalFormer | Unmarked lowest of eight | 1/500 | — | 0.0 | — |
| CrystalFormer | Unmarked single candidate | 0/500 | — | 471.9 | — |

Producer energies are the released corrected CHGNet hull estimates. The penalty is the selected ranking energy minus the minimum in its own batch. The bound applies to that ranking field. The DFT quality metric counts emissions with an available reference energy at or below 0.1 eV/atom above the reference hull, divided by all emissions. Coverage and validity appear beside each audit in the trial files. Structural validity uses pymatgen’s minimum-distance check.

**Processing order.** Start from the released neural-generator CHGNet geometries. Rank and select on these producer geometries. Audit the selected geometry directly, after primitive standardization and a five-decimal CIF write/read, and after an independent CHGNet-r2SCAN relaxation followed by the same CIF deposition. The latter starts from the exact producer coordinates. The upstream producer CHGNet model and the r2SCAN-transfer checkpoint belong to the same model family. The downstream run uses ASE FIRE with a variable cell, force threshold 0.05 eV/Å and a 150-step cap. A nonconverged or failed structure contributes no deposited fingerprint. Every attempted collection stays in the denominator. Precomputing a deterministic, key-independent treatment for each candidate allows exact replay after selection without repeating the physics for each key.

The MatterGen DFT column uses a separately archived same-generation-ID reference plus a reduced-formula check. That archive has no per-record starting coordinates for the DFT job. Only the r2SCAN column verifies processing that starts from the selected producer geometry. Exact intervals, distinct class counts, green fractions, stage coverage, valid yields and trial records are in the [MatterGen results](trials/MatterGen-confirm.json.gz) and [CrystalFormer results](trials/CrystalFormer-confirm.json.gz).

Of the 384 confirmation candidates per generator, 381 MatterGen and 330 CrystalFormer structures converged under the downstream relaxation budget. One CrystalFormer structure also lacked a usable producer fingerprint. All selection policies use the same eligible producer pool for each generator. Ranking compares candidates across the released composition distribution. A composition-constrained production task would require its own comparison.

At the 100 meV/atom window, 9.8% of MatterGen and 17.8% of CrystalFormer emissions fell back to the batch minimum at 32 outputs. The bounded selector exceeded the 90% observed detection target on MatterGen and reached 68.8% on CrystalFormer. At the separately evaluated size of 40 outputs, the corresponding counts were 496/500 and 448/500. On CrystalFormer the unconstrained no-cell selector retained higher detection and incurred a larger energy penalty. The 32-output comparison remains the primary result.

**Development choice.** At 32 emissions, the MatterGen 50, 100 and 150 meV/atom windows gave 188, 280 and 296 detections in 300 trials on paired DFT structures. The 100 meV window was the smallest tested value to reach 90% observed detection on both MatterGen downstream channels. The confirmation uses that same window and fresh keys. Development and confirmation results for both generators are retained.

## Material discrimination and class reuse

The MP-20 test split contains 377 reduced-formula groups with multiple records. Under the specified StructureMatcher settings, 684 within-formula pairs were structurally distinct. The full fingerprint merged 18 of them (2.6%) and the no-cell fingerprint merged 101 (14.8%). Same-formula pair denominators include every eligible pair, so formulas with more structures contribute more pairs.

A fixed 27,000-entry OQMD sample from the public JARVIS distribution shared 165 full fingerprint classes and 267 no-cell classes with the 9,046 MP-20 test structures. The source sample excludes all occurrences of repeated OQMD IDs. In each of 200 trials, an observer sees 1,130 genuine marked MP-20 emissions and resamples 200 matching OQMD records with replacement. Both fingerprints gave 200/200 rejections. Excluding OQMD records that match any exposed structure under StructureMatcher yielded 0/200 for both. An empty pool counts as a non-rejection. Matcher equivalence uses primitive-cell reduction, scale normalization and tolerances of 0.2 for lattice lengths, 0.3 for sites and 5 degrees for angles. The result depends on the two catalogues, sample, and equivalence rule. [Catalogue trial results](trials/catalogues.json.gz)

The within-generator observer also sees the actual sequence of emissions, including non-green fallbacks. The earlier plot used aggregate statistics. The published historical helper prefiltered the observer’s inputs by green status, and the historical archive lacks the per-trial records needed to certify how the plot handled fallbacks. The following results come from the new implementation. At ten observations, the full-fingerprint replay-allowed attack rejected in 194/200 MatterGen trials and 197/200 CrystalFormer trials. With observed row identities excluded, the corresponding rates were 0/200 and 0/200. Sampling probabilities follow the frequency of matching records. Trial files store the seed, emitted row indices, fallback count, eligible pool size, distinct and green counts, and audit decision. [MatterGen reuse](trials/MatterGen-reuse.json.gz) · [CrystalFormer reuse](trials/CrystalFormer-reuse.json.gz)

## Null calibration and effective sample size

On WyCryst, the unmarked minimum-energy control gave 15/20,000 rejections (0.075%) and an independently keyed audit of marked outputs gave 10/20,000 (0.050%). Each trial uses two fresh keys on 123 fixed batches. The attainable conditional null rate averaged 0.0542% in both arms. The earlier 2/200 rates at m = 8 remain in the historical archive. Those trials used different keys and populations. A pooled multiple-testing correction answers a separate question. [New controls, intervals and per-key counts](trials/WyCryst-controls.json.gz) · [Historical controls](../generator_controls.json)

Fix the fingerprint rule, sample rule, threshold and key before observing the collection. For a real provenance claim, register a binding commitment to the secret high-entropy key and these rules with an independent timestamp before inspection. The software’s public benchmark seeds exist for replay. If L prespecified keys are tried, a Bonferroni per-key threshold of α/L controls the familywise false-positive probability under the applicable key-independent null. Searching an unrecorded set of keys or choosing matching records adaptively defeats the single-audit interpretation. A positive audit provides evidence of class association.

Current `n_eff` uses the observed class weights and equals the inverse sum of their squares. `n_eff_population_estimate` preserves the separate population collision estimator. The audit rejection uses the actual distinct-class count K. For seven equally populated classes of nine, the fields are 7, 7.75 and 7, respectively. In the original no-cell populations the observed effective sizes are 1,490.1 (MatterGen) and 981.1 (CrystalFormer). The corresponding K values are 1,495 and 985. [Archived ablation](../fingerprint_ablation.json)

## Constrained strain and computational cost

Apply a key-blind, volume-preserving anisotropic strain to archived MatterGen producer structures. Multiply two lattice vectors by 1.01 and its reciprocal, holding fractional coordinates fixed. Accept the modification only if the geometry is valid and the absolute CHGNet energy change is at most 10 meV/atom. Otherwise deposit the original. All 256 proposals passed this budget. With 32 emissions, full, no-cell and bounded no-cell selection gave 56/300, 300/300 and 299/300 detections, respectively. This evaluates one physically bounded edit. [Strain trial records](trials/strain.json.gz)

A separate CPU profile relaxed 128 raw MatterGen structures in 16 batches of eight. CHGNet converged for 126 within the step budget. All attempted timings are included. Mean relaxation work per emission was 0.867 worker-seconds for the predesignated single candidate and 6.022 for both lowest-of-eight and marked-of-eight. Fingerprinting eight structures and selecting one emission added 0.048 seconds on average. Worker-seconds sum per-candidate durations under concurrent execution. Generator inference and initial model loading are outside the measurement. [Per-structure observations and hardware metadata](measurements/processing.json.gz)

## Input identity and energy checks

The structure and energy tables in the public archives require explicit joins. The CrystalFormer geometry ID is `material_id` and its energy ID is `str_ind`. The corresponding WyCryst columns are `material_id` and `folder_ind`. Pairing by row position gives different formulas on 974/992 CrystalFormer rows and 993/994 WyCryst rows. Joining by ID and checking reduced composition retains 988 and 991 geometry–energy pairs before fingerprint exclusions. The signed energy recorded inside the source error string also preserves a scientific exponent omitted by the numeric column on 10 and 14 energy records. The loader restores only verified mantissa/exponent cases and rejects an unexplained discrepancy. All 10,000 MatterGen producer formulas agree with their source filenames. Of these, 9,764 have finite ranking energy.

## Replay

From the source checkout, install the base package and run the replay driver. It reconstructs each confirmation trial, the two within-generator reuse studies and the WyCryst controls and compares all records with the published outputs. `--all` includes the development and strain trials. Compact gzip files contain ordinary JSON. The catalogue study can be independently rerun from its public source files using the command below.

```bash
python -m pip install -e '.[test]'
python scripts/replay_benchmarks.py --all
python -m pytest
(cd results/v0.3 && shasum -a 256 -c SHA256SUMS)
```

For numerical replay the measured environment used Python 3.12.10, NumPy 2.5.3 and SciPy 1.18.1. For structure regeneration it used pymatgen 2026.9.24, spglib 2.7.0, ASE 3.29.0, CHGNet 0.4.2, torch 2.15.0 and ijson 3.6.0. Potential checkpoint hashes are stored with the measurements. Wall times depend on hardware, scheduling and load.

```bash
python -m pip install -e '.[benchmark]'
mkdir -p data/archives data/runs
for id in 65498772 65498775 65498808 65498610 65498871 65498655; do
  curl -fL "https://ndownloader.figshare.com/files/$id" -o "data/archives/$id.csv.gz"
done
curl -fL https://ndownloader.figshare.com/files/26790182 -o data/archives/26790182.zip

python scripts/prepare_collection.py --dataset MatterGen --directory data/archives \
  --limit 512 --seed 1729 --engine chgnet-r2scan --paired-dft --workers 3 --out data/runs/mg-dev.json
python scripts/prepare_collection.py --dataset MatterGen --directory data/archives \
  --limit 384 --seed 7321 --exclude-cohort data/runs/mg-dev.json --engine chgnet-r2scan \
  --paired-dft --workers 2 --out data/runs/mg-confirm.json
python -m basinmark collection-benchmark --corpus data/runs/mg-confirm.json --sizes 16 32 40 \
  --energy-window 0.10 --seed 98123 --trials 500 --out data/runs/mg-audits.json
python scripts/catalogue_overlap.py --directory data/archives --out data/runs/catalogues.json
python scripts/profile_processing.py --directory data/archives --out data/runs/processing.json
```

For CrystalFormer use the same cohort sizes, seeds and procedure with `--dataset CrystalFormer`, omit `--paired-dft`, and use new output paths. For full producer corpora omit the limit, engine and excluded cohort. The `class-reuse` and `key-controls` commands consume those prepared files directly. The replay driver contains the exact audit sizes and settings via their recorded input results. The public source snapshots and fixed IDs make reconstruction independent of live database queries.

Derived structure measurements use Kazeev and collaborators’ generated-crystal benchmark and MP-20 split, and Choudhary’s JARVIS distribution of OQMD. Both source archives specify CC BY 4.0. This release adds fingerprints, identifier-based joins, independent processing and aggregate and trial statistics. [Generated-crystal data and attribution](https://doi.org/10.6084/m9.figshare.29145101.v6) · [OQMD snapshot and attribution](https://doi.org/10.6084/m9.figshare.14206169.v1) · [Data license](https://creativecommons.org/licenses/by/4.0/)
