# Auditing deformed crystal collections

BasinMark 0.5 uses the same 32 structures selected by the 0.4 energy-budget rule from each pool of 256 candidates. It changes the audit of those structures after processing. The original audit reads the deposited no-cell fingerprint at a symmetry tolerance of 0.10 Å. The recovery audit also reads the same coordinates at 0.30 Å and counts secondary classes with a space-group number of at least three. The original view covers all resolved groups. [Processing order, energy constraints and material-quality definitions](../v0.4/README.md)

## Paired results

The original CrystalFormer strain cohort was used for development. The first fresh cohort independently evaluated the geometry-only change and subsequently informed the energy-ranked extension. The two reserved groups evaluate the frozen ranked method. Each contains 384 source records, disjoint from the other group and all earlier external source IDs and producer classes. All 500 trials count. Intervals reflect keys and candidate draws within a fixed cohort. Every column audits the same selected IDs. Each energy, stability label, validity and diversity outcome is therefore identical within a paired comparison.

| Cohort | Treatment | Original | Geometry recovery | With deposited energy | 95% interval, ranked | Independent-key rejections |
|---|---|---:|---:|---:|---:|---:|
| Original CrystalFormer cohort | Axial strain + relaxation | 306/500 | 479/500 | 487/500 | 95.6–98.6% | — |
| First fresh CrystalFormer cohort | Axial strain + relaxation | 299/500 | 414/500 | 447/500 | 86.4–92.0% | 11/10,000 |
| Reserved CrystalFormer A | Axial strain + relaxation | 269/500 | 434/500 | 447/500 | 86.4–92.0% | 10/10,000 |
| Reserved CrystalFormer B | Axial strain + relaxation | 305/500 | 433/500 | 462/500 | 89.7–94.6% | 10/10,000 |
| Reserved CrystalFormer A | Ordinary relaxation | 467/500 | 468/500 | 489/500 | 96.1–98.9% | 7/10,000 |
| Reserved CrystalFormer B | Ordinary relaxation | 434/500 | 445/500 | 489/500 | 96.1–98.9% | 5/10,000 |
| Fresh MatterGen | Ordinary relaxation | 500/500 | 500/500 | 500/500 | 99.3–100.0% | 13/10,000 |
| Fresh MatterGen | Archived DFT | 492/500 | 495/500 | 495/500 | 97.7–99.7% | 7/10,000 |
| First fresh CrystalFormer | Shear + relaxation | 379/500 | 435/500 | 468/500 | 91.1–95.6% | 10/10,000 |

The prespecified target was 90% observed detection on strained CrystalFormer. The two reserved groups yielded 89.4% and 92.4%. The records give the observed independent-key count and the attainable finite-sample null probability for each channel. The 1.5% shear was applied without consulting the marking key and was allowed only within the 10 meV/atom energy budget. Trials vary keys and candidate draws within an archived cohort. [Frozen design](protocol.md) · [First reserved group](CF-A-audits.json.gz) · [Second reserved group](CF-B-audits.json.gz)

## Secondary 40-output operating points

The 32-output comparison above is primary. These geometry-only 40-output releases use different pools. Six candidates per output cost 240 candidates, and eight cost 320. The table reports strain plus independent relaxation and the stable share of emitted structures. Energy excess is relative to the globally cheapest 40 of the same trial candidates, in meV/atom. Each release size is set before auditing.

| Cohort | Candidates | Original | Recovery | Phonon-stable share | Mean excess (meV/atom) |
|---|---:|---:|---:|---:|---:|
| confirm | 240 | 347/500 | 436/500 | 29.0% | 57.3 |
| confirm | 320 | 398/500 | 460/500 | 34.1% | 53.0 |
| A | 240 | 367/500 | 456/500 | 30.3% | 49.4 |
| A | 320 | 385/500 | 471/500 | 36.3% | 41.5 |
| B | 240 | 381/500 | 464/500 | 22.6% | 48.8 |
| B | 320 | 403/500 | 474/500 | 26.9% | 44.4 |

## Why the strain test lost power

In the original CrystalFormer comparison, ordinary and strained deposition each yielded 15,995 distinct observations summed over 500 collections. Among paired observations, 3,007 changed both the space-group and site fields and 209 changed only the group. Mean green count fell from 30.526 to 28.880. At 32 distinct classes the original test needs 29 green classes. The second tolerance recovered a fraction of the producer assignments. Symmetry resolution in spglib depends on geometric tolerance. [Measured counts, transitions and tolerance sweep](strain-diagnosis.json.gz) · [Spglib symmetry-search algorithm](https://doi.org/10.1080/27660400.2024.2384822)

## Decisions and calibration

Each audit key and procedure are fixed before inspecting the data. Within each view, deduplicate the strings. The primary and secondary sets overlap. Give every distinct string in their union one keyed label shared across all its occurrences. Under a key-independent collection and ideal random-function labels at green probability 5/8, counts in disjoint sets are independent binomials. Use their joint distribution to protect the original rejection region, then add outcomes in descending secondary green count and primary green count until the next outcome exceeds the specified size. This discrete calibration works even when the same class occurs in both views. It preserves every original rejection and keeps the combined level at most 0.001.

The optional third view orders distinct primary classes by a single-point MACE-OMAT-0-medium energy of the deposited structure relative to a composition-matched public LeMat hull. For each primary class choose the lowest scored representative, breaking ties by its recovery class. Drop the lowest quarter, rounded down, and count the eligible recovery classes in the remainder. Low-energy fallbacks often contribute less evidence because the producer's budget can force their selection across keys. On the first fresh CrystalFormer strain cohort, the upper group was green on 95.63% of its distinct recovery-class observations, compared with 86.44% among the remaining classes, aggregated across 500 trials. The audit obtains its own scores from the deposited coordinates. CrystalFormer uses the same pinned potential for both rankings, applied to different geometries. The MatterGen producer uses a CHGNet ranking and the optional auditor uses MACE-OMAT. It needs no producer-energy list, release order or phonon labels. Rank is measured in eV/atom.

The ranked subset is contained in the second view. Partition the union into the five membership patterns induced by these three sets. Convolve the independent Bernoulli labels and retain the full preceding two-view rejection region. Add outcomes ordered by ranked green count, secondary count and primary count while the joint probability stays within the level. The rank extension applies to collections with up to 32 distinct classes in each view and available energy references. Larger collections and missing reference compositions retain the two-view decision. Exact enumeration tests check small joint nulls, duplicates, missing views and identical views. Key-shopping and adaptive collection construction require separate handling. [Key commitment and reuse controls](../v0.3/README.md)

## Use

For a geometry-only audit install the structures extra. The key file contains raw marking-key bytes. Supply the same green probability used at production. Aim to audit the full predetermined collection. Each CLI input is an already-deposited structure. For deposited energy ranking install the phononbench extra and supply both pinned files. The reported decisions and attained joint null probabilities are conditional on the fixed views and their overlap.

```bash
python -m pip install -e '.[structures]'
basinmark audit-deposit --structures deposited/*.cif --key-file key.bin --out geometry.json
python -m pip install -e '.[phononbench]'
basinmark audit-deposit --structures deposited/*.cif --key-file key.bin \
  --weights data/archives/omat.model --hull data/archives/hull.parquet --out ranked.json
python scripts/replay_symmetry_audit.py
```

## Structural resolution and cost

In the repeated-formula MP-20 test pairs, the primary view merged 101 of 684 matcher-distinct pairs. Allowing either direction of recovery matched 105. Among the 55 matcher-equivalent pairs, matches rose from 2 to 23. On validation the corresponding counts rose from 93 to 97 of 752 distinct pairs and from 1 to 30 of 67 equivalent pairs. These are candidate class matches, a different endpoint from collection rejection. A positive audit supports class association. Reproducing the primary, secondary and ranked class sets reproduces all three counts for every key. The earlier single-view reuse experiment directly assesses the primary set. [Test pairs](MP20-test-resolution.json.gz) · [Validation pairs](MP20-validation-resolution.json.gz)

On 64 source structures using one CPU thread, median time per structure was 2.57 ms for the primary fingerprint, 2.04 ms for the additional symmetry view and 0.189 s for deposited-energy ranking. Initializing the pinned model, loading the reference and completing one warm-up structure took 8.18 seconds on this host. The score was checked against an independent source calculation for each structure. Generator and producer relaxation costs are shared across the three audit comparisons. [Profile and hardware](deposit-cost.json.gz) · [Physics and statistical replay](reproduce.md)
