# BasinMark 0.4 collection selection

The producer selects 32 interchangeable outputs from the same 256 measured candidates in each trial. The nominal audit level is 0.001. Every attempted collection counts. The 500 trials vary the marking key and candidate draw within each archived cohort.

## Independent follow-up

| Generator | Policy and green probability | Detected after re-relaxation | 95% interval | Mean excess over unmarked pool minimum (meV/atom) | Independently stable yield |
|---|---|---:|---:|---:|---:|
| MatterGen | Unmarked single | 0/500 | 0.0–0.7% | 107.9 | 41.5% |
| MatterGen | Unmarked fixed minimum | 2/500 | 0.0–1.4% | 7.2 | 82.7% |
| MatterGen | Unmarked pool minimum | 0/500 | 0.0–0.7% | 0.0 | 82.2% |
| MatterGen | Unbounded fixed, 0.5 | 500/500 | 99.3–100.0% | 32.1 | 74.0% |
| MatterGen | Fixed capped, 0.5 | 491/500 | 96.6–99.2% | 21.8 | 77.0% |
| MatterGen | Fixed capped, 0.625 | 475/500 | 92.7–96.7% | 17.8 | 79.1% |
| MatterGen | Strict ranks, 0.5 | 500/500 | 99.3–100.0% | 11.4 | 86.3% |
| MatterGen | Strict ranks, 0.625 | 489/500 | 96.1–98.9% | 8.9 | 86.5% |
| MatterGen | Slot and mean caps, 0.5 | 500/500 | 99.3–100.0% | 16.2 | 86.3% |
| MatterGen | Slot and mean caps, 0.625 | 499/500 | 98.9–100.0% | 11.9 | 86.6% |
| CrystalFormer | Unmarked single | 0/500 | 0.0–0.7% | 3386.3 | 13.3% |
| CrystalFormer | Unmarked fixed minimum | 1/500 | 0.0–1.1% | 49.1 | 32.8% |
| CrystalFormer | Unmarked pool minimum | 0/500 | 0.0–0.7% | 0.0 | 38.8% |
| CrystalFormer | Unbounded fixed, 0.5 | 500/500 | 99.3–100.0% | 400.0 | 22.7% |
| CrystalFormer | Fixed capped, 0.5 | 12/500 | 1.2–4.2% | 56.7 | 31.8% |
| CrystalFormer | Fixed capped, 0.625 | 15/500 | 1.7–4.9% | 55.8 | 31.8% |
| CrystalFormer | Strict ranks, 0.5 | 212/500 | 38.0–46.9% | 22.9 | 35.4% |
| CrystalFormer | Strict ranks, 0.625 | 198/500 | 35.3–44.0% | 19.2 | 36.0% |
| CrystalFormer | Slot and mean caps, 0.5 | 455/500 | 88.1–93.4% | 49.3 | 30.7% |
| CrystalFormer | Slot and mean caps, 0.625 | 476/500 | 92.9–96.9% | 37.3 | 31.9% |

MatterGen stable yield counts emissions with an archived DFT hull energy at or below 0.1 eV/atom, divided by all emissions. CrystalFormer stable yield uses the original PhononBench dynamical stability labels. The two labels describe different physical properties. Both are withheld during selection. Each trial file records absolute ranking energies, coverage, validity, formula diversity, distinct fingerprint counts and false rejections under a separate key.

## How the method changed

Fixed batches sometimes contain no affordable green. Starting with four candidates per slot and spending the reserve on unresolved slots recovered detection on development data but exceeded the CrystalFormer quality allowance. The first pooled rule capped every energy rank against the globally cheapest records and the fixed capped selection. It missed the 90% target on the first external cohort. The protocol records the criteria and both decisions. [Design and decision record](protocol.md)

| First evaluation | Fixed capped | Reserve | Strict rank caps |
|---|---:|---:|---:|
| MatterGen | 499/500 | 500/500 | 500/500 |
| CrystalFormer, previous cohort | 359/500 | 499/500 | 465/500 |
| CrystalFormer, external | 8/500 | 374/500 | 158/500 |
| WyCryst | 499/500 | 500/500 | 500/500 |

The follow-up bounds each sorted output energy by the corresponding sorted original fixed-slot minimum plus 100 meV/atom. It also bounds the mean by the lower of the globally cheapest 32 candidates’ mean plus 100 meV/atom and the mean of the fixed capped selector at the same green fraction. The old strict rule and the follow-up use different per-rank references. The follow-up uses 0.625 at production and audit. The table includes both green probabilities to separate the effects of the two changes.

For a proposed green count, take the cheapest required records of each color. Replacing a record with a cheaper one of the same color can only decrease each sorted energy. Scan the possible green counts from largest to smallest and return the first feasible set. This maximizes the number of green emissions and minimizes energy within that count. Separately requested compositions or properties require separate pools.

## Source records, safeguards and checks

Each fresh physics cohort contains 384 sampled structures. The two new MatterGen cohorts share no source IDs with each other or the previous two physics cohorts. Both PhononBench CrystalFormer cohorts exclude all earlier Figshare producer classes under the no-cell fingerprint. The second cohort also excludes the first external cohort’s IDs and classes. The Figshare and PhononBench populations have different preparation and ranking energies. Policy comparisons are paired within each source.

The downstream channel starts from the exact producer geometry, applies CHGNet-r2SCAN variable-cell FIRE at 0.05 eV/Å for at most 150 steps, standardizes the primitive cell at 0.1 Å and writes and rereads a five-decimal CIF. The physics is key-independent and can be precomputed for each candidate. Missing and nonconverged endpoints supply no fingerprint. Fresh independent audit keys test the actual deposited collections.

- **MatterGen.** The follow-up detected 499/500 collections with 100.0% deposited-record coverage. Dropping the lowest-energy record and repeating selection gave 500/500 versus 497/500 for the original fixed comparator. Independent keys rejected 4/10,000 times. The mean attainable null rate was 0.0373%.

- **CrystalFormer.** The follow-up detected 476/500 collections with 100.0% deposited-record coverage. Dropping the lowest-energy record and repeating selection gave 470/500 versus 15/500 for the original fixed comparator. Independent keys rejected 2/10,000 times. The mean attainable null rate was 0.0375%.

On the second external CrystalFormer cohort, key-blind bounded strain followed by independent re-relaxation gave 306/500 detections with the follow-up and 10/500 with the original fixed method. Two lattice vectors are scaled by 1.01 and its reciprocal. Deformations must be geometrically valid and within 10 meV/atom of the starting MACE single-point energy. The first strict exchange rule gave 56/500 detections on the first combined-treatment cohort.

One first-cohort producer scored 54.91 eV/atom below its model hull. Independent single-point and supercell calculations reproduced the extreme prediction. The primary analysis retains it. Dropping it gave 132/500 detections for the first strict rule. A separate composition-mixture linear program reproduced all 384 first-cohort phase-diagram hull references within 2×10⁻¹³ eV/atom. Interpret model energy alongside the independent phonon outcome and sensitivity analyses.

## Geometry and cost

Each selector uses the same no-cell fingerprint. As an optional reference comparison, average the eight nearest-neighbor log distances by rank, normalized by the cube root of volume per atom. A coarse-matched pair passes at RMS profile distance 0.03 or less. On MP-20 test, the number of matcher-distinct same-formula pairs merged fell from 101 to 92 of 684. On the validation split, it fell from 93 to 85 of 752. The check retained both coarse-matched equivalent test pairs and rejected the sole such validation pair. It is a reference-assisted diagnostic with limited discrimination. It is excluded from all watermark tests.

Median CPU selection time per 32 outputs from prepared inputs was 0.33 ms for fixed, 0.62 ms for reserve, 0.35 ms for strict exchange, and 0.35 ms for budget exchange at γ=0.625. All use the same 256 candidate scores. Generating and relaxing 256 candidates has a different cost from generating one per output. Generator inference is outside these measurements. [Upstream relaxation cost](../v0.3/README.md)

## Use and replay

The select-collection CLI reads a JSON object with a candidates list. Each item has a unique string id, a producer fingerprint and a finite energy in eV/atom. Consecutive blocks of eight define the fixed comparator. The key file stores raw bytes. The command returns the selected IDs and energy references.

```bash
python -m pip install -e '.[test]'
basinmark select-collection --candidates candidates.json --key-file key.bin --m 8 --energy-window 0.1 --gamma 0.625 --policy budget --out selected.json
python scripts/replay_allocation.py
python -m pytest
(cd results/v0.4 && shasum -a 256 -c SHA256SUMS)
```

The replay reconstructs each statistical audit from saved measurements, including development, confirmation and exploratory comparisons. Each result stores its own seed, protocol settings, per-trial decisions and summaries. [Regenerate the measurements](reproduce.md)
