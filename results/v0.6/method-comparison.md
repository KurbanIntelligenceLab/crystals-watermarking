# Original and current BasinMark on the same CrystalFormer candidates

The stronger original baseline already omits cell shape. Under 1% strain followed by CHGNet-r2SCAN relaxation, the current default recovered 37 additional collections on E and 83 on F, net, while choosing crystals with lower predicted energy. The optional preview added four more detections on each. Under ordinary relaxation, both current producers detected two fewer collections per cohort. Each percentage below uses 500 attempted collections. Paired entries list **E / F**.

## Same candidates and keys, each method chooses its own releases

| Pipeline | Ordinary CHGNet | Strain + CHGNet | Strain + MACE | Excess selection energy (meV/atom) | Extra preview relaxations per collection |
|---|---:|---:|---:|---:|---:|
| Original, with cell shape | 4.8% / 7.6% | 1.8% / 4.0% | 1.0% / 1.8% | 298.5 / 396.2 | 0 |
| Original, without cell shape | 100.0% / 100.0% | 91.0% / 82.0% | 98.2% / 96.0% | 293.0 / 393.8 | 0 |
| Current default | 99.6% / 99.6% | 98.4% / 98.6% | 99.4% / 98.0% | 38.3 / 47.5 | 0 |
| Current with preview | 99.6% / 99.6% | 99.2% / 99.4% | 99.2% / 98.6% | 41.2 / 51.9 | 256 |

Each trial gives each producer the same ordered 256 candidates, key and predicted energies. Each emits 32 original geometries. Both original rows implement the submitted rule with 32 fixed batches of eight, no energy cap, green probability 1/2 and a deduplicated binomial audit. The second row implements the paper’s cell-shape ablation and is our principal comparison. The current rows pool interchangeable slots, impose energy limits, use green probability 5/8 and apply the survival-ranked auditor. The optional producer also stress-tests a copy of each candidate. All audits use the same nominal 0.1% level. The table measures the combined effect of the specified changes. [Producer design](producer-protocol.md) · [Audit design](protocol.md)

E and F are the two previously evaluated groups of 384 CrystalFormer structures from PhononBench. Trials vary candidate draws and keys within each group. The head-to-head analysis reuses them. The physical benchmark starts with the archive’s relaxed candidates. The ordinary treatment independently relaxes, standardizes and writes a five-decimal CIF. The strained treatment first applies a volume-preserving 1% strain when it passes the original geometry and energy check, then relaxes and deposits the crystal. Unconverged records contribute no fingerprint. Every attempted collection stays in the denominator. MACE here means MACE-MPA-0 medium. [Source dataset](https://zenodo.org/records/18185662) · [Processing protocol](cross-model-protocol.md)

## Same released crystals, change only the auditor

For each original no-cell release, both auditors get the **identical deposited files and original marking key**. Each uses the original green probability of 1/2. The current auditor reads two symmetry tolerances and uses a single-point deposited-energy score. This isolates the benefit of changing the reader. The producer chooses precisely the same crystals, so its material-energy and compute differences are both zero.

| Processing | Original auditor on original emissions | Current auditor on those exact emissions | Wrong keys firing with current auditor (per 10,000) |
|---|---:|---:|---:|
| Ordinary CHGNet | 100.0% / 100.0% | 100.0% / 100.0% | 14 / 16 |
| Strain + CHGNet | 91.0% / 82.0% | 100.0% / 99.2% | 9 / 12 |
| Strain + MACE | 98.2% / 96.0% | 100.0% / 99.8% | 17 / 9 |

For example, trial 5 on E (seed `6630270238399065`) contains 32 converged CHGNet strain deposits from the original no-cell producer. The original reading found 25 green classes out of 32, below the detection threshold. At the second symmetry tolerance, 30 of 32 classes were green. The calibrated combined audit detected that same collection. This is the first trial on E where the new reader recovered a miss.

We used the 20 previously specified independent control keys per collection. The trial archive records false detections for both readers under these same keys. The theoretical calibration assumes collections independent of the verification key.

## Energy and compute

Selection energy in the first table is the mean excess **over the cheapest 32 candidates in each trial**. Lower is better. The underlying score is MACE-OMAT energy relative to a composition-matched reference hull. Both current selectors had lower mean producer energy than the original no-cell rule in all 500 paired trials on each cohort. For preview versus original, the median paired savings were 156.4 meV/atom on E and 256.8 on F. Endpoint energies below use that same scorer on converged deposits with available scores. The archive records the number scored for each release. Phonon labels refer to source crystals and are withheld during selection. For the CHGNet strain treatment, the original no-cell arm contributed 15,544 and 15,388 converged emission observations out of 16,000 attempts on E and F. The default contributed 15,808 and 15,994. Repeated appearances of a crystal count here in proportion to its selection frequency.

| Producer | CHGNet strained endpoint energy (meV/atom) | MACE strained endpoint energy (meV/atom) | Source phonon pass rate |
|---|---:|---:|---:|
| Original, with cell shape | 346.9 / 431.6 | 352.0 / 446.4 | 18.3% / 19.4% |
| Original, without cell shape | 344.6 / 434.6 | 345.2 / 438.0 | 18.5% / 19.1% |
| Current default | 142.6 / 168.7 | 121.7 / 144.9 | 22.7% / 25.4% |
| Current with preview | 145.4 / 174.0 | 124.5 / 150.2 | 22.7% / 25.2% |

The workload counts below are the primary compute comparison. Every producer starts with 256 candidates and every collection sends 32 emitted crystals for downstream processing. The default producer adds zero relaxation calls. The optional preview adds eight per emitted crystal. For the MACE downstream column we sum the actual optimizer steps of every selected endpoint, including attempts that reached the step limit, and divide by 16,000 attempted emissions. Steps describe the work under this model and optimizer. Structure sizes and model choice also affect cost.

| Pipeline | Extra producer relaxations / 32 outputs | Extra symmetry readings / converged output | Deposited energy scores / converged output | Downstream MACE optimizer steps / output (E / F) |
|---|---:|---:|---:|---:|
| Original, with cell shape | 0 | 0 | 0 | 26.6 / 27.4 |
| Original, without cell shape | 0 | 0 | 0 | 26.5 / 27.2 |
| Current default | 0 | 1 | 1 | 20.9 / 20.8 |
| Current with preview | 256 | 1 | 1 | 21.1 / 21.0 |
| Original outputs, current auditor | 0 | 1 | 1 | 26.5 / 27.2 |

### Shared-laptop timing observations

The laptop was also in use for other work, and machine load was uncontrolled. Use these observations for local context. A reproducible comparison of latency, throughput or speedup requires a controlled machine. Detection, model energies, call counts and recorded optimizer steps are calculated independently of the timers.

The timing table reports seconds **per collection of 32**. Selection, external CHGNet strain relaxation plus deposition, and audit are measured process CPU times on the first registered draw in each cohort. Selection and audit use medians over three interleaved runs after warm-up. Relaxation is timed once per distinct emitted structure and charged to each release that uses it. The preview column instead gives mean measured worker-seconds for all 256 preview attempts over the full 500 draws per cohort. These timings use different clocks and samples. Each column stands on its own.

| Pipeline | Selection, including primary fingerprints | Stress preview | Relaxation and CIF deposition | Reading the deposited collection |
|---|---:|---:|---:|---:|
| Original, with cell shape | 0.88 / 1.03 | 0 / 0 | 45.2 / 29.7 | 0.08 / 0.14 |
| Original, without cell shape | 0.88 / 0.99 | 0 / 0 | 35.4 / 18.0 | 0.16 / 0.17 |
| Current default | 0.92 / 1.00 | 0 / 0 | 21.4 / 16.6 | 2.06 / 2.48 |
| Current with preview | 0.92 / 0.88 | 547 / 479 | 21.4 / 15.4 | 2.23 / 2.27 |
| Original outputs, current auditor | 0.94 / 0.97 | 0 / 0 | 35.4 / 18.0 | 2.32 / 2.40 |

On the identical original no-cell deposits, the extra audit cost was 2.17 CPU seconds on E and 2.23 on F for the profiled collections. These are two first-draw timing examples.

All methods start with the same 256 archive candidates and their ranking energies. Generator inference and upstream candidate preparation were outside this timing experiment. Selection includes reading primary fingerprints from in-memory candidate geometries. The audit starts with parsed deposited geometries and preloaded models and references. The current auditor evaluates deposited energies and symmetry on the selected outputs. The external processing column excludes the common strain-acceptance check. Each preview uses 256 additional candidate relaxations per collection. Times vary with hardware, structure sizes, batching and model initialization. [Timing observations and environment](method-profile.json.gz)

## Paired changes in detection

Each cell gives **collections newly detected / previously detected collections lost**, relative to the original no-cell pipeline on the same draws. For preview versus original no-cell after CHGNet strain, the net gains were 8.2 percentage points on E and 17.4 on F. Conservative 95% paired intervals conditional on these pools were 4.1–12.1 and 12.5–21.9 percentage points. We separately bounded the gain and loss probabilities with 97.5% Clopper–Pearson intervals and subtracted the bounds.

| Processing | Default recovered / lost, E | Default recovered / lost, F | Preview recovered / lost, E | Preview recovered / lost, F |
|---|---:|---:|---:|---:|
| Ordinary CHGNet | 0 / 2 | 0 / 2 | 0 / 2 | 0 / 2 |
| Strain + CHGNet | 45 / 8 | 87 / 4 | 45 / 4 | 89 / 2 |
| Strain + MACE | 9 / 3 | 20 / 10 | 9 / 4 | 20 / 7 |

[Trial records for E](E-method-comparison.json.gz) · [Trial records for F](F-method-comparison.json.gz) · [All strained MACE endpoints for E](E-strain-pool.json.gz) · [All strained MACE endpoints for F](F-strain-pool.json.gz) · [Replay instructions](reproduce.md#original-and-current-methods-on-the-same-candidates)
