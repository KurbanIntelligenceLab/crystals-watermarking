# Statistical follow-ups on archived measurements

Three additions that use only archived measurements, public catalogues and the shipped method. None changes a claim or a contribution. Scripts live in `scripts/`, outputs in this directory.

## 1. Reuse-exposure law

**Question.** Does removing cell shape make class-reuse attacks easier, and by how much?

**Result.** Cross-catalogue transfer (observer sees N genuine marked MP-20 outputs, resamples 200 matching OQMD records) follows a sigmoid in N for both fingerprints. The no-cell fingerprint shifts it left by a factor of about 1.6, the ratio of the two matched shares. Both still need hundreds of observations. Excluding external records geometrically equivalent to an exposed structure gives 0/200 at every N for both fingerprints.

| Observed outputs N | Full: observed | Full: law | No cell: observed | No cell: law |
|---:|---:|---:|---:|---:|
| 128 | 0.0% | 0.1% | 1.0% | 3.5% |
| 192 | 0.5% | 2.3% | 19.0% | 25.0% |
| 256 | 11.5% | 10.8% | 43.0% | 58.7% |
| 384 | 35.0% | 48.5% | 93.0% | 94.8% |
| 512 | 76.5% | 81.9% | 99.5% | 99.7% |
| 768 | 99.0% | 99.2% | 100% | 100% |

200 trials per point; exact 95% intervals in the archive. "Law" is P{Bin(N, s) >= 10} with s the catalogue-matched share among the lowest-energy record of each key-free batch of eight (s = 0.0246 with cell shape, 0.0404 without). Ten is the smallest all-green class count that rejects at alpha = 0.001 and gamma = 1/2. The population share (0.0182 and 0.0295) gives a lower curve; the empirical curve sits between the two within sampling error. The gap is duplicate classes among matched records.

**Why the emitted share exceeds the population share.** Low-energy structures are more often catalogued. By energy decile of MP-20 test, the matched share falls from 4.2% in the cheapest decile to 1.9% in the most expensive (no cell). Energy-preferential selection therefore enriches exposure by about 1.35x. This is a property of the ranking, not of the key.

Files: `scripts/reuse_exposure.py`, `scripts/plot_reuse_exposure.py`, `results/stats/reuse-exposure.json.gz`, `results/stats/reuse-exposure.svg`, cache `results/stats/catalogue-cache.json.gz`.

## 2. Selection power for the shipped selectors

**Question.** The paper's power formula covers independent fixed batches. What governs the fixed capped selector and the pooled budget selector?

**Result.** For the fixed capped selector, the green count is a Poisson-binomial sum over blocks with success probability 1 - (1 - gamma)^{a_i}, where a_i is the number of distinct classes within the window of block i's minimum. This predicts producer-stage power without fitting:

| Cohort | Blocks with a_i = 1 | Predicted power | Observed power |
|---|---:|---:|---:|
| CrystalFormer E | 46.9% | 3.8% | 2.8% |
| CrystalFormer F | 55.1% | 1.9% | 1.6% |
| CrystalFormer second (0.4) | 48.4% | 3.2% | 4.4% |
| MatterGen second (0.4) | 4.8% | 98.3% | 99.0% |

For the pooled budget selector, slot j may draw any class below the j-th sorted block minimum plus the window. The supply condition (at least j green classes under the j-th ceiling for every j) held in all 2,000 trials. In every trial where the budget rule emitted fewer than 32 greens, the binding constraint was the mean-energy ceiling, never a missing green class or a per-slot ceiling. Observed producer-stage power 96.8%, 98.0%, 98.4%, 100%.

Files: `scripts/exchange_power.py`, `results/stats/exchange-power-{E,F,CF2,MG2}.json.gz`, .

## 3. Two-view calibration as a stated proposition

The recovery audit's level control is implemented in `basinmark.symmetry_audit` and checked against exhaustive enumeration of the joint null in `tests/test_symmetry_audit.py`.

