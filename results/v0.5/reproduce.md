# Reproducing the symmetry audit

The statistics replay uses the base Python package and the archived measurements. Every attempt remains in the trial denominator. Producer selection, quality labels and energy limits are inherited from the 0.4 budget policy. The producer sees no downstream observation. [Design, cohorts and locked rules](protocol.md)

```bash
python -m pip install -e '.[test]'
python scripts/replay_symmetry_audit.py
python scripts/replay_allocation.py
python scripts/replay_benchmarks.py --all
python -m pytest
(cd results/v0.5 && shasum -a 256 -c SHA256SUMS)
```

## Regenerating physics

Install the benchmark and phononbench extras. Download the pinned generated-crystal, PhononBench, MACE-OMAT and LeMat archives from the 0.4 input guide. Use Bash, an empty output directory and the following commands. Worker timings are machine-dependent, so regenerated physics files can have different byte digests. The statistics replay verifies the distributed inputs and outcomes exactly. The energy model and the convex-hull reference must match. [Public source URLs and file identities](../v0.4/reproduce.md)

```bash
src=data/archives
out=data/runs/v05
old=results/v0.4
past=results/v0.3/measurements
mkdir -p "$out"
python -m pip install -e '.[benchmark,phononbench,test]'
phonon=(--archive "$src/phonon.zip" --hull "$src/hull.parquet" --weights "$src/omat.model"
        --exclude "$past/CrystalFormer-corpus.json.gz" --exclude "$old/CrystalFormer-external.json.gz"
        --exclude "$old/CrystalFormer-second.json.gz" --hull-solver linear)
python scripts/prepare_phononbench.py "${phonon[@]}" --seed 119233 --limit 384 --workers 5 \
  --out "$out/CrystalFormer-confirm.json.gz"
python scripts/prepare_phononbench.py "${phonon[@]}" --seed 201773 --limit 768 --workers 6 \
  --exclude "$out/CrystalFormer-confirm.json.gz" --out "$out/CF-reserved.json.gz"
python scripts/split_recovery_cohort.py --source "$out/CF-reserved.json.gz" \
  --out "$out/CF-A.json.gz" "$out/CF-B.json.gz"
python scripts/prepare_collection.py --dataset MatterGen --directory "$src" --limit 384 \
  --seed 119233 --workers 3 --paired-dft --engine chgnet-r2scan \
  --exclude-cohort "$past/MatterGen-512.json.gz" --exclude-cohort "$past/MatterGen-384.json.gz" \
  --exclude-cohort "$old/MatterGen-fresh.json.gz" --exclude-cohort "$old/MatterGen-second.json.gz" \
  --out "$out/MatterGen-confirm.json.gz"
for half in A B; do
  seed=931153; nullseed=486211
  if [ "$half" = B ]; then seed=957203; nullseed=486223; fi
  python scripts/measure_symmetry_views.py --corpus "$out/CF-$half.json.gz" --directory "$src" \
    --archive "$src/phonon.zip" --weights "$src/omat.model" --tolerances .3 --workers 5 \
    --out "$out/CF-$half-scores.json.gz"
  python scripts/benchmark_symmetry_audit.py --corpus "$out/CF-$half.json.gz" \
    --views "$out/CF-$half-scores.json.gz" --seed "$seed" --control-seed "$nullseed" \
    --out "$out/CF-$half-audits.json.gz"
done
python scripts/measure_shear.py --corpus "$out/CrystalFormer-confirm.json.gz" --archive "$src/phonon.zip" \
  --weights "$src/omat.model" --score --out "$out/CF-shear-scored.json.gz"
python scripts/benchmark_symmetry_audit.py --corpus "$out/CF-shear-scored.json.gz" \
  --out "$out/CF-shear-ranked.json.gz"
```

The first confirmation uses the same measurement and audit commands with seed 734201 and control seed 486209. For fresh MatterGen supply the matching hull with `--hull` when ranking deposited structures. For geometry-only audits omit the energy model, specify the recorded tolerances and pass the resulting view file. The design sweep measured 0.15, 0.20, 0.30, 0.40, 0.60 and 0.80 Å on the two earlier external cohorts. For structural resolution, pass the corresponding 0.4 matcher pairs to the command below. Timing is measured over 64 source structures following a separate initialization and warm-up. The one-thread CPU profile includes the deposited-energy single point and the composition hull query.

```bash
for split in test validation; do
  python scripts/compare_symmetry_resolution.py --directory "$src" --split "$split" \
    --reference "$old/MP20-$split-reference.json.gz" --out "$out/MP20-$split-resolution.json.gz"
done
python scripts/profile_deposit.py --corpus "$old/CrystalFormer-second.json.gz" \
  --archive "$src/phonon.zip" --weights "$src/omat.model" --hull "$src/hull.parquet" \
  --out "$out/deposit-cost.json.gz"
```

The secondary release size is 40. Set the size and candidates per output separately for each run. Use seed 734201 for the first fresh cohort, 931153 for A and 957203 for B. This size comparison uses the geometry-only test and records the release identities and physical quality values.

```bash
for m in 6 8; do
  python scripts/benchmark_symmetry_audit.py --corpus "$out/CF-A.json.gz" \
    --views "$out/CF-A-scores.json.gz" --seed 931153 --emissions 40 --m "$m" \
    --controls 0 --geometry-only --out "$out/size-A-$m.json.gz"
done
```
