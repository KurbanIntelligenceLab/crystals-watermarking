# Reproducing the allocation experiments

The statistical replay uses Python, NumPy and SciPy. The archive includes measured fingerprints and source IDs. The second PhononBench corpus also records each MACE total and composition-matched reference energy separately. For physical regeneration install the benchmark and phononbench extras. Each measurement and calibration output is written to a separate file. Use an empty output directory for new runs.

## Public inputs

The generated-crystal files and MP-20 splits come from the generated-crystal benchmark by Kazeev and collaborators. The PhononBench CrystalFormer-MP20 relaxed structures and dynamical stability labels come from Xiao-Qi Han. Both archives are distributed under CC BY 4.0. The external ranking uses a pinned MACE-OMAT-0-medium checkpoint and the matching energy column in the LeMat-Bulk MACE-OMAT reference. Check these providers' terms for redistributing the raw sources. The release distributes derived observations. [Generated-crystal benchmark](https://doi.org/10.6084/m9.figshare.29145101.v6) · [PhononBench](https://doi.org/10.5281/zenodo.18185662) · [LeMat reference snapshot](https://huggingface.co/datasets/LeMaterial/LeMat-Bulk-MLIP-Hull/tree/70d505bb294c16658a9551e2b81ace69a28d9790) · [MACE checkpoint release](https://github.com/ACEsuit/mace-mp/releases/tag/mace_omat_0)

```bash
python -m pip install -e '.[benchmark,phononbench,test]'
mkdir -p data/archives data/runs/v04
for id in 65498775 65498808 65498871 65498655 65498670; do
  curl -fL "https://ndownloader.figshare.com/files/$id" -o "data/archives/$id.csv.gz"
done
curl -fL https://zenodo.org/api/records/18185662/files/PhononBench_relaxed_all.zip/content -o data/archives/phonon.zip
curl -fL https://huggingface.co/datasets/LeMaterial/LeMat-Bulk-MLIP-Hull/resolve/70d505bb294c16658a9551e2b81ace69a28d9790/data/mace_omat-00000-of-00001.parquet -o data/archives/hull.parquet
curl -fL https://github.com/ACEsuit/mace-mp/releases/download/mace_omat_0/mace-omat-0-medium.model -o data/archives/omat.model
```

## Producer and downstream observations

The first MatterGen sample excludes the two earlier physics cohorts. The second also excludes the first. The second external cohort excludes all previous external IDs and classes. The second run evaluates each target hull composition with a linear program to avoid enumerating high-dimensional phase-diagram facets. The physics measurements include source and checkpoint hashes, dependency versions, failures and attempt counts. The source archives contain generator outputs. This procedure starts from those coordinates. Timing measurements vary by machine.

```bash
old=results/v0.3/measurements
out=data/runs/v04
for pass in fresh second; do
  seed=47017
  extra=()
  if [ "$pass" = second ]; then
    seed=63593
    extra=(--exclude-cohort "$out/MatterGen-fresh.json.gz")
  fi
  python scripts/prepare_collection.py --dataset MatterGen --directory data/archives \
    --limit 384 --seed "$seed" --paired-dft --engine chgnet-r2scan --workers 3 \
    --exclude-cohort "$old/MatterGen-512.json.gz" --exclude-cohort "$old/MatterGen-384.json.gz" \
    "${extra[@]}" --out "$out/MatterGen-$pass.json.gz"
done
python scripts/prepare_collection.py --dataset WyCryst --directory data/archives --limit 384 \
  --seed 47017 --engine chgnet-r2scan --workers 3 --out "$out/WyCryst-transfer.json.gz"

for pass in external second; do
  seed=47017
  extra=()
  if [ "$pass" = second ]; then
    seed=63593
    extra=(--exclude "$out/CrystalFormer-external.json.gz" --hull-solver linear)
  fi
  python scripts/prepare_phononbench.py --archive data/archives/phonon.zip \
    --hull data/archives/hull.parquet --weights data/archives/omat.model \
    --exclude "$old/CrystalFormer-corpus.json.gz" "${extra[@]}" --seed "$seed" \
    --workers 5 --out "$out/CrystalFormer-$pass.json.gz"
done
```

## Audits and optional geometry comparison

The following Bash commands reproduce the second-cohort audits from the regenerated observations. The replay driver specifies each development, retrospective, sensitivity and calibration input. Passing an existing trial file to the benchmark verification option checks every result field. Raw physical timings are outside exact replay.

```bash
for gen in MatterGen CrystalFormer; do
  python scripts/benchmark_allocation.py --corpus "$out/$gen-second.json.gz" \
    --trials 500 --seed 91571 --followup --out "$out/$gen-second-audits.json.gz"
  python scripts/calibrate_allocation.py --corpus "$out/$gen-second.json.gz" \
    --audits "$out/$gen-second-audits.json.gz" --seed 81689 --out "$out/$gen-second-null.json.gz"
done
for split in test validation; do
  python scripts/benchmark_reference.py --split "$split" --directory data/archives \
    --out "$out/MP20-$split-reference.json.gz"
done
python scripts/profile_allocation.py --corpus "$out/MatterGen-fresh.json.gz" \
  --out "$out/selection-cost.json.gz"
python scripts/replay_allocation.py
```
