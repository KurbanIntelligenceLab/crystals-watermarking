# Replay and physical inputs

Run the commands from the repository root. The base package replays statistics using measured fingerprints, original uncertainties, deposited energy scores and the fixed producer preview pairs. Each comparison retains all 500 attempted collections and includes the selected record indices. The development groups and the two prospective comparisons have different roles. [Audit protocol](protocol.md) · [Producer protocol](producer-protocol.md)

~~~bash
python -m pip install -e '.[test]'
python scripts/fit_survival.py --verify
for cohort in A B first-confirm dev shear mattergen C D; do
  python scripts/benchmark_reliability.py --case "$cohort"
done
for cohort in A B first-confirm dev shear mattergen C D E F; do
  python scripts/benchmark_producer_preview.py --case "$cohort"
done
python scripts/check_preview_controls.py
for cohort in E F; do python scripts/benchmark_cross_model.py --case "$cohort"; done
python scripts/check_preview_controls.py --cross-model
python scripts/benchmark_bounded_removal.py --verify results/v0.6/bounded-removal.json.gz
for rule in both-pass preserve; do
  python scripts/benchmark_bounded_removal.py --phonons results/v0.6/emitted-removal-phonons.json.gz \
    --rule "$rule" --verify "results/v0.6/removal-$rule.json.gz"
done
(cd results/v0.6 && shasum -a 256 -c SHA256SUMS)
~~~

## Regenerate the measurements

Use the pinned source archives and weights in the [physical input guide](../v0.4/reproduce.md). As in the [earlier audit workflow](../v0.5/reproduce.md), PhononBench supplies the clean CrystalFormer inputs and original phonon labels. MACE-OMAT evaluates energies. CHGNet-r2SCAN supplies the separate variable-cell relaxation. Runtime fields depend on the machine. Statistical replay uses the distributed measurements and their file hashes. [PhononBench source](https://doi.org/10.5281/zenodo.18185662)

~~~bash
python -m pip install -e '.[benchmark,phononbench,test]'
src=data/archives
out=data/runs/fresh-comparison
mkdir -p "$out"
common=(--archive "$src/phonon.zip" --hull "$src/hull.parquet" --weights "$src/omat.model"
  --hull-solver linear --limit 768 --workers 5
  --exclude results/v0.3/measurements/CrystalFormer-corpus.json.gz
  --exclude results/v0.4/CrystalFormer-external.json.gz --exclude results/v0.4/CrystalFormer-second.json.gz
  --exclude results/v0.5/CrystalFormer-confirm.json.gz --exclude results/v0.5/CF-reserved.json.gz)
python scripts/prepare_phononbench.py "${common[@]}" --seed 826391 --out "$out/CF-reserved.json.gz"
python scripts/split_recovery_cohort.py --source "$out/CF-reserved.json.gz" --out "$out/CF-C.json.gz" "$out/CF-D.json.gz"
python scripts/prepare_phononbench.py "${common[@]}" --exclude "$out/CF-reserved.json.gz" --seed 589231 --out "$out/EF-reserved.json.gz"
python scripts/split_recovery_cohort.py --source "$out/EF-reserved.json.gz" --out "$out/CF-E.json.gz" "$out/CF-F.json.gz"
for cohort in C D E F; do
  python scripts/measure_symmetry_views.py --corpus "$out/CF-$cohort.json.gz" --directory "$src" \
    --archive "$src/phonon.zip" --weights "$src/omat.model" --tolerances .30 --workers 3 \
    --out "$out/CF-$cohort-scores.json.gz"
  python scripts/measure_producer_preview.py --corpus "$out/CF-$cohort.json.gz" --directory "$src" \
    --archive "$src/phonon.zip" --workers 2 --out "$out/$cohort-preview.json.gz"
done
~~~

The same preview command accepts the older archived input cohorts. Use the input paths in the replay script for A, B, the first confirmation, the original development cohort and MatterGen. The shear comparison uses the first confirmation producer geometries. Each preview deforms a copy of a producer candidate by 0.8%, relaxes it and reads its two symmetry fingerprints. Only the original candidate can be emitted.

## Changing the downstream relaxation model

The archived measurement files contain the deposited fingerprints, common energy, validity and convergence for all emitted IDs under both treatments. The base-package replay compares the full trial records. To repeat the physical work, use the original PhononBench archive and MACE-OMAT ranking weights from the input guide, plus the pinned MACE-MPA-0 medium relaxation weights. The measurement command resumes its single output file when the source files, model hashes, recorded software versions and identities match. Reuse of previously computed endpoints across the 500 releases is intentional. [Frozen comparison](cross-model-protocol.md) · [Foundation model source](https://github.com/ACEsuit/mace-mp/releases/tag/mace_mpa_0)

~~~bash
curl -fL https://github.com/ACEsuit/mace-mp/releases/download/mace_mpa_0/mace-mpa-0-medium.model -o "$src/mpa.model"
for cohort in E F; do
  python scripts/measure_cross_model.py --case "$cohort" --archive "$src/phonon.zip" \
    --relax-weights "$src/mpa.model" --rank-weights "$src/omat.model" --workers 2 \
    --out "$out/$cohort-cross-model-physics.json.gz"
  python scripts/benchmark_cross_model.py --case "$cohort" --physics "$out/$cohort-cross-model-physics.json.gz" \
    --out "$out/$cohort-cross-model.json.gz"
done
~~~

## Active removal and endpoint quality

The attacker receives selected producer geometries and can test up to eight complete collections against a binary detector. The first query tests ordinary relaxation, the second applies all eligible shear edits, and the others apply prespecified subsets (changed primary, changed recovery, or one of four reproducible subsets). Queries stop at the first negative decision. No marking key enters proposal construction. Each edited record must converge, pass the minimum-distance check, and stay within 10 meV/atom both immediately after deformation and relative to its independently relaxed ordinary endpoint at deposition. Endpoints unavailable to the ordinary audit are kept identical across proposals. Every selected record stays in its collection.

To regenerate the paired harmonic screen in a compatible physics environment, use Phonopy 2.43.6, symfc 1.5.3 and spglib 2.6.0. CHGNet-r2SCAN supplies forces on the reconstructed deposited endpoints, checked against both archived fingerprint readings. Use central displacements of 0.01 Å in a 2 by 2 by 2 supercell, remove drift, symmetrize the force constants and evaluate a Gamma-centered 4 by 4 by 4 grid. A passing endpoint has minimum signed frequency at least −0.3 THz. This is a finite-model, finite-grid screen. The two guarded attacks either require both endpoints to pass or require an originally passing endpoint to keep passing. [Phonopy formulation](https://phonopy.github.io/phonopy/formulation.html)

~~~bash
python scripts/measure_removal_phonons.py --archive "$src/phonon.zip" \
  --removal results/v0.6/bounded-removal.json.gz --out "$out/phonons.json.gz" --workers 2
~~~

## Cross-catalogue reuse

For the paired 256-observation comparison, place the pinned MP-20 and OQMD files in the source directory. The command records the selected exposure, original source hashes, matcher settings, and all 200 trials for each fingerprint. The earlier command without an exposure argument reproduces the 1,130-observation setting. [Catalogue sources](../v0.3/README.md)

~~~bash
python scripts/catalogue_overlap.py --directory "$src" --observations 256 --out "$out/catalogue-256.json.gz"
~~~

## Reproduce the class-reuse panel

Each producer emits one record from every batch of eight candidates. It chooses the lowest-energy green candidate when one is available and otherwise emits the batch minimum. The observer receives each emitted ID in order and can read its full structural fingerprint. It samples 200 matching records with replacement from the same generator corpus. One arm allows copying an exposed record, and the second excludes exposed IDs. Empty pools count as unsuccessful trials. The producer uses the key for selection, and the auditor uses it for detection. The observer runs without the key or fallback labels.

Run the statistical replay to reconstruct the archived trials from the measured corpora. Then render the panel into a new directory. The export checks every emitted index and fallback count against the trials before writing the ordered IDs and one-based fallback positions. Look up an ID in the corresponding generator’s corpus to read its full fingerprint. The reference SVG was drawn with Matplotlib 3.11.2. [Panel](class-reuse.svg) · [Observation records](reuse-observations.json.gz) · [MatterGen corpus](../v0.3/measurements/MatterGen-corpus.json.gz) · [CrystalFormer corpus](../v0.3/measurements/CrystalFormer-corpus.json.gz)

~~~bash
python -m pip install matplotlib==3.11.2
python scripts/replay_benchmarks.py
out=$(mktemp -d)
python scripts/replot_reuse.py --out "$out/class-reuse.svg" --observations-out "$out/reuse-observations.json.gz"
~~~

At ten observations, each generator contributed 2,000 emissions in 200 trials. MatterGen had six fallbacks and CrystalFormer had three. Each occurred in a different trial within its generator. Copying an exposed record was allowed in the arm that detected 194/200 MatterGen and 197/200 CrystalFormer collections. Excluding exposed IDs gave 0/200 for each. For a concrete example, in MatterGen trial `1251385316631358` the fourth emitted record was ID `578`. It was a fallback. The observer saw this record along with the other nine. [MatterGen trial outcomes](../v0.3/trials/MatterGen-reuse.json.gz) · [CrystalFormer trial outcomes](../v0.3/trials/CrystalFormer-reuse.json.gz)
