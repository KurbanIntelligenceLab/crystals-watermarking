import json

import numpy as np
import pytest

from basinmark import audit, is_green, select_green
from basinmark.analyses import main
from basinmark.collections import matching_indices
from basinmark.core import critical_value


def test_bounded_selection_and_effective_sample_size():
    key = b"test"
    green = next(str(i) for i in range(100) if is_green(str(i), key))
    red = next(str(i) for i in range(100) if not is_green(str(i), key))
    fps = [red, green, green]
    assert select_green(fps, [0.0, 0.2, 0.1], key) == (2, True)
    assert select_green(fps, [0.0, 0.2, 0.1], key, energy_window=0.1) == (2, True)
    assert select_green(fps, [0.0, 0.2, 0.1], key, energy_window=0.05) == (0, False)
    result = audit([str(i) for i in range(7)] * 9, key)
    assert result.n_eff == pytest.approx(7)
    assert result.n_eff_population_estimate == pytest.approx(7.75)
    assert result.n_units == 7
    assert audit(["a", "a", "a", "b"], key).n_eff == pytest.approx(1.6)
    with pytest.raises(ValueError, match="energy_window"):
        select_green(fps, [0, 1, 2], key, energy_window=np.nan)


def test_reuse_exposes_all_fallbacks_and_separates_record_identities(tmp_path):
    fingerprints = tmp_path / "fps.txt"
    fingerprints.write_text("one\n" * 12)
    out = tmp_path / "audit.json"
    main(
        [
            "class-reuse",
            "--fingerprints",
            str(fingerprints),
            "--m",
            "2",
            "--observations",
            "3",
            "--trials",
            "20",
            "--out",
            str(out),
        ]
    )
    records = json.loads(out.read_text())["rows"]["3"]["trial_records"]
    assert {r[2] for r in records} == {0, 3}
    for _, selected, _, replay, _, _, _, disjoint, _, _, _ in records:
        assert len(set(selected)) == 3
        assert (replay, disjoint) == (12, 9)
    assert matching_indices(["f", "g", "f"], [0], True) == [2]


def test_collection_and_null_trials_replay(tmp_path):
    corpus = tmp_path / "corpus.json"
    rows = []
    for i in range(32):
        fp = {"full": f"{i}", "no_cell": f"{i // 2}"}
        rows.append(
            {
                "id": str(i),
                "energy": i / 100,
                "valid": True,
                "producer": fp,
                "stages": {"cif": fp} if i else {},
            }
        )
    corpus.write_text(
        json.dumps({"metadata": {"engine": None, "paired_dft_source": None}, "rows": rows})
    )
    path = tmp_path / "power.json"
    main(
        [
            "collection-benchmark",
            "--corpus",
            str(corpus),
            "--sizes",
            "8",
            "--energy-window",
            "0.05",
            "--m",
            "4",
            "--trials",
            "5",
            "--out",
            str(path),
        ]
    )
    report = json.loads(path.read_text())
    s = report["summary"]["8"]
    assert s["no_cell/bounded"]["max_penalty_meV_per_atom"] <= 50 + 1e-12
    assert s["full/minimum"]["max_penalty_meV_per_atom"] == 0
    assert s["full/minimum"]["stages"]["cif"]["record_coverage"] < 1
    null = tmp_path / "null.json"
    main(["key-controls", "--corpus", str(corpus), "--m", "4", "--keys", "12", "--out", str(null)])
    report = json.loads(null.read_text())
    batches = np.random.default_rng(report["seed"]).permutation(32).reshape(-1, 4)
    other_hits = 0
    for seed, k0, g0, k1, g1 in report["trial_records"]:
        rng = np.random.default_rng(seed)
        k, other = rng.bytes(16), rng.bytes(16)
        selected, unmarked = [], []
        for batch in batches:
            fps = [rows[i]["producer"]["full"] for i in batch]
            energies = [rows[i]["energy"] for i in batch]
            selected.append(fps[select_green(fps, energies, k)[0]])
            unmarked.append(fps[np.argmin(energies)])
        a, b = audit(unmarked, k), audit(selected, other)
        assert (k0, g0, k1, g1) == (a.n_units, a.green_count, b.n_units, b.green_count)
        other_hits += g1 >= critical_value(k1, 0.5, 0.001)
    assert other_hits == report["summary"]["independent_key"]["rejections"]


def test_prepared_energy_join_uses_source_id_and_restores_exponent():
    import runpy
    from pathlib import Path

    core = pytest.importorskip("pymatgen.core")
    path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_collection.py"
    module = runpy.run_path(str(path))

    def row(rid, energy_id, structure_species, formula, score, message):
        s = core.Structure(core.Lattice.cubic(4), [structure_species], [[0, 0, 0]])
        return {
            "material_id": rid,
            "str_ind": energy_id,
            "structure": json.dumps(s.as_dict()),
            "formula": formula,
            "chgnet_e_above_hull_corrected": score,
            "error_str": message,
        }

    raw = [
        row("0", "1", "Si", "Ge", "-3.4", "(e_h: -3.4e-05)"),
        row("1", "0", "Ge", "Si", "0.2", ""),
    ]
    joined, _, info = module["joined_candidates"]("CrystalFormer", raw)
    assert [(rid, e) for rid, e, _ in joined] == [("0", 0.2), ("1", -3.4e-05)]
    assert info["restored_scientific_exponents"] == 1
    raw[1]["formula"] = "C"
    with pytest.raises(ValueError, match="Composition mismatch"):
        module["joined_candidates"]("CrystalFormer", raw)


def test_compressed_measurement_output_round_trip(tmp_path):
    import gzip

    from basinmark.io import JsonStore

    path = tmp_path / "sample.json.gz"
    store = JsonStore(path)
    store.write({"k": [1, 2]})
    assert json.loads(gzip.decompress(path.read_bytes())) == store.read_object()
    before = path.read_bytes()
    with pytest.raises(ValueError):
        store.write({"k": np.nan})
    assert path.read_bytes() == before
