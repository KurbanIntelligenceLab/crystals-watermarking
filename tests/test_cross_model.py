from pathlib import Path

import pytest

from basinmark.io import JsonStore


def test_cross_model_preserves_failed_endpoints_and_release_denominators(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    from benchmark_cross_model import MPA, benchmark
    from benchmark_producer_preview import sha

    def save(stem, value):
        p = tmp_path / (stem + ".json.gz")
        JsonStore(p).write(value)
        return p

    stages = ("rerelaxed", "strain_rerelaxed")
    rows = [{"id": str(i), "energy": float(i), "producer": {"no_cell": str(i)}} for i in range(2)]
    source = save("CF-E", {"metadata": {"source_sha256": "archive"}, "rows": rows})
    scores = save(
        "CF-E-scores",
        {
            "input_sha256": sha(source),
            "weights_sha256": "rank",
            "rows": [{"id": str(i), "scores": {s: 1.0 for s in stages}} for i in range(2)],
        },
    )
    trials = [
        {
            "seed": i,
            "control_seed": i + 10,
            "selected": [[0], [1]],
            "stages": {s: [0, 0, 0, 0] for s in stages},
        }
        for i in range(2)
    ]
    releases = save(
        "E-producer",
        {
            "source_sha256": sha(source),
            "views_sha256": sha(scores),
            "emissions": 1,
            "controls": 2,
            "records": trials,
        },
    )
    measured = {
        "metadata": {
            "case": "E",
            "archive_sha256": "archive",
            "rank_weights_sha256": "rank",
            "relax_weights_sha256": MPA,
            "fmax": 0.05,
            "steps": 150,
            "source_sha256": sha(source),
            "producer_sha256": sha(releases),
            "indices": [0, 1],
        },
        "records": [
            {"id": "0", "index": 0, "stages": {s: {"converged": False} for s in stages}},
            {
                "id": "1",
                "index": 1,
                "stages": {
                    s: {
                        "converged": True,
                        "primary": "SG10|test",
                        "secondary": None,
                        "energy": 1.002,
                        "valid": True,
                    }
                    for s in stages
                },
            },
        ],
    }
    physics = save("measurements", measured)
    result = benchmark("E", tmp_path, physics)
    for st in stages:
        base, new = (result["summary"][st][p] for p in ("budget", "preview"))
        assert base["alternate"]["trials"] == new["alternate"]["trials"] == 2
        assert base["wrong_tests"] == new["wrong_tests"] == 4
        assert base["attempted_emissions"] == 2 and base["converged_emissions"] == 0
        assert base["mean_paired_energy_change_meV_per_atom"] is None
        assert new["valid_converged_emissions"] == new["jointly_scored_emissions"] == 2
        assert new["mean_paired_energy_change_meV_per_atom"] == pytest.approx(2.0)
    measured["records"].pop(0)
    JsonStore(physics).write(measured)
    with pytest.raises(ValueError, match="coverage"):
        benchmark("E", tmp_path, physics)


def test_method_comparison_aligns_original_fingerprints_and_keeps_missing_deposits(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    from compare_method_versions import endpoint_tables

    first = {"full": "with-cell", "no_cell": "SG15|x"}
    rows = [
        {"id": "a", "stages": {"strain_rerelaxed": first}},
        {"id": "b", "stages": {}},
    ]
    scores = {
        "tolerances": [0.3],
        "rows": [
            {"id": "a", "stages": {"strain_rerelaxed": ["SG2|x"]},
             "scores": {"strain_rerelaxed": 0.20}},
            {"id": "b", "stages": {}, "scores": {}},
        ],
    }
    tables = endpoint_tables(None, "E", rows, scores, None)
    assert tables["CHGNet_rerelaxed"] == [None, None]
    assert tables["CHGNet_strain_rerelaxed"] == [("with-cell", "SG15|x", None, 0.20), None]
    scores["rows"].reverse()
    with pytest.raises(ValueError, match="identities differ"):
        endpoint_tables(None, "E", rows, scores, None)
