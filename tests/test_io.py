import json
import math

import pytest

from basinmark.io import FingerprintCorpus, JsonStore, PairManifest


def test_json_store_round_trip(tmp_path):
    store = JsonStore(tmp_path / "nested" / "result.json")
    store.write({"value": 3})
    assert store.read_object() == {"value": 3}


def test_pair_manifest_resolves_relative_paths(tmp_path):
    structure = tmp_path / "structure.cif"
    structure.write_text("placeholder", encoding="utf-8")
    manifest = tmp_path / "pairs.json"
    manifest.write_text(
        json.dumps([{"raw": "structure.cif", "dft": "structure.cif"}]), encoding="utf-8"
    )
    records = PairManifest(manifest).load(("raw", "dft"))
    assert records[0]["raw"] == str(structure)


def test_fingerprint_corpus_rejects_blank_lines(tmp_path):
    corpus = tmp_path / "fingerprints.txt"
    corpus.write_text("F0\n\nF1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="nonempty fingerprint"):
        FingerprintCorpus(corpus).load()


def test_json_store_rejects_nonstandard_numbers(tmp_path):
    store = JsonStore(tmp_path / "result.json")
    with pytest.raises(ValueError):
        store.write({"value": math.nan})
    assert not store.path.exists()
