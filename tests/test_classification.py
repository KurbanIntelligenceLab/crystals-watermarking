from pathlib import Path

import pytest

from basinmark.analyses import main
from basinmark.classification import ClassificationConfig, FlipComparison


def test_flip_comparison_partitions_changes():
    discrete = FlipComparison("A", "B", "A|1", "B|1").as_record()
    continuous = FlipComparison("A", "A", "A|1", "A|2").as_record()
    assert discrete["discrete_flip"] == 1
    assert discrete["continuous_only_flip"] == 0
    assert continuous["discrete_flip"] == 0
    assert continuous["continuous_only_flip"] == 1


def test_flip_comparison_rejects_inconsistent_fingerprints():
    with pytest.raises(ValueError, match="discrete fingerprint"):
        FlipComparison("A", "B", "same", "same").as_record()


def test_classification_config_rejects_nonpositive_settings():
    with pytest.raises(ValueError, match="must be positive"):
        ClassificationConfig(
            "chgnet",
            Path("in.json"),
            Path("pairs.json"),
            Path("cifs"),
            Path("out.json"),
            fmax=0,
        )


def test_cli_reports_version(capsys):
    with pytest.raises(SystemExit) as raised:
        main(["--version"])
    assert raised.value.code == 0
    assert capsys.readouterr().out.strip() == "basinmark 0.2.0"
