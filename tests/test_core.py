import math

import pytest

import basinmark


def test_public_api_audits_distinct_fingerprints():
    key = b"test-key"
    fingerprints = [f"F{i}" for i in range(40)]
    result = basinmark.audit(fingerprints, key)
    assert result.n_structures == 40
    assert result.n_distinct == 40
    assert result.n_units == 40
    assert math.isfinite(result.p_value)


def test_empty_audit_is_well_defined():
    result = basinmark.audit([], b"test-key")
    assert result.n_units == 0
    assert result.n_eff == 0.0
    assert result.design_effect == 0.0
    assert result.p_value == 1.0
    assert not result.fires


def test_selection_rejects_misaligned_inputs():
    with pytest.raises(ValueError, match="equal length"):
        basinmark.select_green(["F0"], [], b"test-key")


def test_known_design_points():
    assert basinmark.audit_size(m=8) == 10
    assert basinmark.n_forge() == 10
