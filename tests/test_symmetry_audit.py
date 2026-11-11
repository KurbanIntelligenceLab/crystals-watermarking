import itertools

import numpy as np
import pytest

from basinmark.collections import cutoff
from basinmark.symmetry_audit import _region, paired_audit


@pytest.mark.parametrize("gamma", [0.375, 0.625])
@pytest.mark.parametrize("alpha", [0.001, 0.08, 0.25])
def test_paired_region_against_exhaustive_null(gamma, alpha):
    for k, j, overlap in [(0, 0, 0), (1, 3, 0), (4, 3, 1), (4, 5, 3), (5, 5, 5), (5, 0, 0)]:
        c, d, t, size = _region(k, j, overlap, gamma, alpha)
        assert c == cutoff(k, gamma, alpha)
        observed = 0.0
        complement = {}
        for labels in itertools.product((0, 1), repeat=k + j - overlap):
            g, h = sum(labels[:k]), sum(labels[k - overlap :])
            p = gamma ** sum(labels) * (1 - gamma) ** (len(labels) - sum(labels))
            if g >= c or (h, g) >= (d, t):
                observed += p
            else:
                complement[h, g] = complement.get((h, g), 0.0) + p
        assert observed == pytest.approx(size)
        assert size <= alpha
        if complement:
            assert observed + complement[max(complement)] > alpha - 1e-12


def test_duplicate_and_missing_views_keep_legacy_detections():
    a = [str(i) for i in range(32)]
    for key in [np.random.default_rng(i).bytes(16) for i in range(12)]:
        x = paired_audit(a, a, key)
        assert x == paired_audit(a * 2, reversed(a), key)
        assert x.primary_fires == x.fires
        assert paired_audit(a, [], key).fires == x.primary_fires
    assert not paired_audit([], [], b"k").fires


def test_deposit_command_uses_frozen_views_without_storing_key(tmp_path):
    import json
    from dataclasses import asdict

    from basinmark.analyses import main
    from basinmark.symmetry_audit import audit_deposited, eligible_recovery

    pmg = pytest.importorskip("pymatgen.core")
    s = pmg.Structure(pmg.Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    cif = tmp_path / "deposited.cif"
    s.to(filename=cif)
    keyfile, out = tmp_path / "secret", tmp_path / "audit.json"
    keyfile.write_bytes(b"private-test-key")
    main(
        [
            "audit-deposit",
            "--structures",
            str(cif),
            str(cif),
            "--key-file",
            str(keyfile),
            "--out",
            str(out),
        ]
    )
    expected = audit_deposited([pmg.Structure.from_file(cif)] * 2, keyfile.read_bytes())
    result = json.loads(out.read_text())
    assert result["structures"] == 2 and result["shared_classes"] == 1
    assert {k: result[k] for k in asdict(expected)} == asdict(expected)
    assert keyfile.read_bytes() not in out.read_bytes()
    assert eligible_recovery(None) is None
    assert eligible_recovery("SG2|A|A:-1:1") is None
    with pytest.raises(ValueError):
        eligible_recovery("bad")


@pytest.mark.parametrize("geometry", [(5, 4, 3, 2, 1), (3, 4, 0, 3, 0), (4, 4, 4, 2, 2)])
def test_energy_subset_joint_size_by_exhaustive_enumeration(geometry):
    from basinmark.symmetry_audit import _rank_region

    k, j, overlap, top, shared_top = geometry
    gamma, alpha = 0.625, 0.16
    c, d, tie, _ = _region(k, j, overlap, gamma, alpha)
    threshold, size = _rank_region(*geometry, gamma, alpha)
    mass = 0.0
    b = list(range(overlap)) + list(range(k, k + j - overlap))
    subset = b[:shared_top] + b[overlap : overlap + top - shared_top]
    for labels in itertools.product((0, 1), repeat=k + j - overlap):
        g = sum(labels[:k])
        h = sum(labels[i] for i in b)
        u = sum(labels[i] for i in subset)
        if g >= c or (h, g) >= (d, tie) or (u, h, g) >= threshold:
            mass += gamma ** sum(labels) * (1 - gamma) ** (len(labels) - sum(labels))
    assert size == pytest.approx(mass)
    assert size <= alpha


def test_ranked_classes_ignore_order_and_duplicate_records():
    from basinmark.symmetry_audit import ranked_audit, upper_energy_classes

    rows = [(str(i), str(i) if i % 3 else None, float(i)) for i in range(12)]
    assert upper_energy_classes(rows) == {"4", "5", "7", "8", "10", "11"}
    a = ranked_audit(rows, b"key")
    assert a == ranked_audit(rows[::-1] + rows, b"key")
    assert ranked_audit([(p, s, None) for p, s, _ in rows], b"key").fires == a.geometric.fires
    many = [(str(i), str(i), float(i)) for i in range(45)]
    assert (
        ranked_audit(many, b"key").fires
        == paired_audit([p for p, _, _ in many], [s for _, s, _ in many], b"key").fires
    )
