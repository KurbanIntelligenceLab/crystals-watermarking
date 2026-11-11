from itertools import product

import pytest

from basinmark.reliability import membership, reliability_region, reliable_audit, reliable_classes
from basinmark.symmetry_audit import _rank_region, _region, ranked_audit


@pytest.mark.parametrize("gamma,alpha", [(0.625, 0.001), (0.625, 0.16), (0.37, 0.11)])
def test_reliability_region_against_all_label_assignments(gamma, alpha):
    for primary, secondary, ranked, reliable in [
        ({0, 1, 2, 3, 4}, {2, 3, 4, 5, 6}, {2, 5}, {3, 4, 5}),
        ({0, 1, 2}, {3, 4}, set(), {3, 4}),
        ({0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2}, {1, 3}),
        (set(), set(), set(), set()),
    ]:
        k, j, shared = len(primary), len(secondary), len(primary & secondary)
        c, d, t, _ = _region(k, j, shared, gamma, alpha)
        r, _ = _rank_region(k, j, shared, len(ranked), len(primary & ranked), gamma, alpha)
        cut, size = reliability_region(
            membership(primary, secondary, ranked, reliable), gamma, alpha
        )
        mass = 0.0
        for labels in product((0, 1), repeat=len(primary | secondary)):
            g, h, u, v = (sum(labels[i] for i in s) for s in [primary, secondary, ranked, reliable])
            if g >= c or (h, g) >= (d, t) or (u, h, g) >= r or (v, u, h, g) >= cut:
                mass += gamma ** sum(labels) * (1 - gamma) ** (len(labels) - sum(labels))
        assert size == pytest.approx(mass, abs=1e-12)
        assert mass <= alpha + 1e-12


def test_duplicate_rank_and_missing_energy_compatibility():
    rows = [(f"SG40|{i}", f"SG40|{i}" if i % 5 else None, float(i)) for i in range(32)]
    key = b"test"
    result = reliable_audit(rows, key)
    assert reliable_audit(rows + list(reversed(rows)), key) == result
    assert result.reliable_classes == len(reliable_classes(rows)) == 18
    assert result.ranked_fires == ranked_audit(rows, key).fires
    for sample in [rows + [("SG41|X", None, 0.0)], [("SG4|X", "SG4|X", None)]]:
        protected = ranked_audit(sample, key)
        decision = reliable_audit(sample, key)
        assert decision.fires == protected.fires and decision.null_rate == protected.null_rate
