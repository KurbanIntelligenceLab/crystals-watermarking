import itertools
import json

import numpy as np
import pytest

from basinmark import is_green, select_collection
from basinmark.analyses import main


def test_pool_optimizes_under_both_ceilings():
    rng = np.random.default_rng(722)
    for _ in range(40):
        fps = [f"class-{i}" for i in rng.integers(0, 7, 12)]
        e = rng.integers(-3, 6, 12) / 10
        key = rng.bytes(16)
        green = np.array([is_green(f, key) for f in fps])
        selected = select_collection(fps, e, key, m=4, policy="pool")
        baseline = select_collection(fps, e, key, m=4, policy="fixed")
        ceilings = np.minimum(np.sort(e)[:3] + 0.1, np.sort(e[baseline.indices]))
        candidates = (np.array(c) for c in itertools.combinations(range(12), 3))
        scores = [
            (sum(green[c]), -sum(e[c])) for c in candidates if np.all(np.sort(e[c]) <= ceilings)
        ]
        expected = max(scores)
        i = selected.indices
        assert (sum(green[i]), -sum(e[i])) == pytest.approx(expected)
        assert len(set(i)) == 3 and np.all(e[i] <= ceilings)


def test_reserve_counts_every_candidate_and_honors_final_minimum():
    fps = [str(i) for i in range(40)]
    energy = np.random.default_rng(198).random(40)
    for key in (b"a", b"b"):
        s = select_collection(fps, energy, key, policy="reserve")
        assert len(set(s.indices)) == 5
        assert min(s.batch_sizes) >= 4 and sum(s.batch_sizes) == 40
        assert np.all(energy[s.indices] <= s.reference_energies + 0.1)
    with pytest.raises(ValueError):
        select_collection(fps[:-1], energy[:-1], b"a")
    with pytest.raises(ValueError):
        select_collection(fps, np.full(40, np.nan), b"a")


def test_collection_cli_keeps_key_out_of_output(tmp_path):
    records = [{"id": str(i), "energy": i / 10, "fingerprint": str(i)} for i in range(16)]
    candidates, key, out = [tmp_path / s for s in ("c.json", "key", "result.json")]
    candidates.write_text(json.dumps({"candidates": records}))
    key.write_bytes(b"private-example-key")
    main(
        [
            "select-collection",
            "--key-file",
            str(key),
            "--candidates",
            str(candidates),
            "--out",
            str(out),
        ]
    )
    assert len(json.loads(out.read_text())["selected"]) == 2
    assert b"private-example-key" not in out.read_bytes()


def test_reference_profile_survives_equivalent_descriptions():
    pmg = pytest.importorskip("pymatgen.core")
    from basinmark.reference import neighbor_profile, profile_distance

    s = pmg.Structure(pmg.Lattice.cubic(4), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    reference = neighbor_profile(s)
    s.make_supercell([2, 1, 1])
    s.translate_sites(range(len(s)), [0.12, 0.03, -0.09])
    s.scale_lattice(s.volume * 1.1)
    assert profile_distance(reference, neighbor_profile(s)) < 1e-10


def test_budget_default_maximizes_green_under_slot_and_mean_limits():
    rng = np.random.default_rng(123)
    for _ in range(30):
        e = rng.integers(-6, 10, 12) / 16
        fps = list(map(str, rng.integers(0, 10, 12)))
        key = rng.bytes(16)
        marked = np.array([is_green(f, key, 0.625) for f in fps])
        fixed = select_collection(fps, e, key, m=4, gamma=0.625, policy="fixed")
        chosen = select_collection(fps, e, key, m=4)
        cap = np.sort(e.reshape(3, 4).min(axis=1)) + 0.1
        budget = min(e[fixed.indices].mean(), np.sort(e)[:3].mean() + 0.1)
        feasible = [
            np.array(c)
            for c in itertools.combinations(range(12), 3)
            if np.all(np.sort(e[list(c)]) <= cap) and e[list(c)].mean() <= budget
        ]
        score = max((sum(marked[c]), -sum(e[c])) for c in feasible)
        assert (sum(marked[chosen.indices]), -sum(e[chosen.indices])) == pytest.approx(score)
        assert np.all(e[chosen.indices] <= cap)
        assert e[chosen.indices].mean() <= budget
