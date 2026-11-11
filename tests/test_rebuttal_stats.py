import sys
from pathlib import Path

import numpy as np
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from exchange_power import _poisson_binomial_tail  # noqa: E402
from reuse_exposure import exposure_law  # noqa: E402


def test_exposure_law_needs_ten_classes_at_half_and_one_per_mille():
    probability, needed = exposure_law(256, 0.04)
    assert needed == 10
    assert abs(probability - binom.sf(9, 256, 0.04)) < 1e-12


def test_exposure_law_is_monotone_in_observations_and_share():
    values = [exposure_law(n, 0.03)[0] for n in (64, 128, 256, 512)]
    assert values == sorted(values)
    assert exposure_law(256, 0.02)[0] < exposure_law(256, 0.04)[0]


def test_poisson_binomial_tail_matches_binomial_for_equal_probabilities():
    tail = _poisson_binomial_tail(np.full(32, 0.625), 29)
    assert abs(tail - binom.sf(28, 32, 0.625)) < 1e-12


def test_poisson_binomial_tail_handles_certain_and_impossible_blocks():
    assert _poisson_binomial_tail(np.array([1.0, 1.0, 0.0]), 2) == 1.0
    assert _poisson_binomial_tail(np.array([1.0, 1.0, 0.0]), 3) == 0.0
