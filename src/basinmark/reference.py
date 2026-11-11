"""Scale-free local geometry summaries for comparing a deposited record to a reference."""

import numpy as np


def neighbor_profile(structure):
    """Mean log distances of the eight nearest periodic neighbors, by distance rank.

    Distances are normalized by the cube root of volume per atom. The result is
    invariant to uniform scale, atom order and equivalent supercell expansion.
    """
    length = (structure.volume / len(structure)) ** (1 / 3)
    for multiplier in (2, 3, 4, 6):
        centers, _, _, distances = structure.get_neighbor_list(multiplier * length)
        order = np.lexsort((distances, centers))
        counts = np.bincount(centers, minlength=len(structure))
        if counts.min() >= 8:
            starts = counts.cumsum() - counts
            local = distances[order[starts[:, None] + np.arange(8)]]
            if np.any(local <= 0):
                raise ValueError("neighbor distances must be positive")
            return np.log(local / length).mean(axis=0)
    raise ValueError("insufficient neighbors within the search radius")


def profile_distance(first, second):
    """RMS displacement of the eight log-distance ranks."""
    a, b = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if a.shape != (8,) or b.shape != (8,) or not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError("expected two finite eight-rank profiles")
    return float(np.sqrt(np.mean((a - b) ** 2)))
