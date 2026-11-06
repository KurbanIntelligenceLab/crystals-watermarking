"""A paired, deduplicated class audit for two views of a deposited collection."""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.stats import binom

from .collections import cutoff
from .core import is_green


@dataclass(frozen=True)
class PairedAudit:
    primary_classes: int
    secondary_classes: int
    shared_classes: int
    primary_green: int
    secondary_green: int
    primary_cutoff: int
    secondary_cutoff: int
    primary_tiebreak: int
    null_rate: float
    primary_fires: bool
    fires: bool


@lru_cache(maxsize=512)
def _region(k, j, shared, gamma, alpha):
    """Protect the primary rejection region, then fill by secondary count, primary count."""
    first = cutoff(k, gamma, alpha)
    # Conditional on h shared greens, the two exclusive counts are independent.
    # A one-dimensional mixture avoids allocating a joint table for large audits.
    h = np.arange(shared + 1)
    mass = binom.pmf(h, shared, gamma)
    old = float(binom.sf(first - 1, k, gamma))
    miss = binom.cdf(first - 1 - h, k - shared, gamma)

    def first_feasible(upper, null_rate):
        low, high = 0, upper
        while low < high:
            mid = (low + high) // 2
            if null_rate(mid) <= alpha:
                high = mid
            else:
                low = mid + 1
        return low

    def above(boundary):
        return old + float(mass @ (miss * binom.sf(boundary - h, j - shared, gamma)))

    boundary = first_feasible(j, above)
    base = above(boundary)
    tie_mass = mass * binom.pmf(boundary - h, j - shared, gamma)
    previous_tail = binom.sf(first - 1 - h, k - shared, gamma)

    def with_tie(t):
        delta = binom.sf(t - 1 - h, k - shared, gamma) - previous_tail
        return base + float(tie_mass @ np.maximum(delta, 0.0))

    tie = first_feasible(first, with_tie)
    return first, boundary, tie, with_tie(tie)


def paired_audit(primary, secondary, key, *, gamma=0.625, alpha=0.001):
    """Test two prespecified views. Each view contains one string per usable structure.

    Repeated strings share one label, including across views. Supply the marking
    probability used by the producer. Both views must be fixed independently
    of the audit key under the null.
    """
    if not 0 < gamma < 1 or not 0 < alpha < 1:
        raise ValueError("gamma and alpha must lie strictly between zero and one")
    a, b = set(primary), set(secondary)
    labels = {fp: is_green(fp, key, gamma) for fp in a | b}
    g, h = sum(labels[fp] for fp in a), sum(labels[fp] for fp in b)
    k, j, shared = len(a), len(b), len(a & b)
    c, d, tie, rate = _region(k, j, shared, gamma, alpha)
    old = g >= c
    return PairedAudit(k, j, shared, g, h, c, d, tie, rate, old, old or (h, g) >= (d, tie))


RECOVERY_SYMPREC = 0.30


def eligible_recovery(fp):
    """Omit unresolved and triclinic assignments from the secondary view."""
    if fp is None:
        return None
    group, _, _ = fp.partition("|")
    if not group.startswith("SG") or not group[2:].isdigit() or not 1 <= int(group[2:]) <= 230:
        raise ValueError("expected a structural no-cell fingerprint with a space-group number")
    return fp if int(group[2:]) >= 3 else None


def recovery_class(structure):
    """Secondary class of a deposited structure using the frozen symmetry rule."""
    from .core import fingerprint_structure

    try:
        fp = fingerprint_structure(structure, symprec=RECOVERY_SYMPREC, use_cell_shape=False)
    except (ValueError, RuntimeError, TypeError):
        return None
    return eligible_recovery(fp)


def audit_deposited(structures, key, *, gamma=0.625, alpha=0.001, ranker=None):
    """Audit already-deposited pymatgen structures with the frozen pair of symmetry views."""
    from .core import fingerprint_structure

    primary, secondary, ranked = [], [], []
    for structure in structures:
        fp = fingerprint_structure(structure, use_cell_shape=False)
        primary.append(fp)
        alt = recovery_class(structure)
        if alt is not None:
            secondary.append(alt)
        if ranker is not None:
            ranked.append((fp, alt, ranker(structure)))
    if ranker is not None:
        return ranked_audit(ranked, key, gamma=gamma, alpha=alpha)
    return paired_audit(primary, secondary, key, gamma=gamma, alpha=alpha)


@lru_cache(maxsize=512)
def _rank_region(k, j, shared, top, shared_top, gamma, alpha):
    """Extend the fixed two-view rejection region by the high-energy subset of view two."""
    c, d, tie, old = _region(k, j, shared, gamma, alpha)
    shape = (top + 1, j + 1, k + 1)
    joint = np.zeros(shape)
    joint[0, 0, 0] = 1.0
    groups = [
        ((1, 1, 1), shared_top),
        ((1, 1, 0), top - shared_top),
        ((0, 1, 1), shared - shared_top),
        ((0, 1, 0), j - top - shared + shared_top),
        ((0, 0, 1), k - shared),
    ]
    for step, count in groups:
        dest = tuple(slice(s, None) for s in step)
        src = tuple(slice(None, -s or None) for s in step)
        for _ in range(count):
            next_mass = joint * (1 - gamma)
            next_mass[dest] += gamma * joint[src]
            joint = next_mass
    primary = np.arange(k + 1)
    secondary = np.arange(j + 1)[:, None]
    covered = (primary >= c) | (secondary > d) | ((secondary == d) & (primary >= tie))
    joint[:, covered] = 0.0
    # Fill the remaining size in decreasing (high-energy, secondary, primary) order.
    ordered = joint[::-1, ::-1, ::-1].ravel().cumsum()
    n = int(np.searchsorted(ordered, alpha - old, side="right"))
    if n == 0:
        return (top + 1, 0, 0), old
    position = np.unravel_index(n - 1, shape)
    return tuple(int(s - 1 - i) for s, i in zip(shape, position, strict=True)), old + ordered[n - 1]


@dataclass(frozen=True)
class RankedAudit:
    geometric: PairedAudit
    ranked_classes: int
    ranked_green: int
    null_rate: float
    fires: bool


def upper_energy_classes(observations):
    """Take the highest three quarters of distinct primary classes by deposited energy."""
    representatives = {}
    for primary, secondary, score in observations:
        if not np.isfinite(score):
            raise ValueError("deposited ranking energies must be finite")
        value = (score, secondary or "")
        if primary not in representatives or value < representatives[primary]:
            representatives[primary] = value
    ranked = sorted((e, fp, alt) for fp, (e, alt) in representatives.items())
    return {alt for _, _, alt in ranked[len(ranked) // 4 :] if alt}


def ranked_audit(observations, key, *, gamma=0.625, alpha=0.001):
    """Extend the calibrated geometry decision for a release of up to 32 distinct classes.

    Supply (primary, eligible recovery class or None, deposited energy) triples.
    Every view is key-independent under the null. Larger collections retain the
    two-view decision to bound the memory used by exact joint calibration.
    """
    observations = list(observations)
    primary = {p for p, _, _ in observations}
    secondary = {s for _, s, _ in observations if s is not None}
    geometric = paired_audit(primary, secondary, key, gamma=gamma, alpha=alpha)
    if max(len(primary), len(secondary)) > 32 or any(e is None for _, _, e in observations):
        return RankedAudit(geometric, 0, 0, geometric.null_rate, geometric.fires)
    top = upper_energy_classes(observations)
    h = sum(is_green(fp, key, gamma) for fp in top)
    threshold, size = _rank_region(
        len(primary),
        len(secondary),
        len(primary & secondary),
        len(top),
        len(primary & top),
        gamma,
        alpha,
    )
    return RankedAudit(
        geometric,
        len(top),
        h,
        float(size),
        geometric.fires or (h, geometric.secondary_green, geometric.primary_green) >= threshold,
    )
