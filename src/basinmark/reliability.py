"""A calibrated reliability-ranked extension of the deposited-energy audit."""

import json
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from math import comb

import numpy as np

from .core import is_green
from .symmetry_audit import RankedAudit, _rank_region, _region, ranked_audit, upper_energy_classes


@lru_cache(maxsize=1)
def survival_counts():
    table = files("basinmark").joinpath("symmetry_survival.json").read_text(encoding="utf-8")
    return json.loads(table)


@lru_cache(maxsize=460)
def _survival(group, agree):
    counts = survival_counts()
    matched = [v for k, v in counts.items() if k.endswith(f"/{int(agree)}")]
    base = (sum(v[1] for v in matched) + 2) / (sum(map(sum, matched)) + 3)
    family = [counts.get(f"{group}/{i}", (0, 0)) for i in (0, 1)]
    by_group = (sum(v[1] for v in family) + 8 * base) / (sum(map(sum, family)) + 8)
    failures, successes = family[int(agree)]
    return (successes + 5 * by_group) / (failures + successes + 5)


def reliable_classes(observations):
    """Prioritize probable producer fingerprints using only the deposited collection."""
    representatives = {}
    for primary, secondary, energy in observations:
        value = (energy, secondary or "")
        if primary not in representatives or value < representatives[primary]:
            representatives[primary] = value
    ordered = sorted((e, p, s) for p, (e, s) in representatives.items())
    scores = {}
    for rank, (_, primary, secondary) in enumerate(ordered, start=1):
        if secondary:
            survival = _survival(secondary.partition("|")[0], primary == secondary)
            priority = survival * min(1.0, 4 * rank / len(ordered))
            scores[secondary] = max(scores.get(secondary, 0.0), priority)
    return set(sorted(scores, key=lambda fp: (scores[fp], fp), reverse=True)[:18])


def membership(primary, secondary, ranked, reliable):
    groups = Counter(
        (fp in ranked, fp in secondary, fp in primary, fp in reliable) for fp in primary | secondary
    )
    return tuple(sorted(groups.items()))


def _allocations(sizes, remaining):
    """Distribute a fixed number of red labels among occupied membership groups."""
    if not sizes:
        if remaining == 0:
            yield ()
        return
    lo = max(0, remaining - sum(sizes[1:]))
    for count in range(lo, min(remaining, sizes[0]) + 1):
        for tail in _allocations(sizes[1:], remaining - count):
            yield (count, *tail)


@lru_cache(maxsize=2048)
def reliability_region(groups, gamma, alpha):
    """Extend the preceding critical region using the additional green count first."""
    k = sum(n for (_, _, primary, _), n in groups if primary)
    j = sum(n for (_, secondary, _, _), n in groups if secondary)
    top = sum(n for (ranked, _, _, _), n in groups if ranked)
    shared = sum(n for (_, secondary, primary, _), n in groups if secondary and primary)
    shared_top = sum(n for (ranked, _, primary, _), n in groups if ranked and primary)
    first, second, tie, _ = _region(k, j, shared, gamma, alpha)
    third, size = _rank_region(k, j, shared, top, shared_top, gamma, alpha)
    inside = [(tuple(map(int, pattern[:3])), n) for pattern, n in groups if pattern[3]]
    outside = [(tuple(map(int, pattern[:3])), n) for pattern, n in groups if not pattern[3]]
    lengths = tuple(1 + sum(step[i] * n for step, n in outside) for i in range(3))
    rest = np.zeros(lengths)
    rest[0, 0, 0] = 1
    for step, n in outside:
        dest = tuple(slice(i, None) for i in step)
        src = tuple(slice(None, -i or None) for i in step)
        for _ in range(n):
            nxt = rest * (1 - gamma)
            nxt[dest] += gamma * rest[src]
            rest = nxt
    c, b, a = np.ogrid[: top + 1, : j + 1, : k + 1]
    covered = (
        (a >= first)
        | (b > second)
        | ((b == second) & (a >= tie))
        | (c > third[0])
        | ((c == third[0]) & ((b > third[1]) | ((b == third[1]) & (a >= third[2]))))
    )
    sizes = tuple(n for _, n in inside)
    total = sum(sizes)
    all_green = tuple(sum(step[i] * n for step, n in inside) for i in range(3))
    boundary = (total + 1, 0, 0, 0)
    # Integrate one slice for each reliable green count. Leave all previously covered
    # outcomes in place, including attainable outcomes interleaved with the new tail.
    for red in range(total + 1):
        layer = np.zeros(covered.shape)
        mass = gamma ** (total - red) * (1 - gamma) ** red
        for counts in _allocations(sizes, red):
            start = list(all_green)
            probability = mass
            for (step, n), d in zip(inside, counts, strict=True):
                probability *= comb(n, d)
                for i in range(3):
                    start[i] -= step[i] * d
            window = tuple(slice(i, i + length) for i, length in zip(start, lengths, strict=True))
            layer[window] += probability * rest
        layer[covered] = 0
        tail = layer[::-1, ::-1, ::-1].ravel().cumsum()
        upto = int(np.searchsorted(tail, alpha - size, side="right"))
        if upto:
            loc = np.unravel_index(upto - 1, layer.shape)
            point = tuple(int(n - 1 - i) for n, i in zip(layer.shape, loc, strict=True))
            boundary = (total - red, *point)
            size += float(tail[upto - 1])
        if upto < len(tail):
            break
    return boundary, float(size)


@dataclass(frozen=True)
class ReliableAudit(RankedAudit):
    ranked_fires: bool
    reliable_classes: int
    reliable_green: int


def reliable_audit(observations, key, *, gamma=0.625, alpha=0.001):
    """Audit deposited (primary, eligible recovery class, energy) triples."""
    observations = list(observations)
    ranked = ranked_audit(observations, key, gamma=gamma, alpha=alpha)
    geometric = ranked.geometric
    if max(geometric.primary_classes, geometric.secondary_classes) > 32 or any(
        e is None for _, _, e in observations
    ):
        return ReliableAudit(
            **vars(ranked), ranked_fires=ranked.fires, reliable_classes=0, reliable_green=0
        )
    primary = {p for p, _, _ in observations}
    secondary = {s for _, s, _ in observations if s is not None}
    reliable = reliable_classes(observations)
    groups = membership(primary, secondary, upper_energy_classes(observations), reliable)
    boundary, size = reliability_region(groups, gamma, alpha)
    green = sum(is_green(fp, key, gamma) for fp in reliable)
    counts = (green, ranked.ranked_green, geometric.secondary_green, geometric.primary_green)
    return ReliableAudit(
        geometric,
        ranked.ranked_classes,
        ranked.ranked_green,
        size,
        ranked.fires or counts >= boundary,
        ranked.fires,
        len(reliable),
        green,
    )
