"""Budget-matched selection of interchangeable outputs from one candidate collection."""

from dataclasses import dataclass

import numpy as np

from .core import is_green


@dataclass
class Selection:
    indices: np.ndarray
    reference_energies: np.ndarray
    batch_sizes: np.ndarray | None = None


def _fixed(energy, green, batches, window):
    values = energy[batches]
    minimum = values.min(axis=1)
    eligible = green[batches] & (values <= minimum[:, None] + window)
    preferred = np.argmin(np.where(eligible, values, np.inf), axis=1)
    positions = np.where(eligible.any(axis=1), preferred, np.argmin(values, axis=1))
    return Selection(batches[np.arange(len(batches)), positions], minimum)


def _reserve(energy, green, batches, window):
    """Spend reserve candidates on unmarked slots in rotation, then on costly slots."""
    n, m = batches.shape
    start = max(1, m // 2)
    initial = batches[:, :start]
    best = initial[np.arange(n), np.argmin(energy[initial], axis=1)].copy()
    gi = np.argmin(np.where(green[initial], energy[initial], np.inf), axis=1)
    best_green = initial[np.arange(n), gi].copy()
    green_energy = np.where(green[initial].any(axis=1), energy[best_green], np.inf)
    sizes = np.full(n, start)
    cursor = 0
    for index in batches[:, start:].flat:
        fallback = green_energy > energy[best] + window
        unresolved = np.flatnonzero(fallback)
        if unresolved.size:
            later = unresolved[unresolved >= cursor]
            slot = int(later[0] if later.size else unresolved[0])
        else:
            slot = int(np.argmax(green_energy))
        cursor = (slot + 1) % n
        sizes[slot] += 1
        if energy[index] < energy[best[slot]]:
            best[slot] = index
        if green[index] and energy[index] < green_energy[slot]:
            best_green[slot] = index
            green_energy[slot] = energy[index]
    chosen = np.where(green_energy <= energy[best] + window, best_green, best)
    return Selection(chosen, energy[best], sizes)


def _pool(energy, green, order, n, window, fixed, *, mean_budget):
    """For each possible green count the cheapest set is its two color prefixes."""
    reference = energy[order[:n]]
    fixed_sorted = np.sort(energy[fixed.indices])
    limit = np.inf
    if mean_budget:
        limit = min(float(fixed_sorted.mean()), float(reference.mean()) + window)
        reference = np.sort(fixed.reference_energies)
        ceiling = reference + window
    else:
        ceiling = np.minimum(reference + window, fixed_sorted)
    greens = order[green[order]]
    reds = order[~green[order]]
    for count in range(min(n, len(greens)), max(0, n - len(reds)) - 1, -1):
        candidate = np.concatenate((greens[:count], reds[: n - count]))
        ranked = candidate[np.lexsort((candidate, energy[candidate]))]
        if np.all(energy[ranked] <= ceiling) and energy[ranked].mean() <= limit:
            return Selection(ranked, reference)
    raise AssertionError("the globally lowest-energy set must satisfy the two ceilings")


def _select(energy, green, m, window, policies):
    batches = np.arange(len(energy)).reshape(-1, m)
    order = np.argsort(energy, kind="stable")
    n = len(batches)
    fixed = _fixed(energy, green, batches, window)
    selected = {}
    for policy in policies:
        if policy == "fixed":
            value = fixed
        elif policy == "reserve":
            value = _reserve(energy, green, batches, window)
        elif policy in {"pool", "budget"}:
            value = _pool(energy, green, order, n, window, fixed, mean_budget=policy == "budget")
        elif policy == "unbounded":
            value = _fixed(energy, green, batches, np.inf)
        elif policy == "minimum":
            value = _fixed(energy, np.zeros_like(green), batches, 0)
        elif policy == "pooled_minimum":
            value = Selection(order[:n], energy[order[:n]])
        elif policy == "single":
            value = Selection(batches[:, 0], fixed.reference_energies)
        else:
            raise ValueError(f"unknown selection policy: {policy}")
        selected[policy] = value
    return selected


def _preview_selection(energy, green, batches, window, scores):
    """Search energy/preview trade-offs while retaining the budget selector's color counts."""
    order = np.argsort(energy, kind="stable")
    n = len(batches)
    fixed = _fixed(energy, green, batches, window)
    baseline = _pool(energy, green, order, n, window, fixed, mean_budget=True)
    limit = min(np.sort(energy[fixed.indices]).mean(), energy[order[:n]].mean() + window)
    ceiling = baseline.reference_energies + window
    accepted = baseline.indices
    best = scores[accepted].sum()
    count = int(green[accepted].sum())
    greens, reds = np.flatnonzero(green), np.flatnonzero(~green)
    for weight in (0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28):
        cost = energy - weight * scores
        combined = np.concatenate(
            (
                greens[np.argsort(cost[greens], kind="stable")[:count]],
                reds[np.argsort(cost[reds], kind="stable")[: n - count]],
            )
        )
        ranked = combined[np.lexsort((combined, energy[combined]))]
        value = scores[ranked].sum()
        if (
            value > best + 1e-6
            and np.all(energy[ranked] <= ceiling)
            and energy[ranked].mean() <= limit
        ):
            accepted, best = ranked, value
    return Selection(accepted, baseline.reference_energies)


def select_collection(
    fingerprints,
    energies,
    key,
    *,
    policy="budget",
    m=8,
    energy_window=0.1,
    gamma=None,
    previews=None,
):
    """Return distinct row indices and energy references for one exchangeable request.

    Consecutive blocks of m define the fixed-batch comparator. The pool selector
    returns selections in energy order. Its ith reference is the ith cheapest input.
    The budget rule references sorted fixed minima individually and caps mean
    energy against both fixed selected and pooled unmarked energies.
    Supply eV/atom for the default window. A collection with multiple different
    composition or property requests must be split by request before selection.
    """
    if gamma is None:
        gamma = 0.625 if policy in ("budget", "preview") else 0.5
    energy = np.asarray(energies, dtype=float)
    fps = list(fingerprints)
    if not isinstance(m, (int, np.integer)) or m < 1 or not len(fps) or len(fps) % m:
        raise ValueError("expected a positive whole number of complete candidate batches")
    if energy.shape != (len(fps),) or not np.isfinite(energy).all():
        raise ValueError("one finite energy is required per fingerprint")
    if not np.isfinite(energy_window) or energy_window < 0 or not 0 < gamma < 1:
        raise ValueError("invalid energy window or green probability")
    green = np.array([is_green(fp, key, gamma) for fp in fps])
    if policy == "preview":
        if previews is None:
            raise ValueError("the preview policy requires a measured preview pair per candidate")
        values = list(previews)
        if len(values) != len(fps):
            raise ValueError("expected one preview pair per candidate")
        scores = []
        for i, pair in enumerate(values):
            if pair is None:
                pair = (None, None)
            if (
                not isinstance(pair, (list, tuple))
                or len(pair) != 2
                or any(x is not None and not isinstance(x, str) for x in pair)
            ):
                raise ValueError("each preview must be a primary/secondary string pair or None")
            primary = pair[0] or fps[i]
            other = pair[1] or primary
            scores.append(
                0.35 * green[i]
                + 0.20 * is_green(primary, key, gamma)
                + 0.45 * is_green(other, key, gamma)
            )
        return _preview_selection(
            energy, green, np.arange(len(energy)).reshape(-1, m), energy_window, np.asarray(scores)
        )
    return _select(energy, green, m, energy_window, [policy])[policy]
