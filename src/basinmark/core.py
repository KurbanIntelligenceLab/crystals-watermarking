"""Fingerprint construction, keyed selection, and statistical auditing."""

from __future__ import annotations

import hashlib
import hmac
import math
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import beta as _beta
from scipy.stats import binom as _binom


def binom_sf(k: int, n: int, p: float) -> float:
    """P[X > k] for X ~ Bin(n, p).  Exact tail; never the normal approximation."""
    if n < 0 or not 0.0 <= p <= 1.0:
        raise ValueError("n must be nonnegative and p must be in [0, 1]")
    return float(_binom.sf(k, n, p))


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(_beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(_beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def holm_bonferroni(pvals, alpha: float = 0.05):
    """Return a boolean array: True where the hypothesis is rejected."""
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    m = p.size
    out = np.zeros(m, dtype=bool)
    for rank, idx in enumerate(order):
        if p[idx] <= alpha / (m - rank):
            out[idx] = True
        else:
            break
    return out


def critical_value(n: int, gamma: float, alpha: float) -> int:
    """Smallest c with P_H0[G >= c] <= alpha."""
    c = 0
    while c <= n and binom_sf(c - 1, n, gamma) > alpha:
        c += 1
    return c


def q_green(gamma: float, m: int) -> float:
    """Probability an emitted structure is green: q = 1 - (1-gamma)^m."""
    return 1.0 - (1.0 - gamma) ** m


def q_flipped(gamma: float, m: int, p: float) -> float:
    """Audited green probability when the fingerprint flips with probability p."""
    return (1.0 - p) * q_green(gamma, m) + p * gamma


def power(n: int, gamma: float, q: float, alpha: float) -> float:
    return binom_sf(critical_value(n, gamma, alpha) - 1, n, q)


def audit_size(
    gamma: float = 0.5,
    m: int = 8,
    p: float = 0.0,
    target: float = 0.9,
    alpha: float = 1e-3,
    cap: int = 500_000,
) -> int | str | None:
    """Return the first sample size that reaches the target power."""
    q = q_flipped(gamma, m, p)
    if q <= gamma:
        return None
    for n in range(1, cap + 1):
        if power(n, gamma, q, alpha) >= target:
            return n
    return f">{cap}"


def n_forge(gamma: float = 0.5, alpha: float = 1e-3) -> int:
    """Return the smallest n such that gamma**n is at most alpha."""
    return max(1, math.ceil(math.log(alpha) / math.log(gamma) - 1e-12))


def n_eff_plugin(counts) -> float:
    c = np.asarray(counts, dtype=float)
    if c.size == 0 or np.any(c < 0) or c.sum() <= 0:
        return 0.0
    w = c / c.sum()
    return float(1.0 / np.sum(w**2))


def n_eff_debiased(counts) -> float:
    """Return the bias-corrected effective number of classes."""
    c = np.asarray(counts, dtype=float)
    if c.size == 0 or np.any(c < 0) or c.sum() <= 0:
        return 0.0
    n = float(c.sum())
    s2_hat = float(np.sum((c / n) ** 2))
    if n <= 1:
        return 1.0
    s2 = (s2_hat - 1.0 / n) / (1.0 - 1.0 / n)
    if s2 <= 0:
        return n
    return float(min(max(1.0 / s2, 1.0), n))


def design_effect(counts) -> float:
    """Return the plug-in class-frequency design effect."""
    c = np.asarray(counts, dtype=float)
    if c.size == 0 or np.any(c < 0) or c.sum() <= 0:
        return 0.0
    return float(c.sum() * np.sum((c / c.sum()) ** 2))


def prf_uniform(fingerprint: str, key: bytes) -> float:
    """Keyed PRF into [0,1).  HMAC-SHA256, first 8 bytes as a big-endian int."""
    if not isinstance(key, bytes):
        raise TypeError("key must be bytes")
    d = hmac.new(key, fingerprint.encode("utf-8"), hashlib.sha256).digest()
    return int.from_bytes(d[:8], "big") / 2**64


def is_green(fingerprint: str, key: bytes, gamma: float = 0.5) -> bool:
    return prf_uniform(fingerprint, key) < gamma


def select_green(fingerprints, energies, key: bytes, gamma: float = 0.5):
    """Emit the lowest-energy green candidate; the lowest-energy one if none is
    green.  Returns (index, was_green)."""
    fps = list(fingerprints)
    e = np.asarray(energies, dtype=float)
    if not fps or e.ndim != 1 or len(fps) != e.size:
        raise ValueError("fingerprints and energies must be nonempty and have equal length")
    if not np.isfinite(e).all():
        raise ValueError("energies must be finite")
    green = [i for i, f in enumerate(fps) if is_green(f, key, gamma)]
    if green:
        return int(min(green, key=lambda i: e[i])), True
    return int(np.argmin(e)), False


@dataclass
class AuditResult:
    n_structures: int
    n_distinct: int
    n_eff: float
    design_effect: float
    green_count: int
    n_units: int
    p_value: float
    fires: bool
    green_fraction: float
    deduplicated: bool
    ci95: tuple = field(default=(float("nan"), float("nan")))


def audit(
    fingerprints, key: bytes, gamma: float = 0.5, alpha: float = 1e-3, deduplicate: bool = True
) -> AuditResult:
    """The audit of Section 4.  Deduplicated by default, because the naive form
    over-rejects by the design effect and is only exposed for demonstration."""
    fps = list(fingerprints)
    counts = Counter(fps)
    units = list(counts) if deduplicate else fps
    g = sum(1 for f in units if is_green(f, key, gamma))
    n = len(units)
    pval = binom_sf(g - 1, n, gamma) if n else 1.0
    lo, hi = clopper_pearson(g, n) if n else (float("nan"), float("nan"))
    return AuditResult(
        n_structures=len(fps),
        n_distinct=len(counts),
        n_eff=n_eff_debiased(list(counts.values())),
        design_effect=design_effect(list(counts.values())),
        green_count=g,
        n_units=n,
        p_value=pval,
        fires=bool(n and g >= critical_value(n, gamma, alpha)),
        green_fraction=(g / n if n else float("nan")),
        deduplicated=deduplicate,
        ci95=(lo, hi),
    )


def fingerprint_structure(
    structure,
    symprec: float = 0.1,
    cell_bin: float = 0.05,
    use_cell_shape: bool = True,
    composition_only: bool = False,
) -> str:
    """Return the canonical structural fingerprint for a pymatgen structure."""
    if any(not math.isfinite(value) or value <= 0 for value in (symprec, cell_bin)):
        raise ValueError("symprec and cell_bin must be positive")
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "fingerprint_structure needs pymatgen; use the "
            "synthetic backend for statistics-only work"
        ) from exc

    reduced = structure.composition.reduced_formula
    if composition_only:
        return f"COMP|{reduced}"

    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    prim = sga.get_primitive_standard_structure()
    sga = SpacegroupAnalyzer(prim, symprec=symprec)
    sym = sga.get_symmetry_dataset()
    sg = int(sym.number)

    site_sym = sym.site_symmetry_symbols
    equiv = sym.equivalent_atoms
    species = [str(s.specie) for s in prim]
    orbit_size = Counter(equiv)
    triples = sorted({(species[i], site_sym[i], orbit_size[equiv[i]]) for i in range(len(prim))})
    gcd = 0
    for _, _, mult in triples:
        gcd = math.gcd(gcd, int(mult))
    gcd = gcd or 1
    orbits = ";".join(f"{el}:{ss}:{mult // gcd}" for el, ss, mult in triples)

    parts = [f"SG{sg}", reduced, orbits]
    if use_cell_shape:
        a, b, c = prim.lattice.abc
        al, be, ga = prim.lattice.angles
        ratios = [b / a, c / a, al, be, ga]  # scale-free
        parts.append("|".join(f"{round(r / cell_bin) * cell_bin:.4f}" for r in ratios))
    return "|".join(parts)


def synthetic_fingerprints(
    n: int, n_classes: int, concentration: float = 1.0, rng=None, tag: str = "S"
):
    """Draw n fingerprints over n_classes with Dirichlet(concentration) weights.
    Small concentration concentrates the mass and drives n_eff down."""
    rng = rng or np.random.default_rng(0)
    w = rng.dirichlet(np.ones(n_classes) * concentration)
    idx = rng.choice(n_classes, size=n, p=w)
    return [f"{tag}-{i}" for i in idx]


def synthetic_energies(n: int, alpha: float = 0.89, rng=None):
    """Energies with F(t) ~ c t^alpha near the lower endpoint, matching the tail
    model of the cost proposition."""
    rng = rng or np.random.default_rng(0)
    return rng.random(n) ** (1.0 / alpha)


def self_test() -> int:
    rng = np.random.default_rng(11)
    ok = []

    def chk(name, cond):
        ok.append((name, bool(cond)))

    chk("q_green(1/2, 8) == 255/256", abs(q_green(0.5, 8) - 255 / 256) < 1e-12)
    chk("q_flipped(p=1) == gamma", abs(q_flipped(0.5, 8, 1.0) - 0.5) < 1e-12)
    chk("audit_size(m=8) == 10", audit_size(m=8) == 10)
    chk("audit_size(m=2) == 71", audit_size(m=2) == 71)
    high_flip = audit_size(p=0.7)
    low_flip = audit_size(p=0.1)
    chk(
        "audit_size grows with flip rate",
        isinstance(high_flip, int) and isinstance(low_flip, int) and high_flip > low_flip,
    )
    chk("n_forge(1/2, 1e-3) == 10", n_forge(0.5, 1e-3) == 10)
    chk("n_forge == audit_size at large m", n_forge(0.5) == audit_size(m=32))
    chk("n_forge(0.1, 1e-2) == 2 (log-space boundary)", n_forge(0.1, 1e-2) == 2)

    chk("n_eff of k equal classes is k", abs(n_eff_plugin([9] * 7) - 7) < 1e-9)
    chk("debiased n_eff clamps at n", n_eff_debiased([1] * 70) <= 70 + 1e-9)
    chk(
        "deff == n / n_eff_plugin",
        abs(design_effect([5, 3, 2]) - 10 / n_eff_plugin([5, 3, 2])) < 1e-9,
    )

    key = b"k0"
    fps = [f"F{i}" for i in range(4000)]
    frac = sum(is_green(f, key, 0.5) for f in fps) / len(fps)
    chk("PRF green fraction ~ gamma", abs(frac - 0.5) < 0.03)
    chk("PRF is deterministic", is_green("F1", key) == is_green("F1", key))
    agree = sum(is_green(f, b"k0") == is_green(f, b"k1") for f in fps) / len(fps)
    chk("PRF depends on the key (labels agree ~half the time, not always)", abs(agree - 0.5) < 0.03)

    e = synthetic_energies(8, rng=rng)
    i, was = select_green([f"G{j}" for j in range(8)], e, key)
    chk("selection returns a valid index", 0 <= i < 8)
    chk("selection emits a green candidate when one exists", (not was) or is_green(f"G{i}", key))

    marked = [f"M{i}" for i in range(40)]
    marked = [f for f in marked if is_green(f, key)]
    chk(
        "all-green collection above the frontier fires",
        audit(marked, key).fires if len(marked) >= 10 else True,
    )
    unmarked = synthetic_fingerprints(400, 380, 5.0, rng, tag="U")
    chk("unmarked collection does not fire", not audit(unmarked, key).fires)

    low = synthetic_fingerprints(400, 12, 0.08, rng, tag="L")
    a_naive = audit(low, key, deduplicate=False)
    a_dedup = audit(low, key, deduplicate=True)
    chk("deduplication reduces the unit count", a_dedup.n_units < a_naive.n_units)
    chk("low diversity gives a large design effect", a_dedup.design_effect > 5)

    chk("Holm rejects nothing at the null", not holm_bonferroni([0.4, 0.6, 0.9], 0.05).any())
    chk("Holm rejects a clear signal", holm_bonferroni([1e-9, 0.6, 0.9], 0.05)[0])
    _, hi = clopper_pearson(0, 200)
    chk("Clopper-Pearson 0/200 upper bound ~ 0.018", 0.017 < hi < 0.019)

    bad = [n for n, v in ok if not v]
    for n, v in ok:
        print(f"  [{'PASS' if v else '**FAIL**'}] {n}")
    print(f"\n  {len(ok) - len(bad)}/{len(ok)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    import sys

    sys.exit(self_test())
