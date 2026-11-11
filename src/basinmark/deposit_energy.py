"""Source-pinned single-point ranking for deposited crystal structures."""

import hashlib
import math
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations

import pyarrow.parquet as pq
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.optimize import linprog

HULL = "a16020541166e42c3fb0b1c523311996a003c0b939d01b4437186e3355b95319"
WEIGHTS = "d4b14be9afa294eebdbe31a0280b26a0fa29715771e978cbfd3ca24e0d90307a"


def verify(path, digest):
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"SHA-256 mismatch for {path.name}")


def references(path):
    cols = ["species_at_sites", "energy", "mace_omat_energy"]
    data = pq.read_table(path, columns=cols).to_pydict()
    by_system = defaultdict(dict)
    for atoms, energy, omat in zip(*(data[col] for col in cols), strict=True):
        if not math.isclose(energy, omat, abs_tol=1e-7):
            raise ValueError("mixed potential energies in hull reference")
        counts = Counter(atoms)
        divisor = math.gcd(*counts.values())
        formula = tuple(sorted((el, n // divisor) for el, n in counts.items()))
        group = by_system[frozenset(counts)]
        group[formula] = min(group.get(formula, math.inf), energy / divisor)
    return by_system


def hull_at(composition, groups, solver="diagram"):
    entries = []
    elements = [el.symbol for el in composition.elements]
    for n in range(1, len(elements) + 1):
        for subset in combinations(elements, n):
            for formula, energy in groups.get(frozenset(subset), {}).items():
                entries.append(PDEntry(dict(formula), energy))
    if solver == "diagram":
        return float(PhaseDiagram(entries).get_hull_energy_per_atom(composition))
    # At one target composition the hull is the cheapest mixture of reference phases.
    # Solving that point avoids enumerating all facets of high-dimensional diagrams.
    fractions = [
        [entry.composition.get_atomic_fraction(el) for entry in entries] for el in elements
    ]
    target = [composition.get_atomic_fraction(el) for el in elements]
    result = linprog(
        [entry.energy_per_atom for entry in entries],
        A_eq=fractions,
        b_eq=target,
        bounds=(0, None),
        method="highs",
    )
    if not result.success:
        raise ValueError("reference hull does not span the target composition")
    return float(result.fun)


@lru_cache(maxsize=1)
def calculator(weights):
    from mace.calculators import MACECalculator

    return MACECalculator(model_paths=str(weights), device="cpu", default_dtype="float64")


def energy_per_atom(structure, weights):
    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.calc = calculator(weights)
    return float(atoms.get_potential_energy() / len(atoms))


class DepositedEnergy:
    """Score a geometry against a composition-matched public convex hull."""

    def __init__(self, weights, hull):
        verify(weights, WEIGHTS)
        verify(hull, HULL)
        self.weights = weights
        self.groups = references(hull)
        self.hulls = {}

    def __call__(self, structure):
        composition = structure.composition.fractional_composition
        if composition not in self.hulls:
            try:
                self.hulls[composition] = hull_at(composition, self.groups, "linear")
            except ValueError:
                self.hulls[composition] = None
        if self.hulls[composition] is None:
            return None
        return energy_per_atom(structure, self.weights) - self.hulls[composition]
