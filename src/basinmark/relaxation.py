"""Lazy potential construction and structure relaxation shared by analysis commands."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymatgen.core import Structure


class CalculatorRegistry:
    """Construct each requested potential once and reuse it within a process."""

    def __init__(self):
        self._calculators = {}

    def get(self, engine: str):
        key = engine.lower()
        if key not in self._calculators:
            self._calculators[key] = self._build(key)
        return self._calculators[key]

    @staticmethod
    def _build(engine: str):
        e = engine.lower()
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"
        if e == "chgnet":
            from chgnet.model.dynamics import CHGNetCalculator

            return CHGNetCalculator(use_device=device)
        elif e.startswith("mace"):
            from mace.calculators import mace_mp

            return mace_mp(model="medium-mpa-0", default_dtype="float64", device=device)
        elif e.startswith("orb"):
            from orb_models.forcefield import pretrained
            from orb_models.forcefield.calculator import ORBCalculator

            model = pretrained.orb_v3_conservative_inf_omat(device=device)
            return ORBCalculator(model, device=device)
        elif e.startswith("sevennet"):
            from sevenn.calculator import SevenNetCalculator

            return SevenNetCalculator(model="7net-mf-ompa", modal="mpa", device=device)
        elif e.startswith("mattersim"):
            from mattersim.forcefield import MatterSimCalculator

            return MatterSimCalculator(device=device)
        else:
            raise ValueError(f"unknown engine {engine!r}")


@dataclass(frozen=True)
class RelaxationResult:
    structure: Structure
    energy_per_atom: float
    steps: int


@dataclass
class RelaxationEngine:
    """Apply a potential with an explicit force threshold and step budget."""

    registry: CalculatorRegistry
    engine: str
    fmax: float = 0.02
    steps: int = 500

    def run(self, structure: Structure) -> RelaxationResult:
        from ase.filters import FrechetCellFilter
        from ase.optimize import FIRE
        from pymatgen.io.ase import AseAtomsAdaptor

        atoms = AseAtomsAdaptor.get_atoms(structure)
        atoms.calc = self.registry.get(self.engine)
        optimizer = FIRE(FrechetCellFilter(atoms), logfile=None)
        optimizer.run(fmax=self.fmax, steps=self.steps)
        return RelaxationResult(
            AseAtomsAdaptor.get_structure(atoms),
            float(atoms.get_potential_energy() / len(atoms)),
            int(optimizer.get_number_of_steps()),
        )
