"""Regression checks for release input validation and safe checkpoint reuse."""

import json
import sys
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import pytest

from basinmark import analyses
from basinmark.classification import ClassificationConfig
from basinmark.core import select_green


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1.0, 0.0])
def test_classification_rejects_invalid_thresholds(value):
    with pytest.raises(ValueError, match="positive"):
        ClassificationConfig("chgnet", Path("a"), Path("b"), Path("c"), Path("d"), fmax=value)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_selection_rejects_nonfinite_energies(value):
    with pytest.raises(ValueError, match="finite"):
        select_green(["a", "b"], [value, 0.0], b"test")


@pytest.fixture
def potential_run(tmp_path, monkeypatch):
    raw, dft = tmp_path / "raw.cif", tmp_path / "dft.cif"
    raw.write_text("raw input")
    dft.write_text("reference input")
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps([{"pair_id": "p0", "raw": str(raw), "dft": str(dft)}]))
    structure = Namespace(composition=Namespace(reduced_formula="Si"))
    monkeypatch.setitem(sys.modules, "pymatgen.core", SimpleNamespace(
        Structure=SimpleNamespace(from_file=lambda *_: structure)
    ))
    monkeypatch.setitem(sys.modules, "pymatgen.symmetry.analyzer", SimpleNamespace(
        SpacegroupAnalyzer=lambda *a, **k: SimpleNamespace(get_crystal_system=lambda: "cubic")
    ))
    monkeypatch.setattr(analyses, "relax", lambda *a, **k: (structure, 0.0, 1))
    monkeypatch.setattr(analyses, "fp", lambda *a, **k: "Si")
    monkeypatch.setattr(analyses, "rmsd_to_reference", lambda *a: 0.0)
    monkeypatch.setattr(analyses, "audit_size", lambda *a, **k: 10)
    return Namespace(
        pairs=str(pairs), out=str(tmp_path / "out.json"), engines=["chgnet"],
        gamma=0.5, alpha=0.001, m=8, symprec=0.1, cell_bin=0.05, fmax=0.02,
    ), raw


@pytest.mark.parametrize("change", ["settings", "structure"])
def test_completed_run_rejects_changed_inputs(potential_run, change):
    args, raw = potential_run
    analyses.potential_comparison(args)
    before = Path(args.out).read_bytes()
    if change == "settings":
        args.fmax = 0.1
    else:
        raw.write_text("modified structure, same filename")
    with pytest.raises(ValueError, match="inputs or settings differ"):
        analyses.potential_comparison(args)
    assert Path(args.out).read_bytes() == before


def test_completed_run_resumes_identical_inputs(potential_run):
    args, _ = potential_run
    analyses.potential_comparison(args)
    before = json.loads(Path(args.out).read_text())
    analyses.potential_comparison(args)
    assert json.loads(Path(args.out).read_text()) == before


@pytest.mark.parametrize("option,value", [("--symprec", "nan"), ("--fmax", "inf")])
def test_cli_rejects_invalid_tolerances_before_loading_inputs(tmp_path, capsys, option, value):
    with pytest.raises(SystemExit) as exc:
        analyses.main([
            "potential-comparison", "--pairs", "absent.json", "--out", str(tmp_path / "out.json"),
            option, value,
        ])
    assert exc.value.code == 2
    assert "finite and positive" in capsys.readouterr().err


def test_partial_run_rejects_changed_structure(potential_run):
    args, raw = potential_run
    analyses.potential_comparison(args)
    output = Path(args.out)
    completed = json.loads(output.read_text())["chgnet"]
    completed["done"] = 1
    output.write_text(json.dumps({"chgnet__partial": completed}))
    before = output.read_bytes()
    raw.write_text("changed structure")
    with pytest.raises(ValueError, match="inputs or settings differ"):
        analyses.potential_comparison(args)
    assert output.read_bytes() == before


def test_real_structure_fingerprint_is_translation_and_supercell_invariant():
    core = pytest.importorskip("pymatgen.core")
    from basinmark.core import fingerprint_structure

    structure = core.Structure(core.Lattice.cubic(5.43), ["Si", "Si"], [[0, 0, 0], [.25]*3])
    translated = structure.copy()
    translated.translate_sites(range(len(translated)), [.1, .2, .3])
    supercell = structure.copy()
    supercell.make_supercell([2, 2, 2])
    expected = fingerprint_structure(structure)
    assert fingerprint_structure(translated) == expected
    assert fingerprint_structure(supercell) == expected
