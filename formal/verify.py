#!/usr/bin/env python3
"""Build and audit every project theorem. Run from any directory."""
from pathlib import Path
import argparse, re, shutil, subprocess, sys
ROOT = Path(__file__).resolve().parent
lake = shutil.which("lake") or str(Path.home() / ".elan/bin/lake")
def run(args):
    result = subprocess.run([lake, *args], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    if result.returncode:
        sys.exit(result.returncode)
    return result.stdout
args = argparse.ArgumentParser(description=__doc__)
args.add_argument("--clean", action="store_true", help="Rebuild project modules from source; keep dependency cache")
options = args.parse_args()
if options.clean and (ROOT / ".lake/build").exists():
    shutil.rmtree(ROOT / ".lake/build")
if not (ROOT / ".lake/packages/mathlib/.lake/build/lib/lean/Mathlib/Probability/Distributions/Binomial.olean").exists():
    run(["exe", "cache", "get"])
run(["build"])
text = run(["env", "lean", "AxiomAudit.lean"])
found = {}
for name, kind, axioms in re.findall(r"BASINMARK_AUDIT (\S+) (\S+) \[([^\]]*)\]", text):
    if name in found:
        sys.exit(f"Duplicate audit record: {name}")
    found[name] = (kind, set(filter(None, axioms.split(","))))
counts = re.findall(r"BASINMARK_AUDIT_COUNT (\d+)", text)
if len(counts) != 1 or int(counts[0]) != len(found) or not found:
    sys.exit("Missing or incomplete compiled-declaration audit")
# Independently check that every current source module and explicit declaration
# is present in the imported build. The compiled audit is the primary inventory.
imports = (ROOT / "BasinMark.lean").read_text()
for source in (ROOT / "BasinMark").rglob("*.lean"):
    module = ".".join(source.relative_to(ROOT).with_suffix("").parts)
    if not re.search(r"^import " + re.escape(module) + r"\s*$", imports, re.M):
        sys.exit(f"Source module missing from aggregate import: {module}")
    for name in re.findall(r"^\s*(?:theorem|lemma|def|axiom) (\w+)", source.read_text(), re.M):
        if "BasinMark." + name not in found:
            sys.exit(f"Source declaration missing from compiled audit: {source.name}:{name}")
allowed = {"propext", "Classical.choice", "Quot.sound"}
for name, (kind, axioms) in found.items():
    if kind == "axiom" or axioms - allowed:
        sys.exit(f"Unapproved proof dependencies for {name}: {axioms-allowed}")
theorems = sum(kind == "theorem" for kind, _ in found.values())
print(f"Verified {theorems} project theorems and {len(found)-theorems} other declarations; only standard Lean axioms occur.")
