# Imported geometric inequalities

Adapted from [hojonathanho/isoperimetric](https://github.com/hojonathanho/isoperimetric) at commit `29768f8beeaf17295cdf3853d37da35d7e2b0a5f`, under Apache-2.0 (see LICENSE).

The four Lean modules contain the upstream Prékopa–Leindler, Brunn–Minkowski and outer-expansion isoperimetric proofs. Adaptations add the `BasinMark` namespace, change import paths, replace deprecated names, update generated case names and use the current measure-preserving integration API. The inequality statements retain their mathematical hypotheses and conclusions.

These sources are compiled against the same pinned mathlib as the rest of the project. Their declarations and transitive proof dependencies are included in `AxiomAudit.lean`; no upstream precompiled binaries are trusted or shipped.

The outer-expansion inequality is not automatically Euclidean finite-perimeter isoperimetry. `GeometricConnection.lean` states the expansion-limit and approximation conditions used to transfer it. The existence of those approximations for the manuscript’s full class of cells remains outside the formalization.
