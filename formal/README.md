# BasinMark formal proofs

This Lean 4 project checks mathematical components of the manuscript. It contains 67 project-authored theorem declarations and 32 adapted upstream theorem/lemma declarations across thirteen modules. The latest geometric connection adds eleven project theorems and builds the external inequalities from source. Full coverage of the finite-perimeter geometric theorem remains unfinished. The map below distinguishes verified conclusions from the analytic inputs still required. The project does not verify the Python implementation, numerical relaxation, empirical results or PRF security.

## Reproduce the verification

Install [elan](https://github.com/leanprover/elan), the Lean toolchain manager, and Python 3. Run from the companion repository:

```bash
python3 formal/verify.py --clean
```

The script resolves the project directory automatically, obtains missing mathlib build artifacts, rebuilds the project modules, and runs `AxiomAudit.lean`. It enumerates every compiled declaration in the `BasinMark` namespace, including the adapted external proofs, definitions and declarations introduced with `lemma`. It requires an axiom report for each one and rejects any dependency outside `propext`, `Classical.choice` and `Quot.sound`, including `sorryAx`. It also checks that each source module is imported and that explicit source declarations occur in the compiled inventory. `--clean` removes only this project's generated `.lake/build` directory and preserves dependency caches. Internet access is needed for the initial toolchain/dependency setup; subsequent checks use the local files.

The pinned environment is Lean `v4.33.1`, mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474`, and the transitive revisions in `lake-manifest.json`. Direct commands after setup are `lake build` and `lake env lean AxiomAudit.lean`. Include `lean-toolchain`, `lakefile.toml`, `lake-manifest.json`, the Lean sources, this coverage map and `verify.py` in the reproducibility artifact; exclude `.lake/`.

## Manuscript-to-Lean coverage

All declaration names below are in namespace `BasinMark`.

| Manuscript claim | Lean declarations | Verified scope |
|---|---|---|
| Exact-relaxation annihilation, `thm:annihilation` | `annihilation`, `observation_laws_eq`, `detector_eq_of_observation_eq`, `randomized_detector_eq` | A basin-label map, relaxation constant on its fibers, and a basin-preserving mark imply identical channel outputs. Pointwise equality yields equal pushforward laws, and any common detector or common random seed preserves equality. Admissible inputs can be represented by the domain type. No claim that a numerical relaxer satisfies these hypotheses. |
| Fingerprint-preserving substitution, `thm:duality` | `fingerprint_multiset_substitution`, `statistic_laws_eq`, `rejection_probability_eq`, `random_key_rejection_eq` | Multiset equality preserves a factored decision; equality of measurable statistic laws preserves output laws and rejection probabilities. The random-key lemma integrates an almost-everywhere equality of conditional rejection probabilities against the common key law. |
| Deduplicated set specialization | `distinct_set_substitution` | Equal finite fingerprint sets give the same set-based audit output, regardless of record multiplicity. |
| Collision moments, `prop:collision` | `bernoulli_mean`, `bernoulli_memLp`, `bernoulli_variance`, `weighted_mean`, `weighted_variance`, `collision_calibration` | The complete weighted mean/variance calculation follows from Bernoulli marginal laws and pairwise independence. Square integrability is derived from the Bernoulli law, not assumed in the final combined theorem. |
| Effective sample size | `effective_size_bounds` | Nonnegative weights summing to one give `1 ≤ 1 / sum(w_i²) ≤ K`. The record-to-class bound `K ≤ n` is the finite-image cardinality fact used by `coarsening_cardinality`; the equal-occupancy specialization and variance-equivalent rewrite are not separately packaged as declarations. |
| Binomial class count | `independent_green_set_law`, `independent_green_count_binomial`, `deduplicated_binomial`, `binomial_no_green` | Derives the random-subset law and its binomial cardinal from an arbitrary finite mutually independent Bernoulli predicate family. The joint law is proved from marginal laws and `iIndepFun`; it is no longer assumed as the starting model. Predicate truth is the green status. |
| Exact null test | `upper_tail_valid`, `binomial_audit_valid`, `independent_green_audit_valid`, `average_null_bound` | Proves upper-tail superuniformity and transfers it to the observed independent class count. This includes empty rejection regions and degenerate binomial parameters. The integration lemma averages supplied conditional bounds. |
| Candidate selection, `eq:candidate-green` | `selection_green_iff`, `independent_candidate_probability`, `candidate_green_probability`, `average_candidate_probability` | Derives the probability that at least one class is green from arbitrary independent Bernoulli predicates. The selection membership theorem characterizes green output, and the averaging lemma supplies the random-distinct-count algebra. Energy minimization and the construction of a conditional batch law are not verified. |
| Occupancy identity | `iid_class_appearance`, `iid_countable_occupancy`, `occupancy_indicator_count`, `iid_distinct_count_expectation` | Derives class-appearance probabilities from independent draws with a common distribution, identifies the sum of indicators with the actual distinct-class count, and proves the countable occupancy formula using nonnegative integration. This no longer assumes binomial hit-count laws. It includes zero samples and countably many possible labels. |
| Coarsening, `eq:coarsening` | Original event/cardinality lemmas; `coarse_mass_fiber_sum`, `coarse_collision_mass`, `collision_mass_bounds`, `coarse_collision_diversity`, `coarse_real_collision_diversity` | Derives grouped masses from the actual measurable pushforward distribution, proves squared-mass growth for countable labels, proves collision mass lies in `(0,1]` for probability measures, and obtains the manuscript’s real reciprocal inequality. No fiber-mass identification is assumed. |
| Local stability–diversity bound, `thm:tradeoff` | `GeometricConnection.lean`, `Geometry.lean`, and the four `External/` modules | Proves the translated-cell half-symmetric-difference identity, dominated averaging, finite-cell slope, and little-o remainder. Compiles external Prékopa–Leindler, Brunn–Minkowski and expansion isoperimetry; transfers the latter through supplied volume/perimeter approximations into the manuscript’s Jensen bound. Existence of appropriate finite-perimeter approximations, the cellwise directional-variation characterization and the spherical constant are still unproved inputs. |
| Finite-force bound | `gradient_inner_bound_from_hessian`, `displacement_from_hessian`, `energy_displacement_from_hessian`, `block_force_norm_bound` | Derives the gradient inner-product inequality from an actual Fréchet derivative/Hessian lower bound on a convex neighbourhood using the mean-value theorem, then proves displacement control at a stationary point. A separate block-norm theorem proves `norm(g) ≤ sqrt(N) f`. |

## Statement correspondence checks

The formalization uses abstract types for the manuscript's measurable structure, fingerprint and observation spaces. The following substitutions connect the Lean parameters to the proof notation; they do not assert that the physical implementation satisfies the hypotheses.

| Lean parameters | Manuscript interpretation and remaining connection |
|---|---|
| `basin`, `relax`, `canonical`, `quantize`, `mark` | Basin index, exact relaxation, canonicalisation, deposition precision map and keyed mark. The same input law and detector randomness must be used when transferring pointwise equality to hypothesis laws. |
| `s`, `w`, `B`, `p` in calibration | Distinct class indices, multiplicities divided by nonzero record count, class green indicators and `γ`. Pairwise independence suffices for moments; `ProbabilityBridges.lean` separately establishes the binomial law from mutual independence of the green predicates. |
| `ν`, `α` in `upper_tail_valid` | The count null law and rejection threshold. `independent_green_audit_valid` now proves and uses that law for the observed count. Empty index types and endpoint green probabilities are included. |
| `v`, `D`, `κ`, `c`, `slope` in the geometry assembly | Cell volumes, dimension, `κ_D`, `D ω_D^(1/D)` and the limit of `p(ε)/ε`. The Lean conclusion has `(sum v_i²)^(-1/D)`; identifying this with `D₂^(1/D)` remains a mathematical substitution. `GeometricConnection.lean` now instantiates the dimension and unit-ball coefficient through the imported expansion theorem; the spherical `κ_D` identity remains outside the formalization. |
| `d`, `g`, `μ` in the force implication | `y - x*`, the gradient at `y`, and the positive curvature lower bound. `Force.lean` establishes the formerly assumed inner-product inequality and the Cartesian block-norm conversion. The energy-specialized theorem uses `HasFDerivAt (gradient V) (H z) z` to identify `H` as the actual Hessian. It does not separately translate the manuscript’s `C²` notation into this derivative witness or combine the two final norm bounds into one declaration. |
| `F : C → Type*` in `grouped_collision_mass` | Fine-label fibers of a coarse-label map. `CoarseningMeasures.lean` now constructs their masses from `Measure.map`, justifies positive finite collision mass and proves the real-valued `D₂` comparison. |

The random-key and random-batch averaging declarations prove integration steps. They do not construct conditional distributions or derive independence from a keyed implementation. The iid occupancy theorem now starts from independent draws and proves the actual count expectation. No claim in the manuscript has been weakened to fit a Lean statement.

## Geometric connection and remaining analytic work

The following are actual connections between proved components; they do not establish that every finite-perimeter cell satisfies the remaining hypotheses.

| Connection | Verified result | Remaining premise |
|---|---|---|
| Translated-cell loss | `translated_cell_exit_half_symmDiff` and `translated_partition_loss_quotient` prove the factor `1/2` and the quotient formula for Euclidean translated sets. | Finite measurable cells. Interpreting the sum as a probability uses the manuscript’s disjoint unit-volume partition. |
| Directional limits to first order | `averaged_directional_limit`, `partition_slope_from_directional_limits`, `translated_partition_asymptotic`, and `first_order_remainder_of_slope` prove dominated averaging, the finite sum and the little-o remainder. | The directional quotient limit, domination, measurability and its spherical integral must still be established from finite perimeter. |
| External isoperimetry to a perimeter limit | `perimeter_bound_from_expansion_limit` and its real-valued version pass the imported inequality to an actual outer-expansion limit. | That expansion limit must be supplied. It is not identified with arbitrary BV perimeter. |
| Approximation to the manuscript coefficient | `perimeter_bound_from_approximation` and `partition_tradeoff_from_approximations` prove the transfer when approximating volumes and expansion perimeters converge. The final expression contains the Euclidean unit-ball volume and the required dimension exponent. | Existence of approximating sets with both convergences to the manuscript’s cell volume and BV perimeter. |

The finite-perimeter directional translation characterization, strict perimeter approximation and spherical average with the gamma-function constant remain the essential gaps. No finite-perimeter definition has been replaced by a convenient expansion-limit assumption. The normalised volume-weighted Jensen step is already proved. The code’s `κ` still needs to be identified with the manuscript’s spherical constant, and the inverse-square-mass notation remains a mathematical substitution for `D₂`.

The obstruction to a direct identification is mathematical: BV perimeter is invariant under null changes, whereas outer Minkowski content can react to null lower-dimensional pieces. [Galerne, pages 39–40 and 48](https://www.ias-iss.org/ojs/IAS/article/download/22/10/45) explicitly distinguishes these notions. Therefore an unqualified equality between outer-expansion perimeter and the manuscript perimeter would be invalid. A justified approximation route is required.

## External source provenance

The [planar formalization by Samarakkody](https://arxiv.org/abs/2603.14663) uses simple closed `C¹` curves and was not imported.

The four modules under `BasinMark/External/` adapt [hojonathanho/isoperimetric](https://github.com/hojonathanho/isoperimetric) at commit `29768f8beeaf17295cdf3853d37da35d7e2b0a5f`. Their Apache-2.0 license and adaptation notes are included. The upstream theorem statements were compared after whitespace normalization; the updates change namespace/import paths and proof compatibility with the pinned mathlib, while preserving the mathematical statements. The entire adapted dependency chain is rebuilt and axiom-audited with the project. This supersedes the earlier source-inspection-only status.

## Verification result

A clean rebuild and a second verification run passed without Lean warnings. The audit checked 164 compiled declarations, including the imported proofs and generated declarations, and found only `propext`, `Classical.choice` and `Quot.sound`. There are 67 project-authored and 32 adapted upstream theorem/lemma declarations; the compiled count also includes generated declarations and definitions. Kernel acceptance certifies those propositions. The unproved geometric premises above prevent a claim that the manuscript’s full geometric theorem is machine-checked.
