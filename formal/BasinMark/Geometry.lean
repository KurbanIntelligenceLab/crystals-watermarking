import Mathlib.Analysis.Convex.Deriv
import Mathlib.Analysis.Convex.Jensen
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-! The finite-perimeter translation formula, spherical averaging and Euclidean
isoperimetry are NOT proved here. `geometric_bound_from_perimeter_inputs` states
that boundary explicitly. The convexity step is proved, not assumed. -/
namespace BasinMark
open Finset Set
open scoped BigOperators

theorem convex_negative_power (a : ℝ) (ha : a ≤ 0) :
    ConvexOn ℝ (Ioi 0) (fun t : ℝ ↦ t ^ a) := by
  apply convexOn_of_hasDerivWithinAt2_nonneg (f' := fun t ↦ a * t ^ (a - 1))
    (f'' := fun t ↦ a * ((a - 1) * t ^ (a - 1 - 1))) (convex_Ioi 0)
  · intro x hx
    exact (Real.hasDerivAt_rpow_const (Or.inl (ne_of_gt hx))).continuousAt.continuousWithinAt
  · intro x hx
    have hx0 : x ≠ 0 := ne_of_gt (by simpa using hx)
    exact (Real.hasDerivAt_rpow_const (Or.inl hx0)).hasDerivWithinAt
  · intro x hx
    have hx0 : x ≠ 0 := ne_of_gt (by simpa using hx)
    exact ((Real.hasDerivAt_rpow_const (p := a - 1) (Or.inl hx0)).const_mul a).hasDerivWithinAt
  · intro x hx
    have hx0 : 0 < x := by simpa using hx
    have hp : 0 ≤ x ^ (a - 1 - 1) := (Real.rpow_pos_of_pos hx0 _).le
    have hprod : 0 ≤ a * (a - 1) := mul_nonneg_of_nonpos_of_nonpos ha (by linarith)
    nlinarith [mul_nonneg hprod hp]

/-- Jensen's inequality with class-volume weights. -/
theorem collision_diversity_power_bound {ι : Type*} (s : Finset ι) (v : ι → ℝ)
    (hv : ∀ i ∈ s, 0 < v i) (hsum : ∑ i ∈ s, v i = 1)
    (D : ℝ) (hD : 0 < D) :
    (∑ i ∈ s, v i ^ 2) ^ (-(1 / D)) ≤ ∑ i ∈ s, (v i) ^ (1 - 1 / D) := by
  have hc := convex_negative_power (-(1 / D)) (neg_nonpos.mpr (le_of_lt (one_div_pos.mpr hD)))
  have hj := hc.map_sum_le (t := s) (w := v) (p := v)
    (fun i hi ↦ (hv i hi).le) hsum (fun i hi ↦ hv i hi)
  simp only [smul_eq_mul] at hj
  have heq : (∑ i ∈ s, v i * v i) = ∑ i ∈ s, v i ^ 2 := by
    apply Finset.sum_congr rfl; intro i hi; ring
  rw [heq] at hj
  calc
    _ ≤ ∑ i ∈ s, v i * v i ^ (-(1 / D)) := hj
    _ = _ := by
      apply Finset.sum_congr rfl
      intro i hi
      rw [show 1 - 1 / D = 1 + -(1 / D) by ring, Real.rpow_add (hv i hi), Real.rpow_one]

/-- Conditional assembly of the geometric lower bound. `limitFormula` and
`isoperimetry` are the analytic inputs still outside the formalization. -/
theorem geometric_bound_from_perimeter_inputs {ι : Type*} (s : Finset ι)
    (v perimeter : ι → ℝ) (hv : ∀ i ∈ s, 0 < v i)
    (hsum : ∑ i ∈ s, v i = 1) (D κ c slope : ℝ) (hD : 0 < D)
    (hκ : 0 ≤ κ) (hc : 0 ≤ c)
    (limitFormula : slope = κ / 2 * ∑ i ∈ s, perimeter i)
    (isoperimetry : ∀ i ∈ s, c * (v i) ^ (1 - 1 / D) ≤ perimeter i) :
    κ / 2 * c * (∑ i ∈ s, v i ^ 2) ^ (-(1 / D)) ≤ slope := by
  have hsumP : c * ∑ i ∈ s, (v i) ^ (1 - 1 / D) ≤ ∑ i ∈ s, perimeter i := by
    rw [Finset.mul_sum]
    exact Finset.sum_le_sum isoperimetry
  have hj := collision_diversity_power_bound s v hv hsum D hD
  have hcj := mul_le_mul_of_nonneg_left hj hc
  have htot := le_trans hcj hsumP
  have hfinal := mul_le_mul_of_nonneg_left htot (div_nonneg hκ (by norm_num : (0:ℝ) ≤ 2))
  rw [limitFormula]
  nlinarith [hfinal]

/-- The final local force/displacement implication. The Hessian-to-inner-product
step is an explicit input and is not certified by this theorem. -/
theorem displacement_from_gradient_bound {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] (d g : E) (μ : ℝ) (hμ : 0 < μ)
    (hstrong : μ * ‖d‖ ^ 2 ≤ inner (𝕜 := ℝ) d g) :
    ‖d‖ ≤ ‖g‖ / μ := by
  have hc := real_inner_le_norm d g
  apply (le_div_iff₀ hμ).2
  by_cases hd : ‖d‖ = 0
  · simp [hd]
  · have hdn : 0 < ‖d‖ := lt_of_le_of_ne (norm_nonneg d) (Ne.symm hd)
    nlinarith

/-- Conversion of a maximum per-atom force norm to a squared total-force bound. -/
theorem per_atom_force_bound {ι E : Type*} [NormedAddCommGroup E]
    (atoms : Finset ι) (force : ι → E) (f : ℝ)
    (hforce : ∀ i ∈ atoms, ‖force i‖ ≤ f) :
    (∑ i ∈ atoms, ‖force i‖ ^ 2) ≤ (atoms.card : ℝ) * f ^ 2 := by
  calc
    (∑ i ∈ atoms, ‖force i‖ ^ 2) ≤ ∑ _i ∈ atoms, f ^ 2 := by
      apply Finset.sum_le_sum
      intro i hi
      exact pow_le_pow_left₀ (norm_nonneg _) (hforce i hi) 2
    _ = _ := by simp

end BasinMark
