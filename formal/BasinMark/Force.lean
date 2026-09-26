import BasinMark.Geometry
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Analysis.Calculus.Gradient.Basic
import Mathlib.Analysis.InnerProductSpace.PiL2

namespace BasinMark
open Set

/-- A positive lower bound on the derivative of the gradient implies strong
monotonicity along a convex neighbourhood. The mean-value theorem supplies the
same inequality as integrating the Hessian in the supplementary proof. -/
theorem gradient_inner_bound_from_hessian {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] (U : Set E) (hU : Convex ℝ U)
    (g : E → E) (H : E → E →L[ℝ] E)
    (hderiv : ∀ z ∈ U, HasFDerivAt g (H z) z)
    (μ : ℝ) (hH : ∀ z ∈ U, ∀ d : E, μ * ‖d‖ ^ 2 ≤ inner (𝕜 := ℝ) d (H z d))
    (x y : E) (hx : x ∈ U) (hy : y ∈ U) :
    μ * ‖y - x‖ ^ 2 ≤ inner (𝕜 := ℝ) (y - x) (g y - g x) := by
  let d := y - x
  let f : ℝ → ℝ := fun t ↦ inner (𝕜 := ℝ) d (g (x + t • d))
  have hd : ∀ t ∈ Icc (0 : ℝ) 1,
      HasDerivAt f (inner (𝕜 := ℝ) d (H (x + t • d) d)) t := by
    intro t ht
    have hp : HasDerivAt (fun t : ℝ ↦ x + t • d) d t := by
      simpa using ((hasDerivAt_id t).smul_const d).const_add x
    have hg := (hderiv _ (hU.add_smul_sub_mem hx hy ht)).comp_hasDerivAt t hp
    simpa [f] using (hasDerivAt_const t d).inner ℝ hg
  have hcont : ContinuousOn f (Icc 0 1) := fun t ht ↦ (hd t ht).continuousAt.continuousWithinAt
  have hdiff : DifferentiableOn ℝ f (interior (Icc 0 1)) :=
    fun t ht ↦ (hd t (interior_subset ht)).differentiableAt.differentiableWithinAt
  have hb : ∀ t ∈ interior (Icc (0 : ℝ) 1), μ * ‖d‖ ^ 2 ≤ deriv f t := by
    intro t ht
    rw [(hd t (interior_subset ht)).deriv]
    exact hH _ (hU.add_smul_sub_mem hx hy (interior_subset ht)) d
  have hm := (convex_Icc (0 : ℝ) 1).mul_sub_le_image_sub_of_le_deriv
    hcont hdiff hb 0 (by simp) 1 (by simp) (by norm_num)
  simpa [f, d, inner_sub_right] using hm

/-- The local residual-force bound now follows from actual derivative/Hessian
hypotheses, with a stationary reference point. -/
theorem displacement_from_hessian {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] (U : Set E) (hU : Convex ℝ U)
    (g : E → E) (H : E → E →L[ℝ] E)
    (hderiv : ∀ z ∈ U, HasFDerivAt g (H z) z)
    (μ : ℝ) (hμ : 0 < μ)
    (hH : ∀ z ∈ U, ∀ d : E, μ * ‖d‖ ^ 2 ≤ inner (𝕜 := ℝ) d (H z d))
    (x y : E) (hx : x ∈ U) (hy : y ∈ U) (hstationary : g x = 0) :
    ‖y - x‖ ≤ ‖g y‖ / μ := by
  apply displacement_from_gradient_bound (y - x) (g y) μ hμ
  simpa [hstationary] using gradient_inner_bound_from_hessian U hU g H hderiv μ hH x y hx hy

/-- Energy-specialized statement: the vector field is the actual gradient of V. -/
theorem energy_displacement_from_hessian {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] [CompleteSpace E]
    (V : E → ℝ) (U : Set E) (hU : Convex ℝ U) (H : E → E →L[ℝ] E)
    (hderiv : ∀ z ∈ U, HasFDerivAt (gradient V) (H z) z)
    (μ : ℝ) (hμ : 0 < μ)
    (hH : ∀ z ∈ U, ∀ d : E, μ * ‖d‖ ^ 2 ≤ inner (𝕜 := ℝ) d (H z d))
    (x y : E) (hx : x ∈ U) (hy : y ∈ U) (hstationary : HasGradientAt V 0 x) :
    ‖y - x‖ ≤ ‖gradient V y‖ / μ :=
  displacement_from_hessian U hU (gradient V) H hderiv μ hμ hH x y hx hy hstationary.gradient

/-- Cartesian block norms give the stated square-root-of-atom-count factor. -/
theorem block_force_norm_bound {ι E : Type*} [Fintype ι] [NormedAddCommGroup E]
    (g : PiLp 2 (fun _ : ι ↦ E)) (f : ℝ) (hf : 0 ≤ f)
    (hforce : ∀ i, ‖g.ofLp i‖ ≤ f) :
    ‖g‖ ≤ Real.sqrt (Fintype.card ι) * f := by
  have hs : ‖g‖ ^ 2 ≤ (Fintype.card ι : ℝ) * f ^ 2 := by
    rw [PiLp.norm_sq_eq_of_L2]
    simpa using per_atom_force_bound Finset.univ (fun i ↦ g.ofLp i) f (fun i _ ↦ hforce i)
  have h := Real.sqrt_le_sqrt hs
  simpa [Real.sqrt_sq (norm_nonneg g), Real.sqrt_mul (Nat.cast_nonneg (Fintype.card ι)),
    Real.sqrt_sq hf] using h

end BasinMark
