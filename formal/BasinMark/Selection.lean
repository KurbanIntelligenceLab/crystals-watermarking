import BasinMark.Calibration
import Mathlib.MeasureTheory.Integral.Bochner.Set

namespace BasinMark
open MeasureTheory ProbabilityTheory Finset
open scoped BigOperators unitInterval

/-- The selected candidate is green exactly when the candidate batch contains
a green candidate. The two membership hypotheses encode the algorithm's fallback. -/
theorem selection_green_iff {ι : Type*} [DecidableEq ι] (batch : Finset ι)
    (green : ι → Prop) [DecidablePred green] (chosen : ι)
    (hchosen : chosen ∈ batch)
    (heligible : (batch.filter green).Nonempty → chosen ∈ batch.filter green) :
    green chosen ↔ (batch.filter green).Nonempty := by
  constructor
  · intro hg
    exact ⟨chosen, Finset.mem_filter.mpr ⟨hchosen, hg⟩⟩
  · intro h
    exact (Finset.mem_filter.mp (heligible h)).2

/-- Averaging the distinct-count formula over a random batch size. -/
theorem average_candidate_probability {Ω : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (D : Ω → ℕ) (γ : ℝ)
    (h : Integrable (fun ω ↦ (1 - γ) ^ D ω) μ) :
    (∫ ω, 1 - (1 - γ) ^ D ω ∂μ) = 1 - ∫ ω, (1 - γ) ^ D ω ∂μ := by
  rw [integral_sub (integrable_const 1) h]
  simp

/-- Finite-class occupancy expectation from the class-appearance indicators.
The absent-class probabilities must be established from the iid draw model. -/
theorem occupancy_expectation {Ω ι : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (s : Finset ι)
    (seen : ι → Set Ω) (hseen : ∀ i ∈ s, MeasurableSet (seen i))
    (v : ι → ℝ) (n : ℕ)
    (hprob : ∀ i ∈ s, μ.real (seen i) = 1 - (1 - v i) ^ n) :
    (∫ ω, ∑ i ∈ s, (seen i).indicator (fun _ ↦ (1 : ℝ)) ω ∂μ) =
      ∑ i ∈ s, (1 - (1 - v i) ^ n) := by
  rw [integral_finsetSum s (fun i hi ↦ (integrable_const (1 : ℝ)).indicator (hseen i hi))]
  apply Finset.sum_congr rfl
  intro i hi
  rw [integral_indicator (hseen i hi)]
  simpa using hprob i hi

/-- Finite occupancy formula from binomial marginal hit counts. No independence
between different classes is required for this expectation. -/
theorem occupancy_from_binomial_counts {Ω ι : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (s : Finset ι)
    (hits : ι → Ω → ℕ) (hhits : ∀ i ∈ s, Measurable (hits i))
    (v : ι → I) (n : ℕ)
    (hlaw : ∀ i ∈ s, HasLaw (hits i) (binomial n (v i)) μ) :
    (∫ ω, ∑ i ∈ s, {ω | hits i ω ≠ 0}.indicator (fun _ ↦ (1 : ℝ)) ω ∂μ) =
      ∑ i ∈ s, (1 - (1 - (v i : ℝ)) ^ n) := by
  apply occupancy_expectation μ s (fun i ↦ {ω | hits i ω ≠ 0})
    (fun i hi ↦ ((hhits i hi) (measurableSet_singleton 0)).compl) (fun i ↦ (v i : ℝ)) n
  intro i hi
  rw [(hlaw i hi).measureReal_eq (p := fun k ↦ k ≠ 0) (by measurability)]
  exact candidate_green_probability n (v i)

/-- Countable occupancy is the sum of individual appearance probabilities,
using nonnegative integration so the sum/expectation interchange is justified. -/
theorem countable_occupancy_expectation {Ω ι : Type*} [MeasurableSpace Ω] [Countable ι]
    (μ : Measure Ω) (seen : ι → Set Ω) (hseen : ∀ i, MeasurableSet (seen i)) :
    (∫⁻ ω, ∑' i, (seen i).indicator (fun _ ↦ (1 : ENNReal)) ω ∂μ) =
      ∑' i, μ (seen i) := by
  rw [lintegral_tsum (fun i ↦ (measurable_const.indicator (hseen i)).aemeasurable)]
  apply tsum_congr
  intro i
  rw [lintegral_indicator (hseen i)]
  simp

/-- Averaging conditional FPR bounds against a common probability law. -/
theorem average_null_bound {Ω : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (p : Ω → ENNReal) (α : ENNReal)
    (h : ∀ᵐ ω ∂μ, p ω ≤ α) : (∫⁻ ω, p ω ∂μ) ≤ α := by
  calc
    (∫⁻ ω, p ω ∂μ) ≤ ∫⁻ _, α ∂μ := lintegral_mono_ae h
    _ = α := by simp

end BasinMark
