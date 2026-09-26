import Mathlib.Algebra.Order.Chebyshev
import Mathlib.MeasureTheory.Measure.MeasureSpace
import Mathlib.Topology.Algebra.InfiniteSum.ENNReal
import Mathlib.Tactic

namespace BasinMark
open MeasureTheory Finset
open scoped BigOperators

variable {X F C : Type*}

theorem coarse_mismatch_subset (fine : X → F) (merge : F → C) :
    {p : X × X | merge (fine p.1) ≠ merge (fine p.2)} ⊆
      {p : X × X | fine p.1 ≠ fine p.2} := by
  intro p hp heq
  exact hp (congrArg merge heq)

theorem coarsening_mismatch_probability [MeasurableSpace (X × X)]
    (μ : Measure (X × X)) (fine : X → F) (merge : F → C) :
    μ {p | merge (fine p.1) ≠ merge (fine p.2)} ≤
      μ {p | fine p.1 ≠ fine p.2} :=
  measure_mono (coarse_mismatch_subset fine merge)

theorem coarsening_cardinality [DecidableEq F] [DecidableEq C]
    (S : Finset X) (fine : X → F) (merge : F → C) :
    (S.image (merge ∘ fine)).card ≤ (S.image fine).card := by
  rw [← Finset.image_image]
  exact Finset.card_image_le

/-- Collision events increase under merging, for any common pair law. In
particular this applies to two independent draws from the same population. -/
theorem coarsening_collision_probability [MeasurableSpace (X × X)]
    (μ : Measure (X × X)) (fine : X → F) (merge : F → C) :
    μ {p | fine p.1 = fine p.2} ≤ μ {p | merge (fine p.1) = merge (fine p.2)} := by
  apply measure_mono
  intro p h
  exact congrArg merge h

/-- Real collision diversity, once collision probabilities are identified. -/
theorem inverse_collision_antitone {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) :
    1 / b ≤ 1 / a := one_div_le_one_div_of_le ha hab

/-- Squared masses within one merged class include nonnegative cross terms. -/
theorem merged_mass_square {ι : Type*} (s : Finset ι) (v : ι → ℝ)
    (hv : ∀ i ∈ s, 0 ≤ v i) :
    (∑ i ∈ s, v i ^ 2) ≤ (∑ i ∈ s, v i) ^ 2 :=
  Finset.sum_sq_le_sq_sum_of_nonneg hv

/-- The variance-equivalent effective size is between one and class count. -/
theorem effective_size_bounds {ι : Type*} (s : Finset ι) (w : ι → ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i) (hn : ∑ i ∈ s, w i = 1) :
    1 ≤ 1 / (∑ i ∈ s, w i ^ 2) ∧
      1 / (∑ i ∈ s, w i ^ 2) ≤ (s.card : ℝ) := by
  have hupper : (∑ i ∈ s, w i ^ 2) ≤ 1 := by
    simpa [hn] using Finset.sum_sq_le_sq_sum_of_nonneg hw
  have hcs : 1 ≤ (s.card : ℝ) * ∑ i ∈ s, w i ^ 2 := by
    simpa [hn] using (sq_sum_le_card_mul_sum_sq (s := s) (f := w))
  have hpos : 0 < ∑ i ∈ s, w i ^ 2 := by
    have hnonneg : 0 ≤ ∑ i ∈ s, w i ^ 2 := Finset.sum_nonneg (by intros; positivity)
    by_contra h
    have hz : (∑ i ∈ s, w i ^ 2) = 0 := le_antisymm (le_of_not_gt h) hnonneg
    rw [hz] at hcs
    norm_num at hcs
  constructor
  · exact (le_div_iff₀ hpos).2 (by simpa using hupper)
  · exact (div_le_iff₀ hpos).2 (by nlinarith [hcs])

/-- Countable (indeed arbitrary nonnegative) grouped masses: merging cannot
reduce collision probability. Fibers are represented as a dependent family. -/
theorem grouped_collision_mass {C : Type*} {F : C → Type*}
    (v : (c : C) → F c → ENNReal) :
    (∑' c, ∑' i, (v c i) ^ 2) ≤ ∑' c, (∑' i, v c i) ^ 2 := by
  apply ENNReal.tsum_le_tsum
  intro c
  calc
    (∑' i, (v c i) ^ 2) ≤ ∑' i, v c i * (∑' j, v c j) := by
      apply ENNReal.tsum_le_tsum
      intro i
      rw [pow_two]
      gcongr
      exact ENNReal.le_tsum i
    _ = (∑' i, v c i) ^ 2 := by rw [ENNReal.tsum_mul_right, pow_two]

end BasinMark
