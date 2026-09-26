import BasinMark.Coarsening
import Mathlib.MeasureTheory.Measure.Map
import Mathlib.MeasureTheory.Measure.Prod

namespace BasinMark
open MeasureTheory Set
open scoped BigOperators

/-- Coarse singleton masses are derived from the actual pushforward measure. -/
theorem coarse_mass_fiber_sum {F C : Type*} [Countable F]
    [MeasurableSpace F] [MeasurableSingletonClass F]
    [MeasurableSpace C] [MeasurableSingletonClass C]
    (ν : Measure F) (merge : F → C) (hm : Measurable merge) (c : C) :
    (ν.map merge) {c} = ∑' i : merge ⁻¹' {c}, ν {i.val} := by
  rw [Measure.map_apply hm (measurableSet_singleton c)]
  simpa only [tsum_subtype, Set.indicator_univ, Set.preimage_id, Set.preimage_univ, measure_univ] using (tsum_measure_preimage_singleton (μ := ν) (f := id)
    (to_countable (merge ⁻¹' {c})) (fun i _ ↦ measurableSet_singleton i)).symm

/-- Countable coarsening inequality for the actual class-mass distributions. -/
theorem coarse_collision_mass {F C : Type*} [Countable F]
    [MeasurableSpace F] [MeasurableSingletonClass F]
    [MeasurableSpace C] [MeasurableSingletonClass C]
    (ν : Measure F) (merge : F → C) (hm : Measurable merge) :
    (∑' i, (ν {i}) ^ 2) ≤ ∑' c, ((ν.map merge) {c}) ^ 2 := by
  simp_rw [coarse_mass_fiber_sum ν merge hm]
  have h := grouped_collision_mass (fun (c : C) (i : merge ⁻¹' {c}) ↦ ν {i.val})
  rw [ENNReal.tsum_fiberwise (fun i ↦ (ν {i}) ^ 2) merge] at h
  exact h

/-- Inverse collision diversity decreases under a measurable deterministic merge.
Extended nonnegative reals also handle zero-mass and infinite-measure cases. -/
theorem coarse_collision_diversity {F C : Type*} [Countable F]
    [MeasurableSpace F] [MeasurableSingletonClass F]
    [MeasurableSpace C] [MeasurableSingletonClass C]
    (ν : Measure F) (merge : F → C) (hm : Measurable merge) :
    (∑' c, ((ν.map merge) {c}) ^ 2)⁻¹ ≤ (∑' i, (ν {i}) ^ 2)⁻¹ :=
  ENNReal.inv_le_inv.mpr (coarse_collision_mass ν merge hm)

/-- A countable probability distribution has strictly positive collision mass,
bounded by one. This justifies taking real reciprocals. -/
theorem collision_mass_bounds {F : Type*} [Countable F]
    [MeasurableSpace F] [MeasurableSingletonClass F]
    (ν : Measure F) [IsProbabilityMeasure ν] :
    0 < (∑' i, (ν {i}) ^ 2) ∧ (∑' i, (ν {i}) ^ 2) ≤ 1 := by
  have hs : (∑' i, ν {i}) = 1 := by
    have hh := tsum_measure_preimage_singleton (μ := ν) (f := id)
      (to_countable (univ : Set F)) (fun i _ ↦ measurableSet_singleton i)
    change (∑' b : (univ : Set F), ν {b.val}) = ν univ at hh
    rw [tsum_subtype (univ : Set F) (fun i ↦ ν {i})] at hh
    simpa using hh
  constructor
  · by_contra h
    have hz : (∑' i, (ν {i}) ^ 2) = 0 := le_antisymm (le_of_not_gt h) bot_le
    have heach : ∀ i, ν {i} = 0 := by
      intro i
      have hb := ENNReal.le_tsum i (f := fun j ↦ (ν {j}) ^ 2)
      rw [hz] at hb
      simpa using hb
    simp [heach] at hs
  · calc
      (∑' i, (ν {i}) ^ 2) ≤ ∑' i, ν {i} := by
        apply ENNReal.tsum_le_tsum
        intro i
        have hb : ν {i} ≤ 1 := by simpa using (measure_mono (μ := ν) (Set.subset_univ {i}))
        calc
          (ν {i}) ^ 2 = ν {i} * ν {i} := pow_two _
          _ ≤ ν {i} * 1 := by gcongr
          _ = ν {i} := mul_one _
      _ = 1 := hs

/-- The manuscript's real-valued inverse squared-mass diversity inequality. -/
theorem coarse_real_collision_diversity {F C : Type*} [Countable F] [Countable C]
    [MeasurableSpace F] [MeasurableSingletonClass F]
    [MeasurableSpace C] [MeasurableSingletonClass C]
    (ν : Measure F) [IsProbabilityMeasure ν] (merge : F → C) (hm : Measurable merge) :
    (∑' c, ((ν.map merge).real {c}) ^ 2)⁻¹ ≤ (∑' i, (ν.real {i}) ^ 2)⁻¹ := by
  have : IsProbabilityMeasure (ν.map merge) := Measure.isProbabilityMeasure_map hm.aemeasurable
  have hf := collision_mass_bounds ν
  have hc := collision_mass_bounds (ν.map merge)
  have hft : (∑' i, (ν {i}) ^ 2) ≠ ⊤ := ne_top_of_le_ne_top (by simp) hf.2
  have hct : (∑' c, ((ν.map merge) {c}) ^ 2) ≠ ⊤ := ne_top_of_le_ne_top (by simp) hc.2
  have hreal := (ENNReal.toReal_le_toReal hft hct).mpr (coarse_collision_mass ν merge hm)
  have hp : 0 < (∑' i, (ν {i}) ^ 2).toReal := ENNReal.toReal_pos (ne_of_gt hf.1) hft
  have hi := one_div_le_one_div_of_le hp hreal
  simp only [one_div] at hi
  simpa only [ENNReal.tsum_toReal_eq (fun i ↦ ENNReal.pow_ne_top (measure_ne_top ν {i})),
    ENNReal.tsum_toReal_eq (fun c ↦ ENNReal.pow_ne_top (measure_ne_top (ν.map merge) {c})),
    ENNReal.toReal_pow, measureReal_def] using hi

end BasinMark
