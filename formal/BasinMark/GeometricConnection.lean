import BasinMark.Geometry
import BasinMark.External.Isoperimetric
import Mathlib.Analysis.Asymptotics.Lemmas
import Mathlib.MeasureTheory.Measure.Real
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.MeasureTheory.Integral.DominatedConvergence

namespace BasinMark
open MeasureTheory Set Filter
open scoped BigOperators Topology symmDiff Pointwise

/-- Equal-volume sets have equal one-sided difference masses. -/
theorem equal_volume_exit_half_symmDiff {E : Type*} [MeasurableSpace E]
    (μ : Measure E) (A B : Set E) (hA : MeasurableSet A) (hB : MeasurableSet B)
    (hfinA : μ A ≠ ⊤) (hfinB : μ B ≠ ⊤) (hvol : μ.real A = μ.real B) :
    μ.real (A \ B) = μ.real (A ∆ B) / 2 := by
  have ha := measureReal_sdiff_add_inter (μ := μ) (s := A) hB hfinA
  have hb := measureReal_sdiff_add_inter (μ := μ) (s := B) hA hfinB
  rw [Set.inter_comm B A] at hb
  rw [measureReal_symmDiff_eq hA hB hfinA hfinB]
  linarith

/-- The factor one-half is proved for the actual Euclidean translation event. -/
theorem translated_cell_exit_half_symmDiff {d : ℕ}
    (A : Set (EuclideanSpace ℝ (Fin d))) (hA : MeasurableSet A)
    (hfin : volume A ≠ ⊤) (h : EuclideanSpace ℝ (Fin d)) :
    volume.real (A \ ((fun x ↦ x + h) ⁻¹' A)) =
      volume.real (A ∆ ((fun x ↦ x + h) ⁻¹' A)) / 2 := by
  have hv : volume ((fun x ↦ x + h) ⁻¹' A) = volume A :=
    measure_preimage_add_right volume h A
  apply equal_volume_exit_half_symmDiff volume A _ hA (hA.preimage (by fun_prop)) hfin
  · simpa [hv] using hfin
  · simp only [measureReal_def, hv]

/-- Dominated directional limits are transferred through the direction law.
The cellwise finite-perimeter translation result remains an explicit input. -/
theorem averaged_directional_limit {S : Type*} [MeasurableSpace S]
    (σ : Measure S) [IsProbabilityMeasure σ] (q : ℝ → S → ℝ) (variation : S → ℝ)
    (P : ℝ)
    (hm : ∀ᶠ ε in 𝓝[>] (0 : ℝ), AEStronglyMeasurable (q ε) σ)
    (hb : ∀ᶠ ε in 𝓝[>] (0 : ℝ), ∀ᵐ u ∂σ, ‖q ε u‖ ≤ P)
    (hl : ∀ᵐ u ∂σ, Tendsto (fun ε ↦ q ε u) (𝓝[>] (0 : ℝ)) (𝓝 (variation u))) :
    Tendsto (fun ε ↦ ∫ u, q ε u ∂σ) (𝓝[>] (0 : ℝ)) (𝓝 (∫ u, variation u ∂σ)) :=
  tendsto_integral_filter_of_dominated_convergence (fun _ ↦ P) hm hb (integrable_const P) hl

/-- Finite cellwise directional limits yield the partition's first-order slope.
The integral identity supplies the precise spherical normalisation. -/
theorem partition_slope_from_directional_limits {S ι : Type*} [MeasurableSpace S]
    (σ : Measure S) [IsProbabilityMeasure σ] (s : Finset ι)
    (q : ι → ℝ → S → ℝ) (variation : ι → S → ℝ) (P : ι → ℝ) (κ : ℝ)
    (hm : ∀ i ∈ s, ∀ᶠ ε in 𝓝[>] (0 : ℝ), AEStronglyMeasurable (q i ε) σ)
    (hb : ∀ i ∈ s, ∀ᶠ ε in 𝓝[>] (0 : ℝ), ∀ᵐ u ∂σ, ‖q i ε u‖ ≤ P i)
    (hl : ∀ i ∈ s, ∀ᵐ u ∂σ,
      Tendsto (fun ε ↦ q i ε u) (𝓝[>] (0 : ℝ)) (𝓝 (variation i u)))
    (haverage : ∀ i ∈ s, (∫ u, variation i u ∂σ) = κ * P i) :
    Tendsto (fun ε ↦ (1 / 2 : ℝ) * ∑ i ∈ s, ∫ u, q i ε u ∂σ)
      (𝓝[>] (0 : ℝ)) (𝓝 (κ / 2 * ∑ i ∈ s, P i)) := by
  have h := Filter.Tendsto.const_mul (1 / 2 : ℝ)
    (tendsto_finsetSum s (fun i hi ↦ averaged_directional_limit σ (q i) (variation i) (P i)
      (hm i hi) (hb i hi) (hl i hi)))
  have heq : (1 / 2 : ℝ) * ∑ i ∈ s, ∫ u, variation i u ∂σ = κ / 2 * ∑ i ∈ s, P i := by
    rw [Finset.sum_congr rfl haverage, ← Finset.mul_sum]
    ring
  simpa only [heq] using h

/-- The averaged translated-cell loss is connected to the normalised symmetric
-difference integrals used in the dominated-convergence theorem. -/
theorem translated_partition_loss_quotient {d : ℕ} {ι S : Type*} [MeasurableSpace S]
    (σ : Measure S) (s : Finset ι) (A : ι → Set (EuclideanSpace ℝ (Fin d)))
    (hm : ∀ i ∈ s, MeasurableSet (A i)) (hf : ∀ i ∈ s, volume (A i) ≠ ⊤)
    (direction : S → EuclideanSpace ℝ (Fin d)) (ε : ℝ) :
    (∑ i ∈ s, ∫ u, volume.real (A i \ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) ∂σ) / ε =
      (1 / 2 : ℝ) * ∑ i ∈ s, ∫ u,
        volume.real (A i ∆ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) / ε ∂σ := by
  rw [Finset.sum_div, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i hi
  rw [← integral_div, ← integral_const_mul]
  apply integral_congr_ae
  exact Eventually.of_forall (fun u ↦ by
    dsimp only
    rw [translated_cell_exit_half_symmDiff (A i) (hm i hi) (hf i hi)]
    ring)

/-- The imported expansion inequality passes to a supplied expansion limit.
This does not identify that limit with BV perimeter for arbitrary sets. -/
theorem perimeter_bound_from_expansion_limit {d : ℕ}
    (A : Set (EuclideanSpace ℝ (Fin (d + 1)))) (hA : MeasurableSet A)
    (hne : A.Nonempty) (hfin : volume A ≠ ⊤) (P : ENNReal)
    (hlimit : Tendsto
      (fun ε : ℝ ↦ (volume (A + Metric.ball 0 ε) - volume A) / ENNReal.ofReal ε)
      (𝓝[>] (0 : ℝ)) (𝓝 P)) :
    (d + 1 : ENNReal) * (volume A) ^ (1 - ((d : ℝ) + 1)⁻¹) *
      volume (Metric.ball (0 : EuclideanSpace ℝ (Fin (d + 1))) 1) ^ ((d : ℝ) + 1)⁻¹ ≤ P := by
  apply ge_of_tendsto hlimit
  filter_upwards [self_mem_nhdsWithin] with ε hε
  exact isoperimetric_inequality hε hne hA hfin

/-- The imported bound in the manuscript's real volume/perimeter notation. -/
theorem real_perimeter_bound_from_expansion_limit {d : ℕ}
    (A : Set (EuclideanSpace ℝ (Fin (d + 1)))) (hA : MeasurableSet A)
    (hne : A.Nonempty) (hfin : volume A ≠ ⊤) (P : ℝ) (hP : 0 ≤ P)
    (hlimit : Tendsto
      (fun ε : ℝ ↦ (volume (A + Metric.ball 0 ε) - volume A) / ENNReal.ofReal ε)
      (𝓝[>] (0 : ℝ)) (𝓝 (ENNReal.ofReal P))) :
    ((d : ℝ) + 1) * (volume.real A) ^ (1 - ((d : ℝ) + 1)⁻¹) *
      (volume.real (Metric.ball (0 : EuclideanSpace ℝ (Fin (d + 1))) 1)) ^ ((d : ℝ) + 1)⁻¹ ≤ P := by
  have h := perimeter_bound_from_expansion_limit A hA hne hfin (ENNReal.ofReal P) hlimit
  have hr := ENNReal.toReal_mono ENNReal.ofReal_ne_top h
  have hcast : ((d : ENNReal) + 1).toReal = (d : ℝ) + 1 := by
    rw [ENNReal.toReal_add (by simp) (by simp)]
    simp
  simpa [hcast, ENNReal.toReal_mul, ENNReal.toReal_rpow, ENNReal.toReal_ofReal hP, measureReal_def] using hr

/-- Isoperimetry transfers from approximating sets when BOTH their volumes and
expansion perimeters converge. Existence of such an approximation is not assumed
silently or proved by this transfer lemma. -/
theorem perimeter_bound_from_approximation {d : ℕ}
    (A : ℕ → Set (EuclideanSpace ℝ (Fin (d + 1))))
    (hm : ∀ n, MeasurableSet (A n)) (hne : ∀ n, (A n).Nonempty)
    (hfin : ∀ n, volume (A n) ≠ ⊤) (P : ℕ → ℝ) (hP : ∀ n, 0 ≤ P n)
    (hexp : ∀ n, Tendsto
      (fun ε : ℝ ↦ (volume (A n + Metric.ball 0 ε) - volume (A n)) / ENNReal.ofReal ε)
      (𝓝[>] (0 : ℝ)) (𝓝 (ENNReal.ofReal (P n))))
    (v perimeter : ℝ) (hv : 0 < v)
    (hvol : Tendsto (fun n ↦ volume.real (A n)) atTop (𝓝 v))
    (hperim : Tendsto P atTop (𝓝 perimeter)) :
    ((d : ℝ) + 1) * v ^ (1 - ((d : ℝ) + 1)⁻¹) *
      (volume.real (Metric.ball (0 : EuclideanSpace ℝ (Fin (d + 1))) 1)) ^ ((d : ℝ) + 1)⁻¹ ≤ perimeter := by
  have hc := ((Real.continuousAt_rpow_const (x := v) (q := 1 - ((d : ℝ) + 1)⁻¹)
    (Or.inl (ne_of_gt hv))).tendsto.comp hvol).const_mul ((d : ℝ) + 1)
  apply le_of_tendsto_of_tendsto (hc.mul_const _) hperim
  exact Eventually.of_forall (fun n ↦ real_perimeter_bound_from_expansion_limit
    (A n) (hm n) (hne n) (hfin n) (P n) (hP n) (hexp n))

/-- The imported isoperimetry, approximation transfer and volume-weighted Jensen
step are connected to the manuscript's stability-diversity coefficient. -/
theorem partition_tradeoff_from_approximations {d : ℕ} {ι : Type*}
    (s : Finset ι) (v P : ι → ℝ) (hv : ∀ i ∈ s, 0 < v i)
    (hvolume : ∑ i ∈ s, v i = 1) (κ : ℝ) (hκ : 0 ≤ κ)
    (A : ι → ℕ → Set (EuclideanSpace ℝ (Fin (d + 1))))
    (hm : ∀ i ∈ s, ∀ n, MeasurableSet (A i n))
    (hne : ∀ i ∈ s, ∀ n, (A i n).Nonempty)
    (hfin : ∀ i ∈ s, ∀ n, volume (A i n) ≠ ⊤)
    (Pn : ι → ℕ → ℝ) (hPn : ∀ i ∈ s, ∀ n, 0 ≤ Pn i n)
    (hexp : ∀ i ∈ s, ∀ n, Tendsto
      (fun ε : ℝ ↦ (volume (A i n + Metric.ball 0 ε) - volume (A i n)) / ENNReal.ofReal ε)
      (𝓝[>] (0 : ℝ)) (𝓝 (ENNReal.ofReal (Pn i n))))
    (hvol : ∀ i ∈ s, Tendsto (fun n ↦ volume.real (A i n)) atTop (𝓝 (v i)))
    (hperim : ∀ i ∈ s, Tendsto (Pn i) atTop (𝓝 (P i))) :
    κ / 2 * (((d : ℝ) + 1) *
      (volume.real (Metric.ball (0 : EuclideanSpace ℝ (Fin (d + 1))) 1)) ^ ((d : ℝ) + 1)⁻¹) *
      (∑ i ∈ s, v i ^ 2) ^ (-(1 / ((d : ℝ) + 1))) ≤ κ / 2 * ∑ i ∈ s, P i := by
  apply geometric_bound_from_perimeter_inputs s v P hv hvolume ((d : ℝ) + 1) κ _ _
    (by positivity) hκ (by positivity) rfl
  intro i hi
  have h := perimeter_bound_from_approximation (A i) (hm i hi) (hne i hi) (hfin i hi)
    (Pn i) (hPn i hi) (hexp i hi) (v i) (P i) (hv i hi) (hvol i hi) (hperim i hi)
  simpa only [one_div, mul_assoc, mul_left_comm, mul_comm] using h

/-- A one-sided slope limit gives the actual little-o remainder. -/
theorem first_order_remainder_of_slope (p : ℝ → ℝ) (slope : ℝ)
    (h : Tendsto (fun ε ↦ p ε / ε) (𝓝[>] (0 : ℝ)) (𝓝 slope)) :
    (fun ε ↦ p ε - slope * ε) =o[𝓝[>] (0 : ℝ)] (fun ε ↦ ε) := by
  apply (Asymptotics.isLittleO_iff_tendsto' ?_).mpr
  · have hzero := h.sub_const slope
    simp only [sub_self] at hzero
    apply hzero.congr'
    filter_upwards [self_mem_nhdsWithin] with ε hε
    change 0 < ε at hε
    field_simp
  · filter_upwards [self_mem_nhdsWithin] with ε hε
    change 0 < ε at hε
    exact fun hz ↦ False.elim (ne_of_gt hε hz)

/-- The actual translated-cell loss has the stated slope and first-order
remainder once the cellwise BV and spherical-average inputs are supplied. -/
theorem translated_partition_asymptotic {d : ℕ} {ι S : Type*} [MeasurableSpace S]
    (σ : Measure S) [IsProbabilityMeasure σ] (s : Finset ι)
    (A : ι → Set (EuclideanSpace ℝ (Fin d)))
    (hA : ∀ i ∈ s, MeasurableSet (A i)) (hfin : ∀ i ∈ s, volume (A i) ≠ ⊤)
    (direction : S → EuclideanSpace ℝ (Fin d)) (variation : ι → S → ℝ)
    (P : ι → ℝ) (κ : ℝ)
    (hm : ∀ i ∈ s, ∀ᶠ ε in 𝓝[>] (0 : ℝ), AEStronglyMeasurable
      (fun u ↦ volume.real (A i ∆ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) / ε) σ)
    (hb : ∀ i ∈ s, ∀ᶠ ε in 𝓝[>] (0 : ℝ), ∀ᵐ u ∂σ,
      ‖volume.real (A i ∆ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) / ε‖ ≤ P i)
    (hl : ∀ i ∈ s, ∀ᵐ u ∂σ, Tendsto
      (fun ε ↦ volume.real (A i ∆ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) / ε)
      (𝓝[>] (0 : ℝ)) (𝓝 (variation i u)))
    (haverage : ∀ i ∈ s, (∫ u, variation i u ∂σ) = κ * P i) :
    let p := fun ε ↦ ∑ i ∈ s, ∫ u,
      volume.real (A i \ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) ∂σ
    Tendsto (fun ε ↦ p ε / ε) (𝓝[>] (0 : ℝ)) (𝓝 (κ / 2 * ∑ i ∈ s, P i)) ∧
      (fun ε ↦ p ε - (κ / 2 * ∑ i ∈ s, P i) * ε) =o[𝓝[>] (0 : ℝ)] (fun ε ↦ ε) := by
  dsimp only
  have hs := partition_slope_from_directional_limits σ s
    (fun i ε u ↦ volume.real (A i ∆ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) / ε)
    variation P κ hm hb hl haverage
  have h : Tendsto
      (fun ε ↦ (∑ i ∈ s, ∫ u, volume.real
        (A i \ ((fun x ↦ x + ε • direction u) ⁻¹' A i)) ∂σ) / ε)
      (𝓝[>] (0 : ℝ)) (𝓝 (κ / 2 * ∑ i ∈ s, P i)) := by
    simpa only [translated_partition_loss_quotient σ s A hA hfin direction] using hs
  exact ⟨h, first_order_remainder_of_slope _ _ h⟩

end BasinMark
