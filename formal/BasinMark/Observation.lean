import Mathlib.MeasureTheory.Integral.Lebesgue.Basic
import Mathlib.MeasureTheory.Measure.Map
import Mathlib.Data.Finset.Card
import Mathlib.Data.Multiset.MapFold

/-! Observation-channel and statistic-resolution results. Physical basin preservation
is an explicit hypothesis, not a claim about a numerical relaxation implementation. -/
namespace BasinMark
open MeasureTheory

variable {X Y Z K F O : Type*}

def channel (relax : X → Y) (canonical : Y → Z) (quantize : Z → O) : X → O :=
  quantize ∘ canonical ∘ relax

theorem annihilation (basin : X → F) (relax : X → Y)
    (constantOnBasin : ∀ x y, basin x = basin y → relax x = relax y)
    (mark : K → X → X) (preservesBasin : ∀ k x, basin (mark k x) = basin x)
    (canonical : Y → Z) (quantize : Z → O) (k : K) (x : X) :
    channel relax canonical quantize (mark k x) = channel relax canonical quantize x := by
  simp only [channel, Function.comp_apply]
  rw [constantOnBasin _ _ (preservesBasin k x)]

theorem detector_eq_of_observation_eq (observe : X → O) (detect : O → Y)
    {x y : X} (h : observe x = observe y) : detect (observe x) = detect (observe y) :=
  congrArg detect h

/-- Covers randomized detection pointwise for each common random seed. -/
theorem randomized_detector_eq (observe : X → O) (detect : O → Z → Y)
    {x y : X} (h : observe x = observe y) (seed : Z) :
    detect (observe x) seed = detect (observe y) seed := by rw [h]

theorem observation_laws_eq [MeasurableSpace X] [MeasurableSpace O]
    (μ : Measure X) (f g : X → O) (h : ∀ x, f x = g x) :
    Measure.map f μ = Measure.map g μ := by
  have : f = g := funext h
  rw [this]

/-- A measurable decision cannot separate equal statistic laws. -/
theorem statistic_laws_eq [MeasurableSpace X] [MeasurableSpace F] [MeasurableSpace O]
    (μ ν : Measure X) (statistic : X → F) (audit : F → O)
    (hs : Measurable statistic) (ha : Measurable audit)
    (h : Measure.map statistic μ = Measure.map statistic ν) :
    Measure.map (audit ∘ statistic) μ = Measure.map (audit ∘ statistic) ν := by
  rw [← Measure.map_map ha hs, ← Measure.map_map ha hs, h]

theorem rejection_probability_eq [MeasurableSpace X] [MeasurableSpace F]
    (μ ν : Measure X) (statistic : X → F) (hs : Measurable statistic)
    (h : Measure.map statistic μ = Measure.map statistic ν)
    (reject : Set F) (hr : MeasurableSet reject) :
    μ (statistic ⁻¹' reject) = ν (statistic ⁻¹' reject) := by
  rw [← Measure.map_apply hs hr, ← Measure.map_apply hs hr, h]

theorem fingerprint_multiset_substitution (fingerprint : X → F)
    (audit : K → Multiset F → O) (k : K) (xs ys : Multiset X)
    (h : xs.map fingerprint = ys.map fingerprint) :
    audit k (xs.map fingerprint) = audit k (ys.map fingerprint) := by rw [h]

theorem distinct_set_substitution [DecidableEq F] (fingerprint : X → F)
    (audit : K → Finset F → O) (k : K) (xs ys : Finset X)
    (h : xs.image fingerprint = ys.image fingerprint) :
    audit k (xs.image fingerprint) = audit k (ys.image fingerprint) := by rw [h]

/-- If conditional rejection probabilities agree almost everywhere for the
common key law, their unconditional probabilities agree. -/
theorem random_key_rejection_eq [MeasurableSpace K] (ν : Measure K)
    (p₀ p₁ : K → ENNReal) (h : p₀ =ᵐ[ν] p₁) :
    (∫⁻ k, p₀ k ∂ν) = ∫⁻ k, p₁ k ∂ν := lintegral_congr_ae h

end BasinMark
