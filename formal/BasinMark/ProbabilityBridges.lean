import BasinMark.Selection
import Mathlib.Probability.Independence.InfinitePi

namespace BasinMark
open MeasureTheory ProbabilityTheory Set
open scoped BigOperators unitInterval

/-- Arbitrary mutually independent Bernoulli predicates induce the random-subset
law, rather than merely assuming that law for the full collection. -/
theorem independent_green_set_law {Ω ι : Type*} [MeasurableSpace Ω] [Fintype ι]
    (μ : Measure Ω) (B : ι → Ω → Prop) (p : I)
    (hB : ∀ i, HasLaw (B i) (bernoulliMeasure True False p) μ)
    (hind : iIndepFun B μ) :
    HasLaw (fun ω ↦ {i | B i ω}) (setBernoulli (univ : Set ι) p) μ := by
  have hpi := hind.hasLaw_infinitePi hB (aemeasurable_pi_iff.mpr (fun i ↦ (hB i).aemeasurable))
  have heq : setBernoulli (univ : Set ι) p =
      (Measure.infinitePi (fun _ : ι ↦ bernoulliMeasure True False p)).map
        (fun b : ι → Prop ↦ {i | b i}) := by
    simp [setBernoulli_eq_map, bernoulliMeasure_def]
  rw [heq]
  exact (show HasLaw (fun b : ι → Prop ↦ {i | b i})
    ((Measure.infinitePi (fun _ : ι ↦ bernoulliMeasure True False p)).map
      (fun b : ι → Prop ↦ {i | b i})) _ from ⟨by fun_prop, rfl⟩).fun_comp hpi

/-- The number of green classes has the binomial law for an arbitrary finite
mutually independent family, including an empty index type. -/
theorem independent_green_count_binomial {Ω ι : Type*} [MeasurableSpace Ω] [Fintype ι]
    (μ : Measure Ω) (B : ι → Ω → Prop) (p : I)
    (hB : ∀ i, HasLaw (B i) (bernoulliMeasure True False p) μ)
    (hind : iIndepFun B μ) :
    HasLaw (fun ω ↦ ({i | B i ω}).ncard) (binomial (Fintype.card ι) p) μ := by
  have hset := independent_green_set_law μ B p hB hind
  have heq : (setBernoulli (univ : Set ι) p).map Set.ncard = binomial (Fintype.card ι) p := by
    apply Measure.ext_of_singleton
    intro k
    simp [map_ncard_setBernoulli_singleton (finite_univ : (univ : Set ι).Finite), binomial_singleton]
  rw [← heq]
  exact (show HasLaw Set.ncard ((setBernoulli (univ : Set ι) p).map Set.ncard) _ from
    ⟨by fun_prop, rfl⟩).fun_comp hset

/-- Exact audit validity transferred to the observed independent class count. -/
theorem independent_green_audit_valid {Ω ι : Type*} [MeasurableSpace Ω] [Fintype ι]
    (μ : Measure Ω) (B : ι → Ω → Prop) (p : I)
    (hB : ∀ i, HasLaw (B i) (bernoulliMeasure True False p) μ)
    (hind : iIndepFun B μ) (α : ENNReal) :
    μ {ω | (binomial (Fintype.card ι) p) {z | ({i | B i ω}).ncard ≤ z} ≤ α} ≤ α := by
  have h := independent_green_count_binomial μ B p hB hind
  rw [h.measure_eq (p := fun g ↦ (binomial (Fintype.card ι) p) {z | g ≤ z} ≤ α)
    (by measurability)]
  exact binomial_audit_valid _ p α

/-- Candidate success probability from arbitrary independent class labels. -/
theorem independent_candidate_probability {Ω ι : Type*} [MeasurableSpace Ω] [Fintype ι]
    (μ : Measure Ω) (B : ι → Ω → Prop) (p : I)
    (hB : ∀ i, HasLaw (B i) (bernoulliMeasure True False p) μ)
    (hind : iIndepFun B μ) :
    μ.real {ω | ∃ i, B i ω} = 1 - (1 - (p : ℝ)) ^ Fintype.card ι := by
  have h := independent_green_count_binomial μ B p hB hind
  have hevent : {ω | ∃ i, B i ω} = {ω | ({i | B i ω}).ncard ≠ 0} := by
    ext ω
    exact (Set.ncard_pos (Set.toFinite _)).symm.trans Nat.pos_iff_ne_zero
  rw [hevent, h.measureReal_eq (p := fun k ↦ k ≠ 0) (by measurability)]
  exact candidate_green_probability _ p

/-- Appearance probability is derived directly from independent draws with a
common distribution. It is not supplied as an occupancy hypothesis. -/
theorem iid_class_appearance {Ω A : Type*} [MeasurableSpace Ω] [MeasurableSpace A]
    [MeasurableSingletonClass A] (μ : Measure Ω) [IsProbabilityMeasure μ]
    (ν : Measure A) [IsProbabilityMeasure ν] (n : ℕ) (X : Fin n → Ω → A)
    (hX : ∀ j, HasLaw (X j) ν μ) (hind : iIndepFun X μ) (a : A) :
    μ.real {ω | ∃ j, X j ω = a} = 1 - (1 - ν.real {a}) ^ n := by
  have hno := hind.meas_iInter (s := fun j ↦ {ω | X j ω ≠ a})
    (fun j ↦ ⟨{a}ᶜ, (measurableSet_singleton a).compl, rfl⟩)
  have hnoReal : μ.real (⋂ j, {ω | X j ω ≠ a}) = (1 - ν.real {a}) ^ n := by
    have h := congrArg ENNReal.toReal hno
    change μ.real (⋂ j, {ω | X j ω ≠ a}) = _ at h
    rw [ENNReal.toReal_prod] at h
    have hm : ∀ j, (μ {ω | X j ω ≠ a}).toReal = 1 - ν.real {a} := by
      intro j
      change μ.real {ω | X j ω ≠ a} = _
      rw [(hX j).measureReal_eq (p := fun z ↦ z ≠ a) (by measurability)]
      have hs : {z : A | z ≠ a} = ({a} : Set A)ᶜ := by ext z; simp [eq_comm]
      rw [hs]
      simpa using (measureReal_compl (μ := ν) (measurableSet_singleton a))
    simp_rw [hm] at h
    simpa using h
  have hevent : {ω | ∃ j, X j ω = a} = (⋂ j, {ω | X j ω ≠ a})ᶜ := by
    ext ω
    simp
  have hm : NullMeasurableSet (⋂ j, {ω | X j ω ≠ a}) μ :=
    NullMeasurableSet.iInter (fun j ↦ ((hX j).aemeasurable.nullMeasurable (measurableSet_singleton a)).compl)
  rw [hevent, measureReal_compl₀ hm, hnoReal]
  simp

/-- Countable occupancy expectation from an iid sample process. The class
appearance probabilities and the countable interchange are both proved. -/
theorem iid_countable_occupancy {Ω A : Type*} [MeasurableSpace Ω] [MeasurableSpace A]
    [MeasurableSingletonClass A] [Countable A]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (ν : Measure A) [IsProbabilityMeasure ν]
    (n : ℕ) (X : Fin n → Ω → A) (hm : ∀ j, Measurable (X j))
    (hX : ∀ j, HasLaw (X j) ν μ) (hind : iIndepFun X μ) :
    (∫⁻ ω, ∑' a, {ω | ∃ j, X j ω = a}.indicator (fun _ ↦ (1 : ENNReal)) ω ∂μ) =
      ∑' a, ENNReal.ofReal (1 - (1 - ν.real {a}) ^ n) := by
  have hseen : ∀ a, MeasurableSet {ω | ∃ j, X j ω = a} := by
    intro a
    simpa only [← Set.iUnion_ofPred, Set.preimage, Set.mem_singleton_iff, eq_comm] using
      MeasurableSet.iUnion (fun j ↦ (hm j) (measurableSet_singleton a))
  rw [countable_occupancy_expectation μ _ hseen]
  apply tsum_congr
  intro a
  rw [← iid_class_appearance μ ν n X hX hind a]
  exact (ENNReal.ofReal_toReal (measure_ne_top μ _)).symm

/-- The sum of appearance indicators is the number of distinct observed classes. -/
theorem occupancy_indicator_count {Ω A : Type*} (n : ℕ) (X : Fin n → Ω → A) (ω : Ω) :
    (∑' a, {ω | ∃ j, X j ω = a}.indicator (fun _ ↦ (1 : ENNReal)) ω) =
      ((Set.range (fun j ↦ X j ω)).ncard : ENNReal) := by
  classical
  have heq : (fun a ↦ {ω | ∃ j, X j ω = a}.indicator (fun _ ↦ (1 : ENNReal)) ω) =
      (Set.range (fun j ↦ X j ω)).indicator (fun _ ↦ (1 : ENNReal)) := by
    funext a
    simp only [Set.indicator_apply, Set.mem_ofPred_eq, Set.mem_range]
  rw [heq, ← tsum_subtype, ENNReal.tsum_set_one]
  exact_mod_cast (Set.finite_range (fun j ↦ X j ω)).cast_ncard_eq.symm

/-- The occupancy identity in terms of the actual distinct-class count. -/
theorem iid_distinct_count_expectation {Ω A : Type*} [MeasurableSpace Ω] [MeasurableSpace A]
    [MeasurableSingletonClass A] [Countable A]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (ν : Measure A) [IsProbabilityMeasure ν]
    (n : ℕ) (X : Fin n → Ω → A) (hm : ∀ j, Measurable (X j))
    (hX : ∀ j, HasLaw (X j) ν μ) (hind : iIndepFun X μ) :
    (∫⁻ ω, ((Set.range (fun j ↦ X j ω)).ncard : ENNReal) ∂μ) =
      ∑' a, ENNReal.ofReal (1 - (1 - ν.real {a}) ^ n) := by
  simpa only [occupancy_indicator_count] using iid_countable_occupancy μ ν n X hm hX hind

end BasinMark
