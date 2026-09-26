import Mathlib.Probability.Distributions.Binomial
import Mathlib.Probability.Moments.Variance
import Mathlib.Tactic

namespace BasinMark
open MeasureTheory ProbabilityTheory Finset
open scoped BigOperators ProbabilityTheory unitInterval

variable {Ω ι : Type*} [MeasurableSpace Ω] {μ : Measure Ω}

/-- Bernoulli mean, proved from the two-atom measure. -/
theorem bernoulli_mean (p : I) (B : Ω → ℝ)
    (hB : HasLaw B (bernoulliMeasure (1 : ℝ) 0 p) μ) :
    (∫ ω, B ω ∂μ) = (p : ℝ) := by
  rw [hB.integral_eq, integral_bernoulliMeasure]
  simp

theorem bernoulli_memLp [IsProbabilityMeasure μ] (p : I) (B : Ω → ℝ)
    (hB : HasLaw B (bernoulliMeasure (1 : ℝ) 0 p) μ) : MemLp B 2 μ := by
  apply memLp_of_bounded (a := 0) (b := 1) ?_ hB.aemeasurable.aestronglyMeasurable
  apply (hB.ae_iff (p := fun x ↦ x ∈ Set.Icc (0 : ℝ) 1) (by measurability)).2
  rw [bernoulliMeasure_def, ae_add_measure_iff]
  constructor <;> exact Measure.ae_smul_measure (by simp) _

/-- Bernoulli variance, without relying on a placeholder library result. -/
theorem bernoulli_variance (p : I) (B : Ω → ℝ)
    (hB : HasLaw B (bernoulliMeasure (1 : ℝ) 0 p) μ) :
    variance B μ = (p : ℝ) * (1 - p) := by
  rw [hB.variance_eq, variance_eq_integral aemeasurable_id]
  simp only [integral_bernoulliMeasure, id_eq, smul_eq_mul, mul_one, mul_zero, add_zero]
  ring

theorem weighted_mean (s : Finset ι) (w : ι → ℝ) (B : ι → Ω → ℝ) (γ : ℝ)
    (hB : ∀ i ∈ s, Integrable (B i) μ)
    (hmean : ∀ i ∈ s, (∫ ω, B i ω ∂μ) = γ)
    (hn : ∑ i ∈ s, w i = 1) :
    (∫ ω, ∑ i ∈ s, w i * B i ω ∂μ) = γ := by
  rw [integral_finsetSum s (fun i hi ↦ (hB i hi).const_mul (w i))]
  simp_rw [integral_const_mul]
  calc
    ∑ i ∈ s, w i * ∫ ω, B i ω ∂μ = ∑ i ∈ s, w i * γ := by
      apply Finset.sum_congr rfl
      intro i hi
      rw [hmean i hi]
    _ = γ := by rw [← Finset.sum_mul, hn, one_mul]

theorem weighted_variance (s : Finset ι) (w : ι → ℝ) (B : ι → Ω → ℝ) (γ : ℝ)
    (hB : ∀ i ∈ s, MemLp (B i) 2 μ)
    (hind : Set.Pairwise (↑s : Set ι) (fun i j ↦ IndepFun (B i) (B j) μ))
    (hvar : ∀ i ∈ s, variance (B i) μ = γ * (1 - γ)) :
    variance (fun ω ↦ ∑ i ∈ s, w i * B i ω) μ =
      γ * (1 - γ) * ∑ i ∈ s, w i ^ 2 := by
  have hscaled : Set.Pairwise (↑s : Set ι)
      (fun i j ↦ IndepFun (fun ω ↦ w i * B i ω) (fun ω ↦ w j * B j ω) μ) := by
    intro i hi j hj hij
    exact (hind hi hj hij).comp (by fun_prop) (by fun_prop)
  have hv := IndepFun.variance_sum (s := s)
    (fun i hi ↦ (hB i hi).const_mul (w i)) hscaled
  rw [Finset.sum_fn] at hv
  rw [hv]
  simp_rw [variance_const_mul]
  calc
    ∑ i ∈ s, w i ^ 2 * variance (B i) μ = ∑ i ∈ s, w i ^ 2 * (γ * (1 - γ)) := by
      apply Finset.sum_congr rfl
      intro i hi
      rw [hvar i hi]
    _ = _ := by rw [← Finset.sum_mul]; ring

/-- Collision calibration directly from independent Bernoulli class labels. -/
theorem collision_calibration [IsProbabilityMeasure μ]
    (s : Finset ι) (w : ι → ℝ) (B : ι → Ω → ℝ) (p : I)
    (hB : ∀ i ∈ s, HasLaw (B i) (bernoulliMeasure (1 : ℝ) 0 p) μ)
    (hind : Set.Pairwise (↑s : Set ι) (fun i j ↦ IndepFun (B i) (B j) μ))
    (hn : ∑ i ∈ s, w i = 1) :
    (∫ ω, ∑ i ∈ s, w i * B i ω ∂μ) = (p : ℝ) ∧
    variance (fun ω ↦ ∑ i ∈ s, w i * B i ω) μ =
      (p : ℝ) * (1 - p) * ∑ i ∈ s, w i ^ 2 := by
  have hmem : ∀ i ∈ s, MemLp (B i) 2 μ := fun i hi ↦ bernoulli_memLp p (B i) (hB i hi)
  constructor
  · exact weighted_mean s w B p (fun i hi ↦ (hmem i hi).integrable (by norm_num))
      (fun i hi ↦ bernoulli_mean p (B i) (hB i hi)) hn
  · exact weighted_variance s w B p hmem hind (fun i hi ↦ bernoulli_variance p (B i) (hB i hi))

/-- Deduplicated green classes are modeled directly as an independent Bernoulli
random subset of the K class indices. Its cardinal is binomial. -/
theorem deduplicated_binomial (K : ℕ) (p : I) :
    Measure.map Set.ncard (setBernoulli (Set.Iio K) p) = binomial K p := rfl

theorem binomial_no_green (K : ℕ) (p : I) :
    (binomial K p).real {0} = (1 - (p : ℝ)) ^ K := binomial_real_zero K p

theorem candidate_green_probability (d : ℕ) (p : I) :
    (binomial d p).real {k | k ≠ 0} = 1 - (1 - (p : ℝ)) ^ d := by
  have hset : {k : ℕ | k ≠ 0} = ({0} : Set ℕ)ᶜ := rfl
  rw [hset, measureReal_compl (measurableSet_singleton 0)]
  simp

/-- A discrete upper-tail p-value is superuniform. This proves the actual
p-value rejection event, including the case when rejection is impossible. -/
theorem upper_tail_valid (ν : Measure ℕ) [IsProbabilityMeasure ν] (α : ENNReal) :
    ν {g | ν {z | g ≤ z} ≤ α} ≤ α := by
  by_cases hex : ∃ g : ℕ, ν {z | g ≤ z} ≤ α
  · let c := Nat.find hex
    have hc : ν {z | c ≤ z} ≤ α := Nat.find_spec hex
    apply le_trans (measure_mono ?_) hc
    intro g hg
    exact Nat.find_min' hex hg
  · have hempty : {g : ℕ | ν {z | g ≤ z} ≤ α} = ∅ := by
      ext g
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false]
      exact fun hg ↦ hex ⟨g, hg⟩
    rw [hempty, measure_empty]
    exact bot_le

theorem binomial_audit_valid (K : ℕ) (p : I) (α : ENNReal) :
    (binomial K p) {g | (binomial K p) {z | g ≤ z} ≤ α} ≤ α :=
  upper_tail_valid (binomial K p) α

end BasinMark
