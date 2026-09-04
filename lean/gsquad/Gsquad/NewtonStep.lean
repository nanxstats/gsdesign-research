import Mathlib.Tactic

/-!
Formal checks for the guarded Newton-Raphson update used in gsDesign's boundary
search (src/gsbound.c, src/gsbound1.c after the EXTREMEZ change).

The C code updates an iterate `b` towards the root of `phi(b) = target` with
`d = target - phi(b)` and derivative `dphi = phi'(b) ≤ 0` (upper bound: the
crossing probability decreases in the bound) as follows:

    if      d < dphi      then b + 1          (step capped at +1)
    else if d > -dphi     then b - 1          (step capped at -1)
    else if d = 0         then b              (exact hit; guards 0/0)
    else                  b + d / dphi        (Newton step)

and symmetrically for the lower bound with `dplo ≥ 0`. We prove, over the real
numbers, that the division branch is only reached when the derivative is
non-zero, so the update never divides by zero, and that every branch moves the
iterate by at most one unit; together with the iteration cap this bounds the
iterate after `n` steps by `|b₀| + n`.
-/

namespace GsDesign

/-- The guarded update for the upper bound (`dphi ≤ 0`). -/
noncomputable def upperStep (b d dphi : ℝ) : ℝ :=
  if d < dphi then b + 1
  else if d > -dphi then b - 1
  else if d = 0 then b
  else b + d / dphi

/-- The guarded update for the lower bound (`dplo ≥ 0`). -/
noncomputable def lowerStep (a d dplo : ℝ) : ℝ :=
  if d > dplo then a + 1
  else if d < -dplo then a - 1
  else if d = 0 then a
  else a + d / dplo

/-- In the Newton branch of the upper update the derivative is non-zero. -/
theorem upperStep_newton_branch_deriv_ne_zero (d dphi : ℝ)
    (h1 : ¬ d < dphi) (h2 : ¬ d > -dphi) (h3 : ¬ d = 0) : dphi ≠ 0 := by
  intro h0
  subst h0
  push Not at h1 h2
  simp at h1 h2
  exact h3 (le_antisymm h2 h1)

/-- In the Newton branch of the lower update the derivative is non-zero. -/
theorem lowerStep_newton_branch_deriv_ne_zero (d dplo : ℝ)
    (h1 : ¬ d > dplo) (h2 : ¬ d < -dplo) (h3 : ¬ d = 0) : dplo ≠ 0 := by
  intro h0
  subst h0
  push Not at h1 h2
  simp at h1 h2
  exact h3 (le_antisymm h1 h2)

/-- Every branch of the upper update moves the iterate by at most one unit,
provided the derivative has the sign the algorithm assumes (`dphi ≤ 0`). -/
theorem abs_upperStep_sub_le_one (b d dphi : ℝ) (hneg : dphi ≤ 0) :
    |upperStep b d dphi - b| ≤ 1 := by
  unfold upperStep
  split_ifs with h1 h2 h3
  · simp
  · simp
  · simp
  · -- Newton branch: dphi ≤ d ≤ -dphi with d ≠ 0 forces dphi < 0 and |d/dphi| ≤ 1
    push Not at h1 h2
    have hlt : dphi < 0 := by
      rcases lt_or_eq_of_le hneg with h | h
      · exact h
      · exfalso; subst h; simp at h1 h2; exact h3 (le_antisymm h2 h1)
    have hpos : 0 < -dphi := by linarith
    rw [add_sub_cancel_left, abs_div, abs_of_neg hlt, div_le_one hpos, abs_le]
    constructor <;> linarith

/-- Every branch of the lower update moves the iterate by at most one unit,
provided `dplo ≥ 0`. -/
theorem abs_lowerStep_sub_le_one (a d dplo : ℝ) (hpos : 0 ≤ dplo) :
    |lowerStep a d dplo - a| ≤ 1 := by
  unfold lowerStep
  split_ifs with h1 h2 h3
  · simp
  · simp
  · simp
  · push Not at h1 h2
    have hlt : 0 < dplo := by
      rcases lt_or_eq_of_le hpos with h | h
      · exact h
      · exfalso; subst h; simp at h1 h2; exact h3 (le_antisymm h1 h2)
    rw [add_sub_cancel_left, abs_div, abs_of_pos hlt, div_le_one hlt, abs_le]
    constructor <;> linarith

/-- Iterating any update that moves by at most one unit `n` times keeps the
iterate within `n` of its starting value (the iteration cap `GS_MAXITER = 20`
therefore bounds the excursion of the Newton search). -/
theorem abs_iterate_sub_le (f : ℝ → ℝ) (hf : ∀ x, |f x - x| ≤ 1) (b₀ : ℝ) :
    ∀ n : ℕ, |f^[n] b₀ - b₀| ≤ n := by
  intro n
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Function.iterate_succ_apply']
    calc |f (f^[n] b₀) - b₀|
        = |(f (f^[n] b₀) - f^[n] b₀) + (f^[n] b₀ - b₀)| := by ring_nf
      _ ≤ |f (f^[n] b₀) - f^[n] b₀| + |f^[n] b₀ - b₀| := abs_add_le _ _
      _ ≤ 1 + n := add_le_add (hf _) ih
      _ = (n + 1 : ℕ) := by push_cast; ring

end GsDesign
