import HodgeComputer.Verify

/-!
# Kronecker-substitution checker for polynomial identities in `MvPolynomial (Fin 6) ℂ`

Explicit polynomials with integer coefficients are given as balanced binary trees of monomials (`PT`).
An identity `L = R` between expressions built from such polynomials (products, sums, natural scalars,
powers and partial derivatives of the explicit polynomials) is verified by the kernel in three steps:

* a homogeneity/degree check (`homE`, `degE`) that bounds every exponent of `L - R` by `n < m`;
* an ℓ¹ bound (`l1E`) that bounds every coefficient of `L - R` by `2^W`;
* the equality of the two Kronecker values `psiE m W L = psiE m W R`, i.e. the evaluations at
  `X i ↦ 2^(W·mⁱ)` for `i < 5` and `X 5 ↦ 1`.

Together these imply `L = R` (`sound`), first over `ℤ` and then over `ℂ` (`soundC`).  The Kronecker
values are computed with kernel-accelerated big-integer arithmetic (`Nat.shiftLeft`, `Nat.mul`), with a
divide-and-conquer evaluation (`psiRel`) that keeps every intermediate number proportional to the span of
monomial positions in its subtree, so memory stays `O(size · depth)`.

The only trusted interface is `MvPolynomial`; no new axioms, no `native_decide`.
-/

set_option autoImplicit false
set_option linter.unusedSimpArgs false

namespace HodgeComputer.Kron

open MvPolynomial

abbrev ZP := MvPolynomial (Fin 6) ℤ

/-- Exponent vector of a monomial in `X0 X1 X2 Y0 Y1 Y2`. -/
structure Mono where
  e0 : ℕ
  e1 : ℕ
  e2 : ℕ
  e3 : ℕ
  e4 : ℕ
  e5 : ℕ
deriving DecidableEq, Repr

/-- Balanced tree of monomials with integer coefficients (leaves sorted by `Mono.pos`). -/
inductive PT where
  | leaf (m : Mono) (c : ℤ)
  | node (l r : PT)

namespace Mono

def deg (mo : Mono) : ℕ := mo.e0 + mo.e1 + mo.e2 + mo.e3 + mo.e4 + mo.e5

/-- Kronecker position (base `m`, dehomogenized: the last exponent is ignored). -/
def pos (m : ℕ) (mo : Mono) : ℕ := mo.e0 + m * (mo.e1 + m * (mo.e2 + m * (mo.e3 + m * mo.e4)))

/-- the exponent vector as a finitely supported function -/
noncomputable def vec (mo : Mono) : Fin 6 →₀ ℕ :=
  Finsupp.single 0 mo.e0 + Finsupp.single 1 mo.e1 + Finsupp.single 2 mo.e2 +
    Finsupp.single 3 mo.e3 + Finsupp.single 4 mo.e4 + Finsupp.single 5 mo.e5

end Mono

/-! ## Denotations -/

noncomputable def monoZ (mo : Mono) : ZP :=
  X 0 ^ mo.e0 * X 1 ^ mo.e1 * X 2 ^ mo.e2 * X 3 ^ mo.e3 * X 4 ^ mo.e4 * X 5 ^ mo.e5

noncomputable def toZ : PT → ZP
  | .leaf mo c => C c * monoZ mo
  | .node l r => toZ l + toZ r

/-- the complex denotation, in the shape used by certificate modules -/
noncomputable def toC : PT → Poly
  | .leaf mo c => C (c : ℂ) * (X 0 ^ mo.e0 * X 1 ^ mo.e1 * X 2 ^ mo.e2 * X 3 ^ mo.e3 * X 4 ^ mo.e4 * X 5 ^ mo.e5)
  | .node l r => toC l + toC r

theorem toC_eq_map (t : PT) : toC t = MvPolynomial.map (Int.castRingHom ℂ) (toZ t) := by
  induction t with
  | leaf mo c => simp [toC, toZ, monoZ, MvPolynomial.map_C, MvPolynomial.map_X]
  | node l r ihl ihr => simp [toC, toZ, ihl, ihr]

theorem monoZ_eq_monomial (mo : Mono) : monoZ mo = monomial mo.vec 1 := by
  simp only [monoZ, Mono.vec, X_pow_eq_monomial, monomial_mul, mul_one]

/-! ## Kronecker evaluation -/

/-- `v * 2^s` through a kernel-accelerated shift. -/
def shl (v : ℤ) (s : ℕ) : ℤ := v * ((1 <<< s : ℕ) : ℤ)

theorem shl_eq (v : ℤ) (s : ℕ) : shl v s = v * 2 ^ s := by
  simp [shl, Nat.shiftLeft_eq]

/-- Value relative to the minimal position of the subtree, together with that position. -/
def psiRel (m W : ℕ) : PT → ℤ × ℕ
  | .leaf mo c => (c, mo.pos m)
  | .node l r =>
    let a := psiRel m W l
    let b := psiRel m W r
    let p := min a.2 b.2
    (shl a.1 (W * (a.2 - p)) + shl b.1 (W * (b.2 - p)), p)

def psi (m W : ℕ) (t : PT) : ℤ := shl (psiRel m W t).1 (W * (psiRel m W t).2)

/-- the Kronecker point -/
noncomputable def kpt (m W : ℕ) : Fin 6 → ℤ :=
  ![2 ^ W, 2 ^ (W * m), 2 ^ (W * m ^ 2), 2 ^ (W * m ^ 3), 2 ^ (W * m ^ 4), 1]

theorem kpt_zero (m W : ℕ) : kpt m W 0 = 2 ^ W := rfl
theorem kpt_one (m W : ℕ) : kpt m W 1 = 2 ^ (W * m) := rfl
theorem kpt_two (m W : ℕ) : kpt m W 2 = 2 ^ (W * m ^ 2) := rfl
theorem kpt_three (m W : ℕ) : kpt m W 3 = 2 ^ (W * m ^ 3) := rfl
theorem kpt_four (m W : ℕ) : kpt m W 4 = 2 ^ (W * m ^ 4) := rfl
theorem kpt_five (m W : ℕ) : kpt m W 5 = 1 := rfl

theorem eval_monoZ (m W : ℕ) (mo : Mono) : eval (kpt m W) (monoZ mo) = 2 ^ (W * mo.pos m) := by
  simp only [monoZ, map_mul, map_pow, eval_X, kpt_zero, kpt_one, kpt_two, kpt_three, kpt_four, kpt_five,
    Mono.pos, one_pow, mul_one, ← pow_mul, ← pow_add]
  congr 1
  ring

theorem psiRel_spec (m W : ℕ) (t : PT) :
    (psiRel m W t).1 * 2 ^ (W * (psiRel m W t).2) = eval (kpt m W) (toZ t) := by
  induction t with
  | leaf mo c => simp [psiRel, toZ, eval_monoZ, eval_C]
  | node l r ihl ihr =>
    simp only [psiRel, toZ, map_add]
    rw [← ihl, ← ihr]
    simp only [shl_eq, add_mul, mul_assoc, ← pow_add]
    have h1 : W * ((psiRel m W l).2 - min (psiRel m W l).2 (psiRel m W r).2) +
        W * min (psiRel m W l).2 (psiRel m W r).2 = W * (psiRel m W l).2 := by
      rw [← Nat.mul_add, Nat.sub_add_cancel (Nat.min_le_left _ _)]
    have h2 : W * ((psiRel m W r).2 - min (psiRel m W l).2 (psiRel m W r).2) +
        W * min (psiRel m W l).2 (psiRel m W r).2 = W * (psiRel m W r).2 := by
      rw [← Nat.mul_add, Nat.sub_add_cancel (Nat.min_le_right _ _)]
    rw [h1, h2]

theorem eval_kpt_toZ (m W : ℕ) (t : PT) : eval (kpt m W) (toZ t) = psi m W t := by
  rw [psi, shl_eq, psiRel_spec]

/-! ## Partial derivatives -/

def pdLeaf (i : Fin 6) (mo : Mono) (c : ℤ) : PT :=
  match i.val with
  | 0 => .leaf { mo with e0 := mo.e0 - 1 } (c * mo.e0)
  | 1 => .leaf { mo with e1 := mo.e1 - 1 } (c * mo.e1)
  | 2 => .leaf { mo with e2 := mo.e2 - 1 } (c * mo.e2)
  | 3 => .leaf { mo with e3 := mo.e3 - 1 } (c * mo.e3)
  | 4 => .leaf { mo with e4 := mo.e4 - 1 } (c * mo.e4)
  | _ => .leaf { mo with e5 := mo.e5 - 1 } (c * mo.e5)

def pderivPT (i : Fin 6) : PT → PT
  | .leaf mo c => pdLeaf i mo c
  | .node l r => .node (pderivPT i l) (pderivPT i r)

theorem pderiv_X_pow_self' (i : Fin 6) (n : ℕ) :
    pderiv i ((X i : ZP) ^ n) = (n : ZP) * X i ^ (n - 1) := by
  rw [pderiv_pow, pderiv_X_self, mul_one]

theorem pderiv_X_pow_ne' {i j : Fin 6} (h : j ≠ i) (n : ℕ) : pderiv i ((X j : ZP) ^ n) = 0 := by
  rw [pderiv_pow, pderiv_X_of_ne h, mul_zero]

theorem toZ_pdLeaf_zero (mo : Mono) (c : ℤ) : toZ (pdLeaf 0 mo c) = pderiv 0 (toZ (.leaf mo c)) := by
  rcases mo with ⟨a0, a1, a2, a3, a4, a5⟩
  simp only [pdLeaf, Fin.val_zero, toZ, monoZ, pderiv_mul, pderiv_C, pderiv_X_pow_self',
    pderiv_X_pow_ne' (show (1 : Fin 6) ≠ 0 by decide), pderiv_X_pow_ne' (show (2 : Fin 6) ≠ 0 by decide), pderiv_X_pow_ne' (show (3 : Fin 6) ≠ 0 by decide), pderiv_X_pow_ne' (show (4 : Fin 6) ≠ 0 by decide), pderiv_X_pow_ne' (show (5 : Fin 6) ≠ 0 by decide),
    mul_zero, zero_mul, add_zero, zero_add, C_mul, C_eq_coe_nat]
  ring

theorem toZ_pdLeaf_one (mo : Mono) (c : ℤ) : toZ (pdLeaf 1 mo c) = pderiv 1 (toZ (.leaf mo c)) := by
  rcases mo with ⟨a0, a1, a2, a3, a4, a5⟩
  simp only [pdLeaf, Fin.val_one, toZ, monoZ, pderiv_mul, pderiv_C, pderiv_X_pow_self',
    pderiv_X_pow_ne' (show (0 : Fin 6) ≠ 1 by decide), pderiv_X_pow_ne' (show (2 : Fin 6) ≠ 1 by decide), pderiv_X_pow_ne' (show (3 : Fin 6) ≠ 1 by decide), pderiv_X_pow_ne' (show (4 : Fin 6) ≠ 1 by decide), pderiv_X_pow_ne' (show (5 : Fin 6) ≠ 1 by decide),
    mul_zero, zero_mul, add_zero, zero_add, C_mul, C_eq_coe_nat]
  ring

theorem toZ_pdLeaf_two (mo : Mono) (c : ℤ) : toZ (pdLeaf 2 mo c) = pderiv 2 (toZ (.leaf mo c)) := by
  rcases mo with ⟨a0, a1, a2, a3, a4, a5⟩
  simp only [pdLeaf, show ((2 : Fin 6).val = 2) from rfl, toZ, monoZ, pderiv_mul, pderiv_C, pderiv_X_pow_self',
    pderiv_X_pow_ne' (show (0 : Fin 6) ≠ 2 by decide), pderiv_X_pow_ne' (show (1 : Fin 6) ≠ 2 by decide), pderiv_X_pow_ne' (show (3 : Fin 6) ≠ 2 by decide), pderiv_X_pow_ne' (show (4 : Fin 6) ≠ 2 by decide), pderiv_X_pow_ne' (show (5 : Fin 6) ≠ 2 by decide),
    mul_zero, zero_mul, add_zero, zero_add, C_mul, C_eq_coe_nat]
  ring

theorem toZ_pdLeaf_three (mo : Mono) (c : ℤ) : toZ (pdLeaf 3 mo c) = pderiv 3 (toZ (.leaf mo c)) := by
  rcases mo with ⟨a0, a1, a2, a3, a4, a5⟩
  simp only [pdLeaf, show ((3 : Fin 6).val = 3) from rfl, toZ, monoZ, pderiv_mul, pderiv_C, pderiv_X_pow_self',
    pderiv_X_pow_ne' (show (0 : Fin 6) ≠ 3 by decide), pderiv_X_pow_ne' (show (1 : Fin 6) ≠ 3 by decide), pderiv_X_pow_ne' (show (2 : Fin 6) ≠ 3 by decide), pderiv_X_pow_ne' (show (4 : Fin 6) ≠ 3 by decide), pderiv_X_pow_ne' (show (5 : Fin 6) ≠ 3 by decide),
    mul_zero, zero_mul, add_zero, zero_add, C_mul, C_eq_coe_nat]
  ring

theorem toZ_pdLeaf_four (mo : Mono) (c : ℤ) : toZ (pdLeaf 4 mo c) = pderiv 4 (toZ (.leaf mo c)) := by
  rcases mo with ⟨a0, a1, a2, a3, a4, a5⟩
  simp only [pdLeaf, show ((4 : Fin 6).val = 4) from rfl, toZ, monoZ, pderiv_mul, pderiv_C, pderiv_X_pow_self',
    pderiv_X_pow_ne' (show (0 : Fin 6) ≠ 4 by decide), pderiv_X_pow_ne' (show (1 : Fin 6) ≠ 4 by decide), pderiv_X_pow_ne' (show (2 : Fin 6) ≠ 4 by decide), pderiv_X_pow_ne' (show (3 : Fin 6) ≠ 4 by decide), pderiv_X_pow_ne' (show (5 : Fin 6) ≠ 4 by decide),
    mul_zero, zero_mul, add_zero, zero_add, C_mul, C_eq_coe_nat]
  ring

theorem toZ_pdLeaf_five (mo : Mono) (c : ℤ) : toZ (pdLeaf 5 mo c) = pderiv 5 (toZ (.leaf mo c)) := by
  rcases mo with ⟨a0, a1, a2, a3, a4, a5⟩
  simp only [pdLeaf, show ((5 : Fin 6).val = 5) from rfl, toZ, monoZ, pderiv_mul, pderiv_C, pderiv_X_pow_self',
    pderiv_X_pow_ne' (show (0 : Fin 6) ≠ 5 by decide), pderiv_X_pow_ne' (show (1 : Fin 6) ≠ 5 by decide), pderiv_X_pow_ne' (show (2 : Fin 6) ≠ 5 by decide), pderiv_X_pow_ne' (show (3 : Fin 6) ≠ 5 by decide), pderiv_X_pow_ne' (show (4 : Fin 6) ≠ 5 by decide),
    mul_zero, zero_mul, add_zero, zero_add, C_mul, C_eq_coe_nat]
  ring

theorem toZ_pdLeaf (i : Fin 6) (mo : Mono) (c : ℤ) : toZ (pdLeaf i mo c) = pderiv i (toZ (.leaf mo c)) := by
  fin_cases i
  · exact toZ_pdLeaf_zero mo c
  · exact toZ_pdLeaf_one mo c
  · exact toZ_pdLeaf_two mo c
  · exact toZ_pdLeaf_three mo c
  · exact toZ_pdLeaf_four mo c
  · exact toZ_pdLeaf_five mo c

theorem toZ_pderivPT (i : Fin 6) (t : PT) : toZ (pderivPT i t) = pderiv i (toZ t) := by
  induction t with
  | leaf mo c => exact toZ_pdLeaf i mo c
  | node l r ihl ihr => simp [pderivPT, toZ, ihl, ihr]

/-! ## Size data: ℓ¹ norm, degree, homogeneity -/

def l1PT : PT → ℕ
  | .leaf _ c => c.natAbs
  | .node l r => l1PT l + l1PT r

def degPT : PT → ℕ
  | .leaf mo _ => mo.deg
  | .node l _ => degPT l

def homB (n : ℕ) : PT → Bool
  | .leaf mo c => (c == 0) || (mo.deg == n)
  | .node l r => homB n l && homB n r

theorem monoZ_isHomogeneous (mo : Mono) : (monoZ mo).IsHomogeneous mo.deg := by
  have h := (((((isHomogeneous_X ℤ (0 : Fin 6)).pow mo.e0).mul ((isHomogeneous_X ℤ 1).pow mo.e1)).mul
    ((isHomogeneous_X ℤ 2).pow mo.e2)).mul ((isHomogeneous_X ℤ 3).pow mo.e3)).mul
    ((isHomogeneous_X ℤ 4).pow mo.e4) |>.mul ((isHomogeneous_X ℤ 5).pow mo.e5)
  simpa [monoZ, Mono.deg] using h

theorem toZ_isHomogeneous (n : ℕ) (t : PT) (h : homB n t = true) : (toZ t).IsHomogeneous n := by
  induction t with
  | leaf mo c =>
    simp only [homB, Bool.or_eq_true, beq_iff_eq] at h
    rcases h with h | h
    · subst h; simp only [toZ, C_0, zero_mul]; exact isHomogeneous_zero _ _ _
    · subst h; simpa [toZ] using (isHomogeneous_C _ c).mul (monoZ_isHomogeneous mo)
  | node l r ihl ihr =>
    simp only [homB, Bool.and_eq_true] at h
    exact (ihl h.1).add (ihr h.2)

/-- absolute values of the coefficients -/
def absPT : PT → PT
  | .leaf mo c => .leaf mo (c.natAbs : ℤ)
  | .node l r => .node (absPT l) (absPT r)

theorem coeff_monoZ_nonneg (mo : Mono) (d : Fin 6 →₀ ℕ) : 0 ≤ coeff d (monoZ mo) := by
  classical
  rw [monoZ_eq_monomial, coeff_monomial]
  split_ifs <;> simp

theorem coeff_toZ_absPT_nonneg (t : PT) (d : Fin 6 →₀ ℕ) : 0 ≤ coeff d (toZ (absPT t)) := by
  induction t with
  | leaf mo c =>
    simp only [absPT, toZ, coeff_C_mul]
    exact mul_nonneg (by positivity) (coeff_monoZ_nonneg mo d)
  | node l r ihl ihr => simp only [absPT, toZ, coeff_add]; exact add_nonneg (ihl) (ihr)

theorem abs_coeff_toZ_le (t : PT) (d : Fin 6 →₀ ℕ) : |coeff d (toZ t)| ≤ coeff d (toZ (absPT t)) := by
  induction t with
  | leaf mo c =>
    simp only [absPT, toZ, coeff_C_mul, abs_mul, abs_of_nonneg (coeff_monoZ_nonneg mo d)]
    simp
  | node l r ihl ihr =>
    simp only [absPT, toZ, coeff_add]
    exact (abs_add_le _ _).trans (add_le_add ihl ihr)

theorem eval_one_toZ_absPT (t : PT) : eval (fun _ => (1 : ℤ)) (toZ (absPT t)) = (l1PT t : ℤ) := by
  induction t with
  | leaf mo c => simp [absPT, toZ, monoZ, l1PT]
  | node l r ihl ihr => simp [absPT, toZ, l1PT, ihl, ihr]

/-! ## Expressions -/

inductive E where
  | poly (t : PT)
  | pd (i : Fin 6) (t : PT)
  | mul (a b : E)
  | add (a b : E)
  | sub (a b : E)
  | smulN (c : ℕ) (a : E)
  | pow (a : E) (n : ℕ)

noncomputable def denZ : E → ZP
  | .poly t => toZ t
  | .pd i t => pderiv i (toZ t)
  | .mul a b => denZ a * denZ b
  | .add a b => denZ a + denZ b
  | .sub a b => denZ a - denZ b
  | .smulN c a => (c : ZP) * denZ a
  | .pow a n => denZ a ^ n

noncomputable def absE : E → ZP
  | .poly t => toZ (absPT t)
  | .pd i t => toZ (absPT (pderivPT i t))
  | .mul a b => absE a * absE b
  | .add a b => absE a + absE b
  | .sub a b => absE a + absE b
  | .smulN c a => (c : ZP) * absE a
  | .pow a n => absE a ^ n

def psiE (m W : ℕ) : E → ℤ
  | .poly t => psi m W t
  | .pd i t => psi m W (pderivPT i t)
  | .mul a b => psiE m W a * psiE m W b
  | .add a b => psiE m W a + psiE m W b
  | .sub a b => psiE m W a - psiE m W b
  | .smulN c a => (c : ℤ) * psiE m W a
  | .pow a n => psiE m W a ^ n

def l1E : E → ℕ
  | .poly t => l1PT t
  | .pd i t => l1PT (pderivPT i t)
  | .mul a b => l1E a * l1E b
  | .add a b => l1E a + l1E b
  | .sub a b => l1E a + l1E b
  | .smulN c a => c * l1E a
  | .pow a n => l1E a ^ n

def degE : E → ℕ
  | .poly t => degPT t
  | .pd _ t => degPT t - 1
  | .mul a b => degE a + degE b
  | .add a _ => degE a
  | .sub a _ => degE a
  | .smulN _ a => degE a
  | .pow a n => degE a * n

def homE : E → Bool
  | .poly t => homB (degPT t) t
  | .pd i t => homB (degPT t - 1) (pderivPT i t)
  | .mul a b => homE a && homE b
  | .add a b => homE a && homE b && (degE a == degE b)
  | .sub a b => homE a && homE b && (degE a == degE b)
  | .smulN _ a => homE a
  | .pow a _ => homE a

theorem eval_kpt_denZ (m W : ℕ) (e : E) : eval (kpt m W) (denZ e) = psiE m W e := by
  induction e with
  | poly t => simp [denZ, psiE, eval_kpt_toZ]
  | pd i t => simp [denZ, psiE, ← toZ_pderivPT, eval_kpt_toZ]
  | mul a b iha ihb => simp [denZ, psiE, iha, ihb]
  | add a b iha ihb => simp [denZ, psiE, iha, ihb]
  | sub a b iha ihb => simp [denZ, psiE, iha, ihb]
  | smulN c a iha => simp [denZ, psiE, iha]
  | pow a n iha => simp [denZ, psiE, iha]

theorem denZ_isHomogeneous (e : E) (h : homE e = true) : (denZ e).IsHomogeneous (degE e) := by
  induction e with
  | poly t => exact toZ_isHomogeneous _ _ h
  | pd i t => simp only [denZ, degE]; rw [← toZ_pderivPT]; exact toZ_isHomogeneous _ _ h
  | mul a b iha ihb =>
    simp only [homE, Bool.and_eq_true] at h
    exact (iha h.1).mul (ihb h.2)
  | add a b iha ihb =>
    simp only [homE, Bool.and_eq_true, beq_iff_eq] at h
    exact (iha h.1.1).add (h.2 ▸ ihb h.1.2)
  | sub a b iha ihb =>
    simp only [homE, Bool.and_eq_true, beq_iff_eq] at h
    exact (iha h.1.1).sub (h.2 ▸ ihb h.1.2)
  | smulN c a iha =>
    simp only [denZ, degE]
    rw [← C_eq_coe_nat]
    simpa using (isHomogeneous_C _ (c : ℤ)).mul (iha h)
  | pow a n iha => exact (iha h).pow n

theorem coeff_absE_nonneg (e : E) (d : Fin 6 →₀ ℕ) : 0 ≤ coeff d (absE e) := by
  classical
  induction e generalizing d with
  | poly t => exact coeff_toZ_absPT_nonneg t d
  | pd i t => exact coeff_toZ_absPT_nonneg _ d
  | mul a b iha ihb =>
    simp only [absE, coeff_mul]
    exact Finset.sum_nonneg fun x _ => mul_nonneg (iha _) (ihb _)
  | add a b iha ihb => simp only [absE, coeff_add]; exact add_nonneg (iha d) (ihb d)
  | sub a b iha ihb => simp only [absE, coeff_add]; exact add_nonneg (iha d) (ihb d)
  | smulN c a iha =>
    simp only [absE]; rw [← C_eq_coe_nat, coeff_C_mul]; exact mul_nonneg (by positivity) (iha d)
  | pow a n iha =>
    induction n generalizing d with
    | zero => simp only [absE, pow_zero]; rw [← C_1, coeff_C]; split_ifs <;> simp
    | succ k ih =>
      simp only [absE, pow_succ, coeff_mul] at ih ⊢
      exact Finset.sum_nonneg fun x _ => mul_nonneg (ih _) (iha _)

theorem abs_coeff_mul_le {P P' Q Q' : ZP} (hP' : ∀ d, 0 ≤ coeff d P') (hQ' : ∀ d, 0 ≤ coeff d Q')
    (hP : ∀ d, |coeff d P| ≤ coeff d P') (hQ : ∀ d, |coeff d Q| ≤ coeff d Q') (d : Fin 6 →₀ ℕ) :
    |coeff d (P * Q)| ≤ coeff d (P' * Q') := by
  classical
  simp only [coeff_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun x _ => ?_)
  rw [abs_mul]
  exact mul_le_mul (hP _) (hQ _) (abs_nonneg _) (hP' _)

theorem abs_coeff_denZ_le (e : E) (d : Fin 6 →₀ ℕ) : |coeff d (denZ e)| ≤ coeff d (absE e) := by
  induction e generalizing d with
  | poly t => exact abs_coeff_toZ_le t d
  | pd i t => simp only [denZ, absE]; rw [← toZ_pderivPT]; exact abs_coeff_toZ_le _ d
  | mul a b iha ihb =>
    exact abs_coeff_mul_le (coeff_absE_nonneg a) (coeff_absE_nonneg b) iha ihb d
  | add a b iha ihb =>
    simp only [denZ, absE, coeff_add]
    exact (abs_add_le _ _).trans (add_le_add (iha d) (ihb d))
  | sub a b iha ihb =>
    simp only [denZ, absE, coeff_sub, coeff_add]
    exact (abs_sub _ _).trans (add_le_add (iha d) (ihb d))
  | smulN c a iha =>
    simp only [denZ, absE]
    rw [← C_eq_coe_nat, coeff_C_mul, coeff_C_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℤ) ≤ c)]
    exact mul_le_mul_of_nonneg_left (iha d) (by positivity)
  | pow a n iha =>
    induction n generalizing d with
    | zero =>
      simp only [denZ, absE, pow_zero]; rw [← C_1, coeff_C]; split_ifs <;> simp
    | succ k ih =>
      simp only [denZ, absE, pow_succ] at ih ⊢
      exact abs_coeff_mul_le (fun d => by
          have := coeff_absE_nonneg (.pow a k) d; simpa [absE] using this)
        (coeff_absE_nonneg a) ih iha d

theorem eval_one_absE (e : E) : eval (fun _ => (1 : ℤ)) (absE e) = (l1E e : ℤ) := by
  induction e with
  | poly t => simp [absE, l1E, eval_one_toZ_absPT]
  | pd i t => simp [absE, l1E, eval_one_toZ_absPT]
  | mul a b iha ihb => simp [absE, l1E, iha, ihb]
  | add a b iha ihb => simp [absE, l1E, iha, ihb]
  | sub a b iha ihb => simp [absE, l1E, iha, ihb]
  | smulN c a iha => simp [absE, l1E, iha]
  | pow a n iha => simp [absE, l1E, iha]

theorem coeff_le_eval_one {P : ZP} (h : ∀ d, 0 ≤ coeff d P) (d : Fin 6 →₀ ℕ) :
    coeff d P ≤ eval (fun _ => (1 : ℤ)) P := by
  classical
  rw [eval_eq']
  simp only [one_pow, Finset.prod_const_one, mul_one]
  by_cases hd : d ∈ P.support
  · exact Finset.single_le_sum (fun i _ => h i) hd
  · rw [notMem_support_iff.mp hd]
    exact Finset.sum_nonneg fun i _ => h i

theorem abs_coeff_denZ_le_l1 (e : E) (d : Fin 6 →₀ ℕ) : |coeff d (denZ e)| ≤ (l1E e : ℤ) := by
  rw [← eval_one_absE]
  exact (abs_coeff_denZ_le e d).trans (coeff_le_eval_one (coeff_absE_nonneg e) d)

/-! ## Injectivity of the Kronecker evaluation on bounded homogeneous polynomials -/

/-- the position of an exponent vector -/
def posF (m : ℕ) (d : Fin 6 →₀ ℕ) : ℕ := d 0 + m * (d 1 + m * (d 2 + m * (d 3 + m * d 4)))

theorem prod_kpt_pow (m W : ℕ) (d : Fin 6 →₀ ℕ) : ∏ i, kpt m W i ^ d i = 2 ^ (W * posF m d) := by
  rw [Fin.prod_univ_six]
  simp only [kpt_zero, kpt_one, kpt_two, kpt_three, kpt_four, kpt_five, posF, one_pow, mul_one,
    ← pow_mul, ← pow_add]
  congr 1
  ring

theorem mr_head {m a b x y : ℕ} (ha : a < m) (hb : b < m) (h : a + m * x = b + m * y) : a = b ∧ x = y := by
  have hm : 0 < m := lt_of_le_of_lt (Nat.zero_le _) ha
  have h1 : (a + m * x) % m = a := by rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt ha]
  have h2 : (b + m * y) % m = b := by rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hb]
  have hab : a = b := by rw [← h1, h, h2]
  refine ⟨hab, ?_⟩
  subst hab
  have := Nat.add_left_cancel h
  exact Nat.eq_of_mul_eq_mul_left hm this

theorem posF_injective {m n : ℕ} (hn : n < m) {d d' : Fin 6 →₀ ℕ}
    (hd : d 0 + d 1 + d 2 + d 3 + d 4 + d 5 = n) (hd' : d' 0 + d' 1 + d' 2 + d' 3 + d' 4 + d' 5 = n)
    (h : posF m d = posF m d') : d = d' := by
  simp only [posF] at h
  obtain ⟨h0, h⟩ := mr_head (by omega) (by omega) h
  obtain ⟨h1, h⟩ := mr_head (by omega) (by omega) h
  obtain ⟨h2, h⟩ := mr_head (by omega) (by omega) h
  obtain ⟨h3, h4⟩ := mr_head (by omega) (by omega) h
  have h5 : d 5 = d' 5 := by omega
  ext i
  fin_cases i <;> assumption

theorem sum_six_of_isHomogeneous {P : ZP} {n : ℕ} (hP : P.IsHomogeneous n) {d : Fin 6 →₀ ℕ}
    (hd : coeff d P ≠ 0) : d 0 + d 1 + d 2 + d 3 + d 4 + d 5 = n := by
  have := hP hd
  rw [show (1 : Fin 6 → ℕ) = fun _ => 1 from rfl, ← Finsupp.degree_eq_weight_one, Finsupp.degree_eq_sum,
    Fin.sum_univ_six] at this
  exact this

theorem eq_zero_of_eval_kpt {m W n : ℕ} {P : ZP} (hP : P.IsHomogeneous n) (hn : n < m)
    (hc : ∀ d, |coeff d P| < 2 ^ W) (h : eval (kpt m W) P = 0) : P = 0 := by
  classical
  by_contra hne
  have hsupp : P.support.Nonempty := by
    rw [Finset.nonempty_iff_ne_empty, Ne, support_eq_empty]; exact hne
  obtain ⟨d0, hd0, hmin⟩ := Finset.exists_min_image P.support (posF m) hsupp
  rw [eval_eq'] at h
  simp only [prod_kpt_pow] at h
  rw [← Finset.add_sum_erase _ _ hd0] at h
  have hdvd : (2 : ℤ) ^ (W * (posF m d0 + 1)) ∣ ∑ x ∈ P.support.erase d0, coeff x P * 2 ^ (W * posF m x) := by
    refine Finset.dvd_sum fun x hx => ?_
    have hx' := Finset.mem_of_mem_erase hx
    have hne' : x ≠ d0 := Finset.ne_of_mem_erase hx
    have hle := hmin x hx'
    have hposne : posF m x ≠ posF m d0 := by
      intro heq
      exact hne' (posF_injective hn (sum_six_of_isHomogeneous hP (mem_support_iff.mp hx'))
        (sum_six_of_isHomogeneous hP (mem_support_iff.mp hd0)) heq)
    have hlt : posF m d0 + 1 ≤ posF m x := by omega
    exact Dvd.dvd.mul_left (pow_dvd_pow 2 (Nat.mul_le_mul_left W hlt)) _
  have hmain : (2 : ℤ) ^ (W * (posF m d0 + 1)) ∣ coeff d0 P * 2 ^ (W * posF m d0) := by
    have : coeff d0 P * 2 ^ (W * posF m d0) = -∑ x ∈ P.support.erase d0, coeff x P * 2 ^ (W * posF m x) := by
      linarith
    rw [this]
    exact (dvd_neg).mpr hdvd
  have h2 : (2 : ℤ) ^ (W * (posF m d0 + 1)) = 2 ^ (W * posF m d0) * 2 ^ W := by
    rw [← pow_add, Nat.mul_succ]
  rw [h2, mul_comm (coeff d0 P)] at hmain
  have hdiv : (2 : ℤ) ^ W ∣ coeff d0 P :=
    (mul_dvd_mul_iff_left (pow_ne_zero _ (by norm_num))).mp hmain
  have hz : coeff d0 P = 0 := Int.eq_zero_of_abs_lt_dvd hdiv (hc d0)
  exact (mem_support_iff.mp hd0) hz

/-! ## Soundness -/

theorem sound {m W : ℕ} (L R : E) (hL : homE L = true) (hR : homE R = true) (hdeg : degE L = degE R)
    (hm : degE L < m) (hb : l1E L + l1E R < 2 ^ W) (hpsi : psiE m W L = psiE m W R) :
    denZ L = denZ R := by
  have hhom : (denZ L - denZ R).IsHomogeneous (degE L) :=
    (denZ_isHomogeneous L hL).sub (hdeg ▸ denZ_isHomogeneous R hR)
  have hc : ∀ d, |coeff d (denZ L - denZ R)| < 2 ^ W := by
    intro d
    rw [coeff_sub]
    calc |coeff d (denZ L) - coeff d (denZ R)| ≤ |coeff d (denZ L)| + |coeff d (denZ R)| := abs_sub _ _
      _ ≤ (l1E L : ℤ) + (l1E R : ℤ) := add_le_add (abs_coeff_denZ_le_l1 L d) (abs_coeff_denZ_le_l1 R d)
      _ < 2 ^ W := by exact_mod_cast hb
  have hev : eval (kpt m W) (denZ L - denZ R) = 0 := by
    rw [map_sub, eval_kpt_denZ, eval_kpt_denZ, hpsi, sub_self]
  exact sub_eq_zero.mp (eq_zero_of_eval_kpt hhom hm hc hev)

/-! ## Transfer to complex coefficients -/

noncomputable def denC (e : E) : Poly := MvPolynomial.map (Int.castRingHom ℂ) (denZ e)

@[simp] theorem denC_poly (t : PT) : denC (.poly t) = toC t := by rw [denC, denZ, toC_eq_map]
@[simp] theorem denC_pd (i : Fin 6) (t : PT) : denC (.pd i t) = pderiv i (toC t) := by
  rw [denC, denZ, toC_eq_map, pderiv_map]
@[simp] theorem denC_mul (a b : E) : denC (.mul a b) = denC a * denC b := by simp [denC, denZ]
@[simp] theorem denC_add (a b : E) : denC (.add a b) = denC a + denC b := by simp [denC, denZ]
@[simp] theorem denC_sub (a b : E) : denC (.sub a b) = denC a - denC b := by simp [denC, denZ]
@[simp] theorem denC_smulN (c : ℕ) (a : E) : denC (.smulN c a) = (c : Poly) * denC a := by
  simp [denC, denZ]
@[simp] theorem denC_pow (a : E) (n : ℕ) : denC (.pow a n) = denC a ^ n := by simp [denC, denZ]

theorem soundC {m W : ℕ} (L R : E) (hL : homE L = true) (hR : homE R = true) (hdeg : degE L = degE R)
    (hm : degE L < m) (hb : l1E L + l1E R < 2 ^ W) (hpsi : psiE m W L = psiE m W R) :
    denC L = denC R := by
  rw [denC, denC, sound L R hL hR hdeg hm hb hpsi]

/-! ## Auxiliary facts about explicit trees used by certificate modules -/

theorem toC_isHomogeneous (n : ℕ) (t : PT) (h : homB n t = true) : (toC t).IsHomogeneous n := by
  rw [toC_eq_map]; exact (toZ_isHomogeneous n t h).map _

/-- integer evaluation of a tree -/
def evalPT (x : Fin 6 → ℤ) : PT → ℤ
  | .leaf mo c => c * (x 0 ^ mo.e0 * x 1 ^ mo.e1 * x 2 ^ mo.e2 * x 3 ^ mo.e3 * x 4 ^ mo.e4 * x 5 ^ mo.e5)
  | .node l r => evalPT x l + evalPT x r

theorem eval_toZ (x : Fin 6 → ℤ) (t : PT) : eval x (toZ t) = evalPT x t := by
  induction t with
  | leaf mo c => simp [toZ, monoZ, evalPT]
  | node l r ihl ihr => simp [toZ, evalPT, ihl, ihr]

theorem eval_toC (x : Fin 6 → ℤ) (t : PT) : eval (fun i => (x i : ℂ)) (toC t) = (evalPT x t : ℂ) := by
  induction t with
  | leaf mo c => simp [toC, evalPT]
  | node l r ihl ihr => simp [toC, evalPT, ihl, ihr]

theorem toC_ne_zero_of_evalPT (x : Fin 6 → ℤ) (t : PT) (h : evalPT x t ≠ 0) : toC t ≠ 0 := by
  intro h0
  have := congrArg (eval (fun i => (x i : ℂ))) h0
  rw [eval_toC, map_zero] at this
  exact h (by exact_mod_cast this)

end HodgeComputer.Kron
