import HodgeComputer.Verify

/-!
# Reflection checker for polynomial identities in `MvPolynomial (Fin 6) ℂ`

Computable sparse polynomials over ℤ (`SPoly`), with a denotation `eval : SPoly → Poly` and soundness
lemmas for the operations used by certificate identities (add, mul, scalar, partial derivative).
Two `SPoly` terms that reduce (in the kernel) to the same normalized list denote the same polynomial,
so an identity is proved by `congrArg eval (by decide +kernel)`.

Soundness never depends on sortedness or fuel sufficiency: when fuel runs out `addAux` appends,
which is still sound. Sortedness only matters for the final syntactic comparison to succeed.
-/

set_option autoImplicit false
set_option linter.unusedSimpArgs false

namespace HodgeComputer.Reflect

open MvPolynomial

/-- Exponent vector of a monomial in the six variables `X0 X1 X2 Y0 Y1 Y2`. -/
structure Mono where
  e0 : ℕ
  e1 : ℕ
  e2 : ℕ
  e3 : ℕ
  e4 : ℕ
  e5 : ℕ
  deriving DecidableEq, Repr

namespace Mono

def mul (a b : Mono) : Mono :=
  ⟨a.e0 + b.e0, a.e1 + b.e1, a.e2 + b.e2, a.e3 + b.e3, a.e4 + b.e4, a.e5 + b.e5⟩

/-- Strict lexicographic order (as a `Bool`), using kernel-accelerated `Nat` primitives. -/
def lt (a b : Mono) : Bool :=
  Nat.blt a.e0 b.e0 || (Nat.beq a.e0 b.e0 &&
   (Nat.blt a.e1 b.e1 || (Nat.beq a.e1 b.e1 &&
    (Nat.blt a.e2 b.e2 || (Nat.beq a.e2 b.e2 &&
     (Nat.blt a.e3 b.e3 || (Nat.beq a.e3 b.e3 &&
      (Nat.blt a.e4 b.e4 || (Nat.beq a.e4 b.e4 && Nat.blt a.e5 b.e5)))))))))

end Mono

abbrev Term := Mono × ℤ
abbrev SPoly := List Term

/-! ### Denotation -/

noncomputable def evalMono (m : Mono) : Poly :=
  X 0 ^ m.e0 * X 1 ^ m.e1 * X 2 ^ m.e2 * X 3 ^ m.e3 * X 4 ^ m.e4 * X 5 ^ m.e5

noncomputable def evalTerm (t : Term) : Poly := C (t.2 : ℂ) * evalMono t.1

noncomputable def eval : SPoly → Poly
  | [] => 0
  | t :: p => evalTerm t + eval p

@[simp] theorem eval_nil : eval [] = 0 := rfl
@[simp] theorem eval_cons (t : Term) (p : SPoly) : eval (t :: p) = evalTerm t + eval p := rfl

theorem eval_append (p q : SPoly) : eval (p ++ q) = eval p + eval q := by
  induction p with
  | nil => simp
  | cons t p ih => simp [ih, add_assoc]

theorem evalMono_mul (a b : Mono) : evalMono (Mono.mul a b) = evalMono a * evalMono b := by
  simp only [evalMono, Mono.mul, pow_add]; ring

theorem evalTerm_add (m : Mono) (c c' : ℤ) :
    evalTerm (m, c + c') = evalTerm (m, c) + evalTerm (m, c') := by
  simp only [evalTerm, Int.cast_add, C_add]; ring

theorem evalTerm_zero (m : Mono) : evalTerm (m, 0) = 0 := by
  simp [evalTerm]

theorem evalTerm_mul (t s : Term) :
    evalTerm (Mono.mul t.1 s.1, t.2 * s.2) = evalTerm t * evalTerm s := by
  simp only [evalTerm, Int.cast_mul, C_mul, evalMono_mul]; ring

/-! ### Addition by merging -/

def addAux : ℕ → SPoly → SPoly → SPoly
  | 0, p, q => p ++ q
  | _ + 1, [], q => q
  | _ + 1, p, [] => p
  | n + 1, (m, c) :: p, (m', c') :: q =>
    if Mono.lt m m' then (m, c) :: addAux n p ((m', c') :: q)
    else if Mono.lt m' m then (m', c') :: addAux n ((m, c) :: p) q
    else if c + c' = 0 then addAux n p q
    else (m, c + c') :: addAux n p q

def add (p q : SPoly) : SPoly := addAux (p.length + q.length) p q

theorem eval_addAux (n : ℕ) (p q : SPoly) : eval (addAux n p q) = eval p + eval q := by
  induction n generalizing p q with
  | zero => simp [addAux, eval_append]
  | succ n ih =>
    match p, q with
    | [], q => simp [addAux]
    | (m, c) :: p, [] => simp [addAux]
    | (m, c) :: p, (m', c') :: q =>
      simp only [addAux]
      split_ifs with h1 h2 h3
      · simp only [eval_cons, ih]; ring
      · simp only [eval_cons, ih]; ring
      · have hz : evalTerm (m, c) + evalTerm (m', c') = 0 := by
          -- without lt in either direction the monomials coincide
          have hm : m = m' := by
            rcases m with ⟨a0, a1, a2, a3, a4, a5⟩; rcases m' with ⟨b0, b1, b2, b3, b4, b5⟩
            simp only [Mono.lt, Bool.or_eq_true, Bool.and_eq_true, Nat.blt_eq, Nat.beq_eq] at h1 h2
            simp only [Mono.mk.injEq]; omega
          subst hm
          rw [← evalTerm_add, h3, evalTerm_zero]
        simp only [eval_cons, ih]
        linear_combination -hz
      · have hm : m = m' := by
          rcases m with ⟨a0, a1, a2, a3, a4, a5⟩; rcases m' with ⟨b0, b1, b2, b3, b4, b5⟩
          simp only [Mono.lt, Bool.or_eq_true, Bool.and_eq_true, Nat.blt_eq, Nat.beq_eq] at h1 h2
          simp only [Mono.mk.injEq]; omega
        subst hm
        simp only [eval_cons, ih, evalTerm_add]; ring

theorem eval_add (p q : SPoly) : eval (add p q) = eval p + eval q := eval_addAux _ p q

/-! ### Scalar multiplication, negation, subtraction -/

def smul (k : ℤ) (p : SPoly) : SPoly := p.map fun t => (t.1, k * t.2)

theorem eval_smul (k : ℤ) (p : SPoly) : eval (smul k p) = C (k : ℂ) * eval p := by
  induction p with
  | nil => simp [smul]
  | cons t p ih =>
    simp only [smul, List.map_cons, eval_cons, evalTerm, Int.cast_mul, C_mul] at ih ⊢
    rw [ih]; ring

theorem eval_smul_nat (k : ℕ) (p : SPoly) : eval (smul (k : ℤ) p) = (k : Poly) * eval p := by
  rw [eval_smul, Int.cast_natCast, map_natCast]

def sub (p q : SPoly) : SPoly := add p (smul (-1) q)

theorem eval_sub (p q : SPoly) : eval (sub p q) = eval p - eval q := by
  simp only [sub, eval_add, eval_smul, Int.cast_neg, Int.cast_one, C_neg, C_1]; ring

/-! ### Multiplication -/

def mulTerm (t : Term) (q : SPoly) : SPoly := q.map fun s => (Mono.mul t.1 s.1, t.2 * s.2)

theorem eval_mulTerm (t : Term) (q : SPoly) : eval (mulTerm t q) = evalTerm t * eval q := by
  induction q with
  | nil => simp [mulTerm]
  | cons s q ih =>
    simp only [mulTerm, List.map_cons, eval_cons] at ih ⊢
    rw [ih, evalTerm_mul]; ring

def mul (p q : SPoly) : SPoly := p.foldr (fun t acc => add (mulTerm t q) acc) []

theorem eval_mul (p q : SPoly) : eval (mul p q) = eval p * eval q := by
  induction p with
  | nil => simp [mul]
  | cons t p ih =>
    simp only [mul, List.foldr_cons, eval_add, eval_mulTerm, eval_cons] at ih ⊢
    rw [ih]; ring

def one : SPoly := [(⟨0, 0, 0, 0, 0, 0⟩, 1)]

theorem eval_one : eval one = 1 := by simp [one, evalTerm, evalMono]

def pow (p : SPoly) : ℕ → SPoly
  | 0 => one
  | n + 1 => mul p (pow p n)

theorem eval_pow (p : SPoly) (n : ℕ) : eval (pow p n) = eval p ^ n := by
  induction n with
  | zero => simp [pow, eval_one]
  | succ n ih => simp [pow, eval_mul, ih, pow_succ]; ring

/-! ### Partial derivatives -/

theorem pderiv_X_pow_self (i : Fin 6) (n : ℕ) :
    pderiv i ((X i : Poly) ^ n) = C (n : ℂ) * X i ^ (n - 1) := by
  rw [pderiv_pow, pderiv_X_self, mul_one, map_natCast]

theorem pderiv_X_pow_ne {i j : Fin 6} (h : j ≠ i) (n : ℕ) :
    pderiv i ((X j : Poly) ^ n) = 0 := by
  rw [pderiv_pow, pderiv_X_of_ne h, mul_zero]

def pderivTerm0 (t : Term) : SPoly :=
  if t.1.e0 = 0 then [] else [({ t.1 with e0 := t.1.e0 - 1 }, t.2 * t.1.e0)]
def pderivTerm1 (t : Term) : SPoly :=
  if t.1.e1 = 0 then [] else [({ t.1 with e1 := t.1.e1 - 1 }, t.2 * t.1.e1)]
def pderivTerm2 (t : Term) : SPoly :=
  if t.1.e2 = 0 then [] else [({ t.1 with e2 := t.1.e2 - 1 }, t.2 * t.1.e2)]
def pderivTerm3 (t : Term) : SPoly :=
  if t.1.e3 = 0 then [] else [({ t.1 with e3 := t.1.e3 - 1 }, t.2 * t.1.e3)]
def pderivTerm4 (t : Term) : SPoly :=
  if t.1.e4 = 0 then [] else [({ t.1 with e4 := t.1.e4 - 1 }, t.2 * t.1.e4)]
def pderivTerm5 (t : Term) : SPoly :=
  if t.1.e5 = 0 then [] else [({ t.1 with e5 := t.1.e5 - 1 }, t.2 * t.1.e5)]

def pderivTerm (i : Fin 6) : Term → SPoly :=
  if i.val = 0 then pderivTerm0 else if i.val = 1 then pderivTerm1 else if i.val = 2 then pderivTerm2
  else if i.val = 3 then pderivTerm3 else if i.val = 4 then pderivTerm4 else pderivTerm5

def pderivP (i : Fin 6) (p : SPoly) : SPoly := p.flatMap (pderivTerm i)

theorem pderiv_evalTerm0 (t : Term) : pderiv (0 : Fin 6) (evalTerm t) = eval (pderivTerm0 t) := by
  rcases t with ⟨⟨a, b, c, d, e, f⟩, k⟩
  simp only [evalTerm, evalMono, pderivTerm0, pderiv_mul, pderiv_C, pderiv_X_pow_self,
    pderiv_X_pow_ne (show (1 : Fin 6) ≠ 0 by decide), pderiv_X_pow_ne (show (2 : Fin 6) ≠ 0 by decide),
    pderiv_X_pow_ne (show (3 : Fin 6) ≠ 0 by decide), pderiv_X_pow_ne (show (4 : Fin 6) ≠ 0 by decide),
    pderiv_X_pow_ne (show (5 : Fin 6) ≠ 0 by decide), mul_zero, zero_mul, add_zero, zero_add]
  split_ifs with h
  · subst h; simp [pderiv_C]
  · simp only [eval_cons, eval_nil, evalTerm, evalMono, Int.cast_mul, Int.cast_natCast, C_mul, add_zero]; ring

theorem pderiv_evalTerm1 (t : Term) : pderiv (1 : Fin 6) (evalTerm t) = eval (pderivTerm1 t) := by
  rcases t with ⟨⟨a, b, c, d, e, f⟩, k⟩
  simp only [evalTerm, evalMono, pderivTerm1, pderiv_mul, pderiv_C, pderiv_X_pow_self,
    pderiv_X_pow_ne (show (0 : Fin 6) ≠ 1 by decide), pderiv_X_pow_ne (show (2 : Fin 6) ≠ 1 by decide),
    pderiv_X_pow_ne (show (3 : Fin 6) ≠ 1 by decide), pderiv_X_pow_ne (show (4 : Fin 6) ≠ 1 by decide),
    pderiv_X_pow_ne (show (5 : Fin 6) ≠ 1 by decide), mul_zero, zero_mul, add_zero, zero_add]
  split_ifs with h
  · subst h; simp [pderiv_C]
  · simp only [eval_cons, eval_nil, evalTerm, evalMono, Int.cast_mul, Int.cast_natCast, C_mul, add_zero]; ring

theorem pderiv_evalTerm2 (t : Term) : pderiv (2 : Fin 6) (evalTerm t) = eval (pderivTerm2 t) := by
  rcases t with ⟨⟨a, b, c, d, e, f⟩, k⟩
  simp only [evalTerm, evalMono, pderivTerm2, pderiv_mul, pderiv_C, pderiv_X_pow_self,
    pderiv_X_pow_ne (show (0 : Fin 6) ≠ 2 by decide), pderiv_X_pow_ne (show (1 : Fin 6) ≠ 2 by decide),
    pderiv_X_pow_ne (show (3 : Fin 6) ≠ 2 by decide), pderiv_X_pow_ne (show (4 : Fin 6) ≠ 2 by decide),
    pderiv_X_pow_ne (show (5 : Fin 6) ≠ 2 by decide), mul_zero, zero_mul, add_zero, zero_add]
  split_ifs with h
  · subst h; simp [pderiv_C]
  · simp only [eval_cons, eval_nil, evalTerm, evalMono, Int.cast_mul, Int.cast_natCast, C_mul, add_zero]; ring

theorem pderiv_evalTerm3 (t : Term) : pderiv (3 : Fin 6) (evalTerm t) = eval (pderivTerm3 t) := by
  rcases t with ⟨⟨a, b, c, d, e, f⟩, k⟩
  simp only [evalTerm, evalMono, pderivTerm3, pderiv_mul, pderiv_C, pderiv_X_pow_self,
    pderiv_X_pow_ne (show (0 : Fin 6) ≠ 3 by decide), pderiv_X_pow_ne (show (1 : Fin 6) ≠ 3 by decide),
    pderiv_X_pow_ne (show (2 : Fin 6) ≠ 3 by decide), pderiv_X_pow_ne (show (4 : Fin 6) ≠ 3 by decide),
    pderiv_X_pow_ne (show (5 : Fin 6) ≠ 3 by decide), mul_zero, zero_mul, add_zero, zero_add]
  split_ifs with h
  · subst h; simp [pderiv_C]
  · simp only [eval_cons, eval_nil, evalTerm, evalMono, Int.cast_mul, Int.cast_natCast, C_mul, add_zero]; ring

theorem pderiv_evalTerm4 (t : Term) : pderiv (4 : Fin 6) (evalTerm t) = eval (pderivTerm4 t) := by
  rcases t with ⟨⟨a, b, c, d, e, f⟩, k⟩
  simp only [evalTerm, evalMono, pderivTerm4, pderiv_mul, pderiv_C, pderiv_X_pow_self,
    pderiv_X_pow_ne (show (0 : Fin 6) ≠ 4 by decide), pderiv_X_pow_ne (show (1 : Fin 6) ≠ 4 by decide),
    pderiv_X_pow_ne (show (2 : Fin 6) ≠ 4 by decide), pderiv_X_pow_ne (show (3 : Fin 6) ≠ 4 by decide),
    pderiv_X_pow_ne (show (5 : Fin 6) ≠ 4 by decide), mul_zero, zero_mul, add_zero, zero_add]
  split_ifs with h
  · subst h; simp [pderiv_C]
  · simp only [eval_cons, eval_nil, evalTerm, evalMono, Int.cast_mul, Int.cast_natCast, C_mul, add_zero]; ring

theorem pderiv_evalTerm5 (t : Term) : pderiv (5 : Fin 6) (evalTerm t) = eval (pderivTerm5 t) := by
  rcases t with ⟨⟨a, b, c, d, e, f⟩, k⟩
  simp only [evalTerm, evalMono, pderivTerm5, pderiv_mul, pderiv_C, pderiv_X_pow_self,
    pderiv_X_pow_ne (show (0 : Fin 6) ≠ 5 by decide), pderiv_X_pow_ne (show (1 : Fin 6) ≠ 5 by decide),
    pderiv_X_pow_ne (show (2 : Fin 6) ≠ 5 by decide), pderiv_X_pow_ne (show (3 : Fin 6) ≠ 5 by decide),
    pderiv_X_pow_ne (show (4 : Fin 6) ≠ 5 by decide), mul_zero, zero_mul, add_zero, zero_add]
  split_ifs with h
  · subst h; simp [pderiv_C]
  · simp only [eval_cons, eval_nil, evalTerm, evalMono, Int.cast_mul, Int.cast_natCast, C_mul, add_zero]; ring

theorem pderiv_evalTerm (i : Fin 6) (t : Term) : pderiv i (evalTerm t) = eval (pderivTerm i t) := by
  fin_cases i
  · exact pderiv_evalTerm0 t
  · exact pderiv_evalTerm1 t
  · exact pderiv_evalTerm2 t
  · exact pderiv_evalTerm3 t
  · exact pderiv_evalTerm4 t
  · exact pderiv_evalTerm5 t

theorem eval_pderivP (i : Fin 6) (p : SPoly) : eval (pderivP i p) = pderiv i (eval p) := by
  induction p with
  | nil => simp [pderivP]
  | cons t p ih =>
    simp only [pderivP, List.flatMap_cons, eval_append, eval_cons, map_add] at ih ⊢
    rw [ih, pderiv_evalTerm]

/-! ### The verdict -/

theorem eval_eq_of_eq {p q : SPoly} (h : p = q) : eval p = eval q := congrArg eval h

end HodgeComputer.Reflect
