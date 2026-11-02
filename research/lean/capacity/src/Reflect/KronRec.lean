import Cap.K_m3k100
/-! Experiment: recursor-style definitions (no brecOn tables) for the Kronecker value. -/
namespace HodgeComputer.Kron
open HodgeComputer.Cap.K_m3k100

noncomputable def psiRel' (m W : ℕ) (t : PT) : ℤ × ℕ :=
  PT.rec (motive := fun _ => ℤ × ℕ) (fun mo c => (c, mo.pos m))
    (fun _ _ a b => let p := min a.2 b.2; (shl a.1 (W * (a.2 - p)) + shl b.1 (W * (b.2 - p)), p)) t

noncomputable def psi' (m W : ℕ) (t : PT) : ℤ := let a := psiRel' m W t; shl a.1 (W * a.2)

def pdLeaf' (i : Fin 6) (mo : Mono) (c : ℤ) : PT := pdLeaf i mo c

noncomputable def pderivPT' (i : Fin 6) (t : PT) : PT :=
  PT.rec (motive := fun _ => PT) (fun mo c => pdLeaf i mo c) (fun _ _ l r => .node l r) t

noncomputable def psiE' (m W : ℕ) : E → ℤ
  | .poly t => psi' m W t
  | .pd i t => psi' m W (pderivPT' i t)
  | .mul a b => psiE' m W a * psiE' m W b
  | .add a b => psiE' m W a + psiE' m W b
  | .sub a b => psiE' m W a - psiE' m W b
  | .smulN c a => (c : ℤ) * psiE' m W a
  | .pow a n => psiE' m W a ^ n

example : psiE' 17 27 L2 = psiE' 17 27 R2 := by decide +kernel
end HodgeComputer.Kron
