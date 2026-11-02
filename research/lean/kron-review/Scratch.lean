import HodgeComputer.Reflection.Kronecker
import HodgeComputer.Reflection.KroneckerTest
open HodgeComputer HodgeComputer.Kronecker
noncomputable section
-- exact check.py output for "-5*X0**6 + (3*X0**4*X1**2 + -X1**5*Y2)"
def Dgen : HodgeComputer.Poly := (((-(5 : HodgeComputer.Poly)) * ((MvPolynomial.X (R := ℂ) (0 : Fin 6)) ^ (6 : ℕ))) + ((((3 : HodgeComputer.Poly) * ((MvPolynomial.X (R := ℂ) (0 : Fin 6)) ^ (4 : ℕ))) * ((MvPolynomial.X (R := ℂ) (1 : Fin 6)) ^ (2 : ℕ))) + ((-((MvPolynomial.X (R := ℂ) (1 : Fin 6)) ^ (5 : ℕ))) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)))))
def wD : W := .add (.mul (.neg (.nat 5)) (.pow (.var 0) 6))
  (.add (.mul (.mul (.nat 3) (.pow (.var 0) 4)) (.pow (.var 1) 2)) (.mul (.neg (.pow (.var 1) 5)) (.var 5)))
def tD : PT := .node (.leaf ⟨6,0,0,0,0,0⟩ (-5)) (.node (.leaf ⟨4,2,0,0,0,0⟩ 3) (.leaf ⟨0,5,0,0,0,1⟩ (-1)))
example : Dgen = denW wD := rfl
example : denW wD = toC tD := matchW_sound wD tD (by decide +kernel)
-- numeral 1 / 0 cases
def Done : HodgeComputer.Poly := ((1 : HodgeComputer.Poly) * ((MvPolynomial.X (R := ℂ) (0 : Fin 6)) ^ (6 : ℕ)))
example : Done = denW (.mul (.nat 1) (.pow (.var 0) 6)) := rfl
end
#print axioms HodgeComputer.Kronecker.soundN
#print axioms HodgeComputer.Kronecker.sound4
#print axioms HodgeComputer.Kronecker.matchW_sound
#print axioms HodgeComputer.Kronecker.toC_isHomogeneous
#print axioms HodgeComputer.Kronecker.toC_ne_zero_of_evalPT
#print axioms HodgeComputer.KroneckerTest.small
#print axioms HodgeComputer.KroneckerTest.e1
#print axioms HodgeComputer.KroneckerTest.e2
#print axioms HodgeComputer.KroneckerTest.denominator_ne_zero
#print axioms HodgeComputer.KroneckerTest.denominator_homogeneous
