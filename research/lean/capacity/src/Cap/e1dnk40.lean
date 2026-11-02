import HodgeComputer.Verify

set_option Elab.async false

set_option maxHeartbeats 2000000

set_option maxRecDepth 4000

set_option linter.unusedTactic false
set_option linter.unusedSimpArgs false
set_option linter.unnecessarySeqFocus false
set_option linter.unreachableTactic false

noncomputable section

namespace HodgeComputer.Cap.e1dnk40

open MvPolynomial HodgeTestbed

def F : Poly := (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(6 : ℕ)

def A : Poly := (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (538506 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (37 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (358968 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (500 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-213 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-163 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (308 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-109 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-176 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-131 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-210 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-77 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-77 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-119 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-106 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-142 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-128 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (182 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-334 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-121 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (114 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (292 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (558 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-116 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (110 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-137 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-77 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-91 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-67 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-97 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (93 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (248 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (172 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (281 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (383 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-486 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-170 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-116 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-157 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-226 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-346 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (131 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (125 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (131 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (119 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-344 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-425 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-329 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (45 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-215 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (222 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (213 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (222 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-255 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-23 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-212 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-95 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-93 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (94 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-180 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-22 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-199 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-116 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-86 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-345 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-101 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-49 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-426 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-249 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-220 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (129 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-391 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (109 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-74 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-157 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-250 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-91 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-112 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-119 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-418856 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (70 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (145 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-332 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-149 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-614 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-356 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-66 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-273 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-348 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (74 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-348 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-186 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-334 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-450 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-276 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-122 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-402 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-187 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-402 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (113 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (75 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-266 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-152 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-146 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-228 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-262 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-640 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-178 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-151 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-151 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-146 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-324 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-364 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-516 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-397 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-399 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-559 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-198 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-384 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-256 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-498 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-129 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-260 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (66 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (62 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-120 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-148 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (66 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (58 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (82 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-124 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-148 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (58 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-119 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-703 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-136 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-112 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-49 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-70 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-49 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-91 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-599 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-259 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-338 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-201 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-413 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-549 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-171 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-104 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-243 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-77 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-77 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (43 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (86 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-137 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-110 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (153 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (87 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (66 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-134 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-75 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (52 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (146 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (62 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (103 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (99 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (101 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-263 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (179386 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-239498 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-59873 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-321 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-414 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-302 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-83 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (103 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (100 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-386 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-343 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-420 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-272 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-75 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-263 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-180 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-170 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-166 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-584 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-240 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-272 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-49 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-635 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-475 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-311 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-335 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-321 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-231 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-276 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-306 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-154 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-369 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-290 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-628 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-641 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-651 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-350 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-198 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-303 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-247 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (100 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (101 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (104 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (100 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-43 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-46 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-43 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-478578 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-826 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-269 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-286 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-363 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-156 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-225 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-156 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-294 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-340 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-523 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-226 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-367 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-450 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-456 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-505 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-176 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-568 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-463 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-411 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-89 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-47 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-89 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-103 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-61 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-47 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-19 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-103 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-140 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-208 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-31 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-536 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-139 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-194 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-234 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-134 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-469 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-332 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-206 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-418958 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-496 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-422 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-86 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-116 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-90 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-478839 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-62 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-174 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-716 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-438 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-417 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-263 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-546 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-669 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-118 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-65 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-168 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-231 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-299463 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-406 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-359 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-332 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-311 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-404 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-260 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-315 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-89 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-217 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-174 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-473 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-315 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-316 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-343 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-194 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-293 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-239412 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ)

def D : Poly := (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))

def Q : Poly := 0

def N0 : Poly := (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (81 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-81 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-81 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))

def N1 : Poly := (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(7 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (66 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (21 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ)

def N2 : Poly := (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))

def N3 : Poly := (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (25 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-25 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-25 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (40 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))

def N4 : Poly := (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (21 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-49 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))

def N5 : Poly := (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))

def N : Fin 6 → Poly := fun i => if i.val = 0 then N0 else if i.val = 1 then N1 else if i.val = 2 then N2 else if i.val = 3 then N3 else if i.val = 4 then N4 else N5

private theorem derivative_0_0 : pderiv (0 : Fin 6) (X 0 : Poly) = 1 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_0_1 : pderiv (0 : Fin 6) (X 1 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_0_2 : pderiv (0 : Fin 6) (X 2 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_0_3 : pderiv (0 : Fin 6) (X 3 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_0_4 : pderiv (0 : Fin 6) (X 4 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_0_5 : pderiv (0 : Fin 6) (X 5 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_1_0 : pderiv (1 : Fin 6) (X 0 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_1_1 : pderiv (1 : Fin 6) (X 1 : Poly) = 1 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_1_2 : pderiv (1 : Fin 6) (X 2 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_1_3 : pderiv (1 : Fin 6) (X 3 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_1_4 : pderiv (1 : Fin 6) (X 4 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_1_5 : pderiv (1 : Fin 6) (X 5 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_2_0 : pderiv (2 : Fin 6) (X 0 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_2_1 : pderiv (2 : Fin 6) (X 1 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_2_2 : pderiv (2 : Fin 6) (X 2 : Poly) = 1 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_2_3 : pderiv (2 : Fin 6) (X 3 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_2_4 : pderiv (2 : Fin 6) (X 4 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_2_5 : pderiv (2 : Fin 6) (X 5 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_3_0 : pderiv (3 : Fin 6) (X 0 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_3_1 : pderiv (3 : Fin 6) (X 1 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_3_2 : pderiv (3 : Fin 6) (X 2 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_3_3 : pderiv (3 : Fin 6) (X 3 : Poly) = 1 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_3_4 : pderiv (3 : Fin 6) (X 4 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_3_5 : pderiv (3 : Fin 6) (X 5 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_4_0 : pderiv (4 : Fin 6) (X 0 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_4_1 : pderiv (4 : Fin 6) (X 1 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_4_2 : pderiv (4 : Fin 6) (X 2 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_4_3 : pderiv (4 : Fin 6) (X 3 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_4_4 : pderiv (4 : Fin 6) (X 4 : Poly) = 1 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_4_5 : pderiv (4 : Fin 6) (X 5 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_5_0 : pderiv (5 : Fin 6) (X 0 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_5_1 : pderiv (5 : Fin 6) (X 1 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_5_2 : pderiv (5 : Fin 6) (X 2 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_5_3 : pderiv (5 : Fin 6) (X 3 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_5_4 : pderiv (5 : Fin 6) (X 4 : Poly) = 0 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]
private theorem derivative_5_5 : pderiv (5 : Fin 6) (X 5 : Poly) = 1 := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]

theorem e1 : (1 : Poly) * directional F N = A*D := by
  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]
  change (1:Poly) * (pderiv (0 : Fin 6) F * N0 + (pderiv (1 : Fin 6) F * N1 + (pderiv (2 : Fin 6) F * N2 + (pderiv (3 : Fin 6) F * N3 + (pderiv (4 : Fin 6) F * N4 + (pderiv (5 : Fin 6) F * N5)))))) = A*D
  simp only [F,A,D,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,derivative_0_0, derivative_0_1, derivative_0_2, derivative_0_3, derivative_0_4, derivative_0_5, derivative_1_0, derivative_1_1, derivative_1_2, derivative_1_3, derivative_1_4, derivative_1_5, derivative_2_0, derivative_2_1, derivative_2_2, derivative_2_3, derivative_2_4, derivative_2_5, derivative_3_0, derivative_3_1, derivative_3_2, derivative_3_3, derivative_3_4, derivative_3_5, derivative_4_0, derivative_4_1, derivative_4_2, derivative_4_3, derivative_4_4, derivative_4_5, derivative_5_0, derivative_5_1, derivative_5_2, derivative_5_3, derivative_5_4, derivative_5_5]
  apply MvPolynomial.funext; intro x
  norm_num <;> ring

end HodgeComputer.Cap.e1dnk40
