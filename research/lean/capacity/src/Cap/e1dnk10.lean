import HodgeComputer.Verify

set_option Elab.async false

set_option maxHeartbeats 2000000

set_option maxRecDepth 4000

set_option linter.unusedTactic false
set_option linter.unusedSimpArgs false
set_option linter.unnecessarySeqFocus false
set_option linter.unreachableTactic false

noncomputable section

namespace HodgeComputer.Cap.e1dnk10

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

def A : Poly := (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-718065 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (119757 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1012 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-213 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-464 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-240 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (116 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-120 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-145 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-145 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-252 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-876 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-83 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-122 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (99 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (38 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-159 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-134 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-507 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (200 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-446 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-705 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-558 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-971 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-885 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (203 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (271 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (184 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (181 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (295 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (248 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (27 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-74 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-246 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-428 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (142 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-347 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (100 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (50 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-128 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-126 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-666 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (166 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-57 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-263 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-396 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-544 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-436 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-780 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (122 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (37 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-148 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (94 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-288 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-254 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-275 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-521 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-678 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-324 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (129 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-236 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-514 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-660 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-522 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-250 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-398 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (286 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (253 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (160 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (163 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (163 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (197 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (27 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (480 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (384 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-80 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (152 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-240 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-160 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-408 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-128 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (65 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-29 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (239 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-410 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-441 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-228 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (285 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (66 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-614 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (111 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-89 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-140 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-230 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-379 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-310 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-656 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-159 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (55 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-139 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-206 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-191 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-350 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-471 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-391 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (166 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-135 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-326 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-343 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-255 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-321 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-375 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-572 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-446 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-108 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-299 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-314 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-236 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-308 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-503 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-660 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (396 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (325 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (296 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (108 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (142 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-227 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-50 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-331 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-71 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-822 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (147 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-383 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-100 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-261 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-502 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-620 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (243 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (214 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (127 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (162 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-146 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-153 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-399 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-409 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-275 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-339 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-248 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-289 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-165 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-721 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-672 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-524 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-492 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (179634 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (279 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (283 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (222 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (186 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (156 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (250 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (120 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (180 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (147 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (105 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (27 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-164 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (299156 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-134 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-206 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-250 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-182 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-62 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-86 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (160 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-104 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (64 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-236 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-468 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (168 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-246 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-506 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-140 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-174 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-342 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-194 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (180 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-178 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (166 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-368 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-578 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-216 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-496 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-190 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-174 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-712 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (114 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-210 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-446 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-324 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-576 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-404 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (138 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (165 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-285 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-479 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-106 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-326 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-674 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-870 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (182 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (306 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-164 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (62 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-400 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-160 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-250 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-808 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-416 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-244 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-306 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-248 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-362 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-498 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-566 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-108 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-354 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-46 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-168 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-256 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-292 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-490 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-638 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-120 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-124 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-128 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-152 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-540 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-156 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-160 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-448 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-220 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-404 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-348 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-118 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-104 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-124 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-154 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-158 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-318 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (66 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-416 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-226 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-540 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-602 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-119872 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-448 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-256 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-584 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-880 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-837756 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-538710 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (113 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-456 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-216 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-384 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-648 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-744 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (106 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (43 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (64 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (43 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (85 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (99 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-138 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-124 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-336 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-348 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-232 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-296 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-198 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-244 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-620 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-606 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-454 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-430 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (29 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (49 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (92 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (71 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (50 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (49 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (71 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (99 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (99 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (400 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (440 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (416 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (112 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (216 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-152 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (104 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (192 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-272 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-208 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (64 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-272 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-352 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-560 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-352 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-160 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(6 : ℕ)

def D : Poly := (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def Q : Poly := 0

def N0 : Poly := (MvPolynomial.C (σ := Fin 6) (60 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (60 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N1 : Poly := (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-25 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-25 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (35 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N2 : Poly := (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-81 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N3 : Poly := (MvPolynomial.C (σ := Fin 6) (-15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-13 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (34 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N4 : Poly := (MvPolynomial.C (σ := Fin 6) (10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (40 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N5 : Poly := (MvPolynomial.C (σ := Fin 6) (70 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-126 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-84 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (70 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-98 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

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

end HodgeComputer.Cap.e1dnk10
