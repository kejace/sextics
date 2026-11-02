import HodgeComputer.Verify

set_option Elab.async false

set_option maxHeartbeats 2000000

set_option maxRecDepth 4000

set_option linter.unusedTactic false
set_option linter.unusedSimpArgs false
set_option linter.unnecessarySeqFocus false
set_option linter.unreachableTactic false

noncomputable section

namespace HodgeComputer.Cap.e1dnk20

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

def A : Poly := (MvPolynomial.C (σ := Fin 6) (-538548 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (836973 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-478480 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (359298 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (612 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (936 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-542 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-236 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (756 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-252 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (426 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (180 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-175 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-239 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-372 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-366 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-513 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-53 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (544 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (226 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (336 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (200 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-112 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (534 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (468 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (476 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (597 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (221 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (547 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (967 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-111 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-224 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (366 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (324 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (107 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-82 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-235 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (280 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-170 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (52 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (118 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (364 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (94 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (120 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (522 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (648 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1056 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (864 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (202 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (270 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (504 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (432 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (242 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (414 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (194 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (488 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (599 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (262 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (542 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (198 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (288 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (216 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (244 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (735 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1043 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (518 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (280 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (46 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (232 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (280 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (112 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (380 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (450 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (486 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (522 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (648 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-483 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (358 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (198 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (432 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-100 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-41 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-175 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (52 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (58 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (268 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (634 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (510 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (558 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (342 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (738 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (864 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (308 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (432 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (358 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (402 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (204 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (79 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (67 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (750 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (246 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (85 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-76 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-234 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (203 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-124 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-186 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-101 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (764 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (252 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (158 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (186 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (108 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (215 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (424 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (472 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (178 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (466 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-69 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (376 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (402 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (612 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (883 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (370 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (529 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (144 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (520 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (667 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (164 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (52 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (184 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (238 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (604 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (910 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (768 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (174 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (216 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-79 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (196 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (342 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (540 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (211 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (378 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (331 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (458 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-22 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (346 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (380 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (536 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (114 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (528 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (966 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1080 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (842 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (958 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (597 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (582 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1148 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (804 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (336 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (672 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (636 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (343 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (106 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (472 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (461 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (269 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (240 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-359109 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (594 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (312 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (504 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-148 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (159 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-47 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-29 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-250 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-306 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (170 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (452 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (680 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (128 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (428 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (576 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (300 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (336 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (702 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (792 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (473 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (609 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (520 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (495 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (369 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (268 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (184 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (924 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (327 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (178 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-31 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-45 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (284 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-57 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (114 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (128 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (426 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (188 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (146 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (696 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (248 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (200 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (512 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (410 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (732 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (714 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (312 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (600 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (312 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (864 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1224 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (182 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (538556 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (188 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (482 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (548 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (148 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (432 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (235 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (808 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (484 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1048 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (497 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (764 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (656 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (303 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (138 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-34 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (502 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (400 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (296 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (767 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (714 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (570 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (275 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (294 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (76 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (124 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (190 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (264 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (440 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (739 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (639 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (522 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (852 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (570 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (682 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (687 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (712 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (405 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (608 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (324 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (450 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (181 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (484 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (640 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (544 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (521 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (621 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (246 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (562 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (585 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (90 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (232 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (138 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (286 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (312 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (90 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (370 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (612 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (373 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1008 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (672 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (561 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (668 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (716 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (787 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (749 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (360 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (658 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (519 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (446 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (351 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (126 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-88 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (132 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (153 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (478592 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (462 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (333 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (369 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (87 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (249 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (237 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (376 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (238 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (232 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (184 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (346 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (352 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (633 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (775 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (401 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (672 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (507 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (475 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (408 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (810 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (486 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (381 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (112 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (167 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (467 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (279 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (228 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (104 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (202 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (152 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (116 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (120 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (126 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (98 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (150 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (90 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (138 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (134 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (80 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (538622 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (240010 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (120347 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (860 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (886 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (328 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (808 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (931 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (693 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (245 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (632 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (906 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (987 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (460 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (166 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (478872 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (998 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (821 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (592 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (372 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (255 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (418906 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (570 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (370 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (138 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (196 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (102 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(6 : ℕ)

def D : Poly := (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def Q : Poly := 0

def N0 : Poly := (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-70 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (85 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (84 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (112 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-70 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (166 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-23 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-84 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-98 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-84 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (154 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (112 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N1 : Poly := (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N2 : Poly := (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (30 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-82 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

def N3 : Poly := (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (78 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-66 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-45 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-20 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-27 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (63 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N4 : Poly := (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-35 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-21 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-42 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-49 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-34 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (49 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (64 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (56 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N5 : Poly := (MvPolynomial.C (σ := Fin 6) (-90 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (144 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-90 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (144 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-90 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-54 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-126 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (72 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-108 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (126 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (144 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-18 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ)

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

end HodgeComputer.Cap.e1dnk20
