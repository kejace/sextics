import HodgeComputer.Verify

set_option Elab.async false
set_option maxHeartbeats 2000000
set_option maxRecDepth 4000
set_option linter.unusedTactic false
set_option linter.unusedSimpArgs false
set_option linter.unnecessarySeqFocus false
set_option linter.unreachableTactic false

noncomputable section

namespace HodgeComputer.Instances.S_dense_d4_n126

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

def A : Poly := (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (45 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (22 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (13 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (90 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (64 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (45 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (27 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-52 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (39 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (75 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (51 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-34 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-68 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-64 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-96 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-25 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-28 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-60 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-36 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-45 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-34 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-34 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-78 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-64 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-21 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-42 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-33 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-15 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (56 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (30 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-48 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-14 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-11 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (40 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (65 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (26 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (16 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (17 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (10 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-24 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-32 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-18 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-34 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(6 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-70 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-44 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-20 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-12 / 9973 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(5 : ℕ)

def D : Poly := (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def Q : Poly := 0

def N0 : Poly := (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N1 : Poly := (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N2 : Poly := 0

def N3 : Poly := (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N4 : Poly := (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (0 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (1 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-2 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (7 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-9 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (6 : ℂ)) * (MvPolynomial.X (R := ℂ) (2 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (5 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-3 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-1 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-7 : ℂ)) * (MvPolynomial.X (R := ℂ) (3 : Fin 6)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-8 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(5 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (1 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(4 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6)) +
    (MvPolynomial.C (σ := Fin 6) (8 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(3 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(2 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-5 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6))^(2 : ℕ) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(3 : ℕ) +
    (MvPolynomial.C (σ := Fin 6) (-4 : ℂ)) * (MvPolynomial.X (R := ℂ) (4 : Fin 6)) * (MvPolynomial.X (R := ℂ) (5 : Fin 6))^(4 : ℕ)

def N5 : Poly := 0

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

theorem e2 : divergenceNumerator D N = F^1*Q := by
  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]
  change D * (pderiv (0 : Fin 6) N0 + (pderiv (1 : Fin 6) N1 + (pderiv (2 : Fin 6) N2 + (pderiv (3 : Fin 6) N3 + (pderiv (4 : Fin 6) N4 + (pderiv (5 : Fin 6) N5)))))) - (pderiv (0 : Fin 6) D * N0 + (pderiv (1 : Fin 6) D * N1 + (pderiv (2 : Fin 6) D * N2 + (pderiv (3 : Fin 6) D * N3 + (pderiv (4 : Fin 6) D * N4 + (pderiv (5 : Fin 6) D * N5)))))) = F^1*Q
  simp only [F,D,Q,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,derivative_0_0, derivative_0_1, derivative_0_2, derivative_0_3, derivative_0_4, derivative_0_5, derivative_1_0, derivative_1_1, derivative_1_2, derivative_1_3, derivative_1_4, derivative_1_5, derivative_2_0, derivative_2_1, derivative_2_2, derivative_2_3, derivative_2_4, derivative_2_5, derivative_3_0, derivative_3_1, derivative_3_2, derivative_3_3, derivative_3_4, derivative_3_5, derivative_4_0, derivative_4_1, derivative_4_2, derivative_4_3, derivative_4_4, derivative_4_5, derivative_5_0, derivative_5_1, derivative_5_2, derivative_5_3, derivative_5_4, derivative_5_5]
  apply MvPolynomial.funext; intro x
  norm_num <;> ring

end HodgeComputer.Instances.S_dense_d4_n126
