import Cap.K_m3k100
open HodgeComputer.Kron HodgeComputer.Cap.K_m3k100
-- cost of just computing psi of one big tree (N4, 288 leaves) vs its derivative tree
example : psi 17 27 tN4 ≠ 0 := by decide +kernel
