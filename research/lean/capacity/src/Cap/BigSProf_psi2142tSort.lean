import Cap.BigSorted
open HodgeComputer.Kron
set_option maxRecDepth 4000
example : psi 21 42 tSorted ≠ 0 := by decide +kernel
