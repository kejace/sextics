import Cap.BigTree
open HodgeComputer.Kron
set_option maxRecDepth 4000
example : psi 21 8 tBig ≠ 0 := by decide +kernel
