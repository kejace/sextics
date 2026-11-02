import Cap.BigTree
open HodgeComputer.Kron
set_option maxRecDepth 4000
example : homB 3 tBig = false := by decide +kernel
