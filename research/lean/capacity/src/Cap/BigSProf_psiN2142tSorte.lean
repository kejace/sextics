import Cap.BigSorted
open HodgeComputer.Kron
set_option maxRecDepth 4000
example : (psiN 21 42 tSorted).1 ≠ 0 := by decide +kernel
