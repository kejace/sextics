import Cap.BigSorted
open HodgeComputer.Kron
set_option maxRecDepth 4000
example : psiEq 21 42 (.poly tSorted) (.poly tSorted) = true := by decide +kernel
