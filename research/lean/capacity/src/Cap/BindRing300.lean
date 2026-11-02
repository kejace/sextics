import Cap.Bind300
set_option maxHeartbeats 8000000
set_option maxRecDepth 8000
set_option linter.all false
set_option profiler true
set_option profiler.threshold 500
open HodgeComputer.Kron HodgeComputer.Cap.Bind300
-- accepted-template style binding: evaluate both sides and normalize with ring
theorem b_ring : HW_D = toC tD := by
  apply MvPolynomial.funext; intro x
  norm_num [HW_D, tD, toC] <;> ring
