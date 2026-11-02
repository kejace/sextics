import Mathlib.Data.Nat.Basic
-- pow with exponent above 2^24: is it still accelerated?
example : Nat.pow 2 16777217 = Nat.shiftLeft 1 16777217 := by decide +kernel
