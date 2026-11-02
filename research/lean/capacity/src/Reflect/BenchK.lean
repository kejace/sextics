import Mathlib.Data.Nat.Basic
-- big-number kernel ops: shift, pow, mul on ~2M-bit numbers
def a : Nat := (3 ^ 1000000) + 7
def b : Nat := (5 ^ 900000) + 11
example : (a * b) % 1000003 = ((a % 1000003) * (b % 1000003)) % 1000003 := by decide +kernel
example : (1 <<< 2000000) * (1 <<< 3000000) = 1 <<< 5000000 := by decide +kernel
example : (2 ^ 20000000) = 1 <<< 20000000 := by decide +kernel
