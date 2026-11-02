import Mathlib.Data.Nat.Basic
def a : Nat := Nat.pow 3 1000000 + 7
def b : Nat := Nat.pow 5 900000 + 11
example : (a * b) % 1000003 = ((a % 1000003) * (b % 1000003)) % 1000003 := by decide +kernel
