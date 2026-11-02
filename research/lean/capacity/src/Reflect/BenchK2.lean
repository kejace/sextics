import Mathlib.Data.Nat.Basic
-- shifts: 5M-bit numbers, product of two shifted odd numbers
def c : Nat := Nat.shiftLeft 12345678901234567 5000000 + 1
def d : Nat := Nat.shiftLeft 98765432109876543 4000000 + 3
example : (c * d) % 1000003 = ((c % 1000003) * (d % 1000003)) % 1000003 := by decide +kernel
example : Nat.shiftLeft (Nat.shiftLeft 1 2000000) 3000000 = Nat.shiftLeft 1 5000000 := by decide +kernel
