import Mathlib.Data.List.Basic

/-! Kernel-evaluable sparse polynomials: packed monomial keys (base 256 digits = exponents) and ℤ coefficients. -/
namespace HodgeComputer.Reflect

abbrev Term := Nat × Int
abbrev SP := List Term

def B : Nat := 256

def dig : Nat → Nat → List Nat
  | 0, _ => []
  | n+1, k => k % B :: dig n (k / B)

def expo (k i : Nat) : Nat := (dig 6 k).getD i 0

def digitsOk (p : SP) : Bool := p.all (fun t => (dig 6 t.1).all (fun d => d < 128))

def scaleT (t : Term) (q : SP) : SP := q.map (fun u => (t.1 + u.1, t.2 * u.2))

def addT : Nat → SP → SP → SP
  | 0, _, _ => []
  | _+1, [], q => q
  | _+1, p, [] => p
  | n+1, (k, c) :: p, (k', c') :: q =>
    if k' < k then (k, c) :: addT n p ((k', c') :: q)
    else if k < k' then (k', c') :: addT n ((k, c) :: p) q
    else if c + c' = 0 then addT n p q else (k, c + c') :: addT n p q

def add (p q : SP) : SP := addT (p.length + q.length) p q

def mulT : Nat → SP → SP → SP
  | 0, _, _ => []
  | _+1, [], _ => []
  | _+1, [t], q => scaleT t q
  | n+1, p, q => add (mulT n (p.take (p.length / 2)) q) (mulT n (p.drop (p.length / 2)) q)

def mul (p q : SP) : SP := mulT (p.length + 1) p q

def smul (a : Int) (p : SP) : SP := if a = 0 then [] else p.map (fun t => (t.1, a * t.2))

def pderivT (i : Nat) (p : SP) : SP :=
  p.filterMap (fun t => if expo t.1 i = 0 then none else some (t.1 - B ^ i, t.2 * expo t.1 i))

/-- pseudo-random test polynomials: n terms with keys from a simple LCG -/
def mkPoly (seed n : Nat) : SP :=
  let rec go : Nat → Nat → SP
    | 0, _ => []
    | m+1, s => let s' := (s * 1103515245 + 12345) % 2147483648
                let e0 := s' % 7; let e1 := (s' / 7) % 7; let e2 := (s'/49) % 7; let e3 := (s'/343) % 7; let e4 := (s'/2401) % 7; let e5 := (s'/16807) % 7
                (e0 + 256*e1 + 65536*e2 + 16777216*e3 + 4294967296*e4 + 1099511627776*e5, ((s' / 117649) % 19 : Nat) - 9) :: go m s'
  -- normalize via add with []
  let raw := go n seed
  raw.foldr (fun t acc => add [t] acc) []

end HodgeComputer.Reflect
