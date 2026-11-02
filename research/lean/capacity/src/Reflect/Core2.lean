import Mathlib.Data.List.Basic
/-! Variant: explicit Nat.rec recursors (no brecOn tables). -/
namespace HodgeComputer.Reflect2
abbrev Term := Nat × Int
abbrev SP := List Term

def scaleT (t : Term) (q : SP) : SP := q.map (fun u => (t.1 + u.1, t.2 * u.2))

def mergeStep (ih : SP → SP → SP) (p q : SP) : SP :=
  match p, q with
  | [], q => q
  | p, [] => p
  | (k, c) :: p', (k', c') :: q' =>
    match Nat.blt k' k with
    | true => (k, c) :: ih p' ((k', c') :: q')
    | false =>
      match Nat.blt k k' with
      | true => (k', c') :: ih ((k, c) :: p') q'
      | false =>
        match c + c' with
        | Int.ofNat 0 => ih p' q'
        | s => (k, s) :: ih p' q'

def addT (n : Nat) : SP → SP → SP :=
  Nat.rec (motive := fun _ => SP → SP → SP) (fun _ _ => []) (fun _ ih => mergeStep ih) n

def add (p q : SP) : SP := addT (p.length + q.length) p q

def mulStep (ih : SP → SP → SP) (p q : SP) : SP :=
  match p with
  | [] => []
  | [t] => scaleT t q
  | p => add (ih (p.take (p.length / 2)) q) (ih (p.drop (p.length / 2)) q)

def mulT (n : Nat) : SP → SP → SP :=
  Nat.rec (motive := fun _ => SP → SP → SP) (fun _ _ => []) (fun _ ih => mulStep ih) n

def mul (p q : SP) : SP := mulT (p.length + 1) p q

def mkPoly (seed n : Nat) : SP :=
  let rec go : Nat → Nat → SP
    | 0, _ => []
    | m+1, s => let s' := (s * 1103515245 + 12345) % 2147483648
                let e0 := s' % 7; let e1 := (s' / 7) % 7; let e2 := (s'/49) % 7; let e3 := (s'/343) % 7; let e4 := (s'/2401) % 7; let e5 := (s'/16807) % 7
                (e0 + 256*e1 + 65536*e2 + 16777216*e3 + 4294967296*e4 + 1099511627776*e5, ((s' / 117649) % 19 : Nat) - 9) :: go m s'
  let raw := go n seed
  raw.foldr (fun t acc => add [t] acc) []
end HodgeComputer.Reflect2
