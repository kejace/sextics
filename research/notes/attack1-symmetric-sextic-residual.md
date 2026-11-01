# Attack 1 — `symmetric-sextic-residual` (q = 1, A = X0 X1 X2 Y0 Y1 Y2)

**Status: NO WITNESS FOUND.** (Updated continuously; the session that started this note crashed and was
resumed on 2026-11-01 evening. All results below were recomputed or re-read from `log.jsonl`.)

Work directory `/tmp/hodge/results/symmetric-sextic-residual/` — log `log.jsonl`; scripts `h.py`
(helpers/logging), `invsolve.py` (G-invariant reduced solver), `charsolve.py` (character-twisted reduced
solver, this session), `e1_kernel.py`, `e1_kernel_inv_sweep.py`, `families1.py`, `test_d7.py`, `test_K.py`,
`served_scan*.py`, `mult.py`, `check_served*.py`.

## 1. The obstruction and what a certificate must do

* `F = p6(X) + p6(Y) + (5/7)·Q_X·Q_Y·(Q_X + 2 Q_Y)`, `Q = e_2` (elementary symmetric quadric), so `F` is
  `G = S3(X) × S3(Y)`-invariant; the full projective linear symmetry group is `G × <tau>` (tau: Y -> -Y).
  There is no X<->Y swap (coefficients 5/7 vs 10/7).
* `A = e_3(X) e_3(Y)` is a `G`-invariant polynomial; since `Omega` changes sign under every transposition,
  the residue class `alpha_A in H^{3,1}_prim(X)` transforms by `sgn ⊗ sgn` (and is tau-invariant).
* Sector size: the `G`-invariant Y-odd sextics modulo `J(F)` form a **17-dimensional** space
  (`mult.py`, equivariant Hilbert function; confirmed by the reduced solver: 19 orbit representatives, 2 of them in `J`).
  `[A]` is one line in this 17-dimensional `sgn⊗sgn`-residue sector. `A ∉ J(F)` (D = 1 fails E1).
* Reflections of `G × <tau>` are exactly the six transpositions. By the book's §3.5 argument,
  `U = Σ_s H^4(X)^s` is divisor-supported and `G`-stable, and the `sgn⊗sgn`-isotypic part meets `U` trivially
  (projecting a reflection-invariant vector with the `sgn⊗sgn` projector gives 0). So **no reflection quotient,
  no symmetric-function Piola map (`delta = Δ_X Δ_Y`, weights (1,2,3,1,2,3), serves only `Δ_XΔ_Y`), and no
  sum of such constructions can reach any of the 17 classes**. Every element of the index-2 subgroup
  `G' = {sgn σ = sgn ρ}` (under which `alpha_A` is invariant) has `det = 1` and age ≥ 1 at every fixed point
  (3-cycles: ζ, ζ²; double transpositions: −1, −1), so `X/G'` is a (singular) Calabi–Yau complete intersection
  `{G = 0, w² = disc_X·disc_Y} ⊂ P(1,2,3,1,2,3,6)` and Kollár–Larsen gives nothing either.
* Gysin mechanism: a supporting divisor `Z = X ∩ V(D)` (WLOG `G`-stable, replacing `D` by its orbit product)
  must have a resolution with `H^{2,0}(Z~)^{sgn⊗sgn} ≠ 0`, i.e. a holomorphic 2-form anti-invariant under
  every transposition, and its Gysin image must hit `[A]` among the 17 classes.
  Natural candidates: (a) cone with vertex the plane {X=0} over a symmetric plane curve `C ⊂ P²_X`
  (exceptional divisor `C × E_Y`, `E_Y = {p6(Y)=0}` the Fermat sextic curve of genus 10;
  `H^{2,0}(C×E_Y)^{sgn⊗sgn} = H^{1,0}(C)^{sgn} ⊗ H^{1,0}(E_Y)^{sgn}` = (invariant forms of degree deg C − 3) ⊗ 3);
  (b) the mirror with X<->Y; (c) cone with vertex the line `L = {(x,x,x,y,y,y)}` over a `G`-invariant surface in
  `P(std ⊕ std)` (6 vertex points `L ∩ X`, each exceptional divisor ≅ the surface).
* Hard necessary condition (prescribed-support theorem, q = 1): support on `V(D)` ⇒ a certificate with
  denominator `D_red` ⇒ `A·D_red ∈ J(F)`. The E1 kernel `Ann_R([A])_d` is **zero for d ≤ 6** (all isotypic
  types, mod p = 2147483647), so every supporting divisor of degree ≤ 6 has equation in `J(F)`; this kills
  (a), (b) for curves of degree ≤ 6 and (c) for surfaces of degree ≤ 6 outright, in particular all cubic cones
  (Hesse pencil) and all quartic cones over `L`, **regardless of their 2-forms** — their Gysin images land in
  the other 16 classes (or vanish).
* Sum of supports: for `D = D_1 D_2` without common factor, served(D) = served(D_1) + served(D_2)
  (resolution of a reducible divisor = disjoint union). So the right global experiment is the span of the
  served invariant-target spaces over many candidate `D` (§2.4), not only single-`D` tests.

## 2. Experiments

(Tables filled in below as they complete; all mod p = 2147483647 unless stated.)

### 2.1 Recovered from the crashed session (engine v2 / `invsolve.py`, all in `log.jsonl`)

| family | degrees | #D | result |
|---|---|---|---|
| linear forms (X0, e1X, X0±X1, X0+X1−2X2, Y0, e1Y, Y0−Y1) | 1 | 8 | E1 fails (rank 63/63, aug 64) |
| quadrics: QX, QY, p2X, p2Y, e1X², pencils QX+cQY (c ∈ ±1,±2,±1/2,±3), p2X±p2Y, X0²±Y0², X0X1, X0X1−Y0Y1, X0Y0, e1Xe1Y, (X0−X1)(Y0−Y1), p2−Q | 2 | 24 | E1 fails |
| cubics: Hesse pencils p3X+c e3X, p3Y+c e3Y (c ∈ 0,±1,±2,±3,±6,±3/2), e3, e1Q, Δ_X, Δ_Y, mixed e1X·QY etc., X0Y0Y1, X0³+X1³+Y0²X2 | 3 | 36 | E1 fails |
| quartics: p4, p4X±p4Y, QXQY, p2Xp2Y, QX²±QY², Δ·e1 (4 variants), e3·e1 (4), p3·e1 (4), Fermat-type sums of 4 sixth... fourth powers (3), (X0²−Y0²)(X1²−Y1²), QX²+QXQY+QY², (QX−QY)(QX+2QY), QX(QX+2QY), QY(QX+2QY), e1X⁴+e1Y⁴ | 4 | 32 | E1 fails (nullity 6–9 = J_4 part only) |
| E1-kernel sweep (full, all types) `Ann_R([A])_d` | 1–8 | — | dim 0,0,0,0,0,0 for d=1..6; **d=7: 1 per parity** (even: `Δ_X·e2(X0²,X1²,X2²)`, odd: Y-mirror); d=8: 6 per parity |
| G-invariant E1-kernel sweep (`e1inv_d*_p*.pkl`) | 1–16 | — | 0 up to d=7; (8,ev)=2, (9,·)=1, (10,ev)=28, (10,odd)=4, (11,·)=28, (12,ev)=50, (12,odd)=26, (13,·)=43, (14,ev)=54, (14,odd)=34, (15,·)=41, (16,ev)=43, (16,odd)=28 |
| degree-7 kernel elements, full fixed-D system (compressed rref) | 7 | 2 | **inconsistent** (rank 3991, nullity 506, 25935×4497; 343 s / 593 s) |
| invariant degree-8 even kernel: 2 reps, 4 random kernel combos, 6 random points of (J:A)_8^G, 3 random points of J_8^G | 8 | 15 | all inconsistent (reduced system 1858×421, rank 362) |
| invariant degree-9 even kernel: rep, 4 kernel combos, 5 random points of (J:A)_9^G | 9 | 10 | all inconsistent (2648×673, rank 566) |
| radial derivatives EX(F), EY(F), EX−EY (∈ J_6) | 6 | 3 | inconsistent (838×161, rank 140) |

### 2.2 Served-span scan (this session)

`served_scan.py` (crashed before output) was rerun. **Bug found**: its served-target vectors are in
`Sym(m)`-coordinates while `J`/`A` vectors are in `orbit_sum` coordinates (factor 36/|orbit| per monomial
orbit), which made the output read "2 classes beyond J served by every D". Corrected check
(`check_served4.py`): the 2 invariant classes served by `QX` and by `p3X+e3X` are **in `J(F)`**. So all
invariant D of degree 2–6 in the list (≈85 denominators) serve **nothing** in the 17-dimensional sector.
