# Wave 2 — `polarized-k3-remainder` joint search (remote 128-core box, 2026-11-01)

**STATUS: NO WITNESS FOUND.** (Report regenerated 2026-11-01 19:54 from 88 log records; search jobs listed at the end.)


Lane: q = 1, A = X0X1X2Y0Y1Y2, F = ½(f(U)+f(V)), U = X+Y, V = X−Y. Diagonal symmetry group = {1, τ} only, so the
relevant "character sectors" are the two Y-parities of D (Y-even = swap-symmetric in (U,V)). Wave-1 facts used:
[A] has maximal rank on (R/J)_d in every degree (injective d ≤ 9 ⇒ every E1-admissible D of degree ≤ 9 is
Jacobian-type D = Σ a_i F_i; surjective d ≥ 10), degrees ≤ 7 excluded (thread #115), sum principle ⇒ only D with an
irreducible factor of degree ≥ 8 can serve the target.

**Method (new in wave 2, `/tmp/hodge/results/polarized-k3-remainder/w2lib.py`).** Reduced Koszul E2 system of book
§5.2: for an E1-admissible D with exact lift N0 (dF(N0) = A·D; for Jacobian-type D = Σ a_i F_i take N0_i = A a_i),
E2 ⟺ B_D(N0) + T_D(c) ∈ F·S_{2d−6} with c antisymmetric of degree d−4 (the slack G in dF(N) = AD + FG drops out of
E2, so WLOG G = 0). The Q-block is eliminated exactly by the X0-remainder map (F is monic in X0), the c-block is a
dense least-squares problem, and the denominator coordinates a are optimised by Levenberg–Marquardt on the unit
sphere with the Kaufman (variable-projection) Jacobian; F-multiples are projected out of the a-space (lossless by the
prescribed-support theorem). Validated exactly at d = 6 (identity E = F·B_D(N) and the assembled system, mod p).
*Hessian relaxation* = the same system modulo J(F): A·Hess_F(a,a) ∈ J(F)_{2d} (necessary), solved numerically with
deflation against the known spurious family D = g·h (g = f(U)−f(V) ∈ J(F)_6, h a quadric), and exactly (Singular,
mod p) as a system of quadrics with saturation by the spurious component.


## 1. E1-admissible denominator spaces (engine `mult_kernel`, p = 2^31−1), both parities

| d | parity | dim R_d | dim J(F)_d | dim (R/J)_d | admissible dim {D: AD ∈ J} | kernel of [A] mod J | seconds |
|---|---|---|---|---|---|---|---|
| 4 | even | 66 | 0 | 66 | 0 | 0 | 0.7 |
| 4 | odd | 60 | 0 | 60 | 0 | 0 | 0.6 |
| 5 | even | 126 | 3 | 123 | 3 | 0 | 2.2 |
| 5 | odd | 126 | 3 | 123 | 3 | 0 | 2.2 |
| 6 | even | 236 | 18 | 218 | 18 | 0 | 7.2 |
| 6 | odd | 226 | 18 | 208 | 18 | 0 | 7.2 |
| 7 | even | 396 | 63 | 333 | 63 | 0 | 21.9 |
| 7 | odd | 396 | 63 | 333 | 63 | 0 | 22.1 |
| 8 | even | 651 | 168 | 483 | 168 | 0 | 94.5 |
| 8 | odd | 636 | 168 | 468 | 168 | 0 | 107.4 |
| 9 | even | 1001 | 378 | 623 | 378 | 0 | 248.2 |
| 9 | odd | 1001 | 378 | 623 | 378 | 0 | 264.2 |
| 10 | even | 1512 | 750 | 762 | 1044 | 294 | 703.6 |

For d ≤ 9 the admissible space is exactly J(F)_d (Jacobian type); the extra kernel appears at d = 10.

## 3. Exact Hessian relaxation (quadrics mod p, Singular std + saturation by the spurious component)

- generated d=6 pD=0 p=32003: 17 variables (F-multiples removed), 866 quadrics, spurious component dim 0 cut by 17 linear forms; spurious points satisfy all quadrics: max residue []
- generated d=6 pD=0 p=32003: 17 variables (F-multiples removed), 866 quadrics, spurious component dim 0 cut by 17 linear forms; spurious points satisfy all quadrics: max residue []

## 7. Structured linear families vs E1 (exact mod p: dim W ∩ J(F)_d; for d ≤ 9 E1-admissible ⟺ in J(F)_d)

A zero means that no member of the family can serve the target (rigorous: the mod-p rank bounds the rational rank from below, so dim_Q(W ∩ J) ≤ dim_p(W ∩ J) = 0).

| d | family | parity | dim W | dim W ∩ J(F)_d |
|---|---|---|---|---|
| 8 | (U,V)-bidegree (0,8)+(8,0) swap-symmetric | 0 | 45 | **30** |
| 8 | (U,V)-bidegree (0,8)-(8,0) swap-antisymmetric | 1 | 45 | **30** |
| 8 | (U,V)-bidegree (1,7)+(7,1) swap-symmetric | 0 | 108 | **54** |
| 8 | (U,V)-bidegree (1,7)-(7,1) swap-antisymmetric | 1 | 108 | **54** |
| 8 | (U,V)-bidegree (2,6)+(6,2) swap-symmetric | 0 | 168 | **54** |
| 8 | (U,V)-bidegree (2,6)-(6,2) swap-antisymmetric | 1 | 168 | **54** |
| 8 | (U,V)-bidegree (3,5)+(5,3) swap-symmetric | 0 | 210 | **30** |
| 8 | (U,V)-bidegree (3,5)-(5,3) swap-antisymmetric | 1 | 210 | **30** |
| 8 | (U,V)-bidegree (4,4)+(4,4) swap-symmetric | 0 | 120 | **0** |
| 8 | (U,V)-bidegree (4,4)-(4,4) swap-antisymmetric | 1 | 105 | **0** |
| 8 | (X,Y)-bidegree (0,8) | 0 | 45 | **0** |
| 8 | (X,Y)-bidegree (1,7) | 1 | 108 | **0** |
| 8 | (X,Y)-bidegree (2,6) | 0 | 168 | **0** |
| 8 | (X,Y)-bidegree (3,5) | 1 | 210 | **0** |
| 8 | (X,Y)-bidegree (4,4) | 0 | 225 | **0** |
| 8 | (X,Y)-bidegree (5,3) | 1 | 210 | **0** |
| 8 | (X,Y)-bidegree (6,2) | 0 | 168 | **0** |
| 8 | (X,Y)-bidegree (7,1) | 1 | 108 | **0** |
| 8 | (X,Y)-bidegree (8,0) | 0 | 45 | **0** |
| 8 | UV cone vertex line: D in C[U0,U1,V0,V1] swap-antisymmetric | 1 | 80 | **0** |
| 8 | UV cone vertex line: D in C[U0,U1,V0,V1] swap-symmetric | 0 | 85 | **0** |
| 8 | UV cone: h(U) + h(V), h ternary form (Y-even) | 0 | 45 | **30** |
| 8 | UV cone: h(U) - h(V), h ternary form (Y-odd) | 1 | 45 | **30** |
| 8 | cone vertex line: D in C[X0,X1,X2,Y0] (parity 0) | 0 | 95 | **0** |
| 8 | cone vertex line: D in C[X0,X1,X2,Y0] (parity 1) | 1 | 70 | **0** |
| 8 | cone vertex line: D in C[X0,X1,Y0,Y1] (parity 0) | 0 | 85 | **0** |
| 8 | cone vertex line: D in C[X0,X1,Y0,Y1] (parity 1) | 1 | 80 | **0** |
| 8 | cone vertex line: D in C[X0,X1,Y0,Y2] (parity 0) | 0 | 85 | **0** |
| 8 | cone vertex line: D in C[X0,X1,Y0,Y2] (parity 1) | 1 | 80 | **0** |
| 8 | cone vertex line: D in C[X0,X2,Y0,Y2] (parity 0) | 0 | 85 | **0** |
| 8 | cone vertex line: D in C[X0,X2,Y0,Y2] (parity 1) | 1 | 80 | **0** |
| 8 | cone vertex line: D in C[X0,Y0,Y1,Y2] (parity 0) | 0 | 95 | **0** |
| 8 | cone vertex line: D in C[X0,Y0,Y1,Y2] (parity 1) | 1 | 70 | **0** |
| 8 | cone vertex plane {X=0}: D in C[X0,X1,X2] | 0 | 45 | **0** |
| 8 | cone vertex plane {Y=0}: D in C[Y0,Y1,Y2] (Y-even part) | 0 | 45 | **0** |
| 8 | cone vertex plane {Y=0}: D in C[Y0,Y1,Y2] (Y-odd part) | 1 | 0 | **0** |
| 9 | (U,V)-bidegree (0,9)+(9,0) swap-symmetric | 0 | 55 | **45** |
| 9 | (U,V)-bidegree (0,9)-(9,0) swap-antisymmetric | 1 | 55 | **45** |
| 9 | (U,V)-bidegree (1,8)+(8,1) swap-symmetric | 0 | 135 | **90** |
| 9 | (U,V)-bidegree (1,8)-(8,1) swap-antisymmetric | 1 | 135 | **90** |
| 9 | (U,V)-bidegree (2,7)+(7,2) swap-symmetric | 0 | 216 | **108** |
| 9 | (U,V)-bidegree (2,7)-(7,2) swap-antisymmetric | 1 | 216 | **108** |
| 9 | (U,V)-bidegree (3,6)+(6,3) swap-symmetric | 0 | 280 | **90** |
| 9 | (U,V)-bidegree (3,6)-(6,3) swap-antisymmetric | 1 | 280 | **90** |
| 9 | (U,V)-bidegree (4,5)+(5,4) swap-symmetric | 0 | 315 | **45** |
| 9 | (U,V)-bidegree (4,5)-(5,4) swap-antisymmetric | 1 | 315 | **45** |
| 9 | (X,Y)-bidegree (0,9) | 1 | 55 | **0** |
| 9 | (X,Y)-bidegree (1,8) | 0 | 135 | **0** |
| 9 | (X,Y)-bidegree (2,7) | 1 | 216 | **0** |
| 9 | (X,Y)-bidegree (3,6) | 0 | 280 | **0** |
| 9 | (X,Y)-bidegree (4,5) | 1 | 315 | **0** |
| 9 | (X,Y)-bidegree (5,4) | 0 | 315 | **0** |
| 9 | (X,Y)-bidegree (6,3) | 1 | 280 | **0** |
| 9 | (X,Y)-bidegree (7,2) | 0 | 216 | **0** |
| 9 | (X,Y)-bidegree (8,1) | 1 | 135 | **0** |
| 9 | (X,Y)-bidegree (9,0) | 0 | 55 | **0** |
| 9 | UV cone vertex line: D in C[U0,U1,V0,V1] swap-antisymmetric | 1 | 110 | **0** |
| 9 | UV cone vertex line: D in C[U0,U1,V0,V1] swap-symmetric | 0 | 110 | **0** |
| 9 | UV cone: h(U) + h(V), h ternary form (Y-even) | 0 | 55 | **45** |
| 9 | UV cone: h(U) - h(V), h ternary form (Y-odd) | 1 | 55 | **45** |
| 9 | cone vertex line: D in C[X0,X1,X2,Y0] (parity 0) | 0 | 125 | **0** |
| 9 | cone vertex line: D in C[X0,X1,X2,Y0] (parity 1) | 1 | 95 | **0** |
| 9 | cone vertex line: D in C[X0,X1,Y0,Y1] (parity 0) | 0 | 110 | **0** |
| 9 | cone vertex line: D in C[X0,X1,Y0,Y1] (parity 1) | 1 | 110 | **0** |
| 9 | cone vertex line: D in C[X0,X1,Y0,Y2] (parity 0) | 0 | 110 | **0** |
| 9 | cone vertex line: D in C[X0,X1,Y0,Y2] (parity 1) | 1 | 110 | **0** |
| 9 | cone vertex line: D in C[X0,X2,Y0,Y2] (parity 0) | 0 | 110 | **0** |
| 9 | cone vertex line: D in C[X0,X2,Y0,Y2] (parity 1) | 1 | 110 | **0** |
| 9 | cone vertex line: D in C[X0,Y0,Y1,Y2] (parity 0) | 0 | 95 | **0** |
| 9 | cone vertex line: D in C[X0,Y0,Y1,Y2] (parity 1) | 1 | 125 | **0** |
| 9 | cone vertex plane {X=0}: D in C[X0,X1,X2] | 0 | 55 | **0** |
| 9 | cone vertex plane {Y=0}: D in C[Y0,Y1,Y2] (Y-even part) | 0 | 0 | **0** |
| 9 | cone vertex plane {Y=0}: D in C[Y0,Y1,Y2] (Y-odd part) | 1 | 55 | **0** |

The nonzero intersections are exactly the (U,V)-bidegree pieces of J(F)_d (J(F) is bigraded in the (U,V) frame since F = ½(f(U)+f(V))): at d = 8 the swap-symmetric pieces have dimensions 30, 54, 54, 30 for (U,V)-bidegrees (8,0), (7,1), (6,2), (5,3) (plus mirrors), summing to 168 = dim J(F)_8^{even}; pure (X,Y)-bidegree pieces and all coordinate cones (vertex a plane or a line) contain no E1-admissible denominator at d = 8, 9.

## 8. Structural observations made while setting up the search (exact, (U,V)-frame)

* **J(F) is bigraded in (U,V)**, so the Jacobian-type family at d ≤ 9 splits into swap-symmetric (U,V)-bidegree pieces
  P_a = {Σ_i f_i(U) c_i + swap : c_i of bidegree (a−5, d−a)}: dims 30/54/54/30 at d = 8 (a = 8,7,6,5), 45/90/108/90/45 at
  d = 9. No F-multiple and no spurious g·h lies inside a single piece (both have ≥ 4 bidegree components), so the pieces
  are clean sub-searches; but E2 and the relaxation mix bidegrees (A has (U,V)-bidegrees (6,0),(4,2),(2,4),(0,6)).
* **On the pure-block piece P_8 = {h(U) + h(V), h ∈ J(f)_8} the Hessian relaxation is just 6 quadrics.**  With
  cofactors c_i(U) (h = Σ c_i f_i, c_i cubic), Hess_F(v,v) ∝ H(U) + H(V) with H = Σ_ij c_i c_j f_ij ∈ C[z]_10.  In the
  bigraded ring R_f ⊗ R_f the components of A·(H(U)+H(V)) have (U,V)-bidegrees (16−2s, 2s) (and mirrors), s = |S| for the
  8 terms U_{S^c}²V_S² of 64·A: s = 0, 1 vanish identically (R_{f,k} = 0 for k > 12), s = 2 gives [U_i²H] = 0 ∈ R_{f,12}
  and s = 3 gives [H]⊗[z0²z1²z2²] = 0, i.e. H ∈ J(f)_10 (dim R_{f,10} = 6), which implies the s = 2 conditions.  So the
  relaxation restricted to P_8 is H ∈ J(f)_{10}: 6 quadrics in 30 unknowns, a solution variety of projective dimension
  ≥ 23 consisting of denominators that are *not* of the spurious form g·h.  The degree-8 relaxation therefore cannot be
  decisive on its own, and the full reduced-E2 system must be searched; the P_8 full-E2 runs (30 parameters) do this
  directly.  At d = 9 on P_9 (c_i quartic, H of degree 12) only the socle component survives: a single quadric.
* **Optimizer validation.**  On a synthetic consistent instance (the d = 6 reduced system shifted so that a chosen a* is
  an exact solution) the LM/variable-projection search converged from 3/3 random starts to f ≈ 1e−26 and
  |u − u*| ≈ 1e−16 in 13–18 iterations (quadratic tail).  So a residual plateau at 1e−3…1e−1 in the real search is
  evidence against solutions in the searched family, not an optimizer artefact.


## Jobs / reproduction

Remote job logs are mirrored in `/tmp/hodge/remote/_jobs/`; scripts in `/tmp/hodge/results/polarized-k3-remainder/w2_*.py`; queued/launched commands in `w2_jobs.txt`. Each job appends to its own `w2-*.jsonl` (merged by `w2_merge.py`).
