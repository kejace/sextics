# Wave 2 — `polarized-k3-remainder` joint search (remote 128-core box, 2026-11-01)

**STATUS: NO WITNESS FOUND.** (Report regenerated 2026-11-01 20:57 from 212 log records; search jobs listed at the end.)


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


## 0. Summary so far

* **No witness.** Every exact decision is negative; every numerical landscape is bounded away from zero except on the
  predicted relaxation component (which the exact E2 rejects).
* **Exact (mod p = 2^31−1 unless stated), this wave:**
  - E1-admissible spaces per degree and parity (table 1): at d ≤ 9 exactly J(F)_d; kernel appears at d = 10 (294/261).
  - All coordinate cones (vertex a plane or a line in the (X,Y) or (U,V) frames) and all pure (X,Y)-bidegree octics/nonics
    contain **no** E1-admissible denominator (table 7); the only structured linear families meeting J(F)_d are the
    (U,V)-bidegree pieces of J(F)_d itself.
  - 60 structured candidates (Hessian cofactors/pairings, third polars, Fermat, Plücker, polar families, polar pairings,
    determinantal, Hessian determinants) tested for E1 exactly over Q (table 6b): only the Jacobian-type polar families
    pass; 11 of those (all d = 8 polar octics and three controls) are **E2-inconsistent** by the exact reduced system
    (table 6c); the wave-1 engine result for the control octic is reproduced exactly (nullity 865 = 636 + 229).
  - Hessian relaxation at d = 8 (both sectors share one system): 468 quadrics on 147 parameters; restricted to the
    pure-block piece P_8 it is the 6 conditions H ∈ J(f)_{10}; the other bidegree pieces are over-determined
    (90/234/114 quadrics on 54/54/30 parameters).
* **Numerical:** LM/variable-projection search validated (quadratic convergence to an exact zero on a synthetic consistent
  system); relaxation restarts at d = 8 converge only into the P_8 component; deflated restarts find nothing else.
* **Pending on the saturated remote box:** full reduced-E2 variable-projection search over the 156/159-dimensional Jacobian
  families at d = 8 (both parities) and d = 9, the P_8/P_7/P_6/P_5 piece searches, the non-Jacobian d = 10 family, the
  Singular saturation of the d = 8 relaxation variety, and the remaining structured candidates at d = 9, 10, 11.


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
| 10 | odd | 1491 | 747 | 744 | 1008 | 261 | 743.2 |
| 11 | even | 2184 | 1341 | 843 | 1851 | 510 | 1726.4 |

For d ≤ 9 the admissible space is exactly J(F)_d (Jacobian type); the extra kernel appears at d = 10.

## 2b. Hessian relaxation in the exact (U,V)-tensor frame (`w2tensor.py`): one system for both parity sectors

Cofactor cubics c_i(U,V) (K = 168 at d = 8) modulo the 21 Euler directions c_i → c_i + U_i h (which are exactly the F-multiples and the spurious g·h family of the two sectors combined; the relaxation is invariant under them: verified to 1e−11). Quotient: 147 parameters, 468 independent quadrics (the swap-symmetric part of A·Hess vanishes identically, checked to 1e−13). Theoretical check: on the pure-U piece P_8 the residual span has rank exactly 6 = dim R_{f,10}.

- setup d=8: K=168, dimP=147, coords=468, Euler dirs=21, P8 quotient dim=24, pieces by cofactor V-degree {'0': 30, '1': 54, '2': 54, '3': 30}, build 8.1 s, mean f on the sphere 5.079e+05

| d | deflated vs P8 | restarts | f = 0 hits | inside P8 | off P8 | f quantiles (min/10/50/90/max) |
|---|---|---|---|---|---|---|
| 8 | False | 10 | 8 | 6 | 2 | 1.3e-29 / 1.7e-29 / 1.2e-28 / 1.2e-08 / 1.2e-07 |
| 8 | True | 6 | 0 | 0 | 0 | 4.4e-09 / 4.9e-08 / 6.5e-07 / 5.2e-05 / 7.9e-05 |
| 8 | True | 3 | 0 | 0 | 0 | 5.6e-13 / 8.9e-13 / 2.2e-12 / 5.1e-11 / 6.3e-11 |
| 8 | True | 4 | 0 | 0 | 0 | 1.3e-08 / 4.5e-08 / 4.5e-03 / 4.9e-02 / 6.6e-02 |
| 8 | False | 11 | 8 | 6 | 2 | 6.2e-29 / 8.8e-29 / 5.2e-28 / 3.8e-15 / 7.1e-12 |

Zero-residual points classified off P8 by the 1e−6 distance threshold (all of them have cofactor energy ≥ 1 − 1e−7 in the V-degree-0 piece, non-vanishing gradients and f decreasing along the iterations: they are LM runs still converging into P_8 along the flat (4th-order) V-degree-1 directions, not separate components):
- seed 1 restart 2: f=1.1e-19, grad 4.0e-11, dist to P8 3.2e-04, profile {'0': 1.0, '1': 1.02e-07, '2': 9.18e-14, '3': 6.31e-19}, tangent nullity 16
- seed 1 restart 4: f=2.3e-23, grad 5.2e-11, dist to P8 4.9e-05, profile {'0': 1.0, '1': 2.41e-09, '2': 3.76e-16, '3': 1.58e-22}, tangent nullity 21
- seed 5 restart 3: f=1.4e-26, grad 8.7e-13, dist to P8 2.2e-06, profile {'0': 1.0, '1': 4.63e-12, '2': 6.86e-23, '3': 6.81e-32}, tangent nullity 27
- seed 5 restart 9: f=3.6e-20, grad 1.4e-10, dist to P8 1.6e-04, profile {'0': 1.0, '1': 2.51e-08, '2': 2.57e-14, '3': 1.69e-19}, tangent nullity 15

### 2c. Exact relaxation variety (tensor frame, Singular mod p)

- exported d=8 p=32003: 147 variables (168 cofactor coordinates minus 21 Euler), 468 quadrics; the P8 component lives in 24 of the variables

### 2d. Exact structure of the relaxation by (U,V)-bidegree piece (tensor frame, float rank with 1e−9 cutoff)

d = 8: number of independent relaxation quadrics when the cofactors are restricted to the pieces of V-degree v (P_{8−v}):

| cofactor pieces (V-degree) | parameters | independent quadrics |
|---|---|---|
| (0,) | 30 | 6 |
| (1,) | 54 | 90 |
| (2,) | 54 | 234 |
| (3,) | 30 | 114 |
| (0, 1) | 84 | 126 |
| (0, 2) | 84 | 234 |
| (0, 3) | 60 | 300 |
| (1, 2) | 108 | 414 |
| (1, 3) | 84 | 234 |
| (2, 3) | 84 | 398 |
| (0, 1, 2, 3) | 168 | 468 |

(Exact mod-32003 ranks agree: {'(0,)': 6, '(1,)': 90, '(2,)': 234, '(3,)': 114, '(0, 1)': 126, '(0, 3)': 300, '(0, 1, 2, 3)': 468}.)

Only the pure-block piece P_8 (V-degree 0) is under-determined (6 quadrics = dim R_{f,10}, as predicted); every other piece and every pair of pieces has more quadrics than parameters, so relaxation solutions off P_8 can only come from special loci, not from dimension counting.

d = 9: K=378 cofactor quartics, 56 Euler directions, quotient dim 322, 208 independent quadrics (all independent generically: rank 208); on P_9 the residual span has rank 1 (single quadric, as predicted).  With 208 quadrics on 322 parameters the d = 9 relaxation has a solution variety of dimension ≥ 113 and is useless as a filter: the d = 9 search must be done with the full E2 system.

Tangent analysis at converged P_8 relaxation solutions (d = 8): Jacobian nullity 33, 33, 33; kernel energy by cofactor V-degree {'0': 0.5455, '1': 0.4545, '2': 0.0, '3': 0.0}. The ~28–32 kernel directions leaving P_8 point into the V-degree-1 piece; they are a first-order artefact of the degree truncation of R_f (only ≤ 27 of the 468 quadrics have a nonzero mixed term between P_8 and P_7 because R_{f,k} = 0 for k > 12), and are obstructed at second order (the P_7 self-terms meet ≥ 90 independent quadrics).  Deflated restarts (f/dist(P_8)²) found no zero off P_8.

## 3. Exact Hessian relaxation (quadrics mod p, Singular std + saturation by the spurious component)

- generated d=6 pD=0 p=32003: 17 variables (F-multiples removed), 866 quadrics, spurious component dim 0 cut by 17 linear forms; spurious points satisfy all quadrics: max residue []
- generated d=6 pD=0 p=32003: 17 variables (F-multiples removed), 866 quadrics, spurious component dim 0 cut by 17 linear forms; spurious points satisfy all quadrics: max residue []

## 4. Full reduced-E2 variable-projection search (Jacobian-type degrees), LM on the sphere

- setup d=8 pD=0 seed 7: K=30 (dimP=30 after F-multiples), rows 10152 → 8661 after X0-remainder, c-unknowns 954 (eff 954), Q 1491 eliminated; build 1.5 s, one evaluation 13.04 s; probe rel 5.865e-03

| d | pD | restarts | best rel | rel quantiles (10/50/90%) | median iters | median s/restart | dist-to-spurious at best | #rel<1e-6 | #rel<1e-10 |
|---|---|---|---|---|---|---|---|---|---|

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


## 6b. Structured candidates: exact E1 over Q (tensor frame: NF of A·D in R_f ⊗ R_f), all 60 candidates in 12 s

| family | label | d | parity | E1 (exact, Q) |
|---|---|---|---|---|
| control-product | Q2(U;V) (U.V) + swap | 8 | 0 | fail (678 NF terms) |
| control-product | g * (U.V)  [g = f(U)-f(V)] | 8 | 1 | PASS (Jacobian-type) |
| control-spurious | g*(U0^2-V0^2) = 4 g X0Y0 (F-reduced version of the wave-1 control) | 8 | 0 | PASS (Jacobian-type) |
| control-spurious | g*(U0V0+U1V1) (Y-odd sector spurious point) | 8 | 1 | PASS (Jacobian-type) |
| control-wave1 | fU*U0^2 + fV*V0^2 (wave-1: inconsistent rank 6632 nullity 865) | 8 | 0 | PASS (Jacobian-type) |
| fermat | sum U^8 + sum V^8 | 8 | 0 | fail (72 NF terms) |
| fermat | sum U^8 - sum V^8 | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | cof00(U)+cof00(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | cof00(U)-cof00(V) | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | cof11(U)+cof11(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | cof11(U)-cof11(V) | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | cof22(U)+cof22(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | cof22(U)-cof22(V) | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | cof01(U)+cof01(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | cof01(U)-cof01(V) | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | cof02(U)+cof02(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | cof02(U)-cof02(V) | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | cof12(U)+cof12(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | cof12(U)-cof12(V) | 8 | 1 | fail (72 NF terms) |
| hess-cofactor | trace-cof(U)+(V) | 8 | 0 | fail (72 NF terms) |
| hess-cofactor | trace-cof(U)-(V) | 8 | 1 | fail (72 NF terms) |
| hess-pairing | sum_ij f_ij(U) f_ij(V) | 8 | 0 | fail (690 NF terms) |
| hess-pairing | lap(U) lap(V) | 8 | 0 | fail (624 NF terms) |
| hess-pairing | lap(U)^2 + lap(V)^2 | 8 | 0 | fail (72 NF terms) |
| hess-pairing | lap(U)^2 - lap(V)^2 | 8 | 1 | fail (72 NF terms) |
| hess-pairing | sum_i f_ii(U) f_ii(V) | 8 | 0 | fail (624 NF terms) |
| pluecker | sum w_k^4 | 8 | 0 | fail (666 NF terms) |
| pluecker | sum_k w_k^2 w_{k+1}^2 | 8 | 0 | fail (700 NF terms) |
| pluecker | HP-like: sum_k w_k^2 * lap(U)lap(V)^0 .. (sum w_k^2)*(sum_i f_ii(U)+f_ii(V)) | 8 | 0 | fail (762 NF terms) |
| polar | sum_i f_i(U) V_i^3 + swap | 8 | 0 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) V_i^3 - swap | 8 | 1 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) U_i^3 + swap | 8 | 0 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) U_i^3 - swap | 8 | 1 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) U_i^2 V_i + swap | 8 | 0 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) U_i^2 V_i - swap | 8 | 1 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) V_i^2 U_i + swap | 8 | 0 | PASS (Jacobian-type) |
| polar | sum_i f_i(U) V_i^2 U_i - swap | 8 | 1 | PASS (Jacobian-type) |
| third-polar | sum_ij P_ij(U;V) P_ij(V;U) | 8 | 0 | fail (750 NF terms) |
| third-polar | sum_ij P_ij(U;V) f_ij(V) + swap | 8 | 0 | fail (674 NF terms) |
| third-polar | sum_ij P_ij(U;V) f_ij(V) - swap | 8 | 1 | fail (685 NF terms) |
| third-polar | sum_ij P_ij(U;V)^2 + swap | 8 | 0 | fail (762 NF terms) |
| third-polar | sum_ij P_ij(U;V)^2 - swap | 8 | 1 | fail (762 NF terms) |
| polar-d9 | sum_i f_i(U) V_i^4 + swap | 9 | 0 | PASS (Jacobian-type) |
| polar-d9 | sum_i f_i(U) V_i^4 - swap | 9 | 1 | PASS (Jacobian-type) |
| polar-d9 | sum_i f_i(U) U_i^4 + swap | 9 | 0 | PASS (Jacobian-type) |
| polar-d9 | sum_i f_i(U) U_i^4 - swap | 9 | 1 | PASS (Jacobian-type) |
| polar-d9 | sum_i f_i(U) U_i^3 V_i + swap | 9 | 0 | PASS (Jacobian-type) |
| polar-d9 | sum_i f_i(U) U_i^3 V_i - swap | 9 | 1 | PASS (Jacobian-type) |
| control-product-d10 | HP * (U.V) | 10 | 0 | fail (456 NF terms) |
| polar-pairing-d10 | sum_i f_i(U) f_i(V) | 10 | 0 | PASS (Jacobian-type) |
| polar-pairing-d10 | sum_i f_i(U) f_{i+1}(V) + swap | 10 | 0 | PASS (Jacobian-type) |
| polar-pairing-d10 | sum_i f_i(U) f_{i+1}(V) - swap | 10 | 1 | PASS (Jacobian-type) |
| det-d11 | det[grad f(U), grad f(V), U] + swap | 11 | 0 | PASS (Jacobian-type) |
| det-d11 | det[grad f(U), grad f(V), U] - swap | 11 | 1 | PASS (Jacobian-type) |
| hess-det-d12 | det Hess f(U) + det Hess f(V) | 12 | 0 | fail (2 NF terms) |
| hess-det-d12 | det Hess f(U) - det Hess f(V) | 12 | 1 | fail (2 NF terms) |
| hess-det-d12 | det[Hess f(U) + Hess f(V)] | 12 | 0 | fail (218 NF terms) |
| hess-det-d12 | det[Hess f(U) - Hess f(V)] | 12 | 1 | fail (218 NF terms) |
| hess-det-d12 | det P(U;V) + swap | 12 | 0 | fail (198 NF terms) |
| hess-det-d12 | det P(U;V) - swap | 12 | 1 | fail (208 NF terms) |

Every non-Jacobian-type guess (Hessian cofactors/pairings, third polars, Fermat, Plücker, Hessian determinants at d = 12) fails E1 exactly over Q and is excluded; the polar families Σ f_i(U)·(cubic/quartic) ± swap, the polar pairings Σ f_i(U) f_σ(i)(V) and det[∇f(U), ∇f(V), U] ± swap are Jacobian-type (E1 automatic) and go to the exact E2 decision.

### 6c. Exact reduced-E2 decisions (mod p) computed locally

- control-wave1 | fU*U0^2 + fV*V0^2 (wave-1: inconsistent rank 6632 nullity 865): consistent=False rank=2216 nullity=229 (10152×2445, 73.4 s) — inconsistent, rank 6632, nullity 865 (= 636 G-directions + reduced nullity)
- control-spurious | g*(U0V0+U1V1) (Y-odd sector spurious point): consistent=False rank=2216 nullity=211 (10152×2427, 65.5 s) — inconsistent expected (sum principle); Y-odd code path validation
- polar | sum_i f_i(U) U_i^3 + swap: consistent=False rank=2287 nullity=158 (10152×2445, 66.6 s) 
- polar | sum_i f_i(U) U_i^3 - swap: consistent=False rank=2272 nullity=155 (10152×2427, 73.1 s) 
- polar | sum_i f_i(U) V_i^3 + swap: consistent=False rank=2287 nullity=158 (10152×2445, 67.5 s) 
- polar | sum_i f_i(U) V_i^3 - swap: consistent=False rank=2272 nullity=155 (10152×2427, 66.1 s) 
- polar | sum_i f_i(U) U_i^2 V_i + swap: consistent=False rank=2287 nullity=158 (10152×2445, 70.0 s) 
- polar | sum_i f_i(U) U_i^2 V_i - swap: consistent=False rank=2272 nullity=155 (10152×2427, 70.6 s) 
- polar | sum_i f_i(U) V_i^2 U_i + swap: consistent=False rank=2287 nullity=158 (10152×2445, 63.9 s) 
- polar | sum_i f_i(U) V_i^2 U_i - swap: consistent=False rank=2272 nullity=155 (10152×2427, 64.8 s) 
- control-product | g * (U.V)  [g = f(U)-f(V)]: consistent=False rank=2216 nullity=211 (10152×2427, 70.2 s) 

## Remote capacity note

The remote box has been saturated by other agents since ~19:30 (load 300–1150 on 128 cores, 45–90 jobs); the runner refuses new submissions above load 150, so the heavy jobs below were queued with `--wait` and started only when capacity allowed. Latest runner messages:
- `submit-full-d8p0-s1.out`: REMOTE BUSY: load 1735 (max 150), available memory 100 GB (min 60), 76 jobs running; waiting 60 s
- `submit-hess-d8.out`: REMOTE BUSY: load 1725 (max 150), available memory 100 GB (min 60), 76 jobs running; waiting 60 s
- `submit-sweep-a.out`: REMOTE BUSY: load 1725 (max 150), available memory 101 GB (min 60), 76 jobs running; waiting 60 s

## Jobs / reproduction

Remote job logs are mirrored in `/tmp/hodge/remote/_jobs/`; scripts in `/tmp/hodge/results/polarized-k3-remainder/w2_*.py`; queued/launched commands in `w2_jobs.txt`. Each job appends to its own `w2-*.jsonl` (merged by `w2_merge.py`).
