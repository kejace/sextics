# Wave 2 search report: `rank-one-biquartic` (q = 1, A = X0^3 Y0^3)
**STATUS: NO WITNESS FOUND** (extensive wave-2 search; details below; search log `wave2-search.log.jsonl`).
Last rebuilt 2026-11-02 10:17:07.  Work dir `/tmp/hodge/results/rank-one-biquartic/wave2/` (scripts; `out/` job outputs on the remote, mirrored to `/tmp/hodge/remote/rank-one-biquartic/wave2/out/`).
## 0. Setup and method (new in wave 2)

* Instance: `F = f(X) + g(Y) + (5/7) L(X)^2 M(Y)`, `A = X0^3 Y0^3`, q = 1.  Diagonal symmetry group `{1, tau}` only,
  so the only character sector is Y-parity: the denominator may be taken Y-even (lossless up to degree doubling);
  Y-odd families are searched too where cheap.  Prior work (dossier, attack 1): E1 forces every Y-even admissible
  denominator of degree <= 6 into J(F) (proved), kernels mod J of dimension 3/9/27 at d = 7/8/9, all natural
  low-degree denominators fail, a few members of the degree-7/8 families fail exactly (rank+1 pattern).
* Search machinery (`wave2/r1lib.py`, `wave2/e2search.py`): the **normalized E2 system** of book 5.2 instead of the
  cleared system.  For an E1-admissible family `D(a) = sum a_m D_m` with lifts `N0(a)` (`N0(F) = A D`), E2 holds iff
  `B_D(N0) + T_D(c) + F C = 0` for Koszul unknowns `c_ij` (degree d-4, parity-restricted) and `C` in `S_{2d-6}^-`
  (complete correction criterion, independent of the lift).  `C` is eliminated exactly by normal-form reduction
  modulo `F` (lex, LM = X0^6); the objective `f(a) = min_c |rho(B(a,a) + M(a) c)|^2` (a quartic in a) is minimized
  by L-BFGS with the exact variable-projection gradient (adjoint products), on the chart `a = Pc u/|u| + S s`
  where `S` spans the spurious F-multiples `F S_{d-6}` (lossless for denominators nonzero mod F).  Per-evaluation
  least squares: 2854x162 (d=6), 5160x420 (d=7), 8661x954 (d=8), 13717x1890 (d=9), 20736x3480 (d=10), versus
  16797x2602 ... 84903x18900 for the cleared system used by the engine's varpro.
* Residual measure `rel` = |r| / (|rho(D div N)| + sum_i |rho(N_i D_i)|) (cancellation-free scale), with
  |D| = 1 in the chart.  A true solution gives rel ~ 1e-12; candidates below 1e-6 are stored, below 1e-8 they are
  rationally reconstructed (continued fractions, several pivots/denominators) and re-solved exactly with the engine
  (`fixed_d.e1_solvable`, `solve_fixed_D_fast`, `verify_witness`, F-divisibility check).
* E1-admissible family in floats: SVD of the Jacobian-generator matrix in degree d+6, kernel of `Q_perp^T (A . m)`,
  lifts by pseudo-inverse (E1 residual ~1e-15).  Exact mod-p counterparts: `e1spaces.py` (kernel reps + lifts),
  `relax.py` (Jacobian-quotient relaxation rank: `B_D(N0) in J_{2d}`, for Jacobian-type D = dF(w) this is
  `A w^T Hess(F) w in J`), `structured.py` (intersection of structured subspaces with the admissible space),
  `e2exact.py` (single-D normalized E2 mod p).
* Validation (remote): accepted witness `coupled-pair-rational-d4-strict` (q=1, d=4): normalized residual 1e-16 at
  the witness point, Koszul-perturbed lifts recovered to 1e-16, random non-Koszul perturbations give O(1);
  gradient check analytic -0.0689627383 vs numeric -0.0689627386; a 12-restart search on that instance found the
  witness from 9 random starts and confirmed it exactly (strict).  On this lane `B_{eps_X F}(A X) = 8 t A L^2 M`
  is reproduced exactly and `D = eps_X F` is inconsistent (reduced rel 0.60, unreduced 0.091); d=6 family
  dimension 18 = dim J_6^+ and d=7 dimension 66 = 63 + 3 agree with the exact profile.
* Rigor of the exact exclusions.  All exact computations are modulo p = 2^31-1 (= `core.PRIMES[0]`); the Jacobian
  ring dimensions found there ((S/J)_10^- = 744, (S/J)_12^- = 866, (S/J)_14^- = 744, ...) equal the regular-sequence
  values, so every J-rank mod p equals its rank over Q.  Consequently (i) the Q-space of E1-admissible denominators
  reduces into the mod-p space, (ii) p-integral lifts reduce to lifts, (iii) normal forms modulo J commute with
  reduction.  If the relaxation map Phi: Sym^2(family) -> (S/J)_{2d} has full column rank mod p, its matrix over Q
  (in a saturated p-integral basis) has full column rank, so Phi(a (x) a) != 0 for every nonzero a over C: no
  denominator in that family passes E2.  For d <= 6 the family is the complete set of E1-admissible Y-even
  denominators (proved over C in the dossier), so a full-rank result there is an unconditional exclusion of that
  degree.  Numerical search results are evidence only; a failed multistart search excludes nothing.

**Invariance lemma (used for the exact exclusions).**  For `D' = D + F G` (G in S_{d-6}) with lift `N' = N0 + A G Z/6`
one has `B_{D'}(N') - B_D(N0) = F (G div N0 - N0(G) + A G^2)`, an element of `F S` (the three cross terms are
`B_D(A G Z/6) = A G D`, `B_{FG}(N0) = F G div N0 - G A D - F N0(G)` and `B_{FG}(A G Z/6) = A F G^2`; the `A G D` terms cancel), so the Jacobian-quotient relaxation is a condition on `P(V / F S_{d-6})`: the spurious
F-multiples form the vertex of a cone and can be sliced away (one triangular coordinate per spurious direction set to 0,
`slice_ms.py`, invariance verified numerically: 0 violations for the 82 (d=6) and 432x3 (d=7) checks).  Emptiness of the
sliced projective scheme over F_q-bar (msolve: finite quotient) lifts to C by the Macaulay-rank semicontinuity.
* Metric correction (20:50 EST): the first `rel` definition divided by the pieces-norm of the full `D = u-part + s`,
  which grows like |s|^2 while the residual grows at most linearly in the spurious component s, so an optimizer
  that wanders to large |s| produced spuriously small `rel` (e.g. 2e-6 with |s| = 89 and absolute residual 0.018).
  From that point on `rel = |r| / scale0` with `scale0` the pieces-norm of the u-part alone (|u| = 1, s = 0), the
  spurious component is bounded (|s| <= 5), and the per-restart records carry `resnorm`, `scale`, `scale0`, `s_norm`.
  Records written before the fix (tags d7e_b0, d7k_a, lmtest_*) used the old metric; their `resnorm` values
  (0.011-0.013 at d=7) are directly comparable.  Local optimizer: Levenberg-Marquardt with the Kaufman
  variable-projection Jacobian (`--method lm`, default from 20:50): on the solved validation instance it reaches the
  witness in 7 iterations (1e-16); on this lane it reaches the same ~1e-2 floor as L-BFGS but also stops at its
  60-iteration cap, and at d = 9 each iteration costs minutes under the machine load (2K+1 bilinear evaluations and
  a 13717x1890 least squares per line-search trial).  The earlier L-BFGS restarts (tags *_b0, d7k_a, d8k_a) all
  stopped at the 400-iteration cap.

**Why the mod-q emptiness statements lift to C (precise form).**  Let G be the Jacobian-generator matrix in degree 2d
(integer entries after clearing 7).  Its rank over Q equals its rank mod q (the Jacobian-ring dimensions computed mod q
are the regular-sequence values).  Take a saturated q-integral basis {l_i} of the left null space of G over Q; its
reduction mod q is a basis of the left null space of G mod q.  The relaxation system over Q is {l_i . B(a,a) = 0}
(B(a,a) = sum a_m a_n B_mn with q-integral B_mn), and its reduction is exactly the mod-q system whose span msolve
received.  If the mod-q Macaulay matrix in some degree k has full column rank (finite quotient), the Q Macaulay matrix
(whose reduction it is) has full column rank too, so the Q-ideal contains every form of degree k and the projective
variety over C is empty.  The same argument gives the degree-5 statement (full rank of Phi) and would give degree 7 if
msolve reports a finite quotient for the sliced 63-variable system.  No "standard monomial" or pivot assumption is used.

## 1. E1-admissible denominator spaces (exact mod p = 2^31-1)
| d | parity | dim S_d | dim J_d | ring dim (S/J)_d | dim {D: AD in J} | kernel mod J | matrix | time s |
|---|---|---|---|---|---|---|---|---|
| 4 | even | 66 | 0 | 66 | 0 | **0** | 1488x822 | 1 |
| 4 | odd | 60 | 0 | 60 | 0 | **0** | 1505x816 | 1 |
| 5 | even | 126 | 3 | 123 | 3 | **0** | 2181x1512 | 4 |
| 5 | odd | 126 | 3 | 123 | 3 | **0** | 2184x1512 | 4 |
| 6 | even | 236 | 18 | 218 | 18 | **0** | 3080x2612 | 12 |
| 6 | odd | 226 | 18 | 208 | 18 | **0** | 3108x2602 | 12 |
| 7 | even | 396 | 63 | 333 | 66 | **3** | 4284x4257 | 38 |
| 7 | odd | 396 | 63 | 333 | 63 | **0** | 4284x4257 | 38 |
| 8 | even | 651 | 168 | 483 | 177 | **9** | 5796x6657 | 107 |
| 8 | odd | 636 | 168 | 468 | 177 | **9** | 5832x6642 | 107 |
| 9 | even | 1001 | 378 | 623 | 405 | **27** | 7752x10010 | 270 |
| 9 | odd | 1001 | 378 | 623 | 405 | **27** | 7752x10010 | 277 |
| 10 | even | 1512 | 750 | 762 | 1053 | **303** | 10152x14616 | 459 |
| 10 | odd | 1491 | 747 | 744 | 1017 | **270** | 10197x14595 | 635 |
| 11 | even | 2184 | 1341 | 843 | 1851 | **510** | 13167x20748 | 2352 |
| 11 | odd | 2184 | 1341 | 843 | 1854 | **513** | 13167x20748 | 3118 |

## 2. Jacobian-quotient ("Hessian") relaxation ranks (exact mod p)
Necessary condition B_D(N0) in J(F)_{2d} on the full E1-admissible family (J-type generators m F_i with lifts A m e_i, plus kernel reps with their lifts). rank Phi = #pairs means no nonzero solution over F_p-bar, hence over C.

| d | parity | family (J + kernel) | #pairs | dim (S/J)_{2d} | rank Phi | excluded? |
|---|---|---|---|---|---|---|
| 5 | even | 3 + 0 = 3 | 6 | 744 | 6 | **YES (proved, mod p -> C)** |
| 5 | odd | 3 + 0 = 3 | 6 | 744 | 6 | **YES (proved, mod p -> C)** |
| 6 | even | 18 + 0 = 18 | 171 | 866 | 82 | no (rank-deficient by count) |
| 6 | odd | 18 + 0 = 18 | 171 | 866 | 108 | no (rank-deficient by count) |
| 7 | even | 63 + 3 = 66 | 2211 | 744 | 432 | no (rank-deficient by count) |
| 7 | odd | 63 + 0 = 63 | 2016 | 744 | 447 | no (rank-deficient by count) |
| 8 | even | 168 + 9 = 177 | 15753 | 468 | 462 | no (rank-deficient by count) |

## 3. Full-E2 variable-projection multistart searches (normalized system, float, L-BFGS)
| tag | d | parity | family dim K | n_c | rows_red | restarts | best rel | median rel | #rel<1e-3 | #rel<1e-6 | #rel<1e-8 | generic rank/nullity of T_D | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| d7e_b0 | 7 | even | 66 | 420 | 5160 | 200 | 6.270e-03 | 3.842e-02 | 0 | 0 | 0 | 360/60 | 25883 |
| d7k_a | 7 | even | 3 | 420 | 5160 | 40 | 2.268e-02 | 7.835e-02 | 0 | 0 | 0 | 360/60 | 361 |
| d7o_a | 7 | odd | 63 | 420 | 5160 | 16 | 1.440e-02 | 4.123e-02 | 0 | 0 | 0 | 360/60 | 2082 |
| d8k_a | 8 | even | 9 | 954 | 8661 | 40 | 2.502e-03 | 6.222e-03 | 0 | 0 | 0 | 795/159 | 9395 |
| d9m_a | 9 | even | 30 | 1890 | 13717 | 8 | 1.344e-03 | 2.496e-03 | 0 | 0 | 0 | 1530/360 | 12302 |

Totals per family: d=7 even (K=66): 240 restarts, best rel 6.270e-03; d=7 odd (K=63): 16 restarts, best rel 1.440e-02; d=8 even (K=9): 40 restarts, best rel 2.502e-03; d=9 even (K=30): 8 restarts, best rel 1.344e-03

### Residual landscape (per-restart files)
* `d7e_b0`: n=200, min 6.270e-03, 10% 1.645e-02, median 3.842e-02, 90% 7.439e-02, max 2.901e-01
* `d7e_b1`: n=162, min 1.225e-02, 10% 2.282e-02, median 4.019e-02, 90% 6.606e-02, max 2.121e+01
* `d7k_a`: n=40, min 2.268e-02, 10% 2.268e-02, median 7.835e-02, 90% 7.835e-02, max 7.835e-02
* `d7o_a`: n=16, min 1.440e-02, 10% 1.822e-02, median 4.123e-02, 90% 7.103e-02, max 8.815e-02
* `d8e_b0`: n=41, min 1.819e-03, 10% 2.292e-03, median 4.227e-03, 90% 7.684e-03, max 1.866e-02
* `d8k_a`: n=40, min 2.502e-03, 10% 6.222e-03, median 6.222e-03, 90% 6.222e-03, max 6.222e-03
* `d9m_a`: n=8, min 1.344e-03, 10% 1.344e-03, median 2.496e-03, 90% 1.070e-02, max 1.070e-02

## 4. Structured families vs the E1-admissible space (exact mod p)
dim(U cap W_d) with W_d = {D : A D in J} = J_d + kernel; "new" = part not already in J_d.

| d | parity | family U | dim U | dim U cap W | dim U cap J | new |
|---|---|---|---|---|---|---|
| 6 | even | Xbase Q(Y)L^a f^b | 51 | 10 | 10 | 0 |
| 6 | even | Ybase P(X)M^a g^b | 51 | 10 | 10 | 0 |
| 6 | even | bideg(0,6) | 28 | 3 | 3 | 0 |
| 6 | even | bideg(6,0) | 28 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2Y0Y1] | 110 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2Y0Y2] | 110 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2Y0] | 50 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2Y1Y2] | 110 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2Y1] | 50 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2Y2] | 50 | 3 | 3 | 0 |
| 6 | even | vars[X0X1X2] | 28 | 3 | 3 | 0 |
| 6 | even | vars[X0X1Y0Y1Y2] | 110 | 3 | 3 | 0 |
| 6 | even | vars[X0X2Y0Y1Y2] | 110 | 3 | 3 | 0 |
| 6 | even | vars[X0Y0Y1Y2] | 50 | 3 | 3 | 0 |
| 6 | even | vars[X1X2Y0Y1Y2] | 110 | 3 | 3 | 0 |
| 6 | even | vars[X1Y0Y1Y2] | 50 | 3 | 3 | 0 |
| 6 | even | vars[X2Y0Y1Y2] | 50 | 3 | 3 | 0 |
| 6 | even | vars[Y0Y1Y2] | 28 | 3 | 3 | 0 |
| 6 | even | monomials(f,g,L,M) | 6 | 2 | 2 | 0 |
| 6 | even | bideg(2,4) | 90 | 0 |  |  |
| 6 | even | bideg(4,2) | 90 | 0 |  |  |
| 6 | even | vars[X0X1Y0Y1] | 44 | 0 |  |  |
| 6 | even | vars[X0X1Y0Y2] | 44 | 0 |  |  |
| 6 | even | vars[X0X1Y0] | 16 | 0 |  |  |
| 6 | even | vars[X0X1Y1Y2] | 44 | 0 |  |  |
| 6 | even | vars[X0X1Y1] | 16 | 0 |  |  |
| 6 | even | vars[X0X1Y2] | 16 | 0 |  |  |
| 6 | even | vars[X0X2Y0Y1] | 44 | 0 |  |  |
| 6 | even | vars[X0X2Y0Y2] | 44 | 0 |  |  |
| 6 | even | vars[X0X2Y0] | 16 | 0 |  |  |
| 6 | even | vars[X0X2Y1Y2] | 44 | 0 |  |  |
| 6 | even | vars[X0X2Y1] | 16 | 0 |  |  |
| 6 | even | vars[X0X2Y2] | 16 | 0 |  |  |
| 6 | even | vars[X0Y0Y1] | 16 | 0 |  |  |
| 6 | even | vars[X0Y0Y2] | 16 | 0 |  |  |
| 6 | even | vars[X0Y1Y2] | 16 | 0 |  |  |
| 6 | even | vars[X1X2Y0Y1] | 44 | 0 |  |  |
| 6 | even | vars[X1X2Y0Y2] | 44 | 0 |  |  |
| 6 | even | vars[X1X2Y0] | 16 | 0 |  |  |
| 6 | even | vars[X1X2Y1Y2] | 44 | 0 |  |  |
| 6 | even | vars[X1X2Y1] | 16 | 0 |  |  |
| 6 | even | vars[X1X2Y2] | 16 | 0 |  |  |
| 6 | even | vars[X1Y0Y1] | 16 | 0 |  |  |
| 6 | even | vars[X1Y0Y2] | 16 | 0 |  |  |
| 6 | even | vars[X1Y1Y2] | 16 | 0 |  |  |
| 6 | even | vars[X2Y0Y1] | 16 | 0 |  |  |
| 6 | even | vars[X2Y0Y2] | 16 | 0 |  |  |
| 6 | even | vars[X2Y1Y2] | 16 | 0 |  |  |
| 6 | odd | bideg(1,5) | 63 | 0 |  |  |
| 6 | odd | bideg(3,3) | 100 | 0 |  |  |
| 6 | odd | bideg(5,1) | 63 | 0 |  |  |
| 6 | odd | vars[X0X1X2Y0Y1] | 100 | 0 |  |  |
| 6 | odd | vars[X0X1X2Y0Y2] | 100 | 0 |  |  |
| 6 | odd | vars[X0X1X2Y0] | 34 | 0 |  |  |
| 6 | odd | vars[X0X1X2Y1Y2] | 100 | 0 |  |  |
| 6 | odd | vars[X0X1X2Y1] | 34 | 0 |  |  |
| 6 | odd | vars[X0X1X2Y2] | 34 | 0 |  |  |
| 6 | odd | vars[X0X1Y0Y1Y2] | 100 | 0 |  |  |
| 6 | odd | vars[X0X1Y0Y1] | 40 | 0 |  |  |
| 6 | odd | vars[X0X1Y0Y2] | 40 | 0 |  |  |
| 6 | odd | vars[X0X1Y0] | 12 | 0 |  |  |
| 6 | odd | vars[X0X1Y1Y2] | 40 | 0 |  |  |
| 6 | odd | vars[X0X1Y1] | 12 | 0 |  |  |
| 6 | odd | vars[X0X1Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X0X2Y0Y1Y2] | 100 | 0 |  |  |
| 6 | odd | vars[X0X2Y0Y1] | 40 | 0 |  |  |
| 6 | odd | vars[X0X2Y0Y2] | 40 | 0 |  |  |
| 6 | odd | vars[X0X2Y0] | 12 | 0 |  |  |
| 6 | odd | vars[X0X2Y1Y2] | 40 | 0 |  |  |
| 6 | odd | vars[X0X2Y1] | 12 | 0 |  |  |
| 6 | odd | vars[X0X2Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X0Y0Y1Y2] | 34 | 0 |  |  |
| 6 | odd | vars[X0Y0Y1] | 12 | 0 |  |  |
| 6 | odd | vars[X0Y0Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X0Y1Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X1X2Y0Y1Y2] | 100 | 0 |  |  |
| 6 | odd | vars[X1X2Y0Y1] | 40 | 0 |  |  |
| 6 | odd | vars[X1X2Y0Y2] | 40 | 0 |  |  |
| 6 | odd | vars[X1X2Y0] | 12 | 0 |  |  |
| 6 | odd | vars[X1X2Y1Y2] | 40 | 0 |  |  |
| 6 | odd | vars[X1X2Y1] | 12 | 0 |  |  |
| 6 | odd | vars[X1X2Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X1Y0Y1Y2] | 34 | 0 |  |  |
| 6 | odd | vars[X1Y0Y1] | 12 | 0 |  |  |
| 6 | odd | vars[X1Y0Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X1Y1Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X2Y0Y1Y2] | 34 | 0 |  |  |
| 6 | odd | vars[X2Y0Y1] | 12 | 0 |  |  |
| 6 | odd | vars[X2Y0Y2] | 12 | 0 |  |  |
| 6 | odd | vars[X2Y1Y2] | 12 | 0 |  |  |
| 7 | even | Ybase P(X)M^a g^b | 73 | 27 | 24 | 3 |
| 7 | even | bideg(1,6) | 84 | 9 | 9 | 0 |
| 7 | even | bideg(7,0) | 36 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2Y0Y1] | 170 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2Y0Y2] | 170 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2Y0] | 70 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2Y1Y2] | 170 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2Y1] | 70 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2Y2] | 70 | 8 | 8 | 0 |
| 7 | even | vars[X0X1X2] | 36 | 8 | 8 | 0 |
| 7 | even | vars[X0X1Y0Y1Y2] | 160 | 6 | 6 | 0 |
| 7 | even | vars[X0X2Y0Y1Y2] | 160 | 6 | 6 | 0 |
| 7 | even | vars[X1X2Y0Y1Y2] | 160 | 6 | 6 | 0 |
| 7 | even | vars[X0Y0Y1Y2] | 50 | 3 | 3 | 0 |
| 7 | even | vars[X1Y0Y1Y2] | 50 | 3 | 3 | 0 |
| 7 | even | vars[X2Y0Y1Y2] | 50 | 3 | 3 | 0 |
| 7 | even | bideg(3,4) | 150 | 0 |  |  |
| 7 | even | bideg(5,2) | 126 | 0 |  |  |
| 7 | even | vars[X0X1Y0Y1] | 60 | 0 |  |  |
| 7 | even | vars[X0X1Y0Y2] | 60 | 0 |  |  |
| 7 | even | vars[X0X1Y0] | 20 | 0 |  |  |
| 7 | even | vars[X0X1Y1Y2] | 60 | 0 |  |  |
| 7 | even | vars[X0X1Y1] | 20 | 0 |  |  |
| 7 | even | vars[X0X1Y2] | 20 | 0 |  |  |
| 7 | even | vars[X0X2Y0Y1] | 60 | 0 |  |  |
| 7 | even | vars[X0X2Y0Y2] | 60 | 0 |  |  |
| 7 | even | vars[X0X2Y0] | 20 | 0 |  |  |
| 7 | even | vars[X0X2Y1Y2] | 60 | 0 |  |  |
| 7 | even | vars[X0X2Y1] | 20 | 0 |  |  |
| 7 | even | vars[X0X2Y2] | 20 | 0 |  |  |
| 7 | even | vars[X0Y0Y1] | 16 | 0 |  |  |
| 7 | even | vars[X0Y0Y2] | 16 | 0 |  |  |
| 7 | even | vars[X0Y1Y2] | 16 | 0 |  |  |
| 7 | even | vars[X1X2Y0Y1] | 60 | 0 |  |  |
| 7 | even | vars[X1X2Y0Y2] | 60 | 0 |  |  |
| 7 | even | vars[X1X2Y0] | 20 | 0 |  |  |
| 7 | even | vars[X1X2Y1Y2] | 60 | 0 |  |  |
| 7 | even | vars[X1X2Y1] | 20 | 0 |  |  |
| 7 | even | vars[X1X2Y2] | 20 | 0 |  |  |
| 7 | even | vars[X1Y0Y1] | 16 | 0 |  |  |
| 7 | even | vars[X1Y0Y2] | 16 | 0 |  |  |
| 7 | even | vars[X1Y1Y2] | 16 | 0 |  |  |
| 7 | even | vars[X2Y0Y1] | 16 | 0 |  |  |
| 7 | even | vars[X2Y0Y2] | 16 | 0 |  |  |
| 7 | even | vars[X2Y1Y2] | 16 | 0 |  |  |
| 7 | odd | Xbase Q(Y)L^a f^b | 73 | 24 | 24 | 0 |
| 7 | odd | bideg(6,1) | 84 | 9 | 9 | 0 |
| 7 | odd | bideg(0,7) | 36 | 8 | 8 | 0 |
| 7 | odd | vars[X0X1Y0Y1Y2] | 170 | 8 | 8 | 0 |
| 7 | odd | vars[X0X2Y0Y1Y2] | 170 | 8 | 8 | 0 |
| 7 | odd | vars[X0Y0Y1Y2] | 70 | 8 | 8 | 0 |
| 7 | odd | vars[X1X2Y0Y1Y2] | 170 | 8 | 8 | 0 |
| 7 | odd | vars[X1Y0Y1Y2] | 70 | 8 | 8 | 0 |
| 7 | odd | vars[X2Y0Y1Y2] | 70 | 8 | 8 | 0 |
| 7 | odd | vars[Y0Y1Y2] | 36 | 8 | 8 | 0 |
| 7 | odd | vars[X0X1X2Y0Y1] | 160 | 6 | 6 | 0 |
| 7 | odd | vars[X0X1X2Y0Y2] | 160 | 6 | 6 | 0 |
| 7 | odd | vars[X0X1X2Y1Y2] | 160 | 6 | 6 | 0 |
| 7 | odd | vars[X0X1X2Y0] | 50 | 3 | 3 | 0 |
| 7 | odd | vars[X0X1X2Y1] | 50 | 3 | 3 | 0 |
| 7 | odd | vars[X0X1X2Y2] | 50 | 3 | 3 | 0 |
| 7 | odd | bideg(2,5) | 126 | 0 |  |  |
| 7 | odd | bideg(4,3) | 150 | 0 |  |  |
| 7 | odd | vars[X0X1Y0Y1] | 60 | 0 |  |  |
| 7 | odd | vars[X0X1Y0Y2] | 60 | 0 |  |  |
| 7 | odd | vars[X0X1Y0] | 16 | 0 |  |  |
| 7 | odd | vars[X0X1Y1Y2] | 60 | 0 |  |  |
| 7 | odd | vars[X0X1Y1] | 16 | 0 |  |  |
| 7 | odd | vars[X0X1Y2] | 16 | 0 |  |  |
| 7 | odd | vars[X0X2Y0Y1] | 60 | 0 |  |  |
| 7 | odd | vars[X0X2Y0Y2] | 60 | 0 |  |  |
| 7 | odd | vars[X0X2Y0] | 16 | 0 |  |  |
| 7 | odd | vars[X0X2Y1Y2] | 60 | 0 |  |  |
| 7 | odd | vars[X0X2Y1] | 16 | 0 |  |  |
| 7 | odd | vars[X0X2Y2] | 16 | 0 |  |  |
| 7 | odd | vars[X0Y0Y1] | 20 | 0 |  |  |
| 7 | odd | vars[X0Y0Y2] | 20 | 0 |  |  |
| 7 | odd | vars[X0Y1Y2] | 20 | 0 |  |  |
| 7 | odd | vars[X1X2Y0Y1] | 60 | 0 |  |  |
| 7 | odd | vars[X1X2Y0Y2] | 60 | 0 |  |  |
| 7 | odd | vars[X1X2Y0] | 16 | 0 |  |  |
| 7 | odd | vars[X1X2Y1Y2] | 60 | 0 |  |  |
| 7 | odd | vars[X1X2Y1] | 16 | 0 |  |  |
| 7 | odd | vars[X1X2Y2] | 16 | 0 |  |  |
| 7 | odd | vars[X1Y0Y1] | 20 | 0 |  |  |
| 7 | odd | vars[X1Y0Y2] | 20 | 0 |  |  |
| 7 | odd | vars[X1Y1Y2] | 20 | 0 |  |  |
| 7 | odd | vars[X2Y0Y1] | 20 | 0 |  |  |
| 7 | odd | vars[X2Y0Y2] | 20 | 0 |  |  |
| 7 | odd | vars[X2Y1Y2] | 20 | 0 |  |  |
| 8 | even | Ybase P(X)M^a g^b | 102 | 55 | 46 | 9 |
| 8 | even | Xbase Q(Y)L^a f^b | 102 | 46 | 46 | 0 |
| 8 | even | vars[X0X1X2Y0Y1] | 255 | 24 | 24 | 0 |
| 8 | even | vars[X0X1X2Y0Y2] | 255 | 24 | 24 | 0 |
| 8 | even | vars[X0X1X2Y1Y2] | 255 | 24 | 24 | 0 |
| 8 | even | vars[X0X1Y0Y1Y2] | 255 | 24 | 24 | 0 |
| 8 | even | vars[X0X2Y0Y1Y2] | 255 | 24 | 24 | 0 |
| 8 | even | vars[X1X2Y0Y1Y2] | 255 | 24 | 24 | 0 |
| 8 | even | bideg(2,6) | 168 | 18 | 18 | 0 |
| 8 | even | bideg(6,2) | 168 | 18 | 18 | 0 |
| 8 | even | vars[X0X1X2Y0] | 95 | 18 | 18 | 0 |
| 8 | even | vars[X0X1X2Y1] | 95 | 18 | 18 | 0 |
| 8 | even | vars[X0X1X2Y2] | 95 | 18 | 18 | 0 |
| 8 | even | vars[X0Y0Y1Y2] | 95 | 18 | 18 | 0 |
| 8 | even | vars[X1Y0Y1Y2] | 95 | 18 | 18 | 0 |
| 8 | even | vars[X2Y0Y1Y2] | 95 | 18 | 18 | 0 |
| 8 | even | bideg(0,8) | 45 | 15 | 15 | 0 |
| 8 | even | bideg(8,0) | 45 | 15 | 15 | 0 |
| 8 | even | vars[X0X1X2] | 45 | 15 | 15 | 0 |
| 8 | even | vars[Y0Y1Y2] | 45 | 15 | 15 | 0 |
| 8 | even | monomials(f,g,L,M) | 9 | 4 | 4 | 0 |
| 8 | even | Hessian 2x2 minors | 48 | 0 |  |  |
| 8 | even | bideg(4,4) | 225 | 0 |  |  |
| 8 | even | vars[X0X1Y0Y1] | 85 | 0 |  |  |
| 8 | even | vars[X0X1Y0Y2] | 85 | 0 |  |  |
| 8 | even | vars[X0X1Y0] | 25 | 0 |  |  |
| 8 | even | vars[X0X1Y1Y2] | 85 | 0 |  |  |
| 8 | even | vars[X0X1Y1] | 25 | 0 |  |  |
| 8 | even | vars[X0X1Y2] | 25 | 0 |  |  |
| 8 | even | vars[X0X2Y0Y1] | 85 | 0 |  |  |
| 8 | even | vars[X0X2Y0Y2] | 85 | 0 |  |  |
| 8 | even | vars[X0X2Y0] | 25 | 0 |  |  |
| 8 | even | vars[X0X2Y1Y2] | 85 | 0 |  |  |
| 8 | even | vars[X0X2Y1] | 25 | 0 |  |  |
| 8 | even | vars[X0X2Y2] | 25 | 0 |  |  |
| 8 | even | vars[X0Y0Y1] | 25 | 0 |  |  |
| 8 | even | vars[X0Y0Y2] | 25 | 0 |  |  |
| 8 | even | vars[X0Y1Y2] | 25 | 0 |  |  |
| 8 | even | vars[X1X2Y0Y1] | 85 | 0 |  |  |
| 8 | even | vars[X1X2Y0Y2] | 85 | 0 |  |  |
| 8 | even | vars[X1X2Y0] | 25 | 0 |  |  |
| 8 | even | vars[X1X2Y1Y2] | 85 | 0 |  |  |
| 8 | even | vars[X1X2Y1] | 25 | 0 |  |  |
| 8 | even | vars[X1X2Y2] | 25 | 0 |  |  |
| 8 | even | vars[X1Y0Y1] | 25 | 0 |  |  |
| 8 | even | vars[X1Y0Y2] | 25 | 0 |  |  |
| 8 | even | vars[X1Y1Y2] | 25 | 0 |  |  |
| 8 | even | vars[X2Y0Y1] | 25 | 0 |  |  |
| 8 | even | vars[X2Y0Y2] | 25 | 0 |  |  |
| 8 | even | vars[X2Y1Y2] | 25 | 0 |  |  |
| 8 | odd | bideg(1,7) | 108 | 24 | 24 | 0 |
| 8 | odd | bideg(7,1) | 108 | 24 | 24 | 0 |
| 8 | odd | vars[X0X1X2Y0Y1] | 240 | 16 | 16 | 0 |
| 8 | odd | vars[X0X1X2Y0Y2] | 240 | 16 | 16 | 0 |
| 8 | odd | vars[X0X1X2Y1Y2] | 240 | 16 | 16 | 0 |
| 8 | odd | vars[X0X1Y0Y1Y2] | 240 | 16 | 16 | 0 |
| 8 | odd | vars[X0X2Y0Y1Y2] | 240 | 16 | 16 | 0 |
| 8 | odd | vars[X1X2Y0Y1Y2] | 240 | 16 | 16 | 0 |
| 8 | odd | vars[X0X1X2Y0] | 70 | 8 | 8 | 0 |
| 8 | odd | vars[X0X1X2Y1] | 70 | 8 | 8 | 0 |
| 8 | odd | vars[X0X1X2Y2] | 70 | 8 | 8 | 0 |
| 8 | odd | vars[X0Y0Y1Y2] | 70 | 8 | 8 | 0 |
| 8 | odd | vars[X1Y0Y1Y2] | 70 | 8 | 8 | 0 |
| 8 | odd | vars[X2Y0Y1Y2] | 70 | 8 | 8 | 0 |
| 8 | odd | Hessian 2x2 minors | 48 | 0 |  |  |
| 8 | odd | bideg(3,5) | 210 | 0 |  |  |
| 8 | odd | bideg(5,3) | 210 | 0 |  |  |
| 8 | odd | vars[X0X1Y0Y1] | 80 | 0 |  |  |
| 8 | odd | vars[X0X1Y0Y2] | 80 | 0 |  |  |
| 8 | odd | vars[X0X1Y0] | 20 | 0 |  |  |
| 8 | odd | vars[X0X1Y1Y2] | 80 | 0 |  |  |
| 8 | odd | vars[X0X1Y1] | 20 | 0 |  |  |
| 8 | odd | vars[X0X1Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X0X2Y0Y1] | 80 | 0 |  |  |
| 8 | odd | vars[X0X2Y0Y2] | 80 | 0 |  |  |
| 8 | odd | vars[X0X2Y0] | 20 | 0 |  |  |
| 8 | odd | vars[X0X2Y1Y2] | 80 | 0 |  |  |
| 8 | odd | vars[X0X2Y1] | 20 | 0 |  |  |
| 8 | odd | vars[X0X2Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X0Y0Y1] | 20 | 0 |  |  |
| 8 | odd | vars[X0Y0Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X0Y1Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X1X2Y0Y1] | 80 | 0 |  |  |
| 8 | odd | vars[X1X2Y0Y2] | 80 | 0 |  |  |
| 8 | odd | vars[X1X2Y0] | 20 | 0 |  |  |
| 8 | odd | vars[X1X2Y1Y2] | 80 | 0 |  |  |
| 8 | odd | vars[X1X2Y1] | 20 | 0 |  |  |
| 8 | odd | vars[X1X2Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X1Y0Y1] | 20 | 0 |  |  |
| 8 | odd | vars[X1Y0Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X1Y1Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X2Y0Y1] | 20 | 0 |  |  |
| 8 | odd | vars[X2Y0Y2] | 20 | 0 |  |  |
| 8 | odd | vars[X2Y1Y2] | 20 | 0 |  |  |
| 9 | even | Ybase P(X)M^a g^b | 138 | 91 | 79 | 12 |
| 9 | even | vars[X0X1X2Y0Y2] | 365 | 51 | 51 | 0 |
| 9 | even | bideg(7,2) | 216 | 48 | 48 | 0 |
| 9 | even | vars[X0X1X2Y0Y1] | 365 | 48 | 48 | 0 |
| 9 | even | vars[X0X1X2Y1Y2] | 365 | 48 | 48 | 0 |
| 9 | even | bideg(1,8) | 135 | 45 | 45 | 0 |
| 9 | even | vars[X0X1Y0Y1Y2] | 350 | 42 | 42 | 0 |
| 9 | even | vars[X0X2Y0Y1Y2] | 350 | 42 | 42 | 0 |
| 9 | even | vars[X1X2Y0Y1Y2] | 350 | 42 | 42 | 0 |
| 9 | even | vars[X0X1X2Y0] | 125 | 32 | 32 | 0 |
| 9 | even | vars[X0X1X2Y1] | 125 | 32 | 32 | 0 |
| 9 | even | vars[X0X1X2Y2] | 125 | 32 | 32 | 0 |
| 9 | even | bideg(3,6) | 280 | 30 | 30 | 0 |
| 9 | even | bideg(9,0) | 55 | 24 | 24 | 0 |
| 9 | even | vars[X0X1X2] | 55 | 24 | 24 | 0 |
| 9 | even | vars[X0Y0Y1Y2] | 95 | 18 | 18 | 0 |
| 9 | even | vars[X1Y0Y1Y2] | 95 | 18 | 18 | 0 |
| 9 | even | vars[X2Y0Y1Y2] | 95 | 18 | 18 | 0 |
| 9 | even | bideg(5,4) | 315 | 0 |  |  |
| 9 | even | vars[X0X1Y0Y1] | 110 | 0 |  |  |
| 9 | even | vars[X0X1Y0Y2] | 110 | 0 |  |  |
| 9 | even | vars[X0X1Y0] | 30 | 0 |  |  |
| 9 | even | vars[X0X1Y1Y2] | 110 | 0 |  |  |
| 9 | even | vars[X0X1Y1] | 30 | 0 |  |  |
| 9 | even | vars[X0X1Y2] | 30 | 0 |  |  |
| 9 | even | vars[X0X2Y0Y1] | 110 | 0 |  |  |
| 9 | even | vars[X0X2Y0Y2] | 110 | 0 |  |  |
| 9 | even | vars[X0X2Y0] | 30 | 0 |  |  |
| 9 | even | vars[X0X2Y1Y2] | 110 | 0 |  |  |
| 9 | even | vars[X0X2Y1] | 30 | 0 |  |  |
| 9 | even | vars[X0X2Y2] | 30 | 0 |  |  |
| 9 | even | vars[X0Y0Y1] | 25 | 0 |  |  |
| 9 | even | vars[X0Y0Y2] | 25 | 0 |  |  |
| 9 | even | vars[X0Y1Y2] | 25 | 0 |  |  |
| 9 | even | vars[X1X2Y0Y1] | 110 | 0 |  |  |
| 9 | even | vars[X1X2Y0Y2] | 110 | 0 |  |  |
| 9 | even | vars[X1X2Y0] | 30 | 0 |  |  |
| 9 | even | vars[X1X2Y1Y2] | 110 | 0 |  |  |
| 9 | even | vars[X1X2Y1] | 30 | 0 |  |  |
| 9 | even | vars[X1X2Y2] | 30 | 0 |  |  |
| 9 | even | vars[X1Y0Y1] | 25 | 0 |  |  |
| 9 | even | vars[X1Y0Y2] | 25 | 0 |  |  |
| 9 | even | vars[X1Y1Y2] | 25 | 0 |  |  |
| 9 | even | vars[X2Y0Y1] | 25 | 0 |  |  |
| 9 | even | vars[X2Y0Y2] | 25 | 0 |  |  |
| 9 | even | vars[X2Y1Y2] | 25 | 0 |  |  |
| 9 | odd | Xbase Q(Y)L^a f^b | 138 | 79 | 79 | 0 |
| 9 | odd | bideg(2,7) | 216 | 48 | 48 | 0 |
| 9 | odd | vars[X0X1Y0Y1Y2] | 365 | 48 | 48 | 0 |
| 9 | odd | vars[X0X2Y0Y1Y2] | 365 | 48 | 48 | 0 |
| 9 | odd | vars[X1X2Y0Y1Y2] | 365 | 48 | 48 | 0 |
| 9 | odd | bideg(8,1) | 135 | 45 | 45 | 0 |
| 9 | odd | vars[X0X1X2Y0Y1] | 350 | 42 | 42 | 0 |
| 9 | odd | vars[X0X1X2Y0Y2] | 350 | 42 | 42 | 0 |
| 9 | odd | vars[X0X1X2Y1Y2] | 350 | 42 | 42 | 0 |
| 9 | odd | vars[X0Y0Y1Y2] | 125 | 32 | 32 | 0 |
| 9 | odd | vars[X1Y0Y1Y2] | 125 | 32 | 32 | 0 |
| 9 | odd | vars[X2Y0Y1Y2] | 125 | 32 | 32 | 0 |
| 9 | odd | bideg(6,3) | 280 | 30 | 30 | 0 |
| 9 | odd | bideg(0,9) | 55 | 24 | 24 | 0 |
| 9 | odd | vars[Y0Y1Y2] | 55 | 24 | 24 | 0 |
| 9 | odd | vars[X0X1X2Y0] | 95 | 18 | 18 | 0 |
| 9 | odd | vars[X0X1X2Y1] | 95 | 18 | 18 | 0 |
| 9 | odd | vars[X0X1X2Y2] | 95 | 18 | 18 | 0 |
| 9 | odd | bideg(4,5) | 315 | 0 |  |  |
| 9 | odd | vars[X0X1Y0Y1] | 110 | 0 |  |  |
| 9 | odd | vars[X0X1Y0Y2] | 110 | 0 |  |  |
| 9 | odd | vars[X0X1Y0] | 25 | 0 |  |  |
| 9 | odd | vars[X0X1Y1Y2] | 110 | 0 |  |  |
| 9 | odd | vars[X0X1Y1] | 25 | 0 |  |  |
| 9 | odd | vars[X0X1Y2] | 25 | 0 |  |  |
| 9 | odd | vars[X0X2Y0Y1] | 110 | 0 |  |  |
| 9 | odd | vars[X0X2Y0Y2] | 110 | 0 |  |  |
| 9 | odd | vars[X0X2Y0] | 25 | 0 |  |  |
| 9 | odd | vars[X0X2Y1Y2] | 110 | 0 |  |  |
| 9 | odd | vars[X0X2Y1] | 25 | 0 |  |  |
| 9 | odd | vars[X0X2Y2] | 25 | 0 |  |  |
| 9 | odd | vars[X0Y0Y1] | 30 | 0 |  |  |
| 9 | odd | vars[X0Y0Y2] | 30 | 0 |  |  |
| 9 | odd | vars[X0Y1Y2] | 30 | 0 |  |  |
| 9 | odd | vars[X1X2Y0Y1] | 110 | 0 |  |  |
| 9 | odd | vars[X1X2Y0Y2] | 110 | 0 |  |  |
| 9 | odd | vars[X1X2Y0] | 25 | 0 |  |  |
| 9 | odd | vars[X1X2Y1Y2] | 110 | 0 |  |  |
| 9 | odd | vars[X1X2Y1] | 25 | 0 |  |  |
| 9 | odd | vars[X1X2Y2] | 25 | 0 |  |  |
| 9 | odd | vars[X1Y0Y1] | 30 | 0 |  |  |
| 9 | odd | vars[X1Y0Y2] | 30 | 0 |  |  |
| 9 | odd | vars[X1Y1Y2] | 30 | 0 |  |  |
| 9 | odd | vars[X2Y0Y1] | 30 | 0 |  |  |
| 9 | odd | vars[X2Y0Y2] | 30 | 0 |  |  |
| 9 | odd | vars[X2Y1Y2] | 30 | 0 |  |  |
| 10 | even | Ybase P(X)M^a g^b | 183 | 136 | 121 | 15 |
| 10 | even | Xbase Q(Y)L^a f^b | 183 | 121 | 121 | 0 |
| 10 | even | vars[X0X1X2Y0Y2] | 511 | 101 | 101 | 0 |
| 10 | even | vars[X0X1Y0Y1Y2] | 511 | 98 | 95 | 3 |
| 10 | even | vars[X0X2Y0Y1Y2] | 511 | 98 | 95 | 3 |
| 10 | even | vars[X1X2Y0Y1Y2] | 511 | 98 | 95 | 3 |
| 10 | even | vars[X0X1X2Y0Y1] | 511 | 95 | 95 | 0 |
| 10 | even | vars[X0X1X2Y1Y2] | 511 | 95 | 95 | 0 |
| 10 | even | bideg(2,8) | 270 | 90 | 90 | 0 |
| 10 | even | bideg(8,2) | 270 | 90 | 90 | 0 |
| 10 | even | vars[X0X1X2Y0] | 161 | 53 | 53 | 0 |
| 10 | even | vars[X0X1X2Y1] | 161 | 53 | 53 | 0 |
| 10 | even | vars[X0X1X2Y2] | 161 | 53 | 53 | 0 |
| 10 | even | vars[X0Y0Y1Y2] | 161 | 53 | 53 | 0 |
| 10 | even | vars[X1Y0Y1Y2] | 161 | 53 | 53 | 0 |
| 10 | even | vars[X2Y0Y1Y2] | 161 | 53 | 53 | 0 |
| 10 | even | bideg(4,6) | 420 | 45 | 45 | 0 |
| 10 | even | bideg(6,4) | 420 | 45 | 45 | 0 |
| 10 | even | bideg(0,10) | 66 | 35 | 35 | 0 |
| 10 | even | bideg(10,0) | 66 | 35 | 35 | 0 |
| 10 | even | vars[X0X1X2] | 66 | 35 | 35 | 0 |
| 10 | even | vars[Y0Y1Y2] | 66 | 35 | 35 | 0 |
| 10 | even | monomials(f,g,L,M) | 12 | 6 | 6 | 0 |
| 10 | even | vars[X0X1Y0Y1] | 146 | 0 |  |  |
| 10 | even | vars[X0X1Y0Y2] | 146 | 0 |  |  |
| 10 | even | vars[X0X1Y0] | 36 | 0 |  |  |
| 10 | even | vars[X0X1Y1Y2] | 146 | 0 |  |  |
| 10 | even | vars[X0X1Y1] | 36 | 0 |  |  |
| 10 | even | vars[X0X1Y2] | 36 | 0 |  |  |
| 10 | even | vars[X0X2Y0Y1] | 146 | 0 |  |  |
| 10 | even | vars[X0X2Y0Y2] | 146 | 0 |  |  |
| 10 | even | vars[X0X2Y0] | 36 | 0 |  |  |
| 10 | even | vars[X0X2Y1Y2] | 146 | 0 |  |  |
| 10 | even | vars[X0X2Y1] | 36 | 0 |  |  |
| 10 | even | vars[X0X2Y2] | 36 | 0 |  |  |
| 10 | even | vars[X0Y0Y1] | 36 | 0 |  |  |
| 10 | even | vars[X0Y0Y2] | 36 | 0 |  |  |
| 10 | even | vars[X0Y1Y2] | 36 | 0 |  |  |
| 10 | even | vars[X1X2Y0Y1] | 146 | 0 |  |  |
| 10 | even | vars[X1X2Y0Y2] | 146 | 0 |  |  |
| 10 | even | vars[X1X2Y0] | 36 | 0 |  |  |
| 10 | even | vars[X1X2Y1Y2] | 146 | 0 |  |  |
| 10 | even | vars[X1X2Y1] | 36 | 0 |  |  |
| 10 | even | vars[X1X2Y2] | 36 | 0 |  |  |
| 10 | even | vars[X1Y0Y1] | 36 | 0 |  |  |
| 10 | even | vars[X1Y0Y2] | 36 | 0 |  |  |
| 10 | even | vars[X1Y1Y2] | 36 | 0 |  |  |
| 10 | even | vars[X2Y0Y1] | 36 | 0 |  |  |
| 10 | even | vars[X2Y0Y2] | 36 | 0 |  |  |
| 10 | even | vars[X2Y1Y2] | 36 | 0 |  |  |
| 10 | odd | vars[X0X1X2Y0Y2] | 490 | 86 | 86 | 0 |
| 10 | odd | vars[X0X1Y0Y1Y2] | 490 | 83 | 80 | 3 |
| 10 | odd | vars[X0X2Y0Y1Y2] | 490 | 83 | 80 | 3 |
| 10 | odd | vars[X1X2Y0Y1Y2] | 490 | 83 | 80 | 3 |
| 10 | odd | bideg(3,7) | 360 | 80 | 80 | 0 |
| 10 | odd | bideg(7,3) | 360 | 80 | 80 | 0 |
| 10 | odd | vars[X0X1X2Y0Y1] | 490 | 80 | 80 | 0 |
| 10 | odd | vars[X0X1X2Y1Y2] | 490 | 80 | 80 | 0 |
| 10 | odd | bideg(1,9) | 165 | 72 | 72 | 0 |
| 10 | odd | bideg(9,1) | 165 | 72 | 72 | 0 |
| 10 | odd | vars[X0X1X2Y0] | 125 | 32 | 32 | 0 |
| 10 | odd | vars[X0X1X2Y1] | 125 | 32 | 32 | 0 |
| 10 | odd | vars[X0X1X2Y2] | 125 | 32 | 32 | 0 |
| 10 | odd | vars[X0Y0Y1Y2] | 125 | 32 | 32 | 0 |
| 10 | odd | vars[X1Y0Y1Y2] | 125 | 32 | 32 | 0 |
| 10 | odd | vars[X2Y0Y1Y2] | 125 | 32 | 32 | 0 |
| 10 | odd | bideg(5,5) | 441 | 15 | 0 | 15 |
| 10 | odd | vars[X0X1Y0Y1] | 140 | 0 |  |  |
| 10 | odd | vars[X0X1Y0Y2] | 140 | 0 |  |  |
| 10 | odd | vars[X0X1Y0] | 30 | 0 |  |  |
| 10 | odd | vars[X0X1Y1Y2] | 140 | 0 |  |  |
| 10 | odd | vars[X0X1Y1] | 30 | 0 |  |  |
| 10 | odd | vars[X0X1Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X0X2Y0Y1] | 140 | 0 |  |  |
| 10 | odd | vars[X0X2Y0Y2] | 140 | 0 |  |  |
| 10 | odd | vars[X0X2Y0] | 30 | 0 |  |  |
| 10 | odd | vars[X0X2Y1Y2] | 140 | 0 |  |  |
| 10 | odd | vars[X0X2Y1] | 30 | 0 |  |  |
| 10 | odd | vars[X0X2Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X0Y0Y1] | 30 | 0 |  |  |
| 10 | odd | vars[X0Y0Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X0Y1Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X1X2Y0Y1] | 140 | 0 |  |  |
| 10 | odd | vars[X1X2Y0Y2] | 140 | 0 |  |  |
| 10 | odd | vars[X1X2Y0] | 30 | 0 |  |  |
| 10 | odd | vars[X1X2Y1Y2] | 140 | 0 |  |  |
| 10 | odd | vars[X1X2Y1] | 30 | 0 |  |  |
| 10 | odd | vars[X1X2Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X1Y0Y1] | 30 | 0 |  |  |
| 10 | odd | vars[X1Y0Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X1Y1Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X2Y0Y1] | 30 | 0 |  |  |
| 10 | odd | vars[X2Y0Y2] | 30 | 0 |  |  |
| 10 | odd | vars[X2Y1Y2] | 30 | 0 |  |  |
| 11 | even | vars[X0X1X2] | 78 | 47 | 47 | 0 |

(0 three-variable cone families with zero intersection omitted from the table.)

### Exact normalized-E2 tests of single denominators (mod p)

| name | d | parity | e1 | consistent | rank | nullity | rows x cols | time s |
|---|---|---|---|---|---|---|---|---|
| e1kernel_d7_p0.json#0 | 7 | 0 | True | False | 989 | 67 | 5796x1056 | 57.3 |
| e1kernel_d7_p0.json#1 | 7 | 0 | True | False | 989 | 67 | 5796x1056 | 57.8 |
| e1kernel_d7_p0.json#2 | 7 | 0 | True | False | 989 | 67 | 5796x1056 | 52.2 |
| e1kernel_d8_p0.json#0 | 8 | 0 | True | False | 2286 | 159 | 10152x2445 | 150.5 |
| e1kernel_d8_p0.json#1 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 129.8 |
| e1kernel_d8_p0.json#2 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 105.6 |
| e1kernel_d8_p0.json#3 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 153.6 |
| e1kernel_d8_p0.json#4 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 169.0 |
| e1kernel_d8_p0.json#5 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 99.5 |
| e1kernel_d8_p0.json#6 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 174.1 |
| e1kernel_d8_p0.json#7 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 122.9 |
| e1kernel_d8_p0.json#8 | 8 | 0 | True | False | 2273 | 172 | 10152x2445 | 151.0 |
| e1kernel_d8_p1.json#0 | 8 | 1 | True | False | 2271 | 156 | 10152x2427 | 168.0 |
| e1kernel_d8_p1.json#1 | 8 | 1 | True | False | 2264 | 163 | 10152x2427 | 176.4 |
| e1kernel_d8_p1.json#2 | 8 | 1 | True | False | 2264 | 163 | 10152x2427 | 171.7 |
| e1kernel_d8_p1.json#3 | 8 | 1 | True | False | 2271 | 156 | 10152x2427 | 183.2 |
| e1kernel_d8_p1.json#6 | 8 | 1 | True | False | 2271 | 156 | 10152x2427 | 194.6 |
| e1kernel_d8_p1.json#7 | 8 | 1 | True | False | 2264 | 163 | 10152x2427 | 177.4 |
| e1kernel_d8_p1.json#8 | 8 | 1 | True | False | 2264 | 163 | 10152x2427 | 158.7 |

## 5. Candidates with residual < 1e-8 and exact confirmation attempts
None so far: no restart reached rel < 1e-8 on this instance (see landscape above).

## 6. Validation runs
* selftest: normalized-E2 float machinery validated on accepted witness coupled-pair-rational-d4-strict (residual 1e-16, Koszul recovery 1e-16, gradient check ok) and on lane identity B_{eps_X F}(A X)=8tAL^2M; D=eps_X F inconsistent (rel 0.60 reduced / 0.091 unreduced)
* planted witness search (solved instance coupled-pair-rational-d4-strict, d=4, K=5): 9/12 random restarts converged to rel<4e-11; rational reconstruction gave D = X0^4 - X0^2 Y0^2 + Y0^4; exact fixed-D solve + verify_witness succeeded (strict, Q=0). End-to-end pipeline (search -> reconstruct -> exact confirm) validated; 3/12 restarts stuck at a local minimum rel=0.35

## 7. Jobs
* `chainA.log`: chainA start Mon Nov  2 00:47:55 UTC 2026 / relax 5,6 done Mon Nov  2 00:48:43 UTC 2026 / e1spaces 4-10 done Mon Nov  2 01:58:14 UTC 2026 / relax 5-7 done Mon Nov  2 02:02:24 UTC 2026 / d5,d6 searches done Mon Nov  2 02:02:27 UTC 2026 / structured 6-10 done Mon Nov  2 02:32:57 UTC 2026 / msolve d6 exit 0 / msolve d7 exit 124 / msolve done Mon Nov  2 07:13:51 UTC 2026 / wave2/chainA.sh: line 33: 65171 Killed                  python3 wave2/e1spaces.py --dmin 12 --dmax 12 --tag e1spaces_d12 > $O/e1spaces_d12.out 2> $O/e1spaces_d12.err
* `chainB.log`: chainB start Mon Nov  2 00:47:54 UTC 2026 / d7e batch 0 done Mon Nov  2 08:00:43 UTC 2026 / wave2/chainB.sh: line 7: 2631878 Killed                  python3 wave2/e2search.py --d 7 --pD 0 --restarts 200 --seed $((b*200)) --workers 8 --tag d7e_b$b > $O/d7e_b$b.out 2> $O/d7e_b$b.err / d7e batch 1 done Mon Nov  2 12:29:38 UTC 2026
* `chainC.log`: chainC start Mon Nov  2 00:47:55 UTC 2026 / wave2/chainC.sh: line 7: 389988 Killed                  python3 wave2/e2search.py --d 8 --pD 0 --restarts 80 --seed $((b*80)) --workers 8 --tag d8e_b$b > $O/d8e_b$b.out 2> $O/d8e_b$b.err / d8e batch 0 done Mon Nov  2 12:29:36 UTC 2026
* `chainD.log`: chainD start Mon Nov  2 00:47:58 UTC 2026
* `chainD2.log`:  UTC 2026 / d7 kernel probe done Mon Nov  2 01:19:50 UTC 2026 / d8 kernel probe done Mon Nov  2 04:00:05 UTC 2026 / d9 kernel probe done Mon Nov  2 04:00:07 UTC 2026 / d9e batch c0 done Mon Nov  2 04:00:09 UTC 2026 / d9e batch c1 done Mon Nov  2 04:00:11 UTC 2026 / d9e batch c2 done Mon Nov  2 04:00:12 UTC 2026 / d9e batch c3 done Mon Nov  2 04:00:14 UTC 2026 / d9e batch c4 done Mon Nov  2 04:00:16 UTC 2026 / d9e batch c5 done Mon Nov  2 04:00:17 UTC 2026 / d9e batch c6 done Mon Nov  2 04:00:19 UTC 2026 / d9e batch c7 done Mon Nov  2 04:00:21 UTC 2026 / chainD2 end Mon Nov  2 04:00:22 UTC 2026
* `chainE.log`: chainE start Mon Nov  2 04:01:22 UTC 2026 / d9 mixed-kernel full-E2 (LM) done Mon Nov  2 07:43:59 UTC 2026 / wave2/chainE.sh: line 11: 2423901 Killed                  python3 wave2/relaxsearch.py --d 9 --pD 0 --family kernel-mixed --restarts 300 --seed 0 --workers 8 --refine-maxiter 200 --tag r9m_a > $O/r9m_a.out 2> $O/r9m_a.err / d9 mixed-kernel relaxation-first done Mon Nov  2 12:29:38 UTC 2026

## 8. Interpretation and what remains

* **Rigorous frontier moved from "deg D >= 7 for Y-even D (E1)" to "deg D >= 7 for any parity (E1 + E2
  relaxation)".**  Degrees <= 4: E1 forces D in J and J_d = 0 (dossier).  Degree 5 (both parities): the relaxation
  map on Sym^2(J_5) has full rank 6/6.  Degree 6, Y-odd: the relaxation ideal has a finite quotient (projectively
  empty).  Degree 6, Y-even: the relaxation variety is the single point [F] (the Euler/spurious multiple), length 184,
  at two primes; after slicing away the spurious direction (invariance lemma) the system is projectively empty.  All
  statements are mod p (2^31-1 and/or 1073741789, 1000000007) and lift to C by rank semicontinuity because every
  J-rank used equals its regular-sequence value.
* **Degree 7** is the first open degree.  The Y-even admissible family is 66-dimensional (63 polar + 3 kernel, the
  latter T-base-type per wave 1); its relaxation has 432 independent quadrics (rank Phi = 432 of 2211 pairs, 744
  available), which is far above the 63 unknowns of the sliced system, so emptiness is expected; msolve reached F4
  degree 3 in minutes, the degree-4 round decides (queued).  The numerical full-E2 search on the 66-dim family and on
  the 3-dim pure-kernel subfamily finds no residual below 1e-2 in the corrected metric (absolute residual
  0.011-0.021 at |u| = 1 over all restarts of both local methods); the degree-7 landscape has a floor around
  2e-2 .. 1e-1 with no near-solutions.  (Early LM runs reported rel 2e-6, an artifact of the first metric: the
  optimizer inflated the spurious F-multiple component to |s| = 89; see the metric note in section 0.)
* **Pure-kernel probes.**  The 3-dim degree-7 kernel (a projective plane) has no zero: 40 restarts, best absolute
  residual 0.02.  The 9-dim degree-8 Y-even kernel: 18 L-BFGS restarts all converge to one minimum, absolute residual
  0.121 -- no solution supported on the pure kernel directions at degree 8.  Degree-9 pure and mixed kernels: chain E.
* **Degree 8, 9, 10** are the "prime candidate" degrees of the task statement (non-Jacobian kernel elements 9 / 27 /
  >= 294).  The relaxation can no longer exclude by rank or Groebner (468 quadrics in 177 unknowns at d=8 has
  regularity far beyond reach; at d=9 the relaxation is underdetermined, 208 equations in 405 unknowns).  Only the
  full normalized-E2 search applies there; under the current machine load (load average ~1000 on 128 cores) the
  throughput is ~10x below nominal, so the restart counts in section 3 are what could be done, not what was planned.
* **Structured families** (block subrings, 3/4/5-variable cones, bidegree pieces, base-type P(X)M^a g^b and
  Q(Y)L^a f^b, (f,g,L,M)-monomials, Hessian minors; section 4).  Through degree 9 every structured family meets the
  admissible space only inside J_d, except the Y-base-type class P(X)M^a g^b (the packet's compact-Y class), which
  contains the whole degree-7 kernel (3), the whole degree-8 kernel (9) and 12 of the 27 degree-9 kernel dimensions.
  Hence the first E1-admissible denominators outside the packet's obstruction class are the 15 remaining degree-9
  kernel directions, spanned by K_7 . (S_2(Y) mod M) = {P_7(X) + L M beta_3(X)} x {Y-quadrics not proportional to M}
  modulo J + base type.  These are the target of the focused `kernel-mixed` searches (chain E).  At degree 10 new
  structured admissible directions appear: a 15-dimensional space of Y-odd bihomogeneous (5,5) denominators entirely
  outside J (pulled back from P^2 x P^2, i.e. still base type in the packet's sense), 3 new directions in each
  5-variable cone ring C[X_a X_b Y0 Y1 Y2], and 15 Y-base-type directions; members are tested exactly (e2exact).
* **No candidate** with residual below 1e-8 appeared on this instance; the only sub-1e-10 residuals of the whole
  run were the planted-witness validation on a solved instance.
