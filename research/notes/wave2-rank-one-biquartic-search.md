# Wave 2 search report: `rank-one-biquartic` (q = 1, A = X0^3 Y0^3)
**STATUS: NO WITNESS FOUND** (extensive wave-2 search; details below; search log `wave2-search.log.jsonl`).
Last rebuilt 2026-11-01 20:28:59.  Work dir `/tmp/hodge/results/rank-one-biquartic/wave2/` (scripts; `out/` job outputs on the remote, mirrored to `/tmp/hodge/remote/rank-one-biquartic/wave2/out/`).
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
one has `B_{D'}(N') - B_D(N0) = 6 A G D ... in J(F)` exactly (the cross terms are `A G D(a)`, `F G div N0 - G A D(a) - F N0(G)`
and `A F G^2`, all in J), so the Jacobian-quotient relaxation is a condition on `P(V / F S_{d-6})`: the spurious
F-multiples form the vertex of a cone and can be sliced away (one triangular coordinate per spurious direction set to 0,
`slice_ms.py`, invariance verified numerically: 0 violations for the 82 (d=6) and 432x3 (d=7) checks).  Emptiness of the
sliced projective scheme over F_q-bar (msolve: finite quotient) lifts to C by the Macaulay-rank semicontinuity.

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

## 2. Jacobian-quotient ("Hessian") relaxation ranks (exact mod p)
Necessary condition B_D(N0) in J(F)_{2d} on the full E1-admissible family (J-type generators m F_i with lifts A m e_i, plus kernel reps with their lifts). rank Phi = #pairs means no nonzero solution over F_p-bar, hence over C.

| d | parity | family (J + kernel) | #pairs | dim (S/J)_{2d} | rank Phi | excluded? |
|---|---|---|---|---|---|---|
| 5 | even | ? + ? = 3 | 6 | 744 | 6 | **YES (proved, mod p -> C)** |
| 5 | odd | 3 + 0 = 3 | 6 | 744 | 6 | **YES (proved, mod p -> C)** |
| 6 | even | 18 + 0 = 18 | 171 | 866 | 82 | no (rank-deficient by count) |
| 6 | odd | 18 + 0 = 18 | 171 | 866 | 108 | no (rank-deficient by count) |
| 7 | even | 63 + 3 = 66 | 2211 | 744 | 432 | no (rank-deficient by count) |

## 3. Full-E2 variable-projection multistart searches (normalized system, float, L-BFGS)
(pending)

### Residual landscape (per-restart files)
* `d7e_b0`: n=9, min 1.652e-02, 10% 1.652e-02, median 3.312e-02, 90% 1.228e-01, max 1.228e-01

## 4. Structured families vs the E1-admissible space (exact mod p)
(pending)

## 5. Candidates with residual < 1e-8 and exact confirmation attempts
None so far: no restart reached rel < 1e-8 on this instance (see landscape above).

## 6. Validation runs
* selftest: normalized-E2 float machinery validated on accepted witness coupled-pair-rational-d4-strict (residual 1e-16, Koszul recovery 1e-16, gradient check ok) and on lane identity B_{eps_X F}(A X)=8tAL^2M; D=eps_X F inconsistent (rel 0.60 reduced / 0.091 unreduced)
* planted witness search (solved instance coupled-pair-rational-d4-strict, d=4, K=5): 9/12 random restarts converged to rel<4e-11; rational reconstruction gave D = X0^4 - X0^2 Y0^2 + Y0^4; exact fixed-D solve + verify_witness succeeded (strict, Q=0). End-to-end pipeline (search -> reconstruct -> exact confirm) validated; 3/12 restarts stuck at a local minimum rel=0.35

## 7. Jobs
* `chainA.log`: chainA start Mon Nov  2 00:47:55 UTC 2026 / relax 5,6 done Mon Nov  2 00:48:43 UTC 2026
* `chainB.log`: chainB start Mon Nov  2 00:47:54 UTC 2026
* `chainC.log`: chainC start Mon Nov  2 00:47:55 UTC 2026
* `chainD.log`: chainD start Mon Nov  2 00:47:58 UTC 2026
