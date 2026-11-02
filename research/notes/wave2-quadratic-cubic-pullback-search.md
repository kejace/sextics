# Wave 2 search report: quadratic-cubic-pullback

STATUS: NO WITNESS FOUND.  (Report regenerated 2026-11-01 23:26 EST by w2_report.py from wave2-search.log.jsonl; 1705 records.)

Lane: F = sum_{i=0}^5 Q_i^3 (six Y-even quadrics), q = 1, A = X0^3 Y0^3.  Diagonal symmetry group of F: only Z/2 = {1, tau}
(Y-parity), so the only character sectors are the two Y-parities of D (a mixed-parity D can be replaced by D*tau(D)).
Certificate convention: A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q, deg N = d+1, deg Q = 2d-6, gcd(D, F) = 1.

## Headline results
1. Exact E1-admissible spaces K_d = {D : A D in J(F)} for d = 4..12, both parities (mod 2^31-1): K_d = J(F)_d for every
   d <= 9 (both parities); dim K_d/J_d = 294/261 (d=10), 510/510 (d=11), 677/648 (d=12) = dim R_d - dim R_{d+6} exactly
   (multiplication by [A] is surjective in those degrees).  Rank/nullity of the stacked E1 matrices in the table below.
2. Hessian form of the E2 relaxation (verified exactly): for D = sum_l w_l F_l, N0 = A w,
   D div N0 - N0(D) == -A (w^T Hess(F) w) (mod J(F)), so E2 mod J <=> w^T Hess(F) w in K_{2d-6}^even = (J(F):A)_{2d-6}.
3. Rigorous closure of degrees 5 and 6 (both parities): the stage-1 quadratic forms in the w-coordinates have no common
   projective zero (degree-4 Macaulay matrices of full column rank 4845/4845 and 5985/5985 exactly mod p, hence over Q;
   independently msolve's homogeneous reduced Groebner basis has v_i^2 in its leading ideal for every variable).
   Together with K_d = J_d (d <= 9): EVERY certificate for this target has deg D >= 7, in either parity.
   Closed exactly: d=5 even, d=5 odd, d=6 even, d=6 odd.
4. Structured families at d = 10 (pullbacks C[Q]_5 of quintic sections of the cubic, the 20 triples C[Q_S]_5, bihomogeneous
   products, delta * S_4, A * S_4, cones with vertex a plane or a line, all bigraded pieces, products of polars, unions,
   rho-/sigma-invariant subspaces): the only elements of K_10 they contain are F-multiples, except the generic-dimensional
   intersections with the big unions and the rho-/sigma-invariant subspaces (36 / 288 dimensional non-J subfamilies).
5. Numerical joint searches (stage 1 = Hessian relaxation with Levenberg-Marquardt, stage 2 = full E2 with the Koszul
   corrections T_D(c) and the quotient Q eliminated by variable projection; Fermat control gives residual 1e-33):
   d = 6: landscape minimum rel 1.9e-2 (361 restarts, 23 distinct minima), stage-2 minimum 1.6e-2 compressed (x ~8 unscaled).
   d = 7: every stage-1 restart converges to ONE minimizer (rank-1 set of all best points) at rel 1.32e-4, with
   w_0 ~ 0.46 X0^2 + 0.53 X1 X2 + 0.52 Y0 Y1 dominant (D ~ Q_0-like quadric times d_{X0}F, leading terms X0 Q_0^3); its full
   E2 (stage-2) residual is 3.7e-3 compressed (~3.5e-2 unscaled) and random stage-2 restarts give 2.6e-3 -- no certificate.
   d = 8 Y-odd (K_8 = J_8, 168 w's, 159 after gauge removal, 468 equations): the stage-1 landscape is a flat valley at
   rel ~ 1e-6 (18 restarts, all between 1.0e-6 and 2.0e-6, most stopping at the LM iteration cap); long LM refinement of the
   four best points converges (status 3) to genuine local minima at 9.57e-7, 9.57e-7, 1.03e-6, 1.08e-6 -- i.e. two orders of
   magnitude below d = 7 but NOT zeros (||P_J B||/||B|| ~ 1e-3).  Full-E2 stage 2 from these points is still running
   (job d8oddseed; chain job jvp-d89 continues with d = 8 odd stage 2 and d = 9 both parities); results append to
   w2_jvp_d8_odd_seed*.jsonl / w2_jvp_d8_odd_r1_s*.jsonl / w2_jvp_d9_*.jsonl on the remote lane dir.
   d = 8 Y-even, 9, 10: the stage-2 machinery was built (d = 10: K = 1050 parameters, 978 after gauge removal, 26532 rows,
   3480 Koszul + 5796 quotient columns) but the shared box ran at load 300-1750 on 128 cores; the d = 10 job was killed
   before completing a restart and d = 8 even / 9 are still queued in the chains.  No candidate reached 1e-10, so no
   rational reconstruction was attempted.

## Batch 1 (item 1): exact E1-admissible spaces K_d (job search-quadratic-cubic-pullback-kernels, mod p = 2^31-1)
Nullspace of [A*m | J(F)_{d+6} generators]; d <= 8 Y-even and d <= 6 Y-odd were closed in wave 1 (K_d = J_d).

| d | parity of D | dim S_d | dim J_d | dim R_d | dim K_d | dim K_d/J_d | dim R_{d+6} (target) | matrix | rank | s |
|---|---|---|---|---|---|---|---|---|---|---|
| 7 | odd | 396 | 63 | 333 | 63 | 0 | 843 | 4284x4257 | 3774 | 18.8 |
| 8 | odd | 636 | 168 | 468 | 168 | 0 | 762 | 5832x6642 | 5538 | 55.4 |
| 9 | even | 1001 | 378 | 623 | 378 | 0 | 623 | 7752x10010 | 7752 | 151.1 |
| 9 | odd | 1001 | 378 | 623 | 378 | 0 | 623 | 7752x10010 | 7752 | 150.9 |
| 10 | even | 1512 | 750 | 762 | 1044 | 294 | 468 | 10152x14616 | 10152 | 445.5 |
| 10 | odd | 1491 | 747 | 744 | 1008 | 261 | 483 | 10197x14595 | 10197 | 435.3 |
| 11 | even | 2184 | 1341 | 843 | 1851 | 510 | 333 | 13167x20748 | 13167 | 1418.5 |
| 11 | odd | 2184 | 1341 | 843 | 1851 | 510 | 333 | 13167x20748 | 13167 | 1481.2 |
| 12 | even | 3108 | 2223 | 885 | 2900 | 677 | 208 | 16797x28812 | 16797 | 4338.7 |
| 12 | odd | 3080 | 2214 | 866 | 2862 | 648 | 218 | 16852x28784 | 16852 | 4559.4 |

## Batch 2 + 5 (item 4, linear part): structured subspaces W of S_d^even against the exact K_d (jobs struct10, struct1112)
E1 for D is exactly D in K_d; "generic" = dim W + dim K_d - dim S_d (what a random subspace of that dimension would meet).

d = 10:

| family W | dim W | dim W cap K_10 | dim mod J | generic |
|---|---|---|---|---|
| UNION all cones with vertex a line | 1056 | 588 | 294 | 588 |
| sigma-invariant (X_i <-> Y_i) Y-even polynomials | 756 | 288 | 279 | 288 |
| UNION cones vertex a plane + extreme bidegrees | 672 | 204 | 204 | 204 |
| rho-invariant (cyclic shift of both blocks) Y-even polynomials | 504 | 36 | 36 | 36 |
| UNION C[Q]_k + delta*S + A*S + F_iF_j*S | 384 | 33 | 0 | 0 |
| C[Q]_5 (pullbacks of degree-5 hypersurface sections of the cubic) | 252 | 21 | 0 | 0 |
| UNION C[Q]_k + bihomogeneous + all C[Q_S]_k | 252 | 21 | 0 | 0 |
| F_i F_j * S_{d-10} (products of two polars) | 12 | 12 | 0 | 0 |
| C[Q_0,Q_1,Q_2]_a * C[Q_3,Q_4,Q_5]_b bihomogeneous products (a+b=5) | 210 | 9 | 0 | 0 |
| (rho,sigma)-invariant Y-even polynomials | 252 | 0 | 0 | 0 |
| A * S_4^odd | 60 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 2) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 3) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 4) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 2, 3) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 2, 4) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 2, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 3, 4) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 3, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(0, 4, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(1, 2, 3) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(1, 2, 4) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(1, 2, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(1, 3, 4) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(1, 3, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(1, 4, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(2, 3, 4) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(2, 3, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(2, 4, 5) | 21 | 0 | 0 | 0 |
| C[Q_S]_5 S=(3, 4, 5) | 21 | 0 | 0 | 0 |
| C[X0,X1,X2]_d (cones, vertex the Y-plane) | 66 | 0 | 0 | 0 |
| C[Y0,Y1,Y2]_d (cones, vertex the X-plane) | 66 | 0 | 0 | 0 |
| C[z_(0, 1, 2, 3)]_d (cone, vertex line z_4=z_5=0) | 161 | 0 | 0 | 0 |
| C[z_(0, 1, 2, 4)]_d (cone, vertex line z_3=z_5=0) | 161 | 0 | 0 | 0 |
| C[z_(0, 1, 2, 5)]_d (cone, vertex line z_3=z_4=0) | 161 | 0 | 0 | 0 |
| C[z_(0, 1, 3, 4)]_d (cone, vertex line z_2=z_5=0) | 146 | 0 | 0 | 0 |
| C[z_(0, 1, 3, 5)]_d (cone, vertex line z_2=z_4=0) | 146 | 0 | 0 | 0 |
| C[z_(0, 1, 4, 5)]_d (cone, vertex line z_2=z_3=0) | 146 | 0 | 0 | 0 |
| C[z_(0, 2, 3, 4)]_d (cone, vertex line z_1=z_5=0) | 146 | 0 | 0 | 0 |
| C[z_(0, 2, 3, 5)]_d (cone, vertex line z_1=z_4=0) | 146 | 0 | 0 | 0 |
| C[z_(0, 2, 4, 5)]_d (cone, vertex line z_1=z_3=0) | 146 | 0 | 0 | 0 |
| C[z_(0, 3, 4, 5)]_d (cone, vertex line z_1=z_2=0) | 161 | 0 | 0 | 0 |
| C[z_(1, 2, 3, 4)]_d (cone, vertex line z_0=z_5=0) | 146 | 0 | 0 | 0 |
| C[z_(1, 2, 3, 5)]_d (cone, vertex line z_0=z_4=0) | 146 | 0 | 0 | 0 |
| C[z_(1, 2, 4, 5)]_d (cone, vertex line z_0=z_3=0) | 146 | 0 | 0 | 0 |
| C[z_(1, 3, 4, 5)]_d (cone, vertex line z_0=z_2=0) | 161 | 0 | 0 | 0 |
| C[z_(2, 3, 4, 5)]_d (cone, vertex line z_0=z_1=0) | 161 | 0 | 0 | 0 |
| S3 x (X<->Y) shape-symmetric Y-even polynomials | 140 | 0 | 0 | 0 |
| bidegree (0,10) monomials | 66 | 0 | 0 | 0 |
| bidegree (10,0) monomials | 66 | 0 | 0 | 0 |
| bidegree (2,8) monomials | 270 | 0 | 0 | 0 |
| bidegree (4,6) monomials | 420 | 0 | 0 | 0 |
| bidegree (6,4) monomials | 420 | 0 | 0 | 0 |
| bidegree (8,2) monomials | 270 | 0 | 0 | 0 |
| delta * S_4^odd | 60 | 0 | 0 | 0 |

d = 12 (if present):

| family W | dim W | dim W cap K_12 | dim mod J | generic |
|---|---|---|---|---|
| C[Q]_6 (pullbacks of degree-6 hypersurface sections of the cubic) | 462 | 254 | 198 | 254 |
| C[Q_0,Q_1,Q_2]_a * C[Q_3,Q_4,Q_5]_b bihomogeneous products (a+b=6) | 406 | 198 | 162 | 198 |
| delta * S_6^odd | 226 | 54 | 0 | 18 |
| C[z_(0, 1, 2, 3)]_d (cone, vertex line z_4=z_5=0) | 252 | 44 | 44 | 44 |
| C[z_(0, 1, 2, 4)]_d (cone, vertex line z_3=z_5=0) | 252 | 44 | 44 | 44 |
| C[z_(0, 1, 2, 5)]_d (cone, vertex line z_3=z_4=0) | 252 | 44 | 44 | 44 |
| C[z_(0, 3, 4, 5)]_d (cone, vertex line z_1=z_2=0) | 252 | 44 | 44 | 44 |
| C[z_(0, 1, 3, 4)]_d (cone, vertex line z_2=z_5=0) | 231 | 23 | 23 | 23 |
| C[z_(0, 1, 3, 5)]_d (cone, vertex line z_2=z_4=0) | 231 | 23 | 23 | 23 |
| C[z_(0, 1, 4, 5)]_d (cone, vertex line z_2=z_3=0) | 231 | 23 | 23 | 23 |
| C[z_(0, 2, 3, 4)]_d (cone, vertex line z_1=z_5=0) | 231 | 23 | 23 | 23 |
| C[z_(0, 2, 3, 5)]_d (cone, vertex line z_1=z_4=0) | 231 | 23 | 23 | 23 |
| C[z_(0, 2, 4, 5)]_d (cone, vertex line z_1=z_3=0) | 231 | 23 | 23 | 23 |
| C[z_(1, 2, 3, 4)]_d (cone, vertex line z_0=z_5=0) | 231 | 23 | 23 | 23 |
| C[Q_S]_6 S=(0, 1, 2) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 1, 3) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 1, 4) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 1, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 2, 3) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 2, 4) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 2, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 3, 4) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 3, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(0, 4, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(1, 2, 3) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(1, 2, 4) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(1, 2, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(1, 3, 4) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(1, 3, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(1, 4, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(2, 3, 4) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(2, 3, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(2, 4, 5) | 28 | 0 | 0 | 0 |
| C[Q_S]_6 S=(3, 4, 5) | 28 | 0 | 0 | 0 |
| C[X0,X1,X2]_d (cones, vertex the Y-plane) | 91 | 0 | 0 | 0 |
| C[Y0,Y1,Y2]_d (cones, vertex the X-plane) | 91 | 0 | 0 | 0 |
| delta^2 * S_{d-12}^even | 1 | 0 | 0 | 0 |

d = 11 (if present):

| family W | dim W | dim W cap K_11 | dim mod J | generic |
|---|---|---|---|---|
| (pending) | | | | |

## Hessian form of stage 1 (verified exactly, both parities; record item "hessian-identity")
For D = sum_l w_l F_l and N0 = A w: D div N0 - N0(D) == -A (w^T Hess(F) w) (mod J(F)).  Hence the Jacobian-quotient
relaxation of E2 ("Hessian relaxation") is exactly w^T Hess(F) w in K_{2d-6}^even; for d <= 7 that is membership in
J(F)_{2d-6} (independent of A), and the number of independent quadratic equations is dim R_{2d-6}^even (218, 483, 468, ...)
rather than dim R_{2d}^odd.

## Batch 3 (item 2, exact closure of Jacobian-type degrees; w2_closure.py, w2_points.py, job closure678)
At degrees with K_d = J_d every E1-admissible D is sum_l w_l F_l; E2 implies the stage-1 system q(w) = 0, rJ = dim R_{2d}^odd
quadratic forms in the w-coordinates, invariant under the spurious translations D -> D + F G (quotient coordinates v, dim n_q).
Full column rank of the degree-k Macaulay matrix (exact mod p => over Q) certifies that the forms have no common zero over C,
i.e. no denominator of that degree and parity admits a certificate.

| d | parity | n_w | spurious | n_q | rJ | deg-2 rank | higher Macaulay ranks | closed at k |
|---|---|---|---|---|---|---|---|---|
| 5 | even | 3 | 0 | 3 | 744 | 6/6 |  | 2 |
| 5 | odd | 3 | 0 | 3 | 744 | 6/6 |  | 2 |
| 6 | even | 18 | 1 | 17 | 866 | 135/153 | deg3: 967/969 deg4: 4845/4845 | 4 |
| 6 | odd | 18 | 0 | 18 | 866 | 153/171 | deg3: 1138/1140 deg4: 5985/5985 | 4 |

msolve (F4 mod p = 2147483629):

| d | parity | char | system | mode | verdict | time |
|---|---|---|---|---|---|---|
| 6 | even | 2147483629 | 135 forms in 17 vars | affine chart l(v)=1 | no solution in the random affine chart (GB = [1]) | 0.0 s |
| 6 | odd | 2147483629 | 153 forms in 18 vars | affine chart l(v)=1 | no solution in the random affine chart (GB = [1]) | 0.0 s |
| 6 | odd | 2147483629 | 153 forms in 18 vars | homogeneous, -g 2 | EMPTY projective zero set (pure power of every variable in the leading ideal) | 0.1 s |
| 7 | even | 2147483629 | 483 forms in 60 vars | affine chart l(v)=1 | msolve killed (SIGKILL, exit -9) after 6728 s at ~22 GB RSS on the saturated box: no verdict; degree | 6728.4 s |
| 7 | even | 2147483629 | 483 forms in 60 vars | affine chart l(v)=1 | see file | 6728.4 s |

## Batch 4 + 6 (items 2/3, numerical joint search; jobs jvp-d678, d10, jvp-d89): residual landscape
rel = ||residual||^2 / ||B_D(N0)||^2 at the converged point (0 would be a certificate).  F-multiples and the Koszul gauge
are projected out; |u| = 1.  Stage-2 values use a random row compression S (k = n_c + 200 rows), which scales norms by
about k/n_out (0.12 at d = 6); exact stage-2 residuals are recomputed only below 1e-6 (never triggered).
d = 10 real kernel (job d10): dim K = 1044, dim J = 750, non-J part 294; singular values of the Jacobian block at the cut 0.563 / 5.9e-14; E1-lift residual max 1.5e-14.

| d | parity | objective | relative residual statistics | distinct minima | cost/restart |
|---|---|---|---|---|---|
| 6 | even | stage 1 (Hessian relaxation, LM) | n=361 min=1.890e-02 q10=2.303e-02 median=3.463e-02 max=4.660e-02 | 23 | 27 s |
| 7 | even | stage 1 (Hessian relaxation, LM) | n=80 min=1.318e-04 q10=1.318e-04 median=1.318e-04 max=1.843e-04 | 2 | 316 s |
| 8 | even | stage 1 (Hessian relaxation, LM) | n=17 min=1.140e-06 q10=1.179e-06 median=1.555e-06 max=2.468e-06 | 1 | 1092 s |
| 8 | odd | stage 1 (Hessian relaxation, LM) | n=28 min=9.571e-07 q10=1.002e-06 median=1.647e-06 max=3.130e-06 | 1 | 985 s |
| 6 | even | stage 2 (full E2 varpro, L-BFGS; compressed rows) | n=34 min=1.635e-02 q10=1.972e-02 median=2.418e-02 max=3.052e-02 | 27 | 303 s |
| 7 | even | stage 2 (full E2 varpro, L-BFGS; compressed rows) | n=14 min=2.412e-03 q10=2.626e-03 median=3.053e-03 max=3.987e-03 | 10 | 2395 s |

## What remains open after wave 2
* Degree 7 (both parities): stage 1 is 483 quadrics in 60 variables; the degree-3 Macaulay matrix cannot be full (28980 < 37820
  columns) and the degree-4 one (884k x 596k) exceeds dense exact linear algebra; msolve F4 on it is the pending exact route.
* Degrees 8, 9 (Jacobian-type, both parities) and 10 (K_10 = J_10 + 294 non-J): only numerical evidence; the stage-2 systems are
  overdetermined by factors 10-30, and the remote box was saturated, so restart counts stayed far below the planned thousands.
* Structured non-J subfamilies at d = 10 worth finishing first: rho-invariant (36-dim) and sigma-invariant (288-dim) parts of K_10
  (w2_jvp.py --restrict rho|sigma with the real kernel npz; built but not run).
* At d = 12 the E1 condition has codimension 208 only, so every family tested meets K_12 in its generic dimension (C[Q]_6: 254,
  of which 198 non-J; cones with vertex a coordinate line: 23-44 each); E1 stops being a useful filter there and the stage-2
  system (59332 rows) needs sparse exact/iterative linear algebra.
* Wave-1 fixed-D rank/nullity statistics for reference (engine 2, compressed systems): degree-6 E1-failing denominators rank 2214 /
  nullity 162; joins with target delta (consistent) 2358 / 244; delta itself 2223 / 153; Fermat control 2322 / 280; degree-5 polars
  1318 / 128; degree-6 polars served space 18 (nothing).  Exact E1 stacks: see the rank column of the kernel table.
* d = 11, 12 kernels are available (w2_kernel_<d>_<par>.json on the remote mirror); the full E2 systems there (40k-59k rows)
  need a sparse/iterative solver before any joint search is realistic.

## Reproduction
Lane dir /tmp/hodge/results/quadratic-cubic-pullback/: w2_kernels.py (exact K_d), w2_struct.py (structured intersections), w2_closure.py
(Macaulay certificates), w2_points.py (msolve export/verdicts), w2_realkernel.py (real non-J kernel + E1 lifts),
w2_jvp.py (stage-1/2 variable-projection search; Fermat control: F0 = sum Z^6, D = X0^6+X1^6+X2^6 gives residual 1e-33),
w2_merge.py / w2_report.py (log merge, this report), w2_queue.sh (remote submissions).  Remote mirror: /tmp/hodge/remote/quadratic-cubic-pullback/.
Log: /tmp/hodge/results/quadratic-cubic-pullback/wave2-search.log.jsonl (one record per batch/restart; search vectors kept in the mirror files).
