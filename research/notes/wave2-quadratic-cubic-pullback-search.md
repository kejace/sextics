# Wave 2 search report: quadratic-cubic-pullback

STATUS: NO WITNESS FOUND (search starting 2026-11-01 ~19:05 EST; this file is rewritten after every batch).

Lane: F = sum_{i=0}^5 Q_i^3 (six Y-even quadrics), q = 1, A = X0^3 Y0^3. Symmetry group of F: only Z/2 = {1, tau} (Y-parity),
so "character sectors" = Y-parity; the target forces D Y-even (Y-odd D handled separately by averaging, see dossier 1.6).

## Plan (wave 2)
1. Exact E1-admissible spaces K_d = {D : A D in J(F)} per degree d = 4..12, both parities (mod 2^31-1), with bases saved.
2. Jacobian-type degrees (K_d = J_d, proven for Y-even d <= 8): search the w-space D = sum w_i F_i directly with
   (a) stage 1 = Jacobian-quotient ("Hessian") relaxation of E2, (b) stage 2 = full E2 with Koszul-correction unknowns
   (variable projection), thousands of random restarts on the remote.
3. First non-Jacobian degree (d = 10 expected, kernel dim 294 Y-even): real basis of the non-J kernel with exact E1 lifts,
   joint stage-2 search over K_10 = J_10 + nonJ.
4. Structured families at degrees 6-12 (linear intersections with K_d; explicit sweep through E1 + full solve).
5. Rational reconstruction + exact re-solve for any candidate with residual < 1e-10.

## Batches (filled in as they finish)

### Batch 1 (item 1): exact E1-admissible spaces (remote job search-quadratic-cubic-pullback-kernels, mod p = 2^31-1)
K_d^par = {D in S_d^par : A D in J(F)} via nullspace of [A*m | J(F)_{d+6} generators]; d <= 8 Y-even and d <= 6 Y-odd were
already closed in wave 1 (K_d = J_d).  New (this batch):

| d | parity of D | dim S_d | dim J_d | dim R_d | dim K_d | dim K_d/J_d | dim R_{d+6} (target) | matrix | s |
|---|---|---|---|---|---|---|---|---|---|
| 7 | odd | 396 | 63 | 333 | 63 | 0 | 843 | 4284x4257 | 18.8 |
| 8 | odd | 636 | 168 | 468 | 168 | 0 | 762 | 5832x6642 | 55.4 |
| 9 | even | 1001 | 378 | 623 | 378 | 0 | 623 | 7752x10010 | 151.1 |
| 9 | odd | 1001 | 378 | 623 | 378 | 0 | 623 | 7752x10010 | 150.9 |
| 10 | even | 1512 | 750 | 762 | 1044 | 294 | 468 | 10152x14616 | 445.5 |
| 10 | odd | 1491 | 747 | 744 | 1008 | 261 | 483 | 10197x14595 | 435.3 |
| 11 | even | 2184 | 1341 | 843 | 1851 | 510 | 333 | 13167x20748 | 1418.5 |
| 11 | odd | 2184 | 1341 | 843 | 1851 | 510 | 333 | 13167x20748 | 1481.2 |
| 12 | even | 3108 | 2223 | 885 | 2900 | 677 | 208 | 16797x28812 | 4338.7 |
| 12 | odd | 3080 | 2214 | 866 | 2862 | 648 | 218 | 16852x28784 | 4559.4 |

Conclusions: every E1-admissible denominator of degree <= 9 (both parities) lies in J(F) = (F_0..F_5); the first
non-Jacobian degree is d = 10 with exactly the Gorenstein-duality numbers (294 Y-even, 261 Y-odd); d = 11: 510 / 510, d = 12: 677 / 648
(all equal to dim R_d - dim R_{d+6}, i.e. multiplication by [A] is surjective onto R_{d+6} in these degrees).

### Batch 2 (item 4, linear part): structured subspaces W of S_10^even against the exact K_10 (job struct10, mod p)
E1 for D is exactly D in K_10, so dim(W cap K_10) counts the family members passing E1; the last column is the part not in J(F).

| family W | dim W | dim W cap K_10 | dim mod J |
|---|---|---|---|
| C[Q]_5 (pullbacks of degree-5 hypersurface sections of the cubic) | 252 | 21 | 0 |
| F_i F_j * S_{d-10} (products of two polars) | 12 | 12 | 0 |
| C[Q_0,Q_1,Q_2]_a * C[Q_3,Q_4,Q_5]_b bihomogeneous products (a+b=5) | 210 | 9 | 0 |
| A * S_4^odd | 60 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 2) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 3) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 4) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 1, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 2, 3) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 2, 4) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 2, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 3, 4) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 3, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(0, 4, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(1, 2, 3) | 21 | 0 | 0 |
| C[Q_S]_5 S=(1, 2, 4) | 21 | 0 | 0 |
| C[Q_S]_5 S=(1, 2, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(1, 3, 4) | 21 | 0 | 0 |
| C[Q_S]_5 S=(1, 3, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(1, 4, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(2, 3, 4) | 21 | 0 | 0 |
| C[Q_S]_5 S=(2, 3, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(2, 4, 5) | 21 | 0 | 0 |
| C[Q_S]_5 S=(3, 4, 5) | 21 | 0 | 0 |
| C[X0,X1,X2]_d (cones, vertex the Y-plane) | 66 | 0 | 0 |
| C[Y0,Y1,Y2]_d (cones, vertex the X-plane) | 66 | 0 | 0 |
| C[z_(0, 1, 2, 3)]_d (cone, vertex line z_4=z_5=0) | 161 | 0 | 0 |
| C[z_(0, 1, 2, 4)]_d (cone, vertex line z_3=z_5=0) | 161 | 0 | 0 |
| C[z_(0, 1, 2, 5)]_d (cone, vertex line z_3=z_4=0) | 161 | 0 | 0 |
| C[z_(0, 1, 3, 4)]_d (cone, vertex line z_2=z_5=0) | 146 | 0 | 0 |
| C[z_(0, 1, 3, 5)]_d (cone, vertex line z_2=z_4=0) | 146 | 0 | 0 |
| C[z_(0, 1, 4, 5)]_d (cone, vertex line z_2=z_3=0) | 146 | 0 | 0 |
| C[z_(0, 2, 3, 4)]_d (cone, vertex line z_1=z_5=0) | 146 | 0 | 0 |
| C[z_(0, 2, 3, 5)]_d (cone, vertex line z_1=z_4=0) | 146 | 0 | 0 |
| C[z_(0, 2, 4, 5)]_d (cone, vertex line z_1=z_3=0) | 146 | 0 | 0 |
| C[z_(0, 3, 4, 5)]_d (cone, vertex line z_1=z_2=0) | 161 | 0 | 0 |
| C[z_(1, 2, 3, 4)]_d (cone, vertex line z_0=z_5=0) | 146 | 0 | 0 |
| C[z_(1, 2, 3, 5)]_d (cone, vertex line z_0=z_4=0) | 146 | 0 | 0 |
| C[z_(1, 2, 4, 5)]_d (cone, vertex line z_0=z_3=0) | 146 | 0 | 0 |
| C[z_(1, 3, 4, 5)]_d (cone, vertex line z_0=z_2=0) | 161 | 0 | 0 |
| C[z_(2, 3, 4, 5)]_d (cone, vertex line z_0=z_1=0) | 161 | 0 | 0 |
| bidegree (0,10) monomials | 66 | 0 | 0 |
| bidegree (10,0) monomials | 66 | 0 | 0 |
| bidegree (2,8) monomials | 270 | 0 | 0 |
| bidegree (4,6) monomials | 420 | 0 | 0 |
| bidegree (6,4) monomials | 420 | 0 | 0 |
| bidegree (8,2) monomials | 270 | 0 | 0 |
| delta * S_4^odd | 60 | 0 | 0 |

Result: the only intersections are F * C[Q]_2 (21) and F * {Q_a Q_b} (9) -- spurious F-multiples.  No pullback of a quintic
section of the cubic, no cone with vertex a plane or a line, no bigraded (X-degree, Y-degree) piece, no delta * quartic,
no A * quartic, and no product of two polars (mod J) passes E1 at degree 10.

### Hessian form of stage 1 (verified exactly, both parities)
For D = sum_l w_l F_l and N0 = A w, one has D div N0 - N0(D) == -A (w^T Hess(F) w) (mod J(F)), so the stage-1 condition is
exactly  w^T Hess(F) w in K_{2d-6}^even = (J(F) : A)_{2d-6}.  For d <= 7 this is w^T Hess(F) w in J(F)_{2d-6} (independent of A);
for d = 8, 9, 10 it is membership in the explicit kernels K_10, K_12, K_14 (codimensions 468, 208, 60).

### Batch 3 (item 2, exact closure of Jacobian-type degrees): stage-1 quadratic forms (w2_closure.py, w2_points.py)
At degrees with K_d = J_d every E1-admissible D is sum_l w_l F_l and E2 implies q(w) := [D div N0 - N0(D)] = 0 in R_{2d}^odd
(N0 = A w), a system of rJ quadratic forms in the w-coordinates, invariant under the spurious translations D -> D + F G.
Full column rank of a degree-k Macaulay matrix (exact mod p) certifies that the forms have no common projective zero over C,
i.e. NO denominator of that degree/parity admits a certificate.  msolve (F4 mod p, random affine chart l(v)=1) gives
the same conclusion when it returns the Groebner basis [1].

| d | parity | n_w | spurious | n_q | rJ | deg-2 rank | higher Macaulay ranks | closed at k |
|---|---|---|---|---|---|---|---|---|
| 5 | even | 3 | 0 | 3 | 744 | 6/6 |  | 2 |
| 5 | odd | 3 | 0 | 3 | 744 | 6/6 |  | 2 |
| 6 | even | 18 | 1 | 17 | 866 | 135/153 | deg3: 967/969 deg4: 4845/4845 | 4 |
| 6 | odd | 18 | 0 | 18 | 866 | 153/171 | deg3: 1138/1140 deg4: 5985/5985 | 4 |

msolve (exact, mod p = 2147483629):

| d | parity | char | system | mode | verdict | time |
|---|---|---|---|---|---|---|
| 6 | even | 2147483629 | 135 forms in 17 vars | affine chart l(v)=1 | no solution in the random affine chart | 0.0 s |
| 6 | odd | 2147483629 | 153 forms in 18 vars | affine chart l(v)=1 | no solution in the random affine chart | 0.0 s |
| 6 | odd | 2147483629 | 153 forms in 18 vars | homogeneous, -g 2 | EMPTY projective zero set (pure power of every variable in the leading ideal) | 0.1 s |

Conclusions so far (rigorous): d = 5 (both parities, rank-2 certificate) and d = 6 (both parities: Macaulay-4 certificates
4845/4845 and 5985/5985, exact mod p hence over Q; confirmed independently by msolve's homogeneous reduced Groebner basis,
whose leading ideal contains v_i^2 for every variable) admit NO certificate for A.  With K_d = J_d for d <= 9 this proves:
every certificate for this target has deg D >= 7 (both parities).  d = 7, 8 msolve/Macaulay running (job closure678).

### Batch 4 (item 2/3, numerical joint search, remote jobs jvp-d678 / d10; running): residual landscape so far
rel = ||residual||^2 / ||B_D(N0)||^2 at the converged point (0 would be a certificate).  Spurious directions (F-multiples,
Koszul gauge) are projected out; |u| = 1.  The remote box is shared with ~40 other jobs (load 300-670 on 128 cores), so
restart counts are far below the planned thousands; the exact closures of Batch 3 supersede the d <= 6 searches.

| d | parity | objective | relative residual statistics | cost |
|---|---|---|---|---|
| 6 | even | stage 1 (Hessian/Jacobian-quotient, LM) | n=361 min=1.890e-02 q10=2.303e-02 median=3.463e-02 max=4.660e-02 | 27.0 s/restart |
