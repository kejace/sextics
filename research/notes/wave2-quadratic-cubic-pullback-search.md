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

Conclusions: every E1-admissible denominator of degree <= 9 (both parities) lies in J(F) = (F_0..F_5); the first
non-Jacobian degree is d = 10 with exactly the Gorenstein-duality numbers (294 Y-even, 261 Y-odd).  d = 11, 12 pending.

