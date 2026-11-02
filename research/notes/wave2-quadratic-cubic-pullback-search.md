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
(none yet)
