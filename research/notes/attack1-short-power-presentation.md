# Attack 1 — `short-power-presentation` (q = 1, A = X0^3 Y0^3, F = eight sixth powers)

**Status: NO WITNESS FOUND (work in progress; resumed after a crash of the first run).**
No certificate, no Lean claim. Engine v2 (post sign-fix). Raw log:
`/tmp/hodge/results/short-power-presentation/log.jsonl`; E1-kernel table `e1kernel.jsonl`; scripts in the same dir.

## 0. Attack plan (written 2026-11-01, resume)
F = eX + eY + c((L+M)^6 + (L-M)^6), c = 5/7^7, L = X0+2X1+3X2, M = Y0+3Y1+5Y2 (verified exactly by exp0.py).
So F = sum of 8 sixth powers of linear forms l_1..l_8 = X0,X1,X2,Y0,Y1,Y2, g(L+M), g(L-M) (g^6 = c): X is the
codimension-2 linear section of the Fermat sextic sixfold V ⊂ P^7 by P^5 = {z7+z8 = 2gL, z7-z8 = 2gM}.
Steps:
 1. E1 kernel of [A] on the Jacobian ring, both Y-parities of D, d = 1..8 (first run did even d ≤ 6: all zero mod J).
 2. If kernel = 0 up to d0: every denominator of degree ≤ d0 is a polar divisor D = w(F); the E2 condition mod J is
    the quadratic system  w^T Hess(F) w ∈ J  (Koszul-syzygy argument, see §4). Solve it exactly for d = 5, 6
    (18 unknowns) with Singular/msolve → nonexistence theorem for deg D ≤ 6, then attempt d = 7.
 3. Structured denominators from the eight-power geometry: subset sums of the 8 powers (mixed parity handled by
    D·τ(D) or a parity-free builder), partial sums of F, E_X F, E_Y F, hyperplane-pencil denominators
    L^{6a} P(eX/L^6, M/L) (benediamond's (e,m) base), cones over the special curves/surfaces of X.
 4. Geometric (2,0)-budget analysis of the supports (as in attack1-dense-rational §1.3–1.4) for this F.
