# Attack 1 — `eight-plane-configuration` (q = 1, A = X0^3 Y0^3)

Status: IN PROGRESS (file is rewritten at the end of the session; this is the skeleton).

## 0. Top line

- No witness found so far.
- Main structural result: the annihilator of the target class in the Jacobian ring vanishes in all
  degrees d <= 8 (even parity) and d <= 7 (odd parity): every denominator of degree <= 8 that passes
  the necessary E1 condition lies in the Jacobian ideal J(F) (polar/tangency divisors), and those
  tested are inconsistent.

(Sections 1-5 follow; see the final version.)

## Attack plan (resumed session, 2026-11-01 18:10)

Facts carried over (engine v2): A not in J(F) (D=1 fails E1); E1-kernel mod J(F) is 0 for d<=8 (both parities);
the six partials dF/dZ_i (d=5, in J) are inconsistent. Parity-graded Hilbert function of the Jacobian ring:
mult by [A] is square at d=9 (623x623 per parity) and has forced kernel 294 (even) / 261 (odd) at d=10.

1. Kernel at d=9 (both parities) -> if nonzero, test every element (fixed-D at d=9: minutes each).
2. Kernel at d=10 (both parities) -> intersect with structured families: C[Q,K,H,diag]_10 (net pullbacks
   G_5(Q0,Q1,Q2) = preimages of quintic curves under the invariant net), cones with vertex at points of X,
   products of conjugate-plane cones; test the survivors with the fixed-D solver (tens of minutes each).
3. J(F)-type denominators: structured linear vector fields (block-radial, Euler, rotations) at d=6, random
   elements of J_5, J_6; a few structured J_7/J_8 (c vanishing on the 8 planes / on Q_i=0).
4. Geometric E1 sanity sweep (Q-products, K_i, H_i, Fermat part, cones) -- predicted to fail by (fact 2).
