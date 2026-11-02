# Draft discussion comment for #115 — polarized-k3-remainder (wave-2 exclusion lens; NOT POSTED)

Status: no certificate; no claim of nonexistence. Independent re-verification of the degree bound, plus a new exact
exclusion for the K3-sector class that repairs a gap in the "reducible denominators" reasoning. Everything below was
computed with an independent implementation ((U,V) frame, python-flint), primes 65521 / 2^31-1 / 2^31-19, exact
rationals where stated. Reproduction: `/tmp/hodge/results/polarized-k3-remainder/wave2-excl/` (README, scripts, JSON receipts).

**1. Re-verification of `1337isnot1337`'s bound (deg D >= 8).**
* E1 (multiplication by [A] on the Jacobian ring) injective in degrees 0..9: recomputed in the (X,Y) frame with a
  different code at p = 2^31-1 [wave 1] and p = 2^31-19 [this run]; full modular rank is a lower bound for the
  rational rank, so injectivity holds over Q/C. Hence every E1-admissible D with deg D <= 9 lies in J(F).
* Degree 6: the L-relaxation (D = sum P_i a_i, B_D(Ta) in (G)+J(G)^2+(R_ijk)) rebuilt from scratch at p = 65521 and
  exactly over Q gives the same block ranks (603/558/396/260/396/558; augmented 713/715/617/470/617/715), constraint
  rank 437, and the same two-step square elimination (30 off-diagonal, then 5 diagonal entries) => D in {0, cG}.
* Degree 7: [pending: block ranks over Q and F_65521, constraint rank 1779, projective emptiness].
* Degree 8, D = f(U)U0^2: exact rank gap 756/757 in the residue-0 block reproduced; the other two blocks are also
  inconsistent (828/829, 516/517).

**2. New: the K3-sector class alone has no denominator of degree <= 6 (exact).**
Let A_K3 = [(X0+Y0)^2(X1+Y1)^2(X2+Y2)^2 - (X0-Y0)^2(X1-Y1)^2(X2-Y2)^2]/64, the K3-sector projection of A (the book's
P_0 A-hat). Served target spaces are additive over the components of X ∩ V(D), so a *reducible* denominator can serve A
even when no factor does; to use "deg >= 8" for products one needs that no low-degree support carries any K3 class.
* [A_K3] is injective on the Jacobian ring in degrees 0..7 (so E1-admissible D of degree <= 7 lie in J(F); degrees <= 4 impossible).
* Degree 5: the relaxation forces all parameters to vanish (exact).
* Degree 6: the relaxation leaves exactly the gauge class of D = f(U) (residual quadrics (c0-c2)^2, (c0-c2)(c1-c2),
  (c1-c2)^2); the exact Koszul-normalized system for D = f(U) is inconsistent (rank 136 vs 137). Moreover the full
  served space of f(U), of U0, of V0 and of tangent hyperplane sections through {V=0} has zero projection onto the K3
  block (a modular computation that is rigorous for vanishing). So X ∩ V(f(U)) (the join of the two copies of C) only
  supports the joining sector, as Chapter 3A predicts.
* Degree 7: [pending].
Consequence, conditional on ρ(S_2) = 1 (via the dichotomy: End_Hdg(T) = Q makes Λ²T irreducible, so a support either
carries no K3 (3,1)-class or all of them): no denominator whose divisor X ∩ V(D) has only components of degree <= 6·[7]
can serve A; at degree 8 only irreducible supports remain. [adjust once degree 7 is done]

**3. Screens (mod p only, no exclusion claimed):** [structured octics / decics sweep, non-Jacobian E1-kernel elements of
A_K3 at degree 8 — pending remote capacity].

Credits: the degree bound, normalization and packet are `1337isnot1337`'s; the sector decomposition is the book's
(Ch. 3A); the sum principle is the prescribed-support theorem (Ch. 1) plus additivity of Gysin images over components.
