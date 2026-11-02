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
* Degree 7: block L-ranks 1635/1602/1635/1464/1287/1464 (exact over Q = mod 65521), constraint ranks 291/318/291/470/688/470,
  combined rank 1779 on 7260 quadratic monomials — all identical to the thread's values at p = 101; [projective emptiness replay
  (charts + square reduction + homogeneous GB) at p = 65521 pending].
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
  served target space of f(U) (exact rational nullspace, 4923 x 832) has zero projection onto the K3 block modulo J,
  and so do the hyperplanes U0, V0. So X ∩ V(f(U)) (the join of the two copies of C) carries no K3-sector class at all,
  as Chapter 3A predicts (it supports the joining sector).
* Degree 7: the single-target relaxation is far too weak for A_K3 (rank 935 on 7260 quadratic monomials, zero set of
  tangent dimension 84 along the reducible family f(U)·ℓ, no squares). Under the dichotomy with ρ = 1 a denominator
  serving A_K3 serves every K3-block class, so all 19 classes can be imposed with the same Koszul vector; this 19-target
  relaxation (zero-class check passed) does not help either: every K3-block class imposes the same constraints (stacked rank 1092 of 7260, nothing forced), so degree 7 for the K3 class is open by relaxation methods; the reducible family f(U)·ℓ is nevertheless excluded by the sum principle (π_{Λ²T} of each factor vanishes).
  Exact fixed-D: f(U)·U0 and f(U)·V0 fail for A_K3 (377/378).
* Degree 8: the reducible family f(U)·h (h = U0², U0V0, V0²) lies in the K3 relaxation too (mod 65521), so the
  relaxation method cannot decide degree 8 for the K3 class either; f(U)·U0² fails exactly for both targets.
* Degree 7 (K3 class): open by relaxation methods (see above); the reducible family f(U)·ℓ and all tested structured members fail.
Consequence, conditional on ρ(S_2) = 1 (via the dichotomy: End_Hdg(T) = Q makes Λ²T irreducible, so a support either
carries no K3 (3,1)-class or all of them): no denominator whose divisor X ∩ V(D) has only components of degree <= 6 can serve A; at degree 8 only irreducible
reduced supports remain, or supports of the form (irreducible degree-7 component carrying the K3 class) ∪ (hyperplane section).

**2b. Block cones are excluded in every degree (Hodge-theoretic).** For D = g(U) (or g(V)), any plane curve g, the
pullback of X ∩ V(D) under the book's Φ: S₆ × S₆ ⇢ X is the product divisor C_Γ × S₆ (plus the joining locus S₆ × C);
as H¹(S₆) = 0, its Gysin image meets the K3 block only in (algebraic class)⊗V₂ ⊂ NS_prim ⊗ V₂, so it carries no Λ²T
class. Hence no block cone serves A unless η ∈ NS_prim⊗C. This is exactly what the degree-8 relaxation restricted to
character-0 octics leaves over: after square elimination only the 30 coefficients of cones g(U), g ∈ J(f)₈, survive
(6 residual quadrics), so the relaxation itself could never have excluded them.

**3. Screens (mod p only, no exclusion claimed).** For the K3-sector target the E1 kernel is NOT Jacobian from degree 8 on:
[A_K3] has a 90-dimensional kernel modulo J at degree 8, described completely in the (U,V) frame: K4⊗K4 (81 dimensions,
V-degree 4; K4 = ker(η: R_{f,4} → R_{f,10}), η = [z0²z1²z2²], dim 9) plus the 9-dimensional block
{ℓ(U)ℓ'(V)((U0U1U2)² + (V0V1V2)²)} (V-degree 1), the latter consisting of reducible supports (closed form verified: E1 holds for these, fails for the
sign-flipped η(U)−η(V) control). Eight elements of K4⊗K4 (the pairing Σ_k g_k(U)g_k(V) over the exact K4 basis, an alternating pairing, a rank-2
tensor, the rank-1 control, random symmetric tensors and random combinations) all fail the full E2 system for A_K3 (μ6-restricted systems 13572×2430, rank 2206, mod 2^31−1); a random element of the 9-dimensional character-1 block (reducible supports) fails as well. Random
μ6-homogeneous members of J_7 and J_8 (one per character) fail exactly over Q for both targets. [structured octic/decic sweep for A
pending remote compute].

Credits: the degree bound, normalization and packet are `1337isnot1337`'s; the sector decomposition is the book's
(Ch. 3A); the sum principle is the prescribed-support theorem (Ch. 1) plus additivity of Gysin images over components.
