# Wave 4 — polarized-k3-remainder — route `amix-component` — Part 1: component specification

(Companion to /tmp/hodge/notes/wave4-polarized-k3-remainder-amix-component.md; all computations logged in
/tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl, data in .../wave4/.)

## 1. The target split and what the component is for

Conventions: U = X + Y, V = X − Y, F = ½(f(U) + f(V)), f(z) = z0^6 + (5/7)(z0^5 z1 + 2 z1^5 z2 + 3 z2^5 z0 + z0^2 z1^3 z2) + z1^6 + z2^6,
h(U,V) = (U, ηV) (η a primitive sixth root of unity), τ(X,Y) = (X, −Y) = block swap. A = X0X1X2Y0Y1Y2 = 2^{-6} ∏(U_i² − V_i²).

**Exact split (proven, exact over Q, wave4/amix_pieces.json):**
  A = A_K3 + A_mix,
  A_K3 = 2^{-6}(U0²U1²U2² − V0²V1²V2²)   (V-degree b ∈ {0, 6}; residue character η³ = −1: the K3 sector Λ²V_2),
  A_mix = A_2 + A_4, A_2 = −2^{-6}(U0²U1²V2² + U0²V1²U2² + V0²U1²U2²) (b = 2, character η⁵),
                     A_4 = +2^{-6}(U0²V1²V2² + V0²U1²V2² + V0²V1²U2²) (b = 4, character η⁷ = η).
In physical coordinates A_mix is the Y-odd 13-term polynomial
  A_mix = ¾·X0X1X2Y0Y1Y2 − (1/16)·Σ (the twelve monomials X_i²X_jX_kY_k-type and X_iY_j²Y_k... listed in amix_pieces.json);
A_K3 = A − A_mix (13 terms). A_4 = −swap(A_2), so A_mix is swap-antisymmetric (Y-odd), as a legal target must be.

**Non-triviality (modular, p = 2^31−1 and 2^31−19, two primes agree):** none of A, A_K3, A_mix lies in J(F)_6
(rank 18 → 19 after adjoining the target); A_2 and A_4 are individually nonzero mod J (36 → 37).

## 2. Exactly which classes the component serves (statement)

Let X = {F = 0}. By Shioda–Katsura/book §3A.4, H⁴(X)_prim = [H²(S_6)_prim ⊗ H²(S_6)_prim]^{diag μ_6} ⊕ [H¹(C)⊗H¹(C)](−1), S_6 = {u^6 = f},
and a polynomial residue of (U,V)-bidegree (6q − b, b) has h-character η^{b+3}. The **sixth-root sector** is the η^{±1} part,
i.e. b ≡ 2, 4 (mod 6); it is H_1⊗H_5 ⊕ H_5⊗H_1 (ranks 21·21·2 = 882 before the swap; the τ-fixed part has rank 441,
Hodge numbers (90, 261, 90)). At q = 1 the sector is spanned by the swap-antisymmetric combinations of the bidegree-(4,2)
and (2,4) monomials: dim R_f,4 ⊗ R_f,2 = 15·6 = 90 = h^{3,1} of the sector. A_mix is the η^{±1}-isotypic projection of A
(no computation needed: V-degrees 2 and 4 only), so the book's polynomial test P_0 Â_mix = 0 holds and **book Thm 3A.5
asserts (nonconstructively, not formalized) that HasCertificate(F, 1, A_mix) holds** for every polarized-split sextic.

The component therefore serves: every residue class of V-degree b ≡ 2, 4 mod 6 of any F = f(U) + f(V) with smooth
ternary f for which W_{23,f} = {f(z) + s² + t³ = 0} ⊂ P(1,1,1,3,2) has a dominant weighted-conic family with the seed
properties of §3 (the seed is f-specific). It does **not** serve b ≡ 0 mod 6 (K3 sector, A_K3) and is not needed for
b ≡ 1, 5 (third-root pair: triple-plane companion, already formalized for the crossed case) or b ≡ 3 (joining divisor:
explicit strict certificate D = f(V), N = (A/3)E_V, wave-2 B2).

**Consequence for this lane (proven modulo the standard support correspondence used by the book):** the h-action on the
resolution of the joining divisor {f(U) = f(V) = 0} ∩ X (a P¹-bundle over C × C, h rotating the fibre coordinate) is trivial
on H², so its Gysin image lies in the η⁰ sector. Hence **no denominator whose zero locus on X is supported on the joining
divisor (D = f(U)^a f(V)^b, or ≡ such mod F) can serve A_mix or A** — the b = 2, 4 chart identities have a nonzero right-hand
side exactly for this reason. Calibration: the degree-6 Jacobian-type D = f(U) − f(V) is E2-inconsistent for A_mix (engine,
see log `e2-fixed-D-degree6-calibration`).

## 3. The uniruled auxiliary variety, the seed, and its lift

W_{23,f} = {H = 0}, H = 7f(x,a,y) + 7S² + 7T³ = 7x⁶ + 5x⁵a + 5x²a³y + 15xy⁵ + 7a⁶ + 10a⁵y + 7y⁶ + 7S² + 7T³ ⊂ P(1,1,1,3,2),
smooth, −K = O(2); W_{23,f} ≅_bir (S_6 × E_23)/diag μ_6 with E_23 = {α² + β³ = −1} (j = 0, μ_6-action, dβ/α of character η).
Profile-A family (square incidence, 13 unknowns c0..c12, parameters k = (k0,k1,k2,S6,T4)):
  x = 4t, a = c0 + c1t + c2t², y = k0 + c3t + c4t², S = k1 + c5t + … + c9t⁵ + S6·t⁶, T = k2 + c10t + c11t² + c12t³ + T4·t⁴,
equations = coefficients of t⁰..t¹² of H(φ(t)) (all of degree 6 in c; msolve input wave4/inc_seed_*.ms).
Seed (wave 2, re-verified exactly here): k = (2,3,1,1,1), c = (1,1,2, 0,3, 4,0,0,0,4, 0,2,0) over F_5.
  * H(φ) ≡ 0 (mod 5); 13×13 incidence Jacobian has det ≡ 3 (mod 5); explicit inverse mod 5 stored;
  * freeness (polynomial-normal) matrices: 12×16 (all five coordinates) and 12×14 (x fixed) both of rank 12 mod 5, explicit
    right inverses stored; the 12×12 submatrix of the first twelve x-free columns has det ≡ 3, inverse stored
    (these are exactly the objects `seedJacobianExplicit/Right`, `seedFreeExplicit/Right` of the accepted ConicSeed.lean);
  * endpoint gradient at t = 0: (H_x,H_a,H_y,H_S,H_T) ≡ (0,2,4,2,1) ≠ 0 (mod 5);
  * 5-adic Newton lift converges to precision 5^64 with zero residual (wave4/w23_seed_exact_data.json, `lift_mod_5^64`).
"Lift" in the sense required by Lean = `ReflectionCoupledFano.square_system_generic_lift`: special-fibre solution + unit
Jacobian + injective Z → generic field ⇒ a finite separable field factor over C(k0,k1,k2,S6,T4) carrying a solution with unit
Jacobian. The data above are precisely its hypotheses (status: **computed and checked exactly over F_5 / Z_5; not yet in Lean**).
The *algebraic degree* of the exact lift (the number of geometric solutions of the square system at the seed parameters)
is being computed with msolve (local p = 65521; remote p = 32003 and random parameters) — see the status note.

Gysin image: the family is dominant (rank of ∂φ/∂(k0,k1,k2,t) at the seed = 4 = dim of the affine cone, wave 2, mod 5), so
Chapter 3's cylinder-map theorem gives a divisor on a resolution W̃ supporting all of H³(W̃_{23,f}) (rank 42, type (0,21,21,0));
quotient transfer gives [(H_1 ⊕ H_5) ⊗ H¹(E_23)]^{diag μ_6}; the translated-diagonal contraction of the two copies
(U-block and V-block, both on the *same* W_{23,f} since the two blocks share f) gives H_1⊗H_5 ⊕ H_5⊗H_1 on X.
