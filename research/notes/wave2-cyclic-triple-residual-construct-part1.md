# Part 1 — Covers, Galois closures and the averaging ideal: where every quotient mechanism fails for V_1

Receipts: /tmp/hodge/results/cyclic-triple-residual/wave2-construct/c2_abstract.{py,json} (exact group theory, 128
subgroups), zt_smooth.sing (Singular smoothness checks, remote job construct-ctr-ztsmooth, pending at the time of writing).

## 1.1 The class and its representation
F = Z^2 + tZC + H, Z = X2^3, t = 5/7, C = X0^3 + X0Y1Y2 + X1Y0Y1, H = X0^6+X1^6+Y0^6+Y1^6+Y2^6 + t(X0^3X1Y0Y1 +
X0^2Y0Y1Y2^2 + X0X1Y1^3Y2).  Residue class ξ = [X0^3Y0^3 Ω/F^2]: μ3-character ζ3 (X2 ↦ ζX2), τ-invariant, g3-character ζ3;
stabiliser in the diagonal group G (order 18) is ⟨s⟩, s = (1,5,5,2,2,0), all elements senior (dossier §1.4, verified by the
referee).  Note (new bookkeeping): the dossier's V_0, V_1, V_2 are the τ-INVARIANT parts of the μ3-eigenspaces
(32 + 107 + 69 = 208 = h^{3,1}(X)^τ); the full ζ3-eigenspace also contains the τ-anti-invariant (4,0)-class Ω/F.
The exact G-sector S_A = (41,119,26) of ξ has no (4,0)-part, so the GHC predicts S_A ⊕ S̄_A ⊂ N^1 — the lane is
"expected solvable", and all obstructions below are obstructions to specific MECHANISMS.

## 1.2 Natural Galois closure Z̃ and its full symmetry group (order 108)
Z̃ = {X2^3 + X̄2^3 + tC = 0, X2^3 X̄2^3 = H} ⊂ P^6 = X ×_{P^4} Y, Y = {u^3 = H} ⊂ P(1^5,2) (u = X2X̄2), a (3,6) c.i.,
K_Z̃ = O(2).  Symmetries: d (X̄2 ↦ ζX̄2, deck of Z̃ → X), e (X2 ↦ ζX2, lift of μ3), σ (X2 ↔ X̄2), r (lift of τ), g̃3 (lift of
g3, fixing X̄2): G108 = ((Z/3)^2 ⋊ Z/2) × Z/2 × Z/3.  ξ pulls back into the 2-dim irrep ρ_full = ρ ⊗ (r ↦ 1, g̃3 ↦ ζ3),
ρ = span(v10: d ↦ 1, e ↦ ζ; v01: d ↦ ζ, e ↦ 1), σ swapping.  Pseudo-reflections of G108 on Z̃ (6-dim fixed space on C^7):
d^{±1}, e^{±1} (hyperplanes X̄2 = 0, X2 = 0, cyclic-cover ramification, stabiliser order 3) and the three involutions
d^k e^{-k} σ (hyperplanes X2 = ζ^k X̄2, order 2) — nothing else (r, g̃3 and their products fix ≤ 4 dimensions).
Riemann–Hurwitz for the quotient map gives, in O(1)-degrees,
   deg K_{Z̃/K} = 2 − 2[d ∈ K] − 2[e ∈ K] − #{k : d^k e^{−k} σ ∈ K},
reproducing K_X = 0, K_W = O(−2) (K = ⟨d,e⟩), K_Y = O(−1) (K = ⟨σ, de^{−1}⟩ ≅ S3), K_{P^4} = O(−5), K_{Z̃/σ} = O(1).

## 1.3 Averaging-ideal test (book §3.1) over ALL 128 subgroups
c2_abstract.py enumerates every subgroup K ≤ G108, computes dim ρ_full^K and deg K_{Z̃/K}.  Result: 38 subgroups serve
the class (ρ_full^K ≠ 0); their canonical degrees are 0, 1 or 2, never negative; the serving quotients of degree 0 are
X, X' = Z̃/⟨e⟩ and their τ-, s^2-, s-quotients (K = ⟨d⟩ × ⟨r⟩ × ⟨g̃3 e^{−1}⟩ and subgroups: Z̃/K = X/⟨s⟩ etc.).  All
quotients with negative degree (P^4, W, Y, Z̃/⟨σ,r⟩-type Fano-ish quotients, …, 13 subgroups) have ρ_full^K = 0.
Since quotient singularities are klt and a klt variety with pseudo-effective K is not uniruled (BDPP), no uniruled
intermediate quotient of Z̃ (even using τ, g3) serves ξ.  This is the complete version of dossier §1.7 / book §4.7, now
including the "twisted" elements (products of deck, exchange, τ and g3 lifts): the only elements fixing a vector of
ρ_full ⊗ (class line) are in ⟨d, r, g̃3 e^{−1}⟩ (quotient X/⟨s⟩, κ = 0) or are σ-conjugates (quotient degree ≥ 1).

## 1.4 Why enlarging the cover by abelian extensions cannot help (lemma)
Let Z' → Z̃ be any further Galois cover whose group Γ' is an extension of G108 by an abelian deck group N acting trivially
on classes pulled back from X.  The class ξ pulls back to ρ_full ⊗ (trivial N-character).  For k = n·g (n ∈ N, g ∈ G108)
to fix a vector of this component, g must fix it; so the stabilisers are N·(stabilisers in G108), and Z'/K → Z̃/(K∩G108)
is a cover of a quotient with deg K ≥ 0, hence not uniruled.  Twisted elements such as g̃3 e^{−1} (= s^2 on X) are already
in G108.  Conclusion: the ζ3-sector cannot be reached by ANY tower of cyclic covers over the natural structure
(adjoining roots of C, H, the discriminant, or of products thereof); a new mechanism must introduce a non-normalising
exchange with large fixed divisor (Part 2 shows none exists in §2.7 form) or a genuinely different geometry (Parts 4–5).

## 1.5 The §4.6 descent criterion with this closure
Hypotheses: p = Z̃ → X, d, g = μ3, e = σdσ ✓, r, σ commuting ✓.  Missing: divisor support of H^4 of Z̃/⟨r,σ⟩.  That
quotient has deg K = 1; its (de = ζ)-eigenspace contains the multiplicity space M_ρ ≅ (full) ζ3-eigenspace of X, which
contains the (4,0)-class, so H^4(Z̃/⟨r,σ⟩) is NOT entirely divisor-supported and only the r-invariant part of M_ρ
(Hodge numbers (107,341,69)) would be needed — but no sector theorem applies to Z̃/⟨r,σ⟩ (K ample-ish, fibres of all
its natural fibrations have p_g ≥ 9: e.g. over P^2_Y the fibres are (3,6) c.i. surfaces in P^4 with K = O(4)).  The
criterion's output would in any case be H^4(X)^τ ⊂ N^1 (the whole τ-sector), a statement expected true (GHC) and
exactly equivalent to the lane's difficulty.

## 1.6 Galois closure of X ⇢ P^2_X × P^2_Y (the cubic-cover structure of an even sextic)
F(sX, tY) = s^6P6 + s^4t^2P4 + s^2t^4P2 + t^6P0 (cubic in u = t^2/s^2); Galois group Z/2 ≀ S3 (order 48), extended by μ3
acting on P^2_X.  The sign changes ε_j are reflections on the closure (fixed divisors over {P6 = 0}×P^2 etc.), and the
uniruled quotients are P^2×P^2 and P(1,1,3)×P^2 only: the double cover branched along the (12,12) discriminant has
K = O(3,3); X/τ has κ = 0; Z̃_3 (S3-closure) K = π^*O(3,3).  Book §4.7's negative for V_std applies verbatim; adding μ3
produces no new uniruled quotient (the twisted quotients by c·μ3, c a 3-cycle, have no fixed divisor).

## 1.7 Other natural finite maps
Projections from points (generic S5/S6), the map (Z : C : H) (fibres: (3,6) c.i. threefolds, K = O(3)), the P^1×P^1 and
P^2 linear fibrations (Part 3): none has a uniruled Galois-closure quotient serving ξ, by the same canonical-degree
bookkeeping (covers of X branched along divisors have K > 0; X is simply connected so étale covers do not exist).
