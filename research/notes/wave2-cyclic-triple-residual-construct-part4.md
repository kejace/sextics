# Part 4 — The one mechanism specific to V_1 that survives: isotrivial j = 0 elliptic exceptional surfaces

## 4.1 Principle (new formulation, W-language)
Write π: X → W, branch divisor B = {Z = 0} ≅ {H = 0} ⊂ P^4.  Let T_W = W ∩ {d = 0} be a divisor (d a weighted form in
(X0,X1,Z,Y), i.e. D = d(X0,X1,X2^3,Y) is μ3-invariant) and T = π^{-1}(T_W) = X ∩ V(D).  Suppose the surface T_W ∩ B has a
curve Γ of ordinary triple points (transversally c_3(n1,n2) = 0 with three distinct lines).  Then T has, along Γ, the
transversal singularity x^3 = c_3(n1, n2): the cone over the plane cubic {x^3 = c_3} — a simple elliptic Ẽ6 singularity
whose exceptional curve is a j = 0 elliptic curve on which μ3 (x ↦ ζx) acts as the CM automorphism.  The resolution T̃
acquires an isotrivial j = 0 elliptic surface E → Γ, and
   H^{2,0}(E) ⊇ (H^{1,0}(Γ') ⊗ H^{1,0}(E_0))^μ,  μ3 acting on H^{1,0}(E_0) by ζ3 = the residue character of ξ.
So the Gysin image of H^2(T̃) lands partly in the ζ3-sector V_1 — this is the only known way a divisor on X acquires
(2,0)-forms with the twisted character (the other components of T̃ pulled back from W contribute to V_0 only; smooth or
rationally singular components contribute no (3,1)-class at all, book §3.1).  Variants: a quadruple-point curve gives
x^3 = c_4(n1,n2), exceptional curve a Picard genus-3 curve with μ3-characters (2,1) on H^{1,0}; a cusp-type point
x^3 = n1^3 + n2^k etc.  In X-language (attack-1 §4.2) the triple-point case is "D vanishes to order 3 along a plane
Π ⊂ {X2 = 0}" with leading cubic c_3(n1,n2) + a X2^3 — automatic for μ3-invariant D.  The two plane types:
  (A) Π ⊂ {X2 = 0}: Γ ⊂ B, exceptional surface j = 0 isotrivial, forms of character ζ3 (serves V_1 directly);
  (B) Π ∋ e_{X2}: Γ = X ∩ Π is a μ3-stable plane sextic, the fibre cubic is untwisted, μ3 acts along Γ (the ζ3-part of
      H^{1,0}(Γ) ⊗ fibre form); also feeds V_1.
In the Fermat limit t = 0 the join D0 = X0^6+X1^6+X2^6 is NOT of this type (order 6 along P^2_Y, transversal cone over
C_X); the deformation of the join dies (attack-1 §4.1).  The elliptic mechanism is the one compatible with t ≠ 0.

## 4.2 Dimension constraints
Each triple-point curve of genus g contributes ≤ g (≈ g − 1 for a non-product twist) forms; a coordinate plane gives a
plane sextic Γ (g = 10), so ≤ 10 per plane.  Lower bound for the needed number: the (3,1)-dimension of the smallest
Q-Hodge substructure containing ξ — 1 ≤ h ≤ 41 (exact G-sector), 67 with the conjugate sector if S_A ⊕ S̄_A is simple
(dossier E4, unproven).  Counting at degree d (order-3 vanishing along a plane costs Σ_{k≤2} C(k+2,2)·C(d−k+2,2)
linear conditions: 406 at d = 9, 721 at d = 12) against the E1-admissible dimension (d = 9: 521 μ3-invariant, 876 all
sectors; d = 12: ≈ 5762): at most one or two planes at degree 9, about seven at degree 12.  So a degree ≤ 10 certificate
of this type can carry at most ~20 twisted forms; if h ≥ 41 the mechanism needs degree ≥ 12 and several planes/curves.

## 4.3 Computation at degree 9 (c4_planes_d9.py; results in wave2-construct/c4_planes_d9_dims.json and the log)
For each of the 20 coordinate planes and each μ3-invariant G-sector c0x (x = 0..5; exact E1 bases of the search lens,
dims 96, 82, 82, 95, 83, 83) the exact dimensions of the order ≥ 3 and ≥ 4 sub-families were computed, with the
j = 0 diagnostics (X2^3-coefficient and c_3 both nonzero) for type-(A) planes, and random small-integer members of the
largest families were tested with the sector solver mod p.  See §4.4 for the outcome (filled when the run finished).

## 4.4 Outcome
(filled below after the run)
Degree 9, μ3-invariant sectors c00–c05 (E1-admissible dims 96, 82, 82, 95, 83, 83): 120 (plane, sector) pairs admit
order ≥ 3 (every pair), with family dimensions 29–66 (order ≥ 4: 10–43) — far larger than the naive count because the
bulk is J-type (D = G(F), tangency divisors of vector fields, 62–77 per sector) and because sector constraints already
kill most low-normal-degree monomials.  For the planes inside the branch locus {X2 = 0} the random members have a
genuine j = 0 leading cubic (both the X2^3-coefficient and c_3(n1,n2) nonzero).  Mod-p sector tests (two primes,
systems 3383–3539 × 1763–1865) of random small-integer members:
  {X0,X1,X2} sector 03 (family dim 62): rank 1796 / aug 1797 — inconsistent;
  {X2,Y0,Y1} sector 00 (dim 56): 1792 / 1793 — inconsistent;
  {X0,X1,X2} sector 01 (dim 55): 1702 / 1703 — inconsistent;
  {X1,X2,Y0} sector 03 (dim 54): 1796 / 1797 — inconsistent.
(All four: the augmented rank exceeds the rank by exactly one at both primes, the same signature as every degree-8
member; receipts in the log, step C4-test.)  Consistent with Part 6: a degree-9 divisor with one triple-point curve
carries ≤ 10 twisted forms against the 67 required generically, so random members cannot work and a parametric search
in these families is not worth running before the Hodge-theoretic bound is circumvented.
