# cycle-q2 — CONSTRUCTION — part 2: the Γ-tower, averaging-ideal test, equivariant Hodge numbers (lens iii / ii)

Script `w2c_gamma.py` (local, ~3 min; output `w2c_gamma.json`). Conventions: an element is (p, w) acting by
z_i ↦ ζ^{w_i} z_{p(i)}, weights mod scalars; residue character of PΩ/F^{q+1} is wt(P) + Σw − (q+1)χ_F.
Equivariant Hodge numbers from the equivariant Hilbert series of the Jacobian ring,
Π_i (1 − (c_h/λ_i) t⁵)/(1 − λ_i t) (the partials transform by c_h/λ_i), and tr(h | H^{5−k,k−1}_prim) = det(h)·c_h^{−k}·tr(h | R_{6k−6}).
Validation: identity → (1,426,1751,426,1); G-invariant (0,6,26,6,0); K-invariant (0,18,78,18,0); ⟨τ⟩ (0,208,866,208,0); G-sector of χ_A (0,6,26,6,0) — all equal to dossier A.2.

## Facts computed
- |Γ| = 216 (72 diagonal); R-translates: A = (2,2,3,2,2,1), RA = (3,2,2,1,2,2), R²A = (2,3,2,2,1,2).
- Kernels K = ker χ_A, K₁ = ker χ_{RA}, K₂ = ker χ_{R²A}: each of order 24, pairwise distinct; χ_A·χ_{RA}·χ_{R²A} = 1 on G; K_∩ = K∩K₁∩K₂ has order 8 and equals
  {0, (0,0,3,3,3,0), (0,3,0,0,3,0), (0,3,3,3,0,0), (3,0,0,3,0,0), (3,0,3,0,3,0), (3,3,0,3,3,0), (3,3,3,0,0,0)} = P × ⟨τ⟩ (paired signs and Y-parity).
- Fibrewise group of the (c1,c2)-fibration (elements fixing c1 = r1/r0, c2 = r2/r0, i.e. w3−w0 ≡ w4−w1 ≡ w5−w2): G_fib of order 24, G_fib ∩ K = K_∩ (order 8). Elements of G_fib∖K, e.g. (1,0,2,5,4,0) (χ_A = ω) and (0,1,2,2,3,4) (χ_A = ω²), act on v = Y0/X0 by ζ^{±2}, i.e. on t = v² by ω^{∓1}: these are the elements that see the cube root.
- Γ-representation of the target: ρ = Ind_G^Γ χ_A is 3-dimensional irreducible (the three characters are distinct); its multiplicity in H^{p,q}_prim is (0,6,26,6,0) (isotypic part = χ_A ⊕ χ_{RA} ⊕ χ_{R²A} sectors).
- Subgroup lattice of Γ: 226 subgroups. Averaging-ideal test (ρ^H ≠ 0 ⇔ H∩G ⊆ one of K, K₁, K₂; for H ⊄ G this forces H∩G ⊆ K_∩): 154 subgroups reach ρ; orders 1,2,3,4,6,8,12,24 only; **maximal order 24**: the three abelian kernels K, K₁, K₂ (H⁴-invariants (0,18,78,18,0)) and nine nonabelian H' = K_∩ ⋊ ⟨gR^{±1}⟩ ≅ (Z/2)³ ⋊ Z/3 with H⁴(X)^{H'} = (0,20,78,20,0) (three of them) or (0,19,78,19,0) (six). 90 of the 154 are non-diagonal (orders 3, 6, 12, 24). No subgroup of order 36, 72, 108, 216 reaches ρ.
- Non-diagonal H' of order 12 have H⁴-invariants (0,41,163,41,0) or (1,44,163,44,1) (the latter contain the (4,0) class, so their quotients have h^{4,0} = 1).

## What this says about mechanisms (ii)/(iii)
1. The nonabelian quotients X/H' (H' = K_∩ ⋊ Z/3) are new objects (never discussed in the thread; the five-point and lattice-width exclusions concern H ⊆ K only). They reach the target's representation: a certificate for A + RA + R²A projects to one for A (dossier §1.5), and A + RA + R²A is H'-invariant.
2. They are still not uniruled: benediamond's Reid–Tai computation covers every H ⊆ Γ (all ages ≥ 1, no quasi-reflections), so every X/H has canonical singularities and torsion canonical class — κ = 0.
3. Their visible surface fibration is the descended (c1,c2)-fibration: R acts on the base by (c1,c2) ↦ (c2/c1, 1/c1) (order 3, three fixed points (ω^k, ω^{2k})), so the generic fibre of X/H' → (c1,c2)/R is S_c/K_∩ = S3 (K_∩ = P×⟨τ⟩ acts fibrewise: P trivially on the cubic, τ as v ↦ −v, giving the t = v² model) with p_g = 2 — unchanged. Any other H'-invariant base map must be built from R-orbit sums (non-monomial); along a generic R-invariant base curve the Weierstrass coefficients s4, s6 have r-degrees 12 and 36, so those elliptic surfaces have p_g ≫ 0 (the v-line is special only because its coefficients are weighted-homogeneous in T = Rv⁶).
4. Galois closures: X → X/Γ is Galois with group Γ, and the tower above is its complete lattice of intermediate quotients; the covariant map of the cubic fibration (C_b → Jac(C_b), degree 9) is multiplication by 3 on the generic fibre, so its intermediate quotients are the 3-isogenies C_b → C_b/⟨Q⟩; but no element of G acts fibrewise on the cubic fibration over the full base (v,c1,c2) except P (checked: fibrewise ⇒ w ∈ P·scalars), so there is no rational 3-torsion translation to quotient by and no "intermediate" model between the cubic and its Jacobian that is defined over the base field and keeps the character. The cube root enters through the base (t = v²), not through the fibre.

Conclusion for (ii)/(iii): the averaging-ideal test is now complete for the whole Γ-tower; every quotient reaching ρ has κ = 0 and its only visible surface fibration has p_g = 2 fibres. No Fano or p_g = 0 companion exists in this tower.
