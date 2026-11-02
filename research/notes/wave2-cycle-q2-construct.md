# Wave 2 — `cycle-q2` — lens CONSTRUCTION (resumed session, 2026-11-02)

STATUS: finished (session 2) — **no witness, no new certificate**; mechanism map complete; the one unexcluded mechanism (M*, part 5) is a two-variable Gauss–Manin problem outside the polynomial-curve library; no finite-field seed can exist for it (failure at dominance/uniruledness).

Log: `/tmp/hodge/results/cycle-q2/wave2-construct.log.jsonl` (8 lines). Scripts: `/tmp/hodge/results/cycle-q2/w2c_{gamma,cover,qraise,factor,s3}.py` (+ `w2c_gamma.json`, `w2c_cover.json`, `w2c_factor.json`). All local, light (≤ 3 min); no remote jobs were needed.

## Table of contents
- part1 `wave2-cycle-q2-construct-part1.md` — reading summary, the five sub-lenses, plan.
- part2 `wave2-cycle-q2-construct-part2.md` — (ii)/(iii) Γ-tower (216): ρ = Ind χ_A, averaging-ideal test on all 226 subgroups, equivariant Hodge numbers, the new nonabelian quotients X/(K_∩ ⋊ Z/3).
- part3 `wave2-cycle-q2-construct-part3.md` — (i) covers creating a reflection: 0 of 60 coordinate choices; non-coordinate forms excluded via Fermat-curve automorphisms.
- part4 `wave2-cycle-q2-construct-part4.md` — (iv) q-raising identity verified (sympy); all 61 degree-6 families classified: character-free ⇔ G' ⊆ K, RES-surviving ⇔ v-residues {0}; the two sets are disjoint.
- part5 `wave2-cycle-q2-construct-part5.md` — (v)/M*: explicit S3 model (s4, s6 in t³; χ = 3, p_g = 2, κ = 1), target in the ω₀ = dt∧dξ/η sector; exactly where the book's incidence recipe fails; vertical supports excluded; the remaining equation.

## Summary (20 lines)
1. Target: q = 2, A = X0²X1²X2³Y0²Y1²Y2 on the cyclic sextic; residue character ω² on G/K; G-sector (0,6,26,6,0); α^{3,1} ≠ 0 (fqlx).
2. (i) No §2.7 cover: all 60 (y;ℓ1,ℓ2) coordinate choices fail, the best by exactly one coupling monomial (the closing term); non-coordinate forms impossible (Fermat sextic curve automorphisms are monomial). §2.8 inapplicable (dossier §1.8).
3. (iii) Γ = G ⋊ ⟨R⟩ (216): ρ = Ind_G^Γ χ_A is 3-dim irreducible, multiplicity (0,6,26,6,0) in H⁴_prim. 154/226 subgroups reach ρ; max order 24: K, K₁, K₂ and nine nonabelian K_∩ ⋊ Z/3 (K_∩ = P×⟨τ⟩), H⁴-invariants (0,20,78,20,0)/(0,19,78,19,0). All have κ = 0 (Reid–Tai, benediamond) and their only visible surface fibration descends from the (c1,c2)-fibration with the same p_g = 2 fibres (R acts on the base by (c1,c2) ↦ (c2/c1, 1/c1)).
4. (ii) No element of G acts fibrewise on the cubic fibration over (v,c1,c2) except P, so no rational 3-torsion quotient between the cubic and its Jacobian; the cube root sits in the base (t = v²), detected by G_fib∖K (e.g. (0,1,2,2,3,4), χ_A = ω², t ↦ ω²t).
5. (iv) The q-raising identity A − 2w(F_s) + F_s div w = M·E1 − (F_s/2)∂_sE1, w = (Mv − F_s v̇)/2, is verified; a q=1 family certificate for (F + sM', A/M') regular at s = 0 gives a q=2 certificate with denominator D_0².
6. (iv) Of the 61 degree-6 Y-even factorizations: 5 keep the RES structure but keep the character on the fibrewise group; 5 are character-free (G' = K) but force the S3 base change (v-residues {2,4}); none does both — the obstruction is one and the same.
7. (v)/M*: the only mechanism reaching the character is the S3-transport (DiamondJM): S3: η² = ξ³ + s4(t³)ξ − s6(t³)/4 over (c1,c2), s4 = c1c2t³(…)/8 explicit, χ = 3, p_g = 2, κ = 1; the target's (3,1)-part is fibrewise in the sector of ω₀ = dt∧dξ/η.
8. Book recipe fails precisely at *dominance*: a surface dominated by polynomial curves is uniruled; S3 (κ=1), the K3 fibres (κ=0) and X/P (κ=0) are not. Hence no finite-field incidence seed exists for any auxiliary surface carrying the ω²-sector, and `PolynomialCurves`/`HoldoutCover` cannot be the Lean engine here.
9. Vertical divisors of the (c1,c2)-fibration cannot support α (Leray + injective KS ⇒ p_g = 0 for their resolutions); removing curves from fibres leaves the fibrewise (0,2)-defect; so any primitive is a non-fibrewise rational current Ψ₁dc1 + Ψ₂dc2.
10. Remaining equation: [σ] = 0 in H²_dR(B°, V_{ω²}), V_{ω²} = H¹(P¹_T, R¹_RES ⊗ L_{ω²}) (Kummer-twisted cohomology of the rational elliptic surface, rank 12, type (1,10,1)); computable from the RES's twisted Picard–Fuchs data; poles along a horizontal non-normal divisor allowed (NL search, dossier P5).
11. Nothing posted; nothing in the repo touched; no git.
