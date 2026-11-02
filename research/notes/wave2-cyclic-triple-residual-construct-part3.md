# Part 3 — Surface-fibration companions of quotients X/K (book §3.2 theorem), K inside the stabiliser of the class

Receipts: /tmp/hodge/results/cyclic-triple-residual/wave2-construct/c3_fibration_pg.{py,json}.

## 3.1 Setting
The §3.2 surface-fibration theorem needs a fourfold M with a dominant rational map to a surface whose generic fibre has
p_g = q = 0; then ALL of H^4(M) is divisor-supported.  To serve [A] through quotient transfer we need M = X/K (or a
quotient of a cover) with K fixing the class.  Within the diagonal group G (order 18) the stabiliser of the residue
character of A is ⟨s⟩, s = (1,5,5,2,2,0) (order 6, s^3 = τ), all of whose elements are senior (dossier §1.4), so X/K
is never uniruled; a p_g = q = 0 surface fibration is the only §3.2 route left.

## 3.2 Linear fibrations (fibres = sextic surfaces in P^3): exact p_g count of the fibre quotients
For a base P^2 given by three coordinates B (fibre coordinates: the scaling λ of B and the other three), the subgroup
K_b ⊂ ⟨s⟩ acting on the fibre is {s^k : the three base weights agree}, and p_g(S_b/K_b) = number of K_b-invariant
residues m·Ω_3/F_b, m quadratic (Griffiths; the volume form carries the weight sum).  All 20 bases: p_g ∈ {3,4,10},
minimum 3 (bases {X0,X1,X2}, {Y0,Y1,Y2}, …, K_b = ⟨τ⟩).  For P^1 × P^1 bases (two coordinate pairs): minimum 2, attained
only by ([X1:X2],[Y0:Y1]) with K_b = ⟨s⟩ (the surviving forms are the residues of λ^2 and μ^2, where λ, μ scale the pairs).
Trivial bound: p_g(S_b/K_b) ≥ 10/|K_b| ≥ 10/6 > 1 for every linear fibration, so no linear surface fibration of X/K has
p_g = 0 fibres, for any K ⊂ ⟨s⟩.

## 3.3 The P^2 × P^2 (passive-ratio) fibration of X/τ used in book §4.2
Fibres are the (6,3)-surfaces in P^2_X × P^1_{(b:r)} with p_g = 20.  An element g of G acts fibrewise iff
g_{3+i} − g_i ≡ g_3 − g_0 (mod 6) for i = 1,2: only τ does.  So the fibre quotient available here has p_g = 20 (versus
20 → 2 → 0 for the paired-sign sextics of §4.2); the sector theorem of §4.2 does not extend.

## 3.4 K3-fibration remark
W = X/μ3 ⊂ P(1,1,1,1,1,3) is fibred over P^2 (any three weight-one coordinates) by K3 surfaces {Z^2 + tZC + H = 0} ⊂
P(1,1,1,3) (double planes branched along plane sextics); the fibres of X are their μ3-covers branched along the plane
sextic {H = 0}: sextic surfaces with p_g = 10 = 1 (K3 part, V_0) + 6 + 3 (the two twisted characters).  The twisted parts
of p_g are the fibrewise shadow of V_1, V_2 and are never killed by the available fibrewise groups (⟨μ3, τ⟩ at most:
p_g(S_b/τ)_ζ = 2, from X0 μ, X1 μ).  So the K3 structure explains the companion (V_0) and quantifies why V_1 escapes.

## 3.5 Verdict
No fibration-based sector theorem (uniruled quotient, or p_g = q = 0 surface fibres) applies to any quotient of X by a
subgroup fixing the class; the margins are p_g = 2 (best case) instead of 0.
