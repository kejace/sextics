# cycle-q2 — CONSTRUCTION — part 4: q-raising (lens iv) — identity verified, all 61 families classified

## The identity (verified symbolically, `w2c_qraise.py`)
Let F_s = F + sM (deg M = 6), B with A = M·B, v = n/D a q=1 field with E1(s) := B − v(F_s) + F_s div v.
With w := (M v − F_s ∂_s v)/2 one has, identically in all data,
    A − 2 w(F_s) + F_s div w  =  M·E1 − (F_s/2)·∂_s E1.
Hence a q=1 certificate E1 = F_s² Q_s for (F_s, B), rational in s and regular at s = 0 (D_0 ≠ 0, gcd(D_0,F) = 1),
gives the q=2 certificate A − 2w(F) + F div w = −½F³ ∂_sQ|_{s=0} with denominator D_0² (dossier §5 P2 / A.7, now
checked by sympy on a random two-variable instance: residual 0, and the consequence −½F_s³Q̇ checked with a symbolic Q).
This is the certificate-level form of α = −½∇_{∂s}β_s|_{s=0} (benediamond 3.9) and needs no geometry.

## Classification of the families (`w2c_factor.py`, `w2c_factor.json`)
Factorizations A = M'·B' with deg M' = 6 and M' Y-even: **61**. For each: residue character of the sibling
(F, q=1, B') on G; symmetry group G' = Sym(F + sM') = {w ∈ G : wt(M') = χ_F(w)}; whether M' is paired-sign
invariant (P-invariant ⇔ the cubic descent survives on the family); and the torus bookkeeping of the cubic's
invariants: a term of s4/s6 containing n' copies of the new coefficient and (n_d,n_e,n_f) coupling coefficients has
v-degree 2(n_d+n_e+n_f) + n'·Ydeg(M') and is admissible iff the u-weights sum to (k,k,k) mod 3; the elliptic surface
over the v-line stays the rational elliptic surface (coefficients in v⁶ = T/R) iff all admissible v-degrees are ≡ 0 mod 6.

Key identity: res-char(B') = χ_A^{res} + (χ_F − wt(M')), so on G' the sibling has exactly the target's character:
**the sibling is character-free on its own symmetry group iff G' ⊆ K.**

Results:
- 2 families keep all of G (M' G-invariant): M' = X1²Y0²Y1² (B' = X0²X2³Y2) and M' = X0X2³Y0Y2 = M1 (B' = X0X1²Y0Y1²); both keep the RES structure (v-residues {0}) and both siblings have character {1,ω,ω²} on G = G' — the obstruction is unchanged.
- 5 P-invariant families keep the RES structure (v-residues {0}): the two above plus M' = X1²X2²Y1², X0²X1²Y0², X0²X1²X2² (G' of order 24, G' ⊄ K): for all five the sibling character on the fibrewise group G'_fib is {1,ω,ω²} — the transported sibling source still carries the cube root.
- 5 P-invariant families have G' ⊆ K (G' = K, sibling character-free on X_s): M' = X1²X2²Y0², X0X2Y0Y1²Y2, X0²Y0²Y1², X0²X2²Y0², X0²X1X2Y1Y2 — for all five the v-residues are {2,4}: the deformation term puts v² into s4, s6, so the elliptic surface over the v-line becomes the base change S3 (χ = 3, p_g = 2) and the "every polynomial source has a rational primitive" theorem of the RES is lost (this is DiamondJM's second-variation setting for M' = X0X1X2²Y0Y1 — which is not P-invariant and does not even descend: it is among the 51 non-P-invariant factorizations).
- **No family is in both lists.** The two conditions are the same obstruction (invariance under v ↦ ζv of the invariants ⇔ triviality of the character on G'_fib), seen on the family side and on the source side.
- 3 siblings have trivial character on all of G: M' = X1²X2²Y0² (B' = X0²X2Y1²Y2), X0²Y0²Y1² (B' = X1²X2³Y2), X0²X1X2Y1Y2 (B' = X1X2²Y0²Y1); the other two G' = K families (M' = X0X2Y0Y1²Y2, X0²X2²Y0²) have character 2χ_A on G (same kernel). All five have v-residues {2,4}.

## Conclusion for (iv)
The q-raising identity is correct and cheap (a ten-line Lean lemma once a q=1 family certificate exists), but every
one of the 61 one-parameter monomial families reproduces the cube-root obstruction either as a character on the
sibling or as the S3 base change on the auxiliary surface. q-raising moves the problem, it does not remove it; the
sibling problems (q=1, deg D ≥ 7 by E1) are the right place for a *search*, not for a geometric construction.
