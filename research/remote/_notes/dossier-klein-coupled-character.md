# Research dossier: `klein-coupled-character` (OPEN)

Compiled 2026-11-01 from: `data/klein-coupled-character/{instance.rec, why.md, provenance.json, checks.json}`, the complete GitHub discussion #116 (one comment, saved at `/tmp/hodge/discussions/klein-coupled-character.md`), the full reproduction gist of that comment (29 files, downloaded to `/tmp/hodge/downloads/klein-gist/files/`), the accepted solutions of the related instances (`klein-rational`, `klein-hesse`, `holdout`, `holdout-ternary-deformation`, `cycle-three`, `mixed-biweight`), the book chapters `docs/01` (certificates), `docs/02` (§2.1, §2.4, §2.5, §2.7, §2.8, §2.9), `docs/03` (§3.1, §3.3, §3.4), `docs/03a` (§3A.1), `docs/04` (§4.1–§4.7), `docs/05`, `scripts/audit-frontier.py`, and my own computations (scripts and logs in `/tmp/hodge/scratch/kcc/`, summarized in Appendix A).

Conventions: variables `(X0,X1,X2,Y0,Y1,Y2)`; throughout I also use fqlx's short names `(x,y,z,t,r,s) = (X0,X1,X2,Y0,Y1,Y2)`. "Rigorous" means proved (by a cited argument or an exact computation I reproduced); "evidence" means a modular or floating-point computation; "claim" means stated in the thread/packet and not independently re-derived.

## 0. Status in one paragraph

The target is `q = 1`, `A = X0^3 Y0^3` on the coupled Klein sextic below. No certificate, no existential proof and no Lean proof exists. The only public prior work is one comment by `fqlx` (2026-10-26) with a large packet; it contains exactly one *rigorous* denominator exclusion (the degree-9 polynomial `Jhat`, with all its powers, see §3.2 and Appendix A.4; I reproduced it exactly over Q), two further *scoped* exclusions of particular construction families inside fqlx's own transport framework (not independently auditable in full), and a number of explicit intermediate objects. Everything geometric in the book (§4.3) serves only the σ‑invariant sector; the target lives in the ζ₇³‑eigenspace, a sector of primitive Hodge numbers (26,123,37) that no quotient, reflection, cover or Piola mechanism in the book reaches. New rigorous facts established here: the projective symmetry group is exactly Z/14 (diagonal), there are no monomial permutation symmetries and no infinitesimal linear symmetries; every legal Y‑even denominator of degree ≤ 8 must lie in the Jacobian ideal J(F) (Appendix A.3); the polar hypersurfaces `F_i = 0`, the products `Z_j F_i` and the six degree‑11 denominators `xyztrs·F_i` are excluded as supports (modular evidence, A.5; the degree‑5/6 cases are exact because the reduced right‑hand side is nonzero while the Koszul column space is visibly too small — see A.5). This dossier was compiled in two sessions; the first crashed while background computations were running, and Appendix A records exactly which of them finished. The main unexplored route is the complete fixed/joint denominator search in low degree organized by σ‑characters (sizes in §1.9 and §5), which nobody has run.

---

## 1. Exact target

### 1.1 Statement

```
q = 1
F = X0^5*X1 + 10/7*X0^3*X2*Y1*Y2 + 5/7*X0*X1^3*Y0*Y1 + X0*X2^5 + X1^5*X2 + 15/7*X1*X2^3*Y0*Y2 + Y0^6 + Y1^6 + Y2^6
A = X0^3*Y0^3
```

Certificate (formal contract, `formal/README.md`): nonzero homogeneous `D` with `gcd(D,F)=1`, numerators `N_i` homogeneous of degree `deg D + 6q − 5 = deg D + 1`, and

```
E := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^{q+1}·Q = F²·Q.
```

Degrees: `deg E = 2d + 6`, `deg Q = 2d − 6` (`d = deg D`). Normalized form (ch. 1): `E1: N(F) = A·D`, `E2: D·div N − N(D) ∈ (F)`; `d ≤ 3` would force a strict certificate (`Q = 0`).

### 1.2 Structure and sparse presentation

With `(x,y,z,t,r,s) = (X0,X1,X2,Y0,Y1,Y2)`:

```
F = K(x,y,z) + G(t,r,s) + M,
K = x^5 y + y^5 z + z^5 x                      (cyclic "Klein" sextic block)
G = t^6 + r^6 + s^6                             (Fermat block)
M = (5/7) * ( x y^3 (t r) + 2 x^3 z (r s) + 3 y z^3 (t s) )   (three mixed quartic×quadric couplings)
```

Monomial support (9 terms, exponent vectors over `(X0,X1,X2,Y0,Y1,Y2)`):

| monomial | coefficient | σ-weight (mod 7) |
|---|---|---|
| X0^5 X1 | 1 | 5+2 = 0 |
| X1^5 X2 | 1 | 10+4 = 0 |
| X0 X2^5 | 1 | 1+20 = 0 |
| X0^3 X2 Y1 Y2 | 10/7 | 3+4 = 0 |
| X0 X1^3 Y0 Y1 | 5/7 | 1+6 = 0 |
| X1 X2^3 Y0 Y2 | 15/7 | 2+12 = 0 |
| Y0^6, Y1^6, Y2^6 | 1 | 0 |

Integral model: `7F` has integer coefficients. There are no pure powers `X_i^6`; hence the three coordinate points `[1:0:0:0:0:0]`, `[0:1:0:0:0:0]`, `[0:0:1:0:0:0]` lie on `X = V(F)`.

Membership in the book's family (docs/04 §4.3): the complete σ‑invariant, Y‑even family is
`F = a x^5y + b y^5z + c z^5x + d x^2y^2z^2 + x y^3 Q0(Y) + x^3 z Q1(Y) + y z^3 Q2(Y) + g(Y)` (50 coefficients). Here `a=b=c=1`, `d=0`, `Q0 = (5/7) t r`, `Q1 = (10/7) r s`, `Q2 = (15/7) t s` (each a product of two distinct linear forms, i.e. rank‑2 quadrics), `g` = Fermat. The three mixed X‑monomials `x y^3, x^3 z, y z^3` are exactly the three σ‑invariant quartic monomials.

The coefficients `(1,2,3)·5/7` of the three couplings are deliberately distinct: with equal coefficients the equation would acquire the order‑3 monomial symmetry `X0→X1→X2→X0, Y0→Y2→Y1→Y0` (see §1.5), which is absent here.

### 1.3 Parity and degrees

`F(X,−Y) = F(X,Y)` (every mixed term has Y‑degree 2); `A` has Y‑degree 3 (odd). Numerator parities: `N_0,N_1,N_2` Y‑odd, `N_3,N_4,N_5` Y‑even when `D` is Y‑even. `Q` is Y‑odd for Y‑even `D`.

### 1.4 Smoothness (rigorous, from `checks.json`/`why.md`)

Reduce modulo 5: the mixed coefficients `10/7, 5/7, 15/7` vanish, the partials of `F mod 5` are `(X2^5, X0^5, X1^5, Y0^5, Y1^5, Y2^5)` (the Klein partial `∂K/∂x = 5x^4y + z^5 ≡ z^5`, etc.), with no common projective zero; so `X mod 5` is smooth, hence `X` is smooth over Q (nonvanishing of the discriminant mod 5). `F` is irreducible (smooth hypersurface of dimension 4). The Jacobian ideal `J(F)` is a complete intersection; `dim_C S/J(F) = 5^6 = 15625` (confirmed by Singular mod 32003, `vdim = 15625`). Also `A = X0^3 Y0^3 ∉ J(F)`: its reduction mod 5 has all exponents < 5, hence is nonzero in `F_5[Z]/(Z_i^5)`; the standard denominator argument (the 18 generators `Z_jF_i` stay independent mod 5) transfers this to Q. So the residue class of `A` is a nonzero primitive (3,1) class (`targetNonzeroInReducedJacobianRing: true`).

### 1.5 Symmetry group (rigorous)

**Diagonal projective symmetries.** Let `L(F)` be the lattice spanned by exponent differences of the support. The projective diagonal group is `Hom(Z^6_0 / L(F), C^*)`. Smith normal form of `L(F) + Z·(1,…,1)` has invariant factors `(1,1,1,1,1,84)` and `84 = 6·14`, so the group has order **14** (the audit's `projectiveDiagonalGroupOrder: 14`; my explicit enumeration lists all 14 elements). Generators:

* `σ : (X0,X1,X2,Y) ↦ (ζX0, ζ²X1, ζ⁴X2, Y)`, `ζ = e^{2πi/7}`, weights `w = (1,2,4,0,0,0) mod 7`;
* `τ : Y ↦ −Y`, projectively equal to `(X0,X1,X2) ↦ −(X0,X1,X2)`.

`G = ⟨σ⟩ × ⟨τ⟩ ≅ Z/14`. Explicit weight vectors mod 126 (last coordinate normalized to 0): `(18k, 36k, 72k, 0,0,0)` for `σ^k` and `(63,63,63,0,0,0) + σ^k` for `τσ^k`.

**Residue characters.** For a monomial numerator `Z^e` the residue character under weights `w` is `Σ_j (e_j+1) w_j`; since `1+2+4 = 7 ≡ 0`, the volume factor is trivial and the residue character of `Z^e` is simply `Σ e_j w_j mod 7` (numerator character). For `A = X0^3Y0^3`: character `3·1 = 3`, i.e. `ζ³`. Under `τ` the residue of `A` is invariant (as for every Y‑odd numerator). Table of residue characters of `A` under the 14 elements: `σ^k ↦ ζ^{3k}`, `τσ^k ↦ ζ^{3k}` (`k = 0..6`). Only `id` and `τ` fix the character.

**Non‑diagonal symmetries.** (i) Monomial (permutation × diagonal) symmetries: among the 36 block‑preserving coordinate permutations only three preserve the support: the identity and the two rotations `X0→X1→X2→X0` combined with `Y0→Y2→Y1→Y0` (resp. the inverse pair). For these, matching coefficients requires diagonal scalings with coefficient ratios `{1/2, 2/3, 3}` (resp. `{1/3, 3/2, 2}`), which are not roots of unity; since a scaling in a finite group must consist of roots of unity, no such symmetry exists. (ii) Infinitesimal linear symmetries: the space of `6×6` matrices `B` with `(BZ·∇)F = λF` is one‑dimensional (Euler only), so the linear automorphism group of `F` is finite. (iii) Block structure: Hessian symmetrizer dimension 1 (audit) ⇒ no linear block splitting even over C; `subsetEulerAtoms = []`; no coordinate reflections (`coordinateReflections = []`): every `X_i` appears to odd powers in the Klein block and every `Y_j` to odd powers in some mixed term. (iv) No paired signs `(X_i,Y_i) ↦ (−X_i,−Y_i)` (the Klein terms are odd in a single `X_i`). A complete classification of finite non‑diagonal linear automorphisms has not been done, but (i)–(iv) rule out every mechanism the book uses.

**Ages (for Kollár–Larsen junior‑element arguments).** Fixed loci of `σ^k` on `X`: the plane curve `{X=0} ∩ X = {g(Y)=0}` (normal weights `(1,2,4)k/7`, age 1 for every `k`), and the three coordinate points. Tangent weights at `[1:0:0:0:0:0]` are `(3,6,6,6)k/7` → ages `3,3,2,2,1,1` for `k=1..6`; at `[0:1:0:…]`: `(6,5,5,5)k/7` → ages `3,2,1,3,2,1`; at `[0:0:1:0:0:0]`: `(5,3,3,3)k/7` → ages `2,3,1,3,1,2`. No junior (age < 1) fixed point exists, and in any case `σ^k` does not fix the target character. `τ` fixes two plane sextic curves with transverse type `½(1,1,1)`, age 3/2 (not junior). This is the audit's `juniorElementsFixingTarget: []`.

### 1.6 Equivariant Hodge data (rigorous, Koszul character series; Appendix A.1)

Coefficient of `t^{6q} u^k v^1` in `Π_i (1 − t^5 u^{−w_i} v^{e_i})/(1 − t u^{w_i} v^{e_i})` (`u^7 = v^2 = 1`) gives the dimension of the Y‑odd Jacobian‑ring sector of σ‑character `k` in degree `6q`, i.e. `h^{4−q,q}` of the primitive τ‑invariant eigenspace:

| character k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | total |
|---|---|---|---|---|---|---|---|---|
| q=1 (3,1) | 19 | 37 | 37 | **26** | 37 | 26 | 26 | 208 |
| q=2 (2,2) | 128 | 123 | 123 | **123** | 123 | 123 | 123 | 866 |
| q=3 (1,3) | 19 | 26 | 26 | **37** | 26 | 37 | 37 | 208 |

So the sector containing `A` (character 3) has Hodge numbers `(26,123,37)`; its complex conjugate is the character‑4 sector `(37,123,26)`. The Galois orbits under `Gal(Q(ζ₇)/Q) = (Z/7)^*` are `{1,2,4}` (quadratic residues) and `{3,5,6}`. The σ‑trivial sector `(19,128,19)` is exactly what §4.3 supports. Dimension count of raw targets: 26 independent (3,1)‑directions share A's character; a certificate for all of them is a priori a different problem from a certificate for `A` alone.

Y‑even Jacobian‑ring dimensions by character (needed for E1 bookkeeping): `(R_F)_6^+ = [38,30,30,30,30,30,30]`, `(R_F)_9^+ = [83,93,93,87,93,87,87]`, `(R_F)_{10}^+ = [96,111,…,111]`, `(R_F)_{12}^+ = [111,129,…,129]`; Y‑odd: `(R_F)_{15}^- = [83,90,…,90]`, `(R_F)_{16}^- = [57,65,65,72,65,72,72]`.

### 1.7 E1‑admissible denominators (new, Appendix A.3)

Any certificate denominator satisfies `A·D ∈ J(F)` (Euler normalization, ch. 1). Modulo `J(F)` this is the kernel of multiplication by `[A]: (R_F)_d^+ → (R_F)_{d+6}^-` (maps character `k` to `k+3`). Computed mod 32003 (and the F_p rank is a lower bound for the Q rank, see A.3 for the semicontinuity argument):

| d | dim S_d^+ | dim J(F)_d^+ | dim (R_F)_d^+ | kernel of ·[A] (by character 0..6) | total kernel |
|---|---|---|---|---|---|
| 6 | 236 | 18 | 208 | 0,0,0,0,0,0,0 | 0 |
| 7 | 396 | 63 | 333 | 0,0,0,0,0,0,0 | 0 |
| 8 | 651 | 168 | 483 | 0,0,0,0,0,0,0 | 0 |
| 9 | 1001 | 378 | 623 | 2,3,3,2,10,1,3 | 24 (upper bound over Q) |
| 10 | 1512 | 750 | 762 | 24,46,39,39,54,46,46 | 294 (upper bound; generic lower bound is also 294) |

**Consequence (rigorous).** Every Y‑even legal denominator of degree `≤ 8` lies in `J(F)`: `D = Σ_i g_i F_i` with `deg g_i = d − 5 ≤ 3`, and then `N^0 = A·g` is a canonical E1 lift. For `d = 9` the admissible space is `J(F)_9 ⊕ (≤ 24‑dimensional complement)`. `J(F)_d^+` by character: d=6: `[12,1,1,1,1,1,1]`; d=7: `[3,12,12,8,12,8,8]`; d=8: `[51,20,20,19,20,19,19]`; d=9: `[33,63,63,52,63,52,52]`; d=10: `[159,101,101,96,101,96,96]`. (Odd denominators `D` can be replaced by the even `D·τ(D)` of double degree; the odd‑part kernel was not computed.)

### 1.8 Character bookkeeping for a search

For a σ‑eigen denominator of character `χ_D` (Y‑even), averaging (ch. 1, `T_γ`) shows one may assume: `char(N_i) = 3 + w_i + χ_D`, `char(Q) = 3 + 2χ_D`, Koszul correction coefficients `p_ij` (`K = Σ_{i<j} p_ij (F_j ∂_i − F_i ∂_j)`, `deg p_ij = d − 4`) have Y‑parity `1 + π_i + π_j` (`π = (0,0,0,1,1,1)`) and character `3 + w_i + w_j + χ_D`. Restricting a *non*‑eigen `D` to a character is not allowed (averaging multiplies degrees by up to 14).

### 1.9 Sizes of the complete fixed‑D linear system (exact identity `E = F²Q`), q = 1

Unknowns = coefficients of the six `N_i` (degree `d+1`, parities and characters as above) and of `Q` (degree `2d−6`); equations = Y‑odd monomials of degree `2d+6` and character `3+2χ_D`:

| d | dim D‑space (χ_D=0 / 1) | unknowns N | unknowns Q | equations |
|---|---|---|---|---|
| 6 | 50 / 31 | 343 / 357 | 29 | 2393 |
| 7 | 39 / 66 | 522 / 543 | 87 | 3783 |
| 8 | 120 / 89 | 868 / 879 | 209 | 5753 |
| 9 | 116 / 156 | 1238 / 1282 | 436 | 8468 |
| 10 | 255 / 212 | 1890 / 1897 | 824 | 12121 |
| 12 | 495 / 440 | 3698 / 3701 | 2393 | 23169 |
| 15 | 1044 / 1133 | 8609 / 8725 | 8468 | 53503 |
| 18 | 2503 / 2402 | 18265 / 18269 | 23169 | 109537 |

Without character restriction (needed for non‑eigen `D`): d=9: N 9009, Q 3080, equations 59332; d=12: N 25704, Q 16797, equations 162248; d=15: N 61047, Q 59332, equations 374604. The correction‑map form (ch. 5.2: fix an E1 lift, unknowns `p_ij` of degree `d−4`, reduce modulo `F` by the `t`‑monic remainder) is much smaller: d=9, χ_D=2: 240 unknowns, 1502 rows (fqlx's numbers, reproduced); d=11 (`xyztrs·F_i`): ~850 unknowns, ~11k rows; d=15 (`xyztrs·Jhat`): 4550 unknowns, 11158 rows.

---

## 2. Companion and related solved instances

`provenance.json` lists no `sameVarietyInstances`. `why.md` points to the uncoupled Klein examples and to §4.3. Instances sharing the Klein block or the target:

### 2.1 `klein-rational` (accepted; benediamond, PR #42; discussion #105)
`F = K(X) + K(Y)`, `A = X0^3 X1^2 Y0`. Method: cyclic monomial map `U = X0^3X1^2, V = X1^3X2^2, W = X2^3X0^2`, Fourier coordinates over `Q(ω)` (ω³ = 1), logarithmic fields transported by the rational monomial Piola formula, descended by normalized trace. Denominator (explicit, degree **30**, strict `Q = 0`, numerators of degree 31, only the three X‑components nonzero):
`D = Q^3 + 9 S^3 P`, `P = U^2W + UV^2 + VW^2`, `Q = U^2+V^2+W^2−UV−UW−VW`, `S = U+V+W`.
benediamond's comment (#105, 2026-10-24) factors it: `D = (Z^3 + ω² S^3)(Z'^3 + ω S^3)` with `Z = U+ωV+ω²W`, `Z' = conj`, two conjugate degree‑15 factors; and `D = Norm_{Q(ζ₉)/Q}(Z + ζ₉² S)`, six degree‑5 factors over `Q(ζ₉)`; on `X0 = 0` it degenerates to `V^6 = X1^18 X2^12`.
**Why it fails here.** (a) The construction is an X‑block construction with a *passive* `Y0`; our target has Y‑degree 3 (`Y0^3` is not a passive multiplier: `v(Y0^3) ≠ 0` for the transported field). (b) It relies on the additive split `K(X)+K(Y)` (independent blocks); the coupling terms kill every block decomposition (§1.5 (iii)). (c) The served target packet (nine directions `U,V,W × Y_j` by cyclic rotation and passive multipliers) does not contain `X0^3Y0^3` even in the split limit.

### 2.2 `klein-hesse` (accepted; benediamond, PR #22; discussion #85, no comments)
`F = K(X) + 2 (X0X1X2)^2 + K(Y)`, `A = (X0^3X1^2 − X1^3X2^2) Y0`. Method: Hesse relation among `U,V,W`, logarithmic and complementary Euler fields, rational monomial Piola transport. Denominator `D_κ = (κ² − 3κ + 9) r^3 f_κ − p^3` (`r = X0X1X2`, `p = U+V+W`), degree **15** at κ = 2 (coefficient 7), strict, numerators of degree 16. Serves a six‑dimensional packet of zero‑sum targets. Fails here for the same reasons as 2.1 (passive `Y0`, split structure, served packet).

### 2.3 Accepted instances with the same target `A = X0^3Y0^3` (all on Fermat‑type blocks)
* `holdout` (fqlx): `F = Σ(X_i^6+Y_i^6) + X0^3X1Y0Y1 + X1^3X2Y1Y2`. Adjoin `T^6 = Y0^6+Y2^6`, `α = Y0/T`, `β = Y2/T`; the equation becomes symmetric under `y = Y1 ↔ T`; quotient coordinates `S = y²+T²`, `P = yT`; moving polynomial conic on the weighted quotient (nonsingular seed over an algebraic number), polynomial integration, normalized trace, deck‑character selection `Tr(T^{−4}) = 0`. **Existential**; no explicit `D`, no degree bound.
* `holdout-ternary-deformation` (hodgepodgeconjecture, PR #151): same cover and source map for `F = Σ(X_i^6+Y_i^6) + X0^3X1Y0Y1 + X1^3X2Y1Y2 + (X0^5X1 + 2X1^5X2 + X2^5X0)/14` — note the *Klein monomials* appear as a deformation of the Fermat X‑block; new modulo‑5 conic seed; existential.
* `cycle-three`, `mixed-biweight`, `off-paired`‑type family (hodgepodgeconjecture/fqlx): paired‑sign cubic `U(x,y,z)` over passive ratios, Weierstrass Jacobian `η² = ξ³ + s_4 ξ − s_6/4`, polynomial curve with nonsingular mod‑5 seed, Jacobian residue primitive, first variation in a coefficient; existential.
* `one-reflection`, `holdout-moving-branch`, `family-hesse-web`, `family-block-cube`: reflection/binary‑block constructions (explicit witnesses for the family entries).
**Why they fail here.** The `holdout` cover needs the Y‑dependence to be `Y1·(…)` couplings with the pair `(Y0,Y2)` entering only through `Y0^6+Y2^6` and products with `Y1`, so that `y ↔ T` is a symmetry; our couplings involve all three pairs `Y0Y1, Y1Y2, Y0Y2`, and the Klein block has no `X_i^6` terms, so no exchange of the required kind exists after any of the three analogous covers. The paired‑sign cubic needs `F` even under each `(X_i,Y_i) ↦ (−X_i,−Y_i)`, false for `x^5y`. Reflection constructions need a coordinate reflection (none).

### 2.4 The book's own construction on this exact family (docs/04 §4.3) — and why it misses `A`
Rational coordinates `U = x²y/z, V = y²z/x, W = z²x/y` (σ‑invariant, exponent matrix determinant 14, kernel `⟨σ, ±1⟩`), `p = aU²W + bV²U + cW²V + dUVW`, `P = p + UVQ0 + UWQ1 + VWQ2 + g`; `F = P` and `𝒲 = {P=0} ⊂ P(2,2,2,1,1,1)` is the quotient `X/⟨σ,τ⟩` (degree 14), quasismooth, terminal `½(1,1,1)` along the cubic `p=0` at `Y=0`, `−K = O(3)`. Projection to `P²_Y` gives smooth cubic surfaces `p + T[UVQ0+UWQ1+VWQ2] + T³g = 0` (`T = Y0²`, passive `α = Y1/Y0, β = Y2/Y0`); the six‑coordinate Jacobian is `28xyz/Y0`. The 19 σ‑trivial, τ‑odd directions `xyz·k[Y]_3 ⊕ span{x³y², y³z², z³x²}·k[Y]_1` are served strictly by the cubic‑surface theorem (ch. 2 §2.4: common active denominator = the degree‑45 tritangent product, plus chart/passive poles). **The density `28xyz/Y0` and all coordinates are σ‑invariant, so every served density has character 0; `A` has character 3.** Averaging over σ kills `A` (why.md). No quotient of `X` by a subgroup containing a nontrivial `σ^k` can serve `A` (§3.1 averaging‑ideal criterion); the only subgroup fixing the character is `⟨τ⟩`, whose quotient is not uniruled (§4.7: `X/τ` has κ = 0) and whose `H^4` is not divisor‑supported by any known mechanism.

---

## 3. Prior work, item by item

### 3.1 The discussion thread (#116)
Exactly one comment, no replies.

**`fqlx`, 2026-10-26T03:27:03Z (node `DC_kwDOVLOOJc4BJMgE`), "Partial update".** Prepared with Codex/GPT‑6. States: target unchanged, "step 23 remains unsolved" (an internal numbering of fqlx's program; steps 1–22 are not public). Two obstruction claims:
1. On the "source‑adapted degree‑21 support at ε = 1", the three fixed `t/r/s` Hamiltonians, individually balanced at the 13 ordinary residue orbits, give distinct compact controls; cofactor elimination forces the only possible rational base‑amplitude candidate; after clearing its denominator the necessary identity has residual `(39,22,26,84,83) ≠ 0 mod 101` at the specified good point. Excludes that three‑control family with arbitrary rational amplitude poles.
2. A single normalized variable‑trigonal trace with arbitrary rational `w(a,b)`, one common seed, cannot produce the boundary source `K_η`: on a proper resolution of the degree‑7 density cover every polar curve is rational except one deck‑stable curve; the seven target norm‑pole curves are Klein curves permuted cyclically by the deck action; properness, residues and equivariance give a contradiction. Infinity is regular.
Explicitly left open: collective residue cancellation before individual balancing, independent seeds/weights, additional compact controls, fiber‑dependent amplitudes. Link: gist `https://gist.github.com/fqlx/a66456607d1bb86dee4294fa2128f692/06d3d8e2bd59cb6496b1a31e9b9399fdebf82d0a` (SHA‑256 manifest included).

### 3.2 The reproduction packet (gist, 29 files, created 2026-10-26T03:25:49Z)
README: "Coupled Klein: source reduction and scoped obstructions (2026-10-25)". Status line: partial research only; no certificate, Lean theorem, PR, or unrestricted no‑go. Scripts need Python 3 + SymPy; recorded JSON included; "No new numerical replay was performed merely to publish this packet."

The packet is built on a private coordinate framework from earlier "rounds" (round14–24 inputs are cited by file name but **not included**: `round24-S10/S08/S04-complete-response.txt`, `round14-native-comparison-proof.txt`, `round14-D27-compact-proof.txt`, `round24-S08-elliptic-boundary-residue-module.txt`). What is public:

**(a) The S08 chart (verified by me, Appendix A.2).** `c = x³z/t⁴, d = xz⁵/t⁶, τ = yt²/(xz²), a = r/t, b = s/t, h = xt²/z³` with `h⁷ = c⁴/d⁵`, and
`F/t⁶ = H := τc² + (10ab/7)c + 1 + a⁶ + b⁶ + (1 + 15bτ/7)d + (5a/7)τ³cd + τ⁵cd²`.
`c,d,τ,a,b` are σ‑invariant and `h` has character 3 (= the target's). So `{H = 0}` is (an affine model of) the quotient `X/⟨σ,τ⟩`, the fibres over `(τ,a,b)` are plane curves of genus 16 in `(c,d)` (the "actual S08 elliptic fiber"/"genus16" of the notes refer to this curve and its elliptic quotient), and the "physical density" `Ω = h·dc/H_d` carries the character of `A` so that `A·vol/F²` is written as `(14x³/t³)`‑multiples of `h`‑twisted forms on the quotient. The target reduces to a divergence equation `u(H) = 1/14`, `div_α u = source` for a five‑component field `u` on `(c,d,τ,a,b)` with weighted divergence w.r.t. `h·dc dd dτ da db`; the original field is recovered by an explicit inverse Piola map (`V_x = 0, V_y = 14x³y/t³(u_τ/τ − u_c/(7c) + 3u_d/(7d)), …`) and the original defect is `F·(14x³/t³)·R(H)` where `R` is the remaining scalar defect. A certificate is obtained iff `R(H)` is divisible by `H`.
**(b) `round25-S10-S08-five-class-Hermite-physical-transport.txt`.** The compact (second‑kind) cohomology of the genus‑16 fibre in the `h`‑character has rank 5 with basis `β_j = Ω·e_j`, `e = (1, 1/c, d/c, d, cd)` (three holomorphic, two second‑kind). Finite "Hermite" recipe for the five Gauss–Manin connection matrices `M_i` (`i = ℓ,τ,a,b`) and their primitives (using the rank‑3 algebra of critical points of `H_d` and the two Kummer points `d = 0`). Reduces the residual to a five‑vector `r5 = (τ/14)(∂_τ + M_τ) M_ℓ[:,0]`. Also a second chart "S10" `(p,q,w,a,b)` with a cubic `f = G w³ + L w + C` in which the two "physical" first‑`f` coefficients `Rlin = r0 + r1 w` are the obstruction; full off‑`F` quotient ledger (`Qcompact`). Claim/derivation; not independently re‑derived here.
**(c) `round25-S08-normal-parameter-scaling-current.txt`.** Scaling vector field `E = τ(3+2s)∂_τ − a(3+s)∂_a + bs∂_b + (−6−2s)c∂_c + (−9−5s)d∂_d`, `s = 15bτ/7`, with `E(H−ℓ) = (−9−2s)(H−ℓ) + J`, `J = (9+2s)(1−ℓ) − (9+4s)a⁶ + (9+8s)b⁶`. Normal field `u_E = E/(14 J_H)` has `u_E(H) = 1/14`, `div_α u_E = η₀(H)/14`, `η₀ = (N0 + sN1 + s²N2)/J²`, `N0 = −27(1−ℓ+5a⁶+b⁶)`, `N1 = [69(1−ℓ) − 867a⁶ − 561b⁶]/7`, `N2 = −8(1−ℓ+19a⁶+46b⁶)/7`. The remaining source is the single scalar `η₀` times the class `β₀`. In original coordinates `J_H = Jhat/(7xz²t⁶)` with `Jhat` as in (g). Derivation; the resulting polynomial identities in (g) I verified.
**(d) `round25-S08-eta-factored-two-physical-coefficients.txt`.** Explicit factored formulas for the two physical coefficients `r0, r1` of the `η₀`‑source in the S10 cubic field, with a small cubic norm `N = Gα³ + Lαβ² − Cβ³` (`α = 9(1−a⁶+b⁶)`, `β = (15b/7q)(2−4a⁶+8b⁶)`). Derivation only.
**(e) `round25-new-support-thirteen-orbit-compact-control.txt`.** New support `P_ε = x¹⁴t⁷ − (yz)⁷(t⁷+r⁷+s⁷) + ε x³t³y⁷z⁸` (degree 21, σ‑invariant; I verified `P_ε = (yz)⁷t⁷·J/(τ⁷d⁵)` with `J = c⁴ + τ⁷d⁵(εc − Z7)`, `Z7 = 1+a⁷+b⁷`). Hamiltonian field `W = t⁵y⁶z⁷(F_t∂_x − F_x∂_t)/(xP_ε)` (character 3, τ‑odd), pushforward `z = K0(T·A0 + X·B0)` with `A0 = 3c∂_c + d∂_d − τ∂_τ`, `B0 = 4c∂_c + 6d∂_d − 2τ∂_τ + a∂_a + b∂_b`, `X = A0(H)`, `T = 6H − B0(H)`, `K0 = τ⁶d⁵/(14c²J)`. `J` has 13 simple zeros on the generic fibre (degree‑13 map `ε = Z7/c − c³/(τ⁷d⁵)`), giving a rank‑13 residue extension; explicit 13×13 (`L(P0+5P1+7P∞)`) and 16×16 (`L(2P0+6P1+8P∞)` vs. five compact classes + 11 exact) matrices, invertible by Riemann–Roch; individually residue‑balanced compact control `wC` with remaining source `η₀e₀/14 − C_ε`. Construction; not re‑derived.
**(f) `round25-three-Hamiltonian-eta-cofactor-current.txt`, `round25-thirteen-coupled-Hamiltonian-current.txt`, `…-construction.py`, `…-algebraic-elimination.py`, result JSONs.** The three Hamiltonians `W_t, W_r, W_s = m0(F_•∂_x − F_x∂_•)`, `m0 = t⁵y⁶z⁷/(xP_ε)`, their compact fluxes (formula (4): `t: (−τ(T0+2B), aB, bB)`, `r: (−τR, −B, 0)`, `s: (−τS, 0, −B)`), amplitude product rule (6)–(7), cofactor row `ω` of `[H0,R,S,B]`, first algebraic pivot `λ0 F + λ_g g + λ_k k = θω₀`, second pivot `ν`, third pivot `χ`, unique forced rational candidate, and the cleared necessary identity (25). Numerical receipt at `p = 101`, `τ=a=b=ε=1, ℓ=0`: `det E13 = 48`, `det M16 = 43`, `ω = (5,67,82,87,28)`, `(λ0,λ_g,λ_k) = (70,63,79)`, `f0 = 48`, `ν = (0,68,83,17,13)`, `σ=42, δ=39, ρ=23`, `χ = (0,0,30,15,46)`, third pivot `λ = 0` but `∇λ = (65,96,49) ≠ 0` (so `λ ≢ 0`), `r_source = 6`, final cleared residual `(39,22,26,84,83) ≠ 0`. Equal‑amplitude subfamily excluded by a rank argument (`Cstar, H0, θe₀` have rank 3). Status: rigorous *within* the framework, if the framework's identities hold; the framework itself (the pushforwards, the finite projections, the Taylor‑4 caches) is not independently auditable from the public files. Scope is one explicitly specified three‑parameter family of currents; it says nothing about other denominators.
**(g) `round25-root-eta-original-degree9-partial.txt` + `round25-eta-D9-Koszul-fixed-family-proof.txt` + Koszul scripts/JSON.** Original‑coordinate denominator
`Jhat = 7xz²[9(t⁶ − r⁶ + s⁶) − 9F] + 15yst[2t⁶ − 4r⁶ + 8s⁶ − 2F]` (degree 9, 14 terms after expansion, σ‑character **2**, `Jhat = 63` at `[1:−1:1:1:0:0] ∈ X`, so `gcd(Jhat,F)=1`), `D = 14·Jhat`, and the explicit degree‑10 E1 lift
`N_x = −x⁴t³(147xz² + 75yst)`, `N_y = −x³yt³(147xz² + 45yst)`, `N_z = −x³zt³(147xz² + 195yst)`, `N_t = 0`, `N_r = −14x³rt³(21xz² + 15yst)`, `N_s = 210x³ys²t⁴`,
with `N(F) = A·D` and cleared defect `E = F·(R0 + F·Q12)`,
`R0 = 84x³t³[−3087x²z⁴(t⁶+5r⁶+s⁶) + 35xz²yst(69t⁶−867r⁶−561s⁶) − 600y²s²t²(t⁶+19r⁶+46s⁶)]`,
`Q12 = 259308x⁵z⁴t³ − 202860x⁴z²yst⁴ + 50400x³y²s²t⁵`.
Exclusion: with `Λ(P) := [x³y²t⁵s⁸] rem_F(P)` (`t`‑monic remainder, `F = t⁶ + …`), `Λ(−R0) = 2268000` while `Λ(L_ij(p)) = 0` for every Koszul column `L_ij(p) = D[(∂_ip)F_j − (∂_jp)F_i] − p(F_jD_i − F_iD_j)`, `deg p = 5`. Modular receipt: 240 character/parity‑compatible unknowns, 1502 rows, rank 197 mod 101, augmented rank 198; exact dual of support one. **Rigorous and reproduced by me over Q with all 3780 degree‑5 monomials (Appendix A.4).** By the regular‑sequence (Koszul syzygy) fact this excludes every polynomial numerator of degree 10 over `D = Jhat`; fqlx leaves "`D²` or higher pole order" open, but by the prescribed‑support theorem of ch. 1 (at `q = 1`, existence over `D` ⇔ existence over any `D^m`), **every certificate whose poles are contained in `V(Jhat)` is excluded**. The note also mentions two older, now "exhausted" degree‑9 candidates `D9 = x³(t⁶+r⁶+s⁶) + t⁴y⁵` and `H9 = αs + βt`, and a "D27 compact proof" — none of those files is public; status unverifiable.
**(h) `round25-Klein-trigonal-three-orbit-source-current.txt`, `round25-variable-trigonal-norm-pole-and-source-boundary-gate.txt`, `round26-variable-trigonal-source-support-equivariance-gate.txt`.** The curve `Γ: y⁷ = t³(t+1)⁵` (a quotient of the Fermat‑7 curve, `JΓ ∼ E_CM(−7)³`, "Klein curve") appears as the inverse‑density curve of moving poles; it is trigonal: `w = y³/(t(t+1)²)` gives `t³ + t² = w⁷`, `ω_Γ = dt/w²`, and the three conjugate points satisfy `Σ ω_Γ(t_i) = 0`. A three‑orbit residue‑null current with seed `ψ_i = λw⁵` has leading compact feedback `C(w,λ) = w⁴[w_a D_b^κ λ − w_b D_a^κ λ]·F(w⁷, Tλ⁷)`, `F(W,z) = [z(1−4W+W²) + 2z²W⁴]/[1 + z(1−7W+7W²) + 7z²W⁴ − z³W⁷]`, `T = −1/(7³·10⁴·a⁴b⁴G)`, `G = 1+a⁶+b⁶`; the boundary source it must match is `K_η = G(1+5a⁶+b⁶)/(9Q(1−a⁶+b⁶)²)`, `Q = 10ab/7`. The fixed `w = a` family fails by a degree count in `b`; a local pair `(w,λ) = (q³, q^{−7})`, `q = V/(T̄D)`, `V = −49a¹²(1−T_D³)/(270(2+7b⁶))` serves the double pole along `D = 1−a⁶+b⁶ = 0` but has a new double pole at `2 + 7b⁶ = 0`. The equivariance gate: on the degree‑7 density cover `η⁷ = −k0 A⁴B⁴F6` of `P²`, resolved, all polar curves of `κ K_η da∧db` are rational except one deck‑stable genus‑100 curve `C_D` over `D=0` (not a seventh power), the seven norm curves `B_ε: w⁷s³ − ε²s − ε³ = 0` are Klein curves with nonzero residue `(ε/49)dt/w²`, permuted by the deck action, so no single trace current can serve `K_η`. Status: plausible geometric argument, not independently verified; scope: one normalized family.
**(i) `round26-Fermat42-source-CM-proof.txt` + script/JSON.** Koblitz–Rohrlich (Can. J. Math. 30 (1978)) character computation on the Fermat‑42 curve: the genus‑100 curve `C_D: b⁶ = a⁶ − 1, κ⁷ ∝ a^{−10}b^{−4}` has **no** `E_CM(−7)` isogeny factor (26 Galois orbits listed), hence no nonconstant map to `Γ`; the genus‑18 curve at infinity `C_∞: y⁷ ∝ r^{−4}(1+r⁶)^{−1}` has an `E_CM(−7)⁶` block and an explicit degree‑4 map to `Γ` (`U = r², t_K = (U+U^{−1}−1)/3, w_K = c/(Uy²)`, `c⁷ = 1/27`). Plausible; standard inputs; not re‑derived.

### 3.3 What in the packet is rigorous vs. scoped vs. claimed
* Rigorous and reproduced: (g) the `Jhat` exclusion (all supports inside `V(Jhat)`); the chart identities in (a) and the support identity in (e).
* Rigorous within a non‑public framework: (f) the three‑control exclusion (exact nonzero of a specified rational expression at a good point), (h)/(i) geometric gates.
* Not established anywhere: any denominator family beyond `V(Jhat)` and the two unpublished degree‑9 candidates; any nonconstructive support statement for the ζ³‑sector; any lower bound on `deg D`.

### 3.4 Other threads
`cycle-q2` (#?) mentions Klein only in passing (a repair note in a different packet). Discussions #85 (klein‑hesse) is empty; #105 (klein‑rational) has benediamond's factorization comment (quoted in §2.1).

---

## 4. What remains unexplored (concrete gaps)

1. **No complete fixed‑D or joint search in low degree has ever been run** (publicly) — only `Jhat` (one degree‑9 eigen‑denominator of character 2). By §1.7 every legal Y‑even denominator of degree ≤ 8 lies in `J(F)`; the whole space is small (18 / 63 / 168 dimensions at d = 6/7/8, before the gauge `D ↦ D + hF`), and the bilinear joint system `−B_D(Ag) ∈ F·S + im T_D` (ch. 5.3) is tiny (§5, item 1). Not done for any degree; degree 9 and 10 (with the ≤24 / 294‑dimensional kernel parts) are also feasible by σ‑character.
2. **No denominator in a nontrivial σ‑character other than 2 was tried**, although `A·D ∈ J(F)` admits E1‑candidates in every character at d = 9 (kernel dims `2,3,3,2,10,1,3`; character 4 is the richest) and d = 10.
3. **The `(1,3)/(3,1)` asymmetry of the sector was never used.** Sector 3 has `h^{3,1} = 26 < h^{1,3} = 37`; a divisor support must carry a ζ³‑eigenspace in `H^{2,0}` of its resolved components; no candidate divisor (σ‑stable or a σ‑orbit) with holomorphic two‑forms in that eigenspace has been proposed. Natural candidates: `X ∩ {D = 0}` for σ‑eigen `D` from `J(F)` (polar‑type divisors), the degree‑21 support `P_ε` (fqlx), `V(Jhat)`.
4. **Non‑eigen denominators** were never considered (every tested or proposed `D` is a σ‑eigenvector). Averaging multiplies degree by up to 14, so an eigen search at fixed degree is strictly weaker than a general one.
5. **Polynomial Piola maps (§3.4 finite‑Fano theorem).** Whether `F = G(Φ)` for a nonlinear polynomial `Φ` with `x³t³ ∈ δ·k[Φ]` was never examined. Pure‑Y coordinate choices fail (the three products `Y_iY_j` and `Y0⁶+Y1⁶+Y2⁶` are not simultaneously in `k[ℓ₁,ℓ₂,q]` etc.), monomial maps fail (the Klein and Fermat exponents are incompatible with `Σ_i m_i = (4,1,1,4,1,1)`), but mixed non‑monomial maps are open.
6. **Covers with descent (§2.7/§4.6).** No auxiliary cover `Z → X` whose deck/exchange relations absorb the character ζ³ has been proposed. Note the obstruction: `H^4(X)^σ ⊃ H^{4,0}` is not supported, so §4.6's hypothesis "`H^4(X/⟨g⟩)` supported" can only hold with `g = στ` (order 14, quotient = §4.3's `𝒲`).
7. **Deformation/continuation from a solved neighbour.** The split limit `Q_i = 0` (`F₀ = K(X) + Fermat(Y)`) has three Y‑reflections, hence geometric support of every odd target (§3.3), but **no explicit certificate for `X0^3Y0^3` on `F₀` exists in the atlas either** (`klein-rational`/`klein-hesse` have Y‑degree‑1 targets). Nobody has computed a certificate at `Q_i = 0` (where the reflection quotient `{K(X) + u³ + Y1⁶ + Y2⁶ = 0} ⊂ P(1,1,1,2,1,1)` is Fano) and attempted numerical continuation in the coupling parameter `s` (`F_s = F₀ + s·M`), or examined how the denominator degree behaves along the family.
8. **The §4.3 cubic‑surface fibration was never exploited for the nontrivial characters at all**: over `k(P²_Y)` the fibre of `X → P²_Y` (not of the quotient) is the degree‑14 cover of the cubic surface `S_Y`; its ζ³‑isotypic `H^2` (a Prym‑like part) is what carries the target over the passive plane; whether it is algebraic on the generic fibre (which would give a relative support theorem as in §3.2 "surface fibers with algebraic second cohomology") has not been studied.
9. **Multi‑prime / exact linear‑algebra infrastructure for d ≥ 12** has not been set up (sizes in §1.9); modular rank evidence now exists for the polar‑type denominators `F_i`, `Z_jF_i`, `xyztrs·F_i` (all inconsistent, A.5) and for `Jhat`, but not for `Hess(F)`, for products of two polars, or for any denominator of degree ≥ 12 other than the `Jhat` multiples of A.5.
10. **Lower bounds on `deg D`** are absent; §1.7 plus a complete search up to d = 8 would give `deg D ≥ 9` (for Y‑even `D`) as a publishable partial result.

---

## 5. Candidate attack plan (experiments and arguments)

Sizes are exact counts from Appendix A unless marked "≈". All linear/bilinear systems are over Q with `F`‑reduction by the `t`‑monic remainder (`F = t⁶ + …`), which is what makes "modulo `F`" finite.

**1. Complete joint search over `J(F)_d`, d = 6,7,8 (rigorous target: "`deg D ≥ 9`").**
Idea: by §1.7 these are *all* Y‑even legal denominators; parametrize `D = Σ_i g_i F_i` (`g_i` monomials of degree `d−5`, parity `π_i`), remove the gauge `D ↦ D + hF` by a chart, take the canonical lift `N⁰ = A·g`, unknown Koszul coefficients `p_ij` (degree `d−4`, legal parity, all characters), equation `−B_D(Ag) ≡ T_D(p) (mod F)`: bilinear in `(g,p)`.
Sizes: d=6: 17 free `g`‑parameters, 162 `p`‑unknowns, 2345 reduced rows; d=7: 60 / 420 / ≈5170; d=8: 156 / ≈945 / ≈8673. Per σ‑character (eigen `D` only, complete by averaging): d=7, χ=0: 3/48/373; d=8, χ=2: ≈20/152/1182 (I ran these; see A.6).
Method: (i) numerical variable projection + Levenberg–Marquardt with many restarts (ch. 5.4) — evidence; (ii) exact: Gröbner basis of the bilinear ideal over `GF(p)` with the chart normalization, then over Q — a *proof* of nonexistence in each degree if the ideal is `(1)`. Chance of a certificate: low (a Y‑degree‑3 target on a coupled sextic with denominator ≤ 8 would be unprecedented), but the negative result is cheap and rigorous. Also run the Y‑odd‑`D` version (lift `N⁰ = A g` with flipped parities).

**2. Degree 9 and 10 by character, including the kernel part.**
Compute the exact Q‑kernel of `·[A]` on `(R_F)_d^+` (Singular over Q or two primes + reconstruction), lift kernel representatives `D_k` to E1 solutions `N_k` with `lift`, and run the joint search over `D ∈ span(D_k) ⊕ J(F)_d[χ]` per character (`J(F)_9^+` by character `[33,63,63,52,63,52,52]`, kernel ≤ `[2,3,3,2,10,1,3]`). Sizes per character at d=9: ≤ 73 `D`‑parameters, ≈240 Koszul unknowns, ≈1500 rows. `Jhat` is one point of the character‑2 slice; character 4 (10‑dimensional kernel) is the most promising unexplored sector. Chance: moderate for new exclusions, small for a witness.

**3. Natural candidate denominators, complete fixed‑D tests (cheap, modular + exact dual).**
(a) `D = xyztrs·F_i` (support = all coordinate hyperplanes ∪ one polar; degree 11, 792–859 unknowns, 2.4k–3.1k reduced rows) — **done, all six inconsistent mod 32003** (A.5); (b) `D = xyztrs·Jhat` (degree 15; 4550 unknowns, 11158 rows) and `D = F_i·Jhat` (degree 14; e.g. `F_x·Jhat`: 3179 unknowns, 10446 rows) — relaunched in the second session, see A.5 for whatever finished; (c) `D = Hess(F)` (degree 24, E1 lift `N⁰ = (5/6)adj(H)∇A` is automatic by ch. 5.2; correction unknowns of degree 20 — ≈57k per character; needs sparse modular elimination); (d) products of two polars `F_iF_j` with coordinate hyperplanes (degree 16). Each inconsistency is convertible into an exact proof by a sparse rational left‑kernel functional (as fqlx did), and each excludes the whole support by the prescribed‑support theorem. Chance of success: very low for (a),(b); but (c) is the canonical "Hessian" support and has never been tested on any coupled instance.

**4. Relative‑support argument on the §4.3 cubic‑surface fibration (nonconstructive).**
Over `K = k(α,β)` (passive `Y`‑ratios) and the chart `T = Y0²`, `X_K := X ×_{P²_Y} Spec K` is a degree‑14 Galois cover (group `⟨σ,τ⟩`) of the smooth cubic surface `S_K = {p + T[UVQ0+UWQ1+VWQ2] + T³g = 0}`. The target class is the ζ³‑part of `H^2(X_K)` pushed through the Leray/Gysin structure (§3.2 "surface fibers with algebraic second cohomology" requires `H^2` of the generic fibre to be algebraic). Compute the genus/Hodge type of the ζ³‑eigenspace of the cover surface `X_K → S_K` (it is a ramified cyclic‑7 cover of a cubic surface along a curve, `h^{2,0}` of the eigenspaces computable from the branch data). If the ζ³‑eigenspace of `H^2(X_K)` is of type (1,1) only, the surface‑fibration theorem gives support of the whole sector; if it has `h^{2,0} > 0` (likely, since the sector has 26 holomorphic‑type directions), this route is closed and the result still pins down the geometric obstruction. Cost: one genus/branch computation. Chance of a support theorem: low‑to‑moderate; chance of a decisive structural statement: high.

**5. Character transfer by an auxiliary cyclic cover (nonconstructive, §4.6 template).**
Seek `Z → X` with deck automorphism `d`, involutions `r, σ'` and `e = σ'dσ'` with `pe = (στ)p`, `pr = τp`, such that `Z/⟨r,σ'⟩` has supported `H^4`. Concretely: adjoin a 7th root of a σ‑semi‑invariant function of character 3 (e.g. `u⁷ = h⁷ = x⁷t¹⁴/z²¹`, i.e. `u = h = xt²/z³` itself is the degree‑0 function of character 3 that appears in fqlx's chart) and look for an exchange symmetry of the cover equation. The criterion requires finiteness and uniruledness/algebraicity of the intermediate quotient; none is in sight, but this is the only known mechanism in the book that can change characters. Chance: speculative.

**6. Numerical continuation from the split limit.**
Family `F_s = K(X) + G(Y) + s·M`, `s ∈ [0,1]`. At `s = 0` compute an explicit certificate for `X0^3Y0^3` via the `Y0`‑reflection quotient `{K(x,y,z) + u³ + r⁶ + s⁶ = 0} ⊂ P(1,1,1,2,1,1)` (Fano; §3.3 bounded incidence with ≤ 41 coefficients) or via an elliptic‑companion construction (3A) — itself a new result for the atlas family — then continue the fixed‑degree certificate system in `s` (ch. 5.4 variable projection; the solution set over `s` is an algebraic family; watch for degree jumps). Chance of reaching `s = 1` at fixed degree: low; but the `s = 0` certificate and the behaviour of `deg D` are informative.

**7. Polynomial Piola classification (small finite computation).**
Enumerate nonlinear polynomial maps `Φ` (degrees `m_i ≤ 5`, `7 ≤ Σm_i ≤ 12`) with `F ∈ k[Φ]` by a term‑semigroup search on the 9‑term support (for monomial `Φ` the condition is lattice/semigroup membership; for non‑monomial `Φ` restrict to σ‑semi‑invariant coordinates). If none has `x³t³ ∈ δ·k[Φ]`, §3.4 is closed for this instance (a clean negative structural fact); if one exists, §3.4 gives a nonconstructive support proof immediately. Chance of existence: low (the Y‑block analysis in §4 item 5 already kills pure‑Y coordinates).

**8. Exact exclusion certificates for the modular results.**
For every inconsistent fixed‑D test (items 3, A.5) extract a rational left‑kernel vector with small support (sparse dual as fqlx's `Λ`), giving proofs of nonexistence per support that could be formalized (polynomial identity checks). Chance: high (mechanical); value: turns evidence into theorems and feeds a Lean‑checkable "denominator exclusion" library.

**9. Hodge‑theoretic structure of the ζ³‑sector (theory).**
The nontrivial‑character part of `H^4(X)^τ` is a `Q(ζ₇)`‑Hodge structure of rank 186 with signature `(26,123,37)` on the ζ³‑eigenspace. Determine whether it decomposes further (e.g. via the degree‑14 cover of `𝒲` and the Prym‑type decomposition of the fibre curves `{H = 0}` of genus 16, or via the Fermat‑7/Klein curve `Γ` that already surfaces in fqlx's packet with `JΓ ∼ E_CM(−7)³`). A decomposition into pieces with `h^{2,0}`‑type contributions coming from curves/surfaces with CM by `Q(√−7)` would point to the right correspondence (Shioda‑style inductive structure for the split limit `K(X) + Fermat(Y)`, deformed). Chance: this is the most plausible route to a *nonconstructive* proof but is open‑ended.

**10. Use the rank‑13 residue extension of `P_ε` as a denominator family in the original coordinates.**
fqlx's objects translate to original denominators `D = P_ε·(coordinate monomials)·(resultant factors)` of degree ≥ 21; a complete fixed‑D test at degree 21 with character restriction has ≈ `(1/14)·C(44,5)·(1 − (C(38,5)/C(44,5)))` ≈ 70k rows and ≈ `(1/14)·15·C(22,5)` ≈ 28k Koszul unknowns — large but within reach of sparse modular linear algebra (Singular/msolve/LinBox). It would settle, independently of fqlx's current framework, whether the degree‑21 support can carry a certificate. Chance of a witness: unknown; this is the one place where the packet suggests a specific divisor.

---

## 6. Feasibility assessment

**Explicit witness.** Everything known is consistent with a certificate existing only in high degree, if at all. Comparable targets with Y‑degree 3 on Fermat‑type blocks are all *existential* (holdout family, cycle‑three family) with no degree bound; the only explicit Klein witnesses have degrees 15 and 30 for passive‑`Y0` targets. A witness of degree ≤ 10 would have to appear in the complete low‑degree searches of §5 items 1–2; those are cheap and decisive and should be run before anything else. A witness of degree 12–24 is at the edge of exact linear algebra without strong structural guidance; degree ≥ 24 without a construction is hopeless.

**Nonconstructive argument.** Needs divisor support of the ζ³‑eigenspace (or of the class of `A` alone). All book mechanisms are character‑blind in the wrong direction (quotients kill ζ³; reflections, paired signs, exchange covers are absent). A new ingredient is required: either an algebraic‑`H^2` statement for the ζ³‑isotypic part of the surface fibres (item 4), a correspondence changing the character (item 5), or a CM/Prym decomposition (item 9). The §4.4 "one generic target suffices" theorem does not apply (this is an isolated sextic).

**What a Lean proof would require.** (a) Constructive route: `proof.json` + `Solution.lean` + `Witness.lean` with the explicit `(D, N, Q)`, checked through `HodgeComputer.Verify` (exact identity `E = F²Q` over `ℂ` via rational coefficients, homogeneity, `D ≠ 0`, `gcd(D,F) = 1` — the last usually via a point of `X` where `D ≠ 0` plus irreducibility of `F`, as in the accepted Klein proofs). Witness sizes at degree 15 are ~60 KB (klein‑rational's `Witness.lean` is 65 KB at degree 30, strict); a non‑strict witness of degree 20 would have `Q` of degree 34 with up to ~10⁵ terms — still kernel‑checkable with the project's `AlgebraicSpecialization`/digit‑check style, but large. (b) Existential route: the shared `HoldoutCover`/finite‑residue‑clearing libraries formalize "primitive in the function field ⇒ homogeneous polynomial certificate"; what is missing is any finite‑algebraic construction of a residue primitive for the ζ³‑character (a moving‑curve/incidence family on an auxiliary variety reaching this character). Nothing in the packet is close to that interface: fqlx's remaining obligation is a rational function identity (`R(H) ≡ 0 mod H`) in a five‑variable rational framework that has no formal counterpart yet.

---

## Appendix A. Computations performed for this dossier (all scripts in `/tmp/hodge/scratch/kcc/`)

A.1 `structure2.py` (SymPy 1.12): parses `instance.rec`; SNF of `L(F)+Z1` = diag(1,1,1,1,1,84); explicit enumeration of the projective diagonal group (14 elements, weights mod 126 listed in the log) and residue characters of `A`; block permutations preserving the support (identity + two rotations, infeasible coefficient ratios); Koszul character series (table in §1.6); monomial counts for the system sizes (§1.9). `liealg.py`: the Lie algebra of linear symmetries is spanned by Euler.

A.2 `chart.py`: verifies `F/t⁶ = H` for fqlx's S08 chart, `h⁷ = c⁴/d⁵`, σ‑characters `(c,d,τ,a,b,h) ↦ (0,0,0,0,0,3)`, and `P_ε = (yz)⁷t⁷J/(τ⁷d⁵)`.

A.3 `jacring.sing` + `jacpost.py` (Singular 4, `GF(32003)`): `std(jacob(F))` (1798 elements, `vdim = 15625`, 1 s), `reduce(A) ≠ 0`, monomial bases `kbase` in degrees 6–10 and the rank of multiplication by `A` on the Y‑even part per character (table §1.7). Semicontinuity: the standard monomials mod p lift to a `Z_(p)`‑basis of the Jacobian ring (regular sequence mod p ⇒ flat), so the Q‑rank of `·[A]` is ≥ the F_p‑rank; hence the Q‑kernel is ≤ the F_p‑kernel (equal to 0 for d ≤ 8).

A.4 `verify_jhat.sing` (Singular over Q, 0.9 s). Caveat on provenance: the first saved log `verify_jhat.log` shows that the two identity checks on lines 15–16 did *not* execute (Singular error: `NF` is a reserved name, `E` undefined); the script was fixed (`NFv`) and rerun in the second session, log `verify_jhat_rerun.log`, with every check printing `1`. Results: `N(F) − A·D = 0` ✓; `E − F(R0 + F·Q12) = 0` ✓; degrees 10 ✓; `Jhat(1,−1,1,1,0,0) = 63`, `F = 0` there ✓; `LM(F) = t⁶` under `lp` with `t` first; `Λ(−R0) = 2268000` ✓; `Λ(F·x³y²s⁸t⁵) = 0` ✓; all **3780** Koszul columns (`15` pairs × `252` degree‑5 monomials, no character restriction) have `Λ = 0` ✓. Together with the Koszul‑syzygy fact (smooth `F` ⇒ regular sequence) and ch. 1's prescribed‑support theorem this proves: *no `q = 1` certificate has all poles inside `V(Jhat)`.*

A.5 `fixedD.py` (Singular builds the reduced Koszul columns with legal parity/character; NumPy Gaussian elimination mod 32003): complete fixed‑D correction‑map tests. Results (modular evidence; "inconsistent" = no numerator exists mod p, which is a proof only once a rational left‑kernel functional is extracted):
* `D = 14·Jhat`: 240 columns, 1502 rows, rank 197, inconsistent (matches fqlx's mod‑101 receipt).
* `D = F_i` (degree 5, each `i`): inconsistent (columns 0–6; the reduced RHS is nonzero, so these are exact exclusions of the polar supports `V(F_i)`).
* `D = Z_j F_i`, same‑block pairs (18 Y‑even denominators of degree 6): all inconsistent (9–30 legal columns, ranks 9–27, 81–199 rows; log `fixedD_small.log`, JSON `fixedD_small_32003.json`). For these and for `D = F_i` the column spaces have full column rank and the augmented rank is one higher, so a rational left‑kernel functional of small support certainly exists; it has not been extracted yet (that is §5 item 8).
* `D = xyztrs·F_i` (degree 11, Koszul degree 7), completed (`fixedD_polar11.log`, `fixedD_polar11_32003.json`), all **inconsistent** mod 32003:
  `F_x`: χ_D=6, 859 cols, 2746 rows, rank 644; `F_y`: χ_D=5, 859 cols, 2836 rows, rank 644; `F_z`: χ_D=3, 792 cols, 2558 rows, rank 605; `F_t`: χ_D=0, 792 cols, 3071 rows, rank 592; `F_r`: χ_D=0, 792 cols, 2400 rows, rank 590; `F_s`: χ_D=0, 792 cols, 2378 rows, rank 592. (Here `χ_D = −w_i`, E1 lift `N⁰ = A·xyztrs·e_i`.) Note the large column‑rank deficiency (≈ 200 of ≈ 800 Koszul coefficients act trivially modulo `F`), as expected from the second Koszul syzygies.
* `D = 14·Jhat` recheck after fixing a script bug (the JSON `fixedD_jhat_32003.json` from the buggy run shows `cols = 0` and must be ignored; the valid numbers are in `fixedD_jhat.log`): 240 cols, 1502 rows, rank 197, inconsistent.
* `D = F_i·14Jhat` (degree 14, Koszul degree 10) and `D = xyztrs·14Jhat` (degree 15, Koszul degree 11; 4550 cols × 11158 rows, build 81 s): the first session died during the dense elimination of the degree‑15 matrix (machine at load ≈ 180, ≈ 1 GB transient per pivot step). Relaunched in the second session as `fixedD_resume.py` → `fixedD_resume.log` / `fixedD_resume_32003.json`; `F_x·14Jhat` has χ_D = 1, 3179 cols, 10446 rows, |RHS| = 279 nnz. Results are appended in A.8 if they finished before this dossier was closed; otherwise read the log.

A.6 `joint.py`, `joint_char.py` (floating point, variable projection + Levenberg–Marquardt, random restarts): joint bilinear searches `−B_D(Ag) ≡ T_D(p) (mod F)` over `D ∈ J(F)_d^+` (all Y‑even legal `D`, modulo the gauge `D ↦ D + hF`), d = 6 for all `D` and d = 7 per σ‑character. These are **heuristic evidence only** (a nonconvex least‑squares search; a relative residual near 1 means the best fit is essentially the zero correction, a small nonzero residual means nothing rigorous). Final results (the processes were killed with the first session):
* d = 6, all characters (`joint6.log`): 18‑dimensional `J(F)_6^+`, gauge removes 1, 17 free `D`‑parameters, 162 Koszul unknowns, 2345 reduced rows. 4 restarts (121–524 s each): relative residual `1.000` in every restart, `rank M = 144/162` at every visited `D`.
* d = 7 per character (`joint_char_78.log`; `J(F)_7^+[χ]` dimensions `3,12,12,8,12,8,8` for χ = 0..6, matching §1.7): χ=0 (3 params, 48 unknowns, 373 rows): best 1.000 over 10 restarts; χ=1 (11 params, 58 unknowns, 682 rows): best 0.861; χ=2 (11, 48, 507): best 1.000; χ=3 (8, 48, 579): best **0.179** (two restarts converged to the same local minimum, the other eight to 0.9975); χ=4 (11, 102, 577): 8 restarts, best 0.996 (killed during restart 8); χ=5 (8, 58, 642 rows): matrix built, no restart finished; χ=6: not reached.
* d = 8 and d = 9: never ran to completion (one foreground probe at d=8, χ=2, 152 unknowns, 1182 rows, gave relative residual 0.24 after a single restart — not a solution). `joint_char_9.log` does not exist.
Conclusion: no certificate with Y‑even `D` of degree 6 or 7 was found numerically, but this is not an exclusion; the rigorous version (Gröbner basis of the bilinear ideal over GF(p), §5 item 1) remains to be run — it is the single cheapest decisive experiment on the list.

A.7 Environment: local Singular 4.x, SymPy 1.12, NumPy/SciPy 1.18; the `ddg-cas` docker image provides Singular and Julia only (no Sage/M2/SymPy). The machine was heavily loaded (load average ≈ 180 on 40 cores, 125 GB RAM nearly full) during the first session, which is the likely cause of the crash; load ≈ 20 during the second.

A.8 Status of the relaunched `Jhat`‑multiple tests at the time this dossier was closed (2026-11-01 18:13). `fixedD_resume.py` (order: `F_x·14Jhat`, `F_y·…`, `F_z·…`, `F_t·…`, `F_r·…`, `F_s·…` at d = 14, then `xyztrs·14Jhat` at d = 15) was running as an orphaned `nohup` process; the first matrix (`F_x·14Jhat`, χ_D = 1, 3179 Koszul columns, 10446 reduced rows, |RHS| = 279) was built in 16 s and was in dense elimination mod 32003 when the machine load reached ≈ 390 on 40 cores, so no rank had been printed yet. Each finished test appends one `mode=dense: rank(M)=…, consistent=…` line to `/tmp/hodge/scratch/kcc/fixedD_resume.log` and rewrites `/tmp/hodge/scratch/kcc/fixedD_resume_32003.json`; `consistent=False` means "no polynomial numerator of degree d+1 exists over this D mod 32003" (evidence), `consistent=True` would mean a modular solution exists and must then be lifted to Q and checked exactly (it would be the first positive signal on this instance). Expected: inconsistent, by analogy with every polar‑type test so far; but the d = 14/15 tests are the only ones that probe a support strictly larger than `V(Jhat)` ∪ coordinate/polar hyperplanes, which is why they are worth finishing. If the process has died, rerun `cd /tmp/hodge/scratch/kcc && nohup python3 fixedD_resume.py 32003 > fixedD_resume.log 2>&1 &` on a quiet machine (≈ 1 GB RAM transient per pivot step, tens of minutes per test when unloaded).

## Appendix B. Links
* Discussion #116: https://github.com/hodge-computer/sextics/discussions/116
* fqlx packet: https://gist.github.com/fqlx/a66456607d1bb86dee4294fa2128f692/06d3d8e2bd59cb6496b1a31e9b9399fdebf82d0a (files mirrored at `/tmp/hodge/downloads/klein-gist/files/`)
* Companion threads: #85 (klein‑hesse, empty), #105 (klein‑rational, benediamond's factorization), #147 (holdout‑ternary‑deformation partials), #44 (holdout, target‑selection identities)
* Book: `docs/04-known-locus-and-frontier.md` §4.3 (the family and the 19‑direction theorem), §4.6–4.7 (descent and scope), `docs/03-geometric-support.md` §3.1 (averaging ideal), §3.3 (reflections, ages), §3.4 (polynomial Piola), `docs/02-explicit-constructions.md` §2.4–2.5 (cubic surfaces), §2.7–2.8 (covers), `docs/01-certificates.md` (prescribed‑support theorem, character averaging), `docs/05-computation.md` §5.2–5.4 (correction map, joint search).
* Audit script: `scripts/audit-frontier.py` (semantics of `checks.json`: `diagonal_screen` enumerates weights `(a,−5a,25a,21b,21c,0) mod 126`; `center_dimension` = Hessian symmetrizer; `euler_atoms` = additive block splits balanced for `A`).
