# Research dossier: `eight-plane-configuration` (OPEN)

Compiled 2026-11-01 from: `data/eight-plane-configuration/{instance.rec, why.md, provenance.json, checks.json}`, the complete GitHub discussion #176 (https://github.com/hodge-computer/sextics/discussions/176 — exactly one top-level comment, no replies, saved at `/tmp/hodge/discussions/eight-plane-configuration.md` and `.json`), `scripts/audit-structural-frontier.py` (the script that produced `checks.json`), the book chapters `docs/02-explicit-constructions.md` (§2.2–2.3 contact constructions), `docs/02a-binary-blocks.md` (§2a.1), `docs/03-geometric-support.md` (§3.1, §3.2, §3.3), `docs/04-known-locus-and-frontier.md` (§4.2, §4.4, §4.7), `formal/README.md`, the `solution.md`/`witness.json` of the accepted contact-type instances (`contact46006`, `family-passive-contact`, `polarized-contact-coupled-d2`, `polarized-contact-deformation-d2`, `completed-square-nambu-holdout`) and of the accepted instances that carry the *same* target `X0^3·Y0^3` (`family-block-cube`, `family-hesse-web`, `holdout`, `one-reflection`, `cycle-three`, `mixed-biweight`, `holdout-ternary-deformation`, `holdout-moving-branch`, plus `fermat-pair-rational-d4-strict`, `coupled-pair-rational-d4-strict`), and the unpublished computational logs of this workspace (`/tmp/hodge/results/eight-plane-configuration/`, engine `~/devel/sextics-work/hodge/`, 2026-11-01).

Conventions. Variables `(X0,X1,X2,Y0,Y1,Y2)`; I write `Q_i = X_i^2 + Y_i^2`, `τ : Y ↦ −Y` (the global involution), `i = sqrt(−1)`, `S = Q[X,Y]`, `J = J(F) = (∂F/∂Z_0,…,∂F/∂Z_5)`, `R = S/J` (Jacobian ring, Gorenstein with socle degree 24), `S_d^±` = degree-`d` forms with even/odd total `Y`-degree. "Rigorous" = proved by a cited argument or an exact computation I reproduced/inspected; "evidence" = modular computation; "claim" = stated in the thread without proof visible to me; "in flight" = a computation running while this was written (see §3.3 and the addendum at the end).

---

## 0. Status in one paragraph

Target: `q = 1`, `A = X0^3·Y0^3` on a dense `τ`-even sextic (all 236 `Y`-even sextic monomials present) whose only structure is `F ∈ (Q_0,Q_1,Q_2)` — it contains the eight planes `Y_i = ε_i·i·X_i`. No certificate, no existential proof and no Lean proof exist. The public discussion consists of one comment (`benediamond`, 2026-10-23) proving (i) the eight planes are infinitesimally rigid (`h^0(N_{Π/X}) = 0`), (ii) projection from one plane gives quintic residual surfaces with `p_g = 4, q = 0`, and (iii) the `τ`-invariant net `[Q_0:Q_1:Q_2]` has `(2,2,4)` complete-intersection fibres with `p_g = 19` of which **9** invariant holomorphic two-forms survive in the target's character — so neither natural fibration gives a "wholly algebraic invariant fibre `H^2`". Nothing in the thread excludes any denominator. Unpublished local computations (this workspace) add the one rigorous exclusion so far: for every `D` of degree `≤ 8` with even `Y`-parity and every `D` of degree `≤ 7` with odd parity, the necessary condition `A·D ∈ J(F)` forces `D ∈ J(F)`; the six partials `∂F/∂Z_i` (degree 5) are then inconsistent as denominators. Hence **every certificate has `deg D ≥ 9`, or has `D ∈ J(F)` of degree 6, 7 or 8 (even) — the latter family is being swept right now.** The book's own analysis (§4.7) shows the `τ`-invariant sector of *any* even sextic is the hard one: `X/τ` is not uniruled, the natural cubic cover carries the class in its standard `S_3` representation, and no quotient trick serves it. Only a genuinely new correspondence or a (probably degree `≥ 10`) denominator found by structured search can close this instance.

---

## 1. Exact target

### 1.1 Statement and certificate shape

```
q = 1,  A = X0^3*Y0^3,  F = the 236-term sextic of instance.rec (sparse presentation in 1.2)
```

Certificate (formal contract `formal/README.md`, book ch. 1): `D ≠ 0` homogeneous with `gcd(D,F) = 1` over `C`, numerators `N_i` homogeneous of degree `deg D + 6q − 5 = d + 1`, and the exact identity

```
E := A·D^2 − q·D·Σ_i N_i ∂_iF + F·( D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD ) = F^(q+1)·Q = F^2·Q .
```

Degrees: `deg E = 2d + 6`, `deg F^2 = 12`, so `deg Q = 2d − 6` (`Q = 0` forced for `d ≤ 2`; `d = 3` gives constant `Q`). Normalized form (ch. 1, used by every accepted solution): since `F` is irreducible and `gcd(D,F)=1`,

```
E1:  N(F) := Σ N_i ∂_iF ≡ A·D   (mod F)   ⇒  A·D ∈ J(F)     [F ∈ J by Euler]
E2:  D·div N − N(D) ≡ 0        (mod F)
```

No denominator-degree bound and no extra parity condition is imposed by the contract.

**Parity reductions (rigorous, elementary).** `F` is `τ`-even, `A` is `τ`-odd, `J(F)` is `τ`-stable (`∂F/∂X_i` even, `∂F/∂Y_i` odd).
* If `(D,N,Q)` is a certificate then so is `(D·τD, N·τD, Q·(τD)^2)` (a rational field `v = N/D` does not change when numerator and denominator are multiplied by the same form; `gcd(τD,F) = 1` since `F` is `τ`-invariant). Hence **WLOG `D` is `τ`-even**, i.e. `D` is a polynomial in `X_0,X_1,X_2` and the six products `Y_aY_b`; equivalently `V(D)` descends to the quotient `X/τ`.
* For fixed `D` of pure parity the identity is linear in `(N,Q)`; applying `τ` and averaging shows one may take `N_{X_i}` of parity `par(A)+par(D)` and `N_{Y_i}` of parity `par(A)+par(D)+1`, `Q` of parity `par(A)`. With `D` even and `A` odd: `N_{X_i}` odd, `N_{Y_i}` even, `Q` odd. (This is what `hodge.fixed_d.numerator_parities` implements.)
* `E1` splits by parity of `D`, so a parity-split E1-kernel computation covers all `D`.

Linear-system sizes for the full fixed-`D` system (unknowns `3·S_{d+1}` numerator coefficients after parity, `|S_{2d−6}^{odd}|` quotient coefficients; equations = odd monomials of degree `2d+6`), and for the E1-kernel matrix (columns: `D`-monomials of one parity + numerators; rows: monomials of degree `d+6`):

| `d` | fixed-`D` unknowns (N + Q) | fixed-`D` rows | E1-kernel cols / rows |
|---:|---:|---:|---:|
| 4 | 756 + 9 = 765 | 5 796 | 822 / 1 491 |
| 5 | 1 386 + 60 = 1 446 | 10 152 | 1 512 / 2 184 |
| 6 | 2 376 + 226 = 2 602 | 16 797 | 2 612 / 3 080 |
| 7 | 3 861 + 636 = 4 497 | 26 532 | 4 257 / 4 284 |
| 8 | 6 006 + 1 491 = 7 497 | 40 326 | 6 657 / 5 796 |
| 9 | 9 009 + 3 080 = 12 089 | 59 332 | 10 010 / 7 752 |
| 10 | 13 104 + 5 796 = 18 900 | 84 903 | 14 616 / 10 152 |
| 11 | 18 564 + 10 152 = 28 716 | 118 608 | 20 748 / 13 167 |
| 12 | 25 704 + 16 797 = 42 501 | 162 248 | 28 812 / 16 797 |
| 14 | 46 512 + 40 326 = 86 838 | 287 793 | 52 344 / 26 532 |

(The local engine solves the `d = 5` system in ≈50 s and the `d = 6` one in a few minutes, mod a 31-bit prime with random row compression, then lifts rationally if consistent.)

### 1.2 Structure and sparse presentation (all identities verified exactly with sympy, 2026-11-01)

```
F = Σ_{i=0}^{2} Q_i · K_i ,          Q_i = X_i^2 + Y_i^2 ,
K_i = X_i^4 − X_i^2 Y_i^2 + Y_i^4 + (5/p_i) · H_i ,     (p_0,p_1,p_2) = (7, 11, 13),
```

where `H_0, H_1, H_2` are `Y`-even quartics in **all six** variables with 66, 54 and 43 monomials respectively (explicit expansions in `/tmp/hodge/results/eight-plane-configuration/K_decomp.json`, keys `K` and `H`; `H_1` has no pure-`X0` monomials and `H_2` none in `X0,X1` alone, reflecting a triangular elimination order — the `K_i` are only defined modulo the Koszul syzygies `K_i ↦ K_i + Σ_j c_{ij} Q_j`, `c_{ij} = −c_{ji}` quadrics). Since `Q_i·(X_i^4 − X_i^2Y_i^2 + Y_i^4) = X_i^6 + Y_i^6`,

```
F = (X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6) + 5 · Σ_i Q_i H_i / p_i  =  Fermat + 5·G,
```

with `G` having 5-integral coefficients (verified). All denominators of coefficients of `F` are products of powers of 7, 11, 13 (e.g. `31213 = 7^4·13`, `343343 = 7^4·11·13`, `99127 = 7^4·…`). Monomial support: **all 236** `Y`-even sextic monomials (`Σ_{k even} C(k+2,2)·C(8−k,2) = 28+90+90+28`), coefficients all nonzero. There is no block/split structure of any kind (no monomial class is missing; the `H_i` couple every pair of blocks).

**Eight planes.** `V(Q_0,Q_1,Q_2)` is the reduced union of the eight planes `Π_ε : Y_i = ε_i·i·X_i` (`ε ∈ {±1}^3`), and `(Q_0,Q_1,Q_2)` is exactly their ideal (complete intersection of three products of distinct linear forms). `F ∈ (Q_0,Q_1,Q_2)` ⇒ `Π_ε ⊂ X` (asserted and checked symbolically in `audit-structural-frontier.py`; re-verified). Incidence of the arrangement: `Π_ε ∩ Π_{ε'}` is a line if `ε, ε'` differ in one sign (12 lines), the point `[e_i + ε_i·i·e_{i+3}]` if they differ in two signs (6 points, each on four planes: exactly the six singular points of `V(Q_0,Q_1)`-type pencils noted in the thread), empty if `ε' = −ε` (4 antipodal pairs of disjoint planes, each pair spanning `P^5`). `τ(Π_ε) = Π_{−ε}`: the involution pairs the planes; no plane is individually `τ`-stable. Over `Q(i)` each `Q_i` splits as `(X_i + iY_i)(X_i − iY_i)`, so `X ∩ V(Q_i)` is the union of two conjugate hyperplane sections each containing four of the planes (hence each singular, in at least finitely many points, since a smooth threefold hypersurface contains no plane by Lefschetz).

**Why "contact arrangement".** The title refers to `F` lying in the ideal of the arrangement; `why.md` is explicit that this does *not* give a contact presentation in the book's sense (`F = Γ(h_1..h_4) + P·K` with first integrals of degrees summing to 6, or `F = c·H^m + D·K`), and that "writing `F` as a sum of products `Q_i·H_i` does not identify `A` as the density of the book's contact constructions".

### 1.3 Symmetry and the character of the target

* **Diagonal symmetries.** Smith normal form of the lattice spanned by differences of the 236 exponent vectors: invariant factors `(1,1,1,1,2)`, rank 5. Hence the group of diagonal `t ∈ (C^*)^6` with `F(tZ) = λF`, modulo scalars, is exactly `Z/2 = ⟨τ⟩`. (Consistent with `checks.json`: `coordinate_reflections = []`, `coordinate_paired_signs = []`; no single coordinate sign and no triple of paired signs preserves `F`.)
* **Non-diagonal symmetries tested (all negative):** the five nontrivial simultaneous block permutations `(X_i,Y_i) ↦ (X_{σi},Y_{σi})`; the global swap `X ↔ Y`; the three block swaps `X_i ↔ Y_i`; the three block rotations `(X_i,Y_i) ↦ (−Y_i, X_i)`. (The coefficients `5/7, 5/11, 5/13` on `H_0,H_1,H_2` and the unequal pure powers `31218/31213·X0^6`, `49589/49049·X1^6`, `74842/57967·X2^6`, … break them.) Nobody has computed the full linear automorphism group; I expect it to be `⟨τ⟩`.
* `checks.json`: `hessian_symmetrizer_dimension_mod101 = 1`, i.e. the space of `6×6` matrices `M` with `M·Hess(F)` symmetric (as polynomial matrices) is one-dimensional (scalars) mod 101 — the audit's test that `F` has no hidden linear direct-sum decomposition (a direct sum `F = F_1(U) + F_2(V)` would give a 2-dimensional symmetrizer space).
* **Character of the target.** `A = X0^3 Y0^3` has `Y`-degree 3 (odd) and the projective volume form `Ω = Σ ± Z_j dZ_0∧…` is `τ`-odd (`det τ = −1`), so the residue `Res(A·Ω/F^2) ∈ H^{3,1}(X)` is **`τ`-invariant**. In Jacobian-ring terms: `H^{3,1}(X) ≅ R_6` (dim 426); `τ`-invariant classes ↔ odd numerators `R_6^- = S_6^-/J_6^-` of dimension `226 − 18 = 208`; anti-invariant ↔ `R_6^+`, dimension `236 − 18 = 218` (the book §4.4 uses exactly these numbers: `(dim G, dim V) = (218, 208)`). The target lies in the 208-dimensional invariant sector — the sector that §4.7 of the book singles out as transcendental for every even sextic (see §2.4 below).
* **Nonvanishing of the class (rigorous).** `A ∉ J(F)`: (a) mod 5, `F ≡ Fermat` and `J ≡ (Z_0^5,…,Z_5^5)`, and `A` has all exponents `< 5`; (b) directly, the local test `D = 1` gives `rank(J_6^- basis) = 18 < rank(augmented) = 19`. So the target class is a nonzero primitive `(3,1)` class.

### 1.4 Smoothness and primitivity (from `checks.json`, method inspected in the audit script)

All coefficients are 5-integral and `F ≡ X0^6 + … + Y2^6 (mod 5)`; the Fermat sextic is smooth in characteristic 5 (`6 ≢ 0`), so by openness of smoothness in the flat `Z_(5)`-family, `X = V(F)` is smooth over `Q` (hence over `C`). "Primitive target": `A` survives in the Jacobian ring mod 5, and smooth complete-intersection base change transports nonvanishing to characteristic zero.

### 1.5 Jacobian-ring numerics (exact)

`dim R_d` for `d = 0..24`: `1, 6, 21, 56, 126, 246, 426, 666, 951, 1246, 1506, 1686, 1751, 1686, 1506, 1246, 951, 666, 426, 246, 126, 56, 21, 6, 1`. For `d ≤ 9`, `J_d ≅ S_{d−5}^{⊕6}` (no syzygies), `dim J_d^+ = dim J_d^- = 3·|S_{d−5}|` (= 3, 18, 63, 168, 378 for `d = 5..9`). First Koszul syzygies in degree 10: `dim J_10 = 1497`, `J_10^+ = 750`, `J_10^- = 747`; `R_10^+ = 762`, `R_10^- = 744`; `R_8^+ = 483`, `R_8^- = 468`.

Consequence for E1 (rigorous, general): `{D ∈ S_d : A·D ∈ J} / J_d = Ann_R(A)_d = ker(·A : R_d → R_{d+6})`. For `d ≤ 8`, `dim R_d < dim R_{d+6}`, so the kernel is zero *iff* multiplication by `A` has maximal rank — this is what the computation in §3.3 confirms for this `F`. For `d = 9`, `R_9 → R_15` is a square map (`1246 × 1246`), self-dual under the Gorenstein pairing; generic nondegeneracy expected, **not yet verified for this `F`** (in flight). For `d = 10`: `R_10 → R_16 ≅ R_8^∨` cannot be injective; `dim Ann(A)_{10}^{even} ≥ R_10^+ − R_8^- = 762 − 468 = 294` and `dim Ann(A)_{10}^{odd} ≥ 744 − 483 = 261`. So **degree 10 is the first degree where non-Jacobian denominators can pass E1, and for every smooth sextic with this `A` there is a ≥ 294-dimensional even E1-kernel there.** E1 stops being a filter at `d = 10`; the real obstruction is E2.

### 1.6 Geometry recorded in the thread (rigorous, `benediamond`)

* Normal bundle of a plane: `0 → N_{Π/X} → O_{P^2}(1)^3 → O_{P^2}(6) → 0`; `H^0(O(1)^3) = 9 → H^0(O(6)) = 28` is injective (shown mod 5 where the derivative of a normal displacement `Y_i = ε_i i X_i + ηL_i(X)` is `Σ 6(ε_i i)^5 X_i^5 L_i(X)`, nine distinct monomials with nonzero coefficients). Hence `h^0(N_{Π/X}) = 0`, `h^1 = 19`: **each plane is infinitesimally rigid**; moving the planes cannot sweep a divisor.
* Projection from a plane `Π` (write `F = Σ u_i G_i`, `u_i` the linear forms cutting `Π`): residual fibres are quintic surfaces `S ⊂ P^3` with `K_S = O_S(1)`, `p_g = 4`, `q = 0`, smooth for general fibre (Bertini along `λ = 0` because `G_i|_Π` have no common zero by smoothness of `X`). `τ` exchanges the projections from `Π` and `τΠ` rather than acting fibrewise.
* The invariant net `φ = [Q_0:Q_1:Q_2] : X ⇢ P^2`, base locus = the eight planes. General fibre `S_r = {r_1Q_0 − r_0Q_1 = r_2Q_0 − r_0Q_2 = 0, Σ r_iK_i = 0}`, a smooth `(2,2,4)` complete intersection in `P^5` (smoothness argued for this actual family: the six singular points of the `(2,2)` part are avoided because `K_i ≠ 0` there by smoothness of `X`, Bertini on the plane/line strata, generic smoothness elsewhere). `K_{S_r} = O(2)`, `p_g = 21 − 2 = 19`, `q = 0`. Parity bookkeeping including the odd volume form: two-forms with odd quadratic numerator `P` (the 9 monomials `X_iY_j`, no relations) are `τ`-invariant, those with even `P` (12 monomials, 2 relations) are anti-invariant: `h^{2,0}(S_r)^+ = 9`, `h^{2,0}(S_r)^- = 10`. (I re-checked all of these adjunction/dimension counts.)

---

## 2. Companion and related solved instances, and why their methods fail here

`provenance.json` lists no `sameVarietyInstances`; the only references are the two book chapters. The relevant comparison classes are (a) the book's contact constructions that `why.md` names, (b) the accepted instances with the *same* target `X0^3·Y0^3`, (c) the general theorems of ch. 3/4.

### 2.1 Contact-type accepted instances (`why.md`'s "book's contact constructions")

| instance | `F` shape | method | `D` (degree, polynomial) | why not here |
|---|---|---|---|---|
| `contact46006` | polarized plane sextic `f(U)+f(V)`, `U = X+Y, V = X−Y`, with `f = h^2 + C·k`, `h = z_0^3+2z_2^3`, `C` a plane cubic | Hamiltonian/cross-product field `W = ∇C × ∇h` (`div W = 0`, `W(C) = W(h) = 0`), `W(f) = C·W(k)`; polarize | `d = 6`, `D = C(U)·C(V)` | needs a perfect-square contact `f = h^2 + C·k` in a 3-variable half and a split `F = f(U)+f(V)`; `F` here is dense with no polarization structure; target is dictated by the presentation (`A = W(k)`), not prescribed |
| `polarized-contact-coupled-d2` / `-deformation-d2` | `f = r^2 + z_2·Q` in each half plus coupling `c·r(U)r(V)` (and `+ s·D·H`) | same Hamiltonian field, coupling killed because `W(r) = 0`; deformation by multiples of the preserved divisor | `d = 2`, `D = U_2 V_2 = X2^2 − Y2^2` | requires common first integrals preserved by the coupling; the `H_i` here are not first integrals of anything divergence-free that also preserves `A`'s presentation |
| `family-passive-contact` | `X0^6 + X1^6 + X0^2X1^2 R(Y) + passive sixth powers` | linear-contact (Hesse singular value) denominator `D_L = X0^6+X1^6+X0^2X1^2R − R^3/27` (product of three conjugate quadrics) | `d = 6` | needs an active binary block in squared coordinates plus passive variables; here nothing is passive |
| `completed-square-nambu-holdout` | `X1^6 + T^2 + G(others)` | weighted radial field in `(X1^3, T)`, `N(G) = 0` | `d = 6`, `D = X1^6 + T^2` | needs a completed square separating two active variables from the rest |

General statement from the book (§2.2): a Nambu contact field `W/P` with `F = Γ(h_1..h_4) + P·K`, `Σ(deg h_i − 1) = 2`, serves exactly `A = W(K)` — "the contact selects a specific target through a residue pairing; it does not serve an arbitrary prescribed `A`". For `F = Σ Q_iK_i` one would need, e.g. with `P = Q_0`, `Q_1K_1 + Q_2K_2 = Γ(h_1,…,h_4)` with `deg h = (1,1,2,2)` or `(1,1,1,3)` — impossible for a form that is dense in the complementary variables. The only divergence-free fields naturally attached to the arrangement are the block rotations `R_i = X_i∂_{Y_i} − Y_i∂_{X_i}` (`div R_i = 0`, `R_i(Q_j) = 0` for all `j`), which give `R_i(F) = Σ_j Q_j R_i(K_j) ∈ (Q_0,Q_1,Q_2)`; on Fermat (`H = 0`) they serve `a_0·X_0Y_0(X_0^2 − Y_0^2)` with `D = Q_0`, never `X0^3Y0^3`; here they merely produce the Jacobian denominators `R_i(F) ∈ J_6^-` listed in §3.3.

### 2.2 Accepted instances with the identical target `A = X0^3·Y0^3`

All are **sparse** even sextics with Fermat diagonal and a few monomial couplings, or binary-block sextics:

| instance | `F` | method | `D` |
|---|---|---|---|
| `family-block-cube` | `X0^6 + X0^4Y0^2 + 2X0^2Y0^4 + 3Y0^6 + (X1^6 + X1^3X2^3 + 2X2^6) + Y1^6 + Y2^6` | binary cubic `p(u,w)` in `u = X0^2, w = Y0^2`; **osculating cube** `C = (2X0^2 + 9Y0^2)^3/243` matching `p` to first order; spectator sixth powers with forced weights | `d = 6`, `D = C + Y1^6 + (1681/1701)Y2^6` |
| `family-hesse-web` | `X0^6 − 2X0^2X1^2Y0^2 + X1^6 + 2X2^6 + Y0^6 + 3Y1^6 + 5Y2^6` | Hesse cubic in `(X0^2, Y0^2, X1^2)`; three contact-conic channels over the cubic algebra `C^3 − tC^2 − 4 = 0`; norm denominator | `d = 6`, `D = Norm(a + C·x) = a^3 + t a^2 x + 4x^3`, `a = X0^2+Y0^2`, `x = X1^2` |
| `holdout`, `holdout-moving-branch`, `holdout-ternary-deformation`, `cycle-three`, `mixed-biweight`, `one-reflection` | Fermat + 2–3 monomial couplings (`X0^3X1Y0Y1`, …) | nonconstructive-constructive hybrids: Fermat-curve cover on which a reflection appears, moving polynomial conic, polynomial-family integration, finite trace, deck character (`HoldoutCover`, `PolynomialCurves` Lean libraries) | existential (`HasCertificate`), no expanded `D` |
| Fermat itself | `Σ Z_i^6` | binary block `x^6+y^6 = (x^2+y^2)(x^4−x^2y^2+y^4)`; Diamond's paper (arXiv 2603.29157) treats all Fermat targets | block-cofactor `d = 4` for `x^3 y·B`; block-cube type for `x^3y^3` |

Why they fail here: every one of them uses either (i) a binary block in `(X0,Y0)` (so that `F = p(X0^2,Y0^2) + (rest)` or `F = block + passive`), (ii) a coordinate reflection `Y_j ↦ −Y_j` (Reflection theorem §3.3 — requires `F` even in a *single* coordinate), (iii) paired signs `(X_i,Y_i) ↦ (−X_i,−Y_i)` (Paired-sign sector theorem §4.2), or (iv) a Fermat cover made possible by sparse monomial couplings. `F` here has **none**: `coordinate_reflections = []`, `coordinate_paired_signs = []`, all 236 monomials present, and the `H_i` make every pair of blocks interact in every monomial class. `why.md` says this explicitly: "those mixed quartics remove the separate coordinate reflections while preserving the entire arrangement … rather than merely adding an unrelated coefficient to a solved example."

### 2.3 Deformation relation to Fermat (observation, not in the thread)

`F_t = Fermat + t·Σ_i Q_iH_i/p_i`, `t ∈ [0,5]`, is a one-parameter family of `τ`-even sextics all containing the eight planes (`F_t ∈ (Q_0,Q_1,Q_2)` for all `t`), with `F_5 = F` and `F_0 = Fermat` where the target *is* served with a degree-6 denominator. But the local E1 computation (§3.3) shows that at `t = 5` **no even degree-6 denominator outside `J(F)_6` passes E1**, so the Fermat degree-6 witness cannot be continued inside degree 6 to `t = 5` (unless it degenerates into a Jacobian denominator). Any continuation argument would have to let `deg D` grow (to ≥ 9/10) along the path.

### 2.4 Book theorems that do *not* apply (so the instance is genuinely at the frontier)

* §3.3 Reflection theorem / reflection-generated groups: need a reflection in `Aut(F)`; here `Aut` (diagonal) `= ⟨τ⟩` and `τ` is not a reflection (fixed locus codim 3). The book: "The global involution has two fixed plane sextic curves … It supplies no branch divisor or Fano quotient".
* §3.3 Kollár–Larsen age criterion: `τ` has tangent eigenvalues `(1,−1,−1,−1)` at the fixed curves, age `3/2 ≥ 1`; §4.7 proves `Q = X/τ` has torsion `K` and is **not uniruled** (Kodaira dimension 0) — so the quotient-transfer + covering-curves route is closed for the `τ`-invariant sector.
* §4.7 natural cubic cover: `Q → B = P^2_X × P^2_Y`, cubic `a + bt + ct^2 + dt^3` with `(a,b,c,d) = (F_{6,0},F_{4,2},F_{2,4},F_{0,6})`, discriminant `Δ ∈ H^0(O_B(12,12))`, monodromy `S_3`; the invariant primitive `H^4` sits in the standard representation `V_std`, and no intermediate quotient of the Galois closure is uniruled with `V_std^K ≠ 0`. Applies verbatim to this `F` (dense even sextic ⇒ general-position hypotheses are plausible but not checked for this `F`).
* §3.2 Surface-fibration theorem (`p_g = q = 0` fibres): the two natural fibrations have `p_g = 4` (quintic residuals) and `p_g^+ = 9` (invariant net) — `benediamond`'s comment is precisely the check that this theorem does not apply directly.
* §4.4 "one generic nonzero target suffices": needs a *dominating* family of supports over the 236-dimensional even-sextic space; an isolated instance gives nothing.
* §3a split/K3 theorem: needs `F = f(X) + g(Y)`-type polarized splitting; absent.

---

## 3. Prior work, comment by comment

### 3.1 Discussion #176 (the entire public record)

**Comment `DC_kwDOVLOOJc4BJDZL`, `benediamond`, 2026-10-23T17:01Z (edited 17:04Z), ChatGPT-assisted. No replies as of 2026-11-01.** Title: "Rigidity of the eight planes, and the Hodge types of two natural surface fibrations". Content (all three items are proofs, reproduced in §1.6):

1. *Rigidity of each plane* — `h^0(N_{Π/X}) = 0`, `h^1 = 19`, via the mod-5 Fermat reduction. Rigorous. Conclusion drawn: "moving these planes themselves does not produce a supporting divisor. This leaves open other Hilbert-scheme components and families of curves or surfaces related to the arrangement."
2. *Projection from one plane* — quintic residual surfaces, `K_S = O_S(1)`, `p_g = 4`, `q = 0`; `τ` swaps the two projections. Rigorous. Conclusion: the whole-fibre `p_g = q = 0` criterion does not apply and there is no immediate fibrewise `τ`-projector repair.
3. *Invariant net* `[Q_0:Q_1:Q_2]` — smooth `(2,2,4)` fibres for general `r`, `p_g = 19`, `q = 0`, `h^{2,0}(S_r)^+ = 9`, `h^{2,0}(S_r)^- = 10` (with the volume-form sign). Rigorous. Conclusion: "the entire invariant relative `H^2` of this natural net is not of type `(1,1)`. A criterion requiring a wholly algebraic invariant fiber summand cannot directly establish support here."

Explicitly stated limits: "They do **not** show that the particular target is unsupported … The remaining question is more specific: does the fixed absolute target avoid, cancel in, or admit a further algebraic projector away from this nine-dimensional relative holomorphic piece? Nothing in these calculations decides that."

Denominators tested in the thread: **none**. Primes/ranks: none. Reproduction packets/URLs: none beyond the link to `why.md`. Conjectures: none stated. Open sub-problem posed: whether the target class can be separated from the 9-dimensional invariant `H^{2,0}` of the net fibres by a further correspondence.

### 3.2 `why.md` (instance rationale, not a discussion post)

Claims: eight planes exist; `F ∈ (Q_i)`; `X_i^6+Y_i^6 = Q_i(X_i^4 − X_i^2Y_i^2 + Y_i^4)`; additional couplings `Q_i·H_i`; 236 monomials; scalar Hessian symmetrizer excludes hidden direct sums; planes have type `(2,2)` so exhibiting them cannot serve the `(3,1)` target; "writing `F` as a sum of products `Q_i·H_i` does not identify `A` as the density of the book's contact constructions". Suggests as the geometric hope: "whether a family moving away from the arrangement sweeps a divisor with the required Gysin image" — which `benediamond`'s item 1 then answers negatively for first-order motions of the planes themselves.

### 3.3 Unpublished local computations (this workspace, 2026-11-01; NOT in the discussion)

Engine: `~/devel/sextics-work/hodge/{core.py, fixed_d.py}` (python-flint, prime `p = 2^31 − 1 = 2147483647`, random row compression for the big systems). A sign bug in the `F·(D div N − N(D))` term was fixed at 17:25; everything below marked "engine 2" postdates the fix; the E1-kernel code (`e1kernel.py`) never used the buggy routine.

* `D = 1`: `A ∉ J_6` (rank 18 vs augmented 19). Rigorous mod `p`, hence over `Q`.
* **E1-kernel sweep** (`e1kernel.jsonl`): `dim{D ∈ S_d^{par} : A·D ∈ J}` computed mod `p` for `d = 1..8` (even parity) and `d = 1..7` (odd parity):

  | `d` | parity | `D`-cols | nullity of full matrix | projected kernel dim | `dim J_d^{par}` | kernel mod `J` |
  |---:|:-:|---:|---:|---:|---:|---:|
  | 1–3 | both | 3–28 | 0 | 0 | 0 | 0 |
  | 4 | even / odd | 66 / 60 | 9 / 6 | 0 / 0 | 0 | 0 |
  | 5 | even / odd | 126 / 126 | 48 / 48 | 3 / 3 | 3 / 3 | 0 |
  | 6 | even / odd | 236 / 226 | 180 / 171 | 18 / 18 | 18 / 18 | 0 |
  | 7 | even / odd | 396 / 396 | 483 / 483 | 63 / 63 | 63 / 63 | 0 |
  | 8 | even | 651 | 1122 | 168 | 168 | 0 |

  (Nullity minus projected kernel is the Koszul freedom in `N`.) **Interpretation (rigorous over `Q`):** `dim ker_Q ≤ dim ker_{F_p}` (rank can only drop mod `p`, and `p` divides no coefficient denominator — they are products of 7, 11, 13), `J_d(Q) ⊆ ker_Q`, and `dim J_d(F_p) = dim J_d(Q) = 3|S_{d−5}|` (the computed values equal the syzygy-free count). Hence `ker_Q = J_d(Q)` exactly: **for even `D` with `d ≤ 8` and odd `D` with `d ≤ 7`, E1 forces `D ∈ J(F)`.** Equivalently `Ann_R(A)_d = 0` in those degrees, i.e. multiplication by `A` is injective `R_d → R_{d+6}`.
* **Fixed-`D` solves, engine 2** (`log.jsonl`): the six partials `D = ∂F/∂Z_i` (`d = 5`, `D ∈ J_5`, E1 trivially solvable): all six **inconsistent**; system `10152 × 1446`, rank 1318, nullity 128, ≈50 s each.
* **In flight at the time of writing** (launched by the parallel attack agent and by me; results, if finished, in the addendum): (a) `run_J6.py` — the fixed-`D` system for ten natural elements of `J_6`: block Euler derivatives `E_i(F) = X_i F_{X_i} + Y_i F_{Y_i}` (`i = 0,1,2`), `X`-Euler `Σ X_iF_{X_i}`, monomial Jacobian elements `X_0F_{X_0}`, `Y_0F_{Y_0}`, `X_0F_{Y_0}`, `Y_0F_{X_0}`, the `X`-rotation `X_1F_{X_0} − X_0F_{X_1}`, and the block rotation `R_0(F) = Y_0F_{X_0} − X_0F_{Y_0}` (odd; its zero locus contains all eight planes since `R_0(Q_j) = 0`); (b) `e1kernel.py 9` — E1-kernel at `d = 9`, both parities (`7752 × 10010`); (c) E1-only checks of `Q_0Q_1Q_2`, `Q_0^3`, `Q_0^2Q_1` (predicted `false` by the kernel theorem unless they lie in `J_6`); (d) the attack agent's `F − Q_0K_0` (`≡ −Q_0K_0 mod F`, i.e. the reducible divisor `(X∩V(Q_0)) ∪ (X∩V(K_0))`): **E1 false** (rank 2214 vs 2215), as predicted.
* Prepared but not yet run: `koszul_e2.py` — for `D ∈ J` with particular E1 solution `N_0`, decide E2 over `N = N_0 + Σ c_{ij}(F_j∂_i − F_i∂_j)` (`deg c_{ij} = d − 4`) and `Q` — the reduced system for sweeping whole `J_d` families.

Rigorous exclusions vs screened failures: the E1-kernel results are rigorous exclusions of *all* denominators of the stated degrees/parities outside `J(F)`; the partials are rigorous exclusions of six specific denominators (modular inconsistency ⇒ rational inconsistency, since an exact rational solution would reduce mod `p`). Everything else is untested, not excluded.

---

## 4. What remains unexplored (concrete gaps)

1. **E1-kernel at `d = 8` odd, `d = 9` (both), `d = 10` (both).** `d = 8` odd and `d = 9` are cheap (minutes–hour). At `d = 10` the even kernel mod `J` is provably ≥ 294-dimensional; computing it exactly (and a basis) is the entry point for all non-Jacobian denominator searches. Nothing of this has been done for this `F`.
2. **E2 over Jacobian denominators `D ∈ J_6, J_7, J_8`.** Only the six partials (`J_5`) have been tested. `J_6^±` are 18-dimensional each (`D = c(F)`, `c` a linear vector field); the ten natural ones are in flight; a *full* sweep of the 18-dimensional even family is a bilinear problem never attempted (see §5, E3). `J_7`, `J_8` (63-, 168-dimensional) untouched. Geometric meaning: `X ∩ V(c(F))` is the tangency locus of the field `c` on `X`; for `c = R_0` it contains the eight planes.
3. **Even denominators in the ideal of the arrangement, degree ≥ 9/10:** `Q_0Q_1Q_2·G_4`, `Q_i^a Q_j^b · G`, `g_5(Q_0,Q_1,Q_2)` (pull-backs of plane quintics under the net), … — nobody has intersected these structured families with the E1 kernel.
4. **Denominators from the fibrations:** (a) `φ^*` of the discriminant curve of the `(2,2,4)` net (degree `2·deg Δ_net`, computable by an Euler-characteristic count — never computed); (b) the degree-24 even form `Δ(F_{6,0},F_{4,2},F_{2,4},F_{0,6})(X,Y)`, the pull-back of the cubic-cover discriminant of §4.7; (c) for the quintic-surface projections from `Π_ε ∪ Π_{−ε}`: the discriminant of the pencil/net of residual quintics. None computed, none tested.
5. **Hodge-theoretic separation question left by `benediamond`:** compute the image of the target class under the Leray decomposition of the net `φ` (does `α` restrict to zero in `H^0(P^2, R^4)`? does it have a nonzero component in `H^2(P^2, R^2φ_*)` landing in the 9-dimensional `H^{2,0}(S_r)^+` Hodge bundle via Griffiths transversality?). Decidable in principle by infinitesimal computation (the `(3,1)` class and the Kodaira–Spencer map of the fibration), never attempted.
6. **Other Hilbert-scheme components:** conics/lines/curves in the planes (each plane contains a 2-dimensional family of lines and 5-dimensional family of conics, far above the expected dimensions 1 and 1 for a sextic fourfold) and their deformations off the planes; surfaces through the 12 lines / 6 quadruple points; the cones over the planes. `benediamond` explicitly leaves this open; nobody has looked.
7. **Nonlinear joint searches** (unknown `D` in a structured family × linear `(N,Q)`): none run.
8. **Continuation from Fermat along `F_t`** (§2.3): never attempted; the degree-6 obstruction at `t = 5` is known (this dossier), the behaviour of `deg D` along the path is not.
9. **Automorphism group over `C`:** only diagonal and coordinate-permutation symmetries checked.
10. **Reduction to the quotient `X/τ`:** since WLOG `D` is even, one could set the whole search up on `Q = X/τ ⊂ P(1,1,1,2,2,2,2,2,2)`-Veronese model (fewer monomials); never done.

---

## 5. Candidate attack plan

Each item: idea · exact family · system size · chance (my estimate) and reasoning.

**E1. Finish the E1-kernel ladder: `d = 8` odd, `d = 9`, `d = 10` (even), with bases.** Linear algebra mod `p`, matrices `5796×6636`, `7752×10010`, `10152×14616` (dense flint rref: minutes, ≈30–60 min, ≈2–4 h; or compressed/sparse). Output: at `d = 9` either `0` (expected generically) or a small kernel (which would be a structural surprise worth testing immediately); at `d = 10` a basis of the ≥294-dimensional even `Ann_R(A)_{10}` (lifted to `S_10` modulo `J_10`). Chance of directly producing a witness: ≈0; necessity for everything below: total.

**E2. Koszul-reduced E2 over all of `J_6^±` by structured enumeration.** `D = c(F)` with `c = Σ a_{ij} X_j∂_{X_i} + b_{ij} Y_j∂_{Y_i}` (even, 18 parameters) or `c = Σ a_{ij} Y_j∂_{X_i} + b_{ij} X_j∂_{Y_i}` (odd, 18). E1 has the particular solution `N_0 = A·c` (since `A·c(F) = (Ac)(F)`), E2 is linear in the Koszul coefficients `c_{ij} ∈ S_2` (15·21 = 315 unknowns before parity) and `Q ∈ S_6^{odd}` (226): ≈540 unknowns, rows `|S_12^{odd}| = 3080` — seconds per `D`. The consistency locus in `P^17` is algebraic; probe it on (i) the ten geometric candidates in flight, (ii) all `c` with `≤ 2` nonzero entries, (iii) random `c` (to confirm generic inconsistency), (iv) a random pencil `c_0 + s c_1` with `s` symbolic mod `p` (polynomial-matrix rank drop) to detect whether the consistency locus is nonempty at all. Chance a `J_6` denominator works: low (≈5%) — the tangency threefold `X ∩ V(c(F))` of a linear field is, for generic `c`, a `(6,6)` complete intersection with isolated singularities (normal, rational ⇒ `H^{2,0}` of resolution `= 0`, so by §3.1 it cannot support a `(3,1)` class); only special `c` with non-normal / non-rational-singular tangency loci can work, and `R_0` (tangency locus ⊃ eight planes) is the one geometrically motivated candidate.

**E3. Same for `J_7`, `J_8`** (`c` quadratic/cubic vector fields; 63/168 dims per parity): per-`D` Koszul-E2 systems of size ≈ `(15·|S_3|+|S_8^{odd}|) ≈ 1476` unknowns × `|S_14^{odd}|=8 398` rows, and ≈ `(15·|S_4| + |S_10^{odd}|) ≈ 3 381` × `13 167`. Candidates: `c = h·R_0` with `h` a quadric (tangency loci containing the planes), `c = R_0 + R_1`, `E_0 − E_1`-type, `c = ∇(Q_i)`-Hamiltonians. Chance: low (≈5–10%), same geometric reason.

**E4. Plane-ideal denominators of degree 10 inside the E1 kernel.** Families: `D = Q_0Q_1Q_2·G`, `G ∈ S_4^{+}` (66 params); `D = Q_i^2 Q_j Q_k·G_2`; `D = Q_0^5`-type powers times quartics; `D = g_5(Q_0,Q_1,Q_2)` (21 params). Step 1: intersect each family with `Ann(A)_{10}` — a linear system in the family's parameters plus `N` (`66 + 3·|S_11| = 66 + 13104` columns, `10152` rows; minutes). Step 2: if the intersection `W` is nonzero and small-dimensional, solve E2 on `W` (bilinear; if `dim W ≤ 3`, parametrize and use a polynomial-matrix consistency test; if `dim W = 1`, a single fixed-`D` solve of size `18 900 × 84 903`, ≈10–30 min compressed). Chance: low-moderate (≈10%). Motivation: `X ∩ V(Q_0Q_1Q_2)` is a union of six hyperplane sections, each singular along/at the planes it contains; a quartic factor `G` can add a component with non-normal singularities along the 12 lines. Against: §3.1 says unions of hyperplane sections do not support `(3,1)` classes; the whole burden falls on `V(G)`.

**E5. Leray/Griffiths computation for the invariant net (decides whether any `φ^*`-pulled-back denominator can work).** For a curve `C ⊂ P^2_r` of degree `k`, `D = g_k(Q_0,Q_1,Q_2)` (degree `2k`, even). The class `α` is supported on `φ^{-1}(C)` iff its restriction to `H^4(φ^{-1}(P^2∖C))` vanishes. Compute: (i) the image of `α` in `H^0(P^2, R^4) = Q` (zero: `α` is primitive, restricts to zero on the generic fibre for degree reasons); (ii) the component of `α` in `Gr_F^3 H^2(P^2, R^2φ_*Q)` — via the Kodaira–Spencer/Higgs map of the family `S_r`, i.e. whether `α` pairs nontrivially with `H^1(P^2, Ω^1 ⊗ F^2)`-type classes; concretely, test whether the Gysin images of the 9 invariant two-forms of `S_r`, moved along the base, span a subspace of `H^{3,1}(X)^+` containing `α`. This is a finite Jacobian-ring computation (restrict `A·Ω/F^2` to the pencil of threefolds `φ^{-1}(line)` and compute residues on the `(2,2,4)` fibres). Outcome either rules out *all* net-pulled-back denominators (likely, if `α` has a nonzero fibrewise-varying `(2,0)` component — this is what `h^{2,0}(S_r)^+ = 9 ≠ 0` makes possible) or shows `α` is supported on `φ^{-1}(Δ_net)`, giving a specific denominator of degree `2·deg Δ_net` (probably ≥ 20) to test. Chance of a witness: ≈5%; value as a structural result: high.

**E6. Fibration over `P^1` by `τ`-invariant pencils through the arrangement.** Pencils `{λQ_0 + μQ_1}` restricted to `X` give threefold fibres `T_{[λ:μ]} = X ∩ V(λQ_0+μQ_1)` containing all eight planes (hence singular along the 12 lines at least where two contained planes meet). `(3,1)` classes of `X` map to `H^1(P^1, R^3)` and `H^0(P^1, R^4)`; if `H^{2,0}` of a resolution of a *special* fibre is nonzero, that fibre is a candidate support: `D = λ_0Q_0 + μ_0Q_1` of degree 2 — but all degree-2 denominators are excluded by E1. So the only use is as a source of *higher* denominators: products over the finitely many special fibres `Π_s (λ_sQ_0 + μ_sQ_1)`, even, degree `2·#special`. Test by computing the discriminant of the pencil (where `T` acquires extra singularities) — a one-variable resultant computation, cheap — then run E1 on the product. Chance ≈5%.

**E7. Quintic-surface projections from a `τ`-pair of planes.** Projection from `Π_ε` makes `X` birational to a family of quintic surfaces over `P^2`; composing with `τ` gives the conjugate family from `Π_{−ε}`. The relevant Hodge bundle has rank `p_g = 4` per fibre. Candidate denominators: the pull-back of the discriminant of the quintic family (where `S_a` becomes singular), and the "cone over the plane" type denominators `D = Π_ε`-ideal elements: products `Π_{a}(linear forms vanishing on Π_ε)·(conjugates)` (degree 8, even if chosen `τ`-symmetrically) — **excluded by E1** (degree 8 even). Hence only the discriminant (high degree) remains. Chance ≈5%, and the degree is likely beyond current linear algebra; mostly a theoretical route.

**E8. Nonlinear continuation from Fermat.** Set up the bilinear system `Φ(D,N,Q;t) = 0` on `F_t`, start at `t = 0` from the known block-cube witness (`D_0 = C + Y_1^6 + c·Y_2^6`), and track a solution branch in `(D,N,Q)` with `D` allowed in `S_10^+` (1512 params) and `(N,Q)` of degree 11/14, by Newton/predictor-corrector mod `p` or numerically. The solution variety at `t = 0` has positive dimension (Koszul freedom, scaling, alternative block constructions), so one needs to pick the branch; at `t = 5` the degree-6 branch must either escape to infinity or land in `J_6` (§2.3). Chance of success: low (≈5%): overdetermined bilinear systems typically have solution sets that do not deform along a generic path; but a failure mode is informative (where along `t` does the degree-6 family die?). Cost: large (each Newton step is an `≈ 18 900 × 84 903` solve).

**E9. Search guided by the Gorenstein duality description of the degree-10 kernel.** `Ann_R(A)_{10}^{+} ≅ (R_{14}^{-} / A·R_8^{-})^∨`: elements of the E1 kernel are dual to odd degree-14 classes not divisible by `A`. Choose `D` as the "dual" of geometrically meaningful degree-14 forms (e.g. `Q_0Q_1Q_2 · (products of the partials)`), then fixed-`D` solve. Purely heuristic; chance ≈3%; cheap once E1 at `d = 10` is done.

**E10. Nonconstructive correspondence (the only route that avoids the high-degree wall).** Target: a smooth projective variety `Z` with a correspondence `Γ ⊂ Z × X` such that `Γ_*` is surjective onto the `τ`-invariant primitive `H^4`-sector containing `α` (or at least hits `α`) and `H^4(Z)` is known to be divisor-supported (Fano/uniruled `Z`, `p_g = q = 0` surface fibration, or abelian-type `Z` with known Hodge conjecture). §4.7 closes the obvious choices (`X/τ`, the `S_3` Galois closure of the cubic cover, discriminant double cover). New candidates specific to this `F`: (a) the eightfold cover of `X` obtained by adjoining `√(X_i + iY_i)`-type roots that make the planes Cartier/split, (b) the blow-up of `X` along the eight planes mapping to `P^2` with the `(2,2,4)` fibres — then apply the *generalized* surface-fibration argument to the 9-dimensional invariant Hodge bundle: if the family of `(2,2,4)` surfaces `S_r` is a family of surfaces whose transcendental lattice is governed by a K3-type or abelian-type structure (e.g. if `S_r` is a quotient or cover of a K3 surface — `(2,2,4)` complete intersections with `K = O(2)` are double covers of `(2,2,2)` K3 surfaces branched along a `(2,2,2,4)` curve when the quartic is a square mod the quadrics; here it is not, generically), the Kuga–Satake/K3 machinery of ch. 3a could be invoked. Chance that a workable `Z` exists and is found: ≈10% over a long effort; chance of a formal Lean proof following it: lower (the Lean libraries support polynomial-curve integration and elliptic contractions, not general correspondences).

Summary ranking by (information gained)/(cost): E1 ≫ E2 > E4 > E5 > E3 > E6 > E9 > E7 > E8 > E10 (E10 is the only one with a plausible path to a *proof* if explicit witnesses do not exist at tractable degree).

---

## 6. Feasibility assessment

**Explicit witness.** Rigorous facts: `deg D ≥ 9`, or `D ∈ J_6 ∪ J_7 ∪ J_8^{even}`; six degree-5 denominators inconsistent. Generic expectation: the first non-Jacobian E1-passing denominators appear at `d = 10` in a `≥ 294`-dimensional family; E2 then is a bilinear problem over `P^{293}` × linear space, infeasible by brute force but feasible on structured subfamilies (E4). The accepted witnesses on dense-ish instances all have `d ≤ 6` and were *designed* (contact, block, norm denominators), never found by search; the only dense even sextic in the catalogue with a certificate found by linear algebra alone is none. The book's dimension argument (§4.4: any support must have `Σ h^{2,0}(E_i) ≥ 208` to serve the whole invariant sector, though a single target needs less) and §3.1 (support components must be non-normal or have non-rational singularities) suggest that for a dense `F` the needed divisor is a *highly singular* threefold, typically of high degree. My estimate: probability that a certificate with `d ≤ 10` exists ≈ 15%; with `d ≤ 14` ≈ 30%; existence at *some* degree is what the Hodge conjecture predicts (and the challenge's §4.4 shows one generic certificate would settle all). If a witness at `d = 10` were found: `N_i` have `≈13 104` rational coefficients, `Q` `≈5 796`; the Witness.lean would be an order of magnitude larger than the largest accepted one (`contact46006`, 243 KB at `d = 6`), and kernel-checking a degree-26 polynomial identity with `~10^5` terms is the main Lean risk — the existing pipeline (`leangen.py` → `Solution.lean` template → sandbox replay, verified on an accepted witness in ≈5 min) would need scaling tests.

**Nonconstructive argument.** Must prove `HodgeComputer.HasCertificate` over `C` without exhibiting `D`. Available formal infrastructure: `CertificateOperations` (sums of certificates), `HoldoutCover.GeneralizedCover`, `PolynomialCurves` (lifting a square incidence system to a finite extension, integrating a dominant polynomial-curve family with polynomial normal), `EllipticCompanions`, `HigherResidues`. All of these were built for sparse even sextics with a reflection or a cover on which a reflection appears; none handles a dense `τ`-even sextic. A proof here would need (i) a new geometric support theorem for the `τ`-invariant `(3,1)` sector of this `F` (route E10) **and** (ii) its translation into a certificate (ch. 1 equivalence "support ⇔ certificate" is proved in the book over `C`, but the Lean contract needs the exact polynomial identity or an existential derivation inside Lean; the book's §3.6 "from a computed support to every `q`" is prose, not Lean). Realistically, a nonconstructive Lean proof is out of reach without a major new library, and the mathematical theorem itself (support of the `τ`-invariant `(3,1)` sector of a dense even sextic) is exactly the open frontier that §4.7 describes.

**Bottom line.** This instance is a genuine frontier case: the only symmetry is `τ`, whose sector the book proves is transcendental for every even sextic; the eight planes are rigid and of the wrong Hodge type; both natural fibrations have holomorphic two-forms in the target's character; all denominators of degree `≤ 8` (even) / `≤ 7` (odd) outside `J(F)` are rigorously excluded. The cheapest decisive next steps are E1 (finish the kernel ladder), E2/E3 (sweep Jacobian denominators with the Koszul-reduced E2 solver) and E5 (decide whether the net can support the class at all). A witness, if one exists at tractable degree, will most likely come from E4-type structured degree-10 families containing the arrangement.

---

## References and links

* Discussion #176: https://github.com/hodge-computer/sextics/discussions/176 (comment `DC_kwDOVLOOJc4BJDZL`, `benediamond`, 2026-10-23).
* Instance files: `data/eight-plane-configuration/{instance.rec, why.md, provenance.json, checks.json}`; audit: `scripts/audit-structural-frontier.py`.
* Book: `docs/02-explicit-constructions.md` §2.2 (Nambu contact, adjunction interpretation), §2.3 (perfect-power contact); `docs/02a-binary-blocks.md` §2a.1; `docs/03-geometric-support.md` §3.1 (support ⇒ non-normal/non-rational component), §3.2 (surface-fibration theorem, covering curves), §3.3 (reflection theorem, reflection groups, Kollár–Larsen age criterion, residue character formula `Σ(e_j+1)w_j`); `docs/04-known-locus-and-frontier.md` §4.2 (paired-sign sector theorem), §4.4 (one generic target suffices; `(218,208,866)`), §4.7 (`X/τ` not uniruled; cubic cover with `S_3` monodromy and `V_std`).
* Kollár–Larsen, https://arxiv.org/pdf/math/0701466 (Theorem 2, Corollary 4). Diamond, Fermat case, https://arxiv.org/abs/2603.29157. Schütt–Shioda https://arxiv.org/abs/0907.0298 (cited in §4.7).
* Local computation artifacts: `/tmp/hodge/results/eight-plane-configuration/{log.jsonl, e1kernel.jsonl, e1kernel_basis_d*_par*.txt, K_decomp.json, lib.py, e1kernel.py, run_J6.py, koszul_e2.py}`; structural verification script `/tmp/hodge/scratch/epc_struct.py`; engine `~/devel/sextics-work/hodge/`.

## Addendum: results of computations that were in flight (as of 2026-11-01 ≈18:25 local)

All entries below are engine-2 modular computations (`p = 2^31 − 1`), logged in `/tmp/hodge/results/eight-plane-configuration/log.jsonl`; inconsistency mod `p` of a system with `p`-integral coefficients implies inconsistency over `Q` (an exact rational solution would reduce mod `p`), so each "inconsistent" line is a rigorous exclusion of that denominator.

**A. Jacobian denominators of degree 6 (`run_J6.py`, in progress, one result so far).**

| `D` | parity | E1 | fixed-`D` system | result |
|---|:-:|:-:|---|---|
| `E_0(F) = X0·F_{X0} + Y0·F_{Y0}` (block-0 Euler derivative) | even | solvable (automatic) | `16 797 × 2 602`, rank 2359, nullity 243, 412 s | **inconsistent** |
| `E_1(F) = X1·F_{X1} + Y1·F_{Y1}` (block-1 Euler derivative) | even | solvable (automatic) | `16 797 × 2 602`, rank 2359, nullity 243, 631 s | **inconsistent** |

(The identical rank/nullity for the two block Euler derivatives suggests the nullity 243 is structural: 18-dimensional Koszul freedom `Σ c_{ij}(F_j∂_i − F_i∂_j)` with `c_{ij}` quadrics of the right parity, plus the `Q`-freedom from `F·S_6^{odd}`-type identities; it would be worth extracting the kernel once to understand it.)

Still queued in that run: `E_2(F)`, `Σ_i X_iF_{X_i}`, `X0F_{X0}`, `Y0F_{Y0}`, `X1F_{X0} − X0F_{X1}`, `R_0(F) = Y0F_{X0} − X0F_{Y0}` (odd, zero locus ⊃ eight planes), `X0F_{Y0}`, `Y0F_{X0}` (≈7 min each). Results will be appended to `log.jsonl` under `family ∈ {block-radial, block-euler, J6-monomial, J6-rotation, J6-rotation-odd, J6-monomial-odd}`.

**B. E1 sanity sweep of geometric denominators (attack agent + this dossier agent).** All of the following fail E1 (`rank_aug = rank + 1`), exactly as the kernel theorem of §3.3 predicts since none lies in `J(F)`: degree 2: `Q0, Q1, Q2`; degree 4: `Q_i^2`, `K_i`, `H_i`, `diag_i := X_i^4 − X_i^2Y_i^2 + Y_i^4` (`i = 0,1,2`); degree 6: `Q_i^3`, `Q_i·K_i`, `F − Q_iK_i` (≡ `−Q_iK_i` mod `F`), `X_i^6 + Y_i^6 = Q_i·diag_i`, `Q_i·diag_i`, and `Q0·Q1·Q2` (the reduced equation of the whole arrangement). In particular **`Q0Q1Q2 ∉ J(F)_6`** and the arrangement's own equation is excluded as a denominator. Also failed: `Q0^2·Q1`, `Q2·K2`, `F − Q2·K2`. (With `Q_i^3`, `Q_i^2Q_j`, `Q0Q1Q2` all excluded, every degree-6 monomial in the `Q_i` — hence, by the kernel theorem, every even sextic form in the ideal `(Q_0,Q_1,Q_2)` outside `J(F)` — is excluded.)

**C. E1-kernel at `d = 9`** (`e1kernel.py 9`, both parities, `7 752 × 10 010` per parity): still running when this dossier was closed (started 18:07, ≈30 min elapsed; dense nullspace of this size may take an hour); output appends to `e1kernel.jsonl` (fields `kerdim`, `dimJ_d = 378`, `kermodJ`) and `e1kernel_run_d9.out`. Expected outcome `kermodJ = 0` (square map `R_9 → R_15`, generically nondegenerate; the attack agent notes it is `623 × 623` per parity); a nonzero value would be a structural surprise and would immediately provide degree-9 denominators to test.

**D. Consequence for the status line.** Combined with §3.3: every certificate for this instance has either `deg D ≥ 9`, or `D ∈ J(F)_d` with `d ∈ {6,7,8}` (`d = 8` only even) and `D` not among the eight Jacobian denominators tested so far (`∂F/∂Z_i`, `E_0(F)`, `E_1(F)`).
