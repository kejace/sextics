# Research dossier: `reflection-coupled-residual` (OPEN)

Prepared 2026-11-01 from the repository snapshot at `/home/kejace/devel/sextics` (commit `83fbf5c`),
the saved discussion thread `/tmp/hodge/discussions/reflection-coupled-residual.md` (GitHub discussion
#154: 2 top-level comments, 2 replies; last activity 2026-10-28), the companion discussion #111 (read
via the GitHub API, read-only), the accepted companion solution, the book chapters cited by `why.md`
(02 §2.7–2.8, 03 §3.1–3.3, 04 §4.4–4.7), the sibling-agent scratch results in
`/tmp/hodge/results/reflection-coupled-residual/` (local, unpublished), `DiamondJM`'s public sweep
repository (`/tmp/hodge/ext/hodgeHoldout`), and a few independent computations of my own (plain
Python/NumPy/SymPy; scripts and raw outputs in `/tmp/hodge/scratch/rcr-dossier/`).

Conventions. Variables are ordered `X0, X1, X2, Y0, Y1, Y2`; `S = C[Z]`; `J = (∂_0F,…,∂_5F)`;
for a six-tuple `N`, `N(P) = Σ N_i ∂_iP`, `div N = Σ ∂_iN_i`, `B_D(N) = D·div N − N(D)`.
A certificate for `(F,q,A)` is `(D,N,Q)` with `D ≠ 0` homogeneous, `gcd(D,F) = 1` over `C`,
`deg N_i = deg D + 6q − 5`, and

    E := A·D² − q·D·N(F) + F·B_D(N) = F^{q+1}·Q .

Degrees: `deg E = 2·deg D + 6q`, `deg F^{q+1} = 6q + 6`, hence **`deg Q = 2·deg D − 6`, independent
of `q`**. At `q = 1` (this instance) and `d = deg D`: `N_i ∈ S_{d+1}`, `Q ∈ S_{2d−6}`, `E ∈ S_{2d+6}`.
Euler normalization (Chapter 1) makes this equivalent, for fixed `D`, to the two linear conditions
`E1: q·N(F) = A·D` and `E2: B_D(N) ∈ (F^q)`. No denominator-degree bound and no `Q = 0` is imposed.

Labels: **[proved]** (proof given or cited and checked by me), **[computed]** (exact or modular
computation by me, this session), **[reported]** (claimed in the thread / accepted solutions /
sibling scratch logs; not replayed by me), **[inference]** (my own reasoning, unchecked by others),
**[conjecture]**, **[failed]**.

---

## 0. Status at a glance

* Target: `q = 1`, `A = X0³Y0³` on the coupled sextic `F = Σ Z_i⁶ + (5/7)·P` (13 monomials). No
  certificate, no existence proof, no Lean theorem exists. `existenceStatus`: "No existence conclusion
  for this target follows from the screened quotient-density constructions."
* Same variety as the **solved** `reflection-coupled-fano` (`fqlx`, PR #150, accepted 2026-10-20,
  existential, no degree exported), whose target `X0³Y2³` is odd in `Y2`. The present target is
  **even in `Y2`**: its residue has the opposite character under the sole reflection `r: Y2 ↦ −Y2`,
  so the Fano-quotient construction cannot reach it "by changing only the quotient source".
* Symmetry: the full monomial/diagonal projective symmetry group is the **Klein four-group**
  `{1, r, ρ, τ}` with `r: Y2 ↦ −Y2`, `ρ: (Y0,Y1) ↦ −(Y0,Y1)`, `τ = rρ: Y ↦ −Y` **[computed]**.
  The residue of `A` is `τ`-invariant, `r`-anti-invariant, `ρ`-anti-invariant. Its sector has
  Hodge numbers `(h^{3,1}, h^{2,2}_prim, h^{1,3}) = (126, 492, 126)` **[computed; agrees with the
  book's "intermediate one-reflection locus", §4.4]**; the companion's sector (served by the Fano
  quotient `W = X/r`) has `(82, 374, 82)`. The challenge's `Y`-odd target space is exactly
  `82 ⊕ 126`; the 126-dimensional part is entirely open on this variety.
* Rigorous exclusions so far (homogeneous problem, all numerators, unrestricted `Q`, over `C`):
  **every denominator of degree ≤ 4** (the `E1`-admissible space is zero in every character)
  **[computed, this session]**; for degrees 5–7 the `E1`-admissible space equals the Jacobian ideal
  `J_d` (fully for the `Y`-even sector **[reported, sibling scratch]**, by character up to the degrees
  recorded in §4 **[computed]**), so a certificate of degree ≤ 7 must have `D ∈ J(F)`, a small
  bilinear problem never closed. All 61 "natural degree-6 supports" of `DiamondJM`'s sweep and 42
  degree-4 cone/binary denominators of the sibling sweep fail **[reported]**.
* The discussion (4 posts, `fqlx` and `benediamond`) builds an elaborate genus-one / incidence-curve
  reduction of the actual source on `fqlx`'s quartic family; it ends in an unsolved rational-integration
  problem on an algebraic curve with a **proved pole bound (65·H_∞) for the zero-residue baseline** but
  **no uniform bound for the three residue functions**, and a direct F4SAT saturation that blew up.
  Nothing in the thread is a certificate; nothing rigorously excludes any homogeneous denominator.
* One concrete new object **[reported, sibling scratch, mod p]**: at degree 8 there is a **unique
  (mod `J_8`) `Y`-even denominator class `D8` with `A·D8 ∈ J(F)`**; its full `(N,Q)` system has not
  been solved yet (its rational form was reconstructed by the sibling agent while this dossier was written; the full `(N,Q)` test is the next step).

---

## 1. The exact target

### 1.1 Data

From `data/reflection-coupled-residual/instance.rec` (SHA256 of the record, per `fqlx`:
`114c5b4f2c4c4ba36456a88fc19925d28e51399b4957f417ea5ca157904d38a2`):

    q = 1
    A = X0**3*Y0**3
    F = X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6
      + (5/7)*( X0**3*X1*Y0*Y1 + X0*X1*Y0**3*Y1 + X0*X1*Y0**2*Y2**2 + 2*X0*X1*Y1**2*Y2**2
              + 3*X0*X2**3*Y0*Y1 + 2*X1**3*X2*Y0*Y1 + X1*X2*Y0*Y1**3 )

Monomial support (13 monomials, exponent vectors in the order `X0 X1 X2 Y0 Y1 Y2`) **[computed]**:

    (6,0,0,0,0,0) (0,6,0,0,0,0) (0,0,6,0,0,0) (0,0,0,6,0,0) (0,0,0,0,6,0) (0,0,0,0,0,6)
    (3,1,0,1,1,0) (1,1,0,3,1,0) (1,1,0,2,0,2) (1,1,0,0,2,2) (1,0,3,1,1,0) (0,3,1,1,1,0) (0,1,1,1,3,0)

Sparse presentation. With `c = 5/7`, `p = Y0Y1`, `e = X0X1`, `u = Y2²`:

    F = Σ X_i⁶ + Σ Y_i⁶ + c·[ p·C(X) + e·G(Y) + X1X2·Y0Y1³ ],
    C(X) = X0³X1 + 3X0X2³ + 2X1³X2,            G(Y) = Y0³Y1 + Y0²Y2² + 2Y1²Y2².

Equivalently (the companion's notation, `x = (X0,X1,X2,Y0,Y1)`, `u = Y2²`):

    F(x,Y2) = H_c(x, Y2²),   H_c = x0⁶+x1⁶+x2⁶+x3⁶+x4⁶+u³ + c·P(x,u),
    P = x0³x1x3x4 + x0x1x3³x4 + x0x1x3²u + 2x0x1x4²u + 3x0x2³x3x4 + 2x1³x2x3x4 + x1x2x3x4³.

The bihomogeneous `(X,Y)`-decomposition `F = F_{6,0} + F_{4,2} + F_{2,4} + F_{0,6}` is

    F_{6,0} = ΣX⁶,   F_{4,2} = c·Y0Y1·C(X),   F_{2,4} = c·[X0X1·G(Y) + X1X2·Y0Y1³],   F_{0,6} = ΣY⁶.

Note that `F_{4,2}` is a *product* `c·p·C(X)`; this is relevant to the cubic-cover presentation of
§4.7 (see §5, E9).

### 1.2 Degrees and parities

* `deg F = 6`, `deg A = 6`, `q = 1`. For `deg D = d`: `N_i ∈ S_{d+1}`, `Q ∈ S_{2d−6}`, `E ∈ S_{2d+6}`.
* `F` is even in `Y2` alone (`Y2`-degrees of the monomials: `0,0,0,0,0,6,0,0,2,2,0,0,0`) and even
  in `Y` jointly (`Y`-degrees: `0,0,0,6,6,6,2,4,4,4,2,2,4`); hence also even under `ρ`.
  `F` is **not** even in any `X_i`, nor in `Y0` or `Y1` individually.
* `A = X0³Y0³` is `Y`-odd (`τ`-odd as a polynomial), `Y2`-even, `ρ`-odd.
* Linear-system sizes at `q = 1` (used in §5):

  | `d` | `dim S_d` | unknowns `N` (`6·S_{d+1}`) | unknowns `Q` (`S_{2d−6}`) | equations (`S_{2d+6}`) | `E1` rows (`S_{d+6}`) |
  |---|---|---|---|---|---|
  | 5 | 252 | 2,772 | 126 | 20,349 | 4,368 |
  | 6 | 462 | 4,752 | 462 | 33,649 | 6,188 |
  | 7 | 792 | 7,722 | 1,287 | 53,130 | 8,568 |
  | 8 | 1,287 | 12,012 | 3,003 | 80,730 | 11,628 |
  | 9 | 2,002 | 18,018 | 6,188 | 118,755 | 15,504 |
  | 10 | 3,003 | 26,208 | 11,628 | 169,911 | 20,349 |
  | 11 | 4,368 | 37,128 | 20,349 | 237,336 | 26,334 |
  | 12 | 6,188 | 51,408 | 33,649 | 324,632 | 33,649 |
  | 13 | 8,568 | 69,768 | 53,130 | 435,897 | 42,504 |
  | 14 | 11,628 | 93,024 | 80,730 | 575,757 | 53,130 |
  | 15 | 15,504 | 122,094 | 118,755 | 749,398 | 65,780 |

### 1.3 Symmetry **[computed]**

* Exponent-difference lattice of the 13 monomials has rank 5 and Smith normal form `(1,1,1,2,2,0)`.
  Hence the projective diagonal symmetry group is `Z/2 × Z/2` (order 4; agrees with
  `checks.json: projectiveDiagonalGroupOrder = 4`). An exhaustive `μ12` search confirms exactly four
  weight vectors (with `X0`-weight normalized to 0): `0`, `(0,0,0,0,0,6)`, `(0,0,0,6,6,0)`,
  `(0,0,0,6,6,6)` (mod 12), i.e.

      r = diag(1,1,1,1,1,−1),   ρ = diag(1,1,1,−1,−1,1),   τ = rρ = diag(1,1,1,−1,−1,−1).

  There is **no `μ3` or `μ6` diagonal action** (the holdout-style weights `(0,0,0,5,1,5)` are broken
  by `X0X1Y0²Y2²`, weight `2 mod 6`, etc.). Consequently the §2.8 "complementary diagonal quotient"
  and the Kollár–Larsen age argument for order-3/6 elements are unavailable.
* Residue characters (book §3.3: a monomial `Z^e` has residue character `Σ(e_j+1)w_j`; for sign
  actions, numerator sign times determinant):

  | element | `F` | det | `A` as polynomial | residue of `A` |
  |---|---|---|---|---|
  | `r: Y2 ↦ −Y2` | `+F` | `−1` | `+1` | **`−1`** |
  | `ρ: (Y0,Y1) ↦ −(Y0,Y1)` | `+F` | `+1` | `−1` | **`−1`** |
  | `τ: Y ↦ −Y` | `+F` | `−1` | `−1` | **`+1`** |

  So the residue of `A` lies in the `τ`-invariant sector (as every challenge target must: the
  holomorphic 4-form `Ω/F` is `τ`-anti-invariant, so `H⁴(X)^τ ⊂ F¹H⁴` and the general Hodge
  conjecture predicts divisor support for it **[inference]**), but it is anti-invariant under both
  `r` and `ρ`. The companion target `X0³Y2³` has residue characters `(r,ρ,τ) = (+1,+1,+1)`.
* Permutation symmetries: the identity is the only permutation of the six variables preserving the
  monomial support; hence the full monomial (signed/rooted permutation) automorphism group modulo
  scalars is exactly the Klein group above. `F(Y,X) ≠ F` (no `X ↔ Y` exchange). Non-monomial linear
  automorphisms have not been excluded (see §5, E5), but `checks.json: hessianSymmetrizerDimension = 1`
  shows the rational Hessian symmetrizer space is scalar, so `F` is not a linear direct sum of blocks
  in any coordinates (`why.md` of the companion: "ternary-split and binary-block formulas do not apply
  after a linear change of coordinates").
* Sector Hodge numbers from the equivariant Hilbert series of the Jacobian ring
  `Π_i (1 − t⁵χ_i)/(1 − tχ_i)` (the partials form a character-homogeneous regular sequence) **[computed]**;
  "numerator character" `(a,b)` = (`Y2`-parity, `Y0+Y1`-parity); residue character = numerator
  character times `det`, i.e. `a ↦ a+1`:

  | numerator char `(a,b)` | residue char | `(S/J)_6` | `(S/J)_12` | `(S/J)_18` | served by |
  |---|---|---|---|---|---|
  | `(1,0)` (`Y2`-odd, `Y`-odd) | `(0,0)` invariant | **82** | 374 | 82 | Fano quotient `W = X/r` (companion) |
  | `(0,1)` (`Y2`-even, `Y`-odd) — **this target** | `(1,1)` | **126** | 492 | 126 | nothing known |
  | `(0,0)` (`Y`-even) | `(1,0)` | 144 | 529 | 144 | not a challenge target (`τ`-odd) |
  | `(1,1)` | `(0,1)` | 74 | 356 | 74 | not a challenge target (`τ`-odd) |

  Totals `426 / 1751 / 426` are the primitive Hodge numbers of a smooth sextic fourfold. The
  `τ`-invariant total `82 + 126 = 208` and the `τ`-anti-invariant total `144 + 74 = 218` match the
  book's `(dim V, dim G) = (208, 218)`; the `(126, 492, 126)` row matches the "intermediate
  one-reflection locus" of §4.4 verbatim. The 126-dimensional `(3,1)` sector of this target is spanned
  by the 138 `Y`-odd, `Y2`-even sextic monomials modulo the 12 Jacobian elements `Z_i∂_jF` with
  `i ∈ {Y0,Y1}, j ∈ X` or vice versa.
* Non-diagonal symmetries noted by anyone: none. `benediamond` (C2) uses "the two block-sign
  symmetries" (`r`, `ρ`) to halve his current system (twelve coefficient ratios, eleven equations).

### 1.4 Block and split structure

* `F = H_c(x, Y2²)`: the quotient by `r` is the smooth sextic `W = {H_c = 0} ⊂ P(2,1,1,1,1,1)`
  (weights: `u` has weight 2), Fano with `−K_W = O_W(1)` (book §3.3). The branch divisor is
  `B_W = W ∩ {u = 0}`, i.e. the smooth sextic threefold `{H_c(x,0) = 0} ⊂ P⁴`, of class `O_W(2)`.
* `F` is also even under `ρ`, i.e. `F ∈ C[X, Y0Y1, Y0², Y1², Y2²]`; the quotient `X/ρ` has
  `½(1,1)`-singularities along the surface `X ∩ {Y0 = Y1 = 0}` (a smooth sextic surface in `P³`,
  tangent eigenvalues `(−1,−1)`, age 1) and is not of the Fano/branch-divisor type. The quotient
  `X/τ` has the two plane sextic curves `X ∩ {X=0}`, `X ∩ {Y=0}` as fixed locus (age `3/2`) and is not
  uniruled (book §4.7, Kollár–Larsen). **The target's residue survives only in `X/τ`** among the
  Klein quotients, and that quotient has no known support mechanism.
* No three disjoint sign pairs (only two independent signs), so the paired-sign cubic / Jacobian
  route of `cycle-three` (§2.9) does not apply. No `T⁶ = Y_i⁶ + Y_j⁶` exchange structure (§2.7): with
  `T⁶ = Y0⁶ + Y1⁶`, `Y0 = αT`, `Y1 = βT`, the `Y`-part becomes
  `T⁶ + Y2⁶ + c[T²αβ·C(X) + X0X1(T⁴α³β + T²Y2²(α²+2β²)) + X1X2T⁴αβ³]`, which is not symmetric in
  `(T, Y2)` (the `T⁴` terms have no `Y2⁴` partner) **[computed by hand]**; the other coordinate
  pairs are worse because `Y2` appears only through `Y2²·(Y0² + 2Y1²)X0X1`.

### 1.5 Smoothness **[proved, cited]**

`checks.json`: `smoothReductionPrime = 5`. Since `7` is invertible mod 5, `c = 5/7 ≡ 0`, so
`F ≡ Σ Z_i⁶ (mod 5)` with reduced partials `Z_i⁵`, whose only common zero is the origin; good reduction
at 5 implies smoothness over `Q`, hence over `C`. `targetNonzeroInReducedJacobianRing = true`
(`X0³Y0³ ∉ (Z_i⁵)`), so the target class is nonzero. Irreducibility of `F` over `C` is proved in the
companion's `CoupledDomain.lean` by an Eisenstein specialization (`x1 = x2 = x4 = 0`, `x3 = T+1`,
`z = i`).

---

## 2. Companion and related solved instances

| instance | status / PR | method | denominator | why it does not transfer |
|---|---|---|---|---|
| **`reflection-coupled-fano`** (same `F`, `A = X0³Y2³`) | solved, `fqlx` PR #150, accepted 2026-10-20 (owner review by `benediamond`); discussion #111 | Reflection quotient `u = Y2²`; modular weighted quartic on `H_0 = Σx_i⁶ + u³` over `F_25` (`j² = 2`): `x0 = j(s⁴+t⁴), x1 = (2+j)(s³t+st³), x2 = js²t², x3 = j(s³t−st³), x4 = s⁴−t⁴, u = s⁸+s⁴t⁴+t⁸`; invertible `24×24` incidence Jacobian (coefficients of `t¹..t⁴` in six coordinates); étale lifting over the integral endpoint ring `Z[1/7][r2,r3,r4,r5,ξ]/(ξ⁶ + 3cξr2³r3r4 + r2⁶+r3⁶+r4⁶+r5³)` with `c = 5/7` fixed; polynomial normal `ψ` with `Σ G_iψ_i = 1`; five commuting derivations; evaluation determinant `g⁻¹v_{1,1} ≠ 0`; polynomial integration in `t`; normalized trace to the quotient function field; legal clearing `D_H = (∂_0H_c)·d`; reflection transport `D_F = φ(D_H)`, `N_{F,i} = zφ(N_{H,i})` (`i<5`), `N_{F,5} = ½φ(N_{H,5})`; homogeneous extraction | **existential**; no `witness.json`, no degree exported ("without evaluating their full coefficient lists") | Transport needs the source identity `A = z·φ(B)` with `B` polynomial on `W`, i.e. `A` **odd in `z = Y2`**. `X0³Y0³` is even: "the map `u = Y2²` has Jacobian density `2Y2`; changing the polynomial source on the quotient changes the even factor multiplying that density, but cannot produce this target" (`why.md`). Structurally: the quotient `W` serves exactly the `r`-invariant residue sector (82 of the 208 `Y`-odd `(3,1)` directions); "neither an additional incidence seed on the same quotient nor a higher-order integration identity enlarges that quotient's averaging ideal" (§4.7). |
| `holdout` (`fqlx` PR #128), `holdout-ternary-deformation` (`hodgepodgeconjecture` PR #151), `one-reflection-y2-sector` (`DiamondJM` PR #144) | solved | Fermat cover `T⁶ = Y0⁶ + Y2⁶`, exchange `y ↔ T`, quotient in `P(1,1,1,2,2)`, moving-conic primitive, deck character kills the correction `T^e` (`e = b − a − c − 1 ≢ 0 mod 6`) | existential | requires the chain-type `Y`-structure (`Y0⁶+Y1⁶+Y2⁶` plus couplings linear in `Y0Y1`, `Y1Y2`); here no coordinate pair admits an exchange (§1.4). The machinery (`HoldoutCover`, `PolynomialCurves`) is reusable once a serving geometry is found. |
| `holdout-bundle` (`fqlx` PR #161; §2.8) | solved for `X0³Y0³ + X0³Y0Y1²` on the chain sextic | complementary diagonal quotient by weights `(5,1,5)` on `Y`; invariant coordinates `p = Y0Y1, v = Y1Y2, λ = Y1³/(Y0³+iY2³)`; polynomial fiber source `ρ = −X0³p/(6λ)`; conic seed over `F_5` with `det Δ = 3`; `HasCertificate.add` | existential | there is no order-3 or order-6 diagonal symmetry of this `F` (§1.3); `fqlx` (C1) cites #161 only for its "source-transport machinery". |
| `cycle-three` (`hodgepodgeconjecture` PR #165), `mixed-biweight` (`fqlx` PR #170) | solved, same `A = X0³Y0³` | three paired signs ⇒ plane cubic in `X_i²`; Hessian covariants ⇒ Weierstrass Jacobian `η² = ξ³ + s_4ξ − s_6/4`; `ResidueTransport`, moving-map correction | existential | no three disjoint sign pairs here. `fqlx` (R1) notes "[PR165] is useful source-transport machinery, but matching its polynomial-source hypothesis is still open here." |
| `one-reflection-skew` (OPEN; same `A`, `F` with the one `Y2` reflection and the same residual character) | open; see `/tmp/hodge/notes/dossier-one-reflection-skew.md` | — | — | closest open sibling: identical sector type (τ-invariant, `r`-anti-invariant), degrees ≤ 6 and most of 7 excluded there; its CM-elliptic-companion route is specific to its skew term. |

Book boundaries relevant here: §3.3 ("A target even in `Y2` has the opposite residue character and does
not follow by changing only the quotient source"); §4.4 (the residual one-reflection sector, dimensions
`(126,492,126)`: "one dominating support carrying a nonzero residual (3,1) class carries this entire
residual sector" — a statement about a *dominating family* on the 158-dimensional coefficient space,
not about an isolated instance; "a support with the required nonzero residual image has not been
constructed"); §4.6 (deck/exchange descent criterion: needs a cover `Z → X` with commuting involutions
`r_Z, σ`, a deck `d` with `pd = p`, `pe = gp` for `e = σdσ`, `pr_Z = τp`, and supported
`Z/⟨r_Z,σ⟩` and `X/⟨g⟩`; with `g = r` the second input is the Fano quotient `W`, so **the missing
input is a double cover `Z` with an exchange involution `σ`**); §4.7 (`X/τ` not uniruled; the natural
`S_3` cubic cover's standard representation cannot be served through uniruled intermediate quotients).

---

## 3. The discussion thread (GitHub discussion #154), comment by comment

All contributors disclose AI assistance (ChatGPT GPT-6 Pro / Codex / "GPT-6 Astra Pro"). Nobody
claims a certificate, a Lean proof, or a nonexistence result. All "exclusions" below are of
*particular currents in particular presentations*, never of homogeneous denominators.

### C1. `fqlx`, 2026-10-20 23:52 — "dominant genus-one construction and explicit partial primitive" [reported]

URL: https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19084541 (id `DC_kwDOVLOOJc4BIzT9`).

* Setting: `H(x,u) = F(x,z)`, `u = z²`, exact coupling `5/7`. A *divided conic* incidence (degree-2
  map `φ(t)` to the quotient, with `U(t)` the weight-2 coordinate) imposing all 13 coefficients of
  `H(φ(t),U(t)) = 0` "has a nonsingular algebraic characteristic-zero lift". Note the companion
  solution proved that a *degree-2 modular* construction on the characteristic-5 Fermat reduction has
  twisted normal rank ≤ 10 < 12; `fqlx`'s conic here is a different, characteristic-zero object with a
  further slice: `Y0 = φ_3` linear in `t`, `U` monic quartic
  `U = t⁴ + u_3t³ + u_2t² + u_1t + u_0` with the four lower coefficients as independent parameters.
  Dominance on "the same smooth component" is claimed via exact incidence/evaluation minors and the
  endpoint-preserving `GL_2` action; the evaluation determinant vanishes at the original seed, but the
  nonzero derivative of `[t]Y0` and the gauge determinant `−4[t]Y0` give generic transversality.
* Physical cover: `w² = U(t)`, squarefree, **genus one**. The normal source on the new slice is
  `p_A/(2w)`. Rational currents remove the polynomial-exact and infinity terms (a "normal-transgression
  identity `ρ = Σ_j D_jK_j`", moving infinity to the source-zero locus `Y0 = 0`).
* The remaining unrestricted equation on this construction, with `E_2 = t² + ½u_3t`:

      Σ_{α∈{u_0,u_1,u_2,u_3,t}} ( 2U·D_αC_α − (D_αU)·C_α ) = U·(c_0 + c_1E_2),   C_α ∈ E(t),

  `E` = the incidence coefficient field, finite separable over `C(u_0,…,u_3)`. "All parameter
  directions, arbitrary rational denominators, and moving time poles remain allowed."
* Computed obstruction: the compact coefficients are nonzero, `5c_0 ≡ 127 (mod 625)`,
  `c_1/5 ≡ 16 (mod 625)`, "with a proved precision budget". **This excludes time-only completion of
  these particular currents (also after finite algebraic extension); it excludes nothing else.**
* Packet: https://github.com/user-attachments/files/33743911/reflection_residual_genus_one_continuation.zip
  (SHA256 `eaeeeee98815a5ded398e420f3afda8f81f752f7311ac295f358bfe51099adc7`; `python3 verify_all.py`,
  SymPy 1.14.0). Credits PRs #161, #150, repository `PolynomialCurves`.

### R1. `fqlx`, 2026-10-21 02:49 (updated 02:51) — "the actual source now reduces to one holomorphic coefficient" [reported]

URL: https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19087146 (id `DC_kwDOVLOOJc4BIz8q`).

* Deform `G_ε = G + εA` at fixed `U` and gauges; for the full moving-normal determinant `𝒥`, a
  cofactor-degree argument proves `𝒥` and `𝒥̇` are time-constant (also over dual numbers); explicit
  cofactor currents leave the scalar source `g = −½𝒥̇`.
* Full identity for the `(G_2,G_4)` normal:
  `p_A/(2w) = div_5((C^{(0)}+C^{(1)})/w) + R_*/w`, all five current components explicit; 64
  substantive checks replayed; the nine coefficients of the cleared degree-≤8 identity vanish
  mod `5^{143}`; on the characteristic-zero branch `v_5(R_*) = 0`, `R_* ≡ 529 (mod 625)`, so "this
  selected combination is incomplete".
* "Rolling wheel" finishing condition: in depressed coordinates `U = x⁴ + px² + qx + r`, `a = u_3`,
  `x = t + a/4`, find rational `h, k ∈ E` with `g = T(h) + V(k)`, where
  `T = ∂_a`, `V = 6q∂_p + (8r − 2p²)∂_q − pq∂_r`; then `C = hT + k(V + (p+2x²)∂_x)` gives
  `div_5(C/w) = g/w`. **Not found; "scaling leaves an Euler resonance".**
* Reproduction: https://gist.github.com/fqlx/ae3df1cc3d0a3857b7b69f2ad80fc983 (field/source
  definitions, currents, proof, raw jets, verifier, receipts).

### R2. `benediamond`, 2026-10-22 18:13 — "Radial audit of the source-root currents" [reported]

URL: https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19101866-ish (id `DC_kwDOVLOOJc4BI-ZY`).

* Same incidence field and compact basis `(dt/w, (x²+p/6)dt/w)`. For the radial field `S` define
  `𝒟_S f = S(f) + (div S)f + Γ_S f`, `Γ_S = diag(1,−1)`, `div S = −10 + Tν`; covariant Euler identity:
  `𝒟_S f = mf, m ≠ 0 ⇒ Div_∇(S^j f/m) = f`.
* For `fqlx`'s published `Y0`-root current `θ_Y = v⁻¹ι_{W̃}μ`: `dθ_Y = 0`, `ℒ_Sμ = 0`, `[S,W̃] = 0`,
  `S(v) = −2v`, hence `ℒ_Sθ_Y = 2θ_Y` and `θ_Y = ½d(ι_Sθ_Y)`: **the `Y0`-root Cartan forcing already has
  an explicit compact primitive `C_Y^j = S^jk_Y/2`** (weight 2). Same for a closed Cartan current from
  the two finite `X0` roots. But the fixed-infinity compact source `h` satisfies `𝒟_S h = 0` (weight
  0): "constant combinations of these weight-two forcings cannot equal the nonzero source".
* A weight correction needs an even rational multiplier `f` with `Sf = −2f`, `W̃(f) = 0`; at the
  `Y0` root this is `Dφ = (Dr/2r)φ` (`φ/√r` an odd first integral); the evident multiplier `√R` has
  the wrong parity. Certified nonzero 5-adic minors exclude constant combinations of `B²` and the four
  ordinary root coordinates (a "five-dimensional exclusion"). The equation `g = T(h) + V(k)` remains
  unsolved.
* Packet: https://github.com/user-attachments/files/33835113/reflection-coupled-residual-radial-audit-2026-10-22.zip

### C2. `benediamond`, 2026-10-27 12:25 (updated 2026-10-28 15:41) — "Rational integration on the incidence curve, with independent compact coefficients and closed residue freedom" [reported]

URL: https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19140000-ish (id `DC_kwDOVLOOJc4BJTkW`).
Packet: https://github.com/user-attachments/files/34087387/reflection-residual-research.zip (notes, programs, data, checks).

* The four logarithmic compact currents have independent rational coefficient functions; relations
  leave two effective compact directions; eliminating one gives `du = ω_g`, `g = (g_Y, g_0, g_1)`
  — rational integration **on the incidence curve** (the curve of compatible incidence coefficients over
  the field of centered quartic parameters; distinct from the genus-one cover `w² = U`). The three
  residue functions `g_Y, g_0, g_1` (the `Y0` residue current and two channels of the quadratic `X0`
  block) are unknowns subject to closure/weight equations. Formal local solvability holds at a verified
  regular germ; **global rationality is the open requirement**; one rational solution passing
  reconstruction would suffice.
* **Baseline `g = 0`** [proved pole bound, reported]: with incidence coordinate `z`, Jacobian minor
  `δ`, tangent derivation `𝒯 = δ∂_z`, source `ω_0 = (N/δ³)dz`, every rational primitive satisfies

      f = δu ∈ L(65·H_∞),   δ𝒯(f) − 𝒯(δ)f = N,

  on the smooth projective normalization (`H_∞` = largest pole order of the affine incidence coordinates
  at each point at infinity). "The number 65 follows from the actual source and determinant degrees. It
  is not an imposed denominator-degree cutoff." The pole-bound proof covers every normalized place
  (key step: `ω_0 = (H − L(V))dz`, `V = B/δ`, `H = C/δ`, change the lift of the base derivation in a
  local uniformizer; "One Jacobian factor suffices").
* Structural results: the two block-sign symmetries (`r`, `ρ`) give a quotient presentation with 12
  coefficient ratios and 11 equations; character averaging preserves solvability and keeps the three
  residue functions; `f` is invariant. A quartic degeneration certifies ≥ 16 places at infinity across
  regular rank-3 components with cubic ramification; parity excludes local logarithmic terms there;
  the baseline quotient pole bound for `f` there is 31. An explicit kernel classification excludes the
  source-zero conic covers from the closure of every rank-3 component. A regular component is selected
  by a rigorously isolated algebraic point; 91 certified points in a quotient fiber; a period
  calculation has an independently checked evaluator and 14 certified simple branch points but no
  established nonzero period.
* **Computational bottleneck [failed]:** direct F4SAT saturation on the fiber `(p,q,r,z) = (4,1,6,1)`
  over `F_65537` was stopped as successive sparse matrices grew to `675,752 → 2,502,378 → 7,780,680`
  columns (validated 8-thread replay). "This experiment did not reach integration; its failure is no
  obstruction to a certificate."
* **Conceptual gap:** the finite primitive bound applies only after a closed residue tuple is
  specified; "no uniform pole bound for all unknown closed residue tuples has been established", so
  failure of the baseline would not decide the full problem.

### 3.1 Consolidated ledger of what is actually excluded

| claim | scope | status |
|---|---|---|
| `deg D ≤ 4` impossible | all characters of `D`, all `N`, all `Q`, over `C` | **rigorous [computed]**: `E1` space is zero mod `p = 32003`, hence over `Q` and `C` |
| `deg D ∈ {5,6,7}` forces `D ∈ J(F)` | `Y`-even `D` [reported, sibling: kernel of `[A]·` on `(S/J)_d^{Y-even}` is `0` for `d ≤ 7`, `1` at `d = 8`, two primes ≈ 2³⁰, 2³¹]; all four characters for the degrees listed in §4 [computed] | rigorous as a statement about `E1` (modular kernel bounds the characteristic-zero kernel); **the `E2` step for `D ∈ J_d` is not closed** |
| 61 degree-6 "natural supports" (`pow··` sums of 2–4 sixth powers, `Fsub·` coordinate sub-equations of `F` with a coupling) fail | full `(N,Q)` system, `p = 268435399` | [reported, `DiamondJM`, public: https://github.com/DiamondJM/hodgeHoldout `sweep/logs/reflection-coupled-residual.log`, commit `af0589e`, 2026-10-23] |
| 42 degree-4 denominators (cones over Fermat/mixed quartics with/without `Y2`, binary quartics in `(X0,Y0)`, `(X0,Y2)`, `(Y0,Y2)`, the `u`-coefficient `X0X1(Y0²+2Y1²)`, the `Y0Y1`-cofactor `C(X) + X0X1Y0² + X1X2Y1²`) fail `E1` | [reported, sibling scratch `log.jsonl`, engine v2] | subsumed by the rigorous `deg ≤ 4` exclusion |
| time-only completion of `fqlx`'s currents; constant combinations of five root/`B²` functions; the selected `(C^{(0)}+C^{(1)})` combination | one presentation on one quartic family | [reported, 5-adic certificates]; **not exclusions of certificates** |

---

## 4. What remains unexplored

1. **Closing degrees 5–7 rigorously.** Since `E1` forces `D = Σ_i c_i∂_iF` with `c_i ∈ S_{d−5}`, the
   general `N` is `N_i = c_iA + Σ_j K_{ij}∂_jF + m_iF` with `K` antisymmetric (`K_{ij} ∈ S_{d−4}`) and
   `m_i ∈ S_{d−5}`; `E2` becomes bilinear in `(c; K, m, Q)`. Parameter counts: `d = 5`: `c` 6,
   `K` 90, `m` 6, `Q` 126, rows 20,349; `d = 6`: 36 / 315 / 36 / 462, rows 33,649; `d = 7`:
   126 / 840 / 126 / 1,287, rows 53,130. Nobody has run this (the sibling's `cands_J.txt` tests only
   ~20 hand-picked `c`; `batch_J` had not run when this was written).
2. **The degree-8 class `D8`.** The unique non-Jacobian `Y`-even `E1`-admissible class at `d = 8`
   (`|S_8^{Y-even}| = 651`, `rank J_8^{even} = 168`, `E1`-solutions `169`). Not yet tested in the full
   system (`12,012 + 3,003` unknowns, `80,730` rows); not yet understood geometrically (what is
   `X ∩ V(D8)`? is it a cone, is it non-normal, does it have `p_g > 0` components, is it `r`/`ρ`-stable,
   does it factor?). A rational representative (two primes `2³¹−1`, `2³¹−19`, CRT + rational reconstruction) is in
   `/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json`; it is dense (≈ 50 terms, all
   `Y`-even and `Y2`-even, denominators powers of 2, 3, 7), beginning
   `−75/196·X0⁸ − 10609/1764·X0⁶X1² − X0⁶Y1² − 625/4116·X0⁵X1Y0Y1 − 25/441·X0⁴X1²Y0² − 25/294·X0⁴X1X2³ + 5/14·X0⁴X2²Y0Y1 − …`. Also untested: the 169-parameter family `D8 + J_8^{even}` (bilinear), and the
   `Y`-odd characters at `d = 8` (my §4 table, if completed).
3. **Degrees ≥ 9**: nothing at all (no `E1` sweep, no fixed-`D` test, no structured ansatz). The sibling
   queued `sweep_fast.py 9..11` (`Y`-even `E1` kernel dimensions) and a control sweep for the companion
   target `X0³Y2³` (`6..12`), outputs not yet present.
4. **Deformation/continuation from Fermat.** `F_c = Σ Z⁶ + c·P` is a pencil through the Fermat
   sextic (`c = 0`), where `A = X0³Y0³` has explicit small certificates (e.g. the degree-6 "Fermat
   join" `D = X0⁶+X1⁶+X2⁶`, used as a positive control in `DiamondJM`'s sweep). No one has tried to
   deform a Fermat certificate along `c` (first-order system: `D_1 ∈ S_6`, `N_1 ∈ S_7⁶`, `Q_1 ∈ S_6`,
   33,649 rows — tiny) and continue to `c = 5/7`.
5. **Non-monomial linear automorphisms** of `F` have not been searched (a hidden reflection in the
   `(Y0,Y1)`-plane would serve half of the target by the reflection theorem). Monomial ones are
   excluded (§1.3).
6. **Covers and exchanges.** No systematic search for a double cover `Z → X` with an exchange
   involution satisfying the §4.6 criterion with `g = r` (the only quotient with known support). The
   Fermat-type covers over `(Y0:Y1)` fail (§1.4); covers over `X`-ratios, over `(Y_i : Y2)` with twisted
   exchanges, or fibre products with the branch divisor `B_W` were never examined.
7. **The branch divisor geometry.** `why.md` explicitly asks "whether a correspondence involving the
   branch locus, a different finite cover, or a construction with nontrivial reflection character
   transfers support to this complementary component". The twisted description
   `H^{3,1}(X)^{r=−1} ≅ H¹(W, Ω³_W(log B_W) ⊗ O_W(−1))` (cyclic-cover decomposition, `L = O_W(1)`,
   `L² = O(B_W)`) **[inference]** has not been used by anyone; in particular no one has looked for
   divisors of `W` swept by rational curves over which the double cover splits (curves meeting `B_W`
   with even multiplicities), whose preimages in `X` are `r`-non-trivially acted-on divisors with
   possibly anti-invariant Gysin image.
8. **Structured denominators adapted to the sparse presentation** (`p = Y0Y1`, `e = X0X1`, `C(X)`,
   `G(Y)`, `u = Y2²`, the partials, products of the quotient's normal derivative `∂_0H_c` with
   low-degree polynomials in `u`) were never swept beyond degree 4 (the `cands_d6.txt` list of 16
   degree-6 cones/binary/`F`-minus-block denominators was queued, not run).
9. **The incidence-curve baseline** (C2) is stuck on effective arithmetic of one component. Alternatives
   never tried: specialize the three quartic parameters to random values over a prime field, compute
   the normalization/integral basis of the specialized incidence curve with Singular (`normal`,
   `integralBasis`) or Magma, and solve the linear equation in `L(65H_∞)` there (a negative answer at a
   random specialization, with a semicontinuity argument, would settle the baseline).
10. **Hodge-theoretic diagnostics of the genus-one route.** Whether the remaining "holomorphic
    coefficient" `g = −½𝒥̇` (R1) is cohomologically exact on the relevant cover was never tested
    numerically; a nonzero transcendental period would show the compact/genus-one presentation cannot
    close without the residue currents, redirecting effort.

---

## 5. Candidate attack plan

Sizes use `q = 1`, `deg N = d + 1`, `deg Q = 2d − 6`. "Chance" is my subjective probability that the
experiment *solves the instance*; most are more valuable as information.

**E1. Close degrees 5–7 (D ∈ J sector).** Idea: `E1` has been reduced to `D ∈ J_d`; finish `E2`.
Formulation: `D = Σ c_i∂_iF`, `N_i = c_iA + Σ_j K_{ij}∂_jF + m_iF`, impose
`B_D(N) − FQ = 0` in `S_{2d+6}`; treat `c` as symbolic and `(K,m,Q)` as linear unknowns; eliminate
linearly to get necessary polynomial conditions on `c`; Gröbner basis over `F_32003` (Singular in the
`ddg-cas` image), then lift (the conditions are homogeneous in `c`, so a modular empty variety proves
emptiness over `C` by the usual specialization argument once the elimination is exact). Sizes as in
§4.1 (largest: `d = 7`, 126 `c`-parameters, 2,253 linear unknowns, 53,130 rows). Expected outcome:
exclusion. Chance of a solve: < 1%. Value: makes "`deg D ≥ 8`" a theorem for this instance.

**E2. The degree-8 class `D8`.** Idea: it is the *only* candidate below degree 9 outside `J`. Steps:
(a) rational reconstruction (done, see §4.2); (b) fixed-`D` solve (`12,012 + 3,003` unknowns,
`80,730` rows; the sibling engine handles degree 6 in 10–40 s, degree 8 in minutes); (c) if
inconsistent, the family `D8 + Σ_i c_i∂_iF`, `c_i ∈ S_3` (168 `Y`-even parameters) via the same
elimination as E1; (d) geometry of `X ∩ V(D8)`: factorization over `Q̄`, singular locus, whether it
is a cone (independent of some linear forms), `r`- and `ρ`-stability (it is `Y`-even and `Y2`-even by
construction). Chance: ≈ 3% (a unique canonical object sometimes is the answer; but the solved
one-reflection companion needed degree 30 per half). Cost: an hour.

**E3. `E1` sweep by character to `d ≈ 12–14`, then geometry of kernel classes.** Idea: the map
`[A]: (S/J)_d → (S/J)_{d+6}` has growing kernel; its elements are the only `E1`-admissible
non-Jacobian denominators. My script (`/tmp/hodge/scratch/rcr-dossier/e1dims.py`, dense NumPy mod
32003, split by Klein character) reaches `d = 9` in ~1 h; Singular or `python-flint` sparse ranks reach
`d = 14` (`33,649 × 51,408` rows/columns at `d = 12`). For each new kernel class test the fixed-`D`
system and inspect `X ∩ V(D)` for the (3,1)-support signature (non-normal, or a component with
non-rational singularities, cf. Chapter 3 §3.1). Chance: ≈ 2%; information value: high (it tells
where certificates *can* live; in the skew instance the kernel is `9, 144, 388, 699, 1029` at
`d = 7..11`).

**E4. Fermat-pencil continuation.** Idea: `F_c = Fermat + c·P`; at `c = 0` the target has explicit
certificates (`D_0 = X0⁶+X1⁶+X2⁶`, strict, per `DiamondJM`'s control; also `D_0 = X0⁶ + Y0⁶`, etc.).
Seek a formal family `(D(c), N(c), Q(c)) = Σ c^k(D_k, N_k, Q_k)` of fixed degree `d = 6` (or 12 after
multiplying `D_0` by a second Fermat certificate denominator). Order-`k` equation:

    A·Σ_{i+j=k} D_iD_j − Σ_{i+j=k} D_i N_j(F_0) − Σ_{i+j=k−1} D_i N_j(P)
      + Σ_{i+j=k} F_0·B_{D_i}(N_j) + Σ_{i+j=k−1} P·B_{D_i}(N_j) = Σ F_0²Q_k + 2F_0P·Q_{k−1} + P²Q_{k−2};

each order is the *same* linear operator (the linearization at `(D_0,N_0,Q_0)`) with a new right-hand
side: `D_k ∈ S_6` (462), `N_k ∈ S_7⁶` (4,752), `Q_k ∈ S_6` (462), rows `S_18` (33,649). If the
obstruction (right-hand side modulo the image) vanishes at every order, rationally reconstruct the
series in `c` (Padé) and evaluate at `c = 5/7`, checking coprimality there. If it fails at some order,
retry with other Fermat certificates for `A` (there are many: Fermat has three `Y`-reflections and the
full `μ6⁶ ⋊ S_6` symmetry) or with the product of two of them. Chance: ≈ 5% (the pencil may have
Picard/support jumping at `c = 0` only; but the experiment is cheap and never tried on any instance of
this family). Cost: hours.

**E5. Full linear automorphism group.** Idea: find hidden non-monomial symmetries. Restricting to the
block form forced by the characters is *not* justified a priori, so instead solve `F∘g = λF` for
`g ∈ GL_6` by Gröbner/numerical algebraic geometry on the 36 + 1 unknowns with the Fermat-part
constraints (`F|_{Y=0}` and `F|_{X=0}` are Fermat, whose automorphism groups are the finite monomial
groups; any `g` must map the fixed-point data of the Klein group to itself, which already forces `g` to
preserve `{Y2 = 0}` and `{Y0 = Y1 = 0}`, hence a block-triangular form with blocks `X ⊕ (Y0,Y1) ⊕ Y2`
**[inference]**). Then the problem reduces to `g_X ∈ Aut(ΣX⁶)` (monomial, 1296 elements),
`g_{Y01} ∈ Aut(Y0⁶+Y1⁶)` (monomial, 72 elements), `g_{Y2} ∈ μ6`, and the monomial check of §1.3 already
gives **only the Klein group**. So E5 is essentially done modulo the block-form lemma; I record it as a
small proof obligation rather than an experiment. Chance of a new symmetry: ≈ 0.

**E6. Specialized incidence-curve baseline (rescue of C2).** Idea: avoid F4SAT. Specialize
`(p,q,r)` of the centered quartic to random elements of `F_ℓ` (large `ℓ`), obtain the specialized
incidence curve as a plane or space curve, compute its normalization and the Riemann–Roch space
`L(65H_∞)` with Singular (`normal`, `integralBasis`, Hess/`brill_noether` for RR spaces) or Magma, and
solve `δ𝒯(f) − 𝒯(δ)f = N` there. If solvable at several random specializations, lift (Hensel +
rational reconstruction in the parameters) to the generic curve; if unsolvable at one good
specialization, the generic baseline is unsolvable (solvability of a linear system over the function
field specializes well away from a proper closed set). Then the three residue functions remain; a
*bounded* ansatz for them (poles only at the known 16+ places at infinity and the source-zero loci,
orders ≤ 31) would be the natural next step, acknowledged by C2 as non-exhaustive. Chance: ≈ 5%
(this route has a lot of machinery already verified; the risk is that the baseline fails and the
residue tuple is genuinely unbounded). Cost: days, needs the packet's exact data.

**E7. Structured degree-6/8/12 denominators from the sparse presentation.** Families (each a small
linear or bilinear system): (i) cones over plane sextics of geometric genus ≥ 1 in coordinate planes
containing `X0` and `Y0` (`X0⁶+Y0⁶+Y2⁶`, `X0⁶+Y0⁶+Y1⁶`, the Hesse-web members in `cands_d6.txt`);
(ii) `D = (∂_0H_c)(x,Y2²)·L(u, …)` mimicking the companion's legal denominator `D_H = (∂_0H_c)·d`
with `d` of low degree in the quotient coordinates — the companion's actual denominator has this
shape, and a residual-sector certificate, if its denominator is `r`-invariant, is also of the form
`φ(D_H)`; (iii) products `D = D_+ · r(D_+)` for `Y2`-odd `D_+` (the Chapter 1 character-averaging
form of any certificate: if `(D,N)` works then so does the averaged `D·r(D)·ρ(D)·τ(D)`, so one may
assume `D` Klein-invariant at the cost of degree ×4); (iv) `D = Y2²·q(X,Y) + u`-type polynomials
reflecting the branch divisor. Sizes: fixed-`D` tests at `d = 6, 8, 12` (table §1.2). Chance: ≈ 2%.

**E8. Exchange covers for the §4.6 criterion with `g = r`.** Idea: the only missing input of the
descent criterion is a double cover `Z → X` with deck `d` and an involution `σ` such that
`e = σdσ` covers `r`, plus support of `Z/⟨r_Z,σ⟩`. Systematic cheap search: adjoin `T` with
`T² = h(X,Y)` for `h` ranging over `Y2`-even quadrics/quartics/sextics built from the invariants
(`Y2²`, `Y0Y1`, `Y0²`, `Y1²`, `X`-monomials) and look for a linear/monomial involution `σ` of
`(X, Y0, Y1, Y2, T)` preserving `{F = 0, T² = h}` and swapping `T ↔ Y2` up to scalars/twists (as in the
holdout exchange `y ↔ T`). Each candidate is a finite exact check (minutes). Then one still needs the
quotient `Z/⟨r_Z,σ⟩` to be Fano/uniruled (weighted adjunction) — the holdout pattern. Chance: ≈ 3%;
this is the *only* route in the book's toolbox that reaches a complementary character on a fixed
variety ("the cover construction in §2.7 can reach further characters because its useful reflection
acts on an auxiliary variety rather than on the original fourfold", §3.3).

**E9. Anti-invariant support from the branch divisor / twisted sheaf on `W` [conjectural,
nonconstructive].** Idea: `H^{3,1}(X)^{r=−1} ≅ H¹(W, Ω³_W(log B_W)(−1))`. Divisors `E ⊂ X` with
`r`-anti-invariant Gysin image come from `r`-stable divisors `E` whose resolution has
`H^{2,0}(Ẽ)^{r=−1} ≠ 0`, e.g. `E = π⁻¹(E_W)` with `E_W ⊂ W` swept by rational curves tangent to `B_W`
to even order (so that `π⁻¹(curve)` splits into two components exchanged by `r`), parametrized by a
surface `S_W` with `p_g(S_W^{twisted}) > 0`. Experiment: on the modular model `H_0 = Σx⁶ + u³`, look
for low-degree weighted rational curves with `u(t)` a perfect square (then `Y2 = √u` is rational on
the curve and the curve lifts to `X`), count the dimension of the family, and check whether it sweeps a
divisor of `W` (dimension 3 family modulo reparametrization); then lift as in the companion. A lifted
curve family on `X` of dimension 3 sweeping a divisor `E` would give the support candidate; the
cubic-cover analysis of §4.7 also suggests computing the discriminant surface
`Δ(a,b,c,d) = 0` of `a + bt + ct² + dt³` with `b = c_0·p·C(X)` (a product!) and `c`-term
`X0X1G(Y) + X1X2Y0Y1³`, whose special shape might make the `S_3` cover's standard summand partly
algebraic here. Chance: unknown (≈ 10% that some such divisor exists; much lower that it is found
and made explicit soon); value: it is the only *existence*-type idea specific to this variety.

**E10. Mechanization path.** If E4, E7 or E2 produce an explicit `(D,N,Q)`: generate
`Witness.lean/Solution.lean` with the accepted template (`hodge/leangen.py`, sandbox check verified
2026-11-01, ~5 min for a regenerated accepted witness); kernel-checkability requires either small
degree (≤ ~12 with sparse coefficients) or a factored presentation as in `one-reflection` (degree 60
through a five-variable identity). If E6/E8/E9 produce a serving geometry with a polynomial curve
family and a nonsingular incidence seed, the shared modules (`HoldoutCover`, `PolynomialCurves`,
`EllipticCompanions`, `ResidueTransport`, `CertificateOperations`) supply lifting, integration, trace,
clearing, addition and homogeneous extraction, exactly as in PR #150; the instance must supply the
serving correspondence, the exact source transport, and coprimality.

---

## 6. Feasibility: explicit witness vs. nonconstructive argument; Lean requirements

* **Explicit witness.** Degrees ≤ 4 are impossible; 5–7 require `D ∈ J` (one small bilinear
  computation from closure); degree 8 has exactly one non-Jacobian `Y`-even candidate class. Beyond
  that nothing is known. The nearest explicit certificate in the catalogue (`one-reflection`) has
  `deg D = 60` and was Lean-checkable only through a factored five-variable identity; the companion on
  *this* variety exports no degree at all. If a low-degree witness exists it will most likely be found
  by E2/E3/E4; a high-degree one only by a geometric construction presented in factored form.
  **Assessment: an explicit witness by bounded search is unlikely (≤ 10% overall for E1–E4, E7).**
* **Nonconstructive existence.** Predicted by the general Hodge conjecture for the `τ`-invariant
  sub-Hodge structure **[inference]**, but no mechanism in the book serves the `r`-anti-invariant,
  `ρ`-anti-invariant sector of this `F`: the Fano quotient `W` serves only the `r`-invariant 82; `X/ρ`
  and `X/τ` are not uniruled; there is no `μ3`/`μ6` action; no exchange cover is known. The §4.4
  generic theorem does not apply to an isolated instance. The thread's genus-one / incidence-curve
  programme is a *computational* route to a primitive on a specific quartic family, not an existence
  argument; it has no uniform bound for the residue unknowns.
* **What a Lean proof must contain** (CONTRIBUTING, `formal/README.md`): a closed theorem
  `Submission.certificate : HodgeComputer.HasCertificate HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q`
  with axioms only `propext, Classical.choice, Quot.sound`; `Verify` requires `D ≠ 0`, homogeneity
  (`deg N_i = deg D + 1`), `IsRelPrime F D` over `C` (not inferred from rational irreducibility; the
  companion proves it by monic division in `x0` plus complex irreducibility of `F` via Eisenstein),
  and the exact identity `A·D² − D·N(F) + F·B_D(N) = F²·Q` with `Q ≠ 0` allowed. Any construction over
  a number field or a function field descends by normalized trace (Chapter 1); existential proofs
  follow the PR #150 pattern (integral incidence seed with invertible Jacobian, étale lifting, polynomial
  integration, trace, clearing, transport, homogeneous extraction) and import the maintained shared
  modules. Files: `proof.json` (`schemaVersion 2`, `theorem`), `Solution.lean`, `solution.md`, optional
  `witness.json`.
* **Bottom line.** The instance is blocked at the geometry-finding step ("cross a character boundary
  after the curve-finding and integration steps are already understood", `why.md`). The cheapest
  informative next steps are E2 (test `D8`), E1 (close degrees 5–7), E4 (Fermat-pencil continuation,
  never tried, tiny systems), and E8 (exchange-cover search); the only existence-type idea specific to
  this variety is E9.

---

## Appendix A. Verification log (this session) **[computed]**

Scripts and outputs in `/tmp/hodge/scratch/rcr-dossier/` (`e1dims.py`, `e1dims.log`); SymPy/NumPy
snippets inline in the session; Python run with `-I`.

* Monomial support (13 monomials), `Y`/`Y2`-degrees, Smith normal form `(1,1,1,2,2,0)` of the
  exponent-difference lattice, `μ12` search (4 projective diagonal symmetries), sign symmetries with
  determinant and residue characters, absence of permutation symmetries (support-preserving
  permutations: identity only), `F(Y,X) ≠ F`.
* Monomial counts: 138 `Y`-odd/`Y2`-even sextic monomials; 158 `Y`-even/`Y2`-even (= the book's
  coefficient dimension of the one-reflection locus).
* Equivariant Hilbert series of `S/J` under the Klein group: sector dimensions in degrees 6/12/18 as in
  §1.3 (`126/492/126` for this target; `82/374/82` for the companion; totals `426/1751/426`).
* `E1`-admissible denominator dimensions `dim{D ∈ S_d^χ : A·D ∈ J_{d+6}}` mod `p = 32003`, by Klein
  character `χ = (Y2-parity, (Y0+Y1)-parity)` of `D` (kernel over `Q` ⊆ kernel mod `p`, so zeros are
  rigorous and positive values are upper bounds):

  | `d` | `χ=(0,0)` | `χ=(0,1)` | `χ=(1,0)` | `χ=(1,1)` | total | `dim J_d` | non-Jacobian |
  |---|---|---|---|---|---|---|---|
  | 4 | 0 (of 48) | 0 (of 38) | 0 (of 22) | 0 (of 18) | **0** | 0 | 0 |
  | 5 | 3 (of 88) | 2 (of 78) | 1 (of 48) | 0 (of 38) | **6** | 6 | 0 |
  | 6 | 14 (of 158) | 12 (of 138) | 6 (of 88) | 4 (of 78) | **36** | 36 | 0 |
  | 7 | (job still running at write time; sibling: `Y`-even kernel mod `J` = 0) | | | | | 126 | 0 in the `Y`-even sector |
  | 8 | (not run) | | | | | 336 | 1 in the `Y`-even sector (sibling, two primes) |

  (`dim J_d = 6·dim S_{d−5}` for `d ≤ 9`, there being no Koszul syzygies below degree 10.)
  Sibling scratch (`sweep_e1.out`, `k8_d8.out`, primes `≈ 2³⁰` and `2³¹−1`): kernel of `[A]` on
  `(S/J)_d^{Y-even}` is `0` for `d = 1..7` (`dim (S/J)_d^{Y-even} = 3, 12, 28, 66, 123, 218, 333`) and
  `1` at `d = 8` (`dim 483`; `E1`-solutions including `J`: `169 = 168 + 1`).

## Appendix B. Links

* Instance: `data/reflection-coupled-residual/{instance.rec, why.md, provenance.json, checks.json}`.
* Discussion #154: https://github.com/hodge-computer/sextics/discussions/154
  (packets: C1 zip `…/33743911/reflection_residual_genus_one_continuation.zip`; R1 gist
  `https://gist.github.com/fqlx/ae3df1cc3d0a3857b7b69f2ad80fc983`; R2 zip
  `…/33835113/reflection-coupled-residual-radial-audit-2026-10-22.zip`; C2 zip
  `…/34087387/reflection-residual-research.zip`; prefix `https://github.com/user-attachments/files/`).
* Companion: `data/reflection-coupled-fano/solution.md`, PR #150
  https://github.com/hodge-computer/sextics/pull/150, discussion #111
  https://github.com/hodge-computer/sextics/discussions/111, CI run
  https://github.com/hodge-computer/sextics/actions/runs/39376055433.
* Machinery PRs cited in the thread: #128 (holdout), #161 (holdout-bundle), #165 (cycle-three), #151
  (holdout-ternary-deformation), #144 (one-reflection-y2-sector).
* `DiamondJM` sweep: https://github.com/DiamondJM/hodgeHoldout (`sweep/RESULTS.md`,
  `sweep/logs/reflection-coupled-residual.log`).
* Book: `docs/03-geometric-support.md` §3.1, §3.3 (incl. "An algebraic realization by weighted
  quartics"), `docs/04-known-locus-and-frontier.md` §4.4 ("The intermediate one-reflection locus"),
  §4.6, §4.7; `docs/02-explicit-constructions.md` §2.7–2.8.
* Sibling scratch (local, unpublished): `/tmp/hodge/results/reflection-coupled-residual/`
  (`log.jsonl`, `sweep_e1.{out,json}`, `k8.py`, `k8_d8.out`, `cands_*.txt`, `calib.out`); engine
  `~/devel/sextics-work/hodge/` (v2 after the 17:25 sign-bug fix; earlier results unreliable).
* Sibling dossier for the closest open instance: `/tmp/hodge/notes/dossier-one-reflection-skew.md`.
