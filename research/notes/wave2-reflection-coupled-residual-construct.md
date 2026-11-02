# Wave 2 — lane `reflection-coupled-residual` — lens CONSTRUCTION

STATUS: no witness. Structural results: every published mechanism fails at a named hypothesis for the residual
`(r,ρ) = (−,−)` sector; the sector is `H¹(P¹, j_*V₁₁)` for the Klein pencil with `V₁₁` of rank 186, Hodge numbers
`(3,90,90,3)`, and the `∇`-span of its `(3,0)`-forms is all of `V₁₁` (computed, rank 186) — so no level-one
(Abel–Jacobi) support exists and any support must be a horizontal divisor with `(−,−)` `(2,0)`-forms (the twisted-slice
mechanism M5, whose explicit integration is period-obstructed). Fermat-pencil continuation fails at first order. The
hidden-reflection search (`X` and the six pencils' generic fibres) is formulated, validated on a toy, and still running
on the remote at session end (job `construct-rcr-all3`; pull with `construct/finalize_pull.py`).

Log (one JSON line per experiment): `/tmp/hodge/results/reflection-coupled-residual/wave2-construct.log.jsonl`
Scripts/outputs: `/tmp/hodge/results/reflection-coupled-residual/construct/`
Conventions as in the dossier (`/tmp/hodge/notes/dossier-reflection-coupled-residual.md`): Klein group
`{1, r, ρ, τ}`, `r: Y2 ↦ −Y2`, `ρ: (Y0,Y1) ↦ −(Y0,Y1)`, `τ = rρ`; residue characters of the target
`A = X0³Y0³`: `(r, ρ, τ) = (−1, −1, +1)`; residual sector `(126, 492, 126)`; companion sector `(82, 374, 82)`.

## 0. Ledger: mechanisms versus hypotheses (status after this session)

| # | mechanism (task item) | hypothesis needed | status for this `(F, A)` |
|---|---|---|---|
| M1 | §4.6 deck+exchange cover, `g = r`, `X/g = W` Fano | a cover `Z → X` with an involution `σ` **not** lifting an automorphism of `X`; coordinate pencils: fibre separately even in two *new* linear forms | **fails for all 6 coordinate pencils** (proved by inspection of the fibre equations, §2.1); non-coordinate fibre reflections: search running in λ-eliminated form (job `construct-rcr-all3`, Exp B; two heavier formulations stalled on the loaded box; no fibre case had finished at session end) |
| M2 | Galois closure of `X → P⁴_x` (item iii) | an intermediate quotient `Z̃/K` that is uniruled and has `(V ⊗ sgn_ρ)^K ≠ 0` | **fails**: Galois group is the full hyperoctahedral `B₃ × Z/2` [computed]; the target lives in the reflection representation `V ⊗ sgn_ρ`; the three maximal `K` reaching it give `X/τ` (κ = 0), a `(4,6)`-CI in `P(1⁶,2)` mod an involution with ≥ 9 invariant canonical forms, a degree-8 hypersurface in `P⁵` mod an involution with ≥ 9 invariant canonical forms; none uniruled (§2.2) |
| M3 | Fano / surface-fibration companions of Klein quotients (item ii) | uniruled quotient, or fibration over a surface with `p_g = q = 0` fibres | **fails**: only `W = X/r` and `X/⟨r,ρ⟩` are Fano (target not invariant); `X/ρ`, `X/τ` have `K ≡ 0`; the natural surface fibrations of `X/τ` over `P²_Y`, `P²_X` have fibres of weighted degree 6 in `P(1,1,1,2)` with `p_g = 3` (§2.3) |
| M4 | q-raising / Gauss–Manin from the companion (item iv) | a *family* `F + tC` on which the companion's certificate deforms, with `C` `r`-odd (to flip the character) | **fails**: `C` `r`-odd leaves the one-reflection locus, where nothing is supported; `C` `r`-even preserves the character (§2.4) |
| M5 | twisted slices of the companion's curve family on `W` (item v) | a slice `S` of the rational-curve family on `W` whose hyperelliptic lift carries anti-invariant 2-forms, i.e. `H⁰(S, Ω¹_S ⊗ λ_hyp) ≠ 0`, with nonzero Gysin image | **the only mechanism on this variety that can produce `(−,−)` classes**; explicit integration is obstructed by `H³(S, R¹p_*L) ≠ 0` (positive-genus lifts) — exactly the period obstruction met in discussion #154; a finite-field incidence seed cannot be used because the lifted curves have genus `e − 1 ≥ 1` for every weighted degree `e ≥ 2`, and `e = 1` (weighted lines) gives `r`-invariant classes only (§3) |
| M6 | Fermat-pencil continuation (item iv) | formal deformation of the degree-6 Fermat certificate along `F_c = Fermat + cP` | **fails at first order** for the strict subset-Euler base `D₀ = X0⁶+X1⁶+X2⁶`, `N₀ = A·X_i/6`: the linearised system (9306 rows, 1664 unknowns in the Klein `(+,+)` sector) has rank 1484 but the augmented rank is 1485 (`p = 2³¹−1`, job `construct-rcr-all`); the `Y`-block base is being tested |
| M8 | nonlinear polynomial Piola / finite Fano maps (§3.4; "`F = ΣQ_i³`", block maps) | `F = G(Φ)` with `Φ` polynomial of equal dimension; served sector = `det(dΦ)·B(Φ)` | **only `Φ = (x, Y2²)` exists** (density `2Y2`, the companion): `F` is not a polynomial in `(X, Y0², Y1², Y2²)` (the `Y0Y1` terms) nor in `(X, Y0Y1, q(Y0,Y1), Y2²)` for a single quadratic `q` (`Y0²Y2²` vs `2Y1²Y2²`, `Y0³Y1·X0X1` vs `Y0Y1³·X1X2` need three invariant quadrics), and `F` is not a sum of cubes of quadrics beyond the Fermat part |
| M7 | pencil `(Y0:Y1)` and the variation `V₁₁` (new) | Abel–Jacobi/normal-function supports would need a level-one sub-variation `V' ⊂ V₁₁` | the residual sector **is** `H¹(P¹, j_*V₁₁)`, `V₁₁` of rank 186 with Hodge numbers `(3,90,90,3)` [computed]; **the `∇`-span of the three `(3,0)`-forms is all of `V₁₁` (rank 186, computed mod `p`)**, so `V₁₁` has no proper sub-local-system containing `F³` and **no level-one sub-variation** (§4) ⇒ Abel–Jacobi-type supports are impossible; supports must be horizontal divisors carrying `(−,−)` `(2,0)`-forms, consistent with M5 and excluding the "irregular surfaces in the fibres" route |

## 1. Exact data used

`F = Σ Z_i⁶ + c·P`, `c = 5/7`, `P = X0³X1Y0Y1 + X0X1Y0³Y1 + X0X1Y0²Y2² + 2X0X1Y1²Y2² + 3X0X2³Y0Y1 + 2X1³X2Y0Y1 + X1X2Y0Y1³`.
As a cubic in `u = Y2²` over `C(x)`, `x = (X0,X1,X2,Y0,Y1)`:

    F = u³ + L(x)·u + M(x),   L = c·X0X1(Y0² + 2Y1²),
    M = ΣX_i⁶ + Y0⁶ + Y1⁶ + c·(Y0Y1·C(X) + X0X1Y0³Y1 + X1X2Y0Y1³),   C(X) = X0³X1 + 3X0X2³ + 2X1³X2.

[computed] `M` is irreducible over `Q`; `disc = −4L³ − 27M²` is irreducible of degree 12 (both checked
on a random line); the Galois group of `Y2⁶ + L·Y2² + M` over `Q(x)` specialised at two random rational
points is the transitive group of order 48 (`C2 ≀ S3`, sympy `galois_group`). Hence the Galois closure
`Z̃` of `X → P⁴_x` has group `G = B₃ = (Z/2)³ ⋊ S₃` (order 48), and `ρ` lifts to a central
involution `ρ̃` commuting with `G`.

## 2. Negative results (proved modulo the stated computations)

### 2.1 No exchange on any coordinate pencil (item i)
Fibre of the pencil `Y0 = λY1` in coordinates `(X, T, Y2)` with `(Y0, Y1) = (α, β)T`:

    Φ = ΣX⁶ + T⁶ + Y2⁶ + c[ αβ·T²·C(X) + (α³β·X0X1 + αβ³·X1X2)·T⁴ + (α² + 2β²)·X0X1·T²Y2² ].

It is separately even in `T` and in `Y2` (these are the inherited `ρ|` and `r|`; `ρ` restricted to the
hyperplane `Y0 = λY1` has a one-dimensional `(−1)`-eigenspace, so it *is* a reflection of the fibre), and
the monomials `T²C(X)`, `T⁴X0X1`, `T⁴X1X2` have no `Y2²`, `Y2⁴` partners, so no linear involution of the
`(T, Y2)`-plane other than the inherited ones preserves `Φ`. The pencils `Y0/Y2`, `Y1/Y2` fail because
`Y0³Y1·X0X1` and `Y0Y1³·X1X2` have different `X`-coefficients and because of the unpaired `Y2⁴`/`T⁴`
terms; the `X`-pencils fail because `C(X)` restricted to `X_i = λX_j` acquires a pure fourth power with
no partner (all seven cases written out in the session; each is a one-line inspection). The quotient of
every fibre by its inherited Klein group is the Fano threefold

    V_λ = { ΣX⁶ + s³ + u³ + c[αβ·s·C(X) + (α³βX0X1 + αβ³X1X2)·s² + (α²+2β²)·X0X1·s·u] = 0 } ⊂ P(1,1,1,2,2),

which only reproduces the `(+,+)` (companion) sector (§4.6 with inherited `σ` gives `g = 1`).

### 2.2 Averaging ideal of the `B₃`-closure (item iii)
`Q[G/H]` for `H = Stab(z₁) ≅ D₄` (so `Z̃/H = X`) decomposes as `(3,∅) ⊕ (21,∅) ⊕ (2,1)`; the `r`-sign on
the `H`-fixed vectors is `+1, +1, −1`. Hence `H⁴(X)^{r=+}` is the `(3,∅) ⊕ (21,∅)` isotypic part (served
by `W = Z̃/⟨H, r̃⟩`, `Q[G/⟨H,r̃⟩] = (3,∅) ⊕ (21,∅)`), and the residual sector is the isotypic part of
`V := (2,1)`, the 3-dimensional reflection representation of `B₃`, tensored with `sgn_ρ`. A subgroup
`K ⊂ G × ⟨ρ̃⟩` reaches it iff `K ⊂ Stab(v) ⊔ Stab(v)·g·ρ̃` for some `v ∈ Q³` and `g` with `gv = −v`.
The maximal such `K` (one per orbit type of `v`):

| `v` | `K` | `Z̃/K` | canonical forms / Kodaira |
|---|---|---|---|
| `e₁` | `⟨H, s₁ρ̃⟩` (16) | `X/τ` | `K ≡ 0`, not uniruled (book §4.7) |
| `e₁+e₂` | `⟨s₃, (12), s₁s₂ρ̃⟩` (8) | `{σ₂² − (σ₁² − 2σ₂)² = L, σ₂²(σ₁² − 2σ₂) = M} ⊂ P(1⁶,2)` mod `(σ₁,Y0,Y1) ↦ −(σ₁,Y0,Y1)` (`σ₁ = z₁+z₂`, `σ₂ = z₁z₂`) | `K_Y = O(2)`; invariant sections ↔ quadrics odd in `(Y0,Y1,σ₁)`: `X_iY_j, X_iσ₁` ⇒ `p_g ≥ 9` if singularities are canonical |
| `e₁+e₂+e₃` | `S₃ × ⟨−ρ̃⟩` (12) | `{w⁸ − 8Lw⁴ + 64Mw² + 16L² = 0} ⊂ P⁵` (`w = z₁+z₂+z₃`) mod `(w,Y0,Y1) ↦ −(w,Y0,Y1)` | `K_Y = O(2)`; `p_g ≥ 9` likewise |

Smaller `K` give covers of these. So no uniruled intermediate quotient of the natural Galois closure
reaches `V ⊗ sgn_ρ`: the §4.7 obstruction for the `S₃`-cover of `X/τ` has an exact analogue for the
`B₃`-cover of `W` (the two covers are different; this one is new).

### 2.3 Klein quotients and their fibrations (item ii)
`X/ρ ⊂ P(1,1,1,1,2,2,2)` is the CI `{F̃ = 0, ac = b²}` of degrees `(6,4)` (`a,b,c = Y0², Y0Y1, Y1²`),
`−K = O(0)`; `X/⟨r,ρ⟩ ⊂ P(1,1,1,2,2,2,2)` is a `(6,4)`-CI with `−K = O(1)` (Fano, invariant sector only);
`X/τ` has `K ≡ 0`. The rational maps `X/τ ⇢ P²_Y` and `X/τ ⇢ P²_X` have fibres `{f + s·a + s²·b + s³·κ = 0}
⊂ P(1,1,1,2)` (`s = t²`), weighted degree 6, `K = O(1)`, `p_g = 3`: the surface-fibration theorem
(§3.2.3, needs `p_g = q = 0`) does not apply. Fibres of `X ⇢ P²_Y` are sextic surfaces (`p_g = 10`).

### 2.4 Gauss–Manin / q-raising (item iv)
For `F_t = F + tC`, `∂_t[A'Ω/F_t²] = −2[A'CΩ/F_t³]`. The residue character of `A'C` is `χ(A')·χ(C)`, so
only an `r`-odd even sextic `C` moves the companion target `A' = X0³Y2³` into the residual sector (at
`q = 2`, target `A'C`); the §4.5 identity then needs the companion certificate to deform to first order
along `C`, i.e. the class `[A'Ω/F_t²]` to stay supported for `t ≠ 0` — but `F_t` is no longer `Y2`-even
and nothing is supported there. Invariant `C` preserves the character. The same holds at Fermat: the
pencil `F_c = Fermat + cP` has only the Klein group at the generic point, so no sector theorem applies
off `c = 0`.

## 3. The twisted-slice mechanism and the exact obstruction (item v) — the single most promising route

Let `π: X → W = X/r`, `B_W = W ∩ {u = 0}` (branch divisor, class `O_W(2)`), `L = (π_*Q)^−` the rank-one
local system on `W ∖ B_W` with monodromy `−1` around `B_W`. Then `H⁴(X)^{r=−1} ≅ H⁴(W, j_!L)`-type
twisted cohomology, and the residual `(3,1)` sector is its `ρ`-anti-invariant part. Divisors `E_W ⊂ W`
swept by a 2-parameter slice `S` of rational curves `C_s` of weighted degree `e` lift to
`E = π⁻¹(E_W) ⊂ X`, fibred over `S` in the double covers `π⁻¹(C_s) → C_s ≅ P¹` branched at the
`2e` zeros of `u|_{C_s}`, i.e. hyperelliptic curves of genus `e − 1`. Anti-invariant 2-forms on `Ẽ`:

    H^{2,0}(Ẽ)^{r=−1} ↪ H⁰(S, Ω¹_S ⊗ f_*ω_{Ẽ/S}),   rank f_*ω_{Ẽ/S} = e − 1

(the `f^*Ω²_S` part of `Ω²_Ẽ` is `r`-invariant). Consequences:
* `e = 1` (weighted lines on `W`, a 2-dimensional family by the count `e + 1`): the lifts are the
  `r`-invariant conics of `X`, rank 0 Hodge bundle, **no** anti-invariant classes — the "line divisor"
  of `X` carries only the companion sector.
* `e = 2` (weighted conics, the "divided conic" family of discussion #154): genus-one lifts, line bundle
  `λ = f_*ω`; anti-invariant classes need `H⁰(S, Ω¹_S ⊗ λ) ≠ 0` — a condition on the *slice*, satisfiable
  for sufficiently positive slices (sections of twisted 1-forms), not for the whole family.
* `e = 4` (the companion's quartic family): genus-three lifts, rank-3 Hodge bundle.
* **Why the shared Lean recipe cannot apply**: the polynomial-curve primitive (§2.7, `PolynomialCurves`,
  `HoldoutCover`) needs the pulled-back residue form to be exact on the affine total space `A¹ × S`;
  with the local system `L`, the Leray term `H³(S, R¹p_*L)` with `R¹p_*L = H¹(P¹ ∖ 2e pts, L)` of rank
  `2e − 2` is the obstruction. It vanishes identically only for `e = 1`, which serves the wrong sector.
  So the best possible explicit route is `e = 2` with a *period* obstruction — precisely the state of
  `fqlx`/`benediamond`'s computation (no established nonzero period, no primitive). An existence proof
  via M5 would have to go through the Gysin/Hodge-theoretic statement (`H⁰(S, Ω¹_S ⊗ λ) ≠ 0` plus
  nonvanishing of the Gysin map), not through termwise integration; and because the §4.4 "one class
  suffices" theorem is a statement at the generic point of the one-reflection locus, it would yield the
  residual sector *generically on the 158-dimensional locus*, not at this isolated `F` without an
  explicit transport.

## 4. New structural description: the residual sector as `H¹` of a weight-3 variation (item iv, new)

Blow up the base surface `Σ₀ = X ∩ {Y0 = Y1 = 0}` (the Fermat sextic surface in `P³_{X,Y2}`) to get
`f: X̃ → P¹_{(Y0:Y1)}` with fibres the Klein-symmetric sextic threefolds `X_λ` of §2.1. Leray gives
`H⁴(X̃) = Q ⊕ H¹(P¹, R³f_*Q) ⊕ Q`, and since `Σ₀`'s classes are `ρ`-invariant,

    H⁴(X)^{(r,ρ) = (−,−)} = H¹(P¹, j_*V₁₁),   V₁₁ := (R³f_*Q)^{(−,−)}.

[computed, `p = 32003`, random `λ`; Jacobian ring of `X_λ` by character of the *residue*]

| residue character `(ρ, r)` | `h^{3,0}` | `h^{2,1}` | `h^{1,2}` | `h^{0,3}` | rank | meaning |
|---|---|---|---|---|---|---|
| `(+,+)` | 0 | 41 | 41 | 0 | 82 | `H³` of the Fano `V_λ`; companion sector (82 = rank!) |
| `(−,+)` | 1 | 62 | 62 | 1 | 126 | `H³(X_λ/r) ⊖ H³(V_λ)`; `X_λ/r` is the CY3 sextic in `P(1,1,1,1,2)` (`h^{2,1} = 41 + 62 = 103`) |
| `(+,−)` | 1 | 62 | 62 | 1 | 126 | likewise for `X_λ/ρ` |
| `(−,−)` | **3** | **90** | **90** | **3** | **186** | residual sector: `V₁₁` |

The three `(3,0)`-forms are `Res(X_i Ω/Φ_λ)`, `i = 0,1,2`.

**Proposition (no level-one sub-variation seen by the nodal fibres).** Let `E ⊂ V₁₁` be the span of the
`(−,−)`-projections `δ_s` of the vanishing cycles of the *nodal* singular fibres (free Klein orbits of four
nodes and `r`-fixed pairs; each such fibre contributes one `(−,−)` class and acts on `V₁₁` by the symplectic
transvection in it). Every sub-variation `V'' ⊂ V₁₁` containing the three `(3,0)`-forms `ω_i = Res(X_iΩ/Φ_λ)`
contains `E`. Consequently a level-one sub-variation `V'` (type `(0,a,a,0)`) of `V₁₁` is contained in
`E^⊥`, the part of `V₁₁` with trivial local monodromy at every nodal fibre (nontrivial monodromy only at the
six degenerate fibres `λ⁶ = −1`, whose singular point `(0:0:0:1:ζ:0)` has local equation
`cλ·X1·(λ²X0 + X2) + 3cλ·X0⁴ + 2X0⁶ + Y2⁶ + …` (Morse in `(X1, λ²X0+X2)`, semi-quasihomogeneous principal part
`X0⁴ + Y2⁶` in the other two variables): corank 2, Milnor number `3·5 = 15`). **Computed (`gm_rank.py`, job `construct-rcr-gm2`, `p = 2²⁶ − 5`, random `λ₀`): the `∇`-span of the three
`(3,0)`-forms has rank 186 = rank `V₁₁`** (the span grows by exactly 3 per derivative order up to order 61); since a
modular rank bounds the characteristic-zero rank from below, every sub-local-system of `V₁₁ ⊗ C` containing `F³` is
everything, so `E = V₁₁` and **`V₁₁` has no level-one sub-variation at all**.
*Proof.* `V''` is stable under the transvection `T_s(x) = x + c⟨x,δ_s⟩δ_s`, so either `δ_s ∈ V''` or
`δ_s ⊥ V''`. The period `∫_{δ_s} ω_i(λ) = c'·X_i(p_s)·(λ − s) + O((λ−s)^{3/2})` with a universal `c' ≠ 0`
(integral of `dw₁dw₂dw₃/∂₄Q` over the vanishing sphere in a Morse chart), and `X(p_s) ≠ 0` at every nodal
critical point (the only critical points on `{X = 0}` are the six Klein-fixed degenerate ones). Hence
`⟨δ_s, ω_i⟩ ≢ 0`, so `δ_s ∈ V''` for every nodal `s`. For the level-one statement decompose
`V₁₁ = V' ⊕ V'^⊥` (polarized, semisimple); `F³V' = 0` forces `ω_i ∈ V'^⊥`, hence `E ⊂ V'^⊥` and
`V' ⊂ E^⊥`. ∎

**Singular-fibre bookkeeping (computed).** The critical scheme `{∂₀F = ∂₁F = ∂₂F = ∂₅F = F = 0}` of the pencil
has `vdim = 3750` in a random affine chart (Bezout `6·5⁴`, no points at infinity), its `Y2 = 0` part has `vdim = 750`,
and at the Klein-fixed point `(0:0:0:1:ζ:0)` the local multiplicities are `m = 15` (confirming the Milnor-number
prediction from the principal part `X0⁴ + Y2⁶`) and `m' = 3` for the `Y2 = 0` system (Singular, local ring `ds`,
`p = 1000033`, jobs `construct-rcr-critpts3/4`). With free Klein orbits of 4 nodes and `r`-fixed pairs:
`4N_free + 2N_r + 6·15 = 3750`, `2N_r + 6·3 = 750` ⇒ **`N_r = 366`, `N_free = 732`** nodal singular fibres. The Leray
sector dimensions of `H¹(P¹, R³f_*Q)` for `X̃` are `744 / 756 / 1132 / 580` for `(−,−)/(+,−)/(−,+)/(+,+)` (the 106
classes of `Σ₀` split 63/43 under `r`), so the naive Euler-characteristic count `dim H¹ = Σ drops − 2·rank` requires
degenerate-fibre drops `18 / 24 / 34 / 12` (total `88 = 6·15 − 2`); `18 = 6·3` in the residual sector is consistent
with six conjugate fibres, while `34` is not divisible by 6, so the naive formula misses a small `R⁴/R²` jump at the
degenerate fibres in the `(−,+)` sector — irrelevant for the proposition, which uses only the 1098 nodal fibres.

**Consequences for supports.** Write the Gysin map of a resolved divisor `Ẽ → X` through the Leray
filtrations of `Ẽ_U → U` and `X_U → U`. The graded piece `H¹(U, R¹f_*Q_Ẽ)(−1) → H¹(U, V₁₁)` has image in
`H¹(U, V'_E)` for the level-one sub-variation `V'_E = im(R¹f_*Q_Ẽ(−1) → R³f_*Q)`, which is zero by the
proposition and the rank computation: **Abel–Jacobi-type contributions (irregular surfaces in the fibres `X_λ`,
normal functions of families of 1-cycles) cannot support any residual class.** What remains is the secondary map on
monodromy-invariant fibre classes, `H⁰(U, R²f_*Q_Ẽ)^{(−,−)}(−1) → H¹(U, V₁₁)` (the Gysin map restricted to
classes whose fibrewise degree vanishes), i.e. the `(2,0)`-classes of `Ẽ` restricted to the fibres `E_λ`.
So a supporting divisor must be *horizontal* with respect to every such pencil and carry `(2,0)`-forms of
character `(−,−)` — exactly the shape of the twisted-slice mechanism M5 (`E → S` with hyperelliptic fibres
transverse to the pencils, `H^{2,0}(Ẽ)^− ⊂ H⁰(S, Ω¹_S ⊗ λ_hyp)`), and the opposite of the
"Fano-quotient" shape, where the fibrewise intermediate Jacobians `J(V_λ)` (level one, rank 82) do the
work for the companion sector.

(An earlier draft of this note claimed that GHC forces the monodromy of `V₁₁` to be reducible; that
argument ignored the secondary Leray map and is withdrawn. The proposition above is the correct statement.
The computation `gm_rank.py` returned rank 186, closing the loophole `E^⊥ ≠ 0`.)

## 5. Where each shared Lean module fails here

| module | hypothesis | failure |
|---|---|---|
| `HoldoutCover` / `GeneralizedCover` | Fermat-type cover with exchange `y ↔ T` and deck character killing the correction | no exchange on any coordinate pencil (§2.1) |
| `PolynomialCurves` (+ companion's quartic machinery) | polynomial curve family on a uniruled quotient to which the target descends | target descends to no uniruled quotient (§2.2–2.3); twisted version has `R¹p_*L ≠ 0` (§3) |
| `EllipticCompanions` (split theorem) | `F = P(U) + Q(V)` block-sign split (Hessian symmetrizer) | `hessianSymmetrizerDimension = 1` (not split) |
| `TransitiveCubic` / cyclic Jacobian transport | three paired signs ⇒ plane cubic in `X_i²` | only one sign pair (`ρ`), no `μ₃`/`μ₆` |
| `HigherResidues` | split ternary, `q = 2` | not split, `q = 1` |
| `ResidueTransport`, `AveragingProjectors`, `CertificateOperations` | generic (transport, averaging, addition) | usable once a serving correspondence exists; none found |

## 6. Experiments (job names; see log)

* **E (local, done)**: Galois group / irreducibility checks of §1; `p_g` counts of §2.2.
* **D (local, done)**: character-wise Hodge numbers of the pencil fibre (§4 table).
* **B (job `construct-rcr-all3`, running at session end)**: all reflections of `X` and of generic fibres of the six
  coordinate pencils over `F_p`, `p = 2³¹ − 1`, via the 0-dimensional system
  `{λ_i = ∂_iΦ(v); D_vΦ, D_v³Φ vanish on ker(Σλ_iz_i)}` in three charts (solved with msolve; decoder
  verifies each solution exactly and maps eigenvectors back). Control: Fermat in 5 variables must give
  `5 + 10·6 = 65` reflections (local toy in 4 variables gave the expected 40 genuine solutions).
* **C (job `construct-rcr-all`, done for base X)**: formal continuation of the strict Fermat certificate
  `D₀ = X0⁶+X1⁶+X2⁶`, `N₀ = A·X_i/6` along `F_c`, degree 6, Klein `(+,+)` sector, orders 1–5 mod `2³¹−1`
  (also the `Y`-block base `D₀ = Y0⁶+Y1⁶+Y2⁶`).
* **G (job `construct-rcr-gm2`, DONE: rank 186)**: `gm_rank.py` — Griffiths–Dwork over `F_p[ε]/ε^71`, `p = 2²⁶ − 5`, rank of the
  `∇`-span of the `(3,0)`-forms in `V₁₁` (186 would close the `E^⊥` loophole of §4). Setup verified at this second
  prime and a second random `λ₀`: `dim R_d^{even,even} = 3, 90, 90, 3, 0` in degrees `1, 7, 13, 19, 25` (the degree-25
  Jacobian decomposition, a `6825 × 14718` system, took 29 min on the loaded box); a first run crashed in the reduction
  bookkeeping and was resubmitted.
* **F**: `critpts.py` — the 3750 critical points of `Y0/Y1` on `X`: msolve segfaults on the degree-3750 system (rc −11 after 52 s, also with a 120 GB cap); replaced by the Singular script `critpts3.sing` (job `construct-rcr-critpts3`): `vdim` of the critical scheme in a random chart, of its `Y2 = 0` part, and the local multiplicities `m`, `m'` at the Klein-fixed point (local ring, `ds` ordering) — exactly the unknowns of the bookkeeping in §4.

## 6b. Cross-references to the sibling lenses (same lane, same session)
* exclusion lens (`wave2-reflection-coupled-residual-excl.md`): rigorous `deg D ≥ 7` for every denominator (E1 in
  all characters at `d ≤ 6`; Hessian relaxation zero-dimensional at `d = 5, 6`), degree 7 excluded in three of four
  pure Klein characters, `K8` (the unique non-Jacobian degree-8 E1 class) certified **not** a denominator (fixed-D
  system inconsistent, full and symmetric solvers). Any construction here must therefore produce a denominator of
  degree ≥ 8 (very likely ≥ 9), consistent with the geometric mechanisms (the companion exports no degree; the solved
  one-reflection instance needed degree 60).
* search lens (`wave2-reflection-coupled-residual-search.md`): E1-kernels at `d = 9, 10` are dimension-forced
  (`21/0/21/0`, `186/126/136/108` by character); Koszul-reduced full-E2 variable projection finds no zero at `d = 7, 8`.

## 6c. State of the remote work at session end (02:40 local / 06:40 UTC)
* `construct-rcr-all3` (6 cpus) is running `refl_all.sh` in λ-eliminated form: `fib_Y0Y1` (112 equations of degree 26 in
  4 unknowns; msolve 22 min so far), then `fib_Y0Y2, fib_Y1Y2, fib_X0X1, fib_X0X2, fib_X1X2, full_X`, each followed by
  `decode.py` (exact verification of every solution, eigenvectors mapped back to the original coordinates), then
  `fermat_cont.py` with the `Y`-block base. Outputs: `construct/<case>.json` (`decode.genuine` = number of genuine
  reflections; the inherited ones are `r|` and `ρ|` for the `Y0/Y1` pencil, `r|` for the others), `fermatcont_Y.json`.
* To collect: `cd /tmp/hodge/results/reflection-coupled-residual/construct && python3 finalize_pull.py` (runs
  `rjob pull` and appends one log line per finished case to `wave2-construct.log.jsonl`). Expected outcome: only the
  inherited reflections (a hidden reflection would be a new Fano quotient of a fibre and would reopen M1).
* Finished on the remote: `fermatcont_X.json` (order-1 inconsistent), `critpts3.out`/`critpts4.out` (counts and local
  multiplicities), `gmrank_sep2.json` (rank 186).

## 7. Open questions recorded for the lane
1. (done) `∇`-span of the `(3,0)`-forms = `V₁₁` (rank 186); critical scheme: 3750 = 4·732 + 2·366 + 6·15.
2. (withdrawn) irregular surfaces in the fibres cannot help (§4).
3. A slice `S` of the conic family on `W` with `H⁰(S, Ω¹_S ⊗ λ) ≠ 0` and nonzero Gysin image would give
   generic existence on the one-reflection locus (§3), not yet a certificate at this `F`.

## 8. Ten-line summary (kept current)
1. No witness; no mechanism in the book or the shared Lean library reaches the residual `(r,ρ) = (−,−)` sector of this `F`; each fails at a named hypothesis (ledger §0).
2. Coordinate pencils: fibres of `Y0 = λY1` are separately even in `T` and `Y2` (inherited `ρ|, r|`), quotient Fano in `P(1,1,1,2,2)`; no new exchange on any of the six coordinate pencils (§2.1).
3. The Galois closure of `X → P⁴_x` is the hyperoctahedral `B₃ × Z/2`; the target is the reflection representation `V ⊗ sgn_ρ`; every intermediate quotient reaching it (`X/τ`, a `(4,6)`-CI mod involution, a degree-8 hypersurface mod involution) is non-uniruled (§2.2).
4. Klein quotients: only `W` and `X/⟨r,ρ⟩` are Fano (wrong character); the surface fibrations of `X/τ` have `p_g = 3` fibres (§2.3); q-raising needs an `r`-odd deformation that leaves the one-reflection locus (§2.4); only one polynomial Piola map exists (M8).
5. Twisted slices of the companion's curve family are the only mechanism on this variety that can produce `(−,−)` classes: `H^{2,0}(Ẽ)^− ↪ H⁰(S, Ω¹_S ⊗ λ_hyp)` with hyperelliptic lifts of genus `e − 1`; termwise (`PolynomialCurves`) integration is obstructed by `H³(S, R¹p_*L)`, which is exactly the period obstruction of discussion #154 (§3).
6. New description: the residual sector is `H¹(P¹, j_*V₁₁)` for the `(−,−)` variation of the Klein pencil, rank 186, Hodge numbers `(3,90,90,3)`; companion sector = `H³` of the Fano fibres (rank 82, level one) (§4).
7. Proposition + computation: every sub-variation containing the `(3,0)`-forms contains all nodal `(−,−)` vanishing classes (period `∼ X_i(p_s)(λ−s)`), and the `∇`-span of the `(3,0)`-forms has full rank 186 (Griffiths–Dwork mod `p`), so `V₁₁` has **no level-one sub-variation**: Abel–Jacobi/irregular-surface supports are impossible; supports must be horizontal divisors with `(−,−)` `(2,0)`-forms (§4).
8. Remote experiments (scripts in `construct/`): Fermat-pencil continuation of the strict degree-6 certificate fails at first order; the pencil has 732 free-orbit and 366 `r`-fixed nodal fibres plus six Milnor-15 degenerate fibres (Singular); Gauss–Manin rank 186 (done); hidden-reflection search on `X` and the six pencils' generic fibres (0-dim msolve system with verified decoder) running in a λ-eliminated form after two slower formulations stalled on the loaded box.
9. Sibling lenses: `deg D ≥ 7` rigorous, `K8` excluded; any construction must yield degree ≥ 8.
10. Recommended next step for the lane: a Hodge-theoretic existence argument for M5 (twisted slice with nonzero Gysin image) at the generic point of the one-reflection locus, then an explicit transport at this `F`; and the two queued computations.
