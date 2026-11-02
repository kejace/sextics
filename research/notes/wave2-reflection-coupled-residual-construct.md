# Wave 2 — lane `reflection-coupled-residual` — lens CONSTRUCTION

STATUS: no witness (in progress). New structural results on *where* a certificate for the residual
(r,ρ) = (−,−) sector can come from; every published mechanism is shown to fail at a precise hypothesis;
two remote experiments (hidden-reflection search, Fermat-pencil continuation) are queued on the saturated
remote (load > 1200 at 01:10 remote time).

Log (one JSON line per experiment): `/tmp/hodge/results/reflection-coupled-residual/wave2-construct.log.jsonl`
Scripts/outputs: `/tmp/hodge/results/reflection-coupled-residual/construct/`
Conventions as in the dossier (`/tmp/hodge/notes/dossier-reflection-coupled-residual.md`): Klein group
`{1, r, ρ, τ}`, `r: Y2 ↦ −Y2`, `ρ: (Y0,Y1) ↦ −(Y0,Y1)`, `τ = rρ`; residue characters of the target
`A = X0³Y0³`: `(r, ρ, τ) = (−1, −1, +1)`; residual sector `(126, 492, 126)`; companion sector `(82, 374, 82)`.

## 0. Ledger: mechanisms versus hypotheses (status after this session)

| # | mechanism (task item) | hypothesis needed | status for this `(F, A)` |
|---|---|---|---|
| M1 | §4.6 deck+exchange cover, `g = r`, `X/g = W` Fano | a cover `Z → X` with an involution `σ` **not** lifting an automorphism of `X`; coordinate pencils: fibre separately even in two *new* linear forms | **fails for all 6 coordinate pencils** (proved by inspection of the fibre equations, §2.1); non-coordinate fibre reflections: search queued (`construct-rcr-reflsearch`, Exp B) |
| M2 | Galois closure of `X → P⁴_x` (item iii) | an intermediate quotient `Z̃/K` that is uniruled and has `(V ⊗ sgn_ρ)^K ≠ 0` | **fails**: Galois group is the full hyperoctahedral `B₃ × Z/2` [computed]; the target lives in the reflection representation `V ⊗ sgn_ρ`; the three maximal `K` reaching it give `X/τ` (κ = 0), a `(4,6)`-CI in `P(1⁶,2)` mod an involution with ≥ 9 invariant canonical forms, a degree-8 hypersurface in `P⁵` mod an involution with ≥ 9 invariant canonical forms; none uniruled (§2.2) |
| M3 | Fano / surface-fibration companions of Klein quotients (item ii) | uniruled quotient, or fibration over a surface with `p_g = q = 0` fibres | **fails**: only `W = X/r` and `X/⟨r,ρ⟩` are Fano (target not invariant); `X/ρ`, `X/τ` have `K ≡ 0`; the natural surface fibrations of `X/τ` over `P²_Y`, `P²_X` have fibres of weighted degree 6 in `P(1,1,1,2)` with `p_g = 3` (§2.3) |
| M4 | q-raising / Gauss–Manin from the companion (item iv) | a *family* `F + tC` on which the companion's certificate deforms, with `C` `r`-odd (to flip the character) | **fails**: `C` `r`-odd leaves the one-reflection locus, where nothing is supported; `C` `r`-even preserves the character (§2.4) |
| M5 | twisted slices of the companion's curve family on `W` (item v) | a slice `S` of the rational-curve family on `W` whose hyperelliptic lift carries anti-invariant 2-forms, i.e. `H⁰(S, Ω¹_S ⊗ λ_hyp) ≠ 0`, with nonzero Gysin image | **the only mechanism on this variety that can produce `(−,−)` classes**; explicit integration is obstructed by `H³(S, R¹p_*L) ≠ 0` (positive-genus lifts) — exactly the period obstruction met in discussion #154; a finite-field incidence seed cannot be used because the lifted curves have genus `e − 1 ≥ 1` for every weighted degree `e ≥ 2`, and `e = 1` (weighted lines) gives `r`-invariant classes only (§3) |
| M6 | Fermat-pencil continuation (item iv) | formal deformation of the degree-6 Fermat certificate along `F_c = Fermat + cP` | queued (`construct-rcr-fermatcont-X`, Exp C) |
| M7 | pencil `(Y0:Y1)` and the variation `V₁₁` (new) | a level-1 sub-VHS `V' ⊂ V₁₁` realised by irregular surfaces in the fibres | the residual sector **is** `H¹(P¹, j_*V₁₁)` with `V₁₁` of rank 186 and Hodge numbers `(3,90,90,3)` [computed]; GHC forces `V₁₁ = V' ⊕ V''` with `H¹(P¹, j_*V'') = 0` (§4) — a sharp, testable structural prediction and the real content of "cross the character boundary" |

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

**Proposition (Leray constraint).** Let `E ⊂ X` be a divisor and `Ẽ` a resolution. The part of the Gysin
image of `H²(Ẽ)(−1)` in the `(−,−)` sector lies in `H¹(P¹, j_*V'_E)`, where
`V'_E = im(R¹f_*Q_{Ẽ}(−1) → R³f_*Q)` is a sub-variation of `V₁₁` of level one (type `(0,a,a,0)`).
*Proof.* Restrict to an open `U ⊂ P¹` over which `X̃` and `Ẽ` are smooth; the Gysin map respects the
Leray filtrations with a shift of two, the `(−,−)` sector of `H⁴(X_U)` is `H¹(U, V₁₁)`, the only Leray
piece of `H²(Ẽ_U)` mapping there is `H¹(U, R¹f_*Q_{Ẽ})`, its image lies in `H¹(U, V'_E)`, and
`H¹(P¹, j_*V₁₁) ↪ H¹(U, V₁₁)` is compatible with the decomposition `V₁₁ = V'_E ⊕ (V'_E)^⊥`. ∎

**Corollary.** Since `h^{3,0}(V₁₁) = 3 ≠ 0`, the level-one part `V' := Σ_E V'_E` is a *proper*
sub-variation. The general Hodge conjecture (which predicts divisor support for the whole `τ`-invariant
`H⁴`) therefore predicts that the monodromy of the pencil on `V₁₁` is **reducible**, `V₁₁ = V' ⊕ V''`,
with `H¹(P¹, j_*V'') = 0` and `V''` containing the `(3,0)`-part; conversely, if the monodromy on `V₁₁`
is irreducible, no certificate exists for this instance (and GHC fails). The vanishing cycles of the
pencil in the `(−,−)` sector come from three kinds of critical points of `Y0/Y1` on `X`
(`∂₀F = ∂₁F = ∂₂F = ∂₅F = F = 0`, 3750 points with multiplicity): free Klein orbits (4 nodes per singular
fibre), `r`-fixed critical points (on `Y2 = 0`: 750 points, orbits of size 2, contributing only to the
`r`-anti-invariant sectors), and the six Klein-fixed points `(0:0:0:1:ζ:0)`, `ζ⁶ = −1`, at the degenerate
fibres `λ⁶ = −1`. A splitting `V' ⊕ V''` with `H¹(j_*V'') = 0` requires the `(−,−)` vanishing cycles to fall
into two mutually orthogonal families with exactly `2·rank V''` fibres feeding `V''`; the different
fixed-point types make this possible in principle (they are not conjugate under monodromy). This is the
precise sense in which the residual sector needs *new* algebraic cycles: level-one pieces realised by
irregular surfaces in the fibres `X_λ` (e.g. ruled surfaces over positive-genus parameter curves,
necessarily singular complete-intersection surfaces), varying over the pencil.

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
* **B (`construct-rcr-reflsearch`, queued)**: all reflections of `X` and of generic fibres of the six
  coordinate pencils over `F_p`, `p = 2³¹ − 1`, via the 0-dimensional system
  `{λ_i = ∂_iΦ(v); D_vΦ, D_v³Φ vanish on ker(Σλ_iz_i)}` in three charts (solved with msolve; decoder
  verifies each solution exactly and maps eigenvectors back). Control: Fermat in 5 variables must give
  `5 + 10·6 = 65` reflections (local toy in 4 variables gave the expected 40 genuine solutions).
* **C (`construct-rcr-fermatcont-X`, queued)**: formal continuation of the strict Fermat certificate
  `D₀ = X0⁶+X1⁶+X2⁶`, `N₀ = A·X_i/6` along `F_c`, degree 6, Klein `(+,+)` sector, orders 1–5 mod `2³¹−1`.

## 7. Open questions recorded for the lane
1. Monodromy of `V₁₁` (rank 186): reducible or not? (Decides existence by §4.) A finite computation:
   classify the 3750 critical points (0-dimensional msolve system), then compute the intersection graph of
   the `(−,−)` vanishing cycles (needs a Picard–Lefschetz/period computation — heavy).
2. Irregular surfaces in the Klein-symmetric sextic threefolds `X_λ` (ruled surfaces over positive-genus
   curves): any such family moving with `λ` is a support candidate for the residual sector.
3. A slice `S` of the conic family on `W` with `H⁰(S, Ω¹_S ⊗ λ) ≠ 0` and nonzero Gysin image would give
   generic existence on the one-reflection locus (§3), not yet a certificate at this `F`.
