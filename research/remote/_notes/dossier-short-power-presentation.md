# Research dossier: `short-power-presentation` (Eight powers · linear-section remainder)

Status: OPEN (no `proof.json`, no `Solution.lean`, no accepted proof as of 2026-11-01).
Discussion: https://github.com/hodge-computer/sextics/discussions/174 (one top-level comment, no replies).
Added to the catalogue in commit `a3c8c1b` (2026-10-21, Benjamin Diamond, "Add structural frontier instances and refine atlas browsing") together with the sibling frontier instances `quadratic-cubic-pullback`, `eight-plane-configuration`, `diagonal-permutation` (same target `X0^3Y0^3`, all still open). All four are screened by `scripts/audit-structural-frontier.py` (writes `checks.json`). **[verified]** by `git log` on `data/short-power-presentation/`.

Sources read in full: `data/short-power-presentation/{instance.rec, why.md, provenance.json, checks.json}`, `scripts/audit-structural-frontier.py` (all 189 lines), the whole thread #174 (`/tmp/hodge/discussions/short-power-presentation.md`), `formal/README.md`, `docs/01` (parity/averaging/lattice passages), `docs/02` §2.2 (Euler/Nambu), §2.7 (Fermat cover), `docs/03` §3.1, §3.4, `docs/03a` headings and §3A.2–3A.4 context, `docs/04` §4.1, §4.5–4.7, the accepted `fermat-q2`, `fermat-pair-rational-d4-strict`, `holdout` solutions (Fermat-type methods and the "Geometric scope of the Fermat cover" section), the unpublished local sweep `/tmp/hodge/ext/hodgeHoldout/sweep/` (RESULTS.md + `logs/short-power-presentation.log`), and the crashed session's own experiment directory `/tmp/hodge/results/short-power-presentation/` (`log.jsonl`, `e1kernel.jsonl`, scripts). `provenance.json` lists **no** `sameVarietyInstances`; `why.md` names no related instance; the book chapters 02/03/03a/04 do **not** mention this instance by name (grep for `short-power`, `eight power`, `Waring`, `catalecticant`, `sixfold` finds only an unrelated sixfold in 03a). `checks.json` cites a `screening-note.md` that does not exist anywhere in the repository or its history.

Tags: **[verified]** = checked exactly or mod a large prime in this session or by the audit script; **[rigorous]** = a proof, not evidence; **[modular evidence]** = single-prime inconsistency, not a proof; **[claim]** = stated by an author, not independently checked here; **[derived here]** = elementary derivation in this dossier, not checked by anyone else; **[conj]** = conjecture.

Notation. `Z=(X0,X1,X2,Y0,Y1,Y2)`, `τ(X,Y)=(X,−Y)`, `S=C[Z]`, `J=J(F)=(∂_0F,…,∂_5F)` (contains `F` by Euler), `R=S/J` the Jacobian ring, `±` = Y-parity, `e(X)=X0^6+X1^6+X2^6`, `e(Y)=Y0^6+Y1^6+Y2^6`, `c=5/7^7=5/823543`. A certificate for `(F,q,A)` is a nonzero homogeneous `D` coprime to `F`, six homogeneous `N_i` of degree `deg D+6q−5`, and `Q` with

`E := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^(q+1)·Q`.

Degrees at `q=1`, `d=deg D`: `deg N_i = d+1`, `deg E = 2d+6`, `deg F² = 12`, so `deg Q = 2d−6`. Normalized form (Chapter 1): `E1: q·N(F) ≡ A·D (mod F)`, `E2: D·div N − N(D) ≡ 0 (mod F^q)`. Necessary condition **[rigorous]**: reducing `E` mod `F` gives `D·(A·D − q·N(F)) ≡ 0`; `F` is irreducible (smooth) and coprime to `D`, hence `A·D ∈ J_{d+6}`. Call `K_d^± := {D ∈ S_d^± : A·D ∈ J}` the *E1 kernel*; `J_d^± ⊆ K_d^±` always, and only `K_d^±/J_d^±` matters (a `D ∈ J` is tested separately, see §3.4).

---

## 1. Exact target

### 1.1 The sextic

`q = 1`. The record `F` has **236 terms** = every Y-even sextic monomial (`28+90+90+28`), all with denominator `823543 = 7^7`. Sparse presentation (**[verified]** exactly by the audit script's `assert s.expand(F-expected)==0` and again by the crashed session's `exp0.py`):

```
F = e(X) + e(Y) + c·((L+M)^6 + (L−M)^6),     c = 5/7^7,
L = X0 + 2·X1 + 3·X2,   M = Y0 + 3·Y1 + 5·Y2.
```

Expanded in `(L,M)` (odd powers of `M` cancel):

```
F = e(X) + e(Y) + 2c·(L^6 + 15·L^4·M^2 + 15·L^2·M^4 + M^6)
  = [e(X) + 2c·L^6] + [e(Y) + 2c·M^6] + 30c·L^2·M^2·(L^2 + M^2).
```

So `F` is a split form `f(X)+g(Y)` (two Waring-rank-4 ternary sextics) plus a coupling `30c·L²M²(L²+M²)` that depends on `X` only through `L` and on `Y` only through `M`. Coefficient sanity **[verified]**: `X0^6: 823553/823543 = 1+2c`, `X1^6: 1+2c·2^6 = 1+640/823543`, `X2^6: 1+2c·3^6`, `Y2^6: 1+2c·5^6 = 979793/823543`.

**Eight-form presentation.** Put `κ = c^{1/6}` and `ℓ_1..ℓ_6 = X0,X1,X2,Y0,Y1,Y2`, `ℓ_7 = κ(L+M)`, `ℓ_8 = κ(L−M)`. Then `F = Σ_{j=1}^8 ℓ_j^6`, and the eight forms satisfy exactly two linear relations

```
ℓ_7 + ℓ_8 = 2κ·(ℓ_1 + 2ℓ_2 + 3ℓ_3),      ℓ_7 − ℓ_8 = 2κ·(ℓ_4 + 3ℓ_5 + 5ℓ_6).
```

Hence `X = V(F) ⊂ P^5` is the codimension-two linear section `Φ_6 ∩ Λ` of the Fermat sextic sixfold `Φ_6 = {Σ_{j=1}^8 ℓ_j^6 = 0} ⊂ P^7` by the `P^5 = Λ` cut out by the two relations (why.md). The eight power points in the dual `P^5` are the six coordinate points and `(1,2,3,±1,±3,±5)`; `checks.json`: every 4 of the 8 points are independent, and the only minimal dependent 5-sets ("circuits") are `{0,1,2,6,7}` and `{3,4,5,6,7}` (i.e. exactly the two relations above). **[verified by the audit script]**: middle catalecticant rank 8 mod 101, quartic power rank 8, 13 apolar quadrics annihilate `F` (so the apolar quadratic base locus is the 8 points and the Waring decomposition of rank 8 is unique: `why.md`, `checks.json`). The Hessian symmetrizer space (`{C : Hess(F)·C symmetric}`) is one-dimensional mod 101 → **no nontrivial linear direct-sum presentation** of `F` in any coordinates.

### 1.2 Target, degrees, parity

```
A = X0^3·Y0^3,  deg A = 6 = 6q,  Y-degree 3 (odd);  F is Y-even.
```

`A` has all exponents `≤ 4`, so it is nonzero in `F_5[Z]/(Z_i^5)`, which is the Jacobian ring of the mod-5 reduction; by flatness of the graded Jacobian quotient under good reduction, `A ∉ J(F)` over `Q` and `res(A·Ω/F²)` is a **nonzero primitive `τ`-invariant (3,1) class** (`checks.json: primitive_target`; independently **[verified]** by the crashed session's `D=1` run: E1 fails, rank 18 = `dim J_6^-`). Hodge numbers: `X` is a Calabi–Yau fourfold (`K_X = O`); invariant primitive Hodge numbers `(h^{3,1},h^{2,2},h^{1,3})^τ = (208, 866, 208)` (§4.1), `dim R_6^- = 208`, `dim J_6^- = 18`, `dim S_6^- = 226`; the `(4,0)` class `Ω/F` is `τ`-anti-invariant.

### 1.3 Smoothness **[verified, rigorous]**

Modulo 5 the coupling coefficient `5/7^7` vanishes, `F ≡ Σ Z_i^6`, the six partials are `Z_i^5`, the Jacobian ideal is zero-dimensional, so the reduction is smooth and smoothness lifts to characteristic zero (`checks.json`: "Good reduction to Fermat modulo 5"; asserted in the audit script).

### 1.4 Symmetry **[verified by the audit script; argument in checks.json]**

* Coordinate reflections: none. Paired sign symmetries: none (`coordinate_reflections: []`, `coordinate_paired_signs: []`).
* Projective linear automorphism group: **`{1, τ}`** (`checks.json: linear_automorphisms = "Only global Y sign, projectively"`). Argument: any linear automorphism permutes the 8 apolar base points and must preserve the two circuits; within each coordinate block the coefficient vectors `(1,2,3)` and `(1,3,5)` have distinct entries, so no block permutation survives; the only surviving element is `ℓ_7 ↔ ℓ_8` with `ℓ_4,ℓ_5,ℓ_6 ↦ −ℓ_j`, which is `τ`. (`why.md`: "their incidence and coefficients leave only the global involution".)
* Diagonal torus / characters: for a full-support Y-even sextic the exponent-difference lattice `L(F)` is `{sum 0, even Y-sum}` (Chapter 1, line ~286), so the only diagonal scalings preserving `F` up to scalar are `±1` on `Y`. **Residue character**: `τ` acts on `Ω = ι_E(dZ_0∧…∧dZ_5)` by `det = (−1)^3 = −1` and on `A` by `−1`, so `A·Ω/F²` is `τ`-invariant, as required; there is **no finer character decomposition** of the odd degree-6 targets (Chapter 1: for full support all Y-odd targets of one degree lie in a single character coset).
* Non-diagonal symmetries: in the `ℓ`-coordinates `τ` is the permutation `(7 8)` combined with signs; the ambient Fermat group `μ_6^8 ⋊ S_8` of `Φ_6` does **not** preserve `Λ` beyond `⟨τ⟩` (same argument). Nobody in the thread claims any further symmetry; benediamond's model uses exactly `τ`.
* Age screening (§3.3 criterion): `τ` has fixed locus two plane sextic curves (`X=0` and `Y=0`), transverse type `½(1,1,1)`, age `3/2 ≥ 1`; `X/τ` has `κ = 0` and is not uniruled (Kollár–Larsen, as in §4.7). No reflection quotient, no Fano quotient.

### 1.5 Structural facts beyond the record (derived here unless marked)

1. **Six Fermat-fourfold hyperplane sections** **[verified]**. With `U = L+M`, `V = L−M`: `F = e(X)+e(Y) + c(U^6+V^6)` and `U^6+V^6 = ∏_{ζ^6=−1}(V−ζU)`. On the hyperplane `H_ζ = {V = ζU}` (equivalently `M = t_ζ L`, `t_ζ = (1−ζ)/(1+ζ) = −i·tan(θ/2)` for `ζ = e^{iθ}`, so `t ∈ {±i, ±i(2−√3), ±i(2+√3)}`), `F|_{H_ζ} = e(X)+e(Y)`. Hence `X ∩ H_ζ = Φ_4 ∩ H_ζ` is a hyperplane section of the Fermat fourfold `Φ_4 = {e(X)+e(Y)=0} ⊂ P^5` (a sum of six sixth powers in five variables), and `X ∩ V(e(X)+e(Y)) = X ∩ V((L+M)^6+(L−M)^6) = ⋃_ζ X ∩ H_ζ`. Equivalently `B(1,t) := 2c(1+15t²+15t⁴+t⁶) = 2c(t²+1)(t⁴+14t²+1)`.
2. **Pencil through a Fermat surface.** The pencil `{M = tL}` has base locus `X ∩ {L=M=0}`, a sextic surface in `P^3` which is a sum of six sixth powers of linear forms (Fermat-type, `p_g = 10` if smooth — smoothness not checked); every member `X_t = X ∩ {M=tL}` is a threefold in `P^4` that is a sum of **seven** sixth powers of linear forms (`e(X)+e(Y)+B(1,t)L^6` with `Y0` eliminated), i.e. a hyperplane section of the Fermat fivefold.
3. **Affine Thom–Sebastiani structure on `L ≠ 0`** (this is the content of benediamond's chart `L=1`, restated). With `X0 = 1−2X1−3X2` and `L=1`:
   `X ∖ {L=0} = { e_X(x) + g_1(Y) = 0 } ⊂ A^2_x × A^3_Y`, `e_X(x) = (1−2X1−3X2)^6 + X1^6 + X2^6`, `g_1(Y) = e(Y) + 2c(1+15M²+15M⁴+M⁶)`.
   Both summands are tame polynomials with smooth top-degree forms (global Milnor numbers `5² = 25` and `5³ = 125`); `τ` acts on the `Y`-block only (`g_1` is even in `M`). Homogenizing the second block with a new variable `λ`: `g(Y,λ) = e(Y) + c((λ+M)^6 + (λ−M)^6)` is a sum of **five** sixth powers in `(Y0,Y1,Y2,λ)`, so `S := {g = 0} ⊂ P^3` is a hyperplane section of the Fermat threefold `Φ_3`, and
   `X = Ĵ ∩ {λ = L}`, `Ĵ := {e(X) + g(Y,λ) = 0} ⊂ P^6` = the join of the Fermat curve `{e(X)=0} ⊂ P^2` and the surface `S`.
   So `X` is a hyperplane section of a split (join) fivefold. The hyperplane section at infinity `X ∩ {L=0}` is itself split: `β(X1,X2) + e(Y) + 2cM^6` with `β` the binary sextic `(2X1+3X2)^6+X1^6+X2^6`.
4. **Shioda–Katsura dominance** **[derived here]**. The SK map `ψ: Φ_3 × Φ_3 ⇢ Φ_6`, `ψ(x,y) = (x_0y_4 : x_1y_4 : x_2y_4 : x_3y_4 : εx_4y_0 : εx_4y_1 : εx_4y_2 : εx_4y_3)`, `ε^6 = −1`, is generically finite of degree 6 (diagonal `μ_6`). The two linear equations of `Λ` pull back to two **bidegree-(1,1)** equations, e.g. `εx_4(y_2+y_3) = 2κ·y_4(x_0+2x_1+3x_2)` and `εx_4(y_2−y_3) = 2κ(x_3y_4 + 3εx_4y_0 + 5εx_4y_1)`. Hence `X' := Φ_3 × Φ_3 ∩ {two (1,1)-divisors}` is a fourfold dominating `X` with degree 6, and the projection `X' → Φ_3` (first factor) has as fibre over `x` the plane section `Φ_3 ∩ Λ_x` (a plane sextic of genus 10 for smooth fibres): **`X` is dominated by a family of plane sections of the Fermat threefold parametrized by the Fermat threefold**. The same works with `Φ_2 × Φ_4` and `Φ_1 × Φ_5`. This is the precise meaning of why.md's "vanishing cycles of those sections and correspondences arising from diagonal covers". `X'` is of general type (`K_{X'} = O(3,3)|`), with `h^{4,0} > 0` in nontrivial `μ_6`-characters.
5. **Lines** **[verified numerically, 1e-9 rank tolerance]**. Each of the `105·6^4 = 136,080` Fermat `P^3`s `{ℓ_a = ζ_1ℓ_b, ℓ_c = ζ_2ℓ_d, ℓ_e = ζ_3ℓ_f, ℓ_g = ζ_4ℓ_h}` (`ζ^6 = −1`) meets `Λ` in a subspace of projective dimension **exactly 1**, giving **136,066 distinct lines on `X`**; none of these `P^3`s meets `Λ` in a plane, so this construction yields **no planes** on `X` (contrast `eight-plane-configuration`). The sextic fourfold's Fano scheme of lines has expected dimension 1; whether these lines lie on positive-dimensional families was not checked.
6. **Fermat-join divisors.** For `I ⊂ {1..8}`, `D_I := Σ_{j∈I} ℓ_j^6` cuts `X` in the linear section of the join `J(Φ_{|I|−2}, Φ_{6−|I|}) ⊂ Φ_6` (a complete intersection of the two cones); on `X`, `D_I ≡ −D_{I^c}`. For `|I| = 4` with `I ⊂ {1..6}`, `X ∩ V(D_I)` is a cone with vertex a line over a Fermat sextic surface (`p_g = 10`), which satisfies the necessary (3,1)-support condition of §3.1 (non-normal along the vertex) — yet all such `D_I` fail (§3.4), so the geometric criterion is far from sufficient here.

---

## 2. Companion / related solved instances and why their methods fail here

There is no companion instance on the same variety. Related by structure:

| Related | Method and denominator | Why it fails here |
|---|---|---|
| `fermat-q2` (accepted; `F` = Fermat, `q=2`, `A=X0^3X1^2X2^2Y0^3Y1^2`) | Subset Euler field `v = A·E_I/(q·E_I(F))`, `D = E_I(F) = 6(X1^6+X2^6)` (degree 6, strict), valid when `E_I²(F) = 6E_I(F)` and `deg_I A + |I| = 6` (§2.2). | `E_I²(F) = 6E_I(F)` fails for every proper subset `I` because `(L±M)^6` are not homogeneous in any proper coordinate subset; and **every** Y-even degree-6 `D` is excluded rigorously (§3.4, E1 kernel), which covers all `E_I(F)`-type and sums-of-powers denominators. |
| `fermat-pair-rational-d4-strict` (accepted; Fermat, `q=1`, `A = X0^3X1X2Y0`) | Quartic cofactor `D = x^4−x²y²+y^4` from `x^6+y^6 = (x²+y²)D`, two numerator components, strict. Fermat case in Diamond, arXiv:2603.29157. | Needs the target to be divisible by a coordinate pair whose sixth powers appear as a pure binary block; here no binary block exists (Hessian symmetrizer scalar) and all `D` of degree `≤ 6` of either parity are excluded (§3.4). |
| Diamond arXiv:2603.29157, Thm 3.16 (`F = f(X) − f(Y)`, `f` of Waring rank 3) and the book's "arbitrary `f(X)+g(Y)`" row of §4.1 (finite elliptic-companion construction, every target, every `q`, Chapter 3a) | Split equation; `X/τ ~ (T_f × T_g)/diag μ_3` with `T_f: u³ = f`, pure-cube relation `r³ = −f/g`; Fano elliptic companions `{f(X)+s³+t³=0} ⊂ P(1,1,1,2,2)` and contraction of the elliptic factors (§3A.2–3A.3). | `F` is not split in any linear coordinates (Hessian symmetrizer). Structurally (benediamond + §4.7): the `τ`-quotient relation between the two pencil parameters is the **general cubic** `(s+2c)n³+30cn²+30cn+(e+2c)=0` with `S_3` monodromy, not the pure cube `n³ = −e/s` with `μ_3`; the diagonal-`μ_3` identification that makes the elliptic contraction work is destroyed by the coupling `30cL²M²(L²+M²)`. |
| `holdout` family (accepted, PRs #128, #135–#137, #140, #151, #161; chain sextic `Fermat + X0^3X1Y0Y1 + X1^3X2Y1Y2`; fixed-variety theorem: every odd target, every `q`) | Fermat cover `T^6 = Y0^6+Y2^6` creating the exchange `y ↔ T`, weighted quotient in `P(1,1,1,2,2)` with Fano fibres, moving conic, polynomial integration, deck-character target selection (§2.7, §4.6). | Requires couplings in which a pair of coordinates enters only through `y·T` so that a reflection appears on a cover; the coupling here is a binary form in two *linear forms* involving all six coordinates; no cover with a reflection quotient has been proposed (my assessment; nobody has claimed applicability). |
| Nonlinear polynomial Piola maps (§3.4; e.g. cubic in six Y-even quadrics) | `F = G(Φ)`, `7 ≤ Σ deg Φ_i ≤ 12`, serves `δ·b(Φ)`. | `F` is not coordinatewise even (odd `X_iX_j`, `Y_iY_j` cross terms from `(L±M)²`), and a presentation through five linear forms plus a quadric is impossible since `X0,X1,X2,L,M` span only four of the six directions and `e(Y)` needs the other two. **[derived here]** No polynomial `Φ` with `M ≥ 7` is known; none was found by the audit (its `quadratic-cubic-pullback` density check is specific to that sibling). |
| Polarized split `f(X+Y)+f(X−Y)` (§3A.6, K3 correspondence; the related `polarized-k3-remainder` is itself open) | Order-three plane family, Kuga–Satake criterion. | Only the coupling part `c((L+M)^6+(L−M)^6)` is polarized, and in the pair `(L,M)`, not in `(X,Y)`; the Fermat part breaks it. |
| `dense-rational`, `eight-plane-configuration` (open; same target, same "no symmetry" class) | — | No transferable result; `dense-rational`'s rigorous `deg D ≥ 7` (fqlx) and the extension to `≥ 10` are the model for §3.4 here. |

`why.md`'s own verdict: "The ambient geometry does not itself supply the required support… restricting ambient cohomology does not exhaust [the primitive middle cohomology]" (indeed `H^4(Φ_6)` is one-dimensional, so all 208 invariant (3,1) classes are vanishing cohomology of the linear section). `provenance.json: existenceStatus = "No certificate or existence conclusion for this target is supplied by the screened constructions."`

---

## 3. Prior work

### 3.1 The discussion thread (#174) — complete record

**Comment DC_kwDOVLOOJc4BJDb9, `benediamond`, 2026-10-23T17:10:38Z (never edited, no replies). Title: "A product-of-curves model, with the actual quotient density."** Links: only the permalink to `instance.rec` at commit `faa9df2`. No reproduction packet, no code, no logs.

Content, verbatim in substance:

* Chart `L = 1`, `m := M`, for `m ≠ 0` put `y_i = Y_i/m`; then `X0+2X1+3X2 = 1`, `y0+3y1+5y2 = 1`. With `e = Σ X_i^6`, `s = Σ y_i^6` the equation is **exactly** `e + m^6 s + 2c(1+15m²+15m⁴+m⁶) = 0` **[claim: verified exactly by the author; consistent with §1.5(3)]**.
* The geometric fibre over `(e,m)` is an open part of `C_e × D_s`, `C_e: Σ X_i^6 = e(X0+2X1+3X2)^6 ⊂ P²`, `D_s: Σ y_i^6 = s(y0+3y1+5y2)^6 ⊂ P²`, `s = (−e − 2c(1+15m²+15m⁴+m⁶))/m^6`. Both pencils contain the smooth Fermat member at parameter 0, so generic members are smooth plane sextics of genus 10; the generic compact surface fibre has `p_g = 100`, `q = 20`, and its `H²` retains `H¹(C_e) ⊗ H¹(D_s)`.
* `τ` sends `m ↦ −m` and fixes the normalized `X_i, y_i`; the quotient has base coordinate `n = m²` and the **same** geometric surface fibre; the map from the base to the two pencil parameters is the **cubic relation** `(s+2c)n³ + 30c n² + 30c n + (e+2c) = 0`. "The sign quotient does not project away the tensor-product cohomology."
* Source without hidden Jacobian: coordinates `(x0,x1,y0,y1,m)`, `X = (x0,x1,(1−x0−2x1)/3)`, `Y = m(y0,y1,(1−y0−3y1)/5)`; for the standard `Ω` in the order `(X0,X1,X2,Y0,Y1,Y2)`, restriction to `L=1` gives `Ω = (m²/15)·dx0∧dx1∧dy0∧dy1∧dm`; with `A = m³x0³y0³` and `dn = 2m dm`,
  `A·Ω/F² = (n² x0³ y0³/30) · dx0∧dx1∧dy0∧dy1∧dn / F_q²`, `F_q = e(x0,x1) + n³ s(y0,y1) + 2c(1+15n+15n²+n³)`.
  "The source carries `n²`, not a constant coefficient on the parameter cover." **[claim: constants `1/3` and `m²/5` verified exactly by the author]**.
* Explicit disclaimers: "does not show that the specified absolute residue has a nonzero component in every possible fiber summand"; "the positive `p_g` of this fibration does not exclude other fibrations or selected correspondences"; "No solve, homogeneous certificate, or Lean proof is claimed. ChatGPT-assisted investigation for `benediamond`, with separate AI-agent checks; no independent human review is claimed."

What this comment establishes: a compact birational model `X/τ ≈ {(x,y,n): (ẽ_Y(y)+2c)n³+30cn²+30cn+(e_X(x)+2c)=0}` (my notation), the identification of the obstruction to the split-case machinery (the cubic is not a pure cube), and the exact transported source. It proposes "a correspondence or Gauss–Manin approach to the coupled tensor-product system" as the open sub-problem: find an actual support or rational primitive for the displayed source with its scale density `n²` and the cubic relation retained. **No computational exclusion, denominator family, degree sweep, or prime appears anywhere in the thread.**

Additional elementary consequence **[derived here]**: writing the cubic as `n³ + (30c/σ)n² + (30c/σ)n + (e+2c)/σ`, `σ = s+2c`, its roots satisfy `e_1 = −e_2`, so `n_3 = −(n_1+n_2+n_1n_2)/(1+n_1+n_2)`; the `S_3`-Galois closure of the base cover `S = {(e,s,n)} → A²_{(e,s)}` is the **rational surface `Ŝ ≅ A²_{(n_1,n_2)}`** with `σ = −30c/(n_1+n_2+n_3)`, `e = −σn_1n_2n_3 − 2c`. On the locus `e = s` the cubic factors as `(1+n)·((e+2c)(1−n+n²) + 30cn)`, the root `n = −1` (`m = ±i`) being the two Fermat-fourfold sections `H_{±i}` of §1.5(1). The Galois closure `Ẑ = C_tot ×_{A¹_e} Ŝ ×_{A¹_s} D_tot` of `X/τ → A⁴_{(x,y)}` is the §4.7 object; by §4.7 neither `Ẑ` nor `Ẑ/C_2 = X/τ` is uniruled and `Ẑ/A_3`, `Ẑ/S_3` carry no `V_std`-invariants, so **the natural `S_3` cover cannot serve the standard summand through the averaging ideal of §3.1**, for this `F` as for every even sextic. The special feature here — the whole cover is pulled back from a rational *surface* — has not been exploited by anyone.

### 3.2 Unpublished local sweep (DiamondJM-style; `/tmp/hodge/ext/hodgeHoldout/sweep`, commits by "Claude" 2026-10-23; not posted to thread #174) — **[modular evidence]**

`fixed_d.py` (prime `p = 268435399`, full numerator space restricted only by the Y-parity forced by a Y-even `D`) tested the full cleared system `A D² − q D N(F) + F(D div N − N(D)) = F² Q` for **106** natural degree-6 Y-even denominators:
* `pow{I}`: `Σ_{k∈I} Z_k^6` for all subsets `I` of the six coordinates with `|I| ∈ {2,3,4}` (15+20+15 = 50);
* `Fsub{I}`: the sub-equation of `F` in the coordinates `I` (other coordinates set to 0), `|I| ∈ {2,3,4,5}` (15+20+15+6 = 56).
Result: **0 solvable**; each run `rows = 16797, cols = 2602`, rank 2304–2352, 41–58 s. Controls (RESULTS.md): on the Fermat sextic `D = X0^6+X1^6+X2^6` serves `X0^3Y0^3` and `X1^6+X2^6+Y1^6` does not. RESULTS.md's own caveat: "On the solved one-reflection the accepted support has degree 30, so these negatives are weak evidence." All 106 candidates are Y-even of degree 6, so they are in fact **rigorously** excluded by §3.4 below.

### 3.3 The crashed session's experiments (`/tmp/hodge/results/short-power-presentation/log.jsonl`, engine v2 after the 17:25 sign-bug fix; entries before that time involve only `D` with `dD = 0` or the E1 test, which the bug did not touch)

| `D` | degree | E1 (`A·D ∈ J`) | full system / served space | time |
|---|---|---|---|---|
| `1` | 0 | fails (rank 18) | — (confirms `A ∉ J`) | 0.01 s |
| `e(X)` | 6 | fails | — | 19 s |
| `e(X) − L^6` (pencil member, cone over a genus-10 plane curve with vertex `P²_Y`) | 6 | fails (`rank_aug = rank+1`) | served space of *all* odd sextic targets = 18-dim = `J_6^-` only | 793 s |
| `E_X(F) = Σ X_i ∂_{X_i}F` (lies in `J`) | 6 | trivially OK | served space = 18-dim = `J_6^-` only → no primitive target served | 915 s |

Planned but not run before the crash (`exp2.py` lists): `e(X)+L^6`, `e(Y)−M^6`, `e(Y)`, `e(X)+M^6`, `e(Y)+L^6`, `E_Y(F)`, `e(X)+L⁴M²`, `e(X)+L²M⁴`, `e(X)−2L⁶+M⁶`, `e(X)+3L⁴M²−M⁶`. All are Y-even of degree 6 and are now excluded by §3.4 (those outside `J`) or are of the `J`-type of `E_X(F)` (not excluded: see §4).

### 3.4 E1-kernel sweep — **[rigorous exclusions]**

Rank computations mod `p = 2147483647` (`core.PRIMES[0]`); since `rank_p ≤ rank_Q`, nullity over `Q` is at most the nullity mod `p`, so a zero kernel mod `p` is a proof over `Q` (the `1/7^7` coefficients are `p`-units). Files: `e1kernel.jsonl` (even `D`, crashed session, 17:37–17:40), `e1kernel2.jsonl` (odd `D`, this session). Dimension bookkeeping from §4.1 Hilbert series `H(t)=(1+…+t⁴)^6`, `Θ(t)=(1+t²+…+t⁸)³`:

| `d` | `dim S_d^+` | `dim J_d^+` | `K_d^+ mod J` | `dim S_d^-` | `dim J_d^-` | `K_d^- mod J` | `dim R_d^±` | `dim R_{d+6}^∓` |
|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 0 | **0** | 3 | 0 | **0** | 3 | 333 |
| 2 | 12 | 0 | **0** | 9 | 0 | **0** | 12/9 | 468/483 |
| 3 | 28 | 0 | **0** | 28 | 0 | **0** | 28 | 623 |
| 4 | 66 | 0 | **0** | 60 | 0 | **0** | 66/60 | 744/762 |
| 5 | 126 | 3 | **0** | 126 | 3 | **0** | 123 | 843 |
| 6 | 236 | 18 | **0** (92 s) | 226 | 18 | **0** (266 s) | 218/208 | 866/885 |
| 7–9 | 396/651/1001 | | pending (parallel agent) | | | pending | 333,483,623 | 843,744,623 |
| 10 | 1512 | | ≥ 294 forced by dimension | 1491 | | ≥ 261 forced | 762/744 | 468/483 |

**Theorem (rigorous, this and the crashed session).** Any certificate for `(F,1,X0^3Y0^3)` whose denominator is Y-even or Y-odd of degree `≤ 6` has `D ∈ J(F)`; in particular **no Y-pure denominator of degree `≤ 6` outside `J` exists**, and (since `D·τ(D)` is Y-even of degree `2d`) **no mixed-parity denominator of degree `≤ 3` exists**. For mixed parity `D = D^+ + D^-` of degree 4–6, `A·D ∈ J` forces `D^+ ∈ J_d^+` and `D^- ∈ J_d^-` separately (`J` is `τ`-stable). **[verified]** `(L+M)^6+(L−M)^6 ∉ J_6^+` and `e(X)+e(Y) ∉ J_6^+`, so in particular the single-power denominators `(L±M)^6` and their sums with coordinate powers are excluded.

For `dense-rational` the analogous kernel is zero through `d = 9` and jumps at `d = 10` to exactly the dimension-forced values `294/261` (memory notes; fqlx proved `≥ 7`). If the pending `d = 7,8,9` runs return zero here as well, the **minimal degree of any Y-pure denominator is 10**, with the E1-admissible space at `d = 10` of dimension exactly `294` (even) and `261` (odd) modulo `J`.

What is *not* excluded by this: `D ∈ J_d` (e.g. `E_X(F)`, which serves nothing; the other 17 directions of `J_6^+` and all of `J_6^-` untested), and all `d ≥ 7`.

---

## 4. What remains unexplored

1. **E1 kernel at `d = 7,8,9` for both parities, and `d = 10` basis** — running now (parallel agent); not yet in the record. A nonzero kernel at `d < 10` would be structurally meaningful (it would single out special denominators of this `F`).
2. **Denominators inside `J`**: the full 18-dim family `J_6^+ = span{X_j∂_{X_i}F, Y_j∂_{Y_i}F}` (only `E_X(F)` tested) and `J_6^- = span{Y_j∂_{X_i}F, X_i∂_{Y_j}F}`; general `D = Σ a_i ∂_iF` with linear `a_i`. (For `D ∈ J`, `N = A·a + Koszul + F·w` solves E1 and only E2 is a condition; the served space is cheap to compute, ~15 min each at `d = 6` with the present engine.)
3. **Degree 7–12 structured denominators built from the eight forms** (never tested at any degree ≥ 7): products of the `ℓ_j` (hyperplane arrangements, `d ≤ 8`); `D_I · D_{I'}` (products of join sections, `d = 12`); `D_I · τ(D_I)` for mixed `I` (e.g. `((L+M)^6 + X0^6)((L−M)^6 + X0^6)`); binary forms in `(L,M)` times coordinate powers; the pencil members `e(X) − λL^6` (λ at the discriminant values where `C_λ` is singular) times Y-factors.
4. **Nonlinear/joint search** (bilinear in `D` and `N`) has never been run at any degree. At `d = 10` the E1-admissible family is (expected) 294-dim; a certificate, if one exists at `d = 10`, has `D` in that family.
5. **Deformation / continuation from the solved neighbour** `c → 0` (Fermat). Never attempted. At `c = 0` all targets have explicit certificates; the pencil `F_c = Fermat + c((L+M)^6+(L−M)^6)` is a one-parameter family, and by semicontinuity the degree-6 E1 kernel is zero for generic `c` (it is zero at `c = 5/7^7`), so no degree-6 Fermat certificate deforms; the question is which *non-reduced* Fermat certificates `(K·e(X), K·N^{Fermat}, K²Q)` of degree `≥ 7` deform to first (and higher) order in `c`.
6. **Geometric support**: no divisor on `X` has been shown to carry the (3,1) class. Untested geometric candidates satisfying the §3.1 necessary condition: cones with vertex a line over surfaces with `p_g > 0` (the Fermat-join sections `D_I`, `|I|=4`, satisfy it and fail algebraically); ruled threefolds swept by one-parameter families of the 136,066 lines (not known to exist); the singular members of the pencil `X_t` (critical values of `B(1,t)`-twisted Fermat; singularity types unknown).
7. **Symmetry sectors**: there are none beyond `τ`; the only "sector" idea is the `S_3`-cover of §4.7 / benediamond, now known to have a rational Galois-closure base `Ŝ` (§3.1), which nobody has used.
8. **Covers creating reflections** (holdout method): no candidate cover has been written down for a coupling that is a binary form in two linear forms.
9. **The SK model** `X' ⊂ Φ_3 × Φ_3` (§1.5(4)) and the join model `X = Ĵ ∩ H` (§1.5(3)) are new here; no one has analysed the Hodge structure of the plane-section family `Φ_3 ∩ Λ_x` over `Φ_3` or the residue sequence `H^5(Ĵ∖X) → H^4(X)_van`.
10. **Lean**: nothing has been attempted; no shared-library mechanism applies directly.

---

## 5. Candidate attack plan

Linear-system sizes at `q = 1` for a Y-even `D` of degree `d` (unknowns: `N_i` of degree `d+1` with the forced parities, `Q` of degree `2d−6` Y-odd; rows: Y-odd monomials of degree `2d+6`):

| `d` | `#N` | `#Q` | unknowns | rows | feasibility (dense `nmod_mat` rank) |
|---|---|---|---|---|---|
| 6 | 2376 | 226 | 2602 | 16797 | 40 s (sweep); 15 min with served-space nullspace |
| 7 | 3861 | 636 | 4497 | 26532 | minutes |
| 8 | 6006 | 1491 | 7497 | 40326 | ~10–30 min |
| 9 | 9009 | 3080 | 12089 | 59332 | ~1–2 h dense; sparse recommended |
| 10 | 13104 | 5796 | 18900 | 84903 | needs sparse/structured elimination (13 GB dense) |
| 12 | 25704 | 16797 | 42501 | 162248 | sparse only; 40 cores / 125 GB available |

**E1 first.** Every candidate `D` should be screened by the E1 test (`A·D ∈ J_{d+6}`, one rank computation, seconds to minutes) before the full system.

### Experiment 1 — finish the low-degree classification (cheap, decisive)
*Idea.* Make the "minimal denominator degree" statement exact.
*Construction.* (a) E1 kernels at `d = 7,8,9` (both parities; running), `d = 10` with an explicit basis of `K_10^±/J`. (b) Served spaces for a basis of `J_6^+` and `J_6^-` (36 runs of the `target_image` routine, ≈ 15 min each; or one joint run with `D = Σ a_{ij} Z_j ∂_iF` treated symbolically — not linear in `a`, so do the 36). (c) For each `D ∈ K_10^±` basis element with few terms in the `(ℓ_1..ℓ_8)`-basis, run the full `d=10` system.
*Size.* As in the table; (a) at `d = 10` is a `10168 × 14.6k` rank mod `p` (≈ 10 min).
*Chance.* Near-certain to complete; expected outcome `K_7 = K_8 = K_9 = 0`, `dim K_10^+/J = 294` (pattern of `dense-rational`), served spaces of `J_6` trivial. Value: a rigorous "`deg D ≥ 10`" for the record; low chance (<5%) of a witness.

### Experiment 2 — structured denominators from the eight forms (degrees 7–12)
*Idea.* The only non-generic structure of `F` is the eight-form configuration; natural divisors are unions of Fermat-join sections and `ℓ`-hyperplane sections.
*Construction.* Screen by E1 (seconds each): `∏_{j∈I} ℓ_j` (`|I| = 7,8`; mixed parity → use `D·τ(D)`), `D_I·D_{I'}` and `D_I·τ(D_I)` (`d = 12`), `ℓ_j · D_I` (`d = 7`), `ℓ_jℓ_k · D_I` (`d = 8`), `β(L,M)·D_I` with `β` the quadratic and quartic factors `(t²+1)`, `(t⁴+14t²+1)` of `B` (i.e. `L²+M²`, `L⁴+14L²M²+M⁴`; degrees 8, 10), `(e(X)+2cL^6)·(e(Y)+2cM^6)` (`d = 12`, the product of the two "block" sextics), `e(X)·e(Y)`. Only E1-passing candidates go to the full solve.
*Size.* E1: rank of `[J_{d+6} | A·D]`, up to `18k × 20k` at `d = 12`. Full solve only if E1 passes.
*Chance.* The E1 pass rate for structured `D` is the real test; I expect most to fail (as all degree-6 ones do). 10–15% that some product denominator passes E1; a passing `D` at `d ≤ 10` would be a strong lead. Overall chance of a witness: ~5%.

### Experiment 3 — first-order deformation from Fermat (`c → 0`)
*Idea.* `F_c = F_0 + c·F_1`, `F_0` Fermat, `F_1 = (L+M)^6+(L−M)^6`. A certificate algebraic in `c` specializes to `c = 5/7^7` (clear denominators in `c`; the specialized `D` must be nonzero and coprime). Determine which Fermat certificates of degree `d` admit a first-order deformation.
*Construction.* Take `(D_0, N_0, Q_0) = (K·e(X), K·N^F, K²·Q^F)` with `(e(X), N^F, Q^F)` the known Fermat certificate for `X0^3Y0^3` (RESULTS.md control) and `K` an even form of degree `d−6`. Order-`c` equation: `∂_D E·D_1 + ∂_N E·N_1 + [−D_0 N_0(F_1) + F_1(D_0 div N_0 − N_0(D_0))] = F_0² Q_1 + 2F_0F_1 Q_0`, linear in `(D_1,N_1,Q_1)` for fixed `K`, quadratic in `K`. Solve for the variety of unobstructed `K` (start with `K ∈ span{L⁴, L²M², M⁴, X_i²X_j², Y_i²Y_j², e(X)-type}` and random `K`), then second order. Use the Fermat `μ_6^6`-multigrading of `F_0` to block-diagonalize the linear part (the `F_1` term couples gradings in a controlled way: `F_1` has all 236 monomials).
*Size.* At `d = 10` the linear part is the `84903 × 18900` system of the table, but graded into ≈ 6^5 character blocks at Fermat — each tiny; the coupling term is dense. Realistic with sparse tools.
*Chance.* Novel and informative: it either produces a degree bound for certificates algebraic in `c` or a concrete `(D_0 + cD_1 + …)` ansatz. Formal solvability at every order does not by itself give an algebraic certificate. ~10% that it leads to a witness; high value as structural information.

### Experiment 4 — joint (bilinear) search in the `d = 10` E1 family
*Idea.* `D ∈ K_10^+` (294-dim) and the served space `T(D) ⊂ R_6^-`; the condition `A ∈ T(D)` is bilinear. Random `D` in the family will serve nothing; look instead for `D` in the family with small support in the `ℓ`-monomial basis or with a factorization.
*Construction.* Compute the basis of `K_10^+/J`; intersect with structured subspaces (`Q[ℓ_1..ℓ_8]_{10}` has `C(17,7) = 19448` monomials but the `ℓ`'s satisfy two relations — use the lattice-reduced representation); for each candidate compute `T(D)` (nullspace projection, as `target_image` does) rather than solving for `A` only, since `T(D)` tells which targets *are* served.
*Size.* `d=10` full systems, sparse; 30–60 min each with a sparse rank engine.
*Chance.* ~5%; mainly diagnostic (tells whether *any* degree-10 denominator serves *any* primitive target).

### Experiment 5 — geometric pre-screen by singularities (cheap necessary condition)
*Idea.* §3.1: a divisor supporting a nonzero (3,1) class must have a component that is non-normal or has non-rational singularities. Use this to discard candidate `D`'s before linear algebra, and to search for *good* ones.
*Construction.* For each candidate from Experiment 2 and for the pencil members `X_t`, compute `Sing(X ∩ V(D))` (Jacobian ideal of `(F, D)` mod a prime, Singular/`msolve`), its dimension, and for isolated points the singularity type (is `x^a+y^b+z^c+w^d`-type rational? `1/a+1/b+1/c+1/d > 1`). Enumerate the critical values `t` of the pencil (discriminant of `e(X)+e(Y)+B(1,t)L^6` on `P^4`) and check which `X_t` have non-isolated or non-rational singularities (the six `H_ζ` sections are smooth Fermat-fourfold sections if `H_ζ` is not tangent to `Φ_4` — check).
*Size.* Gröbner computations in 6 variables of degree ≤ 12 mod `p`: minutes each.
*Chance.* Will not produce a witness by itself; raises the hit rate of Experiments 2 and 4 and documents, for the record, which natural divisors are even geometrically admissible.

### Experiment 6 — the `S_3` cover with rational Galois base `Ŝ`
*Idea.* benediamond's cubic has `e_1 = −e_2`, so the Galois closure of the base is rational (`Ŝ ≅ A²_{(n_1,n_2)}`), and `Ẑ = T_1 ×_Ŝ T_2` with `T_1 = {(x,n_1,n_2): e_X(x) = e(n_1,n_2)}`, `T_2 = {(y,n_1,n_2): ẽ_Y(y) = s(n_1,n_2)}`. The invariant (3,1) classes of `X` sit in the `V_std`-isotypic part of `H^4(Ẑ)` and come from `H^2(Ŝ, H^1(C_e)⊗H^1(D_s))`-type pieces.
*Construction.* (a) Compute the genus of the level curves `{e(n_1,n_2) = const} ⊂ Ŝ` (explicit rational function `e + 2c = 30c·n_1n_2n_3/(n_1+n_2+n_3)` with `n_3 = −(n_1+n_2+n_1n_2)/(1+n_1+n_2)`); if genus 0, `T_1` and `T_2` are ruled over `A²_x`, `A²_y` and the restriction `H²(T_1)⊗H²(T_2) → H^4(Ẑ)` is algebraic. (b) Look for *other* correspondences on `Ẑ` (e.g. from the Fermat curve pencil automorphisms `C_e → C_e` acting fibrewise, `μ_6`-actions on `T_1`) and compute whether the resulting averaging ideal contains the `V_std ⊗ (relevant character)` projector; §4.7's obstruction only concerns the subgroup quotients of `S_3`.
*Size.* Symbolic; no large linear algebra.
*Chance.* Low (≤10%) for support; §4.7 shows the natural quotients are insufficient and `Ẑ` is not uniruled. But (a) is a one-hour calculation that settles whether this route is dead.

### Experiment 7 — the Shioda–Katsura model `X' ⊂ Φ_3 × Φ_3` and the plane-section family
*Idea.* `X' → Φ_3` is a family of plane sections `Γ_x = Φ_3 ∩ Λ_x` of the Fermat threefold. Fermat varieties have abundant algebraic cycles and known Hodge structures (Shioda's inductive structure; the Hodge conjecture is known for Fermat varieties of prime degree and degree 4, and for several composite degrees in low dimension — for degree 6 in dimension 3 I have not verified the literature status). The pullback of the target class to `X'` lies in `H^3(Φ_3, R^1)`-type Leray pieces of the curve family.
*Construction.* (a) Describe `Λ_x` explicitly (two (1,1) forms → a plane in `P^4_y` depending linearly on `x`); determine the discriminant locus in `Φ_3` and the monodromy of `H^1(Γ_x)` (rank 20). (b) Ask whether the relative Jacobian `J(Γ_x)/Φ_3` has a correspondence to a product/abelian scheme with known support (Fermat-curve Jacobians are CM; plane sections of `Φ_3` are *not* Fermat curves, but the family is special). (c) Compare with the §3A.3 contraction: the elliptic companions used there are weighted Fano hypersurfaces `{f(X)+s³+t³=0}`; the analogous "Fano companion" for the `Y`-block is `W_Y = {g(Y,λ)+s³+t³=0} ⊂ P(1,1,1,1,2,2)`, a quasi-smooth degree-6 hypersurface with `−K = O(2)` (Fano fourfold, hence uniruled and divisor-supported by §3.2), and for the `X`-block `W_X = {e(X)+s³+t³=0} ⊂ P(1,1,1,2,2)` (the split-case companion). The obstruction is that the join structure is only affine (chart `L ≠ 0`) and the contraction requires the diagonal `μ_3` identification; find whether a correspondence `X ⇢ W_X × W_Y` exists with the boundary `X ∩ {L=0}` (itself split) handled by the transfer lemma's boundary terms.
*Size.* Theoretical; small verifications in Macaulay2/Singular.
*Chance.* This is the route most faithful to why.md's intent ("vanishing cycles of those sections and correspondences arising from diagonal covers"). 10–15% that a support argument emerges; it would be a genuinely new mechanism ("diagonal methods beyond rewriting a split equation") and would need its own Lean development.

### Experiment 8 — engine upgrade enabling `d = 10–12`
*Idea.* Experiments 2–4 are blocked at `d ≥ 10` by dense rank computations.
*Construction.* Sparse modular elimination (structured Gaussian elimination + Wiedemann, or `msolve`/`LinBox`), exploiting that each column of the certificate matrix has ≤ 236·(#terms of `N_i`) nonzeros and the `τ`-parity block structure; parallelize over 40 cores; verify any candidate solution by rational reconstruction and exact check (`core.py`).
*Size.* `162k × 42.5k` at `d = 12` with ~2% density.
*Chance.* Pure engineering; needed before any degree-12 claim can be made.

---

## 6. Feasibility assessment

**Explicit witness.** Rigorously no Y-pure `D` of degree `≤ 6` outside `J` (both parities) and, pending the running `d = 7–9` checks, very probably none below degree 10. A degree-10 witness would have ~19k unknown coefficients and an identity of ~85k monomials; a degree-12 witness ~42k unknowns / 162k monomials. Nothing in the eight-form structure has yet produced a denominator passing even E1. The accepted solve on the (symmetric) `one-reflection` family has support of degree 30 and is nonconstructive. I estimate the probability that a fixed-denominator search (Experiments 1, 2, 4, 8) finds a witness at `d ≤ 12` at **≈ 10%**; the structured candidates of Experiment 2 are the only ones worth the compute. Experiment 3 (deformation from Fermat) is the most original computational idea and could at least pin the degree.

**Nonconstructive argument.** The variety has only `τ`; `X/τ` has `κ = 0`; the natural `S_3` cover is obstructed (§4.7); the split-case elliptic contraction fails because the pencil relation is a general cubic. What *is* special: (i) `X` is a hyperplane section of a join of a Fermat curve and a Fermat-threefold hyperplane section; (ii) `X` is dominated with degree 6 by a `(1,1)`-complete intersection in `Φ_3 × Φ_3`, a family of plane sections of `Φ_3` over `Φ_3`; (iii) six hyperplane sections of `X` are Fermat-fourfold sections; (iv) `X` carries 136,066 Fermat lines; (v) the Galois base of benediamond's cover is rational. None of (i)–(v) is known to give divisor support for the 208-dimensional invariant (3,1) space; each is a concrete object on which the §3.1 averaging-ideal method or the §3A contraction could be attempted (Experiments 6–7). My estimate for a complete support proof within the current book's toolkit is **≈ 10–15%**; it would require a new theorem (support of the vanishing cohomology of a hyperplane section of a join, or of a family of plane sections of `Φ_3`), not a re-application.

**What a Lean proof would require.**
* *Constructive*: a `Witness` with exact rational `D`, `N_0..N_5`, `Q`; the verifier (`HodgeComputer.Verify`) checks the cleared identity by polynomial normalization (at `d = 10` an identity with ~85k monomials and 7^7-denominators — heavy but the pipeline `hodge/leangen.py → leancheck.py → sandboxcheck.py` is in place and passed on a regenerated accepted witness), homogeneous degrees, `D ≠ 0`, and **coprimality over `C`** (not inferable from rational irreducibility; since `F` is irreducible over `C`, `F ∤ D` suffices, provable by exhibiting a point of `X` where `D ≠ 0` only if irreducibility of `F` is available — check how the accepted degree-≥6 witnesses discharge this, e.g. `AlgebraicSpecialization`).
* *Nonconstructive*: the protected `HasCertificate` existential via the Chapter 1 bridge; a geometric support theorem would need to be reduced to algebraic identities as in `HoldoutCover` (polynomial family + integration + trace) or `EllipticCompanions` (contraction) — i.e. explicit incidence seeds, dominance (`CoordinateIndependence`), source transport, and legal clearing. For Experiment 7 this means: an explicit polynomial-curve family on `W_Y` (or on `Ĵ`), the residue-transport identities for the hyperplane section `λ = L`, and a deck/character argument isolating `X0^3Y0^3`; none of these exist yet.

**Bottom line.** This instance has no usable symmetry and no low-degree denominator; the only leverage is the eight-form/Fermat-section geometry, which so far yields many lines and sections but no supported divisor. The cheapest high-information steps are Experiment 1 (closing the degree classification) and Experiment 6(a)/7(a) (one-hour symbolic checks that decide whether the cover and SK routes are alive), followed by Experiment 3.

---

## Appendix A — reproduction notes

* Engine: `~/devel/sextics-work/hodge/{core.py, fixed_d.py}` (engine v2; prime list `core.PRIMES`), helpers `/tmp/hodge/results/short-power-presentation/hlib.py`, scripts `exp0.py, exp1.py, exp2.py, e1kernel.py, e1kernel2.py (parallel agent's version; my version used keys parity_D/prime), geom1.py`.
* `geom1.py` (this session) outputs: `dim J_6^+ = 18; B in J_6^+? False; e(X)+e(Y) in J_6^+? False; Fermat P^3 ∩ Λ projective dimensions: {1: 136080}; distinct lines 136066; t = M/L on H_ζ: I(−2+√3) with (1+t)^6+(1−t)^6 = 0`.
* Sweep logs: `/tmp/hodge/ext/hodgeHoldout/sweep/logs/short-power-presentation.log` (106 lines, all `no`), prime `268435399`.
* Thread: `/tmp/hodge/discussions/short-power-presentation.md` (4 KB, one comment).
* The E1 test used (`fixed_d.e1_solvable`) solves `q·N(F) = A·D` with `deg N_i = d+1` and forced parities; the `F·G` term of `J` is absorbed because `F = Σ Z_i∂_iF/6`.
