# Research dossier: `diagonal-permutation` (Simultaneous permutations · standard sector)

Status: OPEN (no `proof.json`, no `Solution.lean`, no accepted proof as of 2026-11-01).
Discussion: https://github.com/hodge-computer/sextics/discussions/177 (one top-level comment, no replies).
Added to the catalogue in commit `a3c8c1b` (2026-10-21, Benjamin Diamond, "Add structural frontier instances"), together with the three sibling frontier instances `short-power-presentation`, `quadratic-cubic-pullback`, `eight-plane-configuration`. All four are screened by `scripts/audit-structural-frontier.py` and all four are still open.

Resumed and completed 2026-11-01 after a crash of the first session (placeholders for the modular kernel results filled with the a-priori data and the exact status of the running computation; see §4.1–4.4).

Sources read in full: `data/diagonal-permutation/{instance.rec, why.md, provenance.json, checks.json}`, `scripts/audit-structural-frontier.py`, the whole discussion thread #177 (`/tmp/hodge/discussions/diagonal-permutation.md` and the raw JSON), `data/symmetric-sextic-residual/*` (the only instance named in `why.md`; it is itself open and its thread #117 is empty), `docs/01`, `docs/02` (§§2.1–2.9), `docs/03` (all), `docs/03a` (§3A.1, §3A.6), `docs/04` (all), `docs/05`, `formal/README.md`, `CONTRIBUTING.md`, and the `solution.md` of the accepted instances with the closest techniques (listed in §2). Everything marked **[verified]** below was recomputed by me in exact arithmetic (sympy 1.12) or modulo a large prime (Singular 4.3.2 / python-flint 0.9); scripts are in `/tmp/hodge/work/dp/` (`structure.py`, `equivariant.py`, `kernels.py`, `candidates.py`, `jac.sing`, `nf_all.sing`). Nothing in the repository was modified and nothing was posted.

Notation. `Z = (X0,X1,X2,Y0,Y1,Y2)`, `τ(X,Y) = (X,−Y)`. A certificate for `(F,q,A)` is a nonzero homogeneous `D` coprime to `F`, six homogeneous `N_i` of degree `deg D + 6q − 5`, and a quotient `Q` with
`E := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^(q+1)·Q`.
Normalized form (Chapter 1): `E1: q·N(F) = A·D`, `E2: D·div N − N(D) = F^q·C`. Here `q = 1`, so `deg N_i = d+1`, `deg Q = 2d − 6`, `deg E = 2d + 6` (`d = deg D`). `R_F = S/J(F)` is the Jacobian ring; `(R_F)_6^−` (Y-odd, degree 6) is the primitive τ-invariant `(3,1)` space, of dimension 208 for every smooth even sextic.

---

## 1. Exact target

### 1.1 The sextic

`q = 1`. The record `F` has 24 terms. Its sparse presentation (the one used by the audit script, **[verified]** equal to `instance.rec`):

```
F = Σ_{i=0}^{2} (X_i^6 + Y_i^6)
  + (5/7)  · Σ_{(i,j,k) ∈ S3}  X_i^3 X_j Y_i Y_j        (orbit O1, 6 terms, Y-degree 2)
  + (5/11) · Σ_{(i,j,k) ∈ S3}  X_i^3 X_j Y_k^2          (orbit O2, 6 terms, Y-degree 2)
  + (5/13) · Σ_{(i,j,k) ∈ S3}  Y_i^3 Y_j X_k^2          (orbit O3, 6 terms, Y-degree 4)
```

where the sums run over the six permutations `(i,j,k)` of `{0,1,2}` (so `k` is the index complementary to the ordered pair `i ≠ j`). Monomial support by Y-degree: 3 monomials of Y-degree 0 (`X_i^6`), 12 of Y-degree 2, 6 of Y-degree 4, 3 of Y-degree 6 (`Y_i^6`). Every monomial has even Y-degree (globally Y-even), but `F` is not even in any single coordinate (each `O1`, `O2`, `O3` contains odd powers of every variable), and no pair `(X_i,Y_i)` carries a paired-sign symmetry (`checks.json`: `coordinate_reflections: []`, `coordinate_paired_signs: []`; **[verified]** by enumerating all 64 sign patterns: the only sign symmetries are `±1` and `±τ`).

Structure of the couplings. Each "pair" `(X_i,Y_i)` is Fermat on its own; `O1` couples pairs `i,j` through `X_i^2·(X_iY_i)(X_jY_j)`; `O2`, `O3` couple all three pairs (`X_i^3X_j` times the square of the third pair's `Y`, and the mirror with `X`/`Y` roles exchanged but *different* coefficient, so there is no `X ↔ Y` exchange symmetry). The coefficients `5/7, 5/11, 5/13` are chosen so that `F ≡ Fermat (mod 5)`.

### 1.2 The target and its degrees

```
A = X0^3·Y0^3 − X1^3·Y1^3  =  T_0 − T_1,     T_i := X_i^3 Y_i^3,     deg A = 6 = 6q,  Y-degree 3 (odd).
```

A is a difference of two monomials; each `T_i` has all exponents ≤ 4, so `A` is a nonzero element of the Fermat Jacobian ring `F_5[Z]/(Z_i^5)` and hence (smooth complete-intersection base change) a nonzero primitive `(3,1)` class for the actual `F` (`checks.json`: `primitive_target`). Degrees at `q=1`: `deg N_i = d+1`, `deg Q = 2d−6` (so any certificate with `d ≤ 3` is automatically strict), `deg E = 2d+6`. Component parities for even `D`: `N_{X_i}` Y-odd, `N_{Y_i}` Y-even, `Q` Y-odd.

### 1.3 Smoothness **[verified]**

Modulo 5 all coupling coefficients vanish and `F ≡ Σ Z_i^6`, whose six partials are `6Z_i^5`; the Jacobian ideal is `(Z_0^5,…,Z_5^5)`, zero-dimensional, so the reduction is smooth and smoothness lifts to characteristic zero (`checks.json`: "Good reduction to Fermat modulo 5"). Independently, over `F_p` with `p = 1000003` the Jacobian ideal `J(F)` has `vdim = 15625 = 5^6` (Singular), i.e. the partials form a regular sequence there too, with Hilbert series `(1+t+t^2+t^3+t^4)^6`.

### 1.4 Symmetry **[verified]**

* Linear monomial symmetries. Enumerating all 720 coordinate permutations with all sign patterns: the signed permutations preserving `F` up to scalar are exactly the 12 projective elements of
  `G = S3 × ⟨τ⟩`, where `S3` acts by *simultaneous* permutation `(X_i, Y_i) ↦ (X_{σ(i)}, Y_{σ(i)})`. There is no `X ↔ Y` exchange and no reflection.
* Projective diagonal group: the exponent-difference lattice of `supp F` has Smith normal form `diag(1,1,1,1,2,0)`, so the diagonal substitutions preserving `F` up to scalar form (modulo scalars) a group of order 2: `{1, τ}`. Residue character of `A` under `τ`: trivial (as required). There is no other diagonal character to decompose against.
* Continuous symmetries: the Lie algebra of `{M ∈ gl_6 : (MZ)·∇F = λF}` is one-dimensional (scalars only), so the full linear automorphism group of `F` is finite. (Whether it contains non-monomial elements beyond `G` was not computed; the scalar Hessian symmetrizer below makes a hidden block structure impossible.)
* Hessian symmetrizer: dimension 1 modulo 101 and modulo 1009 (`checks.json`: 1 mod 101), so there is no linear direct-sum decomposition of `F`.
* Transpositions `s = (01)` etc. have determinant `+1` on `C^6` and fixed subspace `{X_0=X_1, Y_0=Y_1}` of codimension 2 — not reflections. Three-cycles have eigenvalues `1, ω, ω²` each with multiplicity 2.

### 1.5 The representation of the residue **[verified]**

Since every element of `S3` has determinant `+1` on the six coordinates, the residue `Res(A·Ω/F²)` transforms exactly like the polynomial `A`. Under `S3`,
`A = T_0 − T_1`, `(01)·A = −A`, `(02)·A = T_2 − T_1`, `(12)·A = T_0 − T_2`, three-cycles give `T_1 − T_2`, `T_2 − T_0`; the `S3`-average of `A` is `0`. So `A` spans, with its translates, the two-dimensional standard representation `V_std = {Σ a_i T_i : Σ a_i = 0}`. Combined with `τ` (determinant `−1`, `A` odd): the residue is `τ`-invariant and lies in the `V_std ⊗ 1` isotypic component of `H^4(X)^τ`.

Decomposition of the τ-invariant primitive cohomology under `S3` (equivariant Koszul series `det(1 − t^5 ρ(g))/det(1 − t ρ(g))`, **[verified]** totals 208/866/218):

| space | dim | trivial | sign | standard (multiplicity, ×2 = dim) |
|---|---:|---:|---:|---:|
| `(R_F)_6^−` = primitive `H^{3,1}` invariant | 208 | 39 | 31 | 69 (138) |
| `(R_F)_12^−` = primitive `H^{2,2}` invariant | 866 | 154 | 136 | 288 (576) |
| `(R_F)_6^+` = even sextic deformations mod `J` | 218 | 46 | 28 | 72 (144) |

Cross-check by hand: the Y-odd sextic polynomials `S_6^−` (226) have `S3`-character `(226, 10, 1)` on (identity, transposition, three-cycle) — 10 transposition-fixed and 1 three-cycle-fixed monomials — hence `43·triv + 33·sign + 75·std`; the 18 odd Jacobian directions `Y_jF_{X_i}, X_iF_{Y_j}` form two copies of `perm ⊗ perm = 2·triv + sign + 3·std`; subtracting gives `39, 31, 69` as in the table.

So `A` lies in a 138-dimensional standard-isotypic subspace of the 208-dimensional invariant `(3,1)` space; the `(01)`-invariant part of that subspace (the "`B`-sector" of the discussion) has dimension 69. The sector that the known invariant-quotient machinery could reach at best (trivial isotypic) has dimension 39; the sign sector has dimension 31.

### 1.6 Age screen (why no finite quotient is uniruled by the Kollár–Larsen criterion)

From `checks.json`, **[verified]** by hand: for every nontrivial `g ∈ G` and every fixed stratum on `X`, the age is ≥ 1.

| element | fixed strata in `P^5` (eigenvalue weight mod 6, dim of eigenspace) | age on `T_X` |
|---|---|---|
| `τ` | `{Y=0}`, `{X=0}` (two `P^2`) | 3/2, 3/2 |
| transposition `s` | weight 0, dim 4 (`P^3 = {X_0=X_1,Y_0=Y_1}`, meets `X` in a sextic surface); weight 3, dim 2 (`P^1`) | 1, 2 |
| `τs` | weight 0, dim 3; weight 3, dim 3 (two `P^2`) | 3/2, 3/2 |
| three-cycle `c` | weights 0, 2, 4, each dim 2 (three `P^1`) | 2, 2, 2 |
| `τc` | no fixed point on `X` (the six isolated fixed points of `P^5` are coordinate points, not on `X`) | — |

Hence Kollár–Larsen ("some fixed point of age `< 1` ⇒ quotient uniruled") never applies; every quotient `X/⟨g⟩` has canonical singularities and torsion canonical class, so none is uniruled.

### 1.7 Immediate consequences for the certificate problem (rigorous, elementary)

* **Normalization WLOG.** For any certificate `(D,N)` of `A` and any `g ∈ G`, `(D∘g, g^{-1}·(N∘g))` is a certificate of `A∘g` (chain rule; permutations and `τ` are linear). Multiplying numerator and denominator by the same polynomial changes nothing, so **one may assume `D` is `G`-invariant** (replace `D` by `∏_{g∈G} D∘g`, degree ×12; each factor is coprime to `F` because `F` is `G`-invariant). Invariance of `D` is therefore not a restriction of the existence problem, only a (12-fold) degree cost.
* **Reduction to one transposition-fixed source (benediamond, §3 below; [verified]).** With `B_i = T_0+T_1+T_2 − 3T_i`: `A = (B_1 − B_0)/3`, `B_2 = T_0 + T_1 − 2T_2` is fixed by `(01)`, and the `B_i` are `S3`-conjugate. A certificate for `B_2` gives certificates for `B_0, B_1` by permutation and hence for `A` by the linear combination (common denominator = product, or the same `D` if `D` is `S3`-invariant). Conversely a certificate for `A` gives one for every `B_i`. So the problem is equivalent to serving `B := B_2` with `(01)`-equivariant data; the `(01)`-average of a certificate of `B` is again one.
* **E1 is a Jacobian-ring condition.** Reducing `E` modulo `F` gives `D·(A·D − N(F)) ∈ (F)`, hence (`D` coprime to the irreducible `F`, `F ∈ J(F)` by Euler) `[A]·[D] = 0` in `R_F`. So every admissible denominator lies in the kernel of multiplication by `[A]: (R_F)_d → (R_F)_{d+6}`. This kernel is computed exactly in §4.1 below.

---

## 2. Companion and related instances

### 2.1 What `provenance.json` and `why.md` point to

* `provenance.json` has **no** `sameVarietyInstances`; references are `docs/03-geometric-support.md` and `docs/04-known-locus-and-frontier.md`; `existenceStatus`: "No certificate or existence conclusion for this target is supplied by the screened constructions."
* `why.md` contrasts the instance with the **separate-permutation example `symmetric-sextic-residual`** (discussion #117). That instance is **also open** and its discussion thread is empty. Its data: `F` = Fermat + `(5/7)`·(products of elementary symmetric quadrics `Q_X, Q_Y`), symmetric under `S3 × S3` acting separately on `X` and on `Y`; `A = X0X1X2Y0Y1Y2`; `q = 1`. There the residue transforms by `sign ⊗ sign` (an `X`-transposition has determinant `−1`, so the polynomially symmetric `A` has anti-invariant residue); the quotient by the elementary symmetric functions (weights `(1,2,3,1,2,3)`, Jacobian `Δ_X Δ_Y`) serves multiples of `Δ_X Δ_Y`, the opposite character; `τ` is not in the permutation group so the reflection-group theorem (§3.3) does not apply; Hessian symmetrizer scalar; projective diagonal group of order 2. So the two instances are "siblings" of the same design (a permutation symmetry whose averaging projector misses the target), but the present instance is harder in one specific way: its missing component is a genuinely two-dimensional irreducible representation, whereas `symmetric-sextic-residual` misses a one-dimensional character (`sign ⊗ sign`) that an appropriate double cover of the quotient (discriminant cover) could in principle carry.
* The three sibling frontier instances introduced in the same commit (`short-power-presentation`, `quadratic-cubic-pullback`, `eight-plane-configuration`) all have target `X0^3Y0^3` and are all open; none shares the variety.

### 2.2 Why every mechanism of the book fails here (from `why.md`, `checks.json`, Chapters 2–4)

| mechanism | requirement | status for this `F`/`A` |
|---|---|---|
| cubic-surface / six-line construction (§2.4–2.5) | `F` even in enough coordinates (coordinatewise-even, or one `X` and one `Y` reflection) | no coordinate is even |
| conic construction, Piola / weighted-Fano density (§2.6, §3.4) | `F = G(Φ)` for a nonlinear polynomial map; serves only `det(dΦ)·b(Φ)` | no such presentation known; target not a density |
| reflection quotient (§3.3) and reflection-generated groups | a reflection in the symmetry group containing `±τ` | transpositions have fixed codimension 2; no reflection |
| Fermat cover / moving conic (§2.7, `holdout`) | a cover `T^6 = Y_a^6 + Y_b^6` making the equation symmetric in `(Y_c, T)` | couplings `Y_iY_j`, `Y_k^2`, `Y_i^3Y_j` are not of that shape |
| paired-sign cubic → Jacobian transport (§2.9, §4.2, `cycle-three`, `mixed-biweight`, `mixed-polar`) | paired signs `(X_i,Y_i) ↦ (−X_i,−Y_i)` and the diagonal weight action | `O2`, `O3` break every paired sign |
| diagonal cyclic quotients (§3.4, §4.3, `klein-*`, `cyclic-triple-fano`) | a nontrivial diagonal symmetry | projective diagonal group is `{1,τ}` |
| split / crossed / polarized blocks (§3A) | `F = f(X)+g(Y)` or block structure | Hessian symmetrizer scalar: no linear block decomposition |
| Kollár–Larsen low-age quotient (§3.3 end) | a fixed point of age `< 1` | all ages ≥ 1 (table in §1.6) |
| invariant quotient `X/S3`, `X/⟨(01)⟩`, `X/A3` via averaging ideal (§3.1) | the served representation must have `ρ^K ≠ 0` and the quotient must have supported `H^4` | `V_std^{S3} = V_std^{A3} = 0`; `V_std^{(01)} ≠ 0` but `X/(01)` is birationally a family of Calabi–Yau threefolds over `P^1` (see §3), not Fano, not uniruled (fixed surface of age 1 ⇒ canonical singularities, `K` torsion) |
| natural cubic cover and its `S3` Galois closure (§4.7) | uniruled intermediate quotients | the subgroup quotients serving `V_std` (`K = 1, C_2`) are not uniruled |
| generic-target theorem (§4.4) | a dominating family over the full even-sextic space | `F` is a single special point of a 50-dimensional invariant family; not applicable |

### 2.3 Accepted instances with the nearest techniques (for transfer), with denominators

None of the 66 accepted solves uses a nonabelian permutation symmetry; every accepted construction relies on a reflection, a paired sign, a diagonal cyclic symmetry, or a block split, all absent here. The closest relatives:

| instance (author, PR) | why it is relevant | method | denominator |
|---|---|---|---|
| `klein-hesse` (κ=2; `A = X0^3X1^2Y0 − X1^3X2^2Y0`) | the **only accepted target that is a zero-sum ("Fourier") combination under a cyclic permutation** (`Z/3` on `U,V,W = X0^3X1^2, X1^3X2^2, X2^3X0^2`) — the abelian analogue of a standard-representation target | monomial map to Hesse coordinates, logarithmic and complementary Euler fields, rational monomial Piola transport, Fourier channels | `D_κ = 7 r^3 f_κ − p^3`, `d = 15`, strict, explicit witness |
| `klein-rational` (`A = X0^3X1^2Y0`) | same cyclic Fourier machinery | Fourier channels + normalized trace | `D = Q^3 + 9(U+V+W)^3 P`, `d = 30`, strict |
| `one-reflection` (`fqlx`, PR 161) | dense coupling, one hidden reflection; binary-cubic Hessian construction | two half-fields, `√δ`-averaging | `d = 60` (two halves of degree 30), strict, over `Q(√3)` |
| `holdout` (`fqlx`, PR 128), `holdout-ternary-deformation`, `holdout-moving-branch` | create a reflection on a finite cover; moving conic; trace | existence only (finite-algebraic), no expanded `D` | unknown degree |
| `cycle-three`, `mixed-biweight`, `mixed-polar`, `holdout-moving-branch` (`hodgepodgeconjecture`, `fqlx`) | paired-sign cubic, Jacobian residue transport, first variation | existence only | unknown degree |
| `cyclic-triple-fano` / open companion `cyclic-triple-residual` | a solved density sector and an open complementary-character sector on the *same* variety: the same "character gap" as here | free weighted line on the cubic cover | existence only |
| `reflection-coupled-fano` / open companion `reflection-coupled-residual` | idem with a reflection quotient | modular quartic, `24×24` incidence, trace | existence only |
| `two-ternary-dense-two-odd-hessian-d48` | dense blocks, explicit large non-strict witness | ternary cubic covariants | `d = 55`, non-strict |
| `hesse-triangle`, `family-hesse-web` | Hesse pencil in a block, norms of contact factors | three-quadric Piola, cubic algebra norm | `d = 6`, strict |
| `fermat-pair-rational-d4-strict`, `fermat-q2` | the `c → 0` limit of the family `F_c = Fermat + c·(couplings)` is Fermat, where everything is solved (coordinatewise-even) | quartic cofactor `X0^4 − X0^2X1^2 + X1^4`; subset Euler | `d = 4`; `d = 6` |

The precise reason the `klein-*` approach does not transfer: those sextics are split (`f(X) + g(Y)`) with a *diagonal* `Z/7` and a monomial map whose Jacobian is a monomial; here the three-cycle acts by a permutation, the natural "Fourier coordinates" `ξ_k = Σ_j ω^{jk} X_j`, `η_k = Σ_j ω^{jk} Y_j` make `F` dense, and there is no monomial map.

---

## 3. The discussion thread, comment by comment

Thread #177 "Simultaneous permutations · standard sector". Exactly **one** comment, no replies, no edits after posting.

### 3.1 `benediamond`, 2026-10-23T17:09:14Z — "A source-preserving reduction to one transposition quotient"

Provenance stated in the comment: "ChatGPT-assisted investigation for `benediamond`, including separate AI-agent checks; no independent human review is claimed"; "no certificate or Lean theorem claimed"; starting framework = `docs/03-geometric-support.md` at commit `faa9df2`. Links: the pinned `instance.rec` and `docs/03` at `faa9df2ec6406e4eacc135f677cf3711315d0e41`. No reproduction packet, no code, no external URL.

Claims, with my verification status:

1. **Representation-theoretic reduction (proven; [verified]).** `T_i = X_i^3Y_i^3`, `B_i = T_0+T_1+T_2 − 3T_i`, `A = (B_1 − B_0)/3`; each `B_i` is fixed by the transposition of the other two indices; the `B_i` are `S3`-conjugate. Group-algebra identity
   `e_std = (2/3) Σ_{s transposition} (1+s)/2 − 2 e_{S3}`
   ([verified]: both sides equal `(1/3)(2e − c − c²)` with `c` a three-cycle). Consequence: it suffices to support `B = B_2 = T_0 + T_1 − 2T_2` on the quotient by `(01)`, pull back, apply the two permutations, and take the linear combination; simultaneous permutations have determinant 1 so the identity holds for residue classes, not only numerators. **This is correct and also holds at the certificate level** (§1.7). It does *not* supply support; it only concentrates the search on one source.
2. **Explicit quotient chart including the density ([verified]).** `X0 = (p+r)/2, X1 = (p−r)/2, X2 = a, Y0 = (q+tr)/2, Y1 = (q−tr)/2, Y2 = b, u = r²`. The transposition negates `r` and fixes `p,q,a,b,t,u`. `F` becomes even in `r` and descends to `F_s(p,q,a,b,u,t)` of weighted degree 6 with weights `(1,1,1,1,2,0)`. The descended numerator is
   `B_s = [(pq+tu)^3 + 3(pq+tu)(pt+q)^2 u]/32 − 2a^3b^3`
   ([verified] symbolically). On `a=1`, `dX0∧dX1∧dY0∧dY1∧db = −(1/8) dp∧dq∧du∧dt∧db` ([verified]: `det ∂(X0,X1,Y0,Y1)/∂(p,q,r,t) = −r/4`, times `dr = du/(2r)`). So the quotient source is the residue of `−B_s dp dq du dt db/(8F_s²)`. The global sign `τ` acts by `(q,b,t) ↦ (−q,−b,−t)`; `B_s` and the five-form are both odd, so the residue stays invariant ([verified] parity).
   Geometric reading (mine): for fixed `t`, `F_s = 0` is a degree-6 hypersurface in `P(1,1,1,1,2)` with `K = O(6−6) = O`, i.e. `X/⟨(01)⟩` is birationally a one-parameter family (parameter `t = (Y_0−Y_1)/(X_0−X_1)`) of **Calabi–Yau threefolds**, not a Fano quotient. This is consistent with `why.md` and the age table.
3. **First surface fibration retains three holomorphic forms (plausible, standard; partially [verified]).** Fixing `t` and `λ = b/a` gives sextic surfaces in `P(1,1,1,2)` with coordinates `(p,q,a,u)`. Modulo 5 the fiber equation is `(p^6+q^6+(1+t^6)u^3)/32 + (1+λ^6)a^6` ([verified], with `1/32 ≡ 3 mod 5`); at `t = λ = 0` the weighted cone is smooth away from the vertex and the surface misses the vertex; openness of smoothness gives a smooth generic characteristic-zero fiber; adjunction gives `K = O(1)` and `p_g = h^0(O(1)) = 3`. Hence the surface-fibration theorem (`p_g = q = 0`, §3.2) does not apply to this fibration. Stated explicitly as "not an obstruction to the target": another fiber component, fibration, or correspondence could still support `B`.
4. **The exceptional divisor over the fixed surface does not help (plausible).** The blowup summand `H^2(S)(−1)` over the transposition's fixed surface `S = X ∩ {X_0=X_1, Y_0=Y_1}` (a sextic surface in `P^3`) is killed by pushforward to `X`; a correspondence mixing it with the original primitive summand "would need an explicit source-action calculation" — not done.

What the thread does **not** contain: no denominator search, no linear-algebra exclusion, no numerical experiment, no proposed denominator family, no Lean work, no reproduction packet, no second contributor. The rigorous content is item 1; items 2–4 are correct scoping calculations that close one geometric route (the obvious `(01)`-quotient fibrations) without opening another.

### 3.2 Related mentions elsewhere

GitHub search for `diagonal-permutation` across discussions returns only #177 (the other hits #44, #61, #115 are fuzzy matches; grepping the saved threads and the full #44 thread finds no mention). No PR or issue mentions the instance.

---

## 4. Computational evidence produced for this dossier

All computations are in `/tmp/hodge/work/dp/`; prime `p = 1000003` where modular.

### 4.1 Rigorous denominator-degree bound from `E1` **[verified modulo p; injectivity transfers to Q and C]**

For the smooth `F`, `J(F)` over `F_p` has `vdim = 15625` with the Koszul Hilbert function (so the computation is on a complete intersection with the same dimensions as in characteristic zero). The map `[A]·: (R_F)_d → (R_F)_{d+6}` was computed in the standard-monomial bases for `d = 0,…,14`, together with its restrictions to the Y-even and Y-odd parts and to the `G = S3×⟨τ⟩`-invariant part. Because the rank modulo `p` is at most the rank over `Q`, **a zero kernel modulo `p` is a zero kernel over `Q` and `C`**, while a nonzero kernel modulo `p` is only an upper-bound statement.

**Status of the exact modular computation (resumed run, 2026-11-01).** The first attempt (`kernels.py`, pure-Python loops over `nf_all.out`) did not finish before the session crashed; `nf_all.out` (normal forms of every monomial, needed for the symmetry sectors) is complete for `d = 0,…,13` and truncated inside `d = 14`, while `jac.out` (standard bases of `(R_F)_d`, `(R_F)_{d+6}` and `NF(A·k)` for every standard monomial `k`, `d = 0,…,14`) is complete. A numpy/flint rewrite `kernels2.py` (streams `nf_all.out`, builds the maps `S_d^+ → (R_F)_{d+6}^-`, `m ↦ [A·m]`, and the orbit-sum rows for the sectors `G`-invariant, `S3`-sign, `A3`-invariant, `(01)`-invariant, `(01)`-anti-invariant, `G`-invariant-Y-odd) was launched; it writes one line per degree to `/tmp/hodge/work/dp/kernels2.log` and a pickle `kernels2.pkl`. **At the time this dossier was written it had not yet emitted degree results**, so the table below is the *rigorous a-priori* information (equivariant Hilbert series, recomputed [verified]) and must be completed from `kernels2.log` by whoever resumes:

| `d` | `dim (R_F)_d` (all / Y-even / `G`-inv) | `dim (R_F)_{d+6}^-` (all / std-multiplicity) | forced kernel, all Y-even | forced kernel, `G`-inv | exact kernel (mod `p`) |
|---:|---|---|---:|---:|---|
| 6 | 426 / 218 / 46 | 866 / 288 | 0 | 0 | pending (`kernels2.log`) |
| 7 | 666 / 333 / 63 | 843 / 281 | 0 | 0 | pending |
| 8 | 951 / 483 / 94 | 744 / 248 | 0 | 0 | pending |
| 9 | 1246 / 623 / 115 | 623 / 207 | 0 | 0 | pending |
| 10 | 1506 / 762 / 144 | 468 / 156 | ≥ 294 | 0 | pending |
| 11 | 1686 / 843 / 153 | 333 / 111 | ≥ 510 | ≥ 42 | pending |
| 12 | 1751 / 885 / 167 | 208 / 69 | ≥ 677 | ≥ 98 | pending |
| 13 | 1686 / 843 / 153 | 123 / 41 | ≥ 720 | ≥ 112 | pending |

Interpretation rules (rigorous). `D` must be Y-even for the component parities to match (a Y-odd `D` is also admissible with swapped parities; its sector is reported separately). For `d ≤ 9` the map `[A]·: (R_F)_d^+ → (R_F)_{d+6}^-` goes from a smaller to a larger space, so a zero kernel modulo `p` there is a *proof* that no denominator of that degree exists (good reduction: `vdim J(F) mod p = 15625` = the characteristic-zero value, so `J(F)_n mod p` has the characteristic-zero dimension and membership over `Q` implies membership mod `p`). For `d ≥ 10` (all Y-even) and `d ≥ 11` (`G`-invariant) the kernel is forced nonzero by dimension count alone, so the **minimal possible denominator degree is ≤ 10 only if the exact kernel at some `d ≤ 9` is nonzero**, and the first degree where the `E1` constraint stops being a complete obstruction is `d_min ∈ [?, 10]` (all `D`) resp. `[?, 11]` (invariant `D`). Note also that the `E1` kernel is a necessary condition only; `E2` (`D·div N − N(D) ∈ (F)`) is the real test and is only run on `E1`-admissible `D` (§4.2–4.4).

Expected lower bounds from the equivariant Hilbert series alone (forced kernels, before any computation): among all Y-even denominators the kernel is forced nonzero from `d = 10` (`762 → 468`, so `≥ 294`); among `G`-invariant denominators it is forced from `d = 11` (`153` invariant classes against standard-multiplicity `111` in `(R_F)_17^−`, so `≥ 42`); at `d = 12` the forced invariant kernel is `≥ 98` out of `167`. The dimensions of `G`-invariant polynomials (Molien series) are `1,1,4,7,16,26,50,76,126,185,280,392,561,…` in degrees `0,…,12`; in particular the `G`-invariant *sextics* form the 50-dimensional family containing `F`.

### 4.2 Natural structured denominators tested for `E1`

`candidates.py`/`kernels2.py` test, for each structured `D`, whether `NF(A·D) = 0` modulo `J(F)` over `F_p` (exact, rigorous exclusion when nonzero, see §4.1). Candidates (all `S3×τ`-invariant and Y-even unless marked Y-odd), by degree:

* `d = 2`: `p_2(X) = ΣX_i²`, `p_2(Y)`, `p_{XY} = ΣX_iY_i`.
* `d = 4`: `p_2(X)²`, `p_2(X)p_2(Y)`, `p_{XY}²`, `ΣX_i⁴`, `ΣX_i²Y_i²`.
* `d = 6`: `Δ_X²`, `Δ_Y²`, `Δ_XΔ_Y` (Y-odd), `∏_{i<j}((X_i−X_j)² ± (Y_i−Y_j)²)`, `O1`, `O2`, `O3`, `ΣX_i⁶`, `ΣY_i⁶`, `ΣX_i⁶−ΣY_i⁶`, `p_2(X)³`, `p_2(X)²p_2(Y)`, `p_2(X)p_{XY}²`, `F − ΣX_i⁶ − ΣY_i⁶`, `ΣT_i` (Y-odd), `p_{XY}³` (Y-odd).
* `d = 8`: `Δ_X²·{p_2(X), p_2(Y), p_{XY}}`, `O_a·p_2`, `(ΣX_i⁶)·p_2`, `qp·p_2(X)`, `p_2(X)⁴`.
* `d = 10`: `O1·{p_2(X)², p_2(X)p_2(Y), p_{XY}²}`, `Δ_X²·{p_2(X)², p_2(Y)², p_{XY}²}`, `Δ_XΔ_Y·p_2(X)p_2(Y)`, `Δ_XΔ_Y p_{XY} p_2(X)`, `qp·p_2(X)²`, `qm·p_2(X)²`, `(ΣX_i⁶)·{p_2(X)², p_{XY}²}`, `O2·p_2(Y)²`, `O3·p_2(X)²`, `F_c·{p_2(X)², p_{XY}²}` (`F_c` = coupling part of `F`), `(ΣX_i⁶+ΣY_i⁶)p_2(X)²`, `p_2(X)⁵`.
* `d = 12`: `Δ_X²Δ_Y²`, `O_aO_b` (all six products incl. squares), `O_a·ΣX_i⁶`, `O_a·ΣY_i⁶`, `(ΣX_i⁶)(ΣY_i⁶)`, `(ΣX_i⁶)²`, `(ΣX_i⁶−ΣY_i⁶)²`, `qp·qm`, `qp²`, `qm²`, `Δ_X²·ΣX_i⁶`, `Δ_X²·ΣY_i⁶`, `Δ_Y²·ΣX_i⁶`, `Δ_XΔ_Y·O1`, `Δ_XΔ_Y·ΣX_i⁶` (Y-odd), `p_2(X)³p_2(Y)³`, `F_c²`, `F_c·ΣX_i⁶`, `F_c·ΣY_i⁶`, `F²` (control: `F² ∈ J(F)` by Euler, so `A·F² ∈ J`, but `F²` is not coprime to `F`), `A²`, `(ΣT_i)²`, `ΣT_i²`, `B_2² = (T_0+T_1−2T_2)²`, `Σ_iB_i²`, `p_2(X)²p_2(Y)²p_{XY}²`, `Δ_X²Δ_Y² − qp·qm`.

**Results: pending in `/tmp/hodge/work/dp/kernels2.log` (lines beginning `cand d=`).** Expected on dimensional grounds: every candidate with `d ≤ 9` fails `E1` (if the exact kernel there is zero, as the forced count suggests), and the floating-point complete screens of §4.4 already show all fifteen degree-6 candidates inconsistent. Any `d = 10, 12` candidate that *passes* `E1` must then be run through the complete `E1+E2` system (`fixedD2.py d10` / `d12`, §4.4) before being excluded or accepted.

### 4.3 Sizes of the complete fixed-denominator linear systems at `q = 1` (parity-reduced)

| `d` | unknowns in `N` (3 Y-odd + 3 Y-even components, degree `d+1`) | unknowns in `Q` (degree `2d−6`, Y-odd) | equations (degree `2d+6`, Y-odd) |
|---:|---:|---:|---:|
| 6 | 2376 | 226 | 16797 |
| 8 | 6006 | 1491 | 40326 |
| 10 | 13104 | 5796 | 84903 |
| 11 | 18564 | 10152 | 118608 |
| 12 | 25704 | 16797 | 162248 |
| 14 | 46512 | 40326 | 287793 |
| 18 | 127512 | 162248 | 766843 |
| 24 | 427518 | 766843 | 2502990 |

Using `S3`-equivariance (target in `V_std`, `D` invariant) divides every count by roughly 3; the `(01)`-reduction to `B` with `(01)`-equivariant numerators divides by 2. Dense rank computations modulo `p` with python-flint are comfortable up to `d ≈ 10` (minutes) and feasible to `d ≈ 12` with the equivariant reduction; beyond that, sparse elimination or Singular/msolve modules are needed.

### 4.4 Complete fixed-denominator screens (numerical, floating point)

Script `fixedD2.py`: the full parity-reduced system of §4.3 for a given `D` (all six numerator components and the quotient `Q`, no symmetry restriction) is assembled as a sparse real matrix and solved in least squares (scipy LSQR, `atol = btol = 1e−14`). The reported number is the relative residual `|Mx − b|/|b|`; a consistent system gives `~1e−15` (positive control: the accepted `fermat-pair-rational-d4-strict` witness gives `2.5e−16`), an inconsistent one gives a residual of order `1e−1` (negative control: `T_0 = X_0^3Y_0^3` on Fermat over `X_0^4 − X_0^2Y_0^2 + Y_0^4`, which is excluded by hand, gives `8.7e−2`). **These are floating-point screens, not proofs**; at `q = 1` an exact failure of the complete system over `D` would exclude the whole support `V(D)` (prescribed-support theorem), and that exactness still has to be established modulo two primes or over `Q` before anything is posted.

Degree 6 (2602 unknowns, ≈16 800 equations each; every candidate is `G`-invariant unless marked):

| `D` (degree 6) | relative residual |
|---|---:|
| `Δ_X²` | 0.46 |
| `Δ_Y²` | 0.42 |
| `Δ_XΔ_Y` (Y-odd `D`; component parities swapped) | 0.39 |
| `∏_{i<j}((X_i−X_j)² + (Y_i−Y_j)²)` | 0.50 |
| `∏_{i<j}((X_i−X_j)² − (Y_i−Y_j)²)` | 0.40 |
| `O1`, `O2`, `O3` | 0.57, 0.58, 0.59 |
| `ΣX_i^6`, `ΣY_i^6`, `ΣX_i^6 − ΣY_i^6` | 0.43, 0.27, 0.51 |
| `(ΣX_i²)³`, `(ΣX_i²)²(ΣY_i²)`, `(ΣX_i²)(ΣX_iY_i)²` | 0.53, 0.42, 0.39 |
| `F − ΣX_i^6 − ΣY_i^6` (the coupling part) | 0.51 |

All fifteen are far from consistent, as the `E1` kernel in degree 6 already predicts (§4.1).

**Degrees 10 and 12.** The batches `fixedD2.py d10` (12 candidates: `O1·p_2(X)²`, `O1·p_2(X)p_2(Y)`, `O1·p_{XY}²`, `Δ_X²p_2(X)²`, `Δ_X²p_2(Y)²`, `Δ_X²p_{XY}²`, `Δ_XΔ_Y p_2(X)p_2(Y)`, `qp·p_2(X)²`, `(ΣX_i⁶)p_2(X)²`, `(ΣX_i⁶)p_{XY}²`, `O2·p_2(Y)²`, `O3·p_2(X)²`; systems `18900 × 84903`) and `fixedD2.py d12` (14 candidates: `Δ_X²Δ_Y²`, `O1O2`, `O1O3`, `O2O3`, `(ΣX_i⁶)(ΣY_i⁶)`, `O1²`, `O1·ΣX_i⁶`, `O1·ΣY_i⁶`, `qp·qm`, `qp²`, `Δ_X²ΣY_i⁶`, `Δ_XΔ_Y·O1`, `p_2(X)³p_2(Y)³`, `(ΣX_i⁶−ΣY_i⁶)²`; systems `42501 × 162248`) were launched by the first session but **did not run** (`screen_d10.log`, `screen_d12.log` are empty; the session crashed and the Python-list COO assembly in `fixedD2.py` needs several GB per candidate at `d = 10`, which was not available on the shared machine). They are superseded by the exact `E1` test of §4.2, which is cheaper and rigorous: a complete screen is only needed for candidates that pass `E1`. **No complete `E1+E2` fixed-denominator solve has been completed for any `D` of degree `≥ 7` on this instance by anyone.**

---

## 5. What remains unexplored (concrete gaps)

1. **No fixed-denominator solve of any kind has been reported** by anyone for this instance, at any degree. The `E1` computation above is the first exclusion data. The complete `E1+E2` fixed-`D` test has not been run for a single candidate (§4.3 gives the sizes).
2. **No joint (bilinear) search in `(D,N,Q)`**, numerical or exact, in the minimal admissible degrees `10–12`, with or without the `S3×τ`-invariant restriction.
3. **The `(01)`-reduction has not been pushed to a source/primitive construction**: nobody has looked for a curve family, a second fibration with `p_g = q = 0` fibers, or a correspondence on the CY3-fibered quotient `W = {F_s = 0}` (weights `(1,1,1,1,2)` over the `t`-line), nor computed the Hodge numbers of the generic CY3 fiber or of the `(3,3)`-type surfaces obtained by fixing two parameters other than `(t, λ)`.
4. **The three-cycle / `A3` Fourier route** (the abelian analogue used in `klein-hesse`) has not been tried: in Fourier coordinates `ξ_k, η_k` the standard representation splits into the `ω` and `ω²` characters of `A3`, and `A` becomes a sum of two character components, each a target for a *cyclic* group — but a non-diagonal one.
5. **Deformation/continuation from a solved neighbour.** The family `F(a,b,c) = Fermat + a·O1 + b·O2 + c·O3` contains three solved loci one or two coefficients away (**[verified]**, see E8): Fermat (`a=b=c=0`), the separately-`Y`-even line `a=c=0` and the separately-`X`-even line `a=b=0`; the instance is `(a,b,c) = (5/7, 5/11, 5/13)`. No explicit certificate for `T_0 − T_1` (or `T_0`) at any of these starting points has been recorded, and the first-order deformation problem (linear in the correction terms once a starting certificate is fixed) has not been set up by anyone.
6. **No `sector theorem` for the 50-dimensional family of `G`-invariant even sextics**: a construction valid on an open of that family would be a new row in the §4.1 coverage table; the §4.4 generic-target argument has not been adapted to the invariant deformation space (where the Higgs field only uses the 46 invariant deformation directions and the generated algebra commutes with `S3`, so at best it could act irreducibly on each isotypic block).
7. **Covers creating a reflection**: no cover analogous to `T^6 = Y_0^6 + Y_2^6` (§2.7) has been proposed; the obvious obstacle is that the couplings are not symmetric functions of `(Y_c, T)`.
8. **Correspondences from the fixed surfaces** `S_{ij} = X ∩ {X_i = X_j, Y_i = Y_j}` (sextic surfaces in `P^3`, `p_g = 10`) and from the fixed curves of the three-cycles have not been examined beyond the remark that the blowup summand pushes forward to zero.
9. **The sibling `symmetric-sextic-residual`** (one-dimensional missing character) is untouched; a solution there (e.g. through the discriminant double cover of the `S3×S3` quotient) would be the natural warm-up.
10. Nothing at pole order 2 is needed for this instance, but the `E1` kernel at `q = 2` (denominators `D²`) has not been computed either.

---

## 6. Candidate attack plan

Degrees below use `q = 1`: `deg N_i = d+1`, `deg Q = 2d−6`, `deg E = 2d+6`.

### E1. Finish the exact `E1` picture and extract explicit admissible denominators (cheap, rigorous)
*Idea.* The kernel of `[A]` on `(R_F)_d^{G}` (and on `(R_F)_d^+`) is linear; its elements, lifted to polynomials, are the only possible denominators of degree `d` (modulo `J(F)_d`, which is also available as denominator freedom since `D` is only constrained through `[D]`… note that adding an element of `J(F)_d` to `D` changes `D` itself, so the true `E1`-admissible set in `S_d` is the preimage of the kernel, of dimension `dim ker + dim J(F)_d`).
*Construction.* Singular `std(J)`, `kbase`, `reduce`; python-flint rank; Reynolds operator over `G` for the invariant sector; `factorize` of kernel elements modulo `p` to look for products of lower-degree factors (branch strata, Fermat-like factors).
*Size.* Trivial (matrices ≤ `2000 × 1000`).
*Chance of a certificate.* Zero by itself; chance of a *rigorous* statement "no certificate has `deg D ≤ d_0`": certain, with `d_0` as in §4.1.

### E2. Complete fixed-`D` solves (`E1 + E2`) for structured invariant denominators at `d = 10,…,12`
*Idea.* At `q = 1` a failure of the complete numerator system over `D` excludes the whole support `V(D)` (prescribed-support theorem, Chapter 1), so each test is a rigorous exclusion of a geometric support, not just of a polynomial.
*Families to test (all `G`-invariant, all Y-even).* (a) `Δ_X²Δ_Y²` (`d = 12`), `Δ_X²Δ_Y²·(quadratic invariant)`, where `Δ_X = ∏_{i<j}(X_i−X_j)`; these vanish on the three transposition-fixed `P^3`'s. (b) Products of the three transposition quadrics `∏_{i<j}((X_i−X_j)^2 ± (Y_i−Y_j)^2)` (`d = 6`) times invariant quartics/sextics. (c) Products of orbit sums `O_a O_b`, `O_a·ΣX_i^6`, `(ΣX_i^6)(ΣY_i^6)` (`d = 12`) and `O_a·(invariant quartic)` (`d = 10`). (d) Random elements of the invariant `E1` kernel at `d = 11, 12` (control group). (e) Any kernel element that factors.
*System size.* From §4.3: `d = 10`: `18900` unknowns, `84903` equations; with `S3`/`(01)`-equivariance ≈ `6300 × 28300` dense modulo `p` — minutes with flint. `d = 12`: `42501 × 162248`, ≈ `14000 × 54000` after equivariance — tens of minutes to an hour.
*Caveat.* Inconsistency modulo one prime is evidence, not proof (a rational solution could have denominators divisible by `p`); repeat with a second prime or finish over `Q` for any exclusion worth stating publicly.
*Chance.* Low (< 5 %) that a structured low-degree `D` works — every accepted explicit witness with `d ≤ 15` came from a visible reflection/contact/cofactor structure that is absent here. High value as rigorous negative data and as a pipeline check for E3.

### E3. Joint bilinear search at the minimal admissible degree, restricted to the invariant `E1` kernel
*Idea.* Parametrize `D = Σ_k c_k D_k` over a basis `D_1,…,D_m` of the `G`-invariant `E1`-admissible polynomials of degree `d` (`d = 11` or `12`; `m` from §4.1 plus the `J(F)_d^G` freedom), lift `N = N_0(D) + K(c_{ij})` with Koszul syzygies `c_{ij}` of degree `d − 4`, and impose `E2: D·div N − N(D) = F·Q`. For fixed `c` the problem is linear in `(c_{ij}, Q)`; the solvable `c` form an algebraic subset of `P^{m−1}`.
*Method.* Variable projection (Chapter 5): minimize the least-squares residual of the linear subproblem over `c` (Gauss–Newton / Levenberg–Marquardt, many random starts), monitor rank drops and `D mod F ≠ 0`, then reconstruct over `Q` and re-solve exactly. Alternatively an exact Gröbner/elimination attack on the bilinear system restricted to a few `c`-directions.
*Size.* `d = 11`: `m ≈ 42–60`; linear unknowns ≈ `15 × 792/2 + 10152 ≈ 16000` (fewer with equivariance), equations `118608`. Each residual evaluation is a sparse least-squares of that size (seconds with LSQR); a search needs `10^3–10^5` evaluations.
*Chance.* Low-to-moderate (10 %) *if* a certificate of degree ≤ 12 exists at all; the a-priori probability that the minimal degree is that small is itself low (`one-reflection` needed 60, the dense `two-ternary` 55). The experiment is still the only generic method that does not presuppose a geometric mechanism.

### E4. Three-cycle Fourier channels (abelian analogue of the `klein-hesse` solve)
*Idea.* Let `ω³ = 1`, `ξ_k = Σ_j ω^{jk} X_j`, `η_k = Σ_j ω^{jk} Y_j` (`k = 0,1,2`). The three-cycle `c` acts diagonally with weights `(0,1,2,0,1,2)` on `(ξ_0,ξ_1,ξ_2,η_0,η_1,η_2)`, and `V_std|_{A3} = χ_ω ⊕ χ_{ω²}`; `A = A_ω + A_{ω²}` with `A_ω = (1/3)Σ_j ω^{-j}… ` the `c`-eigencomponents (each is a combination of `T_0,T_1,T_2` with coefficients `1, ω^{∓1}, ω^{∓2}` and is a legitimate target over `Q(ω)`; a certificate for `A_ω` plus its conjugate gives one for `A` by the trace argument of Chapter 1). In the `ξ,η` coordinates `F` is a `Z/3`-invariant sextic with a *diagonal* weight action, so the cyclic-quotient machinery (§3.3 end: residue character `Σ(e_j+1)w_j mod 3`; weighted quotient `P(…)` by the monomial invariants `ξ_1³, ξ_1ξ_2, ξ_2³, …`) can be written down. The residue character of `A_ω` under `c` is nontrivial by construction, so the invariant quotient does not serve it; but the quotient by `τc` (which acts **freely** on `X`, §1.6) and the age-2 fixed loci of `c` leave only non-uniruled quotients. The real content of `klein-hesse` was an explicit *rational Piola transport through a monomial map plus Fourier channels*; the corresponding attempt here would be: find a rational map in which the `ω`-component becomes a polynomial density times an explicit derivation of a simpler equation.
*Expected size.* Symbolic; `F` in Fourier coordinates has all 236 even monomial types in general (it is dense), so the search is for structure, not for a linear system.
*Chance.* Low (5 %); the diagonal `Z/7` in the Klein examples came with a *split* sextic and a monomial map with monomial Jacobian, both absent here.

### E5. First-order deformation of a certificate from a solved neighbour
*Idea.* At `c = 0` (Fermat) a certificate `(D_0, N_0, Q_0)` for `A` exists (coordinatewise-even construction; by `S3`-symmetry it suffices to serve `T_0 = X_0^3Y_0^3`, for which the §2.5 cubic-surface construction with passive ratios applies); the same holds at `Fermat + (5/11)O2` and `Fermat + (5/13)O3` (reflection quotients, E8), which are only two coefficients away from the instance instead of three. Write `F_c = F_0 + cM`, `M = (5/7)O1 + (5/11)O2 + (5/13)O3` (or `M = (5/7)O1 + (5/13)O3` when starting from the `O2`-line), and seek `D = D_0 + cD_1 + …`, `N = N_0 + cN_1 + …`, `Q = Q_0 + cQ_1 + …` of fixed degree. Order `c^1` is a *linear* system in `(D_1, N_1, Q_1)` whose right-hand side is built from `M` and `(D_0,N_0,Q_0)`: `A·2D_0D_1 − D_1 N_0(F_0) − D_0 N_1(F_0) − D_0 N_0(M) + F_0(…)_1 + M(D_0 div N_0 − N_0(D_0)) = F_0² Q_1 + 2F_0 M Q_0`.
*Steps.* (i) Produce an explicit Fermat certificate for `T_0` of minimal degree by fixed-`D` solves over the natural Fermat denominators (`X_0^4 − X_0^2Y_0^2 + Y_0^4` fails: **[verified]** by hand, the `E2` obstruction `X_0^3Y_0 D/3` is not cancelled by the constant Koszul freedom; try the degree-6/8 cofactors and the §2.5 cubic-surface denominators). (ii) Solve the order-1 system; if solvable, iterate to order 2; if the series terminates or the obstruction space is small, reconstruct a polynomial solution in `c` and specialize `c = 1`.
*Size.* Same as a fixed-`D` solve at `deg D_0` (small if a small Fermat certificate exists).
*Chance.* Low (5 %): the solution variety of certificates jumps at Fermat (all reflections present), and the order-1 obstruction is expected to be nonzero; but the computation is cheap and the *dimension of the obstruction* at order 1 is informative.

### E6. The `(01)`-quotient as a CY3 fibration: look for a second fibration or a correspondence with `p_g = q = 0` fibers
*Idea.* `W = {F_s(p,q,a,b,u;t) = 0} ⊂ P(1,1,1,1,2) × A^1_t`. Divisor support of all of `H^4(W)` would serve all of `H^4(X)^{(01)}` (108-dimensional `(3,1)` part: 39 trivial + 69 standard) and hence the entire standard sector by the averaging ideal. The discussion checked one surface fibration (`t, λ = b/a` fixed; `p_g = 3`). Other choices: fix `t` and the ratio `q/p`, or `u/p²`, or use the `τ`-quotient as well (`τ` acts on `W` by `(q,b,t) ↦ (−q,−b,−t)`; the quotient `W/τ` has the invariant coordinates `q², b², t², qb, qt, bt`); or exploit that the generic CY3 fiber is a sextic in `P(1,1,1,1,2)` whose `H^3` has `h^{3,0} = 1`: the `(3,1)` classes of `W` come from `H^1(P^1_t, R^3)`, and divisor support would require the variation of `H^3` of the fibers to be "Tate" over a divisor — not expected. More promising: look for surfaces inside the fibers (e.g. the sextic surface `{u = 0}` of weighted degree 6 in `P(1,1,1,1)` — a smooth sextic surface in `P^3` with `p_g = 10`, no good) or for a Calabi–Yau *threefold* fibration with a section whose relative Picard is large.
*Chance.* Very low as a route to support (the fiberwise `H^{3,0}` is the obstruction); moderate value as a clean statement of why the quotient geometry fails.

### E7. A reflection on a finite cover (holdout-style), adapted
*Idea.* Choose a cover `Z → X` on which some element of a larger group acts as a reflection, with `pd = p`, `pe = gp`, `pr = τp` as in the descent proposition (§4.6), so that uniruled quotients of `Z` plus the (unsupported!) quotient `X/⟨g⟩` combine. The proposition needs `H^k(X)^{⟨g⟩}` supported for a *physical* `g`; for `g` a transposition that is exactly the open `B`-sector problem, so the proposition cannot be applied directly — a genuinely new descent statement for a two-dimensional representation (the "matrix-valued transfer" of `why.md`) would be needed.
*Chance.* Low; this is a theory task, not an experiment.

### E8. Sector theorem on the invariant family, and the solved sub-families (theory + a concrete bridge)
*Idea.* Adapt §4.4 to the 50-dimensional family `V^G` of `G`-invariant even sextics (46 directions modulo the 4 invariant Jacobian directions, i.e. the trivial-isotypic part of `(R_F)_6^+`): the Higgs fields `M_g` for invariant `g` commute with `S3`; the algebra they generate (together with adjoints) acts on each isotypic component; if it acts irreducibly on the standard-isotypic `(3,1)` block (138-dim), then *one* nonzero supported standard-sector class on a family dominating `V^G` would give the whole standard sector generically on `V^G`. This would convert the problem into finding **one** certificate for **any** standard-type target on **some** member dominating the family. (Caveat: §4.4 needs a family *dominating* the parameter space; a construction on a proper sub-family does not qualify.)

*The sub-families, computed **[verified]** (sign symmetries, coordinate parities, diagonal group):*

| member of the family | extra symmetry | status of `A` there |
|---|---|---|
| `Fermat` (`c = 0`) | coordinatewise even | served (cubic-surface construction, §2.5) |
| `Fermat + (5/11)·O2` (`O1 = O3 = 0`) | **separately even in `Y_0, Y_1, Y_2`** (12 terms; diagonal group = the 8 `Y`-sign patterns) | served: `T_0` is odd in `Y_0`, `T_1` odd in `Y_1`, so each summand is served by one reflection quotient (§3.3, finite covering-curve construction; this is the 118/124/107 row of the §4.1 table) |
| `Fermat + (5/13)·O3` (`O1 = O2 = 0`) | separately even in `X_0, X_1, X_2` | served (the `X`-version of §3.3) |
| `Fermat + (5/7)·O1` (`O2 = O3 = 0`) | the three **paired signs** `(X_i,Y_i) ↦ (−X_i,−Y_i)` and `τ` (diagonal group of order 8), but **no order-three diagonal element** (both edge orientations occur in `O1`), so it is *not* in the 26-coefficient `Γ`-family of §4.2 | not served by any theorem: the paired quotient surface is the bidegree-`(3,3)` surface with `p_g = 2` of §4.7; `A` is `P`-invariant with invariant residue (paired signs have determinant `+1`), so it lies in the unserved "transcendental part of the surface fibers" |
| any two orbits together | only `S3 × ⟨τ⟩` | open |

So the "hard" coefficient directions are `O1` (which destroys the coordinate reflections of the `O2`/`O3` sub-families) together with whichever of `O2`, `O3` destroys the paired signs. The most concrete bridge found in this dossier is therefore: **start from `F_0 = Fermat + (5/11)O2` (every odd target served, with an explicit finite construction available in the style of `reflection-coupled-fano`'s modular weighted quartic) and continue in the two coefficients `(5/7, 5/13)` of `O1, O3`**, or start from `Fermat + (5/13)O3` and continue in `(5/7, 5/11)`. Each continuation step is a linear system (first-order deformation, as in E5) once an explicit certificate is in hand at the starting point.
*Chance.* As a theorem: unknown, requires real Hodge-theoretic work; as a guide for experiments: useful (it says exactly which two coefficients carry the difficulty, and provides solved neighbours one coefficient away in two different directions).

### E9. Numerical certificate hunt without symmetry restriction (control experiment)
*Idea.* Same as E3 but with `D` ranging over the full Y-even `E1` kernel at `d = 10` (`≥ 294`-dimensional); the `S3`-invariant restriction costs a 12-fold degree in the worst case, so a non-invariant solution of degree 10–12 is not excluded by the invariant search.
*Size.* `d = 10`: `m ≈ 294 + dim J(F)_10^+`, linear unknowns `≈ 15 × 462/2 + 5796 ≈ 9300`, equations `84903`.
*Chance.* Low (5 %); mainly to make the exclusion statements in E2/E3 symmetric-independent.

### E10. Warm-up on `symmetric-sextic-residual`
*Idea.* The sibling needs a one-dimensional character (`sign ⊗ sign` of `S3 × S3`); the discriminant double cover of the elementary-symmetric quotient carries exactly `sign`, so a *double cover* of the quotient plus the Fano/curve machinery might serve it, and the resulting technique (serving a nontrivial character of a permutation group through a cover) is the missing primitive for the standard representation here (`V_std` becomes `χ_ω ⊕ χ_{ω²}` on the `A3`-cover, cf. §4.7 of the book).
*Chance.* Moderate for the sibling, indirect for this instance.

---

## 7. Feasibility assessment

**Explicit witness.** Any certificate has `deg D ≥ d_min` with `d_min` from §4.1 (at least 10 by the Hilbert-series count alone, exactly determined by the `E1` kernels). An explicit witness at `d ≈ 10–12` would mean a polynomial identity in degree `26–30` with `~10^4–10^5` coefficients; Lean verification of identities of that size is routine in this project (`one-reflection` verified a strict identity with `deg D = 60`, `two-ternary-dense` with `deg D = 55`), using `witness.json` and `Submission.witness`. The difficulty is entirely in *finding* `D`: there is no visible denominator family, and all accepted explicit witnesses came from a structure (reflection, cofactor, contact, cyclic monomial map) that this `F` lacks by design. Realistic probability of an explicit witness from the experiments above: ≈ 10 % in total, dominated by E3/E9 if the minimal degree happens to be small.

**Nonconstructive / finite-algebraic existence.** Every accepted existence proof is *finite-algebraic*: a polynomial curve family on an auxiliary (quotient or cover) variety with an invertible incidence Jacobian, polynomial integration, normalized trace, Piola transport, and legal clearing, formalized through the shared modules `HoldoutCover`, `PolynomialCurves`, `EllipticCompanions`, `ResidueTransport`, `HigherResidues`, `TransitiveCubic`, `AveragingProjectors`, `CertificateOperations`. All of them start from a quotient or cover whose relevant cohomology is **served by a reflection, a paired-sign Jacobian, a diagonal cyclic density, or a Fano cover**. For this instance no such auxiliary variety is known: the invariant quotients miss `V_std` (averaging ideal), the `(01)`-quotient is a CY3 fibration, the three-cycle quotients have age 2 or act freely, and the natural cubic cover's Galois closure does not help (§4.7). A purely cohomological existence argument (Kollár–Larsen uniruledness, surface fibrations with `p_g = q = 0`, Hodge-theoretic generic-target arguments) has **never** been formalized in Lean in this project; the Lean library only supplies the group-algebra combination step (`AveragingProjectors`) once supports are known. So a nonconstructive route would require both new mathematics (a correspondence reaching a two-dimensional representation) and new formalization infrastructure.

**What a Lean proof would require (either route).** `theorem Submission.certificate : HodgeComputer.HasCertificate HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q`, with `F`, `A`, `q = 1` generated from the base branch; a witness route proves `Verify F A 1 d D N Q` (homogeneity of `D` and of the six `N_i` in degree `d+1`, `D ≠ 0`, `IsRelPrime F D` over `ℂ` — typically via `relPrime_of_specialization` on a monic specialization — and the exact identity `E = F²Q`); a finite-algebraic route proves the existential through the shared constructions and must supply its own incidence seed, dominance, source identities, and legality. Allowed axioms: `propext`, `Classical.choice`, `Quot.sound`; budget 60 minutes, 14 GiB.

**Bottom line.** The instance was built to sit exactly outside every mechanism in the book, and the one discussion comment (by the project owner) is a correct reduction plus a correct demonstration that the obvious quotient geometry fails. The rigorous facts available now are the symmetry/representation data (§1), the equivalence `A ⇔ B = T_0+T_1−2T_2` with `(01)`-equivariant data, the `E1` degree bound (§4.1), and the observation (E8) that the target *is* served (by the reflection theorem, with a finite construction) on the two sub-families `Fermat + (5/11)O2` and `Fermat + (5/13)O3`, each two coefficients away from the instance, while the paired-sign sub-family `Fermat + (5/7)O1` is itself unserved (no order-three diagonal, `p_g = 2` quotient surface). The hard part is therefore the orbit `O1` together with one of `O2`, `O3`. A solve is unlikely without a new idea; the experiments E1–E3 and E5 are cheap enough to run and would at least produce publishable negative results (minimal denominator degree, excluded supports, order-1 deformation obstruction) for the discussion thread.
