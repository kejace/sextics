# Accepted solves, group 1: contact46006, coupled-pair-rational-d4-strict, crossed-ternary-asymmetric, cycle-three, cyclic-triple-fano, f3-block-x0-quartic-strict

Catalogue prepared 2026-11-01 from the repository at `/home/kejace/devel/sextics` (branch `main`, HEAD `83fbf5c`), the GitHub discussion threads (read-only via `gh api graphql`), and independent sympy/Singular checks (scripts in `/tmp/hodge/scratch/`). Nothing in the repository was modified.

Conventions. Variables `Z = (X0,X1,X2,Y0,Y1,Y2)` in this order; `τ(X,Y) = (X,−Y)`. A certificate for `(F,q,A)` is `(d, D, N, Q)` with `D ≠ 0` homogeneous of degree `d`, coprime to `F`, `N_i` homogeneous of degree `d+6q−5`, and the cleared identity

    E := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^(q+1)·Q .

Normalized form (sufficient, used by every Lean proof below via `HodgeComputer.normalized_identity`): `E1: q·Σ N_i ∂_iF = A·D`, `E2: D·div N − N(D) = F^q·Q`. "Strict" means `Q = 0`.

Epistemic labels used below: **[proved]** = checked by the accepted Lean proof (kernel replay, axioms `propext`, `Quot.sound`, `Classical.choice`); **[verified]** = recomputed here with sympy 1.12 / Singular 4 (exact rational arithmetic); **[computed]** = my own exact computation not in the record; **[claimed]** = stated in `solution.md`/discussion and not independently rechecked here; **[failed/bounded]** = negative or partial results on record.

---

## 0. Shared infrastructure the six proofs rely on

### 0.1 Verifier interface (`formal/lean/HodgeComputer/Verify.lean`, 85 lines)

```
structure Verify (F A : Poly) (q d : ℕ) (D : Poly) (N : Fin 6 → Poly) (Q : Poly) : Prop where
  denominator_ne_zero     : D ≠ 0
  denominator_homogeneous : D.IsHomogeneous d
  numerator_homogeneous   : ∀ i, (N i).IsHomogeneous (d + 6*q - 5)
  denominator_legal       : IsRelPrime F D
  identity                : defect F A D q N = F^(q+1) * Q
def HasCertificate (F A : Poly) (q : ℕ) : Prop := ∃ d D N Q, Verify F A q d D N Q
```
with `defect F A D q N = A*D^2 − q*D*directional F N + F*divergenceNumerator D N`, `directional P N = Σ pderiv i P * N i`, `divergenceNumerator D N = D*divergence N − directional D N` (`HodgeTestbed/Algebra.lean`). Helper theorems: `normalized_identity` (E1,E2 ⇒ identity), `relPrime_of_specialization` (coprimality in `R[z]` descends from a specialization `R → S` when `f` has unit leading coefficient), `relPrime_of_retract` (transfer along `toUni : Poly → Polynomial Coeff` with a left inverse), `localization_identity` (bridge to the rational field statement).

### 0.2 The two proof architectures that occur in this group

**(A) Explicit-witness pattern** (contact46006, coupled-pair, f3-block; all by `benediamond`, same generator). Files: `Solution.lean` (≈2.7–3.0 KB), `Witness.lean` (17 KB–243 KB), `witness.json`, `proof.json` with `"theorem": "Submission.certificate"` and `"witnessTheorem": "Submission.witness"`.
- `Witness.lean` (namespace `HodgeComputer.Instances.<Id>`; `set_option maxHeartbeats 2000000`, `maxRecDepth 1000` or `4000`): defines `F A D Q N0..N5` as explicit `MvPolynomial (Fin 6) ℂ` sums; 36 private lemmas `pderiv i (X j) = δ_ij`; `e1 : (1:Poly) * directional F N = A*D` and `e2 : divergenceNumerator D N = F^1*Q` proved by `simp only [defs, pderiv_mul, pderiv_pow, derivative_*]` then `apply MvPolynomial.funext; intro x; norm_num <;> ring` (i.e. the identity is checked as a polynomial function on ℂ⁶, legitimate since ℂ is infinite); homogeneity by chained `IsHomogeneous.add/mul/pow` terms; `D ≠ 0` by evaluation at `(1,2,3,4,5,6)`.
- **Coprimality over ℂ**: choose a variable `z` that does *not* occur in `D` and in which `F` has a unit leading coefficient (a pure power `c·z⁶`). `toUni` sends `z ↦ Polynomial.X` and the other five to coefficients `MvPolynomial (Fin 5) ℂ`; `retract : fromUni.comp toUni = id`; `relPrime_of_retract` + `relPrime_of_specialization` with the specialization "evaluate the five coefficient variables at a point where `D ≠ 0`": the fibre of `D` is a nonzero constant, the Bézout pair is `(0, 1/const)`. This trick works whenever `D` omits a variable in which `F` is monic up to a unit; it never needs a multivariate resultant.
- `Solution.lean` (`maxHeartbeats 8000000`, `maxRecDepth 8000`) identifies the verifier-generated `HodgeSubmission.F/A/q` and `HodgeWitness.d/D/N/Q` (parsed from `instance.rec`/`witness.json`) with the `Witness.lean` definitions by `MvPolynomial.funext; intro x; norm_num [...] <;> ring` (six `fin_cases` for `N`), then `simpa ... using <Id>.certificate`.

**(B) Nonconstructive pattern** (crossed-ternary-asymmetric, cycle-three, cyclic-triple-fano). No `witness.json`; `proof.json` has only `"theorem": "Submission.certificate"`. All three end in the shared theorem

```
HoldoutConstruction.hasCertificate_of_normalized
  (F A D : Poly) (N : Fin 6 → Poly) (Q : Poly)
  (hF : F.IsHomogeneous 6) (hA : A.IsHomogeneous 6) (hprime : Prime F)
  (hlegal : IsRelPrime F D)
  (normal : directional F N = A*D) (divergence : divergenceNumerator D N = F*Q)
  (i : Fin 6) (hchart : ¬ F ∣ pderiv i F) : HasCertificate F A 1
```
(`Constructions/HoldoutCover/HomogeneousCertificate.lean`). Internally the denominator is first multiplied by `∂_iF·(D+F·U)` to get a legal initial homogeneous component and then homogeneous extraction produces the final denominator of degree `(d'+2)·d'` where `d'` is the order of that auxiliary denominator — **no explicit degree bound is ever extracted**, and `Q ≠ 0` is allowed and in general used.
- **Coprimality over ℂ** in pattern (B) is always: `F` irreducible ⇒ `Prime F` ⇒ `IsRelPrime F D ↔ ¬ F ∣ D`, and `¬ F ∣ D` comes from an evaluation `ρ` into a function field with `ρ F = 0`, `ρ D ≠ 0`. Irreducibility of `F` is proved uniformly by: write `F` as a monic degree-6 polynomial in one coordinate over `ℂ[five others]` (`MvPolynomial.finSuccEquiv`, possibly after a `rename`), specialize the five coefficients to a one-parameter family in which the polynomial becomes `X⁶ + g(T)` with `g(0)=0`, `g'(0)≠0`, and apply `HoldoutConstruction.irreducible_power_add_constant` (Eisenstein at the prime `(T)`), then `Polynomial.Monic.irreducible_of_irreducible_map`.
- The chart hypothesis `¬ F ∣ ∂_iF` is proved by an axis specialization (`F ↦ c·X⁶`, `∂_iF ↦ 6c·X⁵`) or by degree.

### 0.3 Shared library modules used in this group
`HodgeComputer.Constructions.HoldoutCover.*` (fqlx, PR #128/#134: polynomial families, generic incidence, polynomial normals, finite derivations/trace, residue clearing, homogeneous certificate), `PolynomialCurves.{ArithmeticLifting, NormalizedCertificate}` (fqlx PR #150, hodgepodgeconjecture PR #151), and the three libraries *extracted from* the solves below: `EllipticCompanions` (from crossed-ternary-asymmetric), `ResidueTransport` (from cycle-three), `TransitiveCubic` (fqlx, mixed-biweight, extending cycle-three's method).

Policy limits relevant to the big proofs (`.github/solve-policy.json`): 64 files, 16 MiB source, 60-minute proof budget, 14 GiB, Lean 4.33.0, pinned mathlib `db584cd6`.

---

## 1. contact46006 — "Contact cone 46006"

**Record.** `q = 1`. `F`: 56 monomials (coefficients up to 201480), `A`: 51 monomials. Discussion #91 (no comments). PR #28 by `benediamond`, accepted 2026-10-15. Files: `Solution.lean` 2.7 KB, `Witness.lean` 243 KB (951 lines), `witness.json` 11 KB.

**Structure of (F,A)** [verified]. Let `z=(z0,z1,z2)`,
```
h = z0³ + 2 z2³,
C = z1³ + 5 z1² z2 − 22 z1 z2² − 56 z2³ = (z1 − 4z2)(z1 + 2z2)(z1 + 7z2)   (three concurrent lines through [1:0:0]),
k = z1³ − 10 z1² z2 + z1 z2² + 120 z2³      (recovered: k = (F(X,0)/2 − h²)/C, exact division),
f = h² + C·k  (plane sextic with square restriction to each of the three lines),
F = f(U) + f(V),  U = X+Y, V = X−Y   (polarization; τ = block swap U↔V).
```
`W := ∇C × ∇h = (6 z2² C_{z1}, 3 z0² C_{z2}, −3 z0² C_{z1})` satisfies `div W = 0`, `W(C) = 0`, `W(h)=0`, `W(f) = C·W(k)`; the plane target is `A_C := W(k) = 3 z0²·{C,k}_{(z1,z2)} = 135 z0²z1⁴ − 414 z0²z1³z2 − 2817 z0²z1²z2² − 720 z0²z1z2³ + 23256 z0²z2⁴`, and `A = A_C(U) − A_C(V)` exactly.

**Method.** Perfect-power contact field (book Chapter 2 §2.3 with `m = 2`, `H = h`, `D = C`, `K = k`, `B` constant), applied to each polarized block and subtracted. In `(U,V)` coordinates the field is `v = (W(U)/C(U), −W(V)/C(V))`: a difference of two independent strict block certificates (`A_C(U)` over `C(U)`, `−A_C(V)` over `C(V)`), combined by the product-denominator rule of `CertificateOperations`. **Constructive, strict.**

**Witness** [verified: E1, E2 and the cleared identity hold exactly; Q = 0].
- `d = 6`, `D = C(U)·C(V)`, 44 monomials, `D ∈ ℚ[X1,X2,Y1,Y2]` (no `X0, Y0`). Over ℚ already a **product of six hyperplanes** [verified, sympy `factor` and Singular `factorize`]:
  `D = ∏_{c∈{−4,2,7}} (X1 + cX2 − (Y1 + cY2))·(X1 + cX2 + (Y1 + cY2))`.
  So the support is six hyperplanes = pullbacks of the three contact lines of the plane sextic under the two polarized projections.
- `N_X = (W(U)C(V) − W(V)C(U))/2`, `N_Y = (W(U)C(V) + W(V)C(U))/2` [verified componentwise], degree 7, 48/84/84/48/84/84 monomials; parities correct (X-components Y-odd, Y-components Y-even).
- `N(F) = A·D`, `div N = 0`, `N(D) = 0` ⇒ strict.

**Lean proof.** Pattern (A). Coprimality: univariate variable `X0` (`toUni (X 0) = Polynomial.X`), `fUni.leadingCoeff = C 2`; specialization at `(X1,X2,Y0,Y1,Y2) = (1,2,3,4,5)` gives `dFiber = C (−45874512)` (constant, since `D` is free of `X0`), Bézout `(0, −1/45874512)`. `Witness.lean`: `maxHeartbeats 2000000`, `maxRecDepth 4000`. Tactics: `simp only` + `MvPolynomial.funext` + `norm_num <;> ring` for `e1/e2`; no `decide`, no `linear_combination`.

**Symmetry data** [computed]. Projective diagonal group (diagonal maps preserving `F` up to scalar, mod scalars) has order **6** (Smith form of the exponent-difference lattice `[1,1,1,1,6]`); no coordinate reflection preserves `F`; the residue `A·Ω/F²` is τ-invariant but its character on the full diagonal group has order 6.

**Scope / caveats.** The record serves exactly the "difference" target `A_C(U) − A_C(V)` attached to the contact presentation; nothing is claimed for other targets on this polarized sextic. Smoothness of `F` is a standing assumption of the challenge, not re-verified here.

**Related instances.** Solved siblings by the same contact-ladder mechanism: `polarized-contact-coupled-d2` (Hamiltonian field with two first integrals, `D = X2² − Y2²`, coupling `1000003·r(U)r(V)` annihilated by the fields), `polarized-contact-deformation-d2`, `cone-ladder-linear` (`f = e³ + P·k`, hyperplane `D = X1 − X2`), `family-passive-contact`. **Open**: `polarized-k3-remainder` — same polarized-swap structure `F = (f(U)+f(V))/2` but with `f` = the split-ternary plane sextic (`X0⁶ + 5X0⁵X1/7 + 5X0²X1³X2/7 + 15X0X2⁵/7 + X1⁶ + 10X1⁵X2/7 + X2⁶`, i.e. the X-block of `split-ternary-q2`/`rank-one-biquartic`) [verified], target `X0X1X2Y0Y1Y2`. Contact fields of the contact46006 type only serve targets of the form `B(U) − B(V)` with `B = W(k)` (the "joining" sector); the K3-remainder target `(B(U) − B(V))` with `B = U0²U1²U2²/64` is shown in its `why.md` to survive into `Λ²H²(w² = f)_prim`, which no contact presentation of `f` has been found to reach. No same-variety open instance (`sameVarietyInstances = []`).

---

## 2. coupled-pair-rational-d4-strict — "Quartic cofactor · coupled deformation"

**Record.** `q = 1`, `F = X0⁶+X1⁶+X2⁶+Y0⁶+Y1⁶+Y2⁶ + Y1²·(X0⁴ − X0²Y0² + Y0⁴)` (9 monomials), `A = X0³X1X2Y0 + X0X1X2Y0Y1²/3`. Discussion #77 (no comments). PR #15 by `benediamond`, accepted 2026-10-15. `Witness.lean` 18 KB.

**Method.** Two-variable Fermat-cofactor logarithmic field with passive multiplier, transported along a deformation of `F` by a multiple of its own denominator. **Constructive, strict.**

**Mechanism** [verified]. With `x = X0`, `y = Y0`, `B = X1X2`: `x⁶ + y⁶ = (x²+y²)·D`, `D = x⁴ − x²y² + y⁴` (the 12th cyclotomic form; irreducible over ℚ, splits over ℚ(ζ₁₂) into the four hyperplanes `X0 = ζY0`, `ζ` a primitive 12th root of unity). Field `N0 = B·y(x²−y²)/6`, `N3 = B·x³/6`, others 0. Then `N(F0) = A0·D` with `A0 = x³yB`, `div N = xyB/3`, `N(D) = D·div N` (logarithmic). Deformation `F = F0 + D·W`, `W = Y1²`, with `N(W) = 0` ⇒ `N(F) = D·(A0 + W·div N) = D·A`; E2 is unchanged. General principle recorded in `solution.md`: any `F0 + D·W` with `N(W)=0` is served with target shifted by `W·div N`.

**Witness** [verified]. `d = 4`, `D = X0⁴ − X0²Y0² + Y0⁴`, `N = (X1X2Y0(X0²−Y0²)/6, 0, 0, X0³X1X2/6, 0, 0)`, degree 5, `Q = 0`; `gcd_ℚ(F,D) = 1`.

**Lean proof.** Pattern (A). Coprimality: univariate variable `X1` (`F` monic in `X1`, `D` free of `X1`), specialization at all-ones gives `dFiber = C 1`, Bézout `(0,1)`. `maxHeartbeats 2000000`, `maxRecDepth 1000`.

**Symmetry data** [computed]. All six coordinate reflections preserve `F` (every exponent even); projective diagonal group order **864** (`[2,2,6,6,6]`). The residue of `A·Ω/F²` is invariant under the reflections of `X0, X1, X2, Y0` (A odd in each) and anti-invariant under those of `Y1, Y2`. Hence nonconstructive existence was already available from the reflection theorem (Chapter 3 §3.3); the solve supplies the explicit quartic.

**Related.** Undeformed control `fermat-pair-rational-d4-strict` (solved; same field with roles `(x,y,B) = (X0,X1,X2Y0)`), `fermat-q2` (same Fermat variety, `q=2`, solved). No same-variety open instance. Method-related open targets: `dense-rational` (a near-Fermat deformation with all 236 even monomials; the deformation principle above needs the perturbation to be `D·W` with `N(W) = 0` for a chosen Fermat-cofactor `D`, which a generic perturbation is not — the `why.md` warns that formal deformation of a field does not produce a homogeneous rational field); `rank-one-biquartic` (`F = f(X)+g(Y) + (5/7)L(X)²M(Y)` with non-Fermat split blocks, so no cyclotomic cofactor is available).

---

## 3. f3-block-x0-quartic-strict — "F₃ · quartic contact"

**Record.** `q = 1`, `F = X0⁶ + X0²X1²Y0² + X1⁶ + X2⁶ + Y0⁶ + Y1⁶ + Y2⁶` (7 monomials), `A = X0³Y0³ − (5/3)X0X1⁴Y0` (a bundle direction in the 2-dimensional odd-coordinate character space). Discussion #72 (no comments). PR #10 by `benediamond`, accepted 2026-10-15 (earliest of the six). `Witness.lean` 17 KB.

**Method.** Contact conic in the square-coordinate plane of a Hesse cubic; logarithmic field. **Constructive, strict.**

**Mechanism** [verified]. In `(u,x,w) = (X0², X1², Y0²)` the active cubic is `u³ + x³ + w³ + uxw` (Hesse pencil parameter `t = 1`), plus passive `X2⁶+Y1⁶+Y2⁶`. On the conic `u² = xw`, parametrized `(u,x,w) = (sr, s², r²)`, the cubic restricts to `(s³+r³)²` — a perfect square (contact). Pulling the conic back gives `D = X0⁴ − X1²Y0² = (X0² − X1Y0)(X0² + X1Y0)` (two rank-3 quadric cones). Field: `N0 = Y0(X0²Y0² − X1⁴)/6`, `N1 = X0X1Y0³/3`, `N3 = −X0³X1²/3`; `N(F) = A·D`, `div N = 2X0Y0³/3`, `N(D) = D·div N`. The ratio `−5/3` is forced by the contact construction.

**Witness** [verified]. `d = 4`, `N` degrees 5, `Q = 0`, `gcd_ℚ(F,D)=1`.

**Lean proof.** Pattern (A). Coprimality: univariate `X2` (`F` monic in `X2`; `D ∈ ℚ[X0,X1,Y0]`), specialization `(X0,X1,Y0,Y1,Y2) = (1,2,3,4,5)` gives `dFiber = C (−35)`, Bézout `(0, −1/35)`. `maxHeartbeats 2000000`, `maxRecDepth 1000`.

**Symmetry data** [computed]. All six reflections preserve `F`; projective diagonal group order **2592** (`[2,6,6,6,6]`). Residue invariant under the `X0` and `Y0` reflections (A odd in both), anti-invariant under `X1, X2, Y1, Y2`. Reflection theorem applies; the solve is the explicit effective version.

**Related.** Solved Hesse-type siblings: `family-hesse-web` (`t = −2` fiber of the same pencil with passive weights `2,3,5`, target `X0³Y0³`, `D` = degree-6 norm `a³ + t a² x + 4x³` over the cubic algebra `C³ − tC² − 4 = 0` tracing three contact channels), `hesse-radial-…`, `fixed-witness-hesse-…`, `twoties-*`, `klein-hesse`, `hesse-triangle`, `collision-h3`, `m3-coupled`. No same-variety open instance. Method-related open targets: `quadric-contact-residual` (quadratic coordinates `S = 2(X0²+Y0²)`, `T = 4X0Y0` reduce to a cubic surface, but `X0³Y0³` has the wrong `X0↔Y0` exchange character for that density), `eight-plane-configuration` (`F ∈ (X_i² + Y_i²)`, a contact arrangement of eight planes; planes have type (2,2) so cannot serve the (3,1) target by themselves).

---

## 4. crossed-ternary-asymmetric — "Crossed ternaries · elliptic transfer"

**Record.** `q = 1`, `F = P(X0,X1,Y0) + P(Y1,Y2,X2)` with `P(u,v,w) = u⁶+v⁶+w⁶ + (5/7)(u⁵v + u³vw² + uvw⁴ + 2uv⁵)` (14 monomials), `A = X0⁴X1Y1` (block degrees (5,1)). Discussion #114: one thread by `1337isnot1337` (F₅ seed → conditional contraction → PR). PR #149 by `1337isnot1337` (OpenAI Codex-assisted), accepted 2026-10-20. 14 Lean files, 9039 lines, 424 KB; `checks.json`: smooth reduction mod 5, target nonzero in the reduced Jacobian ring, coordinate reflections at indices 2 (`X2`) and 3 (`Y0`) serving nothing, Hessian symmetrizer dimension 2, projective diagonal group order 24.

**Why it is not a reflection case** (`why.md`). With `r` = negate `Y0`, `s` = negate `X2`, `ρ` = negate the second block, `τ = ρ·r·s`. `A` is even in `Y0` and `X2`, so the residue is anti-invariant under `r` and `s` and lies in the `ρ`-invariant ("split") sector; existence was known nonconstructively from Chapter 3A §3A.3 (involution-stable ternary decompositions), but the elliptic contraction had to be made algebraic.

**Method.** Elliptic-companion / moving-conic construction with an algebraic (rational-residue) contraction of the elliptic factors; compatible-fraction clearing. **Nonconstructive** (classical choice, finite extensions, trace); denominator unknown, `Q ≠ 0` in general.

**Mechanism** (from `solution.md`; the Lean statements are [proved], the formulas quoted are [claimed] unless marked).
1. Companion cone `P(x,a,y) + S³ + T³ = 0` in weights `(1,1,1,2,2)`; integral form `H = 7(x⁶+a⁶+y⁶+S³+T³) + 5(x⁵a + x³ay² + xay⁴ + 2xa⁵)`. Moving conic: `x = 4t`, `a, y` quadratic, `S, T` quartic in `t`, constant terms `k0,k1,k2` free, 13 unknowns `c_0..c_12`, 13 incidence equations (coefficients of `H` in degrees 0..12). Seed over 𝔽₅: `(k0,k1,k2) = (4,2,3)`, `c = (2,0,4,0,2,0,1,0,3,0,1,0,2)`, i.e. `x=4t, a=2−t², y=4+2t², S=2+t²+3t⁴, T=3+t²+2t⁴` (from the line `Z2 = 2Z0+4Z1, s = 3Z0+3Z1, t = 4Z0+2Z1` on the Fermat cubic threefold pulled back by `Z_i = z_i²`, per the discussion). Explicit right inverses of the 12×14 freeness and 13-column incidence matrices over 𝔽₅ (`ConicSeed.lean`, `decide`). Lifting + polynomial normal + dominance (four independent coordinates) + finite trace give `CrossedFixedConic.companion_residue_solution`: a companion residue primitive for *every* polynomial source.
2. Radial descent to `P(x,a,1) + S³ + T³ = 0`; sources `7x⁴a` and `7xST`.
3. Fermat elliptic curve `a_e³ + b_e³ + 1 = 0`, presented as `y_e² = x_e³ − 432`, `a_e = −(y_e+36)/(6x_e)`, `b_e = (y_e−36)/(6x_e)`, invariant derivation `δa_e = b_e²`, `δb_e = −a_e²`. Substituting `S = u a_e`, `T = u b_e` with `u³ = p := P(x,a,1)` converts companion data into scalar equations `δf + div_{x,a}C = α`, `α = x⁴a/u⁴`, and (after the coboundary `δ(b_e²/a_e) = b_e/a_e² − a_eb_e`) `δg + div_{r,z}G = β·θ`, `β = r/v²`, `v³ = q := P(r,z,1)`, `θ = b_e/a_e²`.
4. Contraction: translate the first primitive by a generic elliptic point; residue functional `R` over finite fibres (excluding the three theta fibres `r³ = 1728`) with `R(δh) = 0`, `R(∂_ih) = ∂_iR(h)`, `R(g) = 0` (pole separation proved, not assumed). Theta kernels `K0 = −36x(x³+3456)/(x³−1728)²`, `K1 = 3888x⁴/(x³−1728)²`; local formula `Res_r(h0K0 + h1K1) = −12h0'(r) + 432h1'(r) + (r²/2)h1(r) = −δh(r,−36)`; infinity correction via a finite rational jet approximation at `1/x_e = 0` (denominators `5+2n`). Result `R(fθ) = 3α − div_{x,a}ξ`, then `div_{x,a}U + div_{r,z}V = αβ` with `U_i = (R(C_ig) + βξ_i)/3`, `V_j = R(fG_j)/3`.
5. Multiply by `u⁴v²`, descend by normalized trace through the cube-root extension: `W` over `ℂ(x,a,r,z)` with `Σ_i(∂_iW_i − (4∂_ip/3p)W_i − (2∂_iq/3q)W_i) = x⁴ar`.
6. Physical chart `X0 = tx, X1 = ta, X2 = tλ, Y0 = t, Y1 = tλr, Y2 = tλz` (`F = t⁶(p + λ⁶q)`, `A = t⁶λ·x⁴ar`), chart vector `w_i = −λW_i/(3p)`, `w_λ = (3ps + P_d + λ⁶Q_d)/(18pλ⁴q)`; `𝒱(F) = A`, `div 𝒱 = F·(−(3ps+P_d)/(9p²qλ⁵t⁶))`.
7. Clearing without introducing a factor of `F`: every fraction is tracked simultaneously under `θ` (injective, into `FractionRing Poly` via a monomial map) and `ρ` (into `RatFunc` of the splitting field of `X⁶ + p/q`, where `ρF = 0`); `normalized_data_of_compatible` yields `D ≠ 0`, `IsRelPrime F D`, E1, E2; then `hasCertificate_of_normalized` with chart `Y0` (index 3).

**Coprimality over ℂ** [proved]. `physicalF_irreducible`: after `normalPermutation` + `finSuccEquiv`, `F` is monic of degree 6 in `Y0`; specialization `X0 = X1 = Y1 = 0`, `X2 = T−1`, `Y2 = i` gives `Y0⁶ + (T−1)⁶ − 1`, Eisenstein at `T` (`irreducible_power_add_constant`). `physical_normal_legal : ¬ F ∣ ∂_{Y0}F` by degrees.

**Lean architecture.** `Solution.lean` (663 B, no `set_option`): `hF`, `hA` by `unfold; simp; ring`, then `exact CrossedPhysicalChart.physical_certificate`. Chain: `CertificateClosure.physical_certificate` ← `normalized_data` ← `CrossedCompatibleFractions.normalized_data_of_compatible` + `CrossedProductField.weighted_divergence` (`ProductDivergence.lean`) ← `ordinary_product_divergence` ← `CrossedEllipticContraction.contraction` (`EllipticContraction.lean`) + `CrossedFanoElliptic.first_scalar_primitive / second_scalar_primitive` (`EllipticKernels.lean`) + `scalar_primitive_of_data` (`ScalarCalculus.lean`) + `CrossedFanoAffineData.firstSource_data / secondSource_data` (`CompanionPrimitives.lean`) ← `CrossedFixedConic.companion_residue_solution` (`MovingConic.lean`) ← `ConicSeed.lean`; `PhysicalTransport.lean` (chart derivations, irreducibility, physical identities); `RadialResidues.lean` (rational residues, radial descent); `EllipticTranslation.lean`, `InfinityPairing.lean`. Shared imports: 16 distinct `HoldoutCover.*` modules (`AlgebraicFields, EvaluationMatrix, FiniteDifferentials, FiniteTrace, GenericIncidence, HomogeneousCertificate, HomogeneousField, HypersurfaceDomain, Piola, PolynomialCalculus, PolynomialFamilyBridge, PolynomialFieldMap, PolynomialNormal, ResidueClearing, ResidueDerivatives, ResiduePiolaPullback`). Heartbeats: sectionwise, mostly `100000–800000` with `backward.isDefEq.respectTransparency false`; `ConicSeed` sections `3000000–8000000` with `maxRecDepth 10000–100000`; `CompanionFields` up to `3000000`; `MovingConic` has sections at `20000` with `profiler true`. Tactic census (14 files): `simp` 529, `ring` 152, `field_simp` 83, `norm_num` 60, `linear_combination` 59, `decide` 44, no `native_decide`/`grind`/`polyrith`.

**Symmetry data** [computed; consistent with `checks.json`]. Projective diagonal group order **24** (`[2,2,6]`); reflections `X2`, `Y0`; residue anti-invariant under each reflection, invariant under `ρ`; character order 6 on the full group.

**Failed / bounded attempts** [failed/bounded]. `solution.md` records earlier bounded searches "concerning specified denominator ansatzes, degree limits, or finite numerical trials" in discussion #114 — none used as a nonexistence claim and none visible in the thread as fetched (the thread has a single comment with two replies, all by the solver).

**Extracted library.** `HodgeComputer.Constructions.EllipticCompanions` (entry point `Contraction` / `CrossedElliptic.product_divergence`; `ScalarEmbedding`, `AffineCompanion`), explicitly stated *not* to extend pole order automatically.

**Related.** `split-ternary-q2` (solved by `fqlx`, 3+3 block-sign split at `q = 2`, `HigherResidues` library). Open: `polarized-k3-remainder` (swap instead of block sign; 3A.4–3A.6 K3 remainder), `four-plus-two-residual` (4+2 split: triple cover of the 4-variable core has trivial canonical class, no Fano companion), `rank-one-biquartic` (`f(X)+g(Y) + (5/7)L(X)²M(Y)`: the branch locus of the cubic cover now varies with both blocks, so the fixed elliptic contraction above does not transfer; the split blocks are the same plane sextics as in `split-ternary-q2`).

---

## 5. cycle-three — "Three couplings · cyclic"

**Record.** `q = 1`, `F = Σ(X_i⁶ + Y_i⁶) + ½(X0³X1Y0Y1 + X1³X2Y1Y2 + X2³X0Y2Y0)` (9 monomials), `A = X0³Y0³`. Discussion #55. PR #165 by `hodgepodgeconjecture` (substantial OpenAI Codex assistance), merged 2026-10-20, accepted 2026-10-21. **60 Lean files, 14625 lines, 2.14 MB** (`CubicIdentityData.lean` 1.59 MB / 8539 lines of expression data; `CubicReflectionData.lean` 193 KB; `CubicCovariants.lean` 71 KB). `sameVarietyInstances = [cycle-q2]`.

**Why the holdout cover fails here** (`why.md`) [claimed]. With `T⁶ = Y0⁶ + Y2⁶` the closing term becomes `½αβX0X2³T²` with no matching `y²` term, so the `y ↔ T` exchange of the holdout Fermat-cover construction is lost; the §2.8 complementary quotient (weights `(5,1,5)` on `Y`) is lost too (closing monomial has weight 4).

**Method.** Paired-sign cubic reduction + classical covariants + Jacobian residue transport with a moving Weierstrass map ("cyclic cubic / Jacobian residue primitive", book §2.9). **Nonconstructive**; no denominator degree or witness.

**Mechanism** (Lean-proved statements; formulas as in `solution.md`).
1. Chart: `p_i` = images of the coordinates in the function field of `F`; `R = p2²`, `v = p3/p0`, `c1 = (p4/p1)/v`, `c2 = (p5/p2)/v`, `x = p0²/R`, `y = p1²/R`. Then `F = R³·U(x,y,1)` with `U = ax³+by³+cz³+dx²y+ey²z+fz²x`, `(a,b,c) = (1+v⁶, 1+c1⁶v⁶, 1+c2⁶v⁶)`, `(d,e,f) = ½v²(c1, c1c2, c2)` (`CyclicCubic.paired_sign_reduction` [proved as a ring identity]).
2. Integral covariants `H, Θ, J` of the six-coefficient cubic; `s4 = ae²f + bdf² + cd²e`, `s6 = 27a²b²c² + 4a²ce³ + 4ab²f³ + 6abcdef + 4bc²d³ − d²e²f²`; syzygy on `U = 0`: `J² = 4Θ³ + 4s4ΘH⁴ − s6H⁶` (`CubicSyzygy.full_syzygy`, universal polynomial identity including the `U`-divisible term). `ξ = Θ/H²`, `η = −J/(2H³)` satisfy `η² = ξ³ + s4ξ − s6/4`; tangent identity `D_yξ/(2η) = 3/g`, `g = U_x` (`CubicFieldDifferential.invariant_differential`). Nonvanishing of `g, H, J, η, D_yξ` proved via an explicit point (`CubicPoint.lean`).
3. Jacobian model: `T = Rv⁶`, `X = R²ξ`, `Y = R³η`; `V = Y² − X³ − a8·X/8 − b256/256 = 0`, with explicit homogeneous `a8(c1,c2,α,r,t)` (degree 4 in `(r,t)`), `b256` (degree 6) in `JacobianModel.lean` (`a8 = 8R⁴s4`, `b256 = −64R⁶s6` on the physical values). Every polynomial source on `V` has a rational residue primitive (`JacobianResidue.family_residue_solution`, `JacobianClearing.polynomial_clearing`): polynomial family of degrees `(2,4,6)` in `(R,X,Y)` with parameter `T`, endpoint parameters `R(0), X(0)`, 13 unknown coefficients; 𝔽₅ seed at `(c1,c2,R(0),X(0)) = (2,3,1,3)`: `R = 1+2T+4T²`, `X = 3+T+3T²+4T³+4T⁴`, `Y = 2+3T+3T²+4T³+4T⁴+T⁵+T⁶` with explicit inverse coefficient Jacobian (`JacobianSeed.inverse_product` by `decide`); shared `ArithmeticLifting`, polynomial normal (`JacobianNormal`), evaluation determinant ≠ 0, polynomial-family integration and normalized trace.
4. First variation: replace `a = 1+v⁶` by `α+v⁶`, differentiate at `α = 1` using dual numbers `K[ε]/(ε²)` (`TrivSqZeroExt`, `CubicJetCalculus.lean`): `ġ/g = (ȧ_Wξ + ḃ_W)/(2η²) − (g/3)·D_y(ξ̇/(2η))`, `ġ = 3x² − x³U_xx/g`. Jacobian source `P_J = −V_α/(144R) = (R⁵/144)(ȧ_Wξ + ḃ_W)`.
5. Transport: tangent derivations `δ_i = ∂_{i+1} − (F_{i+1}/F_0)∂_0`; coordinate matrix `M_{ji} = δ_i(T,c1,c2,R,X)_j` invertible with residue-volume multiplier `F_0·det M/(2Y) = −144Rv³` (`PhysicalJacobianCoordinates.density_multiplier`); dominance by `CoordinateIndependence`; the Jacobian function field embeds (`PhysicalJacobianMap`). Piola pullback (`CycleResidueTransport.piola_pullback`) gives physical residue divergence `−v³(ȧ_Wξ + ḃ_W)/(2η²)`. Moving-map correction: `W = (λhp1, 0, −λh1p3, λhp4, 0)` with `h = p0F0 + p3F3`, `h1 = p1F1 + p4F4`, `λ = p3³/(4p0⁵p1²p2²)` has zero residue divergence and acts as `v³gD_y`; add `ξ̇W/(6η)`. Normal-chart correction `C = −p0³p3⁴/h` in the `p3` slot. The sum has residue divergence `−A_0/F_0 + A·F_00/F_0²` (`PhysicalPrimitive.physical_residue_solution`, closed by `linear_combination ratio^3*source_transport`).
6. `PhysicalClearing.certificate_of_residue_solution`: shared `residue_solution_clearing` (denominator `D = ∂_0F·d` with `ρ(d) ≠ 0`) then `certificate_of_normalized` = `hasCertificate_of_normalized` with chart `X0` (`normal_legal` by the axis specialization `F ↦ X⁶`, `∂_0F ↦ 6X⁵`).

**Coprimality over ℂ** [proved]. `PhysicalField.equation_irreducible`: `F` monic in `X0` over `ℂ[X1,X2,Y0,Y1,Y2]`; `holdoutSpecialization Complex.I` gives `X⁶ + ((X+1)⁶ − 1)`, Eisenstein; `sextic_irreducible` via `finSuccEquiv`. Legality then from `residue_solution_clearing` (`IsRelPrime f D` in its conclusion).

**The big identity** (how it is discharged). `CubicSyzygy.data_identity` is split into **190 coefficient identities** (`CubicCoefficientsA–H.lean`, `coefficient_18_0 … coefficient_0_0`), each proved by
`apply identity_of_encoding_context v Left_i_j Right_i_j (grade) (2^42); all_goals decide +kernel`
with `maxRecDepth 131072`, `maxHeartbeats 16000000`. The soundness chain is fully proved: expressions are `Lean.Grind.CommRing.Expr` trees (`CubicReflectionData`, `SignedExpressions.ringEval_denote`); each coefficient has a fixed diagonal weight (`CubicGradeCheck.grade_sound`); at fixed weight the exponents of `a,b,c` are determined by those of `d,e,f` (< 19), so the Kronecker code `m3 + 19m4 + 361m5` is injective (`CubicKronecker`); signed expressions split into two ℕ-polynomials with a proved coefficient-mass bound, so evaluation at `2⁴²` is injective (`PolynomialDigitCheck.eval_injective_bounded`). Products of the bivariate data are `rfl` (`CubicIdentityProducts`). No `native_decide`; this is kernel `decide` on numerals. Tactic census (60 files): `decide` 249, `rfl` 221, `simp` 285, `ring` 82, `norm_num` 64, `linear_combination` 29, `field_simp` 28, `nlinarith` 1. Other heartbeat settings: `JacobianSeed/Incidence/SeedBinding` `8000000` with `maxRecDepth 100000`; `PhysicalJacobian*`, `PhysicalCubicVolume` `3000000`; most physical files `2000000`; `backward.isDefEq.respectTransparency false` nearly everywhere.

**Lean architecture (top).** `Solution.lean` (588 B): `hF` by `unfold; ring`, `exact CyclePhysicalPrimitive.certificate`. Shared imports: `HoldoutCover.{AlgebraicFields, EvaluationMatrix, FiniteDifferentials, FiniteTrace, HomogeneousCertificate, HypersurfaceDomain, PolynomialCalculus, PolynomialFamilyBridge, PolynomialFieldMap, PolynomialNormal, ResidueClearing, ResidueDerivatives, ResiduePiolaPullback}`, `PolynomialCurves.ArithmeticLifting`.

**Symmetry data** [computed]. Paired signs `(ℤ/2)³`, cyclic ℤ/3 permutation of the pairs (not diagonal), diagonal ℤ/6 with weights `(2,1,0,3,2,1)`. Full projective diagonal group order **72** (`[2,6,6]`); no coordinate reflections. The residue of `A·Ω/F²` is invariant under the paired signs and under the weights-`(2,1,0,3,2,1)` action, but has **order 3** on the full diagonal group: e.g. `γ = (ζ3, ζ3, 1, ζ3², 1, 1)` preserves `F` (check: couplings have weights 6, 3, 3 ≡ 0 mod 3) and multiplies the residue by `ζ3^{9+4} = ζ3`. (This is not used or contradicted by the proof, which only needs the paired signs; it does mean the target is outside the invariant sector of the full diagonal group.)

**Failed / bounded attempts on record** [failed/bounded].
- `rbajaj5` (#55, 2026-10-17): exact-rational exclusion of every denominator of degree ≤ 5 (kernel of `D ↦ [AD]` in `S/J` is zero for `d ≤ 4`, equals `J_5` at `d = 5`; Hessian map rank 21 forces `a = 0`), and of degree 6 (`J_6` of dimension 36; after the Euler/`F` gauge 35 parameters; full second-equation correction matrix 2352 columns of rank 1451; seven quadrics in five parameters with a full-rank degree-four multiplier matrix ⇒ all parameters vanish; "the only pre-gauge denominator is a multiple of `F`"). Full numerator spaces and unrestricted `Q`. Hence the accepted (unknown) denominator has degree ≥ 7.
- `hodgepodgeconjecture` (2026-10-20): Gauss–Manin pilot on the slice `a_0 = t, a_1 = a_2 = 1` (rank-two connection over `ℚ(t)`); a bounded search for a rational `h = (P_0/D, P_1/D)` solving `h' + Ch = b` (degree-42 data) — partial, superseded by the Jacobian route; a surface-quotient proposal (`ℙ²_u × ℙ¹_{(b:r)}` surface, generically smooth at `(c1,c2) = (2,3)` mod 31) also superseded.

**Extracted library.** `HodgeComputer.Constructions.ResidueTransport` (`Density`, `NormalChange`, `MovingMap`, `FirstJets`, `CoordinateIndependence`); book §2.9 "Coefficient and source scope" generalizes the Jacobian model to the 26-monomial equation space with the same symmetries and lists 22 odd degree-6 monomial targets with invariant residue that get polynomial Jacobian sources.

**Related.** **Same variety, open:** `cycle-q2` (`q = 2`, `A = X0²X1²X2³Y0²Y1²Y2`, residue character `ζ6²` under the weights-`(2,1,0,3,2,1)` action, hence erased by the invariant-sector/paired-sign construction; discussion #61 as of 2026-10-27: `benediamond`'s reduction to a rank-two de Rham/connection problem and a proper-K3 rational-current formulation with one unsolved scalar closure equation; `fqlx`'s checkpoints: the target has a nonzero `H^{3,1}` component (so no `H^{2,2}`-only construction), a conditional Lean pole-three clearing adapter, a fixed-current obstruction and unfinished scalar-curl repair, recovered packets; `DiamondJM`: holdout-q2 line construction does not transfer, the target is a mixed second variation, transcendental part is the whole target sector; `benediamond`: five-point obstruction to characterwise monomial surface strategies, the visible 216-element symmetry group has no uniruled quotient. No certificate claimed by anyone.) **Deformation-related open:** `off-paired` (holdout chain + `X0²X1X2Y0Y2/3`: breaks separate paired signs; the term is odd under simultaneous negation of `(X0,Y0)`, so the `u_i = X_i²` cubic substitution of this method is unavailable), `one-reflection-skew`, `quadric-contact-residual` (paired signs survive but the generic auxiliary elliptic surface has `χ(𝒪) = 3`, Kodaira dimension 1 — not uniruled — so a Jacobian curve family with fixed passive parameters cannot be found; `scripts/audit-jacobian-scope.py`). Solved relatives: `holdout` family, `mixed-biweight` (transitive cubic, `fqlx`), `mixed-polar`.

---

## 6. cyclic-triple-fano — "Cyclic triple cover · Fano density"

**Record.** `q = 1`, `F = Σ(X_i⁶+Y_i⁶) + (5/7)(X0³X1Y0Y1 + X0³X2³ + X0²Y0Y1Y2² + X0X1Y1³Y2 + X0X2³Y1Y2 + X1X2³Y0Y1)` (12 monomials), `A = X0X2²Y0³`. Discussion #112 (`rbajaj5` free line 2026-10-17; `benediamond` confirmation 2026-10-19). PR #153 by `benediamond` (OpenAI Codex-assisted; free line credited to `rbajaj5`; machinery to `fqlx`, `HeadOfNaps`), accepted 2026-10-20. 15 Lean files, 2158 lines, 96 KB. `checks.json`: smooth reduction mod 5; no coordinate reflections; Hessian symmetrizer dimension 1; projective diagonal group order 18; junior elements fixing the target `(0,0,2,0,0,0)`, `(0,0,4,0,0,0)` (i.e. `X2 ↦ ζ3^{±1}X2`). `sameVarietyInstances = [cyclic-triple-residual]`.

**Structure.** With `Z = X2³`: `F = G(X0,X1,Z,Y0,Y1,Y2)`, `G = Z² + (5/7)Z(X0³ + X0Y1Y2 + X1Y0Y1) + H` of weighted degree 6 for weights `(1,1,3,1,1,1)` (sum `M = 8`); `Φ = (X0,X1,X2³,Y0,Y1,Y2)`, `det dΦ = 3X2²`; `A = det(dΦ)·Φ*(X0Y0³/3)` with quotient source `B = X0Y0³` of weight `4 = 12 − M`. The quotient `W = {G = 0} ⊂ ℙ(1,1,1,1,1,3)` is Fano of index 2; existence was known nonconstructively from the weighted-Fano density theorem (Chapter 3 §3.4), which serves exactly the targets `det(dΦ)·b(Φ)`.

**Method.** Free weighted line + anchored polynomial family + polynomial-family integration on the Fano quotient + cubic pullback ("cubic-cover line and its density", §3.4). **Nonconstructive**; no denominator degree.

**Mechanism** [proved in Lean; numbers as in `solution.md`].
1. Free line `[X0:X1:Z:Y0:Y1:Y2] = [0:s:0:is:t:it]` on `W` (rbajaj5; normal map `H⁰(𝒪⁵ ⊕ 𝒪(2)) → H⁰(𝒪(5))` of rank 6, minor `−22500/2401`; `MM* = diag(72, 25/49 ×4, 72)` per benediamond). Chart `s = 1`, proof coordinate order `(Y0,Y1,X0,X1,Y2,Z)`, curve `(i, t, 0, 1, it, 0)`.
2. Anchored family through a generic endpoint `(ξ,0,k0,k1,k2,k3)` on `Y1 = 0` with `ξ⁶ + k0⁶ + k1⁶ + k2⁶ + k3² + (5/7)k0³k3 = 0`: `Y0 = ξ + at`, `Y1 = t`, `X0 = k0 + bt`, `X1 = k1`, `Y2 = k2 + ct`, `Z = k3 + dt + et² + ft³`; six coefficient equations (`t¹..t⁶`) with Jacobian determinant `−22500/2401` at the seed `(a,b,c,d,e,f) = (0,0,i,0,0,0)`, endpoint `(i,0,0,1,0,0)` (`LineSeed.lean`: explicit `seedInverse`, `seedMatrix * seedInverse = 1` by `fin_cases`/`norm_num`/`ring_nf`).
3. Endpoint ring integral: monic in `ξ`, specialization `(k0,k1,k2,k3) = (0,T+1,0,i)` gives `ξ⁶ + (T+1)⁶ − 1`, Eisenstein at `T` (`LineEndpoint`, `CyclicDomain`). Shared `square_system_generic_lift` gives a finite extension `E` of the endpoint function field carrying all six coefficients (`LineLifting`). Polynomial normal with constant value `(H_{Y0}(c(0))⁻¹, 0, …, 0)` (`LineNormal`). Dominance: the four endpoint coordinates plus `t` make evaluation injective; `E(t)/K` finite with both of transcendence degree 5 (`LineEmbedding`). Five commuting derivations, nonzero normal-plus-tangent determinant (`LineSource`). Integration + normalized trace + legal clearing via shared `ReflectionCoupledFano.normalized_data_from_polynomial_family` (`LineClosure.quotient_normalized_of_line_family`): `(d, N, Q)` on `W` with `gcd(H,d) = 1`, `ΣN_jH_j = B·d`, `d·ΣN_j' − ΣN_j∂_jd = H·Q`.
4. Cubic pullback `φ : Z ↦ X2³` (`CubicLift.cubeVariable`): `D = φ(d)`, `Ñ_j = X2²φ(N_j)` (`j ≠ Z`), `Ñ_{X2} = φ(N_Z)/3`; chain rule gives `ΣÑ_jF_j = A·D` and `D·ΣÑ_j' − ΣÑ_j∂_jD = F·X2²φ(Q)` (`cubic_normalized_identities`). Legality (`cubeLast_legal`): `φ` is injective on coefficients and `G` is monic in the unaffected coordinate `Y0`, so `F ∣ D` would descend to `H ∣ d`; `physical_irreducible` upgrades non-divisibility to `IsRelPrime`.
5. `hasCertificate_of_normalized` with chart index 0 (`hypersurface_not_dvd_pderiv_zero`), then `certificate_rename` (`CoordinateOrder.lean`, transport of `Verify` under a `Fin 6` permutation) returns to the repository order: `coordinateOrder = (0,1,2,3,4,5) ↦ (3,4,0,1,5,2)`.

**Coprimality over ℂ** [proved]. `CyclicDomain.polynomial_irreducible` (monic `X⁶ + …` over `MvPolynomial (Fin 5) ℂ`, specialization with `η⁶ = −1`), `physical_irreducible coupling`, `weighted_quotient_irreducible`; `cubeLast_legal`.

**Lean architecture.** `Solution.lean` (1.9 KB): `protected_F`, `protected_A` by `simp only [...]; ring` after `rename coordinateOrder`; `Submission.certificate := certificate_rename … physical_certificate`. `Certificate.lean` (61 lines, `maxHeartbeats 1500000`): `quotient_certificate_data` (obtains the generic field `I` from `Lifting.exists_field`, `letI` instances) → `physical_certificate`. Shared imports: `HoldoutCover.{GenericIncidence, AlgebraicFields, HypersurfaceDomain, PolynomialCalculus, ResidueDerivatives, ResidueClearing, FiniteDifferentials, EvaluationMatrix, PolynomialNormal, HomogeneousCertificate}`, `PolynomialCurves.{ArithmeticLifting, NormalizedCertificate}`. Heartbeats `1000000–2000000`; `LineIncidence.lean` and `QuotientGradient.lean` use `maxHeartbeats 0` (unlimited); `maxRecDepth 10000`; `backward.isDefEq.respectTransparency false` (set back to `true` in `LineLifting`). Tactic census (15 files): `simp` 154, `rfl` 44, `decide` 35, `norm_num` 22, `ring` 14, `omega` 5, `linear_combination` 4, no `field_simp`/`native_decide`.

**Symmetry data** [computed; consistent with `checks.json`]. Projective diagonal group order **18** (`[3,6]`); the residue of `A·Ω/F²` is invariant under `X2 ↦ ζ3X2` (character `ζ3²` of `A` cancels `ζ3` of the volume) but has order 6 on the full group.

**Failed / bounded attempts.** None recorded against this target in #112; the thread went directly from the free line to the proof. `benediamond`'s 2026-10-19 note cautions that the free line itself has zero value of the density `X0Y0³/3` and that a nonzero deformation of the density is not by itself a global Gysin statement — the anchored family resolves this.

**Related.** **Same variety, open:** `cyclic-triple-residual` (`A = X0³Y0³`): cubic pullback supplies only sources `X2²·B(…,X2³,…)` (all `X2`-exponents ≡ 2 mod 3), the proposed source has exponent 0; residue character under `X2 ↦ ζ3X2` is `ζ3` (non-invariant) — the Fano density theorem gives nothing; discussion #155: `fqlx` degree-84 model with integral coordinate, rank-five marked connection, research stopped 2026-10-23 without a certificate. **Same mechanism:** `reflection-coupled-fano` (solved; reflection quotient `Y2 ↦ −Y2`, weighted quartic family, PR #150 machinery) vs `reflection-coupled-residual` (open; even target, wrong reflection character — the exact analogue of the residual pair here); `holdout-moving-branch`, `holdout-ternary-deformation` (solved; anchored families). **Other Fano-density open targets:** `quadratic-cubic-pullback` (six quadrics → Fermat cubic fourfold; §3.4 serves the constant-source density `δ = det(∂Q_i/∂Z_j)`, which reduces mod 5 to `4X0X1X2Y0Y1Y2`, but `X0³Y0³` is a different Jacobian-ring class), `symmetric-sextic-residual` (the `S3×S3` quotient density is `Δ_XΔ_Y`, not the symmetric monomial).

---

## 7. Cross-cutting observations

1. **Denominator families seen in this group.** Hyperplane arrangements (contact46006: six hyperplanes over ℚ; coupled-pair: four hyperplanes over ℚ(ζ₁₂), irreducible quartic over ℚ), a product of two quadric cones (f3-block), and three unknown-degree denominators produced by trace/choice (crossed-ternary, cycle-three, cyclic-triple-fano) — for cycle-three the degree is ≥ 7 by rbajaj5's exclusion. The three explicit solves all have strict certificates (`Q = 0`) and degree-`(d+1)` numerators with the right parities.
2. **All three explicit solves are "reflection-served" or "contact-served" sectors**: f3-block and coupled-pair have residues invariant under coordinate reflections preserving `F` (so Chapter 3's reflection theorem already gave existence); contact46006's target is a polarized difference `B(U) − B(V)`. None of the three touches the sectors that make the open companions hard (nontrivial diagonal characters, swap-invariant K3 remainders, even targets under the only available reflection).
3. **Coprimality over ℂ is never a bottleneck.** Pattern (A): `D` omits a variable in which `F` is monic ⇒ specialization to a constant. Pattern (B): irreducibility of `F` by one Eisenstein specialization (the same `X⁶ + ((X+1)⁶ − 1)` / `(T−1)⁶ − 1` trick each time) and an evaluation with `ρF = 0 ≠ ρD`.
4. **Big identities**: the only place where a large polynomial identity is kernel-checked rather than `ring`-normalized is cycle-three's covariant syzygy (190 graded coefficient identities via a proved digit-encoding, `decide +kernel`, 16 M heartbeats, recDepth 131072). Everything else is `simp`/`ring`/`field_simp`/`linear_combination` on function-field identities, plus `decide` on finite-field seed matrices (𝔽₅ in crossed-ternary and cycle-three, ℂ with explicit inverse in cyclic-triple-fano).
5. **Characters** (projective diagonal group orders / residue character orders, [computed]): contact46006 6/6; coupled-pair 864/6; crossed-ternary 24/6; cycle-three 72/3; cycle-q2 72/3; cyclic-triple-fano 18/6; cyclic-triple-residual 18/6; f3-block 2592/6. The pairs (cycle-three, cycle-q2) and (cyclic-triple-fano, -residual) differ in *which* elements act nontrivially, not merely in the order: cycle-q2's residue is nontrivial on the weights-`(2,1,0,3,2,1)` ℤ/6 (character `ζ6²`), cyclic-triple-residual's residue is nontrivial on `X2 ↦ ζ3X2`.

## 8. Computation log (reproducible)

- `/tmp/hodge/scratch/verify_witness.py` (sympy 1.12, `python3 -I`): parses `instance.rec` + `witness.json`, checks E1, E2, the cleared identity, degrees, parities, `gcd_ℚ(F,D)`, `factor_ℚ(D)`, and the structural claims for contact46006 (`f = h² + Ck`, `F = f(U)+f(V)`, `D = C(U)C(V)`, `W(f) = C·A_C`, `A = A_C(U) − A_C(V)`, `N` formulas), coupled-pair (`F = F0 + D·Y1²`, `N(Y1²) = 0`), f3-block (`N(D) = D div N`). All checks returned `True`.
- `/tmp/hodge/scratch/characters.py`: Smith normal form of the exponent-difference lattice of `F` inside the sum-zero lattice (projective diagonal group), reflections, and the order of the residue character `w = exp(A) + (1,…,1) − (q+1)·m_F` in `Λ_0/L(F)`. Orders 24 and 18 agree with the repository's `checks.json` for crossed-ternary and cyclic-triple-fano.
- Singular 4 (`/usr/bin/Singular`, also in the `ddg-cas:latest` docker image together with giac; that image has no sympy/sage): `factorize(x^3+5x^2y−22xy^2−56y^3) = (x−4y)(x+2y)(x+7y)`.
- `/tmp/hodge/scratch/fetch_disc.py N`: read-only `gh api graphql` dump of discussion `N` (used for #55, #61, #72, #77, #91, #112, #114, #155). Nothing was posted.
