# Accepted solves, group 5: holdout-q2, holdout-ternary-deformation, holdout-x04-p, holdout-x0-x14-y0, holdout-x0-y03-y12, holdout-x12-p-y02

Catalogue notes prepared 2026-11-01 from the local checkout of `hodge-computer/sextics` at `83fbf5c` (branch `main`).
Nothing in the repository was modified; no Lean build was run. Every statement below is labelled as one of

- **[Lean]** a kernel-checked fact from the accepted sources (acceptance.json for all six: Lean 4.33.0, `kernelReplay: true`, axioms exactly `propext, Quot.sound, Classical.choice`, `witnessChecked: false`, result `lean-checked-awaiting-owner`, owner review by `benediamond`);
- **[CAS]** an identity I re-derived independently today (sympy 1.12 on the host, run with `python3 -I` from `/tmp/hodge/check/group5_checks.py`; the `ddg-cas:latest` docker image has Singular 4, python-flint and Julia but **no sympy**, so I used Singular in the container only for a second check of the quintic identity, `/tmp/hodge/check/quintic.sing`);
- **[Doc]** a claim taken from `solution.md` / `why.md` / discussion threads that I did not re-verify;
- **[Obs]/[Conj]** my own observation or conjecture.

Conventions. Variables `Z = (X0,X1,X2,Y0,Y1,Y2)`; `∂_i = ∂/∂Z_i`; `N(F) = Σ N_i ∂_iF`; `div N = Σ ∂_iN_i`; `B_D(N) = D·div N − N(D)`.
Protected statement (`formal/lean/HodgeComputer/Verify.lean`): `HasCertificate F A q := ∃ d D N Q, Verify F A q d D N Q` with
`D ≠ 0`, `D` homogeneous of degree `d`, each `N_i` homogeneous of degree `d+6q−5`, `IsRelPrime F D`, and
`A·D² − q·D·N(F) + F·B_D(N) = F^(q+1)·Q`.
Normalized form (sufficient, `normalized_identity`): `q·N(F) = A·D` and `B_D(N) = F^q·Q`.

The six instances fall into three mechanisms:

| instance | q | A | author / PR | mechanism | constructive? |
|---|---|---|---|---|---|
| holdout-q2 | 2 | `X0^3 X1^2 X2^2 Y0^3 Y1^2` | `HeadOfNaps`, PR #124, accepted 2026-10-19T14:03Z | scaled endpoint reflection + cubic reduction + degree-27 line algebra + Euler pole-2 correction | yes (structured witness, strict `Q = 0`, clearing exponent existential) |
| holdout-ternary-deformation | 1 | `X0^3 Y0^3` | `hodgepodgeconjecture`, PR #151, accepted 2026-10-20T12:09Z | Fermat cover + **new** mod-5 conic seed + étale lifting + polynomial integration + deck trace | no (existential `D`, `Q` allowed nonzero) |
| holdout-x04-p | 1 | `X0^4 X2 Y2` | `fqlx` (Codex for `benediamond`), PR #140, 2026-10-19T20:04Z | specialization of the shared Fermat-cover library (`certificate_of_source_pair`) | no |
| holdout-x0-x14-y0 | 1 | `X0 X1^4 Y0` | `fqlx`, PR #138, 2026-10-19T19:53Z | same | no |
| holdout-x0-y03-y12 | 1 | `X0 Y0^3 Y1^2` | `fqlx`, PR #139, 2026-10-19T19:58Z | same | no |
| holdout-x12-p-y02 | 1 | `X1^2 X2 Y0^2 Y2` | `fqlx`, PR #142, 2026-10-19T20:17Z | same | no |

None of the six folders contains a `witness.json`; `proof.json` is in every case `{"schemaVersion": 2, "theorem": "Submission.certificate"}`.

Chronology worth noting [Lean/acceptance.json]: holdout-q2 (q=2) was accepted at 14:03Z on 2026-10-19, about three hours **before** the central holdout q=1 instance (PR #128, 17:08Z). The four q=1 specializations were accepted the same evening, after `fqlx`'s library had been generalized to arbitrary sources (`GeneralizedCover.lean`).

The fixed holdout sextic for five of the six instances:

```
F_hold = X0^6 + X1^6 + X2^6 + Y0^6 + Y1^6 + Y2^6 + X0^3 X1 Y0 Y1 + X1^3 X2 Y1 Y2
```

Symmetries used below [Doc, partly CAS]: the three paired sign changes `(X_i,Y_i) ↦ (−X_i,−Y_i)`; the exchange `X2 ↔ Y2`; the diagonal `μ6` action with weights `(0,0,0,5,1,5)`; and the **scaled endpoint reflection** `(X2,Y2) ↦ (k·Y2, k²·X2)` for `k² + k + 1 = 0` (checked: it preserves `F_hold`).

---

## 1. holdout-q2

**Data.** `F = F_hold`, `q = 2`, `A = X0^3 X1^2 X2^2 Y0^3 Y1^2` (degree 12). The proof writes the protected variables as `(x,a,c,y,b,d) = (X0,X1,X2,Y0,Y1,Y2)` so `F = x^6+a^6+c^6+y^6+b^6+d^6 + x^3ayb + a^3cbd`, `A = x^3a^2c^2y^3b^2`. Same target on Fermat is the accepted control `fermat-q2` (PR #37). Discussion #60.

**Method name.** "Endpoint-reflection cubic reduction with a line-algebra contact solver and an Euler pole-two correction." Author `HeadOfNaps`; the endpoint-reflection observation and preliminary construction were developed with ChatGPT assistance for Jason Long; Codex assisted the formalization [Doc]. It reuses (by copying, not importing) infrastructure from the accepted `quadric-contact-broken` (PR #122): coefficient localization, finite-free `AdjoinRoot` specialization, matrix contact calculus, homogeneous clearing, matrix Bezout identities (namespaces `QuadricContact*`, `RegularSpecialization`, `AdjoinRootSpecialization` live in the local files). The `one-reflection` solve informed the approach.

**Constructive, strict.** The final witness has `Q = 0` (strict certificate), but it is *not exported*: the denominator is a structured product whose homogenizing exponent `n` is obtained existentially (`exists_q2_representatives_with_fibre`). So: constructive in structure, existential in one integer. No `witness.json`.

**Mathematical mechanism (step by step).**

1. *Target split by the scaled endpoint reflection* [CAS]. With `P = x^3a^2y^3b^2` and `Δ_k = c² − k²d²`,
   `P c² = (2+k)/3 · P Δ_k + (2+k²)/3 · P Δ_{k²}` modulo `k²+k+1`. Each endpoint target `P Δ_k` (resp. `k²`) receives its own strict certificate; the two are combined over a product denominator (`HoldoutReflectionCombination.combine_q2`: `D_{A+B} = DE`, `N = l·E·N + m·D·M`, `Q = l E² Q + m D² R`). The reflection `(c,d) ↦ (kd, k²c)` sends `Δ_k ↦ −Δ_k` and fixes `F` (`delta_swap`, `sextic_reflection`) — this is the motivation only; no reflection theorem is used.
2. *Cubic reduction* [CAS]. Active coordinates `U = x², V = a², S = c² + k²d², T = k c d`; passive ratios `α = y/x`, `β = b/a`. For
   `G(𝖺,𝖻,𝖼,𝖽) = S³ − 3ST² + 𝖺U³ + 𝖻V³ + 𝖼U²V + 𝖽V²T` and `(𝖺,𝖻,𝖼,𝖽) = (1+α⁶, 1+β⁶, αβ, βk²)` one has exactly `G = F` (`cubic_ratio_substitution`). "Good fiber" `(α,β) = (k², k)` gives `(𝖺,𝖻,𝖼,𝖽) = (2,2,1,1)`; all unit checks are done in the localization of `ℚ[𝖺,𝖻,𝖼,𝖽]` at polynomials nonvanishing there (`RegularSpecialization.Coefficients point`).
   The endpoint source on the cubic is `U³V²` (passive scalar `α³β²/2` restores `P/2`, hence `P Δ_k` after the horizontal lift; `pullback_source_density`: `2·Δ_k · pullback(source) = reflectionTarget`).
3. *Explicit quintic Jacobian reduction* [CAS, also Singular]. Four explicit cubic polynomials `M_i ∈ ℚ[𝖺,𝖻,𝖼,𝖽][U,V,S,T]` (file `QuinticReduction.lean`) satisfy
   `Σ_i M_i ∂_iG = H · U³V²`, `H = R(R⁴ − 2916 𝖺⁸𝖽⁶)`, `R = 27𝖺²𝖻 + 4𝖼³`, and `H(2,2,1,1) = 515198970880 ≠ 0`.
   I re-derived this identity in sympy and in Singular (both return `lhs − rhs = 0`, same value of `H`). Then `N := M/H` is a polynomial primitive for the pole-one source and `h := div N` is an (explicit) quadratic form in `U,V,S,T` (its 10 coefficients are in `QuadraticDivergence.coefficients`). This replaced a preliminary 56×56 inversion [Doc].
4. *Lines on the cubic serve every quadratic source* [Lean, Doc]. A line `S = eU + fV`, `T = gU + hV` lies on `G = 0` iff `rel0..rel3 = 0` (`LineContact`). With `f = re`, elimination gives a degree-27 polynomial in `r` with leading coefficient `729𝖺¹¹` (`NewEliminant.rawPolynomial`; `quadratic_factor` is the exact factorization identity `qa·num² + qb·num·den + qc·den² = 18 d⁴(3ar−c) κ⁵ · eliminant`). It is normalized monic and `E = AdjoinRoot P` is a free rank-27 algebra; the line coordinates are rational in the root with unit denominators proved by Bezout identities on the good fiber. The elementary line field `(1/ℓ)((Q₂)_V ∂_U − (Q₂)_U ∂_V)` with `ℓ = S − eU − fV`, `G = ℓQ₁ + mQ₂`, has zero divergence and sends `G` to a quadratic (`LineContact.certificate_equations`: `Equations (cubic) (source) (pole e f) (numerator) 0`). Four constant gradient fields + six coefficient functionals of the line source give an invertible 10×10 matrix on the good fiber in the basis `U²,UV,US,UT,V²,VS,VT,S²,ST,T²` (`GoodFiberSources`, `GenericQuadratic`), hence every quadratic source — in particular `h` — is served. The line pole is cleared by the regular representation: `D₀ = det(S·I − U·M_e − V·M_f)` (27×27, degree 27 in `U,V,S`), numerators `P_i` of degree 27 via adjugate/column functional (`QuadraticSolver`). The determinant is never expanded.
5. *Euler pole-two correction* [Lean]. With `W₁ = P/D₀` (serving `h`) and the four-variable Euler field `E`,
   `W₂ = (N + G·W₁)/2 − h·E/6` satisfies `2 W₂(G) = U³V²` and `div W₂ = 0` (`Q2Reduction.q2_equations`, proved over any commutative ring from `E(G)=3G`, `E(h)=2h`, `div E = 4`). Cleared form: denominator `6D₀` (degree 27), numerators `3(D₀N_i + G P_i) − D₀ h X_i` (degree 30). This is the step that "raises the pole order" and is the reusable idea.
6. *Passive substitution, horizontal lift, homogeneous clearing* [Lean]. `PassiveQ2` maps coefficients `𝖺,𝖻,𝖼,𝖽 ↦ (1+α⁶,…)` into the localization `Regular k`; `D_axis` shows `D` evaluates to `6` on the `S`-axis, hence `D ≠ 0` survives the non-injective substitution. Horizontal six-variable fields (with `T = kcd` normalized)
   `C_U = (Δ_k/x)∂_x + (Δ_k y/x²)∂_y`, `C_V = (Δ_k/a)∂_a + (Δ_k b/a²)∂_b`, `C_S = c∂_c − d∂_d`, `C_T = −2kd∂_c + 2k²c∂_d`
   have zero divergence, annihilate `α, β`, and act on `(U,V,S,T)` as `2Δ_k` times the coordinate derivatives (`EndpointHorizontal`, derivations on the genuine fraction field via the quadric-contact `FractionDerivation.extend`). `HoldoutActiveClearing.exists_q2_representatives_with_fibre` clears the finitely many passive-ratio denominators with a multiplier `(xa)^n`: `Q` homogeneous of degree `2(n+27)`, `P_i` of degree `2(n+30)`. `EndpointLift.numerators` then gives the six-variable data
   `D_end = x²a²·Q` (degree `2(n+27)+4`), `N_i` of degree `2(n+30)+5 = deg D_end + 7` (as required for `q=2`), and `lift_strict_q2` proves the cleared defect is identically zero.
7. *Denominator legality over ℂ* [Lean]. On the fiber `(x,a,y,b,d) = (1,1,k²,k,0)`: `F ↦ c⁶ + 5` [CAS] (monic in `c`), and `D_end ↦ 6u·det(c²I − M)` with `M = ρ(e+f)` the 27×27 regular-representation matrix of the line sum (rational entries, cast to ℂ). An exact inverse `K` of `(e+f)³ + 5` in `AdjoinRoot P` gives `(M³ + 5I)K = I` (`DenominatorFiberLink.complex_matrix_inverse`, obtained by `decide +kernel` on padded `List ℚ` coefficient arithmetic). `QuadricContactMatrixBezout.relPrime_cubic_det` turns this into `IsRelPrime (z³+5) (det(zI − M))` with `z = c²`, i.e. `IsRelPrime (c⁶+5) (det(c²I − M))` **without expanding the degree-54 determinant**. Then `HodgeComputer.relPrime_of_specialization` (needs only that `F` has unit leading coefficient in `c`) and `relPrime_of_retract` (`toUni`/`fromUni`, `c = X2` as the univariate variable) lift relative primality to `ℂ[Z]`. `D ≠ 0` follows from legality since `F` is not a unit.

**Denominator family.** Per endpoint `k`: `x²a² · (xa)^n · 6·det(S·I − U·M_e − V·M_f)|_{passive coefficients substituted}`, i.e. a *norm from the degree-27 line algebra over the cubic*, homogenized; product of the two endpoint denominators for `A`. Degree `2(n_k+27)+4 + 2(n_{k²}+27)+4` with `n` unspecified. [Obs] Nothing is proved about minimality; `fermat-q2` (same `A`, Fermat `F`) has the strict degree-12 certificate `v = A(Fε_X − Dε_Y)/(12D²)` with `D = X0⁶+X1⁶+X2⁶` (Chapter 1), so the holdout version is vastly larger than the control.

**Lean architecture.** 42 files, 5414 lines (largest: `GoodFiberSources.lean` 92 KB / 283 lines, `NewSpecialization.lean` 56 KB, `DenominatorFiberLink.lean` 26 KB — these are kernel-checked rational coefficient lists). Imports: only `HodgeTarget`, `HodgeComputer.Verify`, `HodgeTestbed.Algebra`, Mathlib, and the local `Proposal.*` files; **no `HodgeComputer.Constructions.*` import**. Main chain:
`Submission.certificate` (rewrites `HodgeSubmission.F = HoldoutEndpointTarget.F` by `unfold; ring`, `A` by `rfl`)
→ `HoldoutFinalAssembly.certificate` → `HoldoutEndpointTarget.hasCertificate_of_endpoint_solver` (root `k` from `Complex.exists_root` on `X²+X+1`, `compute_degree!`)
→ `hasCertificate_of_endpoints` (`combine_q2` + `target_combination`, `linear_combination`)
→ `endpointCertificate k` → `endpointCertificate_of_norm` (with `PassiveQ2.equations`, `DenominatorFiberLink.complexSumMatrix/complexInverseMatrix`, `passive_q2Denominator_fiber`)
→ `HoldoutSexticAssembly.hasCertificate_endpoint_from_cubic` (uses `exists_q2_representatives_with_fibre`, `homogeneous_clearing_source/divergence`, `denominator_legal_of_scaled_normFiber`, `lift_strict_q2`; produces `⟨2(n+d)+4, x²a²Q, numerators, 0, …⟩` with `d = 27`).
Tactics: `ring`, `linear_combination` (modular identities in `k`), `norm_num`, `field_simp`, `polynomial` (for `toUni_F`), `compute_degree!`, `decide +kernel` (list arithmetic over ℚ padded to length 32–80), `fin_cases`, `MvPolynomial.ringHom_ext`, `Matrix` adjugate/det API. Coprimality: `relPrime_of_specialization` + `relPrime_of_retract` + matrix Bezout (**not** irreducibility of `F`). Heartbeats: `maxHeartbeats 0` in `DenominatorFiberLink`, `GoodFiberMinorLink`, `GoodFiberCoordinates`, `GoodFiberSources`, `NewEliminant`, `NewLineAlgebra`, `NewSpecialization`; `4000000` in `LineContact`; `2000000` in `GenericQuadratic`, `QuadraticDivergence`, `QuadraticSolver`, `QuinticReduction`; `maxRecDepth` up to `1000000`. Sizes of exact data: `P` is monic of degree 27 over ℚ with coefficients like `12129697/1492992`; the line coordinates `e,f,g,h` have ~60-digit rational coefficients.

**Reusable facts / corollaries.**
- [Lean] `Q2Reduction.q2_equations`: from any polynomial Jacobian primitive `N(G) = B` with homogeneous cubic `G` (degree-3 `N_i`) and any strict certificate `(D,P)` for `div N`, the Euler-corrected pair is a strict pole-two certificate for `B` with denominator `6D`. Ring-generic; worth extracting to the shared library if a `q=2` instance on another cubic reduction appears.
- [Doc, discussion #60, `benediamond`] because the endpoint field kills `α, β`, multiplying by `α^{i−3}β^{j−2}` gives strict `q=2` certificates for `A_ij = X0^{6−i}Y0^i X1^{4−j}Y1^j X2²` (`0≤i≤6`, `0≤j≤4`, `i+j` odd) and, via `X2 ↔ Y2`, the same with `Y2²`: 34 monomial targets on `F_hold`. Not submitted to Lean.
- [Obs] The cubic reduction needs: two pairs with passive ratios fixed and a third pair entering `F` only through `c², d², cd` (here `c⁶+d⁶` and `a³·c·b·d`). The `μ3`-scaled reflection is what makes `(S,T) = (c²+k²d², kcd)` the right invariants: `c⁶+d⁶ = S³ − 3ST²` holds **only** with the `k²` weighting (checked symbolically via `cubic_substitution`).

**Relations to open instances.**
- `cycle-q2` (OPEN, `q=2` on the cyclic sextic `F_hold/2-coefficients + X0X2³Y0Y2/2`, `A = X0²X1²X2³Y0²Y1²Y2`). [Obs] The closing coupling `X0X2³Y0Y2 = α x² c³ d` is **not** a polynomial in `(U,V,S,T)` (since `c²` alone is not a rational function of `S = c²+k²d²` and `T = kcd`: `c²` and `k²d²` are the two roots of `z² − Sz + T²`), so the endpoint cubic reduction does not apply verbatim; `why.md` additionally notes the residue has nontrivial `μ6`-character `ζ6²` so invariant-quotient methods miss it. The transferable part is step 5 (Euler pole-two correction from a Jacobian primitive plus a `q=1` certificate for its divergence); the cycle-three solve (PR #165, Jacobian/Weierstrass model) would have to supply the pole-one primitive on *its* model.
- `fermat-q2`, `block-q2`, `split-ternary-q2` (all accepted) are the other `q=2` solves; `split-ternary-q2` (`fqlx`, third-order residues, `HigherResidues` library) is the only other mechanism for raising the pole order on a non-separated equation.
- `quadric-contact-residual`, `one-reflection-skew`: the infrastructure (line algebra, matrix Bezout coprimality on a fiber) is shared with `quadric-contact-broken`; the residual target there is reflection-invariant, which is exactly the situation (`A = P c²` fixed by the endpoint reflection) that holdout-q2 handled by *splitting* `c² = (2+k)/3 Δ_k + (2+k²)/3 Δ_{k²}` into two reflection-odd pieces for two different scalings `k, k²`. [Conj] This "two scaled reflections whose odd parts span the invariant" trick requires an order-3 scaling symmetry of the pair (here `c⁶ + d⁶` with the `a³cbd` coupling); it does not exist for `X0 ↔ Y0` in `quadric-contact-residual` where only the plain exchange is a symmetry.

---

## 2. holdout-ternary-deformation

**Data.** `F = F_hold + (X0⁵X1 + 2X1⁵X2 + X2⁵X0)/14`, `q = 1`, `A = X0³Y0³`. Author `hodgepodgeconjecture` (Codex assistance), PR #151. Discussion #147 (one comment, by the author, announcing the mod-5 seed as a partial result the day before the PR). `provenance.json` records `existenceStatus`: geometric existence was known beforehand (generic weighted-Fano quotient, uniruledness, deck-character trace) but not formalized. `checks.json` (curation, generated by `scripts/audit-cover-frontier.py` with sympy) records: smoothness via maximal-coordinate derivative bounds `[24/7, 67/14, 3/2, 1, 2, 1] < 6`; Hessian symmetrizer dimension 1; no coordinate reflections; projective diagonal group of order 6, target-fixing subgroup of order 2, no serving junior elements; Jacobian ranks `[36, 37]`; `originalCoverSwapDifference = 0` (the `y ↔ T` exchange survives); generic quotient quasismooth at `(α,β) = (1,0)`, weights `(1,1,1,2,2)`, degree 6, anticanonical index 1, three ambient quotient points of type `½(1,1,1)`; `originalConicDefect = k⁵u⁵v⁷/7 + k u¹¹v/14 + u²v¹⁰/14` — [CAS] this is exactly the ternary correction evaluated on the central holdout conic `(X0,X1,X2) = (s², k s t, t²)`, so the old conic is *not* on the deformed quotient.

**Method name.** "Replacement conic for the Fermat cover: exact nonsingular incidence over `F5`, étale/flat lifting to a generic finite field factor, polynomial integration, deck trace, legal clearing." Nonconstructive.

**Shared modules imported** (all under `HodgeComputer.Constructions.HoldoutCover.`): `GenericIncidence` (étale localization of a square polynomial system, `square_system_jacobian_etale`, `square_system_jacobian_nontrivial`, `genericEtale_exists_fieldFactor`, `squareGeneric_exists_fieldFactor`), `PolynomialNormal` (`full_polynomial_combination_one`), `FiniteTrace`, `PolynomialFamilyBridge`, `EvaluationMatrix`, `FiniteDifferentials`, `PolynomialCalculus`, `PolynomialFieldMap`, `AlgebraicFields`, `FermatFields`, `HypersurfaceDomain`, `ResidueClearing`, `ResidueDerivatives`, `PhysicalCoverEmbedding`, `QuotientCoverEmbedding`, `CoverFields`, `CoverCoordinateDerivatives`, `CoverScalars`, `CoverSourceBinding` (`scaledCoverQuotientSource`, `cover_source_identity`, `cover_normal_source_identity`), `FermatCoefficientDerivations`, `CoverDeckField` (`normalizedTrace_inv_fourth_zero`, `exists_nontrivial_sixth_root`), `FiniteResidueClearing` (`finite_residue_solution_clearing_value`), `HomogeneousCertificate` (`hasCertificate_of_normalized`), `Piola`, `ResiduePiolaPullback` (`residue_solution_piola_pullback`). Local: 27 files, 3042 lines.

**Mathematical mechanism.**

1. *Cover and quotient* [CAS]. `T⁶ = Y0⁶ + Y2⁶`, `α = Y0/T`, `β = Y2/T`, `y = Y1`, `S = y² + T²`, `P = yT`; `α⁶ + β⁶ = 1`. Then
   `14·F(X, αT, y, βT) = H(α,β)(X0,X1,X2,S,P) + 14·T⁶(α⁶+β⁶−1)`, with
   `H = 14(X0⁶+X1⁶+X2⁶ + P(αX0³X1 + βX1³X2) + S³ − 3SP²) + X0⁵X1 + 2X1⁵X2 + X2⁵X0` (`TargetIdentities.protected_cover_equation`; Lean's `linear_combination 98*t^6*hbase` coefficient reflects the post-`field_simp` normalization). The ternary correction involves only `X`, hence is trivially `y ↔ T` invariant: this is why the exchange, the quotient `(X0,X1,X2,S,P)` with weights `(1,1,1,2,2)`, and the deck trace all survive unchanged, while the conic does not.
2. *New seed over `F5`* [CAS]. At `(α,β) = (1,0)`:
   `x0 = 1 + 2t²`, `x1 = t`, `x2 = 4 + 2t + 2t²`, `S = 4 + t + t² + 3t³ + 3t⁴`, `P = 4 + t + 4t² + 4t³ + 3t⁴`.
   `H(curve)` has degree 12 in `t` and all 13 coefficients `[56590, 173305, 388710, 577270, 705900, 661860, 523520, 326680, 173960, 66370, 21710, 3040, 1100]` are divisible by 5 (Lean: `seed_incidence` exhibits `H = 5·(…)`). Unknowns: 13 coefficients (the 3 of `x0`, the `t, t²` coefficients of `x2`, the 4 positive-degree coefficients of each of `S, P`); constants `(x2(0), S(0), P(0)) = (4,4,4)` and `x1 = t` are fixed. Equations: `[t^i]H = 0`, `i = 0..12` (unlike the central holdout, the constant equation is included and `x0(0)` varies). The 13×13 Jacobian at the seed, reduced mod 5, has determinant `2 ≠ 0` and equals Lean's `explicitJacobian` (rows = equations, columns = variables) with `explicitJacobian · inverseJacobian = I` mod 5 — all re-checked today. The seed was found by modular experiment; **the proof does not claim the mod-5 point lifts**; it uses it only to show nontriviality of an étale algebra.
3. *Lifting* [Lean]. `integerEquations : Fin 13 → MvPolynomial (Fin 13) ℤ` (with `(α,β) = (1,0)`, parameters `(4,4,4)`); the Jacobian-localized incidence ring is étale over ℤ and nontrivial (it has an `F5`-point with unit Jacobian, `square_system_jacobian_nontrivial`); flatness gives a nonzero generic fiber over ℚ-type fields and a finite field factor; `IsAlgClosed.lift` embeds it into ℂ (`ConicLifting.complex_point_of_etale_nontrivial`, `ComplexConic.exists_complex_conic`). Then over `B = Frac(ℂ[α,β]/(α⁶+β⁶−1))` with three independent endpoint parameters `k1,k2,k3` (`ParameterRing = MvPolynomial (Fin 3) FermatRing`), the complex point shows the generic incidence algebra is nontrivial; `genericEtale_exists_fieldFactor` yields a finite extension `E/B(k1,k2,k3)` and a curve `c(t) ∈ E[t]⁵` with `H(c) = 0`, `x1 = t`, endpoint `(k1,k2,k3)`, unit coefficient Jacobian (`GenericConic.conic_equation`, `conic_jacobian_unit`). Classical choice is used only to pick the maximal ideal.
4. *Normal, dominance, integration, trace* [Lean]. `full_polynomial_combination_one` produces polynomials `n_i(t)` with `Σ H_i(c) n_i = 1`, only `n_0` having a (nonzero) constant term (`ConicNormal`); the 5×5 matrix `(n, ∂_t c, ∂_{k1}c, ∂_{k2}c, ∂_{k3}c)` has nonzero determinant (its value at `t=0` is the constant coefficient of `n_0`; `ConicFlow.evaluationMatrix_det_ne_zero`). `G` is irreducible: monic in `x0`, specialize to the shared Eisenstein example (the ternary correction vanishes under that specialization) — `DeformedDomains.quotientPolynomial_irreducible`. Transcendence degrees: `trdeg_B(FunctionField) = trdeg_B(TargetField) = 4` ⇒ `E(t)/K_q` finite (`ConicDimension`). `ConicResidue.quotient_residue_solution`: for every polynomial source `U`, there is `w` with `𝒟_g(w) = −(ρ(U_0)/g − ρ(U)ρ(G_00)/g²)` where `g = ρ(∂_0G)`.
5. *Transport back* [Lean, CAS]. Common field `L ⊇ K_p, K_q` (`DeformedCoverFields`; `trdeg_ℂ K_p = 5 = trdeg_ℂ L`, `trdeg_B L = 4`, so both extensions are finite). Coordinate Jacobian `J = 2α⁵(T² − y²)/T` (`actualCoverJacobian_det`), with `α ≠ 0`, `T ≠ 0`, `T² − y² ≠ 0` proved in the actual field. Source `U = X0³S/(2α²)` satisfies `J·U = A(1 − y⁴/T⁴)` and the same for the normal-divergence sources. The first and second `X0`-normals of `G` and `F` agree under the cover maps; the extra ternary contributions are `(5X0⁴X1 + X2⁵)/14` and `10X0³X1/7` [CAS]. Deck action `T ↦ ζT`, `(α,β) ↦ ζ⁻¹(α,β)`, `ζ⁶ = 1`, `ζ⁴ ≠ 1`, is an actual `K_p`-algebra automorphism of `L` (`DeformedDeck.physicalDeckAlgHom`) ⇒ `Tr^norm_{L/K_p}(T⁻⁴) = 0` (`cover_trace_inverse_fourth`).
6. *Clearing* [Lean]. `F` irreducible (`PhysicalModel.protected_irreducible`: `MulEquiv.irreducible_iff (finSuccEquiv ℂ 5)` + `physicalPolynomial_irreducible` with `η = I`, `I⁶ = −1`: specialize `X1 = X2 = Y1 = 0`, `Y2 = η` to get `X0⁶ + Y0⁶ − 1`, Eisenstein at `Y0 − 1`); `∂_0F ≠ 0` mod `F` by restriction to the `X0`-axis (`X0⁶ ↦ 6X0⁵`). `finite_residue_solution_clearing_value` gives `D₀, N', Q'` with `N'(F) = AD₀`, `B_{D₀}(N') = FQ'`, `IsRelPrime F D₀` (from `ρ(D₀) ≠ 0` and primality of `F`); `hasCertificate_of_normalized` performs finite homogeneous extraction (`D ↦ h^p`, numerators the degree-`(pe+1)` components of `S_p N`).

**Denominator.** Existential: `D = (∂_0F)·d` where `d` is any common polynomial denominator of the traced residue vector in `Frac(ℂ[Z]/F)`, then homogenized as `h^p` with `h` the legal initial form. No degree information. `Q` may be nonzero.

**Lean architecture.** Chain: `Submission.certificate := DeformedCover.deformed_hasCertificate`, whose body mirrors `GeneralizedCover.certificate_of_source_pair`: `quotient_residue_solution coverQuotientSource` → `residueDivergence_append_zero/map` → `residue_solution_piola_pullback` → `cover_normal_source_binding`, `cover_normal_binding` → `hasCertificate_of_cover_residue coverT cover_trace_inverse_fourth v`. Coprimality: irreducibility + nonvanishing (no `relPrime_of_specialization`). Options: `backward.isDefEq.respectTransparency false` in every file; `maxHeartbeats` between `50000` and `1000000` per file, with `set_option maxHeartbeats 8000000 in` for `gradientJacobian_eq` (`fin_cases i <;> fin_cases j <;> norm_num [...] <;> decide`) and `seedJacobian_inverse` (`decide` on a 13×13 product over `ZMod 5`); `maxRecDepth 100000`; `profiler true` left on in `ConicSeed`, `DeformedConic`, `ConicFlow`, `GenericConic`, `FermatParameters`, `DeformedEquation` (harmless). Tactics for the exact finite-field part: `ring` (for `H(curve) = 5·(...)`), explicit congruence lemmas `(14 : Polynomial K) = 4` etc. via `congrArg Polynomial.C` + `decide`, `fin_cases`, `norm_num`, `decide`. No `native_decide`.

**Facts, evidence, conjectures.**
- [Lean] Existence of a certificate for this exact `(F, 1, X0³Y0³)`.
- [CAS] Seed incidence, Jacobian unit mod 5, cover substitution identity, the two extra normal terms, the old-conic defect.
- [Obs] The only genuinely new mathematics relative to PR #128 is the seed (a *different-shape* curve: `x0` quadratic with free constant, `x2` quadratic, `S,P` quartic, 13 equations including degree 0) and the fact that lifting through `GenericIncidence` with the deformation coefficient fixed at `1/14` works; everything else is the shared library specialized to a different `H`. The `1/14` scaling (`H = 14G`) is chosen so that `H` is integral; `checks.json` and `solution.md` both insist the normal-pairing constant is rescaled by 14 in the residue calculation.
- [Conj] The same pipeline should work for **any** deformation of the pure-`X` ternary part `f(X)` of `F_hold` (keeping the `Y`-Fermat part and the two couplings), provided the weighted quotient `f(X) + P(αX0³X1 + βX1³X2) + S³ − 3SP²` stays irreducible and admits one nonsingular polynomial-curve incidence over some `F_p`; the exchange, deck trace and source identities are literally unchanged because `f` does not involve `y, T`. This is a conjecture about generic smoothness of the incidence, not a theorem in the repository.
- [Doc] The curation audit (`scripts/audit-cover-frontier.py`) is sympy-based and establishes "curation facts, not Lean acceptance."

**Relations.**
- `holdout-moving-branch` (ACCEPTED, `benediamond` PR #172): the complementary deformation `+X0²Y0²Y2²/7` **breaks** the exchange (`α²β²T⁴X0²` is not `y ↔ T` invariant) and was solved by the transitive-cubic Jacobian route instead; `why.md` of the present instance explicitly contrasts the two.
- `off-paired` (OPEN): added term `X0²X1X2Y0Y2/3 ↦ αβ X0²X1X2 T²/3`, not `y ↔ T` invariant, and the paired signs of the first two pairs are broken, so neither the cover nor the cubic-in-squares Jacobian routes apply as-is. [Obs] The ternary-deformation technique helps only after a serving quotient exists; it does not address the lost exchange.
- `one-reflection-skew` (OPEN): same obstruction (`why.md` says so explicitly and cites this instance as supplying reusable lifting/integration machinery only).
- `dense-rational` (OPEN): [Obs] the contrast is instructive — the ternary deformation is "dense in `X`" but has zero mixed deformation; `dense-rational` perturbs every bidegree, so no exchange survives.
- `rank-one-biquartic`, `cyclic-triple-residual`, `reflection-coupled-residual`, `quadric-contact-residual`: no direct relation beyond the shared étale-lifting/trace library.

---

## 3. The four q=1 specializations (shared mechanism)

All four have `F = F_hold`, `q = 1`, a single `Solution.lean` (3.7–4.1 KB, about 60 lines), `import HodgeTarget` and `import HodgeComputer.Constructions.HoldoutCover.GeneralizedCover`, `set_option maxRecDepth 10000`, `set_option maxHeartbeats 1000000`, and the same six local instances (`fermatCurve_adjoinRoot_isDomain`, `fermat_base_charZero`, `holdoutCoordinateRing_isDomain`, `genericCommonCoverRing_isDomain`, `quotient_adjoinRoot_isDomain`, `quotient_charZero`, `commonCoverFunctionField_finite`, plus a `CharZero CommonCoverFunctionField` instance). They were prepared with Codex assistance for `benediamond`, submitted from `fqlx`'s account, and claim no new conic; the library module `GeneralizedCover.lean` (11 KB, `maxHeartbeats 200000`) was added in a maintenance step and is credited to `benediamond`/Codex with `fqlx`'s construction underneath.

**Entry point** [Lean]:

```
HoldoutGeneralization.certificate_of_source_pair
  (A C : Poly) (B : MvPolynomial (Fin 5) (FermatCoefficientField ℂ)) (t : CommonCoverFunctionField)
  (hA : A.IsHomogeneous 6)
  (ht : Algebra.normalizedTrace HoldoutFunctionField CommonCoverFunctionField t = 0)
  (hsource  : actualCoverJacobian.det * quotientCoverFieldMap (ρq B) =
              physicalCoverFieldMap (ρp A) − physicalCoverFieldMap (ρp C) * t)
  (hpartial : same identity with integerPartial 0 applied to A, B, C) :
  HodgeComputer.HasCertificate holdoutSextic A 1
```

Inside: `HoldoutConic.exists_seed` (the central real-algebraic seed, `r ∈ (7,8)` root of the degree-23 `f`), `conic_target_residue_solution p B` (primitive for every polynomial source on the fixed weighted quotient), `residue_solution_piola_pullback`, `cover_normal_binding`/`cover_second_normal_binding`, `trace_difference` (kills the `C·t` term), `finite_residue_solution_clearing_value` with `holdoutSextic_irreducible Complex.I (I⁶ = −1)`, `hasCertificate_of_normalized` with `holdoutSextic_not_dvd_pderiv_zero`. Hence **coprimality over ℂ is by irreducibility of `F_hold` (monic specialization + Eisenstein) and `F ∤ D`**, not by `relPrime_of_specialization`. The denominator is again `(∂_0F)·d` homogenized; `Q` unconstrained; no degree exported.

**Target-to-source rule** [CAS, derived from `solution.md` §11 of `holdout` and verified on all five cases]. Write `A = M(X)·Y0^a Y1^b Y2^c`, `m := a + c + 1`, `e := b − m`. With `J = 2α⁵(T² − y²)/T`:

- correction `C = M(X)·Y0^a·Y2^c·Y1^m`, and `A − C·T^e = A·(1 − (y/T)^{m−b})`;
- the deck character of `T^e` is nontrivial iff `e ≢ 0 (mod 6)`; the library provides `trace_inverse_square` (`e = −2`), `trace_square` (`e = 2`), `commonCover_trace_inverse_fourth` (`e = −4`);
- the quotient source is `B = M(X)·α^{a−5}·β^{c}·P^{b}·(T^{m−b} − y^{m−b})/(2(T² − y²))`, which is a polynomial in `(X0,X1,X2,S,P)` with coefficients in the Fermat field when `m − b` is even: `m − b = 2 ⇒ B = M α^{a−5} β^c P^b/2`, `m − b = 4 ⇒ B = M α^{a−5} β^c P^b S/2`.

| instance | `(a,b,c)` | `e` | `t` | `C` (Lean `reflected`) | `B` (Lean `quotientSource`) | trace lemma |
|---|---|---|---|---|---|---|
| holdout-x04-p | (0,0,1) | −2 | `(coverT^2)⁻¹` | `X0⁴X2·Y1²·Y2` | `β/(2α⁵)·X0⁴X2` | `trace_inverse_square` |
| holdout-x0-x14-y0 | (1,0,0) | −2 | `(coverT^2)⁻¹` | `X0X1⁴·Y0·Y1²` | `1/(2α⁴)·X0X1⁴` | `trace_inverse_square` |
| holdout-x0-y03-y12 | (3,2,0) | −2 | `(coverT^2)⁻¹` | `X0·Y0³·Y1⁴` | `1/(2α²)·X0·P²` | `trace_inverse_square` |
| holdout-x12-p-y02 | (2,0,1) | −4 | `(coverT^4)⁻¹` | `X1²X2·Y0²·Y1⁴·Y2` | `β/(2α³)·X1²X2·S` | `commonCover_trace_inverse_fourth` |
| (holdout, central) | (3,0,0) | −4 | `(coverT^4)⁻¹` | `X0³Y0³Y1⁴` | `1/(2α²)·X0³·S` | `commonCover_trace_inverse_fourth` |

Both identities `J·B = A − C·T^e` and their `∂_{X0}`-derivatives were re-checked in sympy for all five rows (`True` throughout).

**Lean proof shape (identical in the four files).** `target_partial`, `reflected_partial`, `quotient_partial` by `norm_num [integerPartial_apply, <def>, MvPolynomial.pderiv_pow, MvPolynomial.pderiv_mul, Pi.single_apply, Fin.ext_iff] <;> ring`; `source_binding`/`partial_binding` by `rw [actualCoverJacobian_det]; simp only [map_mul, map_pow, …, quotientCoverFieldMap_residue_X, quotientCoverFieldMap_residue_C, physicalCoverFieldMap_residue_X, map_div₀]; change <explicit field identity>; field_simp [coverAlpha_ne_zero, coverT_ne_zero] <;> ring`; `certificate` by `change HasCertificate holdoutSextic target 1; apply certificate_of_source_pair target reflected quotientSource t` with homogeneity discharged by `simpa only [target, pow_one] using (isHomogeneous_X_pow … ).mul …`. For `holdout-x12-p-y02` the target contains no `X0`, so all three partials are `0` and `partial_binding` is `simp only [map_zero, mul_zero, zero_mul, sub_zero]`.

**Context from `why.md` (shared)** [Doc]: among the 58 paired-sign-invariant primitive sextic directions the direct quotient constructions support 48; the ten remaining split 6 + 4 into conjugate residue characters; `x04-p`, `x0-x14-y0`, `x0-y03-y12` lie in the 4-dimensional component, `x12-p-y02` is also outside the direct span; for each, some original-variety exchange (`X2 ↔ Y2`, etc.) fixes the numerator but reverses the volume form, so original-variety reflection quotients cannot serve them. The cover resolves this because the obstruction character is carried by the deck action on `T`, not by a reflection of the original variety. [Obs] This is the key transferable idea for the other "residual"/character-boundary instances.

**Geometric scope (not formalized)** [Doc, `holdout/solution.md` §"Geometric scope"]: for `F_hold`, every `Y`-odd target of degree `6q` for every `q ≥ 1` has a certificate (cover + physical diagonal quotient with weights `(0,0,0,5,1,5)` + deck/exchange descent). Of the 226 odd sextic monomials the cover criterion `e ≢ 0 mod 6` serves 158 (149 of the 208 primitive basis monomials); the complementary diagonal quotient serves the `e = 0` ones (47 basis monomials) by polynomial sources; twelve `e = −6` basis densities have fiber denominators and are only geometrically supported. The counts are reproduced by a short Python snippet in that file.

### 3a. holdout-x04-p

`A = X0⁴X2Y2`, degree 6, "transverse quartic". PR #140, discussion #53 (no comments). `(a,b,c) = (0,0,1)`, `e = −2`; `B = βX0⁴X2/(2α⁵)` (note the `α⁻⁵`: the source lives in the Fermat coefficient field, so negative powers of `α` are harmless). Partials: `∂_0A = 4X0³X2Y2`, etc. Related: `holdout-x12-p-y02` (same endpoint product `p = X2Y2`, conjugate residual character); the central holdout and its specializations; all same-variety instances are now accepted.

### 3b. holdout-x0-x14-y0

`A = X0X1⁴Y0`, "fourth-power X". PR #138, discussion #52. `(a,b,c) = (1,0,0)`, `e = −2`; `B = X0X1⁴/(2α⁴)`. Discussion #52 contains an exact **bounded exclusion** by `1337isnot1337` (Codex-assisted, reproduction code posted): no certificate with homogeneous `deg D ≤ 5` exists over ℂ, even with nonzero quotient — the `E1` matrix in degree 5 has 3024 columns, rank 2921, 13-dimensional denominator space, and the necessary `B_D(N) ∈ J(F)` equations span 65 independent quadrics in 91 monomials forcing all 13 coordinates to zero; degree 6 was left open (68 denominator coordinates, 40 forced to zero, Gröbner attempt timed out at 90 s); the specific denominator `X0X1X2Y0Y1Y2` fails `E1` at degree 6 (ranks 4437 vs 4438). The accepted proof is compatible with this (it exports no degree) and does not use it. Related: `holdout-x0-y03-y12` (same 4-dimensional residual component, degree shifted from even to odd coordinates); `holdout-x04-p`.

### 3c. holdout-x0-y03-y12

`A = X0Y0³Y1²`, "inverse character". PR #139, discussion #50 (no comments). `(a,b,c) = (3,2,0)`, `e = −2`; `B = P²X0/(2α²)` (the `Y1²` factor becomes `P²/T²·…` and is absorbed into the `(T²−y²)` cancellation; the middle odd coordinate is *not* an inert multiplier — it enters `B` through `P = yT`). Related: `holdout-x0-x14-y0`; conjugate character to the central cubic.

### 3d. holdout-x12-p-y02

`A = X1²X2Y0²Y2`, "mixed quadratic", linear in `p = X2Y2`. PR #142, discussion #48 (no comments). `(a,b,c) = (2,0,1)`, `e = −4`; `B = SX1²X2β/(2α³)`; trace lemma is the original `commonCover_trace_inverse_fourth`. Only instance of the four with `m − b = 4` (same as the central holdout), hence the `S` factor. Target has no `X0`, so the normal-derivative binding is vacuous. Related: `holdout-x04-p`; `holdout-all-six` (accepted, the six-variable monomial).

---

## 4. Cross-cutting observations

1. **Two different coprimality technologies.** Constructive/structured (holdout-q2): specialize to a univariate fiber where `F` is monic, prove coprimality there via a matrix Bezout identity `(M³+5I)K = I` in a finite algebra, and lift with `relPrime_of_specialization` + `relPrime_of_retract`. Nonconstructive (the other five): prove `F` irreducible once (monic specialization + Eisenstein, with `I⁶ = −1`), then `F ∤ D` follows from nonvanishing of `D` in `Frac(ℂ[Z]/F)`. The second is cheaper in Lean (no large determinant data) but gives no witness.
2. **The deck character replaces the original-variety reflection.** All four `q=1` specializations are targets fixed by an original reflection that reverses the volume; the cover's `μ6` deck action on `T` has character `ζ^e` on the correction `C·T^e` with `e = b − a − c − 1`, and `e ≢ 0 mod 6` is the exact serving criterion [Doc, consistent with the four cases]. This is the right lens for `reflection-coupled-residual`, `quadric-contact-residual`, `cyclic-triple-residual` (all OPEN, all "wrong-character on the same variety"): look for a finite cover on which the complementary character becomes a nontrivial deck character, rather than a new curve on the already-known quotient.
3. **Pole order two without new geometry.** holdout-q2 raises `q` on an existing contact/line construction using only (i) an explicit Jacobian primitive for the degree-5 source and (ii) a `q=1` solver for *arbitrary quadratic sources* on the cubic. Any future `q=2` attempt on a Jacobian model (e.g. `cycle-q2`) will need the analogue of (ii) on that model, which is the real cost; (i) is `ring`-level.
4. **What is *not* established by these six.** No denominator degree for any of them; no witness; no statement about other targets on the deformed ternary sextic; no formalization of the geometric fixed-variety theorem (the 158/68/47/12 counts are documentation only); the discussion #60 family of 34 `q=2` targets is unformalized; the degree ≤ 5 exclusion in #52 is a sympy computation, not Lean.
5. **Build hygiene notes** [Obs]: holdout-q2 sets `maxHeartbeats 0` in seven files and relies on `decide +kernel` over `List ℚ` with 60-digit rationals; the ternary deformation leaves `profiler true` on and uses `backward.isDefEq.respectTransparency false` globally. Neither affects soundness (kernel replay passed) but both matter if one wants to reuse these files as library code.

## 5. Independent checks performed today (all passed)

File `/tmp/hodge/check/group5_checks.py` (sympy 1.12, host) and `/tmp/hodge/check/quintic.sing` (Singular, `ddg-cas` container):

- source and `∂_{X0}` identities `J·B = A − C·T^e` for the five Fermat-cover targets;
- `cubic_substitution`, `sextic_good_fiber` (`F ↦ c⁶+5`), `sextic_reflection`, `target_decomposition`, good-fiber coefficients `(2,2,1,1)` — all modulo `k²+k+1`;
- `Σ M_i ∂_iG = R(R⁴ − 2916𝖺⁸𝖽⁶)·U³V²` with `R = 27𝖺²𝖻 + 4𝖼³` (sympy and Singular), value `515198970880` at `(2,2,1,1)`, `div M` quadratic;
- mod-5 seed: all 13 coefficients of `H(curve)` divisible by 5, degree 12, Jacobian `det ≡ 2`, equality with Lean's `explicitJacobian`, `explicitJacobian·inverseJacobian ≡ I`;
- `14F(cover) − H = 14T⁶(α⁶+β⁶−1)`; extra normal terms `(5X0⁴X1 + X2⁵)/14`, `10X0³X1/7`; the old-conic defect equals `checks.json`'s `originalConicDefect`.
