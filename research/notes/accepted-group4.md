# Accepted solves, group 4: the holdout family

Instances: `holdout`, `holdout-all-six`, `holdout-bundle`, `holdout-moving-branch`, `holdout-p3`, `holdout-p-y14`.
Repository: `/home/kejace/devel/sextics`, read-only. Variables are ordered `X0,X1,X2,Y0,Y1,Y2` = `X 0 … X 5` in Lean.

Legend used below: **[proved]** = Lean-verified in the accepted sources (axioms `propext`, `Classical.choice`, `Quot.sound`); **[evidence]** = exact computation recorded in the repo but not part of a `HasCertificate` theorem; **[unformalized]** = mathematical claim in prose only; **[failed/limit]** = approach recorded as not working; **[obs]** = my own observation while cataloguing.

---

## 0. Common ground for the five fixed-variety solves

### 0.1 The chain sextic and its symmetries

```
F_H = X0^6 + X1^6 + X2^6 + Y0^6 + Y1^6 + Y2^6 + X0^3 X1 Y0 Y1 + X1^3 X2 Y1 Y2,   q = 1
```

Five of the six instances use exactly this `F` (only `holdout-moving-branch` deforms it). Facts:

- Smooth **[unformalized in Lean; exact estimate in docs §2.7]**: a maximal-modulus coordinate contributes `6M^5` to its own partial, the mixed terms at most `4M^5`.
- Irreducible **[proved]**: `HoldoutConstruction.holdoutSextic_irreducible` (`data/holdout/HypersurfaceCoordinates.lean`). Argument: `F` is monic of degree 6 in `X0`; specialise `X1=X2=Y1=0`, `Y2=η` with `η^6=-1` (Lean uses `Complex.I`, `I^6=-1`) to get `X0^6+Y0^6-1`, Eisenstein at `Y0-1`; monicity preserves factor degrees. Not via `relPrime_of_specialization`/`relPrime_of_retract`.
- Symmetry group (physical): the three paired sign changes `(X_i,Y_i) -> (-X_i,-Y_i)`; the coordinate exchange `X2 <-> Y2` (fixes `A=X0^3Y0^3`, `X2^3Y2^3`, `X2Y1^4Y2`, `X0X1X2Y0Y1Y2`, but has determinant `-1`, so it reverses the residue volume: all these targets lie **outside** the reflection-invariant residue sector **[failed/limit]**); the diagonal action `g` of order 6 with weights `(0,0,0,5,1,5)` mod 6 on `(X,Y)`.
- Prior obstructions recorded in `why.md` **[failed/limit]**: separated-block formulas do not apply; direct reflection and diagonal-quotient support arguments serve 182 of the 208 primitive odd sextic directions, and all five targets here are among the ten missed directions inside the 58-dimensional paired-sign sector. These are limits of specified constructions, not nonexistence results.
- Weighted-fibration remark **[unformalized]**: over `[X0:Y0],[X1:Y1]` the fibre is a degree-6 surface in `P(2,2,1,1)` whose generic resolution is K3; not used in any proof.

### 0.2 The Fermat cover and its quotient (fqlx, PR #128)

Adjoin `T` with `T^6 = Y0^6 + Y2^6`; set `alpha = Y0/T`, `beta = Y2/T`, `y = Y1`. Then `alpha^6 + beta^6 = 1` (the coefficient curve; `K_b = Frac(C[alpha,beta]/(alpha^6+beta^6-1))`, irreducible by Eisenstein at `beta-1`). On the cover

```
G = X0^6+X1^6+X2^6 + y^6 + T^6 + y T (alpha X0^3 X1 + beta X1^3 X2),
```

which is symmetric under `sigma: y <-> T`. With `r: (y,T) -> (-y,-T)` and invariants `S = y^2+T^2`, `P = yT` (`y^6+T^6 = S^3 - 3SP^2`), the quotient is

```
H = X0^6+X1^6+X2^6 + P(alpha X0^3 X1 + beta X1^3 X2) + S^3 - 3 S P^2,
```

coordinates `(X0,X1,X2,S,P)` of weights `(1,1,1,2,2)`, weighted degree 6 (Lean: `weightedHypersurface`, irreducible: monic in `X0`, specialise `X1=X2=0,P=1` to `X0^6+S^3-3S`, Eisenstein at `S`). Normal chart `X0`: `¬ H ∣ ∂_0 H` and `¬ F ∣ ∂_0 F` by restriction to the `X0`-axis (`X0^6` vs `6X0^5`).

Physical-to-quotient Jacobian (base coordinate `beta`, distinguished `X0` fixed):

```
J_c = 2 alpha^5 (T^2 - y^2) / T          [proved: actualCoverJacobian_det]
```

Deck action on the common cover field, fixing all physical coordinates:

```
d: T -> zeta T,  alpha -> zeta^{-1} alpha,  beta -> zeta^{-1} beta,   zeta^6 = 1.
```

`e = sigma d sigma` acts by `y -> zeta y`; `d,e` commute, have order 6, `r = d^3 e^3`; the group generated with `sigma` is a quotient of `(mu_6 x mu_6) ⋊ Z/2` **[unformalized; prose in holdout solution.md]**. Lean only needs one `zeta` with `zeta^4 ≠ 1` (a primitive cube root of unity; `exists_nontrivial_sixth_root`) and the alg-hom `physicalDeckAlgHom`.

### 0.3 Target selection by deck character (the monomial criterion)

For an odd sextic monomial `A = M(X) Y0^a Y1^b Y2^c` put `m = a+c+1`, `e = b - m = b - a - c - 1`, `C = M(X) Y0^a Y2^c Y1^m`, and

```
B = M(X) alpha^{a-5} beta^c / 2 * (y^m T^b - T^m y^b)/(y^2 - T^2)   in K_b[X0,X1,X2,S,P].
```

Then `J_c B = A - C T^e` and (because `J_c`, `t = T^e` are `X0`-independent) `J_c ∂_0 B = ∂_0 A - (∂_0 C) T^e`. The correction `C T^e` is a `zeta^e`-eigenvector of `d`, so if `e ≢ 0 (mod 6)` its normalized trace to the physical field vanishes. Residue character of `A dvol` under the diagonal `g` is `5(a+1)+(b+1)+5(c+1) ≡ e (mod 6)`, so `e ≡ 0` is exactly the `g`-invariant sector.

Exponents realised by the five targets: `holdout` `e=-4`; `holdout-p3` `e=-4`; `holdout-all-six` `e=-2`; `holdout-p-y14` `e=+2`; the bundle's second summand `X0^3Y0Y1^2` has `e=0` (criterion fails: numerator of `B` vanishes). In the primitive basis (exponents ≤ 4) only `e ∈ {-6,-4,-2,0,2}` occur **[obs]**; the shared library supplies exactly the three needed trace lemmas: `commonCover_trace_inverse_fourth` (`T^{-4}`), `trace_inverse_square` (`T^{-2}`), `trace_square` (`T^2`).

Counts **[evidence; python snippet in holdout solution.md]**: 226 odd sextic monomials, 208 primitive (exponents ≤ 4, complement of the 18 Jacobian leading monomials `6 Z_i^5 Z_j`). Cover (`e ≢ 0`) serves 158 / 149 basis; diagonal sector 68 / 59 basis, split 47 with `e = 0` and 12 with `e = -6`. The `e=0` ones are polynomial sources on the complementary quotient of §0.5; the `e=-6` ones have a fibre pole there. Exchange `X2<->Y2` (an automorphism of `F_H`) carries some of them to served monomials, so ≥ 201 of 208 directions have a displayed finite-algebraic route.

### 0.4 Shared pipeline (all nonconstructive; how the denominator arises)

None of the six solves exports `D`, a degree, or `witness.json`. `proof.json` is `{schemaVersion 2, theorem "Submission.certificate"}` in every case. The denominator is **not** a norm from the cover; it arises as follows (`data/holdout/ResidueClearing.lean`, `FiniteResidueClearing.lean`, `HomogeneousIdeal.lean`, `HomogeneousField.lean`, `HomogeneousCertificate.lean`, all now mirrored in `HodgeComputer.Constructions.HoldoutCover.*`):

1. Residue (normal-chart) equation. With `g_i = ∂_i F`, tangent derivations `E_j = ∂_{j} - (g_j/g_0) ∂_0` (`j=1..5`) on `K_F = Frac(C[X,Y]/(F))`, and `r_A = ∂_0A/g_0 - A ∂_0^2F/g_0^2`, a primitive is `w ∈ K_F^5` with `div_{g_0}(w) := Σ_j (E_j w_j - w_j E_j g_0 / g_0) = -r_A`.
2. Trace descent: a primitive over a finite extension `L/K_F` (common cover field, or `E(t)` over the quotient) descends by `Algebra.normalizedTrace`, which commutes with the uniquely extended derivations (`FiniteDifferentials`, `FiniteTrace.residue_solution_descends_finite`) **[proved]**.
3. Legal clearing (`residue_solution_clearing`, `finite_residue_solution_clearing_value`): choose representatives `w_i = n_i/d` with `ρ(d) ≠ 0` (surjectivity of the quotient map; Classical choice), set
   `D = g_0 d`, `N_0 = A d - Σ_{i≥1} g_i n_i`, `N_i = g_0 n_i`. Then `N(F) = A D` exactly and `D div N - N(D) ∈ (F)`, so `A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q` with `Q` generally nonzero.
4. Coprimality over `C` **[proved]**: `legal_of_residue_ne_zero : Irreducible f → (∀p, ρ p = 0 ↔ f ∣ p) → ρ D ≠ 0 → D ≠ 0 ∧ IsRelPrime f D`. Since `ρ(g_0) ≠ 0` and `ρ(d) ≠ 0` in a field, `ρ(D) ≠ 0`; irreducibility of `F` over `C` is from §0.1. Nothing about `D` itself is ever computed.
5. Homogeneous extraction (`exists_legal_initial_representative`, `homogeneous_data_from_initial_denominator`, `homogeneous_data_of_normalized`): strip lowest homogeneous components of `D` that are multiples of `F` (`D_red = D + F U`), re-normalise `D' = g D_red`, `N' = g N + A F U e_0` (`g = ∂_0 F`); if `h` is the legal initial component of degree `e > 0`, pick `p > e+1`, `S_p = Σ_{i<p} h^i (h-D)^{p-1-i}` so `D S_p = h^p - (h-D)^p`, and extract the degree-`pe+1` components `M_i` of `S_p N_i`. Final data: `D = h^p` homogeneous of degree `pe`, `N_i` homogeneous of degree `pe+1`, `IsRelPrime F h^p`. Degree is existential.
6. `hasCertificate_of_normalized (F A D N Q) (hF : F.IsHomogeneous 6) (hA : A.IsHomogeneous 6) (hprime : Prime F) (hlegal : IsRelPrime F D) (normal : directional F N = A*D) (divergence : divergenceNumerator D N = F*Q) (i) (hchart : ¬ F ∣ pderiv i F) : HasCertificate F A 1`.

Denominator family **[obs]**: `D = (∂_0 F · d)`, homogenised to a power `h^p` of its legal initial form; `d` is a common denominator of a traced primitive and is never identified. Any attempt to extract an explicit witness has to make the finite field factor `E`, the trace, and `d` explicit; the repo explicitly declines to do so.

### 0.5 Shared library modules relevant to this group

- `HodgeComputer.Constructions.HoldoutCover.*` (58 modules = verbatim copies of `data/holdout/*.lean` minus `Solution.lean`), entry points `ConicResidue.conic_target_residue_solution` (primitive for **every** polynomial source `B ∈ K_b[X0,X1,X2,S,P]` on the fixed quotient `H`) and `GeneralizedCover.certificate_of_source_pair`:

  ```
  certificate_of_source_pair (A C : Poly) (B : MvPolynomial (Fin 5) (FermatCoefficientField ℂ))
      (t : CommonCoverFunctionField) (hA : A.IsHomogeneous 6)
      (ht : normalizedTrace HoldoutFunctionField CommonCoverFunctionField t = 0)
      (hsource : actualCoverJacobian.det * quotientCoverFieldMap (ρq B)
                  = physicalCoverFieldMap (ρp A) - physicalCoverFieldMap (ρp C) * t)
      (hpartial : same identity with integerPartial 0 applied to A, B, C)
      : HasCertificate holdoutSextic A 1
  ```
  It is hard-wired to `holdoutSextic` (not parametrised by `F`); the family, seed, étale lifting, trace and clearing are inherited from `holdout`.
- `HodgeComputer.CertificateOperations` (`Verify.add`, `HasCertificate.add`) extracted from `holdout-bundle`'s `BundleAssembly.lean`.
- `HodgeComputer.Constructions.PolynomialCurves` (`ArithmeticLifting.square_system_generic_lift`, `ComplexPoint.complex_point_of_etale_nontrivial`, `NormalizedCertificate`) used by `holdout-bundle` and `holdout-moving-branch`.
- `HodgeComputer.Constructions.TransitiveCubic` and `ResidueTransport.{MovingMap,Density,CoordinateIndependence}` used by `holdout-moving-branch`.

### 0.6 Catalogue status of the whole family (from folder inspection, `acceptance.json` presence)

All twelve same-variety instances on `F_H` are **solved**: `holdout` (#128), `holdout-all-six` (#135), `holdout-p-y14` (#136), `holdout-p3` (#137), `holdout-bundle` (#161), `holdout-q2` (q=2, `A = X0^3X1^2X2^2Y0^3Y1^2`, solved by a *strict constructive* endpoint-cubic/lines method, not by the cover), `holdout-x0-x14-y0`, `holdout-x0-y03-y12`, `holdout-x04-p`, `holdout-x1-y04-y1`, `holdout-x12-p-y02`, `holdout-x13-y13`. Nearby solved deformations: `holdout-ternary-deformation`, `mixed-biweight`, `mixed-polar`, `cycle-three`, `one-reflection` (= the equal-coupling `F_eq` of docs §2.7), `completed-square-nambu-holdout`, `holdout-moving-branch`. Open relatives: `off-paired` (deformation of `F_H`), and the non-holdout pairs listed in §7. `fermat` and `split` have no data folder.

---

## 1. `holdout` — PR #128, `fqlx`, accepted 2026-10-19 (discussion #44)

**(F,q,A)**: `F_H`, `q = 1`, `A = X0^3 Y0^3`.
**Method**: Fermat cover + moving weighted conic + polynomial integration + finite trace + deck character. **Nonconstructive**: no `D`, no degree, `Q` allowed nonzero.

### Mechanism (10 lines)
1. On the cover quotient `H` take the polynomial source `B = X0^3 S / (2 alpha^2)` (sign fixed by the `(S,P)` order; `det ∂(S,P)/∂(y,T) = 2(y^2-T^2)`, `det ∂(Y0,Y1,Y2)/∂(y,T,beta) = -T/alpha^5`). Then `J_c B = A (1 - y^4/T^4)`, i.e. `C = X0^3Y0^3Y1^4`, `e = -4`.
2. Exact complex conic seed on one fibre of the Fermat curve: real `r ∈ (7,8)` root of the degree-23 `f` (appendix of solution.md; `f(7) < 0 < f(8)`, IVT), `m = r^4-4r^3+5r^2+r-2`, `R = r(1-r-r^2)/(1-r)`, `M = r^2 m/(r-1)^2`, `a = 3(1-2r)/(1-r)`, `b = 3(1-2r)`; `q_c^4 = R`, `k_c^6 = M/q_c^2`, `p_c = r/q_c`, `alpha_c = a/(k_c q_c)`, `beta_c = b/(k_c^3 q_c)`; `f(r)=0 ⇔ a^6/(MR) + b^6/M^3 = 1 ⇒ alpha_c^6+beta_c^6 = 1`.
3. Conic: `X0 = s^2, X1 = k_c s t, X2 = t^2, S = -(s^4 + q_c^2 s^2 t^2 + t^4), P = p_c s^3 t + q_c s t^3`; `H(conic) ≡ 0`; passes through `(1,0,0,-1,0)` at `[1:0]`.
4. Incidence: on `s=1`, vary the 12 positive coefficients (`X0,X2` through degree 2; `S,P` through degree 4; `X1 = k_c τ` fixed) against the 12 equations `[τ^e] H(φ) = 0`, `e=1..12`. Parity splits the Jacobian into two 6×6 blocks whose determinants are nonzero scalars times `g_even(r)`, `g_odd(r)` (degree 11, positive for `r > 7` since `r = u+7` gives positive coefficients).
5. Étale algebra over the endpoint ring `(x0,0,x2,s0,p0)` with `x0^6+x2^6+s0^3-3s0p0^2 = 0`; localise at the Jacobian determinant; seed ⇒ nontrivial; flat base change to the endpoint function field `K_0`; finite field factor `E/K_0`; actual family `φ ∈ E[τ]^5`, `φ_1 = k_c τ`.
6. `K_H ↪ L = Frac(E[τ])` injective and finite (both trdeg 4 over `K_b`).
7. Polynomial normal `ψ ∈ E[τ]^5`, `g·ψ = 1`, `ψ(0) = e_0/(6x0^5)`, built from the same inverse Jacobian. `M = (ψ, ∂_τ φ, D_1 φ, D_2 φ, D_3 φ)`, `J = det M`, `J(0) = k_c/(6x0^5) ≠ 0`.
8. Tangent chain rule `E_j = Σ_i T_{ij} D_i` with `T = (M^{-1})` blocks; Piola determinant identity gives a polynomial density `P ∈ E[τ]`; integrate coefficientwise `Σ c_n τ^n ↦ Σ c_n τ^{n+1}/(n+1)`; divide by `J`; add the tangent correction. This solves `div_g(w) = -r_B` for every polynomial `B`.
9. Normalized trace `L → K_H`; then quotient → common cover field; Piola pullback with `J_c`; source `A - X0^3Y0^3Y1^4 T^{-4}`; deck trace kills the `zeta^{-4}` eigenvector; result: physical residue solution for `A`.
10. Legal clearing and homogeneous extraction (§0.4).

### Lean architecture
- 59 files, 13,569 lines; largest `CurveJacobianStaged.lean` (1256 lines, 74 KB) and `CurveGradientBlocks.lean` (1094 lines, 62 KB) — the two exact block-determinant calculations, kept separate "to control compilation memory".
- Heartbeats: **six files with `maxHeartbeats 0` (unbounded)**: `CoverSourceBinding`, `CurveJacobianScaling`, `ConicSourceDerivations`, `CurveJacobianStaged`, `CurveGradientBlocks`, `CurveSeedIncidence`; `CoverCertificate` 1,000,000; `ConicCoefficientField` 800,000; `EndpointFields` 600,000; `maxRecDepth 10000` (3 files) and `100000` (1 file).
- `Solution.lean` (27 lines):
  ```lean
  theorem Submission.certificate : HasCertificate HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q := by
    change HasCertificate holdoutSextic holdoutSource 1
    obtain ⟨p⟩ := HoldoutConic.exists_seed
    obtain ⟨w, hw⟩ := conic_target_residue_solution p coverQuotientSource
    exact holdout_hasCertificate_of_quotient_divergence w hw
  ```
- Lemma chain: `exists_seed` (`CurveSeed`: `intermediate_value_Icc`, `IsAlgClosed.exists_pow_nat_eq` for `q_c,k_c`) → `seedJacobian_nonsingular` (`CurveSeedIncidence`; `Matrix.det_fromBlocks_zero₂₁`, `det_fin_three`, `det_succ_row` expansions, `ring`/`linear_combination`) → `fermat_incidence_nontrivial` → `square_system_jacobian_etale`, `genericEtale_exists_fieldFactor` (`GenericIncidence`, uses `Mathlib.RingTheory.Smooth.Flat`, `Etale.Field`) → `conicCoefficientFactor_exists` → `conicCurve_equation/jacobian` → `conicTargetRingMap_injective` → `exists_polynomial_normal` → `conicEvaluationMatrix_det_ne_zero`, `conicSourceDerivations_commute` → `piola_divergence`, `normal_density_formula`, `residue_solution_of_polynomial_family`, `residue_solution_from_polynomial_curve` → `conic_target_residue_solution` → `residue_solution_descends_finite` → `holdout_hasCertificate_of_quotient_divergence` → `holdout_hasCertificate_of_quotient_residue` (`residue_solution_piola_pullback`, `cover_normal_source_binding`) → `holdout_hasCertificate_of_cover_residue` (`commonCover_trace_inverse_fourth`) → `finite_residue_solution_clearing_value` → `holdout_hasCertificate_of_normalized` → `hasCertificate_of_normalized`.
- Tactics for identities: `field_simp [coverAlpha_ne_zero, coverT_ne_zero] <;> ring`, `linear_combination`, `simp only [map_*]`, `MvPolynomial`/`Polynomial.funext` plus `norm_num only [coeff_*]` for coefficient extraction, `fin_cases <;> rfl` for matrix entries. No `decide` on complex data. Coprimality: §0.4 item 4 via `holdoutSextic_irreducible Complex.I hi`.
- Field-valued fraction derivation adapted from `HeadOfNaps`'s `data/quadric-contact-broken/FractionDerivation.lean`.

### Geometric scope section (prose only) **[unformalized]**
Fixed-variety theorem: every `Y`-odd target of degree `6q`, every `q ≥ 1`, has a certificate on `F_H`. Ingredients: generic fibre of `H` over the Fermat curve is a quasismooth degree-6 hypersurface in `P(1,1,1,2,2)` with three `½(1,1,1)` points at `X=0`, `S(S^2-3P^2)=0`, `-K = O(1)` (normal Q-Fano threefold ⇒ covering rational curves); physical diagonal `g` has a fixed point of age `5/6`; descent criterion docs §4.6 (`pd = p`, `pe = gp`, `pr = τp`) ⇒ `H^4(X,Q)^τ ⊆ N^1`. Contrast: this is exactly the `cycle-q2`-type question (higher `q`) for which no Lean route through the cover exists; `holdout-q2` was solved by a different strict method.

---

## 2. `holdout-all-six` — PR #135, `fqlx` (Codex for `benediamond`), accepted 2026-10-19 (discussion #51)

**(F,q,A)**: `F_H`, `q = 1`, `A = X0 X1 X2 Y0 Y1 Y2` (`a=b=c=1`, `m=3`, `e=-2`).
**Method**: target specialisation of the Fermat cover through `certificate_of_source_pair`. Nonconstructive.

Mechanism: `B = beta P X0 X1 X2 / (2 alpha^4)`, `C = X0 X1 X2 Y0 Y1^3 Y2`, `J_c B = A - C T^{-2}`; `∂_0` identity is nontrivial (`∂_0 A = X1X2Y0Y1Y2`, `∂_0 B = beta P X1 X2/(2 alpha^4)`); correction has character `zeta^{-2}`, `trace_inverse_square`. The `why.md` records **[failed/limit]** that the coordinatewise-squaring Jacobian density idea fails because the chain terms are not coordinatewise even, and that the paired-sign quotient and the endpoint reflection do not supply support (this monomial lies in the 4-dimensional character component conjugate to the central cubic's).

Lean (`Solution.lean`, 70 lines, `maxHeartbeats 1000000`, `maxRecDepth 10000`): imports `HodgeTarget` and `HodgeComputer.Constructions.HoldoutCover.GeneralizedCover` only; local instances for domains/char-zero; `def target/reflected/quotientSource`; `target_partial`, `reflected_partial`, `quotient_partial` by `norm_num [integerPartial_apply, pderiv_pow, pderiv_mul, Pi.single_apply, Fin.ext_iff] <;> ring`; `source_binding`/`partial_binding` by `rw [actualCoverJacobian_det]; simp only [map_*, quotientCoverFieldMap_residue_X/C, physicalCoverFieldMap_residue_X, map_div₀]; change …; field_simp [coverAlpha_ne_zero, coverT_ne_zero] <;> ring`; homogeneity by an `isHomogeneous_X_pow … .mul` chain; `certificate := certificate_of_source_pair target reflected quotientSource ((coverT^2)⁻¹) hA trace_inverse_square source_binding partial_binding`. Coprimality inherited from the library (§0.4).

---

## 3. `holdout-p3` — PR #137, `fqlx` (Codex for `benediamond`), accepted 2026-10-19 (discussion #46)

**(F,q,A)**: `F_H`, `q = 1`, `A = X2^3 Y2^3` (`a=0,b=0,c=3`, `m=4`, `e=-4`).
**Method**: same specialisation. Nonconstructive.

Mechanism: `B = beta^3 S X2^3 / (2 alpha^5)`, `C = X2^3 Y1^4 Y2^3`, `J_c B = A - C T^{-4}`, `commonCover_trace_inverse_fourth`. Since `A`, `B`, `C` do not involve `X0`, all three `∂_0` identities are `0 = 0` (`partial_binding` is `rw [quotient_partial, target_partial, reflected_partial]; simp only [map_zero, mul_zero, zero_mul, sub_zero]`). `why.md` notes the endpoint relation `(X2^6 - Y2^6)^2 = s^2 - 4p^6` (`s = X2^6+Y2^6`, `p = X2Y2`) organises the endpoint double cover but gives no primitive as the other coordinates vary **[failed/limit]**; reversing the oriented chain is not a symmetry, so the endpoint target is not a relabelling of the central one (the endpoint factors enter the coefficient curve through `beta`, not `alpha`).

Lean: 65 lines, identical skeleton and options to §2.

---

## 4. `holdout-p-y14` — PR #136, `fqlx` (Codex for `benediamond`), accepted 2026-10-19 (discussion #45)

**(F,q,A)**: `F_H`, `q = 1`, `A = X2 Y1^4 Y2` (`a=0,b=4,c=1`, `m=2`, `e=+2`).
**Method**: same specialisation. Nonconstructive.

Mechanism: `B = -beta P^2 X2 / (2 alpha^5)`, `C = X2 Y1^2 Y2`, `J_c B = A - C T^{2}` (check: `J_c B = -X2 beta y^2 T (T^2-y^2) = X2 Y1^4 Y2 - X2 Y1^2 Y2 T^2`), `trace_square`. `∂_0` identities trivial as in §3. `why.md` records that multiplying a known field by `Y1^4` is not passive because `Y1` sits in both chain terms **[failed/limit]**.

Lean: 65 lines, same skeleton; the correction argument to `certificate_of_source_pair` is `coverT^2` (not inverted).

---

## 5. `holdout-bundle` — PR #161, `fqlx` (ChatGPT collaborating agents), accepted 2026-10-20 (discussion #62)

**(F,q,A)**: `F_H`, `q = 1`, `A = X0^3 Y0^3 + X0^3 Y0 Y1^2`.
**Method**: certificate addition of (i) the accepted cover certificate for `X0^3Y0^3` and (ii) a new *complementary diagonal quotient* primitive for `A_* = X0^3 Y0 Y1^2` (`e = 0`: the cover's antisymmetrised source is identically zero **[failed/limit]**). Nonconstructive; both quotient remainders retained. Geometric existence was already known (docs ch. 4 §4.6) **[unformalized]**; the Lean contribution is the closed certificate.

### Mechanism (second summand)
1. Diagonal `g: (Y0,Y1,Y2) -> (zeta^5 Y0, zeta Y1, zeta^5 Y2)` preserves `F_H`; `A_*` has character `zeta`, `dvol` has `zeta^5`, so `A_* dvol/F^2` is invariant.
2. Invariant coordinates (with `i^2 = -1`): `p = Y0Y1`, `v = Y1Y2`, `lambda = Y1^3/(Y0^3 + i Y2^3)`; `a = lambda + lambda^{-1}`, `b = i(lambda - lambda^{-1})`. Relative hypersurface over `k = C(lambda)`:
   `K = X0^6+X1^6+X2^6 + X0^3X1 p + X1^3X2 v + a p^3 + b v^3`.
   For `theta = (X0,X1,X2,p,v,lambda)`: `K(theta) = F_H`, `J_theta = det(d theta) = -6 lambda Y1`, `rho = -X0^3 p/(6 lambda)`, `J_theta rho(theta) = A_*`. (Intermediate chart `p = Y0Y1, r = Y2/Y0, R = Y0^6`: `G = (1+r^6)R^2 + (S+Cp)R + p^6`, `G = RK`, `det ∂(p,v,lambda)/∂(p,r,R) = -p lambda/R`.)
3. Generic endpoint `(c0,0,c2,p0,v0)`, `c0^6+c2^6+a p0^3+b v0^3 = 0`, over `H = Frac k[c0,c2,p0,v0]/(…)`; curve `Φ0 = c0+u0 t+u1 t^2`, `Φ1 = t`, `Φ2 = c2+u2 t+u3 t^2`, `Φp = p0+u4 t+…+u7 t^4`, `Φv = v0+u8 t+…+u11 t^4`; equations `E_j = [t^j]K(Φ) = 0`, `j=1..12`.
4. Seed over `F_5` **[evidence → lifted, proved]**: `(a,b) = (0,3)`, endpoint `(1,1,4,1)`, `u = (1,0,4,1,4,2,2,4,4,1,1,2)`, `det Δ = 3`. The endpoint relation holds over `Z` at `(a,b) = (0,-2)`, realised by `lambda = i`; neither `a` nor `lambda^2+1` is inverted at the seed. `ComplexPoint.complex_point_of_etale_nontrivial` gives a complex nonsingular conic; étale flatness + injectivity `IntegralEndpoint → EndpointField` give nontriviality at the generic endpoint; a maximal ideal gives a finite separable `E/H` (`BundleIncidence.coefficientFactor`).
5. Polynomial normal `Ψ = Ψ^(0) - Σ (Δ^{-1} η)_ℓ V_ℓ`, `Ψ^(0) = e_0/(6c0^5)`; `M = (Ψ, D_t Φ, D_{c2} Φ, D_{p0} Φ, D_{v0} Φ)`, `J(0) = 1/(6 c0^5)`.
6. `L = Frac(k[X0,X1,X2,p,v]/K) ↪ E(t)` finite (trdeg 4; `Φ1 = t`, monic irreducible `K`).
7. Primitive `w_j = beta Ψ_j - (I/J) D_t Φ_j` with `I` the coefficientwise integral; `resDiv_{K_0}(w) = -r_rho`, `r_rho = p(4X0^3 - X1 p)/(18 lambda (2X0^3+X1 p)^2)`; trace to `L` (`quotient_polynomial_residue_solution` works for **any** polynomial density over `C(lambda)`).
8. Physical embedding `L ↪ K_F`: `Y1^6 = lambda(p^3 + i v^3)`, `Y0 = p/Y1`, `Y2 = v/Y1`; trdeg 5 ⇒ `(lambda,X1,X2,p,v)` algebraically independent (`physicalParameters_independent`, `physical_parameter_evaluation_injective`). Five-tail Jacobian with `lambda` as an active coordinate, determinant `-6 lambda Y1` (`physicalDiagonalJacobian_det`), `W_lambda = 0`, `W_0 = (rho - Σ K_j W_j)/K_0`, full Piola `V = J_theta (dθ)^{-1} (W∘θ)`.
9. Clearing/homogenisation as §0.4 (`physical_certificate_of_diagonal_residue`).

### Assembly **[proved]** (`BundleAssembly.lean`, now `CertificateOperations.Verify.add`)
`D = D0 D*`, `N_j = D* N_j^(0) + D0 N_j^(*)`, `Q = D*^2 Q0 + D0^2 Q*`; degrees add (`d+e`); `defect F (A+B) (D*E) 1 (E N + D M) = E^2 defect(A,D,N) + D^2 defect(B,E,M)`; legality `IsRelPrime.mul_right`.

### Lean architecture
- 29 files, 3,707 lines. Imports from `HodgeComputer.Constructions.HoldoutCover.{ResidueClearing, ResidueDerivatives, PolynomialFamilyBridge, PolynomialCalculus, AlgebraicFields, ResiduePiolaPullback, PolynomialNormal, PolynomialFieldMap, HypersurfaceDomain, FiniteTrace, FiniteDifferentials, FermatCoefficientDerivations, EvaluationMatrix, CoverSourceBinding, CoverCertificate, ConicResidue}` and `PolynomialCurves.{ComplexPoint, ArithmeticLifting}`.
- Heartbeats: 15 × 1,000,000; 7 × 800,000; 2 × 1,500,000; 1 × 1,600,000; **`BundleConicSeed.lean` 8,000,000** (12×12 `ZMod 5` Jacobian equality via `fin_cases <;> fin_cases <;> norm_num [...] <;> decide`, inverse via `explicitJacobian * inverseJacobian = 1 := by decide`); `maxRecDepth 100000` in 6 files; `set_option backward.isDefEq.respectTransparency false` throughout.
- Chain: `Solution.lean` (22 lines) → `HoldoutBundle.protected_certificate_of_bundle` (`BundleTarget`: `unfold; ring` identifications of `HodgeSubmission.F/A/q`) ← `bundle_certificate` ← `bundle_certificate_of_star first_summand_certificate star_certificate`; `first_summand_certificate` re-runs the holdout chain from the library; `star_certificate` ← `certificate_of_quotient_residue physicalDiagonalFieldMap w hw` (`BundleCompatibility`, structure `DiagonalFieldMap` bundling `fieldMap`, `tangent_compatible`, `normal_binding`, `source_binding`) ← `quotient_density_residue_solution` ← `quotient_polynomial_residue_solution quotientDensity`.
- Coprimality: `F` irreducible from the shared `holdoutSextic_irreducible`; `K` and the endpoint polynomial irreducible by monic + specialisation (`quotientAxis_equation : quotientAxis kPoly = X^6`, `endpointPolynomial_irreducible`); legality via `legal_of_residue_ne_zero` through `physical_certificate_of_residue`.

---

## 6. `holdout-moving-branch` — PR #172, `benediamond` (Codex), accepted 2026-10-21 (no discussion URL; `sameVarietyInstances = []`)

**(F,q,A)**: `F = F_H + (1/7) X0^2 Y0^2 Y2^2`, `q = 1`, `A = X0^3 Y0^3`. This is a **deformation** of the holdout variety, not the same variety.
**Method**: transitive-cubic Jacobian residue transport (fqlx's mixed-biweight method) with a new `F_5` incidence seed. Nonconstructive.

### Why the holdout machinery fails here **[failed/limit; `why.md`, `checks.json`]**
- On the Fermat cover the perturbation is `(1/7) X0^2 alpha^2 beta^2 T^4`; `y <-> T` changes `G` by `(1/7) X0^2 alpha^2 beta^2 (y^4 - T^4)` (`originalCoverSwapDifference`), so `(S,P)` no longer describe the equation. Absorbing it into `T` makes the coefficient curve `alpha^6 + beta^6 + (1/7)(X0/T)^2 alpha^2 beta^2 = 1` depend on an active coordinate.
- The `(5,1,5)` diagonal gives the perturbation weight 2 mod 6 ⇒ `g` no longer preserves `F`.
- `checks.json` (`scripts/audit-cover-frontier.py`): smooth by the maximal-coordinate estimate with per-variable bounds `(23/7, 4, 1, 9/7, 2, 9/7)`, strict bound 6; `hessianSymmetrizerDimension 1`; `coordinateReflections []` (the `X2<->Y2` exchange is broken); diagonal group order 72, target-fixing subgroup order 24, no junior elements fixing the target; `JacobianRanks [36, 37]`; no linear direct-sum decomposition. All are boundaries of specified screens, not nonexistence.
- Surviving symmetry: all three paired sign changes.

### Mechanism (10 lines)
1. Paired-sign chart: `R = X2^2`, `v = Y0/X0`, `c_i = (Y_i/X_i)/v`, `x = X0^2/R`, `y = X1^2/R`. Then `F = R^3 U(x,y,1)` with the transitive cubic `U = a x^3 + b y^3 + c z^3 + d x^2 y + e y^2 z + f x^2 z`,
   `(a,b,c) = (1+v^6, 1+c1^6 v^6, 1+c2^6 v^6)`, `(d,e,f) = (c1 v^2, c1 c2 v^2, c2^2 v^4/7)`; target `A = v^3 R^3 x^3`. (`mixed-biweight` has `(d,e,f) = (c1 v^2/2, c1c2 v^2/2, c2 v^4/3)` with the same diagonal.)
2. Library covariants `H, Θ, J`; on `U = 0`: `xi = Θ/H^2`, `eta = -J/(2H^3)`, `eta^2 = xi^3 + s4 xi - s6/4`, invariant differential `D_y xi/(2 eta) = 3/U_x` (`TransitiveCubic.Map`, `Syzygy`); nonvanishing of `H` etc. proved by evaluation at a point on an explicit localisation.
3. Jacobian model: `T = R v^6`, `X = R^2 xi`, `Y = R^3 eta`; `V = Y^2 - X^3 - A4 X - B6`, `A4 = R^4 s4` (degree 4 in `(R,T)`), `B6 = -R^6 s6/4` (degree 6); integral form `111132 V = 0` with `a8 := 4116 A4`, `b256 := 111132 B6` (`JacobianModel.lean`, weights proved by `ring`).
4. Seed over `F_5` at `c1 = c2 = 1` **[evidence → lifted, proved]**: `R = 1 + 2T^2`, `X = 4 + 4T + T^2 + 3T^3 + T^4`, `Y = 1 + 3T^2 + 4T^4 + 2T^6`; fix `R(0) = 1`, `X(0) = 4`; 13 unknowns, 13 incidence equations; explicit 13×13 inverse (`inverse_product … := by decide`); `JacobianSeedBinding` identifies the matrix with `∂(equations)` (`fin_cases <;> norm_num […] <;> decide`).
5. Arithmetic lifting (`PolynomialCurves.ArithmeticLifting.squareGeneric_exists_fieldFactor`): parameters `c1, c2, R(0), X(0)` stay independent; generic field `Frac C[c1,c2,R0,X0]`; maximal ideal ⇒ finite separable factor with invertible Jacobian (`complex_generic_factor_exists`).
6. Polynomial normal, independent endpoint derivations, curve parameter ⇒ dominant family (`CoordinateIndependence`); polynomial integration (`PolynomialFamilyBridge`); normalized trace ⇒ `target_residue_solution (A : Ambient)` for every polynomial source on the Jacobian hypersurface.
7. Source: vary `a = alpha + v^6` at `alpha = 1`; `source = -V_alpha/(144 R)` is polynomial (`alphaModel.derivative.eval 1 = -16003008 · R · source`, i.e. `-27783 · R · (576 source)`); `∂_0 source = 0`.
8. Physical pullback (`PhysicalJacobianMap`): injective `coefficientMap`, `fieldMap`, `tangent_compatible`; residue-volume multiplier `(gradient 0/(2 targetY)) · det = -144 R v^3` (`density_multiplier`).
9. Corrections (`PhysicalPrimitive`): moving-map correction `(xiDot/(2 eta)/3) · volumeVector` (`ResidueTransport.MovingMap`, `PhysicalCubicVolume.correction_divergence`) and normal-chart correction `correctionVector = Pi.single 2 (-X0^3 Y0^4/pairedGradient)` (`ResidueTransport.NormalChange`-style `paired_normal_change`); `w = pulled + correction - correctionVector`; closed by `linear_combination ratio^3 * source_transport`.
10. `certificate_of_residue_solution` (`PhysicalClearing`): shared `residue_solution_clearing` + local `certificate_of_normalized`.

### Lean architecture
- 31 files, 3,298 lines; `Solution.lean` (20 lines) rewrites `HodgeSubmission.F = CyclePhysical.sextic` and `HodgeSubmission.A = CyclePhysical.source` by `unfold …; ring`, then `exact CyclePhysicalPrimitive.certificate`. Namespaces are inherited from `cycle-three`/`mixed-biweight` (`CyclePhysical*`, `CycleJacobian*`).
- Imports: `TransitiveCubic.{Cubic, Functorial, Map, Syzygy, Diagonal}`, `ResidueTransport.{MovingMap, Density, CoordinateIndependence}`, `PolynomialCurves.ArithmeticLifting`, `HoldoutCover.{ResidueClearing, ResiduePiolaPullback, ResidueDerivatives, PolynomialNormal, PolynomialFieldMap, PolynomialFamilyBridge, PolynomialCalculus, HypersurfaceDomain, HomogeneousCertificate, FiniteTrace, FiniteDifferentials, EvaluationMatrix, AlgebraicFields}`; `Mathlib.Algebra.TrivSqZeroExt.Basic` (dual-number jets).
- Heartbeats: **3 × 8,000,000** (`JacobianIncidence`, `JacobianSeed`, `JacobianSeedBinding`), 1 × 4,000,000, 3 × 3,000,000, 7 × 2,000,000, 6 × 1,000,000; `maxRecDepth 100000` ×2; `backward.isDefEq.respectTransparency false` throughout.
- Coprimality: `PhysicalField.sextic_irreducible` (monic in `X0`, `equation_monic`; specialisation `equation_specialization` ⇒ irreducible; `MvPolynomial.finSuccEquiv`), `normal_legal : ¬ sextic ∣ pderiv 0 sextic` (axis restriction `X^6` vs `6X^5`), then `legal_of_residue_ne_zero` inside `residue_solution_clearing`. Coordinate nonvanishing `coordinate_ne_zero` by a total-degree argument (`totalDegree_le_of_dvd_of_isDomain`).
- Identity tactics: `field_simp; linear_combination he` (`cubic_equation`), `simp […]; ring` for the `alpha`-derivative, `decide` for `ZMod 5` matrices.

---

## 7. Relations to open instances

| Open instance | (F,q,A) relation to this group | What transfers / what blocks |
| --- | --- | --- |
| `off-paired` (OPEN, disc. #57) | `F = F_H` with couplings `1/2` plus `(1/3) X0^2 X1 X2 Y0 Y2`; `A = X0^3Y0^3`. Deformation of `holdout`. | **[failed/limit, why.md]** Added term becomes `(1/3) alpha beta X0^2 X1 X2 T^2` on the Fermat cover: not `y<->T` symmetric, so no descent to `(S,P)`; weight 4 under the `(5,1,5)` diagonal, so the complementary quotient is lost too; it is **odd** under the paired sign change of `(X0,Y0)` (degree 3 in that pair), so the paired-sign cubic reduction used by `cycle-three`, `mixed-biweight`, `holdout-moving-branch` does not give a cubic in `X_i^2`. Survivors: product of the first two paired signs and the third paired sign. `holdout-moving-branch` is the closest solved analogue (cover broken, paired signs intact → Jacobian route works); `off-paired` breaks the Jacobian route one step earlier. Needs a replacement cover/quotient or a method insensitive to sign characters. |
| `one-reflection-skew` (OPEN) | deforms `one-reflection` (= `F_eq`, equal couplings `Y0Y1(X0^3X1+X1^3X2)`, served by the same Fermat cover with `DiamondJM`'s seed, PR #144) by `± (1/7) X0X1 Y0^{3}Y1 / X0X1Y0Y1^3`. | **[obs]** On the cover (`Y0 = alpha T`, `Y1 = y`) the skew pair is `(alpha X0X1/7) · yT · (alpha^2 T^2 - y^2)`, not `y<->T` symmetric, so the §2.7 descent fails unchanged; not checked further here. |
| `cycle-q2` (OPEN, q=2 on `cycle-three`'s variety) | analogue of `holdout-q2` on the cyclic deformation. | `holdout-q2` was solved by a strict endpoint-cubic/line construction, **not** by raising the cover construction to `q=2`; the cover method's `q ≥ 2` extension exists only as the unformalized all-pole theorem (docs §3.6/§4.5). No Lean template for `q=2` via trace/clearing exists in this group. |
| `cyclic-triple-residual`, `reflection-coupled-residual`, `quadric-contact-residual` (OPEN; `A = X0^3Y0^3` on varieties whose "fano" sibling is solved) | different varieties; same target character pattern (residual sector not reached by the sibling's quotient). | Structural lesson from `holdout-bundle`: a residual character can be reached by a *second* quotient and added with `HasCertificate.add`; the `e`-character bookkeeping (§0.3) is the holdout-specific instance of that idea. |
| `rank-one-biquartic`, `dense-rational` (OPEN, `A = X0^3Y0^3`) | dense sextics, no paired signs or exchange. | Nothing here transfers directly; all six proofs depend on sparse symmetric structure (exchange on a cover, a diagonal character, or paired signs). |

Same-variety open instances: **none** — all twelve `F_H` targets in the catalogue are solved (§0.6).

---

## 8. Observations and caveats **[obs]**

1. All six are existence proofs; the only explicit objects are the seeds (one real algebraic number of degree ≤ 23 plus 4th/6th roots for `holdout`; `F_5` points for the bundle and moving-branch). The generic conic is a `Classical.choose` of a maximal ideal; the degree `[E : K_0]` and the trace are never computed. Extracting a witness would need elimination of a 12- or 13-variable square system over a function field of trdeg 3–4 plus a trace — the repo explicitly declines this.
2. Heartbeat profile: the hot spots are exact determinant/inverse checks (`maxHeartbeats 0` in six holdout files; `8,000,000` for the three `ZMod 5` matrix files). Any reuse that changes the seed re-pays that cost.
3. `GeneralizedCover.certificate_of_source_pair` is target-linear but `F`-rigid: it binds `holdoutSextic`. A deformation that keeps the `y<->T` exchange (as `holdout-ternary-deformation` did, with a new 13-coefficient seed by `hodgepodgeconjecture`, PR #151) must rebuild the cover fields and the conic; one that breaks the exchange (`moving-branch`, `off-paired`) cannot use it at all.
4. The bundle's `quotient_polynomial_residue_solution` serves every polynomial density over `C(lambda)` on `K`, hence in principle all 47 `e=0` basis directions on `F_H` with only new source bindings; the 12 `e=-6` directions (fibre pole `1/(p^3 + i v^3)`) are outside it. Archival now that the family is closed, but it is the template for "second quotient + add" on other varieties.
5. The `moving-branch` solve is the only one in this group whose variety is not `F_H`; it documents that the transitive-cubic Jacobian route survives an even perturbation `X0^2Y0^2Y2^2/7` as long as the three paired signs survive and `(d,e,f)` stay in the six-coefficient family. The polynomial-source scope (docs §2.9 "Coefficient and source scope") gives 22 paired-sign-invariant, diagonal-invariant odd monomials with polynomial Jacobian sources in that family — a concrete list for other deformations with intact paired signs.
6. No `witness.json`, no `relPrime_of_specialization`, no `decide` on complex data, and no `MvPolynomial.funext`-style brute-force check of the certificate identity appears anywhere in these six; all identity checking is in function fields via `field_simp`/`ring`/`linear_combination`, and all coprimality is "irreducible `F` + nonzero image in `Frac(C[X,Y]/F)`".
