# Accepted solves, group 6: holdout-x13-y13, holdout-x1-y04-y1, klein-hesse, klein-rational, m3-coupled, mixed-biweight

Catalogue notes prepared from the repository at `/home/kejace/devel/sextics` (branch `main`, head `83fbf5c`), read-only. No repository file was modified; nothing was posted or pushed.

## 0. Conventions, sources, and what was verified

**Certificate convention** (from `formal/README.md` and `formal/lean/HodgeComputer/Verify.lean`). Variables `Z = (X0,X1,X2,Y0,Y1,Y2)` = `MvPolynomial.X 0 … X 5`. For `(F,q,A)` a certificate is `d, D, N : Fin 6 → Poly, Q` with

- `D ≠ 0`, `D` homogeneous of degree `d`, each `N_i` homogeneous of degree `d + 6q − 5`, `IsRelPrime F D`,
- `defect F A D q N := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^(q+1)·Q`.

`normalized_identity` reduces this to the pair `E1: q·Σ N_i ∂_iF = A·D`, `E2: D·div N − N(D) = F^q·Q`. The two coprimality tools in `Verify.lean` are

- `relPrime_of_specialization (φ : R →+* S) (f g : Polynomial R) (hf : IsUnit f.leadingCoeff) (h : IsRelPrime (f.map φ) (g.map φ)) : IsRelPrime f g` (univariate; needs only a unit leading coefficient, not disjointness of the hypersurfaces), and
- `relPrime_of_retract (f : R →+* S) (g : S →+* R) (hg : g.comp f = id) (h : IsRelPrime (f a) (f b)) : IsRelPrime a b`.

**Status labels used below.** PROVEN = established by the accepted, kernel-replayed Lean proof (axioms `propext`, `Classical.choice`, `Quot.sound`; Lean 4.33.0). COMPUTED = my independent check with sympy 1.12 (`/tmp/hodge/scratch/verify_group6.py`, `/tmp/hodge/scratch/klein_extra.py`; the `ddg-cas` docker image and `/usr/bin/Singular` are also present but were not needed). OBSERVATION / CONJECTURE = my reading, not proven anywhere. NOT CLAIMED = explicit non-claims in the accepted explanations.

**Files read per instance:** `instance.rec`, `why.md`, `solution.md`, `proof.json`, `provenance.json`, `acceptance.json`, `witness.json` (where present), `Solution.lean`, and all helper `.lean` files (structure/theorem names; full text of the load-bearing modules). Shared library: `formal/lean/HodgeComputer/Constructions/HoldoutCover/{README.md,GeneralizedCover.lean,ResidueClearing.lean,HomogeneousCertificate.lean,FiniteTrace.lean}`, `Constructions/TransitiveCubic/README.md`, and `data/holdout/solution.md`, `data/cycle-three/solution.md` (the two parent constructions).

**Summary table.**

| instance | q | A | method | D | Q | PR / author / accepted | files |
|---|---|---|---|---|---|---|---|
| holdout-x13-y13 | 1 | `X1^3 Y1^3` | nonconstructive, Fermat cover + moving conic (shared `HoldoutCover.GeneralizedCover`) | existential (no degree) | nonzero allowed | #143, `fqlx` (Codex for `benediamond`), 2026-10-19 | `Solution.lean` 65 lines |
| holdout-x1-y04-y1 | 1 | `X1 Y0^4 Y1` | same | existential | nonzero allowed | #141, `fqlx`, 2026-10-19 | `Solution.lean` 65 lines |
| klein-hesse | 1 | `(X0^3X1^2 − X1^3X2^2) Y0` | constructive, Klein→Hesse monomial map + Piola, strict | deg 15, 10 terms, X-only, irreducible/Q | 0 | #22, `benediamond`, 2026-10-15 | `Witness.lean` 385 lines (31.7 KB) |
| klein-rational | 1 | `X0^3X1^2 Y0` | constructive, Fourier (Q(ω)) channels + rational Piola trace, strict | deg 30, 28 terms, X-only, = norm from Q(ω) | 0 | #42, `benediamond`, 2026-10-16 | `Witness.lean` 460 lines (66 KB) |
| m3-coupled | 1 | `X0X1X2Y1(4X0^2+5X1^2−5Y0^2)` | constructive, conic in quadratic coordinates + Cramer + Euler correction, strict | deg 10, 190 terms, = −2·(irreducible/Q) | 0 | #21, `benediamond`, 2026-10-15 | `Witness.lean` 1231 lines (369 KB) |
| mixed-biweight | 1 | `X0^3 Y0^3` | nonconstructive, paired-sign transitive cubic → Weierstrass Jacobian primitive (adapts cycle-three) | existential | nonzero allowed | #170, `fqlx` (ChatGPT/Codex), 2026-10-21 | 59 helper modules + `Solution.lean` (≈2.2 MB, 15.8k lines) |

All six have `theorem: Submission.certificate`; the three constructive ones also export `witnessTheorem: Submission.witness` with `witnessChecked: true`.

---

## 1. holdout-x13-y13 (Holdout · second cubic)

**Statement.** `F = X0^6 + X1^6 + X2^6 + Y0^6 + Y1^6 + Y2^6 + X0^3 X1 Y0 Y1 + X1^3 X2 Y1 Y2` (the holdout chain), `q = 1`, `A = X1^3 Y1^3`. Same variety as `holdout` and ten other holdout-* instances (all solved). Symmetry noted in `why.md`: `A` is invariant under all three paired sign changes `(X_i,Y_i) ↦ (−X_i,−Y_i)` and under the endpoint exchange `X2 ↔ Y2` (which reverses the residue volume, so that quotient serves the opposite sector); the direct paired-sign quotient leaves this target among ten independent directions outside its 48-dimensional supported subspace of the 58-dimensional paired-sign sector (prior work, quoted from `why.md`). Discussion #49.

**Method.** Nonconstructive specialization of `fqlx`'s Fermat-cover / moving-conic construction (PR #128), through the target-parameterized entry point `HoldoutGeneralization.certificate_of_source_pair` in `HodgeComputer.Constructions.HoldoutCover.GeneralizedCover`. No denominator, degree, or quotient is exhibited.

**Mechanism (PROVEN in the shared library; target data PROVEN in `Solution.lean`).**
1. Cover: adjoin `T^6 = Y0^6 + Y2^6`, put `α = Y0/T`, `β = Y2/T`, `y = Y1`; coefficient field `K_b = Frac(C[α,β]/(α^6+β^6−1))` (Fermat curve). On the cover `G = X0^6+X1^6+X2^6+y^6+T^6 + yT(αX0^3X1 + βX1^3X2)` is invariant under `y ↔ T`; invariants `S = y^2+T^2`, `P = yT` give the weighted quotient `H = X0^6+X1^6+X2^6 + P(αX0^3X1+βX1^3X2) + S^3 − 3SP^2`, weights `(1,1,1,2,2)` on `(X0,X1,X2,S,P)`.
2. The moving-conic theorem `HoldoutConicIncidence.conic_target_residue_solution` gives a residue primitive on `H = 0` for every polynomial source `B ∈ K_b[X0,X1,X2,S,P]` (exact nonsingular conic seed, finite algebraic family through the generic point, polynomial integration, normalized finite trace from the family's coefficient extension).
3. Piola transport back to the cover uses the physical-to-quotient Jacobian `J = 2α^5 (T^2 − y^2)/T`. The target-specific inputs are the two exact identities (both discharged in `Solution.lean` by `rw [actualCoverJacobian_det]; simp only […]; change …; field_simp [coverAlpha_ne_zero, coverT_ne_zero] <;> ring`):

   `J·B = A − C·T^2` with `B = −P X1^3/(2α^5)`, `C = X1^3 Y1`, and the same identity for `∂/∂X0` (trivial here: `A`, `B`, `C` are all `X0`-free, so `target_partial`, `reflected_partial`, `quotient_partial` are `0 = 0`).
4. Deck trace: the common cover field has an automorphism fixing the physical coordinates with `T ↦ ζT`, `α ↦ ζ^{-1}α`, `β ↦ ζ^{-1}β` (`ζ^6 = 1`, `ζ^4 ≠ 1`). The correction `C·T^2` is an eigenvector with eigenvalue `ζ^2 ≠ 1`, so its normalized trace vanishes (`HoldoutGeneralization.trace_square`), leaving the residue equation for `A` itself (`trace_difference`).
5. Shared closure: `residue_solution_clearing` (choose polynomial representatives `d, n_i` of the physical residue solution `w_i ∈ Frac(C[Z]/(F))`; set `D = ∂0F·d`, `N0 = A d − Σ_{i≥1} g_i n_i`, `N_i = ∂0F·n_i`; irreducibility of `F` plus `ρ(D) ≠ 0` give `IsRelPrime F D`) then `hasCertificate_of_normalized` (graded reduction `D_red = D + FU`, then the power trick `H_p = h^p`, `S_p = Σ h^i (h−D)^{p−1−i}`, extracting the degree-`pe+1` component of `S_p N`) produce homogeneous legal data with `N(F) = AD` and `B_D(N) = FQ`. This is where the denominator "arises": it is `∂0F` times a common denominator of trace-descended representatives, then homogenized by a power; it is never written down.

**Character bookkeeping (COMPUTED/OBSERVATION, from `data/holdout/solution.md` §"Target characters").** For `A = M(X) Y0^a Y1^b Y2^c`, put `m = a+c+1`, `e = b − m`, `C = M(X) Y0^a Y2^c Y1^m`; the unnormalized invariant sum has density `A − C·T^e` and the trace removes the correction iff `e ≢ 0 (mod 6)`. Here `(a,b,c) = (0,3,0)`, `m = 1`, `e = 2`. I verified `J·B = A − C T^2` symbolically (sympy) and that the exponent of the lemma used (`trace_square`) matches `e = 2`.

**Lean architecture.** `Solution.lean` (65 lines, 3.4 KB): `import HodgeTarget`, `import HodgeComputer.Constructions.HoldoutCover.GeneralizedCover`; `set_option maxRecDepth 10000`, `set_option maxHeartbeats 1000000`; local instances `fermatCurve_adjoinRoot_isDomain`, `fermat_base_charZero`, `holdoutCoordinateRing_isDomain`, `genericCommonCoverRing_isDomain`, `quotient_adjoinRoot_isDomain`, `quotient_charZero`, `commonCoverFunctionField_finite`, plus a `CharZero CommonCoverFunctionField` instance. Declarations: `target := (X 1)^3*(X 4)^3`, `reflected := (X 1)^3*(X 4)`, `quotientSource : MvPolynomial (Fin 5) (FermatCoefficientField ℂ) := C (−β^0/(2α^5)) * (X 1)^3 * (X 4)` (coordinates `(X0,X1,X2,S,P)`, so `X 4 = P`); `target_partial`, `reflected_partial`, `quotient_partial` (`norm_num [integerPartial_apply, …, pderiv_pow, pderiv_mul, Pi.single_apply, Fin.ext_iff] <;> ring`); `source_binding`; `partial_binding` (rewrites with the three zero partials); `certificate` by `change HasCertificate holdoutSextic target 1; apply certificate_of_source_pair target reflected quotientSource (coverT^2)` discharging homogeneity (`isHomogeneous_X_pow … .mul …`), `trace_square`, `source_binding`, `partial_binding`. Coprimality over ℂ is entirely inside the library (irreducibility of the holdout sextic by monic/Eisenstein specialization, `legal_of_residue_ne_zero`). No helper modules; `witnessChecked: false`.

**Relations.** Same variety: `holdout`, `holdout-p3`, `holdout-p-y14`, `holdout-x12-p-y02`, `holdout-x0-y03-y12`, `holdout-all-six`, `holdout-x0-x14-y0`, `holdout-x04-p`, `holdout-bundle`, `holdout-q2` (all solved). Open deformations of the chain: `off-paired` (adds `X0^2X1X2Y0Y2/3` with chain coefficients `1/2`; breaks `y ↔ T`), `one-reflection-skew` (deformation of the different chain `one-reflection`; the same cover `T^6 = Y0^6+Y2^6` loses its exchange). `mixed-biweight` (solved) is the chain-with-`1/2`-coefficients plus the opposite biweight `X0X2Y0^3Y2/3`, where this cover fails and a different method was needed (see §6).

**NOT CLAIMED.** New conic construction, optimized witness, denominator degree, strictness.

---

## 2. holdout-x1-y04-y1 (Holdout · fourth-power companion)

**Statement.** Same `F` as §1, `q = 1`, `A = X1 Y0^4 Y1`. Invariant under each paired sign change and under `X2 ↔ Y2`; again outside the 48-dimensional directly supported subspace (`why.md`). Discussion #47.

**Method / mechanism.** Identical pipeline to §1 with the target data

`B = P·S·X1/(2α)`, `C = X1 Y0^4 Y1^5`, exact identity `J·B = A − C·T^{−4}`,

`(a,b,c) = (4,1,0)`, `m = 5`, `e = −4`, correction character `ζ^{−4} = ζ^2 ≠ 1`, trace lemma `HoldoutConstruction.commonCover_trace_inverse_fourth` (the one used by the original holdout proof, since `X0^3Y0^3` also has `e = −4`). COMPUTED: `J·B = A − C T^{−4}` holds symbolically after substituting `Y0 = αT`, `Y1 = y`. Again `A, B, C` are `X0`-free so the partial identity is `0 = 0`.

**Lean architecture.** `Solution.lean` (65 lines, 3.6 KB) is line-for-line the §1 template with `target := X 1 * (X 3)^4 * X 4`, `reflected := X 1 * (X 3)^4 * (X 4)^5`, `quotientSource := C (β^0/(2α^1)) * X 1 * X 3 * X 4` (i.e. `S·P·X1/(2α)`), `t := (coverT^4)⁻¹`, and `exact commonCover_trace_inverse_fourth`. Same options (`maxHeartbeats 1000000`, `maxRecDepth 10000`), same local instances. PR #141.

**Holdout family cross-check (COMPUTED by grepping every `holdout-*/Solution.lean`).** All nine GeneralizedCover specializations use exactly the trace lemma predicted by `e`:

| e | lemma | instances |
|---|---|---|
| 2 | `trace_square` | holdout-x13-y13, holdout-p-y14 |
| −4 | `commonCover_trace_inverse_fourth` | holdout-x1-y04-y1, holdout-p3, holdout-x12-p-y02 (and `holdout` itself) |
| −2 | `trace_inverse_square` | holdout-all-six, holdout-x04-p, holdout-x0-x14-y0, holdout-x0-y03-y12 |

`holdout-bundle` (`X0^3Y0^3 + X0^3Y0Y1^2`, second summand `e = 0`) instead uses the complementary diagonal quotient plus certificate addition; `holdout-q2`, `holdout-moving-branch`, `holdout-ternary-deformation` use other constructions.

OBSERVATION. For a Y-odd monomial `a+b+c` is odd, so `e = b − a − c − 1` is always even; with `b ≤ 6`, `e ∈ {−6,−4,−2,0,2,4}`. Hence `e ≡ 0 (mod 6)` means `e ∈ {0,−6}` (cover fails; served by the diagonal quotient when `e = 0`, and the twelve `e = −6` densities have a fiber pole), and the only cover-served class not yet covered by a library trace lemma is `e = 4`, i.e. `b = 5, a = c = 0`: the three monomials `X_i Y1^5`. A `trace_fourth` for `coverT^4` would follow from the same deck argument (`ζ^4 ≠ 1`) and is the only missing ingredient for those three targets on this variety. Not verified in Lean.

**Relations.** As in §1. `why.md` points to `holdout-x12-p-y02` (solved, `e = −4`) as the comparison target.

---

## 3. klein-hesse (Klein–Hesse · κ = 2)

**Statement.** With `K(z) = z0^5 z1 + z1^5 z2 + z2^5 z0` and `r = X0X1X2`: `F = K(X) + 2 r^2 + K(Y)`, `q = 1`, `A = (U − V)·Y0` where `U = X0^3X1^2`, `V = X1^3X2^2`, `W = X2^3X0^2`. Symmetries: cyclic rotation of the X block (and of the Y block separately), the global involution; each block also keeps the order-7 diagonal action `(X0,X1,X2) ↦ (ζX0, ζ^2X1, ζ^4X2)` (`ζ^7 = 1`): every term of `K` and `r^2` has weight `≡ 0 (mod 7)`. OBSERVATION: `U, V, W` all have weight `0 (mod 7)` under that action (`3+4`, `6+8`, `12+2`), so unlike `klein-coupled-character` there is no order-7 character obstruction for these targets. No same-variety instance. Discussion #85.

**Method.** Constructive, strict (`Q = 0`), explicit rational witness with `d = 15`, numerators of degree 16; only `N0, N1, N2 ≠ 0` (X-components), each 11 terms, all divisible by the passive factor `Y0`, coefficients in `(1/42)·Z`.

**Denominator (PROVEN as witness; closed form COMPUTED).** `solution.md` gives the pencil formula `D_κ = (κ^2 − 3κ + 9)·r^3·f_κ − p^3`, `f_κ = K(X) + κ r^2`, `p = U+V+W`; at `κ = 2` the coefficient is 7. sympy: the witness `D` equals `7 r^3 f_2 − p^3` exactly (ratio 1, not up to sign). Useful identities (COMPUTED):

- `r^3·(K(X) + κ r^2) = U^2W + V^2U + W^2V + κ·UVW` — the Klein block is a cubic form in the monomials `(U,V,W)` (this is the "Hesse relation" of the monomial map: `UVW = r^5`, `r^3K = U^2W + V^2U + W^2V`);
- so `D_κ = (κ^2−3κ+9)(U^2W+V^2U+W^2V+κUVW) − (U+V+W)^3`, a cubic in `(U,V,W)`; `D_κ(1,1,1) = κ^3` (hence Lean's `dFiber = 8`);
- `κ^2 − 3κ + 9 = (κ + 3ω)(κ + 3ω^2)` is the norm of `κ + 3ω` from `Q(ω)`, consistent with the Fourier-channel description;
- `D_2(U,V,W)` is irreducible over `Q` and over `Q(ω)` as a form in `(U,V,W)`; the witness `D` (degree 15 in X) is irreducible over `Q` (sympy `factor`). Over `Q(ω)` in the X-variables not checked.
- `D_0 = 9(U^2W+V^2U+W^2V) − p^3` (the `κ = 0` member) does not divide the klein-rational denominator.

**Mechanism (as described in `solution.md`; the Lean proof only checks the resulting identities).** The monomial map `(X0,X1,X2) ↦ (U,V,W)` exposes the cubic relation above; logarithmic and "complementary Euler" fields in the Fourier directions of the cyclic triple (`u_1 = U + ωV + ω^2W`, `u_2 = U + ω^2V + ωW`) are transported back by the rational monomial Piola formula (the Jacobian density of the monomial map is essential for the divergence). After simplification the auxiliary monomial denominators cancel, leaving the polynomial numerators over `D_κ`. The served packet is the zero-sum directions `{U−V, V−W}` times the three passive linear factors `Y_j`: six target directions at this `κ`; the sum direction `U+V+W` is not served (that is the point of `klein-rational`). COMPUTED: `E1` and `E2` hold exactly, `D` homogeneous of degree 15, all `N_i` of degree 16.

**Lean architecture.** `Witness.lean` (385 lines, 31.7 KB): `import HodgeComputer.Verify`; `maxHeartbeats 2000000`, `maxRecDepth 1000`; namespace `HodgeComputer.Instances.KleinHesse`; explicit `F, A, D, Q := 0, N0..N5, N`; `e1 : (1:Poly) * directional F N = A*D` and `e2 : divergenceNumerator D N = F^1*Q` by `simp only [defs, map_add, pderiv_C, pderiv_mul, pderiv_pow, derivative_i_j …]; apply MvPolynomial.funext; intro x; norm_num <;> ring` (identity checked by evaluation at every point); `denominator_homogeneous : D.IsHomogeneous 15` and `N*_homogeneous … 16` by explicit `isHomogeneous_C/X/.mul/.pow/.add` term trees; coprimality block: `toUni : Poly →+* Polynomial Coeff` sends `X0,X1,X2 ↦ C(X0),C(X1),C(X2)`, `Y0 ↦ X` (the univariate variable), `Y1 ↦ C(X3) + X` (a shear so that `Y0^5Y1` contributes the monic `X^6`), `Y2 ↦ C(X4)`; `fromUni` is the left inverse (`retract : fromUni.comp toUni = id`); `fUni` is `toUni F` with `fUni_lc = C 1`; `specialize := MvPolynomial.eval point` with `point = (1,1,1,1,1)`; `fFiber = X^6 + 2X^5 + 5X^4 + 10X^3 + 10X^2 + 6X + 6`, `dFiber = 8`, Bézout `bezoutF = 0`, `bezoutD = 1/8`, `fiber_coprime` by the `polynomial` tactic; `denominator_legal := relPrime_of_retract toUni fromUni retract (relPrime_of_specialization specialize …)`; `denominator_ne_zero` by evaluating at `(1,…,1)`; `certificate : Verify F A 1 15 D N Q` with `identity := normalized_identity (by simpa using e1) e2`; `solved`. `Solution.lean` (42 lines): `import HodgeWitness` (verifier-generated from `witness.json`), `import Proposal.Witness`; `maxHeartbeats 8000000`, `maxRecDepth 8000`; `witness` proves `HodgeSubmission.F = …KleinHesse.F` etc. by `MvPolynomial.funext; norm_num; ring` and `simpa … using KleinHesse.certificate`; `certificate := ⟨d, D, N, Q, witness⟩`.

**Relations.** No same-variety instance. Family link: `klein-rational` is the `κ = 0` member of the same X-block pencil `K + κ r^2` with the same `K(Y)` block. Open: `klein-coupled-character` (Klein X-block coupled to Y by `10/7 X0^3X2Y1Y2 + 5/7 X0X1^3Y0Y1 + 15/7 X1X2^3Y0Y2`, pure `Y^6` block, target `X0^3Y0^3` with nontrivial `Z/7` character `ζ^3`; the monomial-map/Fourier machinery here is block-local and would have to survive the coupling). `diagonal-permutation` and `symmetric-sextic-residual` are not Klein-type.

**NOT CLAIMED.** Denominator minimality; service of individual monomials `U·Y_j` (that is `klein-rational`); formulas valid only at `κ = 0`.

---

## 4. klein-rational (Klein · rational Piola trace)

**Statement.** `F = K(X) + K(Y)` (additively split, neither block a polynomial in squares), `q = 1`, `A = U·Y0 = X0^3X1^2Y0`. Symmetry as in §3 with `κ = 0` (cyclic, order-7 diagonal on each block, global involution). Discussion #105.

**Method.** Constructive, strict (`Q = 0`), `d = 30`; `N0,N1,N2` of degree 31 with 31 terms each, all `∝ Y0`, coefficients in `(1/28)·Z`; `N3 = N4 = N5 = 0`.

**Denominator (PROVEN as witness; structure COMPUTED).** With `P = U^2W + UV^2 + VW^2`, `Q_K = U^2+V^2+W^2 − UV − UW − VW`:

`D = Q_K^3 + 9·(U+V+W)^3·P` (sympy: equals the witness exactly; 28 terms; `D(1,1,1) = 729`).

In Fourier coordinates `u_0 = U+V+W`, `u_1 = U+ωV+ω^2W`, `u_2 = U+ω^2V+ωW`: `Q_K = u_1 u_2` and `U^3+V^3+W^3−3UVW = u_0u_1u_2` (COMPUTED). sympy `factor(D, extension=√−3)` on the abstract form in `(U,V,W)` gives `D = δ·δ̄` with

`δ = U^3+V^3+W^3 + 3ω^2·(U^2V+V^2W+W^2U) − 6ω·(U^2W+V^2U+W^2V) + 6UVW`,

a cubic in `(U,V,W)` (degree 15 in X) over `Q(ω)`; `δ(1,1,1) = −27ω`, `|δ(1,1,1)|^2 = 729`. So the degree-30 rational denominator is the norm `N_{Q(ω)/Q}(δ)` of a degree-15 denominator defined over `Q(ω)`. This is consistent with `solution.md` ("two conjugate fields are combined by a normalized trace, giving coefficients in Q") but the explanation does not display `δ`; whether the construction's intermediate field actually has denominator `δ` is an OBSERVATION, not something I verified. Over `Q` the witness `D` is irreducible (sympy).

**Mechanism (from `solution.md`).** Diagonalize the cyclic structure with a primitive cube root of unity; build logarithmic fields in the Fourier coordinates; transport through the monomial map by the rational Piola formula; combine the two conjugate fields by normalized trace to land in `Q`. `Y0` is a passive multiplier (no extra derivative). The descended tuple satisfies `N(F) = A·D` and `D·div N = N(D)` exactly — the strict identities; "the denominator norm by itself would not" establish them. Scope: cyclic permutations and the three passive multipliers give nine directions `{U,V,W} ⊗ {Y0,Y1,Y2}` (plus a mirrored Y-block packet). COMPUTED: `E1`, `E2` exact; degrees 30/31.

**Lean architecture.** `Witness.lean` (460 lines, 66 KB), same template as §3 (namespace `HodgeComputer.Instances.KleinRational`, `maxHeartbeats 2000000`, `maxRecDepth 1000`): brute-force `e1`/`e2` by `simp only […]; MvPolynomial.funext; norm_num <;> ring`; homogeneity term trees; identical `toUni` shear (`Y0 ↦ X`, `Y1 ↦ C(X3)+X`), `fUni_lc = C 1`, `point = (1,1,1,1,1)`, `fFiber = X^6+2X^5+5X^4+10X^3+10X^2+6X+4`, `dFiber = 729`, `bezoutF = 0`, `bezoutD = 1/729`; `certificate : Verify F A 1 30 D N Q`. `Solution.lean` (42 lines) identical binding pattern, `maxHeartbeats 8000000`, `maxRecDepth 8000`. Coprimality is trivial here because `D` is `Y`-free: after specializing the X's to 1 the fiber of `D` is the nonzero constant `D(1,1,1)`.

**Relations.** See §3. Compared with klein-hesse: higher-degree denominator buys the individual monomial instead of the zero-sum packet. Open relative: `klein-coupled-character` (same X block, coupled).

**NOT CLAIMED.** Minimality of degree 30; all odd targets.

---

## 5. m3-coupled (M₃ · coupled conic)

**Statement.** `F` is the 22-term sextic
`−116X0^6 + 132X0^4X1^2 − 171X0^4Y0^2 + 14X0^2X1^4 − 52X0^2X1^2Y0^2 + X0^2X2^2Y2^2 + 33X0^2Y0^4 − 30X1^6 + 62X1^4Y0^2 − 34X1^2Y0^4 + X1^2Y0^2Y1^2 + 22X2^6 − 56X2^4Y1^2 + 80X2^4Y2^2 + 32X2^2Y1^4 + 6X2^2Y1^2Y2^2 + 40X2^2Y2^4 − 3Y0^6 + 2Y1^6 + 40Y1^4Y2^2 + 8Y1^2Y2^4 − 2Y2^6`,
`q = 1`, `A = 4X0^3X1X2Y1 + 5X0X1^3X2Y1 − 5X0X1X2Y0^2Y1 = X0X1X2Y1·(4X0^2 + 5X1^2 − 5Y0^2)`. Every exponent in `F` is even, so `F` is a polynomial in the squares (invariant under all six individual sign changes `Z_i ↦ −Z_i`); `A` is odd in `X0,X1,X2,Y1` and even in `Y0,Y2`. `why.md`: an exact Hessian calculation excludes every nontrivial linear direct-sum decomposition over `C`. Discussion #83.

**Method.** Constructive, strict (`Q = 0`), `d = 10` (190 terms, integer coefficients, `D = −2·(irreducible over Q)`), six nonzero numerators of degree 11 with (88, 87, 83, 35, 83, 35) terms, coefficients in `(1/3)·Z`.

**Mechanism (PROVEN as witness; structure COMPUTED).** Quadratic coordinates `U = X0^2 − X1^2`, `W = X2^2 − Y1^2`, `r = Y0`, `s = Y2`, `x = X1^2`, `y = Y1^2`. With `(U,W,r,s)` fixed, `F = a x^2 + h xy + b y^2 + c x + d y + e` with (sympy)

- `a = −70U − 161r^2`, `h = r^2 + s^2`, `b = −14W + 126s^2` (degree 2),
- `c = −216U^2 − 394Ur^2 + Ws^2 − r^4`, `d = Us^2 + 10W^2 + 166Ws^2 + 48s^4` (degree 4),
- `e = −116U^3 − 171U^2r^2 + UWs^2 + 33Ur^4 + 22W^3 + 80W^2s^2 + 40Ws^4 − 3r^6 − 2s^6` (degree 6).

The target is `A = 16·X0X1X2Y1·(k x + m)` with `k = 9/16`, `m = U/4 − 5r^2/16` (COMPUTED), i.e. Jacobian density of the squaring map times an affine function of the conic coordinates. Seek an affine radial fiber field `w = (ρx + p1)∂x + (ρy + p2)∂y` with `w(f) − f·div w = kx + m` (`div w = 2ρ`); comparing coefficients gives the 3×3 system
`[[2a, h, −c],[h, 2b, −d],[c, d, −2e]]·(p1,p2,ρ)ᵀ = (k, 0, m)ᵀ`,
whose determinant is `−2((4ab − h^2)e − bc^2 − ad^2 + hcd)`; sympy: this equals the witness `D` exactly. Piola transport through the quadratic map supplies the density `16X0X1X2Y1` and hence `A`; since the fiber identity is strict (`V(F) − F·div V = A`), one has `K_A(V) = 0`. Finally the Euler correction `V ↦ V − (div V)·ε/6` (the divergence of a degree-1 homogeneous rational field is homogeneous of degree 0, so `div(φε) = 6φ` kills it) gives a field with `N(F) = AD` and `div(N/D) = 0`. COMPUTED: `E1`, `E2` exact; `D` homogeneous of degree 10; the univariate fiber `gcd(fFiber, dFiber) = 1`.

**Lean architecture.** `Witness.lean` (1231 lines, 369 KB): `import HodgeComputer.Verify`; `set_option Elab.async false`, `maxHeartbeats 2000000`, `maxRecDepth 4000`; namespace `HodgeComputer.Instances.M3Coupled`. `e1` is brute force (`simp only […]; MvPolynomial.funext; norm_num <;> ring`). `e2` is NOT brute force: the file defines `radialT` (degree 10), `transportP`, `transportR` and proves the structural decomposition
`N0 = X0·radialT + X1X2Y1·transportP`, `N1 = X1·radialT + X0X2Y1·transportP`, `N2 = X2·radialT + X0X1Y1·transportR`, `N3 = Y0·radialT`, `N4 = Y1·radialT + X0X1X2·transportR`, `N5 = Y2·radialT`,
i.e. `N = radialT·ε + X2Y1·transportP·(X1∂0 + X0∂1) + X0X1·transportR·(Y1∂2 + X2∂4)`. The two bracketed fields annihilate `D` (`tangent_X_D : X1·∂0D + X0·∂1D = 0`, `tangent_Y_D : Y1·∂2D + X2·∂4D = 0` — because `D` depends on `(X0,X1)` only through `U` and on `(X2,Y1)` only through `W`), `euler_D : ε(D) = 10D`, `divergence_N : div N = 10·radialT`, hence `directional_D : N(D) = 10·radialT·D` by a `calc`/`ring` and `e2` closes by `ring` ("avoids expanding products of the degree-ten denominator and its derivatives"). COMPUTED: `N5 = Y2·(N3/Y0)` and `div N = 10·radialT`. Coprimality: `toUni` sends `X0 ↦ X` and the other five variables to coefficient variables (no shear needed; `fUni_lc = C (−116)`), `point = (1,1,1,1,1)`, `fFiber = −116X^6 − 39X^4 − 4X^2 + 168`, `dFiber = 3572352X^8 − 2288364X^6 − 735534X^4 − 1029524X^2 − 1071854`, explicit Bézout cofactors of degrees 6 and 4 with denominators such as `57446509982104879067201427` (≈ 5.7·10^25), `fiber_coprime` by `polynomial`; then `relPrime_of_retract`/`relPrime_of_specialization` as in §3. `Solution.lean` (42 lines) binds the generated `HodgeWitness` as in §3 (`maxHeartbeats 8000000`, `maxRecDepth 8000`).

**Relations.** No same-variety instance. The target is that of the solved split instance `passive-quadratic-d4-strict-rational` (`4X0^3X1X2Y1 + 5X0X1^3X2Y1`, `F` split) plus `−5X0X1X2Y0^2Y1`. Technique-wise this is the explicit Cramer-determinant analogue of "the degeneration locus of a moving conic organizes the denominator" (Chapter 2 §2.6); it does not transfer directly to any of the 16 open instances, none of which is a polynomial in squares with such a conic fibration (the open `quadric-contact-residual` has odd-exponent couplings).

**NOT CLAIMED.** Denominator minimality; a general even-sextic theorem.

---

## 6. mixed-biweight (Mixed biweights · asymmetric)

**Statement.** `F = X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6 + (1/2)X0^3X1Y0Y1 + (1/2)X1^3X2Y1Y2 + (1/3)X0X2Y0^3Y2`, `q = 1`, `A = X0^3Y0^3`. The chain terms have `(X,Y)`-bidegree `(4,2)`, the new term `(2,4)`. Symmetries: global involution and all three paired sign changes `(X_i,Y_i) ↦ (−X_i,−Y_i)` survive (each mixed term has even degree in each pair). `why.md` shows why the holdout cover fails: on `T^6 = Y0^6+Y2^6` the new term becomes `(1/3)α^3βX0X2T^4` with no `y^4` partner, so `y ↔ T` is not a symmetry and `F` does not descend to `(S,P)`; the complementary diagonal action `(5,1,5)` on `Y` multiplies it by `ζ^2`. Discussion #56.

**Method.** Nonconstructive. Paired-sign reduction to a plane cubic, Hessian covariants to a Weierstrass Jacobian, a polynomial-family primitive on the Jacobian, first-variation (dual-number) source transport, Piola pullback with moving-map and normal-chart corrections, then the shared clearing/homogenization. Adapts `hodgepodgeconjecture`'s `cycle-three` proof (PR #165) with a different (transitive) cubic family. The cubic machinery was afterwards extracted into `HodgeComputer.Constructions.TransitiveCubic`.

**Mechanism (PROVEN).**
1. Chart: `R = X2^2`, `v = Y0/X0`, `c1 = X0Y1/(X1Y0)`, `c2 = X0Y2/(X2Y0)`, `x = X0^2/R`, `y = X1^2/R` (all in `Frac(C[Z]/(F))`, chart denominator `X0X1X2Y0`). Then `F = R^3·U(x,y,1)` with the transitive cubic
   `U = a x^3 + b y^3 + c z^3 + d x^2y + e y^2z + f x^2z`, `a = 1+v^6`, `b = 1+c1^6v^6`, `c = 1+c2^6v^6`, `d = c1v^2/2`, `e = c1c2v^2/2`, `f = c2v^4/3`,
   and `A = v^3R^3x^3`. (cycle-three instead has the directed cycle with `f z^2 x` and `f = c2v^2/2`; here the three edges `0→1, 1→2, 0→2` are "transitive".)
2. Covariants: integral `H, Θ, J` (3, 9, 27 times classical), invariants `s4 = (−9bcdf + 3cd^2e − e^2f^2)/3`, `s6 = (729a^2b^2c^2 + 108a^2ce^3 + 108b^2cf^3 + 108bc^2d^3 − 108bcdef^2 + 72cd^2e^2f + 8e^3f^3)/27`; full polynomial syzygy `27J^2 = 108Θ^3 + 36 s4·Θ H^4 − s6·H^6 + U·(…)` (`CubicSyzygy.full_syzygy`), so `ξ = Θ/H^2`, `η = −J/(2H^3)` satisfy `η^2 = ξ^3 + s4 ξ − s6/4` on `U = 0`; invariant differential `D_yξ/(2η) = 3/U_x` (`CubicFieldDifferential.Physical.invariant_differential`). The 190 syzygy coefficients are proved by an exact encoding (`Kronecker` weights `(1,19,361)`, signed split into natural polynomials, evaluation at `2^42`, `PolynomialDigitCheck`), kernel-checked in `CubicCoefficientsA`–`H` and `CubicIdentityData.lean` (8705 lines, 1.57 MB).
3. Jacobian model (`JacobianModel.lean`): `T = Rv^6`, `X = R^2ξ`, `Y = R^3η`; `V = 5832Y^2 − 5832X^3 − 27·a8·X − b256`, with `a8 = 216a4`, `b256 = 5832b6` (names kept from the cyclic proof), explicitly
   `a8 = 27c1^3c2r^2t(r+c2^6t) − 108c1c2rt(r+c1^6t)(r+c2^6t) − 2c1^2c2^4r^2t^2`,
   `b256 = −39366(αr+t)^2(r+c1^6t)^2(r+c2^6t)^2 − 729c1^3c2^3r^2t(αr+t)^2(r+c2^6t) − 216c2^3rt^2(r+c1^6t)^2(r+c2^6t) − 729c1^3r^2t(r+c1^6t)(r+c2^6t)^2 + 162c1^2c2^3r^2t^2(r+c1^6t)(r+c2^6t) − 81c1^4c2^3r^3t^2(r+c2^6t) − 2c1^3c2^6r^3t^3`
   (weights 4 and 6 in `(r,t)`; `α` is the auxiliary parameter in `a = α + v^6`, set to 1).
4. Primitive on the Jacobian: polynomial-curve family of coordinate degrees `(2,4,6)` in `(R,X,Y)` with parameter `T`, 13 unknown coefficients, 13 incidence equations; nonsingular seed over `F_5` at `(c1,c2,R(0),X(0)) = (1,1,1,3)`: `R = 1+T^2`, `X = 3+4T^2+3T^4`, `Y = 2+3T^6` (`JacobianSeed`, 13×13 inverse checked by `decide`; `JacobianSeedBinding` binds it to the actual derivatives). `PolynomialCurves.ArithmeticLifting` lifts to a nonsingular solution over a finite separable extension of `C(c1,c2,R0,X0)`; `full_polynomial_combination_one` builds a polynomial normal (`conicNormal_pairing = 1`); `PolynomialFamilyBridge` integrates; `FiniteTrace.residue_solution_descends_finite` descends. Result `JacobianResidue.target_residue_solution`: every polynomial source on the Jacobian hypersurface has a rational residue primitive (stronger than needed). The source used is the α-source `P_J = −V_α/(…)`: `sourceNumerator = −(R+T)(R+c2^6T)(54(R+c1^6T)^2(R+c2^6T) + c1^3c2^3R^2T)`, `source = sourceNumerator/576`, with `∂_α V|_{α=1} = −839808·R·source`.
5. First variation: in `K[ε]/(ε^2)` (`TrivSqZeroExt`), `a_ε = a + ε`, `x_ε = x − εx^3/g` (`g = U_x`), other data fixed; yields `ġ/g = (ȧ_Wξ + ḃ_W)/(2η^2) − (g/3)D_y(ξ̇/(2η))`, `ġ = 3x^2 − x^3U_xx/g` (`CubicJetCalculus`, `PhysicalCubicJets`).
6. Transport: tangent derivations `δ_i = ∂_{i+1} − (F_{i+1}/F_0)∂_0`; passive coordinates `(T, c1, c2, R, X)`; matrix `M_{ji} = δ_i(passive_j)` invertible with density multiplier `F_0·det M/(2Y) = −144Rv^3` (`PhysicalJacobianCoordinates.density_multiplier`); `CoordinateIndependence` (characteristic-zero degree descent) gives injectivity; the monic irreducible Jacobian equation embeds `TargetField ↪ FunctionField` (`fieldMap`), intertwining derivations (`tangent_compatible`). Piola pullback (`ResidueTransport.piola_pullback`) of the Jacobian primitive has residue divergence `−v^3(ȧ_Wξ+ḃ_W)/(2η^2)`; add `(ξ̇/(2η))/3·W` with `W = λ(hp1, 0, −h1p3, hp4, 0)`, `λ = p3^3/(4p0^5p1^2p2^2)`, `h = p0F0 + p3F3`, `h1 = p1F1 + p4F4` (zero residue divergence, acts as `v^3gD_y`; `PhysicalCubicVolume`); subtract the normal-chart correction `C = −p0^3p3^4/h` in the `p3` slot, with `A_0/F_0 − AF_{00}/F_0^2 = v^3ġ/g + div(0,0,C,0,0)` (`physical_source_transport`). `PhysicalPrimitive.physical_residue_solution` assembles `w = pulled + correction − correctionVector`.
7. Closure: `PhysicalClearing.certificate_of_residue_solution` = shared `residue_solution_clearing` (needs `sextic_irreducible`, `residueMap_kernel`, `normal_ne_zero`, commuting partials, `residue_common_representatives`) + `hasCertificate_of_normalized`; `Solution.lean` rewrites `HodgeSubmission.F = CyclePhysical.sextic`, `HodgeSubmission.A = CyclePhysical.source` by `unfold; ring` and `exact CyclePhysicalPrimitive.certificate`. The denominator again arises existentially as `∂0F·d` for a common polynomial denominator `d` of trace-descended representatives, then homogenized.

**Coprimality over ℂ.** `PhysicalField.sextic_irreducible`: `F` is monic of degree 6 in `X0` (`finSuccEquiv`), and the specialization `holdoutSpecialization Complex.I` (library lemma family: `X1 = X2 = Y1 = 0`, `Y2 = i` with `i^6 = −1`, `Y0 = t + 1`) maps the equation to `X^6 + C((X+1)^6 − 1)`, irreducible by the Eisenstein-type `irreducible_power_add_constant`; `Monic.irreducible_of_irreducible_map` lifts it. `normal_legal`: restriction to the `X0`-axis sends `F ↦ X0^6`, `∂0F ↦ 6X0^5`. Then `legal_of_residue_ne_zero` inside `residue_solution_clearing`. Nonvanishing of `H, J` (hence `η`, `D_yξ`) in the function field is shown by the explicit point `CubicPoint` (`(a,b,c,d,e,f,x,y) = (1 − i/8, 1,1,0,0,0, 2i, −2i)`, `H = 324 − 81i/2`).

**Lean architecture (sizes, options).** 59 helper modules + `Solution.lean` (20 lines). Largest: `CubicIdentityData.lean` 8705 lines/1.57 MB, `CubicReflectionData.lean` 302 lines/160 KB, `CubicDifferential.lean` 1301 lines/110 KB, `CubicCovariants.lean` 105 lines/68 KB. Heartbeat settings: `CubicCoefficientsA–H`, `CubicIdentityData`, `CubicIdentityProducts`: `maxHeartbeats 32000000`, `maxRecDepth 131072`; `CubicSyzygy` 8000000/131072; `CubicDifferential`, `CubicIdentityEvaluation` 4000000; `JacobianIncidence`, `JacobianSeedBinding` 8000000 (`maxRecDepth 100000`); `JacobianSeed.inverse_product` by `decide` under 8000000/100000; `JacobianNormal` 4000000/10000; `PhysicalJacobianCoordinates`, `PhysicalJacobianMap`, `PhysicalCubicVolume` 3000000; most other `Physical*`/`Cubic*` 2000000; `JacobianFamily`, `JacobianLifting`, `JacobianResidue`, `JacobianHypersurface`, `JacobianSource`, `PhysicalField` 1000000; `JacobianDimension`, `JacobianEmbedding`, `JacobianGeometry` 200000; `JacobianParameters` 100000. Several `Physical*` modules set `backward.isDefEq.respectTransparency false`. Identities are discharged by `simp`/`ring`/`linear_combination`/`field_simp` over the function field; the syzygy data by `rfl`-generated products and kernel reduction; the seed Jacobian by `decide` over `ZMod 5`. Shared imports: `HoldoutCover.{PolynomialCalculus, HypersurfaceDomain, AlgebraicFields, PolynomialNormal, PolynomialFieldMap, EvaluationMatrix, FiniteDifferentials, ResidueDerivatives, PolynomialFamilyBridge, FiniteTrace, ResidueClearing, ResiduePiolaPullback, HomogeneousCertificate}`, `PolynomialCurves.ArithmeticLifting`, `ResidueTransport.{MovingMap, Density}`. `witnessChecked: false`; acceptance hashes all 60 Lean files.

**Relations.** Model: `cycle-three` (solved; cyclic cubic). Open on the cycle-three variety: `cycle-q2` (`q = 2`, `A = X0^2X1^2X2^3Y0^2Y1^2Y2`, nontrivial character `ζ_6^2` under the diagonal `(2,1,0,3,2,1)` action; needs an `F^3` identity, i.e. third-order machinery as in `HigherResidues`/`split-ternary-q2`, not the second-order Jacobian identity). Open deformation of the same chain: `off-paired` (`X0^2X1X2Y0Y2/3`; breaks the individual paired signs in pairs 0,1 — only their product and the third survive — so the `u_i = X_i^2` cubic reduction is not available unchanged; its `why.md` explicitly compares to cycle-three). Other open `X0^3Y0^3` targets where neither the cover nor an obvious paired-sign cubic applies: `dense-rational`, `four-plus-two-residual`, `one-reflection-skew`, `rank-one-biquartic` (`f(X)+g(Y)+tL(X)^2M(Y)`, rank-one `(4,2)` coupling, no paired signs; its `why.md` suggests a binary-cubic fiber `as^3 + bs^2u + du^3`). Solved neighbours using the same lifting/trace toolkit: `mixed-polar`, `holdout`, `cyclic-triple-fano`, `reflection-coupled-fano`, `crossed-ternary-asymmetric`.

**NOT CLAIMED.** Explicit witness; denominator degree; a theorem for every paired-sign sextic (the library README is explicit about this).

---

## 7. Cross-cutting observations

1. **Two coprimality idioms.** Constructive witnesses prove `IsRelPrime F D` by `relPrime_of_retract ∘ relPrime_of_specialization` plus an explicit univariate Bézout identity on one fiber (`fiber_coprime`, tactic `polynomial`); the only requirement is a variable (after an optional shear such as `Y1 ↦ X3 + t`) in which `F` has unit leading coefficient. When `D` omits a variable block (klein-*: `D` is X-only) the fiber of `D` can be made a nonzero constant and the Bézout data is `(0, 1/const)`. Nonconstructive proofs instead prove `Irreducible F` (monic in `X0` + Eisenstein specialization) and use `legal_of_residue_ne_zero` with `D = ∂0F·d`.
2. **Identity discharge.** Constructive: `MvPolynomial.funext` + `norm_num <;> ring` after `simp only` with per-variable derivative tables; m3-coupled shows the scalable pattern for large `D`: decompose `N` into `radialT·ε` plus fields tangent to `D`, so `E2` follows from Euler + tangency without expanding `D·∂D` products. Nonconstructive: all identities live in function fields and are closed by `field_simp`/`ring`/`linear_combination`; large polynomial identities are moved into kernel-checked data modules.
3. **Where existential denominators come from.** In both nonconstructive families: trace from a finite extension of the family's coefficient field, (cover case) a second deck trace to the physical field, choice of polynomial representatives in `Frac(C[Z]/(F))`, `D = ∂0F·d`, graded reduction and a power `h^p`. No norm is ever specialized.
4. **Holdout cover characters.** `e = b − a − c − 1` is always even for Y-odd monomials; library trace lemmas exist for `T^2`, `T^{−2}`, `T^{−4}`; `e = 4` (`X_iY1^5`) has no lemma yet (OBSERVATION).
5. **Klein pencil.** `r^3(K + κr^2) = U^2W + V^2U + W^2V + κUVW`; the klein-hesse denominator is a cubic form in `(U,V,W)` with norm coefficient `κ^2−3κ+9 = N(κ+3ω)`; the klein-rational denominator is `N_{Q(ω)/Q}(δ)` for an explicit cubic `δ` in `(U,V,W)` (COMPUTED). Whether a degree-15 certificate over `Q(ω)` with denominator `δ` exists for the individual monomial is a CONJECTURE suggested by this factorization, not checked.
6. **Failed attempts recorded by the authors (none of mine).** The holdout cover and the `(5,1,5)` diagonal quotient both fail on `mixed-biweight`, `off-paired`, `one-reflection-skew` (`why.md`s); direct paired-sign quotient arguments leave the two holdout targets of §1–2 unsupported (48/58-dimensional statement); the Klein–Hesse difference construction does not serve individual monomials.
