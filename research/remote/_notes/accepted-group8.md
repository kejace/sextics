# Accepted solves, group 8: contact quadrics, broken contact, coupled Fano, split ternary q=2

Catalogue of six accepted instances in `hodge-computer/sextics` (repo at `/home/kejace/devel/sextics`, read 2026-11-01):
`ternary-tangent-residual-h2-d2`, `ternary-collision-cb-d2`, `ternary-collision-cbp-d2`,
`quadric-contact-broken`, `reflection-coupled-fano`, `split-ternary-q2`.

Sources read per instance: `instance.rec`, `why.md`, `solution.md`, `proof.json`, `provenance.json`,
`checks.json`, `acceptance.json`, `witness.json` (where present), `Solution.lean`, all helper `.lean` files
(declaration lists, imports, `set_option`s, and the terminal modules in full). Independent sympy checks
were run locally (`/tmp/hodge/scratch/verify_group8.py`, sympy 1.12); results are marked **[sympy]**.
Nothing in the repository was modified.

Variable order is always `X0, X1, X2, Y0, Y1, Y2` (Lean `Fin 6` indices `0..5`).

## 0. Certificate conventions used below

For `(F, q, A)`, a certificate is a nonzero homogeneous `D` with `IsRelPrime F D` over `ℂ`, and
`N_0..N_5` homogeneous of degree `deg D + 6q − 5`, with the cleared identity

    E := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^(q+1)·Q .

Normalized form (sufficient, `HodgeComputer.normalized_identity`):
`E1: q·Σ N_i ∂_iF = A·D` and `E2: D·div N − N(D) = F^q·Q`. "Strict" means `Q = 0`.

Lean contract (`formal/lean/HodgeComputer/Verify.lean`):

    structure Verify (F A : Poly) (q d : ℕ) (D : Poly) (N : Fin 6 → Poly) (Q : Poly) : Prop where
      denominator_ne_zero : D ≠ 0
      denominator_homogeneous : D.IsHomogeneous d
      numerator_homogeneous : ∀ i, (N i).IsHomogeneous (d + 6*q - 5)
      denominator_legal : IsRelPrime F D
      identity : defect F A D q N = F^(q+1) * Q
    def HasCertificate (F A : Poly) (q : ℕ) : Prop := ∃ d D N Q, Verify F A q d D N Q

Two shared coprimality lemmas appear repeatedly:

    relPrime_of_specialization (φ : R →+* S) (f g : Polynomial R)
        (hf : IsUnit f.leadingCoeff) (h : IsRelPrime (f.map φ) (g.map φ)) : IsRelPrime f g
    relPrime_of_retract (f : R →+* S) (g : S →+* R) (hg : g.comp f = RingHom.id R)
        (h : IsRelPrime (f a) (f b)) : IsRelPrime a b

Pattern: pick one variable `t` in which `F` is monic, map `ℂ[X,Y] → (ℂ[other five])[t]` (`toUni`),
specialize the other five at a point so that the fibre polynomials are visibly coprime, lift. This
**does not** require `V(F) ∩ V(D) = ∅`.

Summary table:

| instance | q | A | contributor / PR / accepted | kind | D |
|---|---|---|---|---|---|
| ternary-tangent-residual-h2-d2 | 1 | `X0·Y0·Q2(u,w,x)` | benediamond, PR #6, 2026-10-15 | constructive, strict | `X0²+X1²−2Y0²`, d=2 |
| ternary-collision-cb-d2 | 1 | `X0·Y0·Q2(u,w,x)` | benediamond, PR #9, 2026-10-15 | constructive, strict | `X0²−X1²−Y0²`, d=2 |
| ternary-collision-cbp-d2 | 1 | `X1·Y0·Q2(u,w,x)` | benediamond, PR #8, 2026-10-15 | constructive, strict | `X1²−X0²−Y0²`, d=2 |
| quadric-contact-broken | 1 | `X1²·X0Y0·(X0²−Y0²)` | HeadOfNaps, PR #122, 2026-10-19 | existential but explicit up to a clearing exponent; strict (Q=0) | `a²c²·Q`, d = 2n+58, n unspecified |
| reflection-coupled-fano | 1 | `X0³·Y2³` | fqlx, PR #150, 2026-10-20 | nonconstructive (quotient cover + trace), Q allowed ≠ 0 | existential, degree unspecified |
| split-ternary-q2 | 2 | `X0⁴X1³Y0³Y1²` | fqlx, PR #159, 2026-10-20 | nonconstructive (elliptic companions + higher residue), strict (Q=0) | existential, degree unspecified |

All six were verified with Lean 4.33.0, kernel replay, axioms exactly `propext, Quot.sound, Classical.choice`.
`witnessChecked` is `true` only for the three ternary instances (the others export no `witness.json`).

---

## 1. `ternary-tangent-residual-h2-d2` ("Ternary tangent · H₂")

**Record.** `q = 1`. With `(u,w,x) = (X0², Y0², X1²)`:

    F = p(u,w,x) + X2⁶ + Y1⁶ + Y2⁶,
    p = u³ + u²w + u²x + uw² + uwx + ux² − 14w³ + (17/2)w²x + wx² − (3/2)x³,
    A = X0⁵Y0 − (1/25)X0X1⁴Y0 − (293/25)X0X1²Y0³ + (699/25)X0Y0⁵
      = X0·Y0·( u² − (1/25)x² − (293/25)wx + (699/25)w² ).

Character: `A` is odd under `X0 ↦ −X0` and `Y0 ↦ −Y0`, even under `X1 ↦ −X1` (and under the sign/μ6 actions
on the passive `X2, Y1, Y2`). The sextic is invariant under every individual coordinate sign change
(`X0, X1, Y0` enter squared; `X2, Y1, Y2` as sixth powers), so its diagonal symmetry group contains
`μ2³ × μ6³`. No `checks.json` exists for the three ternary instances. Discussion #68.

**Method.** Contact geometry of a plane cubic in square coordinates → quadric denominator with a
logarithmic cubic numerator ("tangent-and-residual" configuration). Constructive, strict, `witness.json`:

    d = 2,  D = X0² + X1² − 2·Y0²            (pullback of the line ℓ = u − 2w + x)
    N0 = Y0·(41X0² + 33X1² − 78Y0²)/150
    N1 = −2·X0·X1·Y0/75
    N3 = −X0·(24X0² + 27X1² − 49Y0²)/75
    N2 = N4 = N5 = 0,   Q = 0 .

**Mechanism (proven, record + [sympy]).**
- On `ℓ`, parametrized by `(u,w,x) = (s+2r, s+r, s)`, the cubic restricts to `−(5/2)·s·r²`:
  the line is tangent at `(1:1:1)` and has residual intersection `(2:1:0)`, which lies on the coordinate
  line `x = 0` (the branch coordinate of the *even* active variable `X1` for this character).
  [sympy: `p|_{u=2w−x} = −(5/2)·x·(x−w)²`.]
- The field is logarithmic along `D`: `div N = (137/75)·X0·Y0` and `N(D) = (137/75)·X0·Y0·D`, so
  `E2 = D·div N − N(D) = 0` exactly; and `N(F) = A·D` (E1). Both verified **[sympy]**; `gcd(F,D) = 1` over ℚ.
- `why.md` stresses that tangency alone is not the criterion; the *placement of the residual point on the
  branch line* is what makes the pulled-back line serve the `X0·Y0`-character. This is stated as an
  observation on this fibre, not a theorem for all tangent lines.

**Lean architecture.**
- `Witness.lean` (355 lines, `maxHeartbeats 2000000`, `maxRecDepth 1000`), namespace
  `HodgeComputer.Instances.TernaryTangentResidualH2D2`. Defines `F, A, D, Q, N0..N5` as explicit
  `MvPolynomial (Fin 6) ℂ` sums, `N := fun i => if i.val = 0 then N0 else …`.
  36 private lemmas `derivative_i_j : pderiv i (X j) = …` by `norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]`.
- `e1 : (1:Poly) * directional F N = A*D` and `e2 : divergenceNumerator D N = F^1*Q`:
  `simp only [directional, Fin.sum_univ_succ, …]`, `change` to the explicit 6-term sum,
  `simp only [F, A, D, N0.., map_add, pderiv_C, pderiv_mul, pderiv_pow, derivative_*]`, then
  `apply MvPolynomial.funext; intro x; norm_num <;> ring`.
- Homogeneity: explicit chains `((isHomogeneous_C …).mul ((isHomogeneous_X ℂ i).pow k)).add …`;
  zero components via `isHomogeneous_zero`.
- `denominator_ne_zero`: evaluate at `(1,1,1,1,1,1)`, `norm_num`.
- Coprimality over ℂ: `toUni : Poly →+* Polynomial Coeff` sends `X2 ↦ Polynomial.X` and
  `X0,X1,X3,X4,X5 ↦ C (X0..X4)` (`Coeff = MvPolynomial (Fin 5) ℂ`); `fromUni` is the inverse on the image;
  `retract : fromUni.comp toUni = id` by `MvPolynomial.ringHom_ext` + `fin_cases`.
  `fUni = X^6 + C(rest)` with `fUni_degree = 6` (`compute_degree!`) and `fUni_lc = 1`.
  `point = (1,2,3,4,5)`, `specialize = MvPolynomial.eval point`, fibre polynomials
  `fFiber = X^6 + 12464`, `dFiber = C(−13)`, Bézout `bezoutF = 0, bezoutD = C(−1/13)` closed by the
  `polynomial` tactic; `fiber_coprime : IsCoprime fFiber dFiber`, then
  `denominator_legal := relPrime_of_retract toUni fromUni retract (relPrime_of_specialization specialize … fiber_coprime.isRelPrime)`.
- `certificate : Verify F A 1 2 D N Q` with `identity := normalized_identity (by simpa using e1) e2`;
  `solved : HasCertificate F A 1`.
- `Solution.lean` (`maxHeartbeats 8000000`, `maxRecDepth 8000`): `theorem witness : Verify HodgeSubmission.F … HodgeWitness.d HodgeWitness.D HodgeWitness.N HodgeWitness.Q`
  by proving `HodgeSubmission.F = Instances….F` etc. (`MvPolynomial.funext; norm_num [...] <;> ring`,
  `funext i; fin_cases i` for `N`) and `simpa … using Instances….certificate`; `certificate := ⟨d,D,N,Q,witness⟩`.
  No `HodgeComputer.Constructions.*` imports.

**Related instances.** Same family (all solved by benediamond with the identical template):
`collision-h3` (same line `u − 2w + x`, cubic `p = (u−2w)³ + ℓ·Q2`, triple contact at `(2:1:0)`,
`D = X0² + X1² − 2Y0²`, `N0 = −2Y0(X0² − 2Y0²)`, `N3 = X0X1²`, `div N = −4X0Y0`), and the two
collision characters below. No open catalogue instance is of this "plane cubic in squares plus passive
sixth powers" shape; `quadratic-cubic-pullback` (open) is a different construction (`F = Σ Q_i³`).

---

## 2. `ternary-collision-cb-d2` ("Ternary collision · first character")

**Record.** `q = 1`, `(u,w,x) = (X0², Y0², X1²)`:

    p = −3u³ + 4u²w + u²x − 2uw² + uwx + ux² + w³ + w²x + wx² − x³,   F = p + X2⁶ + Y1⁶ + Y2⁶
    A = −X0⁵Y0 + (7/5)X0X1⁴Y0 + (21/5)X0X1²Y0³ + 4X0Y0⁵ = X0·Y0·( −u² + (7/5)x² + (21/5)wx + 4w² ).

Character `X0·Y0·ℚ[u,w,x]_2` (six-dimensional raw target space; two of those directions,
`X0Y0·p_u` and `X0Y0·p_w`, come from polynomial fields; the quotient is 4-dimensional). Discussion #71, PR #9.

**Witness (constructive, strict).**

    D = X0² − X1² − Y0²     (line ℓ = u − w − x)
    N0 = Y0·(5X0² − X1² − 5Y0²)/10 = Y0·(5D + 4X1²)/10
    N3 = X0·(5X0² − 3X1² − 5Y0²)/5 = X0·(5D + 2X1²)/5
    others 0, Q = 0;   div N = −X0·Y0,  N(D) = −X0·Y0·D.

**Mechanism.** Inflectional collision: `p|_{u = w+x} = −2·x³` **[sympy]**, i.e. the line meets the cubic
three times at `(u:w:x) = (1:1:0)`, on the coordinate line `x = 0` (again the branch coordinate of the even
active variable). Tangency point and residual point have merged. `E1, E2, E` all vanish **[sympy]**;
`gcd(F,D)=1` over ℚ.

**Lean.** Same template as §1 (Witness.lean 355 lines, same options and lemma names). Specialization point
`(1,1,1,1,1)` → `fFiber = X^6 + 6`, `dFiber = C(−1)`, Bézout `(0, −1)`. `F` monic in `X2`.

---

## 3. `ternary-collision-cbp-d2` ("Ternary collision · second character")

**Record.** `q = 1`:

    p = u³ + u²w + u²x + uw² + uwx + ux² − 3w³ + 10w²x − 8wx² + x³,   F = p + X2⁶ + Y1⁶ + Y2⁶
    A = −X0⁴X1Y0 − 3X0²X1Y0³ + 20X1³Y0³ − (265/17)X1Y0⁵ = X1·Y0·( −u² − 3uw + 20wx − (265/17)w² ).

Character `X1·Y0·ℚ[u,w,x]_2` (X0-exponent even; odd under `X1 ↦ −X1`, `Y0 ↦ −Y0`). Discussion #70, PR #8.

**Witness (constructive, strict).**

    D = X1² − X0² − Y0²    (line ℓ = x − u − w; witness.json writes −X0² + X1² − Y0²)
    N1 = Y0·(21X0² − 40X1² + 40Y0²)/34 = Y0·(−19X0² − 40D)/34
    N3 = −X1·(4X0² + 15X1² − 15Y0²)/34 = −X1·(19X0² + 15D)/34
    others 0, Q = 0;   div N = −(25/17)X1·Y0,  N(D) = −(25/17)X1·Y0·D.

**Mechanism.** `p|_{x = u+w} = 4·u³` **[sympy]**: triple contact at `(0:1:1)`, on `u = 0` (the branch
coordinate of the even variable `X0`). The active numerator slots move from `{0,3}` to `{1,3}` with the
character. `E = 0` and `gcd(F,D) = 1` **[sympy]**.

**Lean.** Same template; point `(1,1,1,1,1)` → `fFiber = X^6 + 8`, `dFiber = C(−1)`.

**Observed pattern across the four accepted ternary-contact solves (observation, not a theorem):**
in each case `ℓ` has coefficients in `{±1, ±2}`, `D` is the pullback quadric, the numerator is supported in
the two odd slots of the character plus (for `h2`) one extra slot, every `N_i` is (coordinate)·(quadric), and
`N` is logarithmic along `D` with `D·div N = N(D) = c·(odd monomial)·D`, `c ∈ {−1, −25/17, 137/75, −4}`.
The contact point (collision) or the residual point always lies on the coordinate line of the *third*,
even, square coordinate. Whether this placement is necessary is not established anywhere in the records.

---

## 4. `quadric-contact-broken` ("Quadric contact · broken reflections")

**Record.** `q = 1`.

    F = F* + (5/7)·X0·X1·Y0·Y1·(X2² + 2Y2²),
    F* = X0⁶+X1⁶+X2⁶+Y0⁶+Y1⁶+Y2⁶ + 5X0⁴Y0² + 5X0²Y0⁴ + 5X0²X1²Y1² + 5X0²X1X2Y1Y2 + 5X0²X2²Y2²
         + 5X1²Y0²Y1² + 5X1X2Y0²Y1Y2 + 5X2²Y0²Y2²        (= two-reflection-cubic-surface-d2)
    A = X0³X1²Y0 − X0X1²Y0³ = X1²·X0·Y0·(X0² − Y0²).

[sympy: `F − F* = 5·X0X1Y0Y1(X2²+2Y2²)/7`.] `checks.json`: smooth by reduction mod 5; **no coordinate
reflections**; projective diagonal group of order 24, no junior element fixing the target; Hessian
symmetrizer space scalar (no hidden linear block split); target nonzero in the reduced Jacobian ring.
Provenance: `sameVarietyInstances = [quadric-contact-residual]`; discussion #119; PR #122 (HeadOfNaps,
submitted 2026-10-19 02:02Z, accepted 07:37Z; 56 Lean modules, 5606 lines; `witnessChecked: false`).

**Proven obstruction (why.md, exact linear algebra).** With `D = X0² + Y0²` (the undeformed denominator) or
any nonzero multiple, no numerator tuple works: a legal certificate forces `A·D ∈ (∂0F,…,∂5F)`; the system
in the 336 cubic coefficients has rank 336 and augmented rank 337 modulo 101, and full column rank makes
this a characteristic-zero inconsistency certificate. The old field
`v = X1²/(4(X0²+Y0²))·(−Y0,0,0,X0,0,0)` has `dF*(v) = A*`, `div v = 0`, and on `F` the error is
`dF(v) − A = (5/28)·X1³Y1(X0²−Y0²)(X2²+2Y2²)/(X0²+Y0²)` **[sympy]**. Nothing is claimed about `(X0²+Y0²)^k`,
`k ≥ 2`, or other supports.

**Method name.** "Finite-algebra contact certificate" / cubic-surface reduction with a line over a
degree-27 finite coefficient algebra, regular-representation (norm-type) denominator, horizontal clearing.
Prepared with OpenAI Codex assistance (attribution in solution.md).

**Constructive?** Semi-explicit. The cubic-level denominator `D₀ = det(ρ(L))` is a concrete 27×27
determinant, but the final `D = a²c²·Q` and `N` come from an existential clearing lemma
(`exists_cleared_cubic`) whose exponent `n` is not computed; final degree `d = 2(n+27) + 4 = 2n + 58`,
numerators of degree `2n + 59`. The quotient is `Q = 0` (strict). No `witness.json`.

**Mechanism (proven in Lean; the key identities also [sympy]-checked).**
1. Coordinates `(x,a,c,y,b,d) = (X0,X1,X2,Y0,Y1,Y2)`, ratios `α = b/a`, `β = d/c`, quadratic coordinates
   `S = 2(x²+y²)`, `T = 4xy`, `B = a²`, `C = c²`. Note `∂(S,T)/∂(x,y) = 16(x² − y²) = 16δ`; the factor `δ`
   in `A` is exactly this Jacobian density, which is why the construction serves this character and not
   `X0³Y0³` (see `quadric-contact-residual` below).
2. Cubic: `G = 2S³ + ST² + S(pB² + qBC + sC²) + kTBC + uB³ + vC³` with
   `(p,q,s,k,u,v) = (40α², 40αβ, 40β², (20/7)α(1+2β²), 16(1+α⁶), 16(1+β⁶))`; **`G = 16F`** exactly
   **[sympy]**, and `BT = 4a²xy`, `(δ/4)·BT = δa²xy = A` **[sympy]**.
3. Work over `O = ℚ[α,β]` localized at all polynomials nonzero at `(α,β) = (1,2)` (`RegularSpecialization`,
   `ParameterSpecialization.point`). A line `S = ξB + ηC, T = ζB + ωC` on `G = 0`, with
   `η = rξ`, `ω = rζ + (v − ur³)/(kr)`; eliminating `ζ` leaves a quadratic and a cubic in `ξ`, whose
   resultant is `k¹²r⁹·(degree-27 polynomial in r)` with leading coefficient `u¹⁰v` (a unit in `O`).
   Normalize to monic `P` (`ParameterEliminant.P`, `P_natDegree = 27`); `E = O[r]/(P) = AdjoinRoot P`, free
   with basis `1, r, …, r²⁶`. Subresultant formula gives `ξ`; the chart gives `η, ζ, ω`; the four
   restriction coefficients `rel0..rel3` of `G` vanish (`GenericCubicLines.rel*_zero`). Units
   (`r`, the quadratic leading coefficient, the linear subresultant coefficient) are proved units via the
   specialization at `(1,2)` with explicit rational Bézout identities and lifted through nonzero specialized
   determinants (`isUnit_of_special_isUnit`). **Irreducibility or squarefreeness of `P` is never used.**
4. Strict cubic certificate (`ConcreteLine.lean`): `L = S − ξB − ηC`,
   `U = ηT + (2ξω + ηζ + k)B + 3ηωC`, `V = ξT + 3ξζB + (ξω + 2ηζ + k)C`; the vector
   `(ξU − ηV, ζU − ωV, U, −V)` is divergence-free (`numerator_divergence`), tangent to `L`
   (`numerator_tangency`), and sends `G` to `L·J` with `J` quadratic in `SB, SC, T², TB, TC, B², BC, C²`
   (`source_with_residual`). The 8 source coefficients' first 8 power-basis coordinates give an 8×8 minor
   with an explicit rational inverse on the good fibre (`GenericSpecialization.sourceInverse`,
   `minor_isUnit`), producing an `O`-linear functional `weights` sending `J ↦ TB`.
5. Regular representation `ρ = Algebra.leftMulMatrix powerBasis : E →ₐ[O] Matrix (Fin 27) (Fin 27) O`;
   `matrixPole = ρ(L)` (27×27 matrix of linear forms in `S,B,C`), `D₀ = det ρ(L)`,
   `N₀,i = columnFunctional weights (ρ(numerator_i) · adjugate ρ(L))`. Differentiating the adjugate
   identity yields `Σ N₀,i ∂_iG = D₀·BT` and `D₀ Σ ∂_iN₀,i = Σ N₀,i ∂_iD₀` (`CubicMatrixCertificate`).
   `D₀` and `N₀,i` are homogeneous of degree 27 in `(S,T,B,C)`; `D₀ ≠ 0` because `ρ(L) = I` at
   `(S,T,B,C) = (1,0,0,0)`.
6. Clearing (`HomogeneousClearing`, `ActiveHomogeneousClearing`, `CubicLiftAssembly.exists_cleared_cubic`):
   finitely many coefficients in `O` are cleared simultaneously, `α = b/a, β = d/c` homogenized by `(ac)^n`,
   giving `Q, P_S, P_T, P_B, P_C` homogeneous of degree `h = 2(n+27)`; on `b = a, d = 2c` the clearing factor
   is `κ·(ac)^n`, `κ ∈ ℚ^×`. Final data: `D = a²c²Q`,
   `N = (a²c²(xP_S − yP_T), 2δac²P_B, 2δca²P_C, a²c²(−yP_S + xP_T), 2δbc²P_B, 2δda²P_C)`.
   Horizontal derivations `x∂x − y∂y, −y∂x + x∂y, a∂a + b∂b, c∂c + d∂d` (chain-rule scales `4δ, 4δ, 2a², 2c²`)
   kill the parameter ring (`HorizontalFractionDerivation`, `FractionDerivation`); `lift_strict_equations`
   gives `directional F N = A·D` and `divergenceNumerator D N = 0`.
7. Legality (`TargetDenominator`, `GenericDenominatorLink`, `NumericDenominatorLink`, `MatrixBezout`,
   `TerminalFibre`): specialize `(a,c,y,b,d) = (1,1,0,1,2)`, `x = X0` free. Then `F ↦ x⁶ + 35x² + 67`
   **[sympy]** (monic in `x`), and `D ↦ κ·det(2x²·I − M)` with `M = M_{ξ+η}` the 27×27 multiplication
   matrix in the specialized algebra `AdjoinRoot GenericSpecialization.P`. An explicit polynomial Bézout
   identity shows `M³ + 140M + 536I` is invertible (`RatioMatrixInverse.lean`, 360 KB of explicit matrix data,
   `actual_matrix_inverse`), and `(2x²)³ + 140(2x²) + 536 = 8(x⁶ + 35x² + 67)`. The generic lemma
   `relPrime_cubic_det (z a b) (M K) (hK : (M³ + a•M + b•1)·K = 1) : IsRelPrime (z³ + az + b) (z•1 − M).det`
   (proved by `relPrime_det_of_matrix_bezout`, a quotient-ring argument) gives fibre coprimality;
   `relPrime_of_specialization` (leading coefficient 1) and `relPrime_of_retract` lift it to `ℂ[X,Y]`.
   `denominator_ne_zero_of_legal` derives `D ≠ 0` from coprimality with the non-unit `F`.
8. Base change ℚ → ℂ: `complexMap : MvPolynomial (Fin 6) ℚ →+* Poly` (injective, preserves homogeneity);
   `ComplexTarget.strict_verify` assembles `Verify targetF targetA 1 d (complexMap D) (complexMap ∘ N) 0`.

**Lean proof architecture.**
- `Solution.lean` (15 lines): `rw [← target_F, ← target_A]; exact QuadricContactSextic.certificate`;
  `TargetBinding.lean` (16 lines) identifies `QuadricContactDenominator.targetF = HodgeSubmission.F` and
  `QuadricContactComplex.targetA = HodgeSubmission.A`.
- `SexticCertificate.lean` (45 lines): `D, N := GenericCubicCertificate.D/N` → `hsrc`, `hdiv` →
  `exists_cleared_cubic D N 27 …` → `⟨n,κ,Q,P,…⟩` → `ClearedDenominatorLegal.denominator_legal Q n κ` →
  `lift_strict_equations` → `refine ⟨2*(n+27)+4, complexMap (denominator Q), …, 0, ?_⟩; apply strict_verify`.
- Module sizes (lines): `RatioMatrixInverse` 343 (360 KB), `GenericMonicPolynomial` 281, `MatrixContact` 257,
  `GenericCoordinateData` 248 (73 KB), `GenericSourceReductions` 245 (887 KB: source reductions mod `P`),
  `CubicMatrixCertificate` 189, `PolynomialMatrixMap` 177, `GenericSpecialization` 171 (270 KB: `P`, factors,
  inverses at `(1,2)`), `GenericCubicLines` 152, `ComplexTarget` 152, `GenericCubicCertificate` 146,
  `HomogeneousClearing` 141, `ConcreteLine` 140, `GenericEliminant` 136, `PolynomialLift` 129, … ,
  `MatrixBezout` 72, `SexticCertificate` 45, `ClearedDenominatorLegal` 42, `Solution` 15.
- Options: `maxHeartbeats 0` (unlimited) in `GenericEliminant`, `GenericSourceCoordinates`,
  `GenericMonicPolynomial`, `GenericMinorInverse`, `GenericCoordinateData`, `GenericSourceFamily`,
  `GenericSpecialization`, `GenericSourceReductions`, and `in`-scoped in `GenericCubicLines`,
  `GenericEliminantEvaluation`, `GenericLineSpecialization`, `ParameterEliminant`;
  `16000000` in `RatioMatrixInverse`, `12000000` in `RatioDataLink`, `4000000` in `ConcreteLine`,
  `2000000` in `Generic/NumericDenominatorLink`, `400000` in `RatioDenominatorLink`;
  `maxRecDepth 100000` widely; `backward.isDefEq.respectTransparency false` in `GenericCubicCertificate`.
- Imports: only `HodgeComputer.Verify`, `HodgeTestbed.Algebra`, `HodgeTarget`, mathlib. **No**
  `HodgeComputer.Constructions.*` (the proof predates the shared library; `formal/README.md` later credits its
  `FractionDerivation` idea as adapted into `HoldoutCover`).
- Tactics discharging identities: `simp only` with explicit derivation lemmas, `ring`, `linear_combination`,
  `module`/`abel` (matrix Bézout), `polynomial`/`compute_degree!` for fibres, `decide` for small matrix facts.

**Related instances.**
- `quadric-contact-residual` (OPEN, discussion #125): same `F`, `A = X0³Y0³`. The exchange `X0 ↔ Y0`
  fixes `A` but negates the Jacobian factor `δ = X0² − Y0²` of the quotient map `(x,y) ↦ (S,T)`; so `A` is not
  a polynomial density through this quotient (`δ ∤ A`). `X1 ↔ Y1` also fixes `A`. The paired-sign cubic
  reduction yields a generic auxiliary elliptic surface with `χ(O) = 3`, Kodaira dimension 1 (not uniruled),
  so a new incidence seed with fixed passive parameters cannot supply a dominant rational-curve family
  (`scripts/audit-jacobian-scope.py`). Open directions named in `why.md`: a transfer serving the
  complementary character, varying the passive parameters, or a different supporting geometry.
- `two-reflection-cubic-surface-d2` (SOLVED, benediamond): the undeformed `F*`, `D = 4(X0² + Y0²)`,
  `N = (−Y0X1², 0, 0, X0X1², 0, 0)` (rotation field, strict). The broken instance shows that after the
  deformation the same quadric is provably dead at degree 2 while a (huge) contact-line denominator exists.
- `polarized-contact-coupled-d2` / `polarized-contact-deformation-d2` (both SOLVED) are the other
  "contact + deformation" pairs in the catalogue (not read in detail here).

---

## 5. `reflection-coupled-fano` ("Reflection · coupled Fano quotient")

**Record.** `q = 1`, `c = 5/7`, `(x0,x1,x2,x3,x4,z) = (X0,X1,X2,Y0,Y1,Y2)`, `u = z²`:

    F = H_c(x, z²),   H_c(x,u) = x0⁶+x1⁶+x2⁶+x3⁶+x4⁶ + u³ + c·P(x,u),
    P = x0³x1x3x4 + x0x1x3³x4 + x0x1x3²u + 2x0x1x4²u + 3x0x2³x3x4 + 2x1³x2x3x4 + x1x2x3x4³,
    A = X0³·Y2³ = x0³z³.

Sparse: Fermat part + 7 coupling monomials with coefficients `5/7, 5/7, 5/7, 10/7, 15/7, 10/7, 5/7`.
[sympy: `F − H_{5/7}(x,z²) = 0`; the **only** coordinate reflection is `Y2 ↦ −Y2`.] `checks.json`:
`coordinateReflections = [5]`, `reflectionsServingTarget = [5]`, projective diagonal group order 4, one junior
element `[3,3,3,3,3,0]` fixing the target, Hessian symmetrizer scalar, smooth mod 5.
Provenance: existence "known by the cited geometric theorem" (Chapter 3 §3.3 reflection theorem: the quotient
`H_c = 0 ⊂ ℙ(1,1,1,1,1,2)` is a smooth Fano hypersurface, anticanonical `O(1)`; `A` odd under the reflection
so its residue descends). `sameVarietyInstances = [reflection-coupled-residual]`; discussion #111;
PR #150 (fqlx, prepared with OpenAI ChatGPT collaborating agents; accepted 2026-10-20 11:52Z; 26 modules,
3289 lines).

**Method name.** Modular free quartic on the reflection quotient + étale/arithmetic lifting + polynomial
family integration + finite trace descent + reflection transport + homogeneous extraction
(HoldoutCover machinery). Nonconstructive: no `D`, no degree, `Q` may be nonzero.

**Mechanism (as proved; numbered equations follow solution.md).**
1. Quotient source `B(x,u) = x0³u`; with `φ(p) = p(x,z²)`, `z·φ(B) = A`.
2. Irreducibility over ℂ of both `H_c` and `F_c` (`CoupledDomain.lean`): both are monic of degree 6 in `x0`;
   specialization `x1 = x2 = x4 = 0, x3 = T+1, u = −1` (resp. `z = i`, since `i⁶ = −1`) gives
   `x0⁶ + (T+1)⁶ − 1`, Eisenstein at `T` in characteristic 0
   (`HoldoutConstruction.irreducible_power_add_constant` + `Monic.irreducible_of_irreducible_map`).
   Normal chart `x0`: `H_c ∤ ∂0H_c`, `F_c ∤ ∂0F_c`.
3. Endpoint section `x1 = 0`: `R_end = ℤ[1/7][r2,r3,r4,r5,ξ]/(ξ⁶ + 3cξr2³r3r4 + r2⁶ + r3⁶ + r4⁶ + r5³)`
   (irreducible: `r2=r3=0, r4=1, r5=T−1` → `ξ⁶ + T³ − 3T² + 3T`, Eisenstein). `trdeg_ℂ K_end = 4`,
   `trdeg_ℂ K = 5` where `K = Frac(ℂ[x,u]/(H_c))`.
4. Characteristic-5 seed (`c ≡ 0` since `7` is invertible): `H_0 = Σx_i⁶ + u³`; with `j² = 2` in `𝔽₂₅`,
   `x0 = j(s⁴+t⁴), x1 = (2+j)(s³t+st³), x2 = js²t², x3 = j(s³t−st³), x4 = s⁴−t⁴, u = s⁸+s⁴t⁴+t⁸`
   (weighted quartic, basepoint-free because `x0 ± jx4 = 2js⁴, 2jt⁴`; endpoint `(j,0,0,0,1,1)`).
   `QuarticSeed.seed_equation` is a `linear_combination` over `j² = 2` and `5 = 0`.
   Deformation: 24 unknowns (coefficients of `t¹..t⁴` in each of the six coordinates; `u` keeps `t⁸`-coefficient 1
   and `t⁵..t⁷` zero), 24 equations (coefficients of `t¹..t²⁴` of `H_0(φ(t))`), gradients
   `G0 = −j(1+t²⁰), G1 = (2−j)(t⁵+t¹⁵), G2 = −jt¹⁰, G3 = −j(t⁵−t¹⁵), G4 = 1−t²⁰, G5 = 3+t⁴+4t⁸+t¹²+3t¹⁶`;
   the Jacobian kernel is zero by sparse elimination (`QuarticMatrix.lean`, `maxHeartbeats 8000000`).
5. Arithmetic lifting (`OpenIncidence`, `ArithmeticLifting`, `ArithmeticSeed`, `ArithmeticFamily`): the
   incidence `S = R_end[v_{i,r}]/(E_1..E_24)` localized at the Jacobian determinant `Δ` and at `v_{1,1}`
   (transversality to `x1 = 0`) is étale over `R_end`, nonzero because the mod-5 seed is a point of it
   (`v_{1,1} = 2 + j ≠ 0`); flatness gives a nonzero étale algebra over `K_end`, hence a field factor `E`,
   finite separable over `K_end`. In Lean: `exists_generic_field : ∃ I : Ideal GenericIncidence, I.IsMaximal ∧
   Module.Finite GenericEndpoint (GenericIncidence ⧸ I) ∧ Algebra.IsSeparable … ∧ (equations vanish) ∧
   IsUnit (Jacobian det) ∧ IsUnit (coordinate (1,0))`, instantiated with `AlgebraicClosure (ZMod 5)` to get `j`.
   The coupling `5/7` is fixed throughout (`ArithmeticBase.Base = Localization.Away (7:ℤ)`, `coupling = 5·7⁻¹`).
6. Dominance (`Lifting.lean`, `TargetEmbedding.lean`): substitution `x_i ↦ φ_i(t)` is injective on
   `ℂ[x1,x2,x3,x4,u]` (constant coefficients `r2..r5` independent; `φ1` has zero constant term and
   `t`-coefficient `v_{1,1} ≠ 0`, "anchored substitution" lemma), extends to `K ↪ L = E(t)`; `L/K` finite
   (equal transcendence degree 5 plus essentially-finite-type; `Algebra.trdeg` API).
7. Polynomial normal `ψ ∈ E[t]⁶` with `Σ G_iψ_i = 1`, `ψ(0) = (g⁻¹,0,…,0)`, `g = G0(0) ≠ 0`
   (`QuarticNormal.lean`); five commuting derivations `D0 = ∂t, D1..D4` (endpoint derivations extended through
   the separable extension, coefficientwise to `E[t]`); matrix `M = [ψ | D0φ | … | D4φ]` with
   `(det M)(0) = g⁻¹v_{1,1} ≠ 0` (`QuarticSource.lean`).
8. Integration/trace/clearing (shared `HoldoutCover` library via `CurveClosure.normalized_data_from_polynomial_family`
   and `FamilyClosure.quotient_normalized_of_quartic_family`): the normal-density polynomial `p(t) ∈ E[t]`
   (formula (13)), char-0 antiderivative, inverse tangent matrix → `w_i ∈ L` solving (14); normalized trace
   `tr_{L/K}/[L:K]` commutes with derivations → solution in `K`; representatives `n_i/d` give
   `D_H = (∂0H_c)·d`, `N_{H,0} = Bd − Σ_{i≥1}(∂_iH_c)n_{i−1}`, `N_{H,i} = (∂0H_c)n_{i−1}`, with
   `Σ N_{H,i}∂_iH_c = B·D_H` and divergence numerator `= H_c·Q_H`. `D_H` legal because `(∂0H_c)·d ≠ 0` in `K`
   and `H_c` is irreducible.
9. Reflection transport (`ReflectionLift.lean`): `D_F = φ(D_H)`, `N_{F,i} = zφ(N_{H,i})` (`i<5`),
   `N_{F,5} = ½φ(N_{H,5})`, `Q_F = zφ(Q_H)`; chain rule `∂_zφ(p) = 2zφ(∂_up)`. Coprimality:
   `squareLast_not_dvd`: if `F_c ∣ φ(D_H)` then `H_c ∣ D_H` (monic in `x0`, coefficient map `u ↦ z²` injective,
   `Polynomial.map_dvd_map`); `squareLast_legal` combines it with `physical_irreducible` via
   `Irreducible.isRelPrime_iff_not_dvd`. Then `hasCertificate_of_normalized` (shared
   `HoldoutCover.HomogeneousCertificate`) performs finite homogeneous extraction (needs `F, A` homogeneous of
   degree 6, `F` prime, `F ∤ ∂0F`) and yields `Verify` at some unspecified `d` with a possibly nonzero `Q̃`.

**Lean proof architecture.**
- `Solution.lean` (37 lines; `maxHeartbeats 1000000`, `maxRecDepth 4096`, `respectTransparency true`):
  `obtain ⟨I, hI, hfinite, _⟩ := Arithmetic.exists_generic_field; … have h := physical_certificate_of_quartic_family
  (Arithmetic.arithmeticFamily I); rw [protected_F, protected_A, protected_q]; simpa only [complexMap_coupling] using h`.
  `SourceBinding.lean` (26 lines): `protected_F : HodgeSubmission.F = hypersurface (5/7) 6 2`,
  `protected_A : HodgeSubmission.A = physicalSource`, `protected_q : HodgeSubmission.q = 1 := rfl`.
- Module sizes (lines): `EndpointDomain` 329, `TargetEmbedding` 310, `QuarticIncidence` 291 (`maxHeartbeats 0`),
  `QuarticSource` 248, `QuarticMatrix` 226 (8,000,000), `CoupledDomain` 195 (1,500,000), `ReflectionLift` 192
  (2,000,000), `OpenIncidence` 147, `QuarticFamily` 126, `ArithmeticFamily` 126 (300,000), `FamilyClosure` 125
  (1,500,000), `CurveClosure` 118, `Lifting` 103, `QuarticSeedJacobian` 100 (0), `ArithmeticLifting` 97
  (2,000,000), `QuarticNormal` 90, `ArithmeticEndpoint` 84 (2,000,000), `ArithmeticBase` 70,
  `QuarticSeedIncidence` 56 (0), `QuarticSeed` 52, `QuotientGradient` 45 (0), `Solution` 37, `QuarticGradient` 36
  (4,000,000), `ArithmeticSeed` 36, `SourceBinding` 26, `QuarticGradientCoefficients` 24.
  `backward.isDefEq.respectTransparency false` in most modules.
- Shared imports: `HodgeComputer.Constructions.HoldoutCover.{HypersurfaceDomain, GenericIncidence,
  ResidueDerivatives, ResidueClearing, PolynomialNormal, PolynomialFieldMap, PolynomialFamilyBridge,
  PolynomialCalculus, HomogeneousCertificate, FiniteResidueClearing, FiniteDifferentials, EvaluationMatrix,
  AlgebraicFields}`, plus `HodgeComputer.Verify`, `HodgeTarget`.
- Coprimality over ℂ: via **irreducibility** (`Irreducible.isRelPrime_iff_not_dvd`), not via specialization
  of `D`; irreducibility itself is by monic specialization + Eisenstein (shared lemma).

**Recorded limitation / failed direction.** A degree-2 map to the characteristic-5 Fermat quotient cannot be
free: the fifth powers of five quadratic coordinates span ≤ 3 functions, products with degree-≤1 variations have
rank ≤ 6, the weight-2 coordinate adds ≤ 4 columns, so the twisted normal map has rank ≤ 10 < 12. This only
obstructs that particular modular construction; it says nothing about characteristic-0 free conics.
The quartic "was found after lower-degree exploration"; finite-field and finite-precision experiments guided it
but are not part of the proof.

**Related instances.**
- `reflection-coupled-residual` (OPEN, discussion #154): same `F`, `A = X0³Y0³`, **even** in `Y2`. The residue
  volume changes sign under the reflection; the density of `u = Y2²` is `2Y2`, so no polynomial source on the
  quotient produces this target; the reflection theorem gives no existence statement. Suggested routes in
  `why.md`: a correspondence involving the branch locus, a different finite cover, or a construction with
  nontrivial reflection character.
- `cyclic-triple-fano` (SOLVED, benediamond) / `cyclic-triple-residual` (OPEN): same Fano-quotient vs.
  complementary-character pattern.
- `one-reflection` (SOLVED, fqlx, bundle assembly with `HasCertificate.add`) and `holdout` (SOLVED, fqlx,
  Fermat cover + moving conic): in `one-reflection` the unique reflection is also `Y2` and the solved target
  `X0³Y0³` is *even* in `Y2`, i.e. the "residual" character was reached there by a cover construction in which
  `Y2` only appears as `Y2⁶`. In `reflection-coupled-fano`, `Y2²` is coupled through
  `x0x1x3²u + 2x0x1x4²u`, which is the obstacle to transplanting that route to `reflection-coupled-residual`.
- `one-reflection-skew` (OPEN): one-reflection sextic with an extra skew coupling; same residual character.

---

## 6. `split-ternary-q2` ("Split ternaries · third-order pole")

**Record.** `q = 2`, `c = 5/7`:

    F = f(X0,X1,X2) + g(Y0,Y1,Y2),
    f(x,a,y) = x⁶ + a⁶ + y⁶ + c(x⁵a + x²a³y + 3xy⁵ + 2a⁵y),
    g(r,z,w) = r⁶ + z⁶ + w⁶ + c(2r⁵z + 2r³zw² + rw⁵ + 3z⁵w),
    A = X0⁴·X1³·Y0³·Y1²     (degree 12, block bidegree (7,5)).

[sympy: `F − (f + g) = 0`.] `checks.json`: no coordinate reflections, Hessian symmetrizer dimension 2 (the
`3+3` split is the unique linear decomposition up to swapping blocks; no further binary reduction), projective
diagonal group order 6, smooth mod 5. Numerator degree `d + 7`. Provenance: existence known by the split
theorem + elliptic contraction (Chapter 3A §§3A.1–3A.3); discussion #113; PR #159 (fqlx, with OpenAI Codex;
accepted 2026-10-20 20:26Z; 27 modules, 7054 lines; `witnessChecked: false`).

**Why the elementary route fails (recorded).** The split Euler field is strict only when the block target
degree is 3 in a three-variable block, independently of `q`; here the block degrees are 7 and 5. Neither
ternary has a coordinate reflection. So the pole-order increase is the substantive content: the certificate must
give full divisibility by `F³`, not only the leading Jacobian-ring component.

**Method name.** Polynomial companion curves on weighted cones (moving conic seed mod 5, arithmetic/étale
lifting), **second-normal-residue** calculus for the third-order pole (new), radial descent, translated
elliptic contraction along a Fermat cubic (the crossed-ternary machinery of `1337isnot1337`, PR #149),
physical transport, simultaneous homogeneous clearing. Nonconstructive, but the final identity is **strict**
(`Q = 0`).

**Mechanism.**
1. Chart: `p = f(x,a,1)`, `q0 = g(r,z,1)`, `s = x⁴a³r³z²`. Companion cones `7(f(x,a,y) + S³ + T³) = 0` and
   `7(g(r,z,w) + S³ + T³) = 0` with weights `(1,1,1,2,2)`.
2. Seed (`CompanionSeeds.lean`, over `ZMod 5` where the `5/7` perturbations vanish): moving conic
   `x = 4X, a = 2 − X², y = 4 + 2X², s = 3x² + 3a², t = 4x² + 2a²` with `x⁶ + a⁶ + y⁶ + s³ + t³ = 0`;
   13 coefficient variables (`slotCoordinate = [1,1,1,2,2,3,3,3,3,4,4,4,4]`, `slotDegree = [0,1,2,1,2,1,2,3,4,1,2,3,4]`),
   13 equations with Jacobian determinant `1` in `𝔽₅`; the 12×14 evaluated matrix has an explicit right inverse
   (`right_inverse : explicit * right = 1 := by decide`), free/normal block determinant `4`.
   `CompanionGeometry.lean`: generic incidence `SquareEtaleRing`/`SquareGenericAlgebra` over
   `ParamField = Frac(ℂ[Fin 3])`, `coefficientIdeal := Classical.choose …`, giving `CoefficientField which`
   (finite over the parameter field), polynomial curves `conicCurve`, polynomial normal `conicNormal` with
   `polynomialPairing gradient normal = 1`, nonzero evaluation determinant.
3. Second normal residue (`HigherResidue.lean`, `NormalTaylor`, `NormalJetBinding`, `JacobianJets`,
   `NormalRowJets`, `HigherCurve`): for a hypersurface `H`, normal coordinate with `g0 = ∂0H`,

       R2(B) = ½( B00/g0² − (3B0H00 + BH000)/g0³ + 3BH00²/g0⁴ ).

   Thicken the curve linearly, `c(t) + εψ(t)` with `dH(ψ) = 1`; write `H = ε + aε² + bε³ + …` and
   (numerator × Jacobian) `= K0 + K1ε + K2ε² + …`; the integrand `P = K2 − 3aK1 + (6a² − 3b)K0` satisfies
   `P = J(R2(B) + ResDiv_{g0} Z)` with an explicit tangent correction `Z`. Polynomial integration in the curve
   parameter, Piola transport, and derivative-compatible normalized trace descend the primitive from the curve
   field to the companion field (`higher_residue_solution_of_density`, `second_density_transgression`,
   `weighted_second_transgression`, `linear_thickening_density_identity`, determinant jet lemmas).
4. First block (degree 7): normal coordinate `T`, `h = x⁴a³`, `B = 49hST` on `H = 7(f + S³ + T³)`, `g0 = 21T²`;
   first-order source `B' = ½∂_T(B/g0) = −7hS/(6T²)` has `R1(B') = R2(B)`; tangent equation RHS `−2hS/(9T⁵)`;
   `T ≠ 0` in the companion field (`CompanionPoleTwo`, `CompanionPoleTwoFamily`, `RationalAmbient`).
   Radial descent `(T,x,a,y,S) = (ρ²T0, ρx0, ρa0, ρ, ρ²S0)`, cone weights `(1,1,1,2)`, normal weight 4, source
   weight `−1` → affine coefficient divided by `ρ` → radial residue gives the affine primitive
   (`CompanionTheta`, `CompanionThetaBridge`, `CompanionPrimitives`). Second block (degree 5) uses the ordinary
   companion construction.
5. Elliptic contraction: adjoin `u³ = p`, `v³ = q0`, Fermat curve `a_e³ + b_e³ + 1 = 0` with
   `δa_e = b_e², δb_e = −a_e²`. First companion scalar `−h·a_e/(3u⁵b_e²)`; the identity
   `δ(a_e²/b_e + b_e²/a_e) = b_e/a_e² − a_e/b_e²` converts it (with base-divergence terms) into the theta source
   `(h/u⁵)·(b_e/a_e²)`; second companion gives the ordinary source `r³z²/v⁴`. In Lean (`ProductDivergence.lean`):
   `alpha := x⁴a³/u⁵`, `beta := r³z²/v⁴`,
   `OrdinaryPrimitive : ∃ bf cr cz, delta bf + coefficient (DY 0) cr + coefficient (DY 1) cz = beta ∧ (kill clause)`,
   `ThetaPrimitive : ∃ bg cx ca, delta bg + coefficient (DX 0) cx + coefficient (DX 1) ca = alpha·thetaRatio ∧ (kill clause)`;
   `product_divergence` (shared `EllipticCompanions.Contraction`) combines them into `U : Fin 4 → K` with
   `Σ D_i U_i = alpha·beta`; multiplication by `u⁵v⁴` and normalized trace give `W : Fin 4 → BaseField` with

       Σ_i ∂_iW_i − 5·W(p)/(3p) − 4·W(q0)/(3q0) = s        (`weighted_divergence`).

6. Physical transport (`HigherPoleTransport`, `PhysicalTransport`): `X = t(x,a,1)`, `Y = tλ(r,z,1)`,
   `F = t⁶(p + λ⁶q0)`; `baseCoefficient = (p + l⁶q)·l⁵/(9p²)`, vertical component and Jacobian correction
   `2V_λ/λ`; `transported_strict_identity` is the strict pole-two chart identity.
7. Clearing and legality (`CertificateClearing.lean`): the 6 degree-zero vector components are related
   simultaneously in `M = Frac(ℂ[X,Y])` (via `θ = algebraMap ∘ monomialMap`, injective) and in
   `HypField = RatFunc(SplittingField(X⁶ + p/q))` (via `ρ`, with `ρ(physicalF) = 0`, a function field of `F = 0`)
   through the relation `Rel θ ρ`; `common_denominator` yields homogeneous `D` of degree `d` with `θD ≠ 0`,
   `ρD ≠ 0`, and `N_i = X2⁷·A_i` of degree `d + 7` (`homogeneous_vector_data`).
   `strict_polynomial_identity : defect physicalF physicalA D 2 N = 0` by injectivity of `θ`, Leibniz for
   `physicalD i (θ(N i)/θ D)`, and `linear_combination (θ D)² * transported_strict_identity`.
   Legality: `ρ D ≠ 0 ⇒ F ∤ D`; `physicalF_irreducible` (rename `(1,2,0,3,4,5)` so `F` is monic of degree 6 in
   `X2`; specialize the coefficients at `![0,0,0,1,X + C a]` with `a` a root of `X⁶ + (15/7)X + 1` and
   `6a⁵ + 15/7 ≠ 0`, giving `X⁶ + C((X+a)⁶ + (15/7)(X+a) + 1)`, Eisenstein-type
   `irreducible_power_add_constant`); `Irreducible.isRelPrime_iff_not_dvd` closes `denominator_legal`.
   `hasCertificate_of_weighted_divergence : Irreducible physicalF → ∀ W, (weighted divergence) → HasCertificate physicalF physicalA 2`
   with `Q = 0`.

**Lean proof architecture.**
- `Solution.lean` (20 lines, no options): `hF : HodgeSubmission.F = SplitTernaryPhysical.physicalF` by
  `unfold … ; simp only [show (5/7:ℂ) = 5*(1/7) by ring, map_mul, map_ofNat]; ring`; `hA` by `unfold; ring`;
  `exact SplitTernaryProduct.physical_certificate`.
- `CompanionAssembly.lean` (53 lines): `ordinary_primitive` from `SplitGElliptic.ordinary_scalar_primitive`,
  `theta_primitive` from `SplitFElliptic.theta_scalar_primitive_complete`, bindings
  `firstPolynomial = SplitFCompanion.p`, `secondPolynomial = SplitGCompanion.p`, `alpha = phi h/u⁵`;
  `physical_certificate := certificate_of_companion_primitives ordinary_primitive theta_primitive`.
- Module sizes (lines; `maxHeartbeats`): `CompanionSeeds` 507 (8,000,000; `maxRecDepth 100000`),
  `CompanionGeometry` 504 (3,000,000), `HigherResidue` 491 (4,000,000), `CompanionTheta` 479 (700,000),
  `CompanionPrimitives` 479 (3,000,000), `CompanionResidue` 462 (3,000,000), `CertificateClearing` 462 (1,500,000),
  `CompanionElliptic` 392 (3,000,000), `PhysicalTransport` 374 (500,000–1,000,000 per section),
  `CompanionScalar` 365 (1,000,000), `CompanionThetaScalar` 356 (1,500,000), `CompanionPoleTwoFamily` 323
  (3,000,000), `ProductDivergence` 311 (400,000–800,000), `HigherCurve` 252 (6,000,000), `CompanionThetaBridge`
  227 (3,000,000), `NormalJetBinding` 154 (3,000,000), `NormalRowJets` 127 (4,000,000), `JacobianJets` 125
  (4,000,000), `HigherPoleTransport` 125, `RationalAmbient` 122 (3,000,000), `CompanionPoleTwo` 98 (3,000,000),
  `EllipticSource` 89, `NormalTaylor` 79 (4,000,000), `CompanionAssembly` 53, `CompanionThetaExistence` 48
  (1,000,000), `SourceDensityAlgebra` 30, `Solution` 20. `backward.isDefEq.respectTransparency false` almost
  everywhere; `maxRecDepth 100000` in the jet modules.
- Shared imports: `HodgeComputer.Constructions.EllipticCompanions.{Residues, CompanionAlgebra, Calculus,
  AffineCompanion, ScalarEmbedding, Kernels, Contraction}`,
  `HodgeComputer.Constructions.HoldoutCover.{ResidueClearing, PolynomialFamilyBridge, PolynomialFamily,
  PolynomialCalculus, Piola, HypersurfaceDomain, HomogeneousField, FiniteDifferentials}`,
  `HodgeComputer.Constructions.PolynomialCurves.ArithmeticLifting`. The second-normal-residue modules were
  subsequently extracted into the shared `HodgeComputer.Constructions.HigherResidues` (per `formal/README.md`).
- Coprimality over ℂ: irreducibility of `F` + `ρ D ≠ 0` in a function field of `F = 0` (not a specialization
  of `D`).

**Related instances.**
- `rank-one-biquartic` (OPEN, discussion #118, q=1, `A = X0³Y0³`): **the same `f` and `g`** plus the coupling
  `(5/7)·L(X)²·M(Y)` with `L = X0² + X0X1 + 2X1X2 + 3X2²`, `M = Y0² + Y0Y1 + 2Y1Y2 + 3Y2²` **[sympy]**
  (bidegree `(4,2)`, tensor rank one, no `(2,4)` part). At `t = 0` this is exactly the split-ternary variety where
  the elliptic-companion theorem supplies every target; at `t = 5/7` the fibre cubic
  `a s³ + b s²u + d u³ = 0` (`a = f, b = tL²M, d = g`) has discriminant `−d(4b³ + 27a²d)` varying with both
  blocks, so the fixed elliptic contraction does not transfer. No reflections, scalar Hessian symmetrizer.
- `crossed-ternary-asymmetric` (SOLVED, 1337isnot1337, q=1, `A = X0⁴X1Y1`): the pole-one origin of the
  companion/contraction machinery reused here.
- Other `q = 2` solves: `fermat-q2` (benediamond), `holdout-q2` (HeadOfNaps), `block-q2` (benediamond).
  The only open `q = 2` target is `cycle-q2` (`A = X0²X1²X2³Y0²Y1²Y2` on the `cycle-three` variety), where a
  higher-residue extension of the cyclic-cubic construction would be the analogue of what this proof did for the
  split case.
- `four-plus-two-residual` (OPEN): `4+2` additive split; the core triple cover `T³ = H ⊂ ℙ(1,1,1,1,2)` has
  trivial canonical class (not Fano), so the `3+3` companion argument does not carry over.
- `polarized-k3-remainder` (OPEN): linearly split `F = ½(f(U) + f(V))` but the involution swaps the blocks;
  the target reaches the `Λ²H²(K3)_prim` remainder.
- `dense-rational` (OPEN): all monomials, no structure; the "no mechanism" baseline for every method here.

---

## 7. Cross-cutting observations (separated by epistemic status)

**Proven (in Lean, kernel-replayed):** all six `HasCertificate` statements; the explicit strict witnesses of §§1–3;
`G = 16F` and the degree-27 line algebra of §4; irreducibility over ℂ of the §5 and §6 sextics; the
second-normal-residue identities of §6.

**Independently checked here [sympy]:** the three ternary witnesses (`E1 = E2 = E = 0`, `gcd(F,D) = 1` over ℚ),
the contact factorizations `−(5/2)x(x−w)²`, `−2x³`, `4u³`, the quadric-contact composition `G = 16F`, the
fibre `x⁶ + 35x² + 67`, the deformation error of the old rotation field, `F = H_{5/7}(x,z²)` with `Y2` the sole
reflection, `F = f + g`, and the rank-one coupling `(5/7)L²M`.

**Computational evidence recorded in the repository but not formalized:** the rank-336/337 obstruction for
`D = X0² + Y0²` on `quadric-contact-broken` (mod 101, with a char-0 argument via a nonzero minor); the
rank-≤10 obstruction for degree-2 modular maps in §5; the moving-conic seed determinants in §6 are formalized
(`decide`).

**Open questions / conjectures suggested by the records (none claimed):**
- The actual minimal denominator degree for `quadric-contact-broken` is unknown; the proof gives `2n + 58` with
  `n` existential. A fixed-`D` linear search with `D = (X0² + Y0²)^k`, `k ≥ 2`, is explicitly not excluded.
- For the ternary family, whether "contact/residual point on the even branch line" is necessary for the
  pulled-back line to serve the given character has not been tested in the records; all four accepted fibres
  were *designed* with this placement.
- For `reflection-coupled-residual` and `quadric-contact-residual`, the obstruction in each `why.md` is to a
  specific quotient-density mechanism, not to existence. The solved `one-reflection` (residual character reached
  by a Fermat-cover/bundle route) is the closest precedent.
- For `rank-one-biquartic`, the natural first experiment is to deform the §6 companion construction in `t`:
  the companion cones become `7(f + tL²M(Y) + S³ + T³)`-type objects with `Y` as parameters; nothing in the
  repository reports an attempt.

**Failed or excluded approaches explicitly documented:** old quadric `X0² + Y0²` for §4 (proven dead);
degree-2 modular map for §5 (rank bound); elementary split Euler field for §6 (wrong block degrees);
Jacobian/junior-quotient mechanisms for the two `*-residual` siblings (symmetry screens in `checks.json`).

## 8. File index

- `/home/kejace/devel/sextics/data/ternary-tangent-residual-h2-d2/{instance.rec,why.md,solution.md,witness.json,Solution.lean,Witness.lean}`
- `/home/kejace/devel/sextics/data/ternary-collision-cb-d2/…`, `/home/kejace/devel/sextics/data/ternary-collision-cbp-d2/…` (same layout)
- `/home/kejace/devel/sextics/data/quadric-contact-broken/{solution.md,Solution.lean,SexticCertificate.lean,ClearedDenominatorLegal.lean,TargetDenominator.lean,MatrixBezout.lean,GenericCubicCertificate.lean,CubicMatrixCertificate.lean,ConcreteLine.lean,CubicLiftAssembly.lean,…}` (56 modules)
- `/home/kejace/devel/sextics/data/reflection-coupled-fano/{solution.md,Solution.lean,SourceBinding.lean,CoupledDomain.lean,ReflectionLift.lean,FamilyClosure.lean,CurveClosure.lean,ArithmeticSeed.lean,QuarticSeed.lean,…}` (26 modules)
- `/home/kejace/devel/sextics/data/split-ternary-q2/{solution.md,Solution.lean,CompanionAssembly.lean,CertificateClearing.lean,ProductDivergence.lean,PhysicalTransport.lean,HigherResidue.lean,CompanionSeeds.lean,…}` (27 modules)
- Shared contract: `/home/kejace/devel/sextics/formal/lean/HodgeComputer/Verify.lean`, `/home/kejace/devel/sextics/formal/README.md`, `/home/kejace/devel/sextics/docs/01-certificates.md`
- Local check script: `/tmp/hodge/scratch/verify_group8.py`
