# Accepted solves, group 7 — catalogue notes

Scope: `mixed-polar`, `one-reflection`, `one-reflection-y2-sector`,
`passive-quadratic-d4-strict-rational`, `polarized-contact-coupled-d2`,
`polarized-contact-deformation-d2`. Repo: `/home/kejace/devel/sextics`
(read only; nothing modified, nothing posted).

Conventions. Variables `Z = (X0,X1,X2,Y0,Y1,Y2)`, `τ(X,Y) = (X,−Y)`. A certificate for
`(F,q,A)` is `(d, D, N, Q)` with `D ≠ 0` homogeneous of degree `d`, coprime to `F`
over `ℂ`, `N_i` homogeneous of degree `d+6q−5`, and

    E := A·D² − q·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F^(q+1)·Q .

Normalized form (`HodgeComputer.normalized_identity`): `E1: q·N(F) = A·D`,
`E2: D·div N − N(D) = F^q·Q`. "Strict" means `Q = 0`.

Status key used below. **Proven** = kernel-checked in the accepted Lean source
(axioms `propext, Classical.choice, Quot.sound`). **Evidence** = my own sympy
re-checks (`/tmp/hodge/scratch/checks.py`, `check_or.py`, sympy 1.12), or
finite-field/numeric checks quoted from the write-ups. **Reported** = claims in
`why.md`/`solution.md`/discussions that I did not independently verify.

Shared Lean verifier facts (all six use them):

- `HodgeComputer.Verify` is a 5-field structure (`denominator_ne_zero`,
  `denominator_homogeneous`, `numerator_homogeneous`, `denominator_legal :
  IsRelPrime F D`, `identity`). `HasCertificate F A q := ∃ d D N Q, Verify …`.
- Coprimality lemmas in `Verify.lean` (81 lines):
  `relPrime_of_retract (f : R →+* S) (g : S →+* R) (hg : g.comp f = id)` pulls
  `IsRelPrime` back along a retraction; `relPrime_of_specialization (φ : R →+* S)
  (hf : IsUnit f.leadingCoeff) (h : IsRelPrime (f.map φ) (g.map φ))` lifts
  univariate coprimality from a specialized fiber when `f` has unit leading coefficient.
- Catalogue status (derived from presence of `acceptance.json`): 66 solved / 16 open.
  Open: off-paired, dense-rational, cycle-q2, polarized-k3-remainder,
  klein-coupled-character, symmetric-sextic-residual, rank-one-biquartic,
  four-plus-two-residual, quadric-contact-residual, one-reflection-skew,
  reflection-coupled-residual, cyclic-triple-residual, short-power-presentation,
  quadratic-cubic-pullback, eight-plane-configuration, diagonal-permutation.

---

## 1. `passive-quadratic-d4-strict-rational` (PR #19, `benediamond`, accepted 2026-10-15, discussion #81)

**Instance.** `q = 1`.

    F = X0⁶ + X1⁶ + X2⁶ + Y0⁶ + Y1⁶ + Y2⁶ + X0²X1²(Y0² + Y2²) + X2²Y1²Y2²
    A = X0·X1·X2·Y1·(4X0² + 5X1²)

Sparse: Fermat + a "two-tie" pattern where one passive square is replaced by the
rank-two quadratic `R = Y0² + Y2²`, and `Y2` also appears in the second tie `X2²Y1²Y2²`.
Symmetry: every monomial of `F` has even exponent in every variable, so `F` is
invariant under all `(ℤ/2)⁶` coordinate sign changes and under `X0 ↔ X1`.
Character of `A` under `(ℤ/2)⁶`: `(−,−,−,+,−,+)`; `A` is *not* `X0↔X1`-symmetric
(the point of the `(4,5)` pair is to pick a direction in a 2-dim target space).

**Method.** Hesse-cofactor contact denominator + logarithmic-Hamiltonian field,
rational norm of two conjugate contact factors. Constructive, strict, witness
exported.

**Witness (proven; re-checked).** `d = 4`, `Q = 0`,

    D = X0⁴ − X0²X1² + X1⁴ + (R/3)(X0² + X1²) + R²/9,      R = Y0² + Y2²
      = u² + v² + w² − uv − uw − vw   at (u,v,w) = (X0², X1², −R/3)      [evidence: checked]
      = (u + ωv + ω²w)(u + ω²v + ωw),  ω = e^{2πi/3}                      [evidence: sympy factor over ℚ(√−3)]
    N0 = X1X2Y1(12X0² + 3X1² + 5R)/18,   N1 = −X0X2Y1(3X0² − 15X1² − 4R)/18,   N2 = … = N5 = 0.

`D` is the cofactor of `u+v+w` in `u³+v³+w³−3uvw` ("Hesse cofactor"); its two
conjugate quadratic factors `Q±` over `ℚ(√−3)` are the contact factors.

**Mechanism (reported + evidence).** With active variables `x = X0, y = X1`, passive
`R`, and multiplier `B = X2Y1`, take `h = λ₊ log Q₊ + λ₋ log Q₋` and the
Hamiltonian field `v = B(h_y ∂x − h_x ∂y)`. A 2-variable Hamiltonian field is
divergence-free, and `B` does not involve `x,y`, so `div v = 0`; clearing gives
`D·div N = N(D)` exactly (E2 ≡ 0). Conjugate residues `λ±` make `v` rational.
The residues are tuned so that `v(F) = A`, i.e. `N(F) = A·D`. The write-up reports
that a direct coefficient-system solve for a strict quartic denominator gives
parameter relations `t₂ − t₃²/3 = 0`, `t₃³ + 1 = 0`, with the rational point
`(t₂,t₃,1) = (1/3,−1,1)` and two conjugate points over `ℚ(√−3)`. (Reported.)
Structural lesson claimed by `why.md`: the contact construction treats `Y0²+Y2²`
as one passive quadratic quantity; it need not be a single square. No claim is
made for arbitrary passive quadratics or all targets.

**Lean architecture.** `Witness.lean` (355 lines, 21 KB) + `Solution.lean` (42 lines).
- `Witness.lean`: `set_option maxHeartbeats 2000000`, `maxRecDepth 1000`. Defines
  `F, A, D, Q, N0…N5` as explicit `MvPolynomial (Fin 6) ℂ` sums; 36 private
  `pderiv i (X j)` lemmas by `norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]`;
  `e1 : 1 * directional F N = A*D` and `e2 : divergenceNumerator D N = F^1*Q` by
  `MvPolynomial.funext; intro x; norm_num <;> ring`; homogeneity lemmas per component.
- Coprimality over ℂ: `toUni : Poly →+* Polynomial Coeff` sends `X2 ↦ Polynomial.X`
  and the other five to constants, `fromUni` the inverse, `retract` proved;
  `specialize := MvPolynomial.eval (1,1,1,1,1)`. Since `D` does not involve `X2`,
  `fiber_D : dFiber = C (25/9)` (a nonzero constant) while `fiber_F = X⁶ + X² + 7`;
  Bezout `(0, 9/25)`; `denominator_legal := relPrime_of_retract … (relPrime_of_specialization specialize …)`
  with `fUni_lc = C 1` a unit. `denominator_ne_zero` by evaluating at `(1,…,1)`.
  `certificate : Verify F A 1 4 D N Q` with `identity := normalized_identity e1 e2`.
- `Solution.lean` (`maxHeartbeats 8000000`, `maxRecDepth 8000`): rewrites the protected
  `HodgeSubmission.F/A` and `HodgeWitness.*` to the local definitions by
  `MvPolynomial.funext; norm_num […] <;> ring`, then `simpa using …certificate`;
  exports `Submission.witness` and `Submission.certificate`. `witnessChecked = true`.
- Denominator family: *Hesse-cofactor quartic* in `(X0², X1², −R/3)`; the same
  template (quartic cofactor, `d=4`) recurs across the `*-d4-*` solved instances.

**Related.** Same family (solved, same author/method): `family-passive-contact`
(`F = F0 + X0²X1²(13/17 Y0² + 5/9 Y0Y2 + 7/11 Y2²) + X2²Y1²Y2²`, passive quadratic
with cross term), `coupled-pair-rational-d4-strict`, `family-sylvester-norm`
(same "average the rational fields with their denominators" principle).
No open instance shares this variety. Not related to the cubic/holdout families.

---

## 2. `polarized-contact-coupled-d2` (PR #39, `benediamond`, accepted 2026-10-16, discussion #102)

**Instance.** `q = 1`. In `U = X + Y`, `V = X − Y` (so `τ` swaps `U ↔ V`):

    r(z)  = z0³ + 2z0²z1 + 3z0z1² + 5z1³                    (binary cubic in z0,z1)
    Q(z)  = z0⁵ + 2z0⁴z1 + 3z0³z1² + 5z0²z1³ + 7z0z1⁴ + 11z1⁵ + z2⁵
    f(z)  = r² + z2·Q                                        (ternary sextic, square contact along z2 = 0)
    F = f(U) + f(V) + 1000003·r(U)·r(V),     A = g(U) − g(V),   g := W(Q)
    W := (−r_{z1}, r_{z0}, 0)                                (plane Hamiltonian field of r)

Dense in the `X,Y` presentation (≈120 terms, coefficients up to ~10⁸); `instance.rec`
`F` and `A` equal these expressions (evidence: checked). Symmetry: only the global
`τ` (= `U ↔ V`); `why.md` reports exact coefficient checks exclude any nontrivial
linear direct-sum decomposition. `1000003` is prime (evidence).

**Method.** Contact Hamiltonian field with two first integrals (`r` and `z2`),
transported through the polarization; the coupling `r(U)r(V)` is annihilated
because each field preserves its own `r`. Constructive, strict, witness exported.

**Witness (proven; re-checked).** `d = 2`, `Q = 0`,

    D = U2·V2 = X2² − Y2²
    (N_X, N_Y) = ( (V2·W(U) − U2·W(V))/2 ,  (V2·W(U) + U2·W(V))/2 ),   i.e. (V2·W(U), −U2·W(V)) in (U,V) coordinates;
    N2 = N5 = 0;  N0,N1,N3,N4 are cubics (listed in witness.json; equal to the formula above — evidence).

Identities (proven in Lean, re-checked): `W(f) = z2·g`; `N(F) = A·D`; `div N = 0`;
`N(D) = 0`. Hence E1 exact and E2 ≡ 0.

**Mechanism.** `W` has zero divergence (Hamiltonian), kills `r` and `z2`, so
`W(f) = z2·W(Q) = z2·g`. On `F`: `W(U)` applied to `f(U)` gives `U2·g(U)`, applied to
`f(V)` gives 0, applied to `r(U)r(V)` gives 0; multiply by `V2` and subtract the
`V`-copy to get `N(F) = U2V2(g(U) − g(V)) = D·A`. `N(D) = 0` because `W` has no
`z2` component and `D = U2V2`. The plane contact geometry: on the double plane
`w² = f`, the line `z2 = 0` lifts to the two components `w = ±r`.
`why.md`: the point is that a coupling built from the *same cubic first integral*
preserves the field even though it destroys visible separation; arbitrary mixed
couplings do not (reported). Reduction mod `1000003` recovers the uncoupled control
(used for finite-field smoothness comparison — reported; smoothness is *not* part of
the certificate).

**Lean architecture.** `Witness.lean` (609 lines, 71 KB — large only because `F` is
written out term by term), `Solution.lean` (42 lines, identical template to #1 modulo
namespace; evidence: diff). `maxHeartbeats 2000000`, `maxRecDepth 4000` in Witness;
`8000000 / 8000` in Solution. Same tactic pattern: `e1`, `e2` via
`MvPolynomial.funext; norm_num <;> ring`. Coprimality: `toUni` sends `X0 ↦ Polynomial.X`;
`D = X2² − Y2²` does not involve `X0`; specialize `(X1,X2,Y0,Y1,Y2) = (1,2,3,4,5)` giving
`dFiber = C(−21)`, `fFiber` = explicit degree-6 polynomial with leading coefficient
`1000005` (a unit); Bezout `(0, −1/21)`; `relPrime_of_retract ∘ relPrime_of_specialization`.
`denominator_ne_zero` by evaluation at `(1,…,6)`. Denominator family: *product of
the two polarized contact coordinates* `U2V2`.

**Related.** `polarized-contact-deformation-d2` (below; same `(D,N)`), the open
`polarized-k3-remainder` (`F = (f(U)+f(V))/2` uncoupled, target reaching the
`Λ²H²(K3)_prim` remainder sector — *different* sector, the contact field does not
reach it), `coupled-pair-rational-d4-strict` (solved; comparison for the
deformation principle), and the open `rank-one-biquartic` as the "split + rank-one
coupling" analogue (`F = f(X)+g(Y)+(5/7)L(X)²M(Y)`, no reflections, no shared first
integral — so this mechanism does not transfer).

---

## 3. `polarized-contact-deformation-d2` (PR #40, `benediamond`, accepted 2026-10-16, discussion #103)

**Instance.** `q = 1`. With `F0, A0, D, N` from #2, `H = Y1⁴`, `s = 101000303 = 101·1000003`:

    F = F0 + s·D·H = F0 + 101000303·(X2² − Y2²)·Y1⁴
    A = A0 + s·N(H) = A0 + 101000303·4Y1³·N4      (N4 = Y1-component of N)

(evidence: both equalities re-checked against `instance.rec`). Symmetry: `τ` only.

**Method.** Deformation by a multiple of the preserved divisor ("moving target"):
same `(D, N)`, `Q = 0`, `d = 2`. Constructive, strict, witness exported.

**Mechanism (proven; re-checked).** Product rule:
`N(F) = N(F0) + s(N(D)H + D·N(H)) = D·(A0 + s·N(H)) = A·D` using `N(D) = 0`;
`div N = 0` and `N(D) = 0` are unchanged, so E2 ≡ 0. Degrees: `N` cubic, `H`
quartic, so `N(H)` has degree 6 and the correct odd parity. The coefficient
vanishes mod both `101` and `1000003`, preserving two finite-field comparisons
with #2 (reported design choice). General principle (proven as the identity, not as
a theorem in the library): if `(D,N)` is strict with `N(D) = 0`, then for *any*
`H` of degree `4` the pair serves `(F0 + sDH, A0 + sN(H))`. Whether the deformed
`F` is smooth is a separate question (reported Hessian check excludes a hidden
linear block decomposition).

**Lean architecture.** `Witness.lean` (623 lines, 74 KB), `Solution.lean` (42 lines,
same template). Same options and same coprimality route as #2: `X0 ↦ Polynomial.X`,
specialization `(1,2,3,4,5)`, `dFiber = C(−21)`, `fUni_lc = C 1000005`, Bezout `(0,−1/21)`.
Denominator family: `U2V2` again.

**Related.** #2 (same `(D,N)`); `coupled-pair-rational-d4-strict` (solved; the
"smaller comparison with the same distinction"). No open instance on this variety.

---

## 4. `one-reflection` (PR #123, `fqlx`, accepted 2026-10-19, discussion #54)

**Instance.** `q = 1`.

    F = X0⁶ + X1⁶ + X2⁶ + Y0⁶ + Y1⁶ + Y2⁶ + X0³X1Y0Y1 + X1³X2Y0Y1
    A = X0³Y0³

Fermat + two couplings sharing the odd pair `Y0Y1` (coefficient 1; compare the
holdout's `X0³X1Y0Y1 + X1³X2Y1Y2`). Symmetries: visible reflection `Y2 ↦ −Y2`;
hidden reflections `P ↦ −P`, `M ↦ −M` for `P = Y0+Y1`, `M = Y0−Y1` (`F` is separately
even in `P, M, Y2` — evidence: checked); `σ = diag(1,−1,1,1,−1,1)` preserves `F`;
coordinate-sign group `{ε : ε0ε1ε3ε4 = ε1ε2ε3ε4 = 1}` (order 16, contains `τ`).
Target decomposition (evidence: checked):
`A = ⅛X0³P(P²+3M²) + ⅛X0³M(M²+3P²) = A₊ + A₋`, `A± = X0³(Y0³ ± Y1³)/2`, with
`σ(A₊) = A₋`. The visible `Y2` reflection does **not** serve `A` (it is even in `Y2`);
existence was already known geometrically from the `P`,`M` reflections (`why.md`,
credited to `fringewidth`, discussion #54).

**Method.** Binary-cubic Hessian (Sylvester) construction on a cubic surface, rational
field averaging over `√Δ`, cofactor (Piola) transport back to six variables,
sign-change assembly of two half-targets via `Verify.add_q1`. Constructive, strict
(`Q = 0`), `d = 60 = 30 + 30`, coefficients in `ℚ(√3)`; **no `witness.json`**
(`witnessChecked = false`), the polynomials are defined structurally in Lean.

**Mechanism.**
1. Chart `X0Y2 ≠ 0`, coordinates `Φ = (u = X0², p = Y0Y1, v = Y2², z = (Y0³−Y1³)/Y2, r = X1/X0, w = X2/X0)`.
   Using `Y0⁶+Y1⁶ = (Y0³−Y1³)² + 2(Y0Y1)³`: `F = G_cub∘Φ` with
   `G_cub = f0 + v³ + vz²`, `f0 = a u³ + b u²p + 2p³`, `a = 1+r⁶+w⁶`, `b = r + r³w`
   (`r,w` passive). `det dΦ = 12(Y0³+Y1³)/X0` (evidence: checked), so the density
   `u²/24` pulls back to `A₊`.
2. Binary cubic identities (evidence: all checked): `Δ = 324a² + 24b³`,
   `h0 = −b²u² + 18aup + 6bp²`, `H = (18au+12bp)∂u + (2b²u−18ap)∂p`, `k = H(f0)/3`:
   `H(h0) = 0`, `H(k) = 3Δf0`, `Δf0² − k² = (4/9)h0³`, and
   `Δu² = 108a ∂_u f0 + 12b² ∂_p f0 − 12b h0`  (1).
3. Over `ℚ(√Δ)` the cubic diagonalizes, `f± = ½(f0 ± k/√Δ)`; with `λ² = 3`,
   `v³ + vz² = ½(v + λz/3)³ + ½(v − λz/3)³`, `r0 = v + λz/3`, `g = r0³/2`. Euler fields
   on the `f₊`-coordinate and on `r0` give a divergence-free field sending `G_cub` to
   `h0` with denominator `3(f₊+g)`. **Averaging the two conjugate rational fields with
   their denominators** (not just numerators) removes `√Δ`:
   `d0 = 9Δg² + 9Δf0g + h0³`, `J = Δ(2g+f0)(u,p,r0,λr0) − k(H_u,H_p,0,0)`,
   `β_h = 3h0J/(4d0)`, and with (1): `β = ((36a d0, 4b² d0, 0, 0) − 3b h0 J)/(8Δ d0)`,
   `β(G_cub) = u²/24`, `div β = 0`. Every numerator component is divisible by `Δ`;
   cancelling it is what yields degree 30 (not 42) for a half target.
4. Physical half-target polynomials (degrees): `p=Y0Y1 (2)`, `s=Y0³−Y1³, S=Y0³+Y1³, T=Y2³, R=T+λs/3 (3)`,
   `c = X0³X1 + X1³X2 (4)`, `G0 = X0⁶+X1⁶+X2⁶, f = G0 + cp + 2p³ (6)` (so `F = f + s² + Y2⁶`),
   `h = −c² + 18G0p + 6cp² (8)`, `δ = 324G0² + 24c³ (12)`, `W = R³ + fT (9)`,
   `Z = R³(R³ + 2fT) (18)`, `D₊ = 9δZ + 4h³T² (30)`,
   `p0 = G0c + 12G0p² + c²p/3 + 6cp³ (10)`, `p1 = −3G0p + c²/3 + cp² + 6p⁴ (8)`,
   `U = 81G0Z − 3chWT + 3hp0T² (24)`, `V = 9c²Z − 3chWpT + chp1T² (26)`, `B = −3chWR (24)`,
   `N₊ = X0³·(3SX0U, 3SX1U, 3SX2U, 6Y1²V + Y0(s+2λT)B, 6Y0²V − Y1(s+2λT)B, 3SY2B)` (31).
   Identities `N₊(F) = A₊D₊`, `D₊ div N₊ − N₊(D₊) = 0` (proven in Lean; evidence:
   exact evaluation at 4 random rational points in `ℚ(√3)` gives 0). Parity: `D₊`
   τ-even, `N₊[0..2]` τ-odd, `N₊[3..5]` τ-even (evidence: checked).
5. The 5-variable core (evidence: both checked symbolically): with
   `𝓔 = 6G0∂_{G0} + 4c∂_c`, `𝓚 = (3R+12T)∂_R + 9T∂_T`, `M = 3𝓔(U) + 18U + 6∂_pV + 𝓚(B) + 6B`:
   `D₊ = 6U(6G0+4cp) + 12V(c+6p²) + 36BR²` and `D₊M = 3U𝓔(D₊) + 6V∂_pD₊ + B𝓚(D₊)`  (7).
6. Assembly: `D₋ = σD₊`, `(N₋)_i = ε_i σ((N₊)_i)`, `ε = (1,−1,1,1,−1,1)` (vector
   transport, signs essential); `D = D₊D₋`, `N = D₋N₊ + D₊N₋`, `Q = 0`.

**Lean architecture.** Single file `Solution.lean` (1205 lines, 46 KB), imports
`HodgeTarget`, `Mathlib.Analysis.Complex.Basic`, `Mathlib.Analysis.Real.Sqrt`,
`Mathlib.Tactic.Ring`. **No** `HodgeComputer.Constructions.*` imports (predates them);
its `Verify.add_q1` / `Verify.sign_change` were later extracted into
`HodgeComputer.CertificateOperations` (`formal/README.md`). Sections:
- `CertificateConstruction` (`namespace HodgeComputer.Construction`): `directional_add`,
  `divergence_add`, `defect_add`, `Verify.add_q1` (combine `(D,N,Q)`,`(E,M,R)` into
  `(DE, EN+DM, E²Q+D²R)`), `signChange s : Poly →+* Poly`, `signField`,
  `pderiv_signChange`, `Verify.sign_change`; `OneReflection.target_decomposition`,
  `signChange_F`, `signChange_Aplus`.
- `OneReflectionDivergenceCore` (`maxHeartbeats 12000000`, `maxRecDepth 6000`):
  `divergence_core (g c p r t : ℂ)` — identity (7b) as a *scalar* identity in five ℂ
  variables, proved by `simp only [defs]; ring`.
- `OneReflectionLegality`: `toUni` (`Y2 = X 5 ↦ Polynomial.X`), `fromUni`, `retract`,
  `monic_F`, `specialize := eval (fun _ => 1)`, `fiber_F = X⁶ + 7`,
  `legal_of_fiber D (h : IsRelPrime (X⁶+7) (fiber D)) : IsRelPrime F D` via
  `relPrime_of_retract` + `relPrime_of_specialization`.
- `OneReflectionPhysical` (`maxHeartbeats 8000000`): `tau := (Real.sqrt 3 : ℂ)`,
  `tau_sq`, the definitions above, `F_decomposition`.
- `OneReflectionDenominator`: `fiber_D = 27972X¹⁸ + 391608X¹² + 953312X⁶`
  (`simp […]; ring`), `fiber_D_remainder : fiber D = (X⁶+7)(27972X¹²+195804X⁶−417316) + 2921212`
  (evidence: checked), `fiber_D_coprime` with explicit Bezout coefficients
  `(−C(1/2921212)·(…), C(1/2921212))` proved by `Polynomial.funext; norm_num; ring`;
  `D_legal`, `D_ne_zero` (evaluate fiber at 1).
- `OneReflectionDifferentiation` (`4000000`): block fields `EX = (X0,X1,X2,0,0,0)`,
  `JY = (0,0,0,Y1²,Y0²,0)`, `KY`; `dir_*` rules; `EX_G = 6G`, `EX_c = 4c`, `JY_P = S`,
  `KY_R = S(3R+12T)`, `KY_T = 9ST`, etc.; `normal_derivative_grouped`, `divergence_grouped`.
- `OneReflectionFirstIdentity` (`8000000`): `first_scalar_identity` (= (7a)) by
  `funext; norm_num; ring`; `three_R_sq` (uses `tau_sq`); `directional_grouped`
  (expands `directional F N` with `pderiv` simp and `ring`); `first_identity :
  directional F N = Aplus * D` via `linear_combination (x 0)^3 * eval x S / 2 * hx`.
- `OneReflectionHomogeneity`: `IsHomogeneous` for each named polynomial (degrees as
  in the table), `D_homogeneous : D.IsHomogeneous 30`, `N_homogeneous … 31`.
- `OneReflectionCertificate`: `sector_verify : Verify F Aplus 1 30 D N 0`
  (`identity := normalized_identity first_identity second_identity`), `sigma`,
  `reflectedN`, `target_sum`, `fullD = D * sigma D`, `fullN i = sigma D * N i + D * reflectedN i`,
  `full_verify : Verify F A 1 60 fullD fullN 0` by `Verify.sign_change` + `Verify.add_q1`
  + `simpa`; `certificate_from_second_identity`.
- `OneReflectionSecondIdentity` (`8000000`, `maxRecDepth 6000`): `eval_*` lemmas,
  `second_identity : divergenceNumerator D N = 0` by rewriting to the grouped form,
  `MvPolynomial.funext`, then `linear_combination (x 0)^3 * eval x S * hc` with
  `hc := divergence_core (eval x G) …`.
- Final: `Submission.certificate := certificate_from_second_identity second_identity`.
- Coprimality over ℂ: proven for `D₊` by the specialization `X0=X1=X2=Y0=Y1=1`,
  `Y2 = t` (nonzero constant remainder `2921212`), transported to `D₋ = σD₊` by the
  automorphism, and to the product since `IsRelPrime` is multiplicative. Note
  `fiber D` was computed with `λ` eliminated (`s = 0` on the fiber, so `R = T`).
- Denominator family: *half-target denominators* `D₊ = 9δZ + 4h³T²` — a mix of the
  binary-cubic discriminant `δ` and the Hessian covariant `h` of `f` with the
  "other" cubic's root coordinate `R`; full denominator `D₊·σ(D₊)`, degree 60.

**Related.** Same variety: `one-reflection-y2-sector` (solved, below; the binary
cubic construction *cannot* serve its transverse target — `why.md`). Open deformation:
`one-reflection-skew` (`F + (1/7)X0X1Y0Y1(Y0²−Y1²)`, same `A = X0³Y0³`; breaks the
`P/M` exchanges, keeps `τ` and the `Y2` reflection, which does not serve the target;
`why.md` records that neither the `y↔T` Fermat cover nor the `(5,1,5)` diagonal
quotient survives — reported). The same target `X0³Y0³` is also open on
`dense-rational`, `rank-one-biquartic`, `reflection-coupled-residual`,
`cyclic-triple-residual`, `quadric-contact-residual`, `off-paired` (different
varieties). Technique lineage: `family-sylvester-norm` (PR #30, field averaging).

---

## 5. `one-reflection-y2-sector` (PR #144, `DiamondJM`, accepted 2026-10-19, discussion #126)

**Instance.** `q = 1`, same `F` as #4, `A = X0³Y0Y1Y2`. `A` is odd in `Y2` and
**even** under both hidden reflections `P ↦ −P`, `M ↦ −M` (since `Y0Y1 = (P²−M²)/4`;
evidence: checked), so only the `Y2` reflection serves it. `checks.json`:
`coordinateReflectionsServingTarget = ["Y2"]`, degree-six Jacobian ranks `[36, 37]`,
smoothness by maximal-coordinate derivative estimate (mixed bounds
`3,4,1,2,2,0` vs pure coefficient `6`), `densityObstruction: "Odd under both deck
transformations of the two binary-cubic maps."` (`provenance.json` still says
"no accepted proof for this target" — stale relative to `acceptance.json`.)

**Method.** Nonconstructive: Fermat cover + `y ↔ T` exchange quotient + moving
conic (polynomial-curve family) + polynomial integration + deck trace + legal
clearing + homogeneous extraction — the holdout architecture of `fqlx` (PR #128)
specialized to equal coupling coefficients. No witness, no degree bound exported.

**Mechanism.**
1. Both couplings contain `Y0Y1`. Adjoin `T` with `T⁶ = Y0⁶ + Y2⁶`, put
   `α = Y0/T`, `β = Y2/T` (`α⁶ + β⁶ = 1`, Fermat curve), `y = Y1`. On the cover
   `G = X0⁶ + X1⁶ + X2⁶ + y⁶ + T⁶ + α·yT·(X0³X1 + X1³X2)`, invariant under `y ↔ T`
   (evidence: checked). With `S = y² + T²`, `P = yT`:
   `H = X0⁶ + X1⁶ + X2⁶ + αP(X0³X1 + X1³X2) + S³ − 3SP²` (evidence: checked) — this
   is the holdout quotient `H(α,β)` with **both** coefficients equal to `α`
   (`weightedHypersurface alpha alpha`, `commonCoverPolynomial alpha alpha`).
2. Seed: `fqlx`'s conic family lies on `H(α_c, β_c)` for every parameter `r`; the
   equal-coefficient condition `β_c = α_c ⇔ k_c² = 1 − r`, with `q_c² = M/(1−r)³`,
   `q_c⁴ = R`, becomes `g(r) = r⁸ + 7r⁷ − 58r⁶ + 149r⁵ − 178r⁴ + 115r³ − 44r² + 10r − 1 = 0`.
   `g(−13) = 35869443`, `g(−12) = −35000185` (evidence: checked; `g` is irreducible
   over ℚ, real roots ≈ `−12.62304`, `0.26761`). For `r < −12` the auxiliary
   quantities `m, R, M, a, b` are sign-definite and both Jacobian polynomials
   `gEven, gOdd` (degree 11) are negative, so the seed is nonsingular.
   `exists_seed` uses the IVT in Lean plus `IsAlgClosed.exists_pow_nat_eq` for
   `q, k, baseBeta` (nonconstructive existence of a `Seed` structure).
3. Source/deck trace: pulled-back density symmetrised and divided by `2(y²−T²)`
   gives the weighted-degree-5 source `B = βX0³P/(2α⁴)`, with
   `J_c·B = A − (X0³Y0Y1³Y2)·T⁻²`, `J_c = 2α⁵(T²−y²)/T` (hand check: trivial).
   The companion `A(y/T)²` is a `ζ⁻²`-eigenvector of the deck automorphism `T ↦ ζT`
   (`ζ` a primitive cube root), so its normalized trace vanishes
   (`commonCover_trace_inverse_square`); the closure recovers the exact source `A`.
4. Primitive: polynomial-curve family on `H`, polynomial normal, invertible
   evaluation matrix ⇒ polynomial integration in the curve parameter (characteristic
   zero), finite trace descends through the finite separable extension, Piola
   pullback through the cover Jacobian, legal clearing, homogeneous extraction.

**How the denominator arises (nonconstructive).** `residue_solution_clearing`
(shared `HoldoutCover/ResidueClearing.lean`): given a rational tangent solution
`w` of the residue-divergence equation on the function field of `f = 0`, pick a
common polynomial representative denominator `d` (`residue_common_representatives`,
`ρ d ≠ 0`) and set `D = ∂₀f · d`, `N₀ = A·d − Σ ∂ᵢf·nᵢ`, `Nᵢ = ∂₀f · nᵢ`; coprimality is
`IsRelPrime f D` from irreducibility of `f` and `f ∤ D`. Then
`hasCertificate_of_normalized` / `homogeneous_data_of_normalized`
(`HomogeneousCertificate.lean`) extract homogeneous data from the non-homogeneous
`(D, N, Q)` using a nonzero normal chart `¬ F ∣ ∂ᵢF` ("no bound on the initial
denominator degree is assumed"; the resulting degree is `5 + d'` for an internal `d'`).
The function-field denominators themselves come from the cover/trace: the norm over
the Fermat-curve coefficient field and the deck group `μ₆`, and the inverse cover
Jacobian (`actualCoverJacobian.det`).

**Lean architecture.** 62 files, 13 717 lines (~0.7 MB). `Solution.lean` (28 lines):

    obtain ⟨p⟩ := HoldoutConic.exists_seed
    obtain ⟨w, hw⟩ := conic_target_residue_solution p coverQuotientSource
    exact holdout_hasCertificate_of_quotient_divergence w hw

All modules are **local copies** under `Proposal.*` (110 local imports; the only
library import is `HodgeComputer.Verify` from `HomogeneousCertificate.lean`) —
the submission predates/ignores the extracted `HodgeComputer.Constructions.HoldoutCover`
library. Diff against the library (evidence: file-by-file diff after rewriting import
prefixes): 22 of 61 modules differ, mostly by `alpha alpha` substitutions
(`CommonCoverBaseChange 14`, `CoverCertificate 27`, `CoverDeckField 33`,
`CoverQuotientCompatibility 26`, `CoverScalars 22`, `CoverSourceBinding 76`,
`CurveSeedIncidence 41`, `CurveSeed 244`, `DeckTrace 33`, `HoldoutCertificateClosure 59`,
`HomogeneousCertificate 20`, `PolynomialNormal 39`, … changed lines); the other 39
are verbatim. Key chain: `CurveSeed` (`g`, `exists_parameter`, `gEven_neg`, `gOdd_neg`,
`parameter_relation`, `Seed`, `Seed.beta_eq_alpha`, `Seed.base_relation`, `Seed.curve`)
→ `CurveSeedIncidence.seedJacobian_nonsingular` → `GenericIncidence`/`ConicCoefficientField`
(étale algebra, lifting) → `ConicPolynomialFamily`/`PolynomialNormal`/`EvaluationMatrix`
→ `ConicResidue.conic_target_residue_solution` → `CoverCertificate.holdout_hasCertificate_of_quotient_divergence`
→ `holdout_hasCertificate_of_quotient_residue` (Piola pullback through
`actualCoverJacobian`, `cover_normal_source_binding`, deck trace
`commonCover_trace_inverse_square`) → `HoldoutCertificateClosure.holdout_hasCertificate_of_cover_residue`
→ `FiniteTrace`/`FiniteResidueClearing`/`ResidueClearing.residue_solution_clearing`
→ `HomogeneousCertificate.hasCertificate_of_normalized`. Irreducibility of `F`
over ℂ: `holdoutPolynomial_irreducible` by `Monic.irreducible_of_irreducible_map`
under the Eisenstein-type specialization `holdoutSpecialization η`
(`X⁶ + C((X+1)⁶ − 1)`-shape), lifted through `finSuccEquiv`. Heartbeats: several
files set `maxHeartbeats 0` (unbounded: `CurveSeedIncidence`, `CurveJacobianStaged`
with `maxRecDepth 100000`, `CurveJacobianScaling`, `CurveGradientBlocks`,
`CoverSourceBinding`, `ConicSourceDerivations`); `CoverCertificate 1000000`,
`ConicCoefficientField 800000`, `EndpointFields 600000`.

**Related.** Same variety: `one-reflection` (#4). Open: `one-reflection-skew`
(the `y ↔ T` exchange is destroyed by the skew term, `why.md`, reported);
`off-paired` (chain + `(1/3)X0²X1X2Y0Y2`: cover term `αβX0²X1X2T²` not `y↔T`
invariant, reported); `cycle-q2`. Lineage: `holdout` (PR #128), `holdout-moving-branch`,
`holdout-ternary-deformation` (replacement conics).

---

## 6. `mixed-polar` (PR #182, `benediamond`, accepted 2026-10-22, discussion #58)

**Instance.** `q = 1`.

    F = X0⁶ + X1⁶ + X2⁶ + Y0⁶ + Y1⁶ + Y2⁶ + ½X0³X1Y0Y1 + ½X1³X2Y1Y2 + ⅓X0X1Y0Y1³
    A = X0·X1·X2·Y0·Y1·Y2      ("balanced" monomial)

Holdout chain with coefficients ½ plus a reverse-edge term of bidegree `(2,4)`.
Symmetry group `Γ` = paired sign changes `(X_i,Y_i) ↦ (−X_i,−Y_i)` (each `i`)
together with the diagonal `μ₆` action with weights `(2,1,0,3,2,1)`
(all three couplings have weight `≡ 0 mod 6`; `A` and `Ω` both have weight 9, so
the residue is `Γ`-invariant). Smoothness over ℂ: maximal-coordinate estimate —
the non-Fermat part of each partial has coefficient-absolute-sum ≤ 7/3 < 6
(evidence: per-variable sums `11/6, 7/3, 1/2, 5/6, 2, 1/2`).

**Method.** Nonconstructive: reverse-edge *transitive cubic Jacobian* construction
(the `cycle-three` / `mixed-biweight` architecture): ternary cubic in squared
coordinates over a passive function field, Weierstrass Jacobian, moving-map
(first-variation) correction, polynomial curve on the Jacobian seeded mod 5,
arithmetic lifting, polynomial-family integration, finite trace, Piola pullback,
normal-chart correction, residue clearing + homogeneous extraction. No witness,
no degree bound.

**Mechanism.**
1. Coordinates on the affine cone: `R = X2²`, `v = Y0/X0`, `u = (Y1/X1)/v`,
   `w = (Y2/X2)/v`, `x = X0²/R`, `y = X1²/R`. Then `F = R³·U(x,y,1)`,
   `A = u·w·v³·R³·x·y` (evidence: both checked) with the six-coefficient cubic
   `U = ax³ + by³ + cz³ + dx²y + ey²z + fxy²`,
   `(a,b,c) = (1+v⁶, 1+u⁶v⁶, 1+w⁶v⁶)`, `(d,e,f) = (uv²/2, uwv²/2, u³v⁴/3)`.
   The `fxy²` term is the reversed first edge (vs. `cycle-three`'s coefficient pattern).
2. Weierstrass model `η² = ξ³ + I₄ξ + I₆`, `I₄ = −3acef + cd²e`,
   `I₆ = −(27/4)a²b²c² − a²ce³ + (9/2)abc²df − ac²f³ − bc²d³ + ¼c²d²f²`;
   `ξ = Θ/H²`, `η = −J/(2H³)` from the covariants `H, Θ, J` with the exact syzygy
   `27J² = 108Θ³ + 36sΘH⁴ − tH⁶ + U·(…)` (`CubicSyzygy.full_syzygy`, proven).
   Invariant differential: with `g = U_x` and the relative derivation `D` (`Dy = 1`,
   `Dx = −U_y/g`): `gDξ = 6η`, i.e. `Dξ/(2η) = 3/g`.
3. Balanced source = first variation `U_τ = U + τxyz` (`cubicJet = U + 2kxyz`, `k = τ/2`);
   at `z = 1` move `ẋ = −xy/g`: `ġ = y − xyU_xx/g`, `İ₄ = (9abc − cdf)/2`, `İ₆ = acde²`.
   Moving-map identity: `ġ/g = (İ₄ξ + İ₆)/(2η²) − (g/3)·D(ξ̇/(2η))`.
4. Jacobian curve: `T = Rv⁶`, `X = R²ξ`, `Y = R³η`; `V = Y² − X³ − A₄(R,T)X − B₆(R,T)`,
   integral form `432V = 432Y² − 432X³ − 54·a8·X − b256` (`a8 = 8A₄`, `b256 = 432B₆`).
   Seed over `𝔽₅` at `(u,w) = (2,1)`: `R(T) = 1 + 3T + 2T²`, `X(T) = 3 + T + T² + 2T³ + 2T⁴`,
   `Y(T) = 2 + 4T + T² + 3T⁴ + 4T⁵ + 4T⁶` (evidence: reduced equation vanishes mod 5);
   13 coefficient equations in 13 unknowns `(R₁,R₂,X₁..X₄,Y₀..Y₆)`, Jacobian determinant
   `1` for `V`, `2` for `432V` mod 5 (`inverse_product … := by decide`, proven).
   `ArithmeticLifting` lifts to a finite separable extension of `ℂ(u,w,R(0),X(0))`;
   the coefficient matrix also supplies a polynomial normal; dominance and finiteness
   via `CoordinateIndependence`.
5. Polynomial source `P_J = (uw/144)(P₄X + P₆)`, `P₄ = (R+w⁶T)[(9/2)(R+T)(R+u⁶T) − (u⁴/12)RT]`,
   `P₆ = (u³w²/8)R²T(R+T)(R+w⁶T)` (`P₄ = R³İ₄`, `P₆ = R⁵İ₆`). Transport: residue-volume
   multiplier `(F_{X0}/(2Y))·det M = −144Rv³`; add `(ξ̇/(6η))·𝒱` (divergence-free
   `𝒱` acting as `uwv³gD`) for the moving-map term; subtract the normal-change vector
   `𝒞` with single `Y0`-component `−(X1X2Y1Y2)X0Y0²/(X0F_{X0} + Y0F_{Y0})`
   (`PhysicalNormalTransport`).
6. Geometric support (reported, not in Lean): the paired-sign quotient surface
   `H = 6(r³+b³)(x³+y³+z³) + 3br²(x²y + y²z) + 2b²rxy²` in `ℙ² × ℙ¹` is smooth
   (evidence: the mod-5 Gröbner check from `solution.md` reproduces `[1]` on all
   three `r=1` charts); sector theorem applies since the residue is `Γ`-invariant.

**How the denominator arises.** `PhysicalClearing.certificate_of_residue_solution`
feeds the residue solution `w` into the shared `residue_solution_clearing`
(`D = ∂₀F·d`, `d` = common denominator of polynomial representatives of `w` on the
quotient `ℂ[Z]/(F)`), then `CyclePhysical.certificate_of_normalized :=
hasCertificate_of_normalized … sextic_irreducible.prime … normal_legal` extracts
homogeneous data. The function-field denominators contain `H, J, g, Y, R, v`,
the finite-trace norm, and the lifting-field norm; nothing is expanded.

**Lean architecture.** 52 files, 46 686 lines, ≈ 8.5 MB (dominated by generated
tables: `CubicIdentityTablesA/B/C` 1.7–1.9 MB each (8381 lines each),
`CubicIdentityData` 881 KB, `CubicDifferentialTables` 662 KB,
`CubicCoefficientExpressions` 199 KB). `Solution.lean` (20 lines): `hF : HodgeSubmission.F
= CyclePhysical.sextic` and `hA` by `unfold; ring`, then `CyclePhysicalPrimitive.certificate`.
Shared-library imports (by file):
`HoldoutCover.{AlgebraicFields, EvaluationMatrix, FiniteDifferentials, FiniteTrace,
HomogeneousCertificate, HypersurfaceDomain, PolynomialCalculus, PolynomialFamilyBridge,
PolynomialFieldMap, PolynomialNormal, ResidueClearing, ResidueDerivatives,
ResiduePiolaPullback}`, `PolynomialCurves.ArithmeticLifting`,
`ResidueTransport.{CoordinateIndependence, Density, MovingMap}`,
`TransitiveCubic.{BivariateDifferentiation, BivariateExpressions, PolynomialDigitCheck,
SignedExpressions}`. Three layers:
- `Cubic*` (covariants): `CubicCovariants` (`H, Θ, J, s4Numerator, s6Numerator`),
  `CubicEncoding` (injective encoding, `PolynomialDigitCheck`), `CubicCoefficientsA–D`
  + `CubicIdentityData/TablesA–C` (coefficient tables; `data_identity0/1` by
  `simp only [huge list]`, `maxRecDepth 131072`, `maxHeartbeats 8000000`),
  `CubicSyzygy.full_syzygy`, `full_syzygyJet` (dual numbers `TrivSqZeroExt`, `k² = 0`),
  `CubicDifferential*` (cleared invariant-differential identity), `CubicJetCalculus`,
  `CubicPoint` (Gaussian-rational point `(2, 2i, 1−i, 1+i, 0, 0)` on `F`, with
  `H = 324 − 81i/2 ≠ 0`, `J ≠ 0`, `g ≠ 0`, by `Complex.ext <;> norm_num`).
- `Jacobian*`: `JacobianModel` (`432V`), `JacobianSeed` (`ZMod 5`, `seed_incidence`
  by `dsimp; ring` to a multiple of 5, `inverse_product := by decide`;
  `maxHeartbeats 8000000 in`, `maxRecDepth 100000 in`), `JacobianSeedBinding`,
  `JacobianIncidence` (`8000000`), `JacobianLifting` (`generic_lift`,
  `complex_generic_factor_exists`), `JacobianFamily` (polynomial normal),
  `JacobianGeometry/Embedding/Dimension` (dominance, finiteness),
  `JacobianNormal` (`4000000`), `JacobianResidue` (polynomial-source primitive),
  `JacobianSource`.
- `Physical*`: `PhysicalField` (`equation` monic in `X0` over `Coeff`;
  `equation_specialization` under `holdoutSpecialization Complex.I` gives
  `X⁶ + C((X+1)⁶ − 1)`, `equation_irreducible` by
  `Monic.irreducible_of_irreducible_map` + `irreducible_power_add_constant`
  (Eisenstein-type), `sextic_irreducible` via `finSuccEquiv`; `CoordinateRing :=
  AdjoinRoot equation`, `FunctionField := FractionRing`, `axis` map proving
  `¬ sextic ∣ ∂₀sextic`, `tangent` derivations, `certificate_of_normalized`),
  `PhysicalChart` (`Localization.Away (X0X1X2Y0)`, nonvanishing transferred from the
  point via `nonzero_of_point`), `PhysicalCoordinates`, `PhysicalCubicJets/Tangent/Volume`,
  `PhysicalJacobianCoordinates/Coefficients/Map` (`3000000`, `maxRecDepth 10000`),
  `PhysicalNormalTransport`, `PhysicalClearing`, `PhysicalPrimitive`
  (`physical_residue_solution`, `certificate`; uses `piola_pullback`,
  `correction_divergence`, `source_transport`, `linear_combination`).
- Coprimality over ℂ: **not** by specialization; `IsRelPrime F D` comes from
  `Irreducible F` (Eisenstein) + `F ∤ D` (because `ρ D ≠ 0` in the residue field),
  inside `residue_solution_clearing`. Heartbeats: max `8000000` (seed/incidence/syzygy),
  `maxRecDepth` up to `131072`; `backward.isDefEq.respectTransparency false` in several files.

**Related.** Same target on the holdout: `holdout-all-six` (solved by the Fermat
cover; `why.md` explains why that cover fails here — the extra term lacks the
`y↔T` partner `yT³`, and the `(5,1,5)` diagonal quotient gives it weight 2; reported).
Architecture siblings (solved): `cycle-three` (`hodgepodgeconjecture`, PR #165),
`mixed-biweight` (`fqlx`, PR #170), `holdout-moving-branch`. Open relatives:
`off-paired` (`F = chain(½) + ⅓X0²X1X2Y0Y2`, `A = X0³Y0³`): breaks the first two
paired signs so the squared-coordinate cubic reduction does not exist (`why.md`,
reported); `cycle-q2` (`q = 2` on `cycle-three`'s `F`, `A = X0²X1²X2³Y0²Y1²Y2`,
residue character `ζ₆²`, erased by averaging to the invariant quotient — reported);
`dense-rational` (no structure). Note `mixed-polar` and `cycle-three` share `Γ`
(weights `(2,1,0,3,2,1)`), which is why the Jacobian route works for both.

---

## Cross-cutting observations

**Denominator families seen in this group.**
1. Hesse-cofactor quartic `u²+v²+w²−uv−uw−vw` in `(X0², X1², −R/3)`, `R = Y0²+Y2²`
   (splits over `ℚ(√−3)`); `d = 4`.
2. Product of the two polarized contact coordinates `U2V2 = X2² − Y2²`; `d = 2`;
   stable under the deformation `F ↦ F + sDH` with `A ↦ A + sN(H)` (any `H` of degree 4).
3. Binary-cubic Hessian half denominators `D₊ = 9δZ + 4h³T²` (`δ` = discriminant of
   the binary cubic `f0`, `h` its Hessian, `Z = R³(R³+2fT)` from the `v³+vz²` part),
   full `D = D₊σ(D₊)`, `d = 60`, coefficients in `ℚ(√3)`.
4. Nonconstructive: `D = ∂₀F · d` for a representative denominator `d` of a
   function-field primitive, then homogeneous extraction (degree not exported).

**Coprimality-over-ℂ techniques.** (a) `relPrime_of_retract ∘ relPrime_of_specialization`:
choose a variable absent from `D` (or in which `F` is monic), specialize the others
to constants, exhibit a univariate Bezout identity (#1, #2, #3 with `D` specializing
to a nonzero constant; #4 with the explicit remainder `2921212`). (b) Irreducibility of
`F` by an Eisenstein-type specialization (`Monic.irreducible_of_irreducible_map`) and
`F ∤ D` (#5, #6).

**Identity-discharge tactics.** Explicit witnesses: `MvPolynomial.funext; intro x;
norm_num <;> ring` (#1–#3) — this evaluates at a point and lets `ring` close a
scalar identity, which scales to degree-6 sextics but not to degree-60 data. #4
instead factors the big identities through five-variable scalar identities
(`divergence_core`, `first_scalar_identity`, `simp only; ring`) and `linear_combination`.
#5/#6 never expand polynomials: they prove field identities (`field_simp; ring`,
`linear_combination`) and rely on shared clearing theorems; #6 additionally
certifies the covariant syzygy by table-driven `simp only` and the seed matrix inverse
by `decide` over `ZMod 5`.

**Proven vs evidence vs reported — summary.**
- Proven (Lean, kernel-replayed): all six `Submission.certificate` theorems with the
  stated `(F,q,A)`; the explicit witnesses of #1–#3; the degree-60 strict certificate
  of #4; existence for #5 and #6.
- Evidence (my sympy re-checks): the three witnesses satisfy `E = 0`, `gcd_ℚ(F,D) = 1`,
  the stated degrees; #1's Hesse-cofactor form and `ℚ(√−3)` factorization; #2's
  construction formulas reproduce `instance.rec` and `witness.json`; #3's
  `F = F0 + sDH`, `A = A0 + sN(H)`; #4's identities (1), (6), (7a), (7b), the
  Jacobian determinant, degrees/parities, and 4-point exact evaluations of the
  half-target identities; #5's `g(−13)`, `g(−12)`, cover invariance, quotient
  identity; #6's coordinate identities, mod-5 seed incidence, mod-5 Gröbner
  smoothness of the quotient surface, derivative bounds.
- Reported, not verified here: the coefficient-system parameter relations in #1;
  the linear-decomposition exclusions and finite-field smoothness comparisons in
  #2/#3; the geometric-support and K3 remarks; all "why the other method fails"
  statements in `why.md` files (these are screens of specific mechanisms, not
  nonexistence proofs — the repository's own caveat).
- Conjectural/open directions implied by this group: whether the deformation rule
  of #3 (`N(D) = 0` strict fields) generalizes to a family statement; whether the
  reverse-edge Jacobian architecture of #6 extends to `off-paired` (paired signs
  broken, so no) or to `cycle-q2` (pole order 2, nontrivial `μ₆` character — the
  shared `HigherResidues` module handles third-order poles only for the split-ternary
  family); whether #4's binary-cubic construction can be re-seeded on
  `one-reflection-skew` (the `P/M` reflections are destroyed, so the half-target
  split `A = A₊ + A₋` with `σ(A₊) = A₋` no longer has an obvious analogue).
