# Accepted solves, group 0: block-q2, cofactor-eight, cofactor-mixed, collision-h3, completed-square-nambu-holdout, cone-ladder-linear

Repository state: `/home/kejace/devel/sextics`, branch `main` at `83fbf5c` (2026-11-01). Variables are always ordered
`Z = (X0, X1, X2, Y0, Y1, Y2)` = Lean `X 0 .. X 5`; `tau(X,Y) = (X,-Y)`.

Certificate contract (formal/README.md, HodgeComputer/Verify.lean): a witness for `(F,q,A)` is `(d, D, N_0..N_5, Q)` with
`D != 0`, `D` homogeneous of degree `d`, each `N_i` homogeneous of degree `d + 6q - 5`, `IsRelPrime F D` in `C[Z]`, and

    A*D^2 - q*D*sum_i N_i dF/dZ_i + F*(D*sum_i dN_i/dZ_i - sum_i N_i dD/dZ_i) = F^(q+1) * Q.

Euler-normalized form used by all six solves (`HodgeComputer.normalized_identity`): `E1: q*N(F) = A*D` and
`E2: D*div(N) - N(D) = F^q * Q`. All six have `Q = 0` (strict), so `v = N/D` is a divergence-free rational field with
`q*v(F) = A` exactly.

## 0. Status legend and what was independently checked

* PROVEN (Lean, kernel replay, axioms `propext, Classical.choice, Quot.sound`): existence of a certificate for each of the
  six exact `(F,q,A)`, with the displayed witness (acceptance.json of each folder, `witnessChecked: true`, Lean 4.33.0,
  docker sandbox without network, 2 CPUs, 14336 MB, 3600 s budget from `.github/solve-policy.json`).
* RE-VERIFIED HERE (exact rational arithmetic, sympy 1.12, script `/tmp/hodge/scripts/verify_group0.py`): E1, E2, the
  full cleared defect, homogeneous degrees, Y-parities, `gcd(F,D) = 1` over Q (gcd is field-extension invariant, so this
  is coprimality over C), and every structural identity quoted below from the six `solution.md` files.
* COMPUTATIONAL EVIDENCE (Singular 4, `/tmp/hodge/scripts/smooth.py`, `absfact.py`): all six sextics are smooth
  (`dim std(jacob F) = 0`, `vdim = 15625 = 5^6`); absolute factorizations of the six denominators.
* COMPUTATIONAL EVIDENCE (`/tmp/hodge/scripts/blockq2_lin.py`): the block-q2 "forced coefficient" claim for the fixed
  denominator (see section 1).
* Conjectures / suggestions are marked CONJECTURE.

No open instance shares a variety with any of the six (checked: exact equality of monic `F` across all 82 catalogue
records, and equality up to any permutation of `{X0,X1,X2}`, of `{Y0,Y1,Y2}`, and the X/Y swap; all six have
`sameVarietyInstances: []` in provenance.json). Relations below are therefore by presentation family or method, not
by variety.

## 0a. Common Lean architecture (identical template in all six folders)

Files per folder: `Solution.lean` (42 lines, 2.6-3.0 KB), `Witness.lean` (342-500 lines, 17-52 KB), `witness.json`
(`d, D, N[6], Q` as exact `**`-syntax strings), `proof.json` (`schemaVersion 2`, `theorem: Submission.certificate`,
`witnessTheorem: Submission.witness`), `solution.md`, plus the record files. No helper modules, no imports from
`HodgeComputer.Constructions.*`; the only library imports are `HodgeComputer.Verify` (which pulls
`HodgeTestbed.Algebra` for `directional`, `divergence`, `divergenceNumerator`) and Mathlib.

Verifier wiring (`scripts/proposals/check.py`, `sandbox.py`): `HodgeTarget` defines `HodgeSubmission.F/A/q` by a
literal translation of `instance.rec` (integer constants become `(n : Poly)`, divisions become `* MvPolynomial.C (1/n : C)`);
`HodgeWitness` defines `HodgeWitness.d/D/N0..N5/N/Q` from `witness.json`; `Witness.lean` is copied to module
`Proposal.Witness`; the checker then proves `HodgeSubmission.checked : HasCertificate ...` and
`HodgeSubmission.checkedWitness : Verify ... HodgeWitness.d HodgeWitness.D HodgeWitness.N HodgeWitness.Q` from the
submitted theorems, replays with the kernel, and audits axioms. Banned tokens include `axiom, sorry, admit,
native_decide, unsafe, elab, macro, syntax`.

`Witness.lean` (namespace `HodgeComputer.Instances.<CamelName>`, `set_option maxHeartbeats 2000000`,
`maxRecDepth 1000`, four linters disabled):

1. `def F, A, D, Q, N0..N5 : Poly` with `MvPolynomial.C (c : C)` coefficients (rational numbers written as `(p/q : C)`),
   `def N : Fin 6 -> Poly := fun i => if i.val = 0 then N0 else ...`.
2. 36 private lemmas `derivative_i_j : pderiv (i : Fin 6) (X j : Poly) = [i = j]`, each by
   `norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]`.
3. `theorem e1 : (q : Poly) * directional F N = A*D` and `theorem e2 : divergenceNumerator D N = F^q*Q`:
   `simp only [directional, divergence, Fin.sum_univ_succ, ...]`, a `change` to the explicit six-term sum,
   `simp only [F, A, D, N0, ..., map_add, pderiv_C, pderiv_mul, pderiv_pow, derivative_*]`, then
   `apply MvPolynomial.funext; intro x; norm_num <;> ring`. So the polynomial identities are discharged as ring
   identities in six complex numbers after evaluation (`MvPolynomial.funext` needs an infinite domain; C qualifies).
   This is the heartbeat-heavy step (a degree-18 identity for block-q2; a 30-term `F` for cofactor-eight).
4. Homogeneity: `denominator_homogeneous`, `N0_homogeneous .. N5_homogeneous` by explicit terms built from
   `isHomogeneous_C`, `isHomogeneous_X`, `.mul`, `.pow`, `.add` (zero components use `isHomogeneous_zero`);
   `numerator_homogeneous` by `fin_cases`.
5. Coprimality over C, the only non-mechanical part:
   * `toUni : Poly ->+* Polynomial Coeff` (`Coeff := MvPolynomial (Fin 5) C`) sends one chosen variable `z` to
     `Polynomial.X` and the other five to constants; `fromUni` is a left inverse; `retract : fromUni.comp toUni = id`
     via `MvPolynomial.ringHom_ext` + `fin_cases` + `simp`.
   * `fUni` is the explicit image of `F`, `toUni_F : toUni F = fUni` by `simp only [...]; polynomial` (Mathlib's
     `polynomial` normalization tactic), `fUni_degree : natDegree = 6` by `compute_degree!`,
     `fUni_lc : leadingCoeff = C 1` (so `F` is monic in `z`: `z^6` occurs with coefficient 1).
   * `point : Fin 5 -> C`, `specialize := MvPolynomial.eval point`, explicit univariate `fFiber, dFiber : Polynomial C`
     with `fiber_F`, `fiber_D` by `norm_num [...] <;> polynomial`, and `fiber_coprime : IsCoprime fFiber dFiber` by an
     explicit Bezout pair `bezoutF, bezoutD` and `polynomial`.
   * `denominator_legal : IsRelPrime F D := relPrime_of_retract toUni fromUni retract
     (relPrime_of_specialization specialize _ _ (unit leading coeff) fiber_coprime.isRelPrime)`.
     Mathematics: a common divisor `g` of `F, D` in `Coeff[z]` divides the monic `F`, hence has unit leading
     coefficient; its specialization divides `gcd(fFiber, dFiber) = 1`; Mathlib's
     `Polynomial.isUnit_of_isUnit_leadingCoeff_of_isUnit_map` makes `g` a unit. No resultants, no irreducibility of
     `D`, one rational point suffices. When `D` does not involve `z`, `dFiber` is a nonzero constant and the Bezout pair
     is `(0, 1/const)`.
6. `denominator_ne_zero` by evaluating at `(1,1,1,1,1,1)` or `(1,2,3,4,5,6)` and `norm_num`.
7. `theorem certificate : Verify F A q d D N Q` with `identity := normalized_identity (by simpa using e1) e2`;
   `theorem solved : HasCertificate F A q`; `#print axioms`.

`Solution.lean` (`set_option maxHeartbeats 8000000`, `maxRecDepth 8000`): proves `hF : HodgeSubmission.F =
HodgeComputer.Instances.<Name>.F` (and `hA, hD, hQ`) by `apply MvPolynomial.funext; intro x; norm_num [defs] <;> ring`,
`hN` by `funext i; fin_cases i` + `change` + the same, then
`simpa only [hF,hA,hD,hQ,hN,HodgeSubmission.q,HodgeWitness.d] using HodgeComputer.Instances.<Name>.certificate`, and
`certificate : HasCertificate ... := <d, D, N, Q, witness>`. The bridge is needed because the generated
`HodgeSubmission` terms differ syntactically from the `Witness.lean` terms.

Per-instance coprimality data (chosen variable `z`, specialization point for the other five variables in order,
fibers):

| id | z | point | fFiber | dFiber | Bezout (bF, bD) |
|---|---|---|---|---|---|
| block-q2 | X1 | (1,1,1,1,1) | X^6 + 6 | 3 | (0, 1/3) |
| cofactor-eight | X0 | (1,1,1,1,1) | X^6+3X^5+10X^4+15X^3+15X^2+8X+8 | X^4+2X^3+2X^2+X+1 | ((-X^2-X+1)/3, (X^4+2X^3+6X^2+5X-5)/3) |
| cofactor-mixed | X0 | (1,1,1,1,1) | X^6+3X^5+13X^4+21X^3+35X^2+25X+23 | 4X^4+8X^3+16X^2+12X+10 | ((-12X^2-12X+8)/89, (3X^4+6X^3+22X^2+19X)/89 - 19/178) |
| collision-h3 | X2 | (1,2,3,4,5) | X^6 + c | -13 | (0, -1/13) |
| completed-square-nambu-holdout | X0 | (1,1,1,1,1) | X^6 + c | 29/4 | (0, 4/29) |
| cone-ladder-linear | X0 | (1,2,3,4,5) | X^6 + c | -1 | (0, -1) |

(collision-h3 must use `z = X2` because `X0^6` has coefficient 3 and `X1^6` coefficient 13; `X2^6, Y1^6, Y2^6` have
coefficient 1. In cone-ladder-linear `X2^6` and `Y2^6` have coefficient -1, so `z = X0`.)

All six PRs are by `benediamond`, merged 2026-10-15/16: #11 (completed-square), #24 (collision-h3), #26 (cofactor-eight),
#27 (cofactor-mixed), #36 (cone-ladder-linear), #38 (block-q2). Discussions: #73, #87, #89, #90, #99, #101 respectively.

Common mechanism (PROVEN by the re-verification): every one of the six fields is logarithmic along `D`, i.e.
`N(D) = lambda * D` and `div N = lambda` for a polynomial `lambda`, so `v = N/D` is divergence-free; and
`q*N(F) = A*D` exactly. Moreover the field is logarithmic along each absolutely irreducible component of `D`
separately (six hyperplanes for block-q2; two conjugate cubics for cofactor-mixed and completed-square, with weight
`A/2` on each; `D` is absolutely irreducible for the other three). Script: `/tmp/hodge/scripts/logcomp.py`.

---

## 1. block-q2 (q = 2)

**Instance.** `F = X0^6 + 2 X0^4 Y0^2 - X0^2 Y0^4 + X1^6 + X2^6 + Y0^6 + Y1^6 + Y2^6` (8 terms; the only coupling is
`X0-Y0`). `A = X0 Y0^11 - (610/87) X0 Y0^5 Y2^6` (degree 12, Y-odd). `F = P(X0^2, Y0^2) + X1^6 + X2^6 + Y1^6 + Y2^6`
with the binary cubic `P(u,w) = u^3 + 2u^2 w - u w^2 + w^3` (`(mu,nu,gamma) = (2,-1,1)`), `disc P = -87`,
`P` irreducible over Q. `F` is a cubic in the six squares, hence coordinatewise even; it is also of the Chapter 2a
shape `P0(X1,X2) + P1(Y1,Y2) + P2(X0^2,Y0^2)` (two Fermat binary sextics plus one square-binary cubic), so by the
book's theorems every Y-odd target at `q = 1, 2` is servable on this variety (constructively). Diagonal symmetries
contain `mu_2^2` (signs of X0, Y0) x `mu_6^4` (scalings of X1, X2, Y1, Y2) and `S_4` permuting the four spectator
sixth powers; `A` is `mu_6^4`-invariant and has character `(-1,-1)` under the two sign changes (odd in X0 and in Y0).

**Method.** Fixed-denominator linear solve ("one block, pole order two"): `D = P(X0^2, Y0^2)` (the active block
itself), numerator degree `13 = 6 + 12 - 5`, solving `2 N(F) = A D`, `D div N - N(D) = 0`.

**Constructive data.** `d = 6`, `D = X0^6 + 2X0^4Y0^2 - X0^2Y0^4 + Y0^6`; over C, `D = prod_{k=1}^{6} (a_k X0 - Y0)` with
`a_k` the roots of `a^6 - a^4 + 2a^2 + 1` (Singular absFactorize): six hyperplanes through the codimension-two linear
space `{X0 = Y0 = 0}`. `Q = 0` (strict). Numerator support: `N0` (10 terms, Y-odd), `N1 = 0`, `N2` (6 terms, each
divisible by X2), `N3` (10 terms), `N4 = -(305/261) X0 Y0^5 Y1 (X2^6 - Y2^6)`, `N5 = -(305/261) X0 Y0^5 Y1^6 Y2`;
note `305/261 = (610/87)/6`. The field involves the spectators `X2^6, Y1^6, Y2^6` but not `X1`.
`div N = N(D)/D = X0 Y0^3 (168 X0^2 X2^6 + 478 X2^6 Y0^2 + 87 Y0^8 - 610 Y0^2 Y2^6)/174`,
`N(F)/D = X0 Y0^5 (87 Y0^6 - 610 Y2^6)/174 = A/2`.

**Mechanism (5-10 lines).** In the active block the sextic is `P(u,w)`, `u = X0^2`, `w = Y0^2`. A degree-13 field with
denominator `P(u,w)` can serve `X0 Y0^11` only together with a compensating spectator term: differentiating along the
spectator Fermat blocks (`N2, N4, N5` are of the form `X2 * (...)`, `Y1 * (...)`, `Y2 * (...)`) produces the
`Y2^6`-shifted monomial `X0 Y0^5 Y2^6`, and the ratio of the two target coefficients is pinned by the cubic:
`10 (2 mu^3 - 9 mu nu + 27 gamma)/disc(P) = 10*61/(-87) = -610/87` (RE-VERIFIED). The witness is one point of a
75-dimensional affine solution space (see below); its structure is not "one Euler or Nambu field" but a genuine
linear-algebra solution; the six-hyperplane denominator makes it logarithmic along every hyperplane `a_k X0 = Y0`
separately (RE-VERIFIED by reduction modulo the minimal polynomial of `a`).

**Fixed-denominator compatibility (COMPUTATIONAL EVIDENCE, exact).** Restricting to the `mu_6^4`-equivariant subspace
(lossless: `F, D, A` are invariant, so averaging any solution gives an equivariant one), with `N0, N3` invariant of
degree 13 and `N1, N2, N4, N5` of the form `Z_i *` (invariant of degree 12): 336 unknown coefficients. For
`A = alpha X0 Y0^11 + beta X0 Y0^5 Y2^6`:
* strict (`E2 = 0`): system 303 x 338, nullspace dimension 75, projection to `(alpha, beta)` has rank 1 with
  `beta/alpha = -610/87`;
* non-strict (`E2 = F^2 Q`, `Q` an invariant sextic, 11 unknowns): system 324 x 349, nullspace 79, again rank 1 with
  `beta/alpha = -610/87`.
So with this `D` neither monomial is servable alone, even modulo `F^3`; exactly the recorded bundle is. This says
nothing about other denominators: the Chapter 2a theorem guarantees some certificate for `X0 Y0^11` alone.

**Lean specifics.** `Solution.lean` 42 lines / 2568 B; `Witness.lean` 364 lines / 25041 B; `e1` is a degree-18
identity. `toUni` uses `z = X1` (so `D` is a constant 3 on the fiber). `denominator_ne_zero` at all-ones (`D = 3`).

**Relations.** Solved relatives: `fermat-q2` (q = 2, subset-Euler field `N = (A/2) E_{X1,X2}`, `D = 6(X1^6 + X2^6)`),
`holdout-q2` (q = 2, nonconstructive, 42 Lean files), `xtype-strict-rational-113` / `xtype-sparse-rational-1m142` /
`twoblock-cubic` (fixed binary-cubic denominators serving "compatible bundles" at q = 1). Open relatives: `cycle-q2` is
the only open q = 2 instance; it is a coupled Fermat deformation (`X0^3 X1 Y0 Y1/2 + X0 X2^3 Y0 Y2/2 + X1^3 X2 Y1 Y2/2`),
not a block form, so neither the 2a theorem nor this denominator applies. CONJECTURE/SUGGESTION: the equivariant
fixed-denominator linear algebra used here (a few hundred unknowns at numerator degree 13) is cheap and could be run for
`cycle-q2` with candidate denominators, but the cyclic coupling leaves a much smaller diagonal symmetry group, so the
reduction will be weaker and a successful denominator family is not suggested by this solve.

---

## 2. cofactor-eight (q = 1)

**Instance.** 30-term `F` (couplings among `X0, X1, X2, Y0` and between `Y1, Y2`): with
`C = X0^2 + X0 X1 + 2 X1^2 + X2^2 - 2 Y0^2` (even quadric, rank 4), `r = Y0 (X0^2 + X0 X1 + X1^2 + X2^2 - (3/2) Y0^2)`
(Y-odd cubic), `K = Y1^6 + Y1^4 Y2^2 + Y2^6`, `H = C^3 + C X1^2 X2^2 + X1^6 + X2^6`:

    F = H(C, X1, X2) + 4 r^2 + K(Y1, Y2)            (RE-VERIFIED)

so `F = H4(X0,X1,X2,Y0) + K(Y1,Y2)` is an additive 4+2 split with a separately-even binary block. Target
`A = J * C * Y1` (20 terms) with the Jacobian `J = C_{X0} r_{Y0} - C_{Y0} r_{X0}`, which factors as
`J = (1/2)(2X0 + X1)(2X0^2 + 2X0X1 + 2X1^2 + 2X2^2 - Y0^2)` (RE-VERIFIED). The target lies in the eight-dimensional
packet `J * <C, X1^2, X1 X2, X2^2> * <Y1, Y2>`.

**Method.** Cofactor transport through the quadratic-cubic map `(C, r)` in the active pair `(X0, Y0)` plus the Euler
field of the passive binary block: this is the Chapter 2a.1 field `e_i + E_j/6` with multiplier `M = C_i L_j`
("C*J times linear" packet), applied with a non-coordinate active "block" (`C, r` involve `X0, X1, X2, Y0`).

**Constructive data.** `d = 6`, `D = r^2 + K/4` (17 terms), absolutely irreducible (Singular), `Q = 0`. Nonzero
numerators (degree 7): `N0 = -(1/8) C Y1 r C_{Y0}`, `N3 = (1/8) C Y1 r C_{X0}`, `N4 = A Y1/24`, `N5 = A Y2/24`;
`N1 = N2 = 0`.

**Mechanism.** The Hamiltonian-type pair `(N0, N3) = (C Y1 r/8)(-C_{Y0}, C_{X0})` annihilates `C` and all of
`X1, X2`, so `N(H) = 0`; it sends `r` to `(C Y1 r/8) J = A r/8`, so `N(4r^2) = A r^2`; the Euler part gives
`N(K) = (A/24) * 6K = A K/4`. Hence `N(F) = A (r^2 + K/4) = A D` and `N(D) = A r^2/4 + A K/16 = A D/4`.
Divergence: write the active pair as `g X_C` with `X_C = (-C_{Y0}, C_{X0})` (Hamiltonian field of `C` in `(X0,Y0)`,
divergence-free, `X_C(C) = 0`, `X_C(r) = J`) and `g = C Y1 r/8`; then `div(g X_C) = X_C(g) = (C Y1/8) J = A/8`. The
Euler part `(A/24)(Y1, Y2)` has divergence `(1/24)(2A + Y1 A_{Y1} + Y2 A_{Y2}) = (1/24)(2A + A) = A/8`, since `A` has
joint degree one in `(Y1, Y2)`. Total `div N = A/4 = N(D)/D` (all four identities `N(C) = 0`, `N(r) = A r/8`,
`N(D) = A D/4`, `div N = A/4` RE-VERIFIED). The Jacobian factor in `A` is the Piola density of the coordinate change;
dropping it changes the class.

**Lean specifics.** `Solution.lean` 2700 B; `Witness.lean` 500 lines / 51509 B (largest of the group; `F` has 30 terms,
`N3, N4, N5` have 20 terms each). `z = X0`; `D` involves `X0`, so `dFiber` is a genuine quartic and the Bezout pair has
degrees 2 and 4 (denominator 3). `denominator_ne_zero` at all-ones (`D = 7`).

**Relations.** Solved relatives (same cofactor/quadratic-cubic mechanism): `cofactor-mixed`, the `*-cofactor-d4-*`
family (`fermat-pair-rational-d4-strict`, `coupled-pair-rational-d4-strict`, `hesse-*-cofactor-d4-rational`,
`twoties-crossblock-hesse-cofactor-d4-rational`, `tie-offblock-ytype-d4-strict-rational`), `family-block-cube`,
`family-sylvester-norm`, `twoblock-cubic`. Open relatives: `quadratic-cubic-pullback` (`F = sum Q_i^3`, a finite map to a
Fermat cubic fourfold; the density `det(dQ)` direction is served by the book's Fano/Piola theorem but the recorded
target `X0^3 Y0^3` is outside it, so the cofactor mechanism, which only ever serves Jacobian-density packets, is
predicted not to reach it); `four-plus-two-residual` (same additive 4+2 shape `H(X0,X1,Y0,Y1) + B(X2,Y2)`, target in
the core; the cofactor method would need a `C^3 + r^2`-type presentation of the core and a target of the form
`J * M`, neither of which is visible). CONJECTURE: none of the open instances has a visible quadratic-cubic presentation
with a Jacobian-density target.

---

## 3. cofactor-mixed (q = 1)

**Instance.** 47-term `F`; `X2` appears only as `X2^6` (a "suspension" by a sixth power), everything else is coupled.
With `C = X0^2 + X0 X1 + 2 X1^2 - 2 Y0^2 + Y1^2 + Y2^2`, `r = Y0 (X0^2 + X0 X1 + X1^2 - (3/2) Y0^2 + Y1^2 + Y2^2)`:

    F = H(C, X1, Y1, Y2) + 4 r^2 + X2^6,   H = C^3 + C X1^2 Y1^2 + X1^6 + Y1^6 + Y2^6   (H found and RE-VERIFIED)

Target `A = J X2^2 Y1` (10 terms), `J = C_{X0} r_{Y0} - C_{Y0} r_{X0} = (1/2)(2X0 + X1)(2X0^2 + 2X0X1 + 2X1^2 - Y0^2 +
2Y1^2 + 2Y2^2)`. The packet here is two-dimensional: `J X2^2 <Y1, Y2>`. Diagonal symmetry: `mu_6` on `X2` (character
of `A` is `zeta^2`), the global sign.

**Method.** Same cofactor transport, with the Euler field of the one-variable block `X2^6` (`N2 = A X2/6`) and
multiplier `M = X2^2 Y1`: Chapter 2a.1's "raw J times (quadratic in block j) times (passive linear)" packet, except that
the passive linear factor `Y1` lives inside `C` and `r` rather than in a separate block. The field still annihilates
`Y1, Y2` (`N4 = N5 = 0`), which is why this works ("passivity is a property of the derivation, not of a coordinate
block").

**Constructive data.** `d = 6`, `D = 4 r^2 + X2^6` (21 terms) `= (2r + i X2^3)(2r - i X2^3)`: two conjugate cubics over
`Q(i)` (Singular: minimal polynomial `2a^2 + 2a + 13`, i.e. `a = (-1 +- 5i)/2`). `Q = 0`. Nonzero numerators:
`N0 = -(1/2) X2^2 Y1 r C_{Y0}`, `N2 = (1/6) A X2`, `N3 = (1/2) X2^2 Y1 r C_{X0}`.

**Mechanism.** `N(C) = 0`, `N(X1) = N(Y1) = N(Y2) = 0` so `N(H) = 0`; `N(r) = A r/2`, hence `N(4r^2) = 4 A r^2`;
`N(X2^6) = A X2^6`; so `N(F) = N(D) = A D`. Divergence: the active pair is `(X2^2 Y1 r/2) X_C` with
`div = X_C(X2^2 Y1 r/2) = (X2^2 Y1/2) J = A/2`, and `d_{X2}(A X2/6) = (A + X2 A_{X2})/6 = (A + 2A)/6 = A/2` because `A`
has `X2`-degree two; total `div N = A` (all RE-VERIFIED). Along each conjugate component
`N(2r +- i X2^3) = (A/2)(2r +- i X2^3)` (RE-VERIFIED): the field is logarithmic on both cubic hypersurfaces with equal
weight.

**Lean specifics.** `Solution.lean` 2700 B; `Witness.lean` 481 lines / 39289 B. `z = X0`, quartic `dFiber`, Bezout
denominators 89 and 178. `denominator_ne_zero` at all-ones (`D = 50`).

**Relations.** As for cofactor-eight. The suspension-by-`X2^6` shape also appears in the open `four-plus-two-residual`
("suspension remainder") and in the solved holdout family; those targets (`X0^3 Y0^3`) are not Jacobian densities.

---

## 4. collision-h3 (q = 1)

**Instance.** `F = 3X0^6 + 9X0^4X1^2 - 7X0^4Y0^2 + 20X0^2X1^4 + X0^2X1^2Y0^2 + 11X0^2Y0^4 + 13X1^6 - 14X1^4Y0^2 -
19X1^2Y0^4 + X2^6 - 18Y0^6 + Y1^6 + Y2^6` (13 terms): `F = p(X0^2, Y0^2, X1^2) + X2^6 + Y1^6 + Y2^6`, a cubic in six
squares (coordinatewise even; the cubic-surface all-target theorem applies in principle). In `(u,w,x) = (X0^2, Y0^2,
X1^2)`:

    p = (u - 2w)^3 + (u - 2w + x) * Q2,   Q2 = 2u^2 + 3uw + 5w^2 + 7ux + 12wx + 13x^2     (RE-VERIFIED)

so the line `L = u - 2w + x = 0` meets the plane cubic `p = 0` only at `(2:1:0)` (on the coordinate line `x = 0`), with
multiplicity three (`p|_{L=0} = (u-2w)^3`): an inflectional ("H3") collision of a tangency point with its residual point.
Target `A = -36X0^5Y0 - 50X0^3X1^2Y0 + 56X0^3Y0^3 - 28X0X1^4Y0 + 28X0X1^2Y0^3 - 44X0Y0^5 = X0 Y0 * q2(u,w,x)` with a
specific quadric `q2`; the target lives in the six-dimensional character space `X0 Y0 * Sym^2(u,w,x)` (odd in X0 and
Y0, invariant under all other sign changes and under `mu_6^3` on `X2, Y1, Y2`).

**Method.** Contact construction in square coordinates: pull back the inflectional line to the quadric
`D = X0^2 + X1^2 - 2Y0^2` and use the logarithmic field it carries.

**Constructive data.** `d = 2`, `D = X0^2 + X1^2 - 2Y0^2` (rank-3 quadric cone with vertex plane `{X0 = X1 = Y0 = 0}`,
absolutely irreducible), `Q = 0`. `N = (-2Y0(X0^2 - 2Y0^2), 0, 0, X0 X1^2, 0, 0)` (degree 3).

**Mechanism.** In `(u,w,x)` the field is `w = -2(u-2w) d_u + x d_w`, which satisfies `w(L) = -2L`, `w(u-2w) = -2L`,
`div_{(u,w,x)} w = -2` (RE-VERIFIED), so `L` is an invariant hypersurface and `w/L` is divergence-free; `N` is half the
cofactor transport of `w` through `(u,w) = (X0^2, Y0^2)` with `x = X1^2` passive. Then `div N = -4X0Y0`,
`N(D) = -4X0Y0 D`, and `w(p) = L (-6(u-2w)^2 - 2Q2 + w(Q2))`, giving `N(F)/D = 2X0Y0(-18u^2 + 28uw - 25ux - 22w^2 +
14wx - 14x^2) = A` (RE-VERIFIED). The spectator sixth powers are untouched. The pole is a quadric because the contact
line is linear in the squares; the degree of `A` in `(u,w,x)` is forced to two by homogeneity.

**Lean specifics.** `Solution.lean` 2656 B; `Witness.lean` 354 lines / 18537 B. `z = X2` (the only `X`-variable whose
sixth power has coefficient 1 together with `Y1, Y2`), point `(X0,X1,Y0,Y1,Y2) = (1,2,3,4,5)`, `dFiber = -13`.
`denominator_ne_zero` at `(1,2,3,4,5,6)` (`D = -27`; `D(1,...,1) = 0`, so all-ones would fail).

**Relations.** Solved relatives with the same mechanism and the same or analogous quadric: `ternary-tangent-residual-h2-d2`
(same `D = X0^2 + X1^2 - 2Y0^2`, ordinary tangent + residual configuration, three-component numerator),
`ternary-collision-cb-d2` (`D = X0^2 - X1^2 - Y0^2`), `ternary-collision-cbp-d2` (`D = -X0^2 + X1^2 - Y0^2`, other
character), `f3-block-x0-quartic-strict` (conic contact, `D = X0^4 - X1^2 Y0^2`), `two-reflection-cubic-surface-d2`
(`D = 4X0^2 + 4Y0^2`). Open relatives: none of the sixteen open instances is coordinatewise even or has a visible
cubic-in-squares block; `quadric-contact-residual` shares only the theme "quadric contact divisor" (its `F` is dense
with broken reflections and its solved partner `quadric-contact-broken` needed a 56-file nonconstructive proof).

---

## 5. completed-square-nambu-holdout (q = 1)

**Instance.** `F = X0^6 + X1^6 + 5X2^6 + (49/4)X2^4Y2^2 + X2^2Y1^3Y2 + 43X2^2Y2^4 + Y0^6 + Y1^6 + 2Y1^3Y2^3 + 23Y2^6`
(10 terms; odd individual powers `Y1^3 Y2^3`, so not coordinatewise even). With `T = Y1^3 + Y2^3 + (1/2) X2^2 Y2`:

    F = X1^6 + T^2 + G,   G = X0^6 + Y0^6 + 5X2^6 + 12X2^4Y2^2 + 42X2^2Y2^4 + 22Y2^6   (G free of X1, Y1; RE-VERIFIED)

`F` is a crossed ternary split `p(X0,X1,Y0) + r(X2,Y1,Y2)` with `p = X0^6 + X1^6 + Y0^6` (so the book's crossed
elliptic-companion all-target theorem applies in principle, nonconstructively), and additionally a Fermat suspension in
`X0, X1, Y0` (diagonal `mu_6^3`). Target `A = X1^2 X2 Y1^2 Y2` (one monomial). The identifier's "holdout" is historical;
this is not the central holdout variety.

**Method.** Weighted radial (Euler) field in the two "coordinates" `s = X1^3`, `t = T`, multiplied by the passive
density `X2 Y2`; equivalently (over `Q(i)`) the trace of two conjugate Nambu channels along `X1^3 +- i T`.

**Constructive data.** `d = 6`, `D = X1^6 + T^2 = (X1^3 + iT)(X1^3 - iT)` (Singular: two conjugate cubics over `Q(i)`,
minimal polynomial `a^2 - 2a + 2`), `Q = 0`. `N1 = (1/6) X1^3 X2 Y1^2 Y2`, `N4 = (1/6) X1^2 X2 Y2 T`, others zero.

**Mechanism.** Write `N = (1/6) X1^2 X2 Y2 * V` with `V = X1 Y1^2 d_{X1} + T d_{Y1}`. Since `dT/dY1 = 3Y1^2`,
`V(X1^3) = 3Y1^2 X1^3` and `V(T) = 3Y1^2 T`, so `V = 3Y1^2 * (s d_s + t d_t)` is the Euler field of `(s,t)` up to the
factor `3Y1^2`; hence `V(D) = 6Y1^2 D` and `N(D) = A D`. `N(G) = 0` because `G` is free of `X1, Y1`, so `N(F) = N(D) =
A D`. Divergence: `div N = (1/6) X2 Y2 (d_{X1}(X1^3 Y1^2) + X1^2 d_{Y1} T) = (1/6) X2 Y2 (3 + 3) X1^2 Y1^2 = A`
(RE-VERIFIED: `N(F) = N(D) = A D`, `div N = A`). On each conjugate factor `N(X1^3 +- iT) = (A/2)(X1^3 +- iT)`
(RE-VERIFIED), which is the "trace of conjugate Nambu channels" statement made precise: the rational field is defined
over Q although the channels are not. Changing the passive `G` preserves the identities; changing the active completed
square does not automatically.

**Lean specifics.** `Solution.lean` 3008 B (long namespace name); `Witness.lean` 348 lines / 19368 B. `z = X0`,
`dFiber = 29/4` (constant, `D` free of `X0`). `denominator_ne_zero` at all-ones.

**Relations.** Solved relatives by mechanism (radial/Euler denominators with passive multipliers): `fermat-q2`,
`tie-singleton-3odd-strict` (`D = 7X1^6/8 + 27X2^6/32 + 27Y0^6/32 + 27Y2^6/32`, "active correction + passive radial
motion"), `family-trace-zero`, `hesse-triangle`, `family-shear` (`D = X1^6 + Y0^6 + Y1^6 + Y2^6`); by family,
`crossed-ternary-asymmetric` (crossed split, nonconstructive elliptic companions). Open relatives: none in the crossed
family; `four-plus-two-residual` (another additive split, 4+2, where the book's 3+3 companions do not apply) and
`short-power-presentation` (power-sum presentation: `F = sum of 8 sixth powers`; a completed square is a 2-power
presentation of a summand) are only thematically related.

---

## 6. cone-ladder-linear (q = 1)

**Instance.** `F = f(X) + f(Y)` with `f(z) = z0^6 + (z1 - z2)(z1^5 + z2^5) = z0^6 + z1^6 - z1^5 z2 + z1 z2^5 - z2^6`
(10 terms; couplings `X1-X2`, `Y1-Y2` only). Contact presentation `f = e^3 + P*k`, `e = z0^2`, `P = z1 - z2`,
`k = z1^5 + z2^5` (RE-VERIFIED). The binary part factors as `(z1 - z2)(z1 + z2)(z1^4 - z1^3 z2 + z1^2 z2^2 - z1 z2^3 +
z2^4)` over Q. `F` is a plain split `f(X) + g(Y)` with `f = g`, so the book's split all-target theorem applies in
principle; it has the extra `X <-> Y` swap. Target `A = -10 X0 Y0 (X1^4 + X2^4)`.

**Method.** Hyperplane rung of the "contact ladder" `f = c e^3 + P k`, `1 <= deg P <= 5`: the Nambu/cross-product
field `W = grad P x grad e` in the `X`-slots times a passive odd multiplier `Y0`.

**Constructive data.** `d = 1`, `D = X1 - X2` (a hyperplane; the smallest possible pole), `Q = 0`,
`N = (0, -2X0Y0, -2X0Y0, 0, 0, 0)` (degree 2).

**Mechanism.** `W = grad P x grad e = (0, -2X0, -2X0)` annihilates `P` and `e` and has zero divergence, so
`W(f) = P W(k)`, hence `N(F) = Y0 W(f) = D * Y0 W(k) = A D` with `A = Y0 W(k) = -10 X0 Y0 (X1^4 + X2^4)`;
`N(D) = 0` and `div N = 0` (RE-VERIFIED), so `v = N/D` is divergence-free and tangent to its own polar hyperplane.
Legality: on `X1 = X2`, `F|_{D=0} = X0^6 + f(Y)`, not identically zero, so `gcd(F,D) = 1`. The same identity works for
every rung of the ladder; the degree of the pole equals `deg P`.

**Lean specifics.** `Solution.lean` 2766 B; `Witness.lean` 342 lines / 17317 B (smallest). `z = X0`, point
`(X1,X2,Y0,Y1,Y2) = (1,2,3,4,5)`, `dFiber = -1`. `denominator_ne_zero` at `(1,2,3,4,5,6)` (`D = -1`; all-ones gives 0).

**Relations.** Solved relatives: `contact46006` (cubic rung: `f = h^2 + C k` with `deg C = 3`, polarized to `U = X+Y`,
`V = X-Y`, `D = C(U)C(V)`, `d = 6`), `family-passive-contact` ("passive quadric, linear contact", `d = 6`),
`polarized-contact-coupled-d2` / `polarized-contact-deformation-d2` (`D = X2^2 - Y2^2`, contact fields with two first
integrals), `f3-block-x0-quartic-strict`, `twoties-tie2-x2-quartic-strict`. Open relatives: `polarized-k3-remainder`
(polarized split `f(U) + f(V)` whose involution swaps the summands, the case where the split theorem fails; the target
reaches the K3 remainder sector), `rank-one-biquartic` (a rank-one coupling deformation away from Fermat/split; the
contact presentation and the passive multiplier do not survive generic coupling), `eight-plane-configuration` (planes
`Y_i = +- i X_i`, a contact-type arrangement inside a dense 236-monomial equation). CONJECTURE: a hyperplane or
low-degree contact denominator for the open dense instances would require an exact presentation `F = c e^3 + P k + ...`
with the passive part annihilated by the field; none of the open records displays one.

---

## 7. Cross-cutting summary

| id | q | d | D (structure over C) | #N nonzero | method | Lean (Solution/Witness bytes) |
|---|---|---|---|---|---|---|
| block-q2 | 2 | 6 | 6 hyperplanes through {X0=Y0=0} (sextic field a^6-a^4+2a^2+1) | 5 | fixed-denominator linear solve, binary block | 2568 / 25041 |
| cofactor-eight | 1 | 6 | absolutely irreducible `r^2 + K/4` | 4 | cofactor transport, (C,r) map + Euler on K | 2700 / 51509 |
| cofactor-mixed | 1 | 6 | two conjugate cubics `2r +- i X2^3` | 3 | cofactor transport, (C,r) map + Euler on X2^6 | 2700 / 39289 |
| collision-h3 | 1 | 2 | rank-3 quadric cone | 2 | inflectional contact line in square coordinates | 2656 / 18537 |
| completed-square-nambu-holdout | 1 | 6 | two conjugate cubics `X1^3 +- iT` | 2 | weighted radial in (X1^3, T), passive density X2 Y2 | 3008 / 19368 |
| cone-ladder-linear | 1 | 1 | hyperplane | 2 | Nambu field `grad P x grad e` times passive Y0 | 2766 / 17317 |

Facts common to all six (PROVEN/RE-VERIFIED): strict (`Q = 0`), Euler-normalized (`q N(F) = A D`), logarithmic
(`N(D) = (div N) D`), numerators supported in at most five slots, all data over Q, `gcd(F,D) = 1`, `F` smooth.
Each `F` belongs to at least one of the book's all-target families (coordinatewise-even: block-q2, collision-h3;
binary blocks: block-q2; split: cone-ladder-linear; crossed: completed-square), or is a cofactor/density instance with
a small served packet (cofactor-eight: 8 directions; cofactor-mixed: 2 directions). The six accepted proofs are explicit
witnesses, not the general theorems.

What these solves do not establish: a denominator that fails for a different target ratio (block-q2) says nothing
about other denominators; a packet served by cofactor transport (cofactors) says nothing about targets outside the
Jacobian-density span; a contact field exists only for the displayed exact presentation (collision-h3, cone-ladder).
None of the sixteen open instances shares a variety with these six, and all sixteen have dense couplings (coupling
graph with 14-15 of the 15 possible edges, except `four-plus-two-residual` with 7 edges and `one-reflection-skew` with
9), so the sparse presentations exploited here are not visible in any open record.

Reusable Lean pattern (not library code, but a template): generated `Witness.lean` with `MvPolynomial.funext +
norm_num + ring` for the identities, explicit `isHomogeneous_*` terms for degrees, and the
`relPrime_of_retract (relPrime_of_specialization ...)` chain for coprimality from one rational specialization point.
Requirements for reuse: `F` monic in some variable (every catalogue `F` contains some `Z_i^6`; rescale if the
coefficient is not 1), rational coefficients, and a point where the univariate specializations of `F` and `D` are
coprime (generic; Bezout cofactors from the extended Euclidean algorithm). Heartbeats: 2000000 sufficed for a
degree-18 identity with 50 KB of definitions; 8000000 in the bridging file.

Scripts and logs: `/tmp/hodge/scripts/verify_group0.py`, `blockq2_lin.py`, `absfact.py` (+ `abs_*.sing`),
`smooth.py` (+ `smooth_*.sing`), `logcomp.py`, `samevar.py`.
