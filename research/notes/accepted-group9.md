# Accepted solves, group 9: catalogue notes

Scope: six accepted instances of the hodge.computer sextic challenge, read from
`/home/kejace/devel/sextics/data/<id>/` (instance.rec, why.md, solution.md,
proof.json, witness.json, acceptance.json, Solution.lean, helper .lean files).
Variables are ordered `X0, X1, X2, Y0, Y1, Y2` (Lean indices 0..5). The
certificate contract (formal/README.md, docs/01-certificates.md):

    E := A D^2 - q D sum_i N_i dF/dZ_i + F ( D sum_i dN_i/dZ_i - sum_i N_i dD/dZ_i ) = F^(q+1) Q,

with D nonzero homogeneous of degree d, coprime to F over C, and each N_i
homogeneous of degree d + 6q - 5 (zero allowed). Normalized form (sufficient,
used by every Lean proof here): E1: q N(F) = A D and E2: D div N - N(D) = F^q Q
(`HodgeComputer.normalized_identity`). "Strict" means Q = 0.

Legend used below: **[proven]** = established by the accepted Lean proof or by a
derivation written out here; **[computed]** = re-verified by me with sympy 1.12 /
Singular (scripts in `/tmp/hodge/scripts/`); **[claim]** = asserted in the
solution text and not independently checked; **[observation]** = my own
derived statement, checked computationally where stated.

Common acceptance data: all six were submitted and accepted from
`benediamond`'s PRs on 2026-10-15/16 (toolchain `leanprover/lean4:v4.33.0`,
sandbox `docker-no-network`, kernel replay, axioms `propext`, `Quot.sound`,
`Classical.choice`). None of the six imports any
`HodgeComputer.Constructions.*` module; the five small ones import only
`HodgeComputer.Verify` (plus the generated `HodgeWitness`/`Proposal.Witness`),
the large one imports `HodgeComputer.Verify`, `HodgeTestbed.Algebra` and Mathlib.

---

## 0. The shared constructive template (five of the six)

Five instances (three-block, tie-offblock, tie-singleton, twoblock-cubic,
two-reflection) use one generated proof template, identical modulo the
instance namespace (`diff` of the Solution.lean files is empty after renaming).

**Files.** `Solution.lean` (2.6-3.1 KB) and `Witness.lean` (18-21 KB, 340-350
lines) plus `witness.json`; `proof.json` has `"theorem": "Submission.certificate"`
and `"witnessTheorem": "Submission.witness"` (schemaVersion 2).

**Solution.lean.** `import HodgeWitness`, `import Proposal.Witness`,
`set_option maxHeartbeats 8000000`, `set_option maxRecDepth 8000`.
`Submission.witness : Verify HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q
HodgeWitness.d HodgeWitness.D HodgeWitness.N HodgeWitness.Q` is proved by binding
each protected polynomial to the instance copy:
`apply MvPolynomial.funext; intro x; norm_num [HodgeSubmission.F, Instances.<Name>.F] <;> ring`
(for F, A, D, Q; for N via `funext i; fin_cases i` then `change` + the same),
then `simpa only [hF,hA,hD,hQ,hN,HodgeSubmission.q,HodgeWitness.d] using
Instances.<Name>.certificate`. `Submission.certificate` is the existential
packaging `⟨d, D, N, Q, witness⟩`.

**Witness.lean** (`namespace HodgeComputer.Instances.<Name>`,
`maxHeartbeats 2000000`, `maxRecDepth 1000`, unused-tactic linters off):

1. `def F, A, D, Q, N0..N5 : Poly` as explicit `MvPolynomial (Fin 6) ℂ`
   expressions with `MvPolynomial.C (c : ℂ)` coefficients; `N i` by nested `if`.
2. 36 private lemmas `derivative_i_j : pderiv i (X j) = [i = j]` via
   `norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]`.
3. `e1 : (1 : Poly) * directional F N = A*D` and
   `e2 : divergenceNumerator D N = F^1*Q`:
   `simp only [directional, Fin.sum_univ_succ, ...]; change <explicit sum>;
   simp only [F,A,D,N0..N5, map_zero, map_add, pderiv_C, pderiv_mul, pderiv_pow,
   derivative_*]; apply MvPolynomial.funext; intro x; norm_num <;> ring`.
   So the polynomial identities are discharged by pointwise evaluation over C
   plus `ring`, not by `linear_combination` or `decide`.
4. Homogeneity: `denominator_homogeneous : D.IsHomogeneous d` and
   `Ni_homogeneous : Ni.IsHomogeneous (d+1)` by explicit term-by-term
   combinators `simpa only [D] using ((isHomogeneous_C ..).mul ((isHomogeneous_X ℂ i).pow k)).add ...`;
   zero components by `isHomogeneous_zero`.
5. Coprimality over C (`denominator_legal : IsRelPrime F D`):
   `toUni : Poly →+* Polynomial Coeff` (Coeff = `MvPolynomial (Fin 5) ℂ`) sends one
   chosen variable t to `Polynomial.X` and the other five to constants;
   `fromUni` is its left inverse (`retract : fromUni.comp toUni = RingHom.id`,
   by `ringHom_ext` + `fin_cases`). Then
   `relPrime_of_retract toUni fromUni retract` and
   `relPrime_of_specialization specialize` with `specialize = MvPolynomial.eval point`
   at the all-ones point. Obligations: `fUni_lc` (leading coefficient of toUni F
   is the unit `C 1`, via `compute_degree!`) and `fiber_coprime : IsCoprime fFiber dFiber`
   with explicit Bezout cofactors. In all five instances the chosen t is absent
   from D, so `dFiber` is a nonzero constant c and the Bezout pair is
   `(bezoutF, bezoutD) = (0, C (1/c))`. The variable choices and fibers:

   | instance | t | fFiber (F at others = 1) | dFiber |
   |---|---|---|---|
   | three-block | X0 | t^6 + 4 | 1 |
   | tie-offblock | X2 | t^6 + 6 | 32 |
   | tie-singleton | X0 | t^6 + t^2 + 5 | 109/32 |
   | twoblock-cubic | X2 | t^6 + 7/3 t^4 - 2 t^2 + 29/8 | 205/17 |
   | two-reflection | X1 | t^6 + 10 t^2 + 10 t + 25 | 8 |

   **[proven]** This is a complete proof over C: `relPrime_of_specialization`
   (Verify.lean) shows that a common factor would have unit leading coefficient
   in t and would survive specialization, contradicting coprimality of the fibre.
6. `denominator_ne_zero` by evaluating D at (1,..,1) (`norm_num` on `congrArg eval`).
7. `certificate : Verify F A 1 d D N Q` assembled with
   `identity := normalized_identity (by simpa using e1) e2`; `solved : HasCertificate F A 1`.

**Independent check [computed].** For all five witness.json tuples, sympy
confirms E1 = 0, E2 = 0 with Q = 0 (so E = 0 identically), the degree table,
the component parities (N_0..N_2 Y-odd, N_3..N_5 Y-even, D Y-even), and
gcd(F, D) = 1 over Q (`/tmp/hodge/scripts/verify_witness.py`). In every one of
the five, N(D) = D div N holds separately, i.e. the rational field v = N/D is
divergence-free and satisfies v(F) = A exactly.

---

## 1. three-block-rational-d4-strict

**Data.** q = 1.
F = X0^6 + [6 X1^6 - 15 X1^4 Y1^2 + 9 X1^2 Y1^4 + Y1^6] + X2^6 + Y0^6 + Y2^6
(three binary blocks (X_i, Y_i); first and third diagonal, middle non-diagonal
and even in each variable). A = X0 X1 X2^3 Y1. Discussion #69. PR #7
(benediamond, accepted 2026-10-15T18:13Z).

**Method.** "Cyclotomic cofactor across blocks" (solution.md title). Explicit,
constructive, strict.

**Witness.** d = 4,
D = X1^4 + X1^2 X2^2 - 4 X1^2 Y1^2 + X2^4 - 2 X2^2 Y1^2 + 4 Y1^4
  = L^2 + L X2^2 + X2^4 with L = X1^2 - 2 Y1^2
  = (L - ω X2^2)(L - ω^2 X2^2) = N_{Q(ω)/Q}(L - ω X2^2), ω^2+ω+1 = 0 **[computed]**.
D is irreducible over Q, splits over Q(√-3) into two conjugate quartics.
N = (0, X0 X2^3 Y1 / 6, X0 X1 Y1 (L + X2^2) / 6, 0, X0 X1 X2^3 / 6, 0),
Q = 0. Degrees: D 4, N_i 5. div N = X0 X1 X2 Y1 / 3 and N(D) = D div N **[computed]**.

**Mechanism [proven, derivation].** Put u = X1^2, w = Y1^2. The middle block is
7 (u - w)^3 - L^3 with L = u - 2w **[computed]**, so F restricted to the active
variables is X0^6 + 7(u-w)^3 - (L^3 - X2^6) + (inert). Since L^3 - X2^6 =
(L - X2^2)(L^2 + L X2^2 + X2^4) = (L - X2^2) D, the denominator is the cofactor
of the "cube minus cube" factorization (hence "cyclotomic"). The field
differentiates only X1, X2, Y1:
N(F) = (1/3) X0 X1 X2^3 Y1 (F_u + F_w) + X0 X1 X2^5 Y1 (L + X2^2).
The operator ∂_u + ∂_w kills 7(u-w)^3 and gives (∂_u+∂_w)(-L^3) = 3L^2, so
N(F) = X0 X1 X2^3 Y1 (L^2 + L X2^2 + X2^4) = A D. The divergence identity
N(D) = D div N is a direct computation ((L+X2^2)(L+2X2^2) - X2^2(2L+X2^2) = D).
So v = N/D is a divergence-free field with v(F) = A: strict.
The norm over Q(ω) is what keeps the field rational (why.md asks exactly whether
the "sum of cubes" structure can be realized over Q; answer: yes, through the
norm of the Eisenstein-integer factor).

**Symmetry [computed].** F has all six pure sixth powers; the diagonal symmetry
group modulo scalars is the torsion of Z^6/L(F) = Z ⊕ Z/2 ⊕ (Z/6)^4, order 2592.
All 63 coordinate sign changes fix F (F is a polynomial in the squares).
23 nontrivial variable permutations fix F (permutations of {X0, X2, Y0, Y2},
one block symmetric). A is odd exactly under sign patterns with an odd number of
flips among {X0, X1, Y1} and X2 (as listed by the script); A is odd in X0, X1,
X2, Y1 individually.

**Lean architecture.** Template of §0: Solution.lean 2986 B, Witness.lean
18815 B (341 lines); t = X0; fFiber = X^6 + 4, dFiber = 1. Verify
degree d = 4, numerators 5.

**Related instances.** Same A = X0 X1 X2^3 Y1 appears in
`hesse-radial-hesse-cofactor-d4-rational` (PR 16, d = 4 strict, D a quartic
in X0, X2, Y2) and `fixed-witness-hesse-cofactor-d4-rational` (PR 17, d = 4
strict); same denominator family as `fermat-pair-rational-d4-strict`
(D = X0^4 - X0^2 X1^2 + X1^4 = N(X0^2 - ω X1^2)) and `coupled-pair-rational-d4-strict`.
The three-block family is the subject of docs/02a-binary-blocks.md (Chapter 2a:
all-target theorem for P0(X0,X1) + P1(Y0,Y1) + P2(X2^2, Y2^2); note the present F
is of the form P0(X0,Y0)+P1(X1,Y1)+P2(X2,Y2) with blocks paired as (X_i,Y_i)).
No open instance shares this variety or is a deformation of it;
`sameVarietyInstances` is empty.

---

## 2. tie-offblock-ytype-d4-strict-rational

**Data.** q = 1. F = X0^6 + X0^2 X1^2 Y1^2 + X1^6 + X2^6 + Y0^6 + Y1^6 + Y2^6
(Fermat plus one ternary tie). A = X0^3 Y0^2 Y1 (Y0^2 lies outside the tied
block {X0, X1, Y1}). Discussion #76. PR #14 (accepted 2026-10-15T19:20Z).
`sameVarietyInstances: [tie-singleton-3odd-strict]`.

**Method.** "A Hesse cofactor with a passive multiplier". Explicit, constructive,
strict.

**Witness.** d = 4,
D = 18 X0^4 + 6 X0^2 X1^2 - 18 X0^2 Y1^2 + 2 X1^4 + 6 X1^2 Y1^2 + 18 Y1^4
  = 18 (u^2 + v^2 + w^2 - uv - uw - vw),  u = X0^2, v = Y1^2, w = -X1^2/3
  = 18 (u + ω v + ω^2 w)(u + ω^2 v + ω w) **[computed]** (content 2 over Z;
irreducible over Q, splits over Q(√-3)).
N = (3 Y0^2 Y1 (X0^2 - Y1^2), 0, 0, 0, X0 Y0^2 (3 X0^2 + X1^2), 0), Q = 0.
div N = 6 X0 Y0^2 Y1 = ∂N0/∂X0 (N4 has no Y1); N(D) = D div N **[computed]**.

**Mechanism [proven, derivation].** With u, v, w as above the tie is
X0^2 X1^2 Y1^2 = -3uvw, so the (u,v)-dependence of F is u^3 + v^3 - 3uvw + const,
whose u- and v-partials coincide with those of the Hesse cubic
H = u^3 + v^3 + w^3 - 3uvw (the w^3 = -X1^6/27 term is not the X1^6 of F, but it is
never differentiated: only X0 and Y1 are active) **[computed]**. With
F_u = 3(u^2 - vw), F_v = 3(v^2 - uw) and N0 = 3 m Y1 (u - v), N4 = 3 m X0 (u - w),
m = Y0^2:
N(F) = 2X0 N0 F_u + 2Y1 N4 F_v = 18 m X0 Y1 [(u-v)(u^2-vw) + (u-w)(v^2-uw)]
     = 18 m X0 Y1 · u · (u^2+v^2+w^2-uv-uw-vw) = A D,
using the identity (u-v)(u^2-vw) + (u-w)(v^2-uw) = u (u^2+v^2+w^2-uv-uw-vw).
The quartic D is the cofactor of (u+v+w) in H ("Hesse cofactor"). The passive
multiplier m enters neither active derivative, so div N and N(D) are unchanged
by it: this is the precise sense of "passive".

**Observation [computed].** The same field with m replaced by any Y-even
quadratic monomial in the inert variables serves X0^3 Y1 m with the same
degree-4 D: m ∈ {Y0^2, Y0 Y2, Y2^2, X2^2} all give strict certificates
(sympy, E1 = E2 = 0, gcd 1). In particular m = Y0 Y2 gives a **strict degree-4
certificate for tie-singleton-3odd-strict** (A = X0^3 Y0 Y1 Y2), lower than the
accepted degree 6 there; see §3. (m = X2 Y0 also works algebraically but the
target X0^3 X2 Y0 Y1 is Y-even, outside the problem.)

**Symmetry [computed].** Torsion of Z^6/L(F) = Z/2 ⊕ (Z/6)^4, order 2592
(explicit description: (γ_i) ∈ μ_6^6 with γ0^2 γ1^2 γ4^2 = 1, modulo diagonal
μ_6). All 63 sign changes fix F. 35 nontrivial permutations fix F (S_3 on
{X2, Y0, Y2} times S_2 on {X0, X1}... as polynomial symmetries also Y1 ↔ X0/X1
inside the tie: the script lists e.g. X0 X1 X2 Y2 Y1 Y0 and X0 Y1 X2 Y0 X1 Y2).
A is odd in X0 and Y1 individually.

**Lean architecture.** Template of §0: Solution.lean 3118 B, Witness.lean
18298 B; t = X2 (absent from D), fFiber = X^6 + 6, dFiber = 32,
bezoutD = C (1/32). Verify degree 4, numerators 5.

**Related instances.** Same variety: `tie-singleton-3odd-strict` (solved, §3).
Same tie shape with the tie in a different triple: `f3-block-x0-quartic-strict`
(tie X0^2 X1^2 Y0^2, PR 10, d = 4, D = X0^4 - X1^2 Y0^2), the two-tie instances
`twoties-tie2-x2-quartic-strict`, `twoties-crossblock-hesse-cofactor-d4-rational`
(same Hesse-cofactor quartic in X2, Y1, Y2), `passive-quadratic-d4-strict-rational`,
and `family-passive-contact` (why.md's "mixed passive-quadratic example", solved,
PR 35, d = 6). No open instance is on this variety.

---

## 3. tie-singleton-3odd-strict

**Data.** q = 1. Same F as §2 (single tie X0^2 X1^2 Y1^2).
A = X0^3 Y0 Y1 Y2 (three odd Y coordinates; only Y1 is in the tied block).
Discussion #75. PR #13 (accepted 2026-10-15T19:10Z).
`sameVarietyInstances: [tie-offblock-ytype-d4-strict-rational]`.

**Method.** "Combine an active correction with passive radial motion".
Explicit, constructive, strict.

**Witness.** d = 6,
D = (7/8) X1^6 + (27/32)(X2^6 + Y0^6 + Y2^6)   (= (28 X1^6 + 27 X2^6 + 27 Y0^6 + 27 Y2^6)/32;
no factorization over Q **[computed]**; a weighted diagonal sextic in the four
coordinates {X1, X2, Y0, Y2}, which are exactly the coordinates not
differentiated by the two "active" components).
N0 = -(3/64) X1^2 Y0 Y1^3 Y2,  N4 = (1/64) X0 X1^4 Y0 Y2  (active corrections),
(N1, N2, N3, N5) = (9/64) A · (X1, X2, Y0, Y2)  (radial part),  Q = 0.
Degrees D 6, N_i 7 (all six nonzero). div N = (27/32) A, N(D) = (27/32) A D **[computed]**.

**Mechanism [proven, derivation].** Radial part R = c A (Z_i)_{i∈S},
S = {X1, X2, Y0, Y2}, c = 9/64. Because D is a sum of sixth powers in exactly the
coordinates of S, Euler gives R(D) = 6c A D regardless of the weights; and
div R = c A Σ_{i∈S} (1 + deg_{Z_i} A) = c A (1+1+2+2) = 6c A. The corrections N0, N4
have no X0 resp. Y1, so they contribute nothing to div N, and D has no X0, Y1, so
they contribute nothing to N(D). Hence E2 holds for any weights. The weights of D
are then fixed by E1: the radial part contributes 6c A Σ_{i∈S} (coefficient of
Z_i^6 in F) Z_i^6 = (27/32) A (X1^6 + X2^6 + Y0^6 + Y2^6) plus (2c) A X0^2 X1^2 Y1^2
from the tie; N0 F_{X0} and N4 F_{Y1} cancel the two tie-derived monomials
X0^5 X1^2 Y0 Y1^3 Y2 and X0 X1^4 Y0 Y1^5 Y2 against each other and against the
radial tie term, and N4 · 2 X0^2 X1^2 Y1 adds (1/32) A X1^6, which is why the X1^6
weight is 27/32 + 1/32 = 7/8 **[checked by hand and by sympy]**.
"Passive multiplier" reading: the pair Y0 Y2 enters only through the radial part.

**Observation [computed].** A strict degree-4 certificate also exists:
D' = 18(u^2+v^2+w^2-uv-uw-vw) (the §2 Hesse cofactor) with
N' = (3 Y0 Y2 Y1 (X0^2 - Y1^2), 0, 0, 0, X0 Y0 Y2 (3 X0^2 + X1^2), 0) (sympy: E1 = E2 = 0,
gcd(F, D') = 1). So the accepted degree 6 is not minimal, and the "passive factor
substitution" Y0^2 → Y0 Y2 asked about in why.md works verbatim for the Hesse
field. This is a new observation relative to the solution record, not a claim
about minimality (degree 2 has not been excluded; I did not search).

**Symmetry.** As in §2 (same F). A is odd in X0, Y0, Y1, Y2 individually.

**Lean architecture.** Template of §0: Solution.lean 2898 B, Witness.lean
18891 B; t = X0, fFiber = X^6 + X^2 + 5, dFiber = 109/32, bezoutD = C (32/109).
Verify degree 6, numerators 7.

**Related instances.** Same variety: §2. Radial-sextic denominator family:
`xtype-strict-rational-113` (D = block sextic + X1^6/3 + 27 X2^6/17 + 27 Y1^6/17, d = 6),
`xtype-sparse-rational-1m142` (d = 6), `family-passive-contact` (d = 6). Also
`one-reflection-y2-sector` (A = X0^3 Y0 Y1 Y2 on the one-reflection sextic;
solved by DiamondJM, PR 144, 59 files) is the "three odd Y coordinates" target on
a coupled (non-square) sextic. No open instance on this variety.

---

## 4. twoblock-cubic

**Data.** q = 1.
F = [X0^6 - (11/8) X0^4 Y0^2 + Y0^6] + X1^6 + Y1^6 + [X2^6 + (7/3) X2^4 Y2^2 - 2 X2^2 Y2^4 + Y2^6]
(two binary blocks in squares, both with nonzero mixed coefficients, plus two
spectator pure powers). A = X0^3 Y0^2 Y1 - 3 X0 Y0^4 Y1 = X0 Y0^2 Y1 (X0^2 - 3 Y0^2).
Discussion #88. PR #25 (accepted 2026-10-15T22:43Z). `sameVarietyInstances: []`.

**Method.** "A fixed denominator serving a companion bundle": choose D from the
binary-block theory, then solve the fixed-denominator linear system for N on
degree-7 character supports (solution.md). Explicit, constructive; the solved
completion turned out strict (C = 0), although the method allowed C ≠ 0.

**Witness.** d = 6,
D = Y0^2 (X0^2 - 3 Y0^2)^2 + 8 X1^6 + Y1^6/17 = r^2 + 8 X1^6 + Y1^6/17 with
r = Y0 (X0^2 - 3 Y0^2) the cubic factor of A = X0 Y0 Y1 · r **[computed]**
(irreducible over Q as a polynomial, factored content 1/17).
N0 = -(11/68) X0^4 Y0^2 Y1 + (81/68) X0^2 Y0^4 Y1 + (2/17) X1^6 Y1 - (36/17) Y0^6 Y1 + (2/17) Y1^7,
N1 = -(1/34) X0 X1 Y1 (X0^2 - 8 Y0^2)(4 X0^2 - 17 Y0^2),
N2 = 0,
N3 = -(3/34) X0 Y0 Y1 (X0^2 - 3 Y0^2)(4 X0^2 - 17 Y0^2),
N4 = -(1/34) X0 Y1^2 (2 X0^2 - Y0^2)^2,
N5 = 0, Q = 0. Degrees D 6, N_i 7.
div N = -X0 Y1 (12 X0^4 - 148 X0^2 Y0^2 + 411 Y0^4)/17 and N(D) = D div N **[computed]**.

**Mechanism [proven by the record; structure].** Denominator shape from
docs/02a-binary-blocks.md §2a.1: for a block P_i = C_i^3 + r_i^2 the field
e_i = r_i {C_i, ·}/(2 J_i) has e_i(r_i^2) = r_i^2, and D = r_i^2 + (other block)
gives e(D) = D; here the role of r_i^2 is played by r^2 = Y0^2 (X0^2 - 3Y0^2)^2 and
the "other block" by the spectators 8 X1^6 + Y1^6/17 (not by F's second block). The
first block is not literally of the form C^3 + r^2 with this r (I checked:
P0 - r^2 is not a cube), so the record's description "linear completion with D
fixed" is the honest mechanism: the numerator was found by linear algebra, and
the two identities were then verified exactly. The useful structural fact is that
N involves only X0, X1, Y0, Y1 and D involves only X0, X1, Y0, Y1.

**Observation [proven by inspection, also computed].** Since N2 = N5 = 0, D is
free of X2, Y2 and Q = 0, the identity E1 only uses ∂F/∂{X0, X1, Y0, Y1} and E2 does
not involve F at all. Hence the same (D, N) is a strict certificate for every
F' = P0(X0, Y0) + X1^6 + Y1^6 + g(X2, Y2) with g arbitrary (sympy check with
g = X2^6 + 3 X2^5 Y2 + Y2^6: E1 = E2 = 0, gcd 1). The second block is passive in the
strongest possible sense; the only F'-dependent obligations are smoothness of F'
and gcd(F', D) = 1 (which the all-ones specialization in t = X2 handles whenever g
is monic in X2 with g(X2 = t, Y2 = 1) coprime to the constant dFiber, i.e. always).
This answers why.md's question ("whether the second binary polynomial is
genuinely passive") affirmatively for this field.

**Symmetry [computed].** Torsion of Z^6/L(F) = (Z/2)^2 ⊕ (Z/6)^3, order 864.
All 63 sign changes fix F. Exactly one nontrivial permutation fixes F: X1 ↔ Y1.
A is odd in X0 and Y1 individually (and even in Y0).

**Lean architecture.** Template of §0: Solution.lean 2700 B, Witness.lean
21322 B (slightly longer N definitions); t = X2 (absent from D), fFiber =
X^6 + (7/3) X^4 - 2 X^2 + 29/8, dFiber = 205/17, bezoutD = C (17/205).
Verify degree 6, numerators 7.

**Related instances.** Target shape X0 Y0^2 Y1 (X0^2 - c Y0^2) on a single
deformed block: `ytype-companion-fiberA-d4-strict-rational` (block X0^6 - 2 X0^4 Y0^2 + Y0^6,
A = X0^3 Y0^2 Y1 - (6/7) X0 Y0^4 Y1, d = 4 quartic D in X0, Y0, Y1) and
`ytype-companion-231-d4-strict-rational` (block X0^6 + 2 X0^4 Y0^2 + 3 X0^2 Y0^4 + Y0^6, d = 4);
same block family at q = 2: `block-q2` (block X0^6 + 2 X0^4 Y0^2 - X0^2 Y0^4 + Y0^6,
D = that block + Y0^6 terms, d = 6, strict). The X-type analogues:
`xtype-strict-rational-113`, `xtype-sparse-rational-1m142`. All solved. No open
instance shares this variety or is a listed deformation.

---

## 5. two-reflection-cubic-surface-d2

**Data.** q = 1.
F = X0^6 + 5 X0^4 Y0^2 + 5 X0^2 X1^2 Y1^2 + 5 X0^2 X1 X2 Y1 Y2 + 5 X0^2 X2^2 Y2^2 + 5 X0^2 Y0^4
  + X1^6 + 5 X1^2 Y0^2 Y1^2 + 5 X1 X2 Y0^2 Y1 Y2 + X2^6 + 5 X2^2 Y0^2 Y2^2 + Y0^6 + Y1^6 + Y2^6
  = u^3 + w^3 + 5 (u + w)(u w + B) + X1^6 + X2^6 + Y1^6 + Y2^6,
  u = X0^2, w = Y0^2, B = X1^2 Y1^2 + X1 X2 Y1 Y2 + X2^2 Y2^2 **[computed]**.
A = X0^3 X1^2 Y0 - X0 X1^2 Y0^3 = X1^2 X0 Y0 (u - w). Discussion #80. PR #18
(accepted 2026-10-15T20:00Z). `sameVarietyInstances: []`.

**Method.** "A rotation field with a radial denominator". Explicit,
constructive, strict; the smallest denominator in this group.

**Witness.** d = 2, D = 4 (X0^2 + Y0^2), N = (-X1^2 Y0, 0, 0, X0 X1^2, 0, 0), Q = 0.
Degrees D 2, N 3. div N = 0 and N(D) = 0 **[computed]**.

**Mechanism [proven, derivation].** R = -Y0 ∂_{X0} + X0 ∂_{Y0} is the infinitesimal
rotation in the (X0, Y0) plane: div R = 0 and R(u + w) = 0, so D = 4(u + w) is a
first integral and N = X1^2 R is divergence-free with N(D) = 0. On F,
R(F) = 2 X0 Y0 (F_w - F_u) and for F = u^3 + w^3 + (u+w) Φ one has
F_u - F_w = (u + w)[3(u - w) + Φ_u - Φ_w]; here Φ = 5(uw + B) gives
Φ_u - Φ_w = -5(u - w), so F_u - F_w = -2 (u - w)(u + w) **[computed]** and
N(F) = 4 X0 Y0 X1^2 (u - w)(u + w) = A D. The construction uses only the two
reflections X0 → -X0, Y0 → -Y0 (so that u, w are coordinates) and the rotation
identity; it does not need the other four sign symmetries, which F does not have.
The passive square X1^2 can be replaced by any Y-even quadratic M in
X1, X2, Y1, Y2 (or by u + w), so the rotation trick serves at most the
7-dimensional subspace M · X0 Y0 (u - w) of the 32-dimensional two-reflection
sector mentioned in why.md; the rest of that sector is covered only by the
nonconstructive two-reflection extension of docs/02-explicit-constructions.md §2.5
(degree-45 denominators over k(α, β) before clearing) **[observation; the
dimension 7 is my count, not checked against the sector basis]**.

**Symmetry [computed].** Torsion of Z^6/L(F) = (Z/2)^2 ⊕ (Z/6)^2, order 144.
Sign changes fixing F: the 32-element group {s : s_{X1} s_{X2} s_{Y1} s_{Y2} = 1};
F is even in X0 and in Y0 individually, and only jointly in the other four.
15 nontrivial permutations fix F (X0 ↔ Y0 combined with the pair permutations of
{X1, X2}, {Y1, Y2}, and the swap (X1 X2)(Y1 Y2)-type symmetries listed by the script).
A is odd in X0 and Y0 individually (the two-reflection sector).

**Lean architecture.** Template of §0: Solution.lean 3008 B, Witness.lean
18526 B; t = X1 (absent from D; F is monic in X1), fFiber = X^6 + 10 X^2 + 10 X + 25,
dFiber = 8, bezoutD = C (1/8). Verify degree 2, numerators 3.

**Related instances (the one genuine open relative in this group).**
`quadric-contact-broken` (solved, HeadOfNaps, PR 122, 56 Lean files, strict
"finite-algebra contact certificate" over a cubic-surface model) and
`quadric-contact-residual` (**open**, discussion #125) share one sextic
F_b = F + (5/7) X0 X1 X2^2 Y0 Y1 + (10/7) X0 X1 Y0 Y1 Y2^2 **[computed term diff]**,
a deformation that destroys both individual reflections (the added terms are odd
in X0 and in Y0). `quadric-contact-broken` keeps the same A = X1^2 X0 Y0 (u - w);
`quadric-contact-residual` asks for A = X0^3 Y0^3 on F_b. The rotation certificate
does not survive the deformation: on F_b, E2 still holds but
N(F_b) - A D = (5/7) X1^3 Y1 (X0 - Y0)(X0 + Y0)(X2^2 + 2 Y2^2), a nonzero polynomial of
degree 8 not divisible by D **[computed]**; since a cleared defect must be
divisible by F_b^2 (degree 12), no correction with this (D, N) exists, so the
broken instance genuinely needs a different field (as its 56-file proof shows).
For the residual target X0^3 Y0^3, which is the reflection-invariant sector,
§2.5's endpoint construction is documented as not reaching it even on the
unbroken F (docs/02-explicit-constructions.md line 478: "The endpoint cubic
construction of §2.5 does not reach the reflection-invariant target X0^3 Y0^3").
Also thematically related (small quadric denominators, binary-in-squares
blocks): `ternary-tangent-residual-h2-d2`, `ternary-collision-cb-d2`,
`ternary-collision-cbp-d2`, `collision-h3` (all d = 2, strict, benediamond).

---

## 6. two-ternary-dense-two-odd-hessian-d48

**Data.** q = 1. F = P(X0, X1, Y0) + Q_b(X2, Y1, Y2), each block a dense cubic in
the squares (20 terms, coefficients 14, 7, 21, 23, 12, 23, 4, 39, 73 | 31, 37, 41, 43,
47, 53 | 8 | 59, 61, 67, 71):
P = 14 u^3 + b u^2 + c u + d (u = X0^2) with b = 7 X1^2 + 21 Y0^2,
c = 23 X1^4 + 12 X1^2 Y0^2 + 23 Y0^4, d = 4 X1^6 + 39 X1^4 Y0^2 + 73 X1^2 Y0^4 + 8 Y0^6;
Q_b = 59 w^3 + β w^2 + γ w + δ0 (w = Y1^2) with β = 43 X2^2 + 61 Y2^2,
γ = 37 X2^4 + 47 X2^2 Y2^2 + 67 Y2^4, δ0 = 31 X2^6 + 41 X2^4 Y2^2 + 53 X2^2 Y2^4 + 71 Y2^6.
A = X0^5 Y1 (odd in the two reflected coordinates X0, Y1). Discussion #104.
PR #41 (accepted 2026-10-16T01:37Z). `sameVarietyInstances: []`.
proof.json has no `witnessTheorem`; `witnessChecked: False` (no witness.json).

**Method.** "A rational certificate from ternary cubic covariants": an explicit
witness written in factored form through binary-cubic covariants (Hessian H,
discriminant Δ, the companion covariant t with t^2 = 4 H^3 - 27 Δ f^2), proved
through a four-variable chart, a lifted/homogenized field, and a polynomial
bridge. Constructive (explicit D, N, C as closed formulas), **non-strict**
(C ≠ 0), degree 55.

**Witness (factored; definitions in PExplicitWitnessDefinitions.lean).**
With Δ(a,b,c,d) = b^2 c^2 - 4 a c^3 - 4 b^3 d - 27 a^2 d^2 + 18 a b c d,
H_P = (b^2 - 3ac) u^2 + (bc - 9ad) u + (c^2 - 3bd) (degree 8), H_Q likewise in w,
Δ_P, Δ_Q (degree 12), P' = 3·14 u^2 + 2 b u + c, Q_b' similarly (degree 4):
G = H_P^3 Δ_Q Q_b^2 - H_Q^3 Δ_P P^2 (degree 48),  D = 36 X0 Q_b G (degree 55),
h0 = c^2 - 3bd, T_Q = 3(2(β^2 - 3·59 γ) w + βγ - 9·59 δ0) Q_b - 2 H_Q Q_b',
K_P = (b c^2 - 6 b^2 d + 9·14 c d) u + 2 c^3 - 9 b c d + 27·14 d^2,
L_P = u (b u^2 + 2 c u + 3 d), U = h0 H_P H_Q T_Q Q_b, V = K_P H_Q^3 P^2 - 3 Δ_Q Q_b^2 H_P^2 L_P,
R_P = 6P - 2u P', R_Q = 6 Q_b - 2 w Q_b',
R_0 = P H_Q^3 (R_P K_P - 6 u^3 Δ_P) - 18 u Δ_Q Q_b^2 H_P^2 (u c + 3 d), R_1 = R_0 + 2 u Q_b' U;
N = (0, -6 X1 Y1 Q_b V, X2 Y1 P R_1, -6 Y0 Y1 Q_b V, P (Y1^2 R_0 - X0^2 R_Q U), Y2 Y1 P R_1),
C = -72 X0 Y1 Q_b G R_0 (degree 104).
Expanded sizes **[computed, Singular]**: D 8050 terms, G 4900, N1 = N3 8740,
N2 = N5 14947, N4 15574, C 81016 terms. E1 and E2 hold exactly with this C ≠ 0;
gcd(F, D) = 1 over Q (Singular, 56 s). The solution text's "thousands of terms"
is accurate.

**Mechanism [proven in Lean; summary].** Chart X0 = 1, X1 = x, Y0 = y, Y1 = λ,
X2 = λ r, Y2 = λ z gives F = p(x,y) + λ^6 s(r,z), A = λ, with p = P(1,x,y),
s = Q_b(r,1,z). The core is a four-variable rational field W with
div W - (4/(3p)) dp(W) - (2/(3s)) ds(W) = 1 (weights 4/3, 2/3 come from the two
blocks and the λ-Jacobian). W is built from cubic covariants of the two
"binary cubics in one square" p (in u with coefficients in x, y) and s, and the
cancellation rests on the classical identities t^2 = 4H^3 - 27 Δ f^2,
k H + h0 t = 3 Δ l, and their derivative versions with coefficient weights
0, 2, 4, 6 (`Cubic.t_square`, `Cubic.k_H`, `k_derivative_identity`,
`l_derivative_identity`, `t_derivative_minus/plus`, `mixed_bridge`, all proved by
`simp only [...]; ring`). The master cancellation `Cubic.radial_scalar_identity`
is proved by one large `ring` normalization followed by rewriting with those
identities. The chart field is rescaled by 1/(18 p s), lifted with density λ^2
and a normal component, homogenized, and the divergence error proportional to F
is kept, giving dF(N) = A D and D div N - dD(N) = F C in the fraction field;
injectivity of the polynomial-to-rational-function embeddings and D ≠ 0 turn
these into polynomial identities (`equations_of_rational_field`).

**Legality over C [proven].** Two independent legality proofs appear:
(i) rational level (`PExplicitWitnessLegality`, `PDescentFrozenIrreducible`): F is
irreducible over Q by Eisenstein at p = 31 on the fibre X2 = 1, X1 = Y0 = Y1 = Y2 = 0
(F ↦ 14 X0^6 + 31 **[computed]**), lifted through `finSuccEquiv` with unit
leading coefficient (`irreducible_of_unit_leadingCoeff_map`); then
`DenominatorLegal F D ↔ ¬ F ∣ D`, and F ∤ D by the line X1 = X2 = Y0 = Y1 = Y2 = 1:
lineD (degree 13, odd) = lineF · quotient + remainder with a nonzero remainder
of degree 5 < 6. (ii) complex level (`PDenseP0`): the §0 pattern with t = X0,
`fUni_lc = 14`, fFiber = 14 t^6 + 28 t^4 + 58 t^2 + 634, dFiber the explicit
degree-13 polynomial with ~30-digit coefficients, and explicit Bezout cofactors
bezoutF (degree 12) and bezoutD (degree 5) with ~70-digit rational coefficients
checked by the `polynomial` tactic. Only (ii) feeds the protected theorem; (i)
serves the rational `HodgeTestbed.Certificate` structure (which additionally
records Y-parity of every component). D ≠ 0 from the nonzero X0-coefficient of dFiber.

**Lean architecture.** Solution.lean (723 B): `import HodgeTarget`,
`import Proposal.PDenseP0`, `maxHeartbeats 8000000`; binds `HodgeSubmission.F/A`
to `DenseP0.F/A` by `MvPolynomial.funext; norm_num [.., DenseP0.cm, Frozen.F, Frozen.X0, ..] <;> ring`
then `simpa using DenseP0.solved`. 45 helper files, 3855 lines, ~217 KB, all
named `P*.lean` and imported as `Proposal.*`; layered as:
- `PCertificate` (interface over Q: `SexticPoly`, `reflectY`, `YEven/YOdd`,
  `Certificate` structure extending `HodgeTestbed.Equations`, `Certificate.scale`),
  `PFrozen` (frozen F, A over Q, homogeneity, parity, `P0 := HasCertificate F A`).
- `PConicFields`, `PConicResidueReduction`, `PConicEulerCorrection`,
  `PConicHomogenization`, `PConicReflectionCover`: abstract derivation calculus
  (`fieldAction`, `fieldDivergence`, bracket, residue correction, homogenized
  field with Euler term, cover field) over arbitrary commutative Q-algebras.
- `PFrozenChartPolynomial/Rational/ChartMap/ChartRules/Compatibility/
  SquareCover/NormalDerivative/FrozenIdentities`: `RatFun n = FractionRing (MvPolynomial (Fin n) ℚ)`,
  partial derivatives on fraction fields via dual numbers (`TrivSqZeroExt` jets),
  `Localization.Away X0` chart, injectivity lemmas, square cover, normal derivative.
- `PIncidenceBlocks`, `PEllipticAmbientAmbient/Split/FourVariable/BaseEmbedding/
  Physical/ScalarNumerator/Transfer/Cubic/RadialIdentity/CubicDerivative/
  ActualScalar/ActualField/FactoredField/CompactPhysical/HomogeneousMap/
  HomogeneousRules/CovariantScaling/PolynomialBridge/ActualEquations`:
  the split F = p + λ^6 s, the cubic covariant algebra (`Cubic.Data`,
  `t_square`, ...), the radial identity, the specific coefficient data
  (`pc`, `qc`, `Cubic.Weights`), the four-variable lift (`lifted_divergence_of_primitive`),
  physical Jacobian (`physical_jacobian`, Mathlib determinant), the 30 `hom_*`
  scaling rules (e.g. `homBase (Cubic.H pc) = embed HP / X0^8`,
  `homBase (Cubic.G pc qc) = embed G / (X0^24 Y1^24)`), the bridge
  `equations_of_rational_field`, and `actual_equations : Equations Frozen.F Frozen.A D N C`.
- `PExplicitWitnessDefinitions/Grading/Parity/Legality/Complete`: the factored
  witness, 40 homogeneity lemmas (degrees 0..104), 40 `reflectY` parity lemmas,
  legality over Q, `ExplicitWitness.certificate` and `frozenP0`.
- `PDescentFrozenIrreducible` (Eisenstein/Gauss lemma, `frozen_irreducible`),
  `PDenseP0` (transport to C via `MvPolynomial.map (algebraMap ℚ ℂ)`,
  `directional_cm`, `divergence_cm`, complex legality, `Verify F A 1 55 D N Q`, `solved`).
Heartbeats: only `PDenseP0` (4,000,000, maxRecDepth 1000) and
`PEllipticAmbientFourVariable` (800,000) raise the limit; all other helpers
run at the default. `PDenseP0` ends with `#print axioms`.
Tactics for identities: `ring`, `simp only` + `ring`, `linear_combination`
(PConic modules), `field_simp`, `polynomial`, `compute_degree!`, `norm_num`;
coprimality: `relPrime_of_retract` + `relPrime_of_specialization` + explicit
Bezout; irreducibility: `irreducible_of_eisenstein_criterion`.

**Symmetry [computed].** Torsion of Z^6/L(F) = (Z/2)^4 ⊕ Z/6, order 96 (dense
blocks leave only the sign group and the block-wise μ_6 scalings). All 63 sign
changes fix F (F is coordinatewise even), so §2.5's full cubic-surface theorem
gives nonconstructive existence for every Y-odd target at every q (as why.md
states); no variable permutation fixes F. A is odd in X0 and Y1 individually.

**Related instances.** No open instance shares the variety. why.md points to
`split-ternary-q2` (solved, fqlx, PR 159, 27 files, third-order pole) as the next
step with the square symmetries removed; `crossed-ternary-asymmetric`
(1337isnot1337, PR 149, A = X0^4 X1 Y1, elliptic transfer) and the small
ternary instances `ternary-tangent-residual-h2-d2`, `ternary-collision-*` are the
other "ternary block in squares" solves. Among open instances the only thematic
neighbour is `dense-rational` (full coefficient support, A = X0^3 Y0^3), which has
no coordinate evenness, so the covariant method of this proof does not apply to
it as stated.

---

## 7. Cross-cutting summary

| instance | q | d | D structure | Q | div v = 0? | t for coprimality | Lean size |
|---|---|---|---|---|---|---|---|
| three-block-rational-d4-strict | 1 | 4 | N_{Q(ω)/Q}(L - ω X2^2), L = X1^2 - 2Y1^2 | 0 | yes | X0 | 2 files, 22 KB |
| tie-offblock-ytype-d4-strict-rational | 1 | 4 | 18·N_{Q(ω)/Q}(u + ωv + ω^2 w), Hesse cofactor | 0 | yes | X2 | 2 files, 21 KB |
| tie-singleton-3odd-strict | 1 | 6 | weighted diagonal sextic in {X1,X2,Y0,Y2} | 0 | yes | X0 | 2 files, 22 KB |
| twoblock-cubic | 1 | 6 | r^2 + 8X1^6 + Y1^6/17, r = Y0(X0^2-3Y0^2) | 0 | yes | X2 | 2 files, 24 KB |
| two-reflection-cubic-surface-d2 | 1 | 2 | 4(X0^2 + Y0^2), rotation first integral | 0 | yes | X1 | 2 files, 21 KB |
| two-ternary-dense-two-odd-hessian-d48 | 1 | 55 | 36·X0·Q_b·(H_P^3 Δ_Q Q_b^2 - H_Q^3 Δ_P P^2) | ≠ 0 (deg 104) | no | X0 | 46 files, 217 KB |

Facts worth reusing:
1. Every coprimality proof in this group reduces to `relPrime_of_retract` +
   `relPrime_of_specialization` at the all-ones point with a unit leading
   coefficient in one variable; for the five small solves that variable is absent
   from D, so the fibre Bezout identity is trivial (`0·f + (1/c)·c = 1`).
2. The five small identities are discharged by `MvPolynomial.funext` + `norm_num`
   + `ring` (evaluation at an arbitrary complex point); no polynomial-arithmetic
   tactic beyond `ring` is used. The large proof keeps everything factored and
   proves the cancellation symbolically, then specializes coefficients.
3. Passive multipliers: in §2/§3 the Hesse field admits any Y-even quadratic
   multiplier in the inert variables (four Y-odd targets from one degree-4 D);
   in §4 the whole (X2, Y2) block is inert; in §5 the passive square X1^2 can be
   any Y-even quadratic in the non-reflected variables. These are checked
   computationally, not claimed in the records.
4. Failed/limited approaches recorded by the project for these varieties: the
   rotation field fails on the quadric-contact deformation (explicit remainder
   above); the §2.5 endpoint construction does not reach X0^3 Y0^3, which is the
   open `quadric-contact-residual` target. No failed searches are recorded in
   the six solution.md files themselves.

Scripts: `/tmp/hodge/scripts/verify_witness.py` (sympy, five witnesses),
`/tmp/hodge/scripts/hessian_d48.sing` (Singular, factored witness; 56 s),
plus the inline sympy checks reported above (factorizations over Q(√-3),
Smith normal forms, sign/permutation groups, the three observations).
