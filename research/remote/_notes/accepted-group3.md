# Accepted solves, group 3: family-trace-zero, fermat-pair-rational-d4-strict, fermat-q2, fixed-witness-hesse-cofactor-d4-rational, hesse-radial-hesse-cofactor-d4-rational, hesse-triangle

Catalogue of six accepted constructive proofs in `/home/kejace/devel/sextics/data/<id>/`.
Variables are ordered `X0, X1, X2, Y0, Y1, Y2` (Lean indices 0..5). Conventions of
`formal/README.md`: a certificate for `(F,q,A)` is `(d, D, N_0..N_5, Q)` with `D` nonzero,
homogeneous of degree `d`, `IsRelPrime F D`, `deg N_i = d + 6q - 5`, and

    A*D^2 - q*D*sum_i N_i*dF/dZ_i + F*(D*sum_i dN_i/dZ_i - sum_i N_i*dD/dZ_i) = F^(q+1)*Q.

Notation below: `N(P) := sum_i N_i dP/dZ_i` (directional derivative), `div N := sum_i dN_i/dZ_i`,
`E1 : q*N(F) = A*D`, `E2 : D*div N - N(D) = F^q*Q`.

Status legend used in this file:
- **[proven]** = checked by the accepted Lean proof (kernel replay, axioms `propext, Quot.sound, Classical.choice`, Lean 4.33.0, `docker-no-network` sandbox) and re-checked here in exact arithmetic;
- **[verified-sympy]** = exact rational/algebraic computation done in this session (sympy 1.12, scripts in `/tmp/hodge/scratch/`), not part of any accepted record;
- **[derived]** = hand derivation consistent with the verified computations but not itself machine-checked;
- **[conjecture]** / **[not checked]** = stated as such.

---

## 0. Common facts for all six records

### 0.1 Provenance
| id | PR | accepted (UTC) | discussion | comments |
|---|---|---|---|---|
| fermat-pair-rational-d4-strict | #1 | 2026-10-15 17:04 | #63 | 0 |
| hesse-radial-hesse-cofactor-d4-rational | #16 | 2026-10-15 19:39 | #78 | 0 |
| fixed-witness-hesse-cofactor-d4-rational | #17 | 2026-10-15 19:48 | #79 | 1 (benediamond, 2026-10-28, Codex-assisted; explicit Darboux field, real-section topology) |
| hesse-triangle | #23 | 2026-10-15 22:23 | #86 | 0 |
| family-trace-zero | #32 | 2026-10-15 23:59 | #95 | 0 |
| fermat-q2 | #37 | 2026-10-16 00:52 | #100 | 0 |

All six: `githubLogin: benediamond`, `proof.json` schemaVersion 2 with `theorem: Submission.certificate`,
`witnessTheorem: Submission.witness`; `witnessChecked: true`; `result: lean-checked-awaiting-owner`
then owner-merged. `provenance.json.sameVarietyInstances` is empty for all except the pair
fermat-pair ↔ fermat-q2 (same Fermat sextic).

### 0.2 Shape of the certificates
All six are **constructive, strict (`Q = 0`), rational** witnesses: `witness.json` has `Q: "0"`, every
coefficient is in ℚ, every `N_i` has the involution parity required by Chapter 1 (X-components Y-odd,
Y-components Y-even) **[verified-sympy]**. Summary table:

| id | q | d | D (sparse) | nonzero N_i | deg N | legality variable t |
|---|---|---|---|---|---|---|
| family-trace-zero | 1 | 6 | `2X1^6+3X2^6+5Y1^6+7Y2^6` | all six | 7 | X0 (lc 2) |
| fermat-pair-rational-d4-strict | 1 | 4 | `X0^4-X0^2X1^2+X1^4` | N0, N1 | 5 | X2 (lc 1) |
| fermat-q2 | 2 | 6 | `6X1^6+6X2^6` | N1, N2 | 13 | X0 (lc 1) |
| fixed-witness-hesse-cofactor-d4-rational | 1 | 4 | `X0^4+X0^2X1^2/3-2X0^2X2^2/3+X1^4/9-X1^2X2^2/9+X2^4/9` | N0, N1, N2 | 5 | Y0 (lc 1) |
| hesse-radial-hesse-cofactor-d4-rational | 1 | 4 | `X0^4-X0^2X2^2+X0^2Y2^2/3+X2^4+X2^2Y2^2/3+Y2^4/9` | N0, N2 | 5 | X1 (lc 5) |
| hesse-triangle | 1 | 6 | `X0^6+X1^6+X2^6+2Y2^6` | all six | 7 | Y0 (lc 1) |

In every case `D` omits a variable `t` in which `F` is (up to a nonzero constant) monic of degree 6; this is the
entire content of the Lean coprimality proof (0.3).

### 0.3 The shared Lean template (identical structure in all six; machine-generated, generator not in repo)

Two files per instance, no `HodgeComputer.Constructions.*` imports (these predate the shared library):

**`Witness.lean`** (helper module `Proposal.Witness`, namespace `HodgeComputer.Instances.<CamelCaseId>`), imports only
`HodgeComputer.Verify` (which imports `HodgeTestbed.Algebra` for `directional`, `divergence`, `divergenceNumerator`).
Options: `maxHeartbeats 2000000`, `maxRecDepth 1000`, linters `unusedTactic`, `unusedSimpArgs`, `unnecessarySeqFocus`,
`unreachableTactic` off. Declaration chain:
1. `def F A D Q N0 … N5 : Poly` as explicit `MvPolynomial.C c * (X i)^k` sums; `def N : Fin 6 → Poly` as an `if i.val = k` chain.
2. 36 private lemmas `derivative_i_j : pderiv i (X j) = (0|1)` by `norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]`.
3. `theorem e1 : (q : Poly) * directional F N = A*D` and `theorem e2 : divergenceNumerator D N = F^q*Q`: `simp only [directional/divergence…, Fin.sum_univ_succ, …]`, `change` to the explicit six-term sum, `simp only [F,A,D,N0..N5, map_zero, map_add, pderiv_C, pderiv_mul, pderiv_pow, derivative_*]`, then `apply MvPolynomial.funext; intro x; norm_num <;> ring`. So the polynomial identity is discharged as an identity of evaluation functions on ℂ⁶ (`MvPolynomial.funext`, valid since ℂ is infinite) closed by `ring`. No `linear_combination`, no `decide`.
4. `denominator_homogeneous : D.IsHomogeneous d`, `N{i}_homogeneous : N{i}.IsHomogeneous (d+6q-5)` by explicit `(isHomogeneous_C …).mul ((isHomogeneous_X ℂ i).pow k)` / `.add` chains; `numerator_homogeneous` by `fin_cases`.
5. Coprimality over ℂ: `toUni : Poly →+* Polynomial Coeff` (`Coeff := MvPolynomial (Fin 5) ℂ`, `MvPolynomial.eval₂Hom` sending the chosen variable `t` to `Polynomial.X` and the other five to `Polynomial.C (X k)`), `fromUni` by `Polynomial.eval₂RingHom`, `retract : fromUni.comp toUni = RingHom.id Poly` by `MvPolynomial.ringHom_ext` + `fin_cases` + `simp`; `fUni` = explicit univariate form of `F`, `toUni_F` by `simp only […]; polynomial` (Mathlib's `polynomial` tactic, `Mathlib/Tactic/Polynomial/Basic.lean`, an `algebra`-based normalizer), `fUni_degree : natDegree = 6` by `compute_degree!`, `fUni_lc : leadingCoeff = C c`; `point : Fin 5 → ℂ := (1,1,1,1,1)`, `specialize := MvPolynomial.eval point`; `fFiber, dFiber : Polynomial ℂ` with `fiber_F`, `fiber_D` by `norm_num […] <;> polynomial`; `fiber_coprime : IsCoprime fFiber dFiber := ⟨0, C (1/c'), by polynomial⟩` (Bezout; `dFiber` is the nonzero constant `c' = D(1,…,1)` because `D` omits `t`); finally
   `denominator_legal : IsRelPrime F D := relPrime_of_retract toUni fromUni retract (relPrime_of_specialization specialize (unit leading coeff) fiber_coprime.isRelPrime)`.
   Mathematics: a common divisor of `F` and `D` in `ℂ[others][t]` has unit leading coefficient (it divides a `t`-monic `F`) and unit image under the specialization, hence is a unit (`Polynomial.isUnit_of_isUnit_leadingCoeff_of_isUnit_map`); `relPrime_of_retract` transports back along the coordinate isomorphism. This proves coprimality over ℂ directly, not via ℚ-irreducibility. Reuse condition: `D` must omit some variable in which `F` contains a pure sixth power with nonzero coefficient.
6. `denominator_ne_zero : D ≠ 0` by evaluating at `(1,…,1)`; `certificate : Verify F A q d D N Q` with `identity := normalized_identity (by simpa using e1) e2`; `solved : HasCertificate F A q`; `#print axioms`.

Per-instance fibers used in step 5 **[read from source]**:
| id | t | `fFiber` = F(t at index t, 1 elsewhere) | `dFiber` |
|---|---|---|---|
| family-trace-zero | X0 | `2t^6+2t^4+(62022/10609)t^2+242375/10609` | 17 |
| fermat-pair | X2 | `t^6+5` | 1 |
| fermat-q2 | X0 | `t^6+5` | 12 |
| fixed-witness-hesse | Y0 | `t^6+27` | 7/9 |
| hesse-radial | X1 | `5t^6+21t^4+30t^2+19` | 16/9 |
| hesse-triangle | Y0 | `t^6-3t^2+7` | 5 |

**`Solution.lean`** (42 lines, 2.6–3.2 KB each): `import HodgeWitness` (verifier-generated protected module:
`HodgeSubmission.F/A/q` from `instance.rec`, `HodgeWitness.d/D/N0..N5/N/Q` from `witness.json`; see
`scripts/proposals/check.py::target_source`, `witness_source`) and `import Proposal.Witness`;
`maxHeartbeats 8000000`, `maxRecDepth 8000`.
`theorem witness : HodgeComputer.Verify HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q HodgeWitness.d HodgeWitness.D HodgeWitness.N HodgeWitness.Q`
proves `hF hA hD hQ` and `hN` (per component via `funext i; fin_cases i`) by
`apply MvPolynomial.funext; intro x; norm_num [both defs] <;> ring`, then
`simpa only [hF,hA,hD,hQ,hN,HodgeSubmission.q,HodgeWitness.d] using <Inst>.certificate`.
`theorem certificate : HasCertificate … := ⟨d, D, N, Q, witness⟩`.

File sizes (lines / bytes): family-trace-zero 42/2744 + 362/22407; fermat-pair 42/2986 + 335/17378;
fermat-q2 42/2590 + 333/17120; fixed-witness-hesse 42/3184 + 353/20689; hesse-radial 42/3162 + 347/19022;
hesse-triangle 42/2700 + 338/18677.

### 0.4 Two reusable mechanisms behind the six (Chapter 2 §2.1–2.2 language)

Cofactor (Piola) transport **[proven in docs, verified-sympy on each instance]**: for a polynomial map Φ with
`J = dΦ`, `Δ = det J`, and a target-space field `w`, `V := adj(J)·w(Φ)` satisfies `V(G(Φ)) = Δ·w(G)(Φ)` and
`div V = Δ·(div w)(Φ)`.

(M1) *Radial sextic denominator* (hesse-triangle, family-trace-zero): write `F = G_c(Φ) + D_t(Φ)`; find a polynomial
target field `P` with `P(G_c) = 0`, `P(D_t) = ρ·D_t`, `div P = ρ` (same scalar ρ). Then with `N = adj(J)·P(Φ)`,
`D = D_t(Φ)`, `A = Δ·ρ(Φ)`: `N(F) = A·D`, `div N = A`, `N(D) = A·D`, so `E1` and `E2 = 0` hold strictly.
(This is the §2.1 constructor with `(μ,ν) = (0,1)`: `w(D_0) = D_0`, `div(bw) = b`.)

(M2) *Quartic Hesse cofactor* (fixed-witness-hesse, hesse-radial, and an unrecorded alternative for hesse-triangle):
in quadratic coordinates `(u,v,w)` write `F̃ = S·Q + K` with `S = u+v+w`, `Q = u²+v²+w²-uv-uw-vw`
(`u³+v³+w³-3uvw = S·Q`, `Q = L_1 L_2` the two non-real Hesse lines `L_k = u+ω^k v+ω^{2k} w`). A *linear* field `E` on
`(u,v)` (fixing `w` and all passive coordinates) with `E(Q) = 3Q`, `div E = 3`, `E(K) = 0` makes `W = E/Q`
divergence-free with `W(F̃) = E(S) + 3S`. Piola transport with a passive density `b` gives `D = Q(Φ)` (up to a
constant), `N = adj(J)·(bE)(Φ)`, `A = Δ·b·(E(S)+3S)(Φ)`. The admissible `E` form a 2-dimensional space
(`E = (2u-v-w)∂_u + (u+v-2w)∂_v` with `E(S)+3S = 6u+3v`, and `E' = 3[(v-w)∂_u + (v-u)∂_v]` with
`E'(S)+3S = 9v`), so the served targets are `Δ·b·(αu+βv)` — no `w` term. (Darboux integrating-factor
interpretation, as in benediamond's comment on #79.) **[verified-sympy on the three instances]**

Subset Euler (§2.2) **[proven in docs; verified-sympy]**: `E_I = Σ_{i∈I} Z_i ∂_i`; if every monomial of `F`
has `I`-degree 0 or 6 (equivalently `E_I²(F) = 6E_I(F)`) and `deg_I A + |I| = 6` for a monomial `A`, then
`v = A·E_I/(q·E_I(F))` is strict for every `q`, with `D = E_I(F)`.

---

## 1. family-trace-zero — "Trace zero · line × conic" (PR #32)

**Target.** `q = 1`,
`F = 2X0^6 + 2X0^4Y0^2 + (10201/10609)X0^2X1^4 + (20402/10609)X0^2X1^2Y1^2 + 2X0^2Y0^4 + (10201/10609)X0^2Y1^4 + 2X1^6 + (10201/10609)X1^4Y0^2 + (20402/10609)X1^2Y0^2Y1^2 + 3X2^6 + 2Y0^6 + (10201/10609)Y0^2Y1^4 + 5Y1^6 + 7Y2^6`
(14 monomials; `10201/10609 = (101/103)²`).
`A = 8X0^3X1^2Y0 - 8X0X1^2Y0^3 = 8·X0Y0(X0²-Y0²)·X1²`.
Structure **[verified-sympy]**: with `U = X0²+Y0²`, `V = X0²-Y0²`, `Q = (101/103)(X1²+Y1²)`,
`F = U·(U²+V²+Q²) + D_t`, `D_t = 2X1^6+3X2^6+5Y1^6+7Y2^6`.
The "cubic complement" `U(U²+V²+Q²)` is a line times a smooth conic in the `(U,V,Q)` plane (deformation of the
Hesse triangle of §6).

**Method.** Trace-zero rotation preserving the conic, cofactor/Piola transport through
`Φ = (U, V, X1, X2, Y1, Y2)`, `Δ = det dΦ = -8X0Y0` — mechanism (M1).

**Certificate (constructive, strict).** `d = 6`, `D = D_t = 2X1^6+3X2^6+5Y1^6+7Y2^6` (diagonal sextic in four
variables: irreducible over ℂ, a cone over a smooth diagonal sextic surface with vertex line
`X1=X2=Y1=Y2=0`); `Q = 0`; all six numerators of degree 7:
`N1,N2,N4,N5 = (A/6)·(X1, X2, Y1, Y2)` = `(4/3)X0X1²Y0(X0²-Y0²)·(X1,X2,Y1,Y2)`,
`N0 = -(2/3)·Q²·B·Y0`, `N3 = +(2/3)·Q²·B·X0` with `B = X1²` (so `N0 = -(20402/31827)X1²Y0(X1²+Y1²)²`,
`N3 = (20402/31827)X0X1²(X1²+Y1²)²`).

**Mechanism (5–10 lines) [verified-sympy: Piola transport of the stated flux reproduces the witness exactly].**
Target-space flux `P = B·[(Q²/3)∂_V - (V/6)E_sp]`, `E_sp = X1∂X1+X2∂X2+Y1∂Y1+Y2∂Y2` (Euler on the four
spectators), `B = X1²`. Since `E_sp Q = 2Q`, `P(V²+Q²) = B[(2/3)VQ² - (V/6)·4Q²] = 0` and `U` is untouched, so
`P(U(U²+V²+Q²)) = 0`. `E_sp D_t = 6D_t` gives `P(D_t) = -VB·D_t`; `div P = -VB(3/6 + 1/6+1/6+1/6) = -VB`.
Hence ρ = `-VB` in (M1): `N(F) = Δρ D_t = A·D`, `div N = Δρ = A`, `N(D) = A·D`.
The rotation `Q∂_V - V∂_Q` of the conic is realized as `(Q/3)[Q∂_V - V·E_sp/2]` because `E_sp/2` acts on
functions of `Q` as `Q∂_Q`; the integrating factor `Q/3` and the density `B` make `div P = ρ`.
The target necessarily contains the Jacobian density `Δ = -8X0Y0` and the eigenvalue form `V`; `B` is a free
spectator density with `E_sp B = 2B` **[derived]**: any Y-even quadratic in `(X1,X2,Y1,Y2)` works, so this
construction serves the 6-dimensional space `Δ·V·B`, nothing else. (Chapter 2 §2.1's literal constructor would use
`D_0 = w(G) = -V·D_t` of degree 8; the record divides out `V` because `D_t` itself is an eigenfunction.)

**Symmetry.** Diagonal symmetry group of `F` modulo scalars has order 288 (Smith form of the exponent-difference
lattice `[2,2,2,6,6]`) **[verified-sympy]**; `A` is a two-term eigenvector (its monomials differ by
`2(e_{Y0}-e_{X0}) ∈ L(F)`).

**Lean.** Template of 0.3 with `t = X0`, leading coefficient 2 (`(2:ℂ) ≠ 0`), `dFiber = 17`.
`Witness.lean` 362 lines / 22.4 KB (the largest of the six: 14-term `F`, six nonzero `N_i`).

**Relations.** Solved: hesse-triangle (the degenerate triangle it deforms; same (M1) mechanism),
family-hesse-web (#34), family-sylvester-norm (#30; same weighted-sixth-power spectator style). Open: none on the
same variety; no open instance has a line × conic complement. (`sameVarietyInstances: []`.)

---

## 2. fermat-pair-rational-d4-strict — "Fermat · quartic cofactor" (PR #1)

**Target.** `q = 1`, `F = X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6`, `A = X0^3·X1·X2·Y0`.
Diagonal symmetry group `μ_6^6/μ_6` (order 6⁵ = 7776) **[verified-sympy]**, times coordinate permutations;
character of `A` = exponent vector `(3,1,1,1,0,0) mod 6`.

**Method.** Quartic cofactor of the binary sixth-power sum, passive multiplier `B = X2Y0`:
`x^6+y^6 = (x²+y²)·D`, `x = X0`, `y = X1`.

**Certificate (constructive, strict).** `d = 4`, `D = X0^4 - X0^2X1^2 + X1^4 = X1^4·Φ_12(X0/X1)`;
over ℚ(ω): `(X0²+ωX1²)(X0²+ω²X1²)`; over ℚ(ζ₁₂): four hyperplanes `X0 = ζX1`, ζ primitive 12th root of unity
**[verified-sympy]**. `Q = 0`; `N0 = (B/6)·X1(X0²-X1²)`, `N1 = (B/6)·X0^3`, `N2..N5 = 0`; degree 5.

**Mechanism [proven; decomposition verified-sympy].**
`N(F) = B(x^5·y(x²-y²) + x^3·y^5) = x^3 y B·D = A·D`; `div N = xyB/3`; `N(D) = (xyB/3)·D`.
Decomposition: `v = N/D = v_E + (B/6)·H/(x^6+y^6)` where `v_E = A(x∂_x+y∂_y)/(6(x^6+y^6))` is the §2.2
subset-Euler field (`I = {X0,X1}`, `deg_I A + |I| = 4+2 = 6`) and `H = -y^5∂_x + x^5∂_y` is the Hamiltonian
field of the block (`H(x^6+y^6) = 0`, `div(H/(x^6+y^6)) = 0`). Adding `c·B·H/(x^6+y^6)` preserves the
certificate for every constant `c`; for `c = 1` the numerator acquires the factor `x²+y²`, which cancels against
`x^6+y^6` and lowers the denominator from 6 to 4.
Denominator family **[verified-sympy over ℚ(ω)]**: cancellation requires `c³ = 1` **[derived]**, and the three
values give the three quartic cofactors
`D_c = (x^6+y^6)/(y²+c·x²)`: `D_1 = x^4-x²y²+y^4` (rational, accepted), `D_ω = ω²x^4 - ωx²y² + y^4`,
`D_{ω²} = conj`, each a strict certificate for the same `A`; they are the orbit of the rational one under
`y ↦ ηy`, `η ∈ μ_3 ⊂ μ_6` (Chapter 1 character transport). Only `c = 1` is rational, hence "rational-d4".

**Lean.** Template with `t = X2` (D lives in `ℚ[X0,X1]`), `fFiber = t^6+5`, `dFiber = 1`.
`Witness.lean` 335 lines / 17.4 KB.

**Relations.** Same variety: fermat-q2. Solved deformations of this exact construction:
coupled-pair-rational-d4-strict (#15: `F = Fermat + X0^4Y1^2 - X0^2Y0^2Y1^2 + Y0^4Y1^2`,
`A = X0^3X1X2Y0 + X0X1X2Y0Y1^2/3`), three-block-rational-d4-strict (#7, rational quartic from a binary block).
Open: the Fermat sextic is the base point of every open instance (all sixteen are `Fermat + couplings`, or
perturbed-diagonal variants), see §7.

---

## 3. fermat-q2 — "Fermat · pole order two" (PR #37)

**Target.** `q = 2`, `F` Fermat, `A = X0^3X1^2X2^2Y0^3Y1^2` (degree 12); character `(3,2,2,3,2,0) mod 6`.
Same target polynomial as holdout-q2 (solved, HeadOfNaps #124) on the holdout sextic.

**Method.** Subset Euler field (§2.2) at pole order two, `I = {X1, X2}`: `deg_I A + |I| = 4 + 2 = 6`.

**Certificate (constructive, strict).** `d = 6`, `D = E_I(F) = 6(X1^6+X2^6) = 6(X1²+X2²)(X1^4-X1²X2²+X2^4)
= 6·∏_{ζ^6=-1}(X1-ζX2)` (six hyperplanes over ℚ(ζ₁₂)); `Q = 0`; `N = (A/2)·E_I`: `N1 = A·X1/2`, `N2 = A·X2/2`,
degree 13 = 6 + 12 − 5.

**Mechanism [proven].** `2N(F) = A·E_I(F) = A·D`; `div N = (deg_I A + |I|)·A/2 = 3A`; `N(D) = (A/2)E_I(D) = 3AD`
since `D` is `I`-homogeneous of degree 6. So `D·div N - N(D) = 0` identically and the residual is divisible by `F^3`
trivially. The factor `1/q` is the only `q`-dependence; the same tuple shape works at every `q`.
Legality: `D` omits `X0`, in which `F` is monic.

**Denominator family [verified-sympy, new here].** The accepted field is the ℚ(i)-trace of a certificate with a
*cubic* denominator: `D' = X1^3 + i·X2^3` (a factor of `X1^6+X2^6` over ℚ(i)),
`N'1 = B'·X2²/12`, `N'2 = i·B'·X1²/12`, `B' = A/(X1²X2²) = X0^3Y0^3Y1^2` (degree 10 = 3+12−5): `E1`, `E2` vanish
and `gcd(F, D') = 1`; averaging with the conjugate `X1^3 - iX2^3` returns exactly `(A/12)(X1,X2)/(X1^6+X2^6)`,
the accepted subset-Euler field. Within this one-parameter Hamiltonian gauge the only rational member is the
sextic **[derived]**; whether a rational denominator of degree < 6 exists for this target was not investigated
**[open question]**.

**Lean.** Template with `e1 : (2 : Poly) * directional F N = A*D`, `e2 : … = F^2*Q`, `t = X0`,
`fFiber = t^6+5`, `dFiber = 12`. `Witness.lean` 333 lines / 17.1 KB.

**Relations.** Same variety: fermat-pair. why.md names holdout-q2 (solved) and block-q2 (solved, #38) as the
coupled/intermediate comparisons. Open: **cycle-q2** (`q = 2`, cyclic sextic of cycle-three,
`A = X0^2X1^2X2^3Y0^2Y1^2Y2`) is the open pole-order-two instance for which fermat-q2 is the diagonal control;
the subset-Euler criterion fails there for every subset `I` **[verified-sympy scan]** because the couplings
`X0^3X1Y0Y1/2, X0X2^3Y0Y2/2, X1^3X2Y1Y2/2` have `I`-degrees strictly between 0 and 6 for every proper subset meeting them.

---

## 4. fixed-witness-hesse-cofactor-d4-rational — "Hesse · fixed-target quartic" (PR #17)

**Target.** `q = 1`,
`F = X0^6 + 12X0^4X1^2 + 9X0^4X2^2 + 6X0^2X1^4 - 3X0^2X1^2X2^2 - 6X0^2X2^4 + 2X1^6 + X2^6 + Y0^6 + 2Y1^6 + 3Y2^6`
(11 monomials), `A = X0·X1·X2^3·Y1`.
Hidden structure **[verified-sympy]**: `u = X0²`, `v = X2²-2X0²`, `w = X1²+X0²`, `S = u+v+w`,
`Q = u²+v²+w²-uv-uw-vw`: `F = S·Q + K`, `K = w^3 + Y0^6 + 2Y1^6 + 3Y2^6 = (X0²+X1²)^3 + G(Y)`.
Equivalently `F = (u³+v³+w³-3uvw) + w³ + G(Y)`: the singular Hesse member plus a cube of the third coordinate,
in *mixed* quadratic coordinates (the same shape as hesse-triangle's `C + 2W³`).

**Method.** Hidden Hesse cofactor after a nonlinear quadratic substitution; Darboux integrating-factor field;
Piola transport — mechanism (M2).

**Certificate (constructive, strict).** `d = 4`, `D = Q(u,v,w)/9 = X0^4 + X0^2X1^2/3 - 2X0^2X2^2/3 + X1^4/9 - X1^2X2^2/9 + X2^4/9`;
over ℚ(ω): `9D = (X0² + (1/6 - √3 i/6)X1² - X2²/3)·(conj)·9`, i.e. `D = L_1(Φ)·L_2(Φ)/9`, each factor a rank-3
diagonal quadric in `(X0,X1,X2)` (irreducible) **[verified-sympy]**. `Q = 0`; `N0, N1, N2` of degree 5, `N_Y = 0`:
`N0 = X1X2Y1(3X0²-X1²-X2²)/54`, `N1 = -X0X2Y1(3X0²-X1²-X2²)/54`, `N2 = X0X1Y1(3X0²-4X1²-X2²)/54`.

**Mechanism [verified-sympy: `N = -(Y1/216)·adj(dΦ)·E(Φ)` reproduces the witness exactly].**
`E = (2u-v-w)∂_u + (u+v-2w)∂_v` (fixes `w`, `Y`): `E(Q) = 3Q`, `div E = 3`, `E(K) = 0`; so `W = E/Q` is
divergence-free and `W(F̃) = E(S) + 3S = 6u + 3v = 3(2u+v) = 3X2²`. With `Φ = (u,v,w,Y0,Y1,Y2)`,
`Δ = -8X0X1X2`, density `b = -Y1/24`: `V(F) = Δ·b·3X2² = X0X1X2^3Y1 = A`, `div V = 0`.
Clearing: `D = Q(Φ)/9`, `N = D·V = -(Y1/216)adj(dΦ)E(Φ)`. Identities: `N(F) = A·D`, `div N = A/(9X2²)·…`
precisely `div N = X0X1X2Y1/9` and `N(D) = D·div N`.
Target restriction: this family serves `Δ·b·(αu+βv)` only; the fixed target needs the combination `2u+v = X2²`.
The alternative "complementary sextic" `F - 9SD = K` mentioned in solution.md is an (M1)-type radial denominator;
no witness for it is recorded and none was checked here **[not checked]**.

**Symmetry.** Diagonal group order 864 (Smith form `[2,2,6,6,6]`) **[verified-sympy]**.

**Lean.** Template with `t = Y0`, `fFiber = t^6+27`, `dFiber = 7/9`. `Witness.lean` 353 lines / 20.7 KB.

**Relations.** Same target `X0X1X2^3Y1` as hesse-radial (§5) and three-block-rational-d4-strict (#7). Solved
relatives: hesse-triangle, family-hesse-web (#34), klein-hesse (#22), twoties-crossblock-hesse-cofactor (#20).
Open: none (no open instance has a Hesse-type ternary block).

---

## 5. hesse-radial-hesse-cofactor-d4-rational — "Hesse · radial quartic" (PR #16)

**Target.** `q = 1`,
`F = X0^6 + X0^2X2^2Y2^2 + 5X1^6 + 21X1^4Y1^2 + 30X1^2Y1^4 + X2^6 + Y0^6 + 14Y1^6 + Y2^6` (9 monomials),
`A = X0·X1·X2^3·Y1`. Blocks: Hesse-type `X0^6+X2^6+Y2^6+X0^2X2^2Y2^2` in `(X0,X2,Y2)`; non-diagonal binary block
`5X1^6+21X1^4Y1^2+30X1^2Y1^4+14Y1^6`; spectator `Y0^6`.
Structure **[verified-sympy]**: `(u,v,w) = (X0², X2², -Y2²/3)`, `F = S·Q + K`,
`K = (28/27)Y2^6 + 5X1^6+21X1^4Y1^2+30X1^2Y1^4+14Y1^6 + Y0^6` (because `u³+v³+w³-3uvw = X0^6+X2^6-Y2^6/27+X0²X2²Y2²`).

**Method.** Mechanism (M2) with `w` *passive* (`Y2` is never differentiated), passive multiplier `B = X1Y1`,
`Φ = (X0², X1, X2², Y0, Y1, Y2)`, `Δ = 4X0X2`.

**Certificate (constructive, strict).** `d = 4`,
`D = Q(u,v,w) = X0^4 - X0^2X2^2 + X0^2Y2^2/3 + X2^4 + X2^2Y2^2/3 + Y2^4/9`
(= `(9X0^4-9X0^2X2^2+3X0^2Y2^2+9X2^4+3X2^2Y2^2+Y2^4)/9`), product over ℚ(ω) of the two conjugate Hesse lines
`u+ωv+ω²w`, `u+ω²v+ωw` (diagonal rank-3 quadrics in `X0,X2,Y2`); `Q = 0`;
`N0 = (B·X2/18)(3X2²+Y2²)`, `N2 = (X0·B/6)(X2²-X0²)`, others 0; degree 5.

**Mechanism [verified-sympy: `N = (B/36)·adj(dΦ)·E'(Φ)` reproduces the witness exactly].**
`E' = 3[(v-w)∂_u + (v-u)∂_v]`: `E'(Q) = 3Q`, `div E' = 3`, `E'(K) = 0` (`K` has no `u, v`),
`E'(S)+3S = 9v = 9X2²`. `V = (B/36)·adj(dΦ)·(E'/Q)(Φ)`: `V(F) = Δ·(B/36)·9X2² = 4X0X2·X2²·X1Y1/4 = A`, `div V = 0`.
Why the ugly binary block is harmless: the field differentiates only `X0, X2`; `X1, Y0, Y1, Y2` are first
integrals, so the coefficients 5, 21, 30, 14 never enter `N(F)`. `D` omits `X1`, in which `F = 5X1^6 + …`.

**Symmetry.** Diagonal group order 864 (`[2,2,6,6,6]`) **[verified-sympy]**.

**Lean.** Template with `t = X1`, leading coefficient 5 (`(5:ℂ) ≠ 0`), `fFiber = 5t^6+21t^4+30t^2+19`,
`dFiber = 16/9`. `Witness.lean` 347 lines / 19.0 KB.

**Relations.** Same target as fixed-witness-hesse (§4) and three-block-rational-d4-strict (#7). Solved relatives
as in §4. Open: none.

---

## 6. hesse-triangle — "Hesse · three quadrics" (PR #23)

**Target.** `q = 1`, `F = X0^6 + X1^6 + X2^6 + Y0^6 - 3Y0^2Y1^2Y2^2 + Y1^6 + 3Y2^6` (7 monomials),
`A = X0·Y0·Y1^3·Y2`. With `U = Y0²`, `V = Y1²`, `W = Y2²`: `F = C + D`,
`C = U³+V³+W³-3UVW = (U+V+W)(U+ωV+ω²W)(U+ω²V+ωW)` (triangle of lines in the `(U:V:W)` plane = three quadrics in
ℙ⁵), `D = X0^6+X1^6+X2^6+2Y2^6` (radial part).

**Method.** Three-quadric Piola formula — mechanism (M1) with `Φ = (X0,X1,X2,Y0²,Y1²,Y2²)`, `Δ = 8Y0Y1Y2`,
`adj(dΦ) = diag(Δ,Δ,Δ,4Y1Y2,4Y0Y2,4Y0Y1)`.

**Certificate (constructive, strict).** `d = 6`, `D = X0^6+X1^6+X2^6+2Y2^6` (diagonal sextic in four variables,
irreducible over ℂ; cone with vertex line `X0=X1=X2=Y2=0`); `Q = 0`; six numerators of degree 7:
`N_{Xi} = (A/6)·Xi`, `N3 = X0Y1Y2^5/6`, `N4 = X0Y0^3Y2^3/6`, `N5 = X0Y0Y1^3Y2^2/6`.

**Mechanism [verified-sympy: Piola of the field below reproduces the witness exactly].**
Target field `P = (X0V/48)·E_X + (X0W/24)·Z`, `E_X = X0∂X0+X1∂X1+X2∂X2`, `Z = W∂_U + U∂_V + V∂_W`.
`Z(L_k) = ω^k L_k` for the three Hesse lines `L_k = U+ω^kV+ω^{2k}W` **[verified-sympy]**, so
`Z(C) = (1+ω+ω²)C = 0`: the triangle is an eigenfunction decomposition with eigenvalues of trace zero (this is the
"trace-zero" structure that family-trace-zero keeps when the triangle becomes line × conic).
`P(D_t) = (X0V/8)(X0^6+X1^6+X2^6) + (X0W/24)·6W²V = (X0V/8)·D_t`, `div P = X0V(2+1+1+2)/48 = X0V/8`;
so ρ = `b = X0V/8` and (M1) gives `N(F) = A·D`, `div N = A`, `N(D) = A·D` with `A = Δ·b = X0Y0Y1^3Y2`.
In target space `w = P/b = E_X/6 + (W/(3V))·Z` is rational; the density `b` cancels its pole
(`b` has weight `12 - Σ m_i = 3`, consistent with §2.1).

**Unrecorded quartic alternative [verified-sympy, new here].** The (M2) field `E' = 3[(V-W)∂_U+(V-U)∂_V]` with
density `b' = X0/72` gives a second strict rational certificate for the same `(F, A)`:
`D' = Q(Y0²,Y1²,Y2²) = Y0^4+Y1^4+Y2^4 - Y0²Y1² - Y0²Y2² - Y1²Y2²` (degree 4, product of the two non-real Hesse
lines), `N' = (X0/6)·(0, 0, 0, Y1Y2(Y1²-Y2²), -Y0Y2(Y0²-Y1²), 0)` (degree 5), `div N' = X0Y0Y1Y2/3`,
`gcd(F, D') = 1` (`D'` omits `X0`). solution.md mentions "the quartic alternative also discussed in the linked
research record"; that record is not in the repo and I could not confirm it is this tuple. The two witnesses are
not interchangeable (different numerators).

**Symmetry.** Diagonal group order 2592 (`[2,6,6,6,6]`) **[verified-sympy]**; also `S_3` on `(X0,X1,X2)` and
`Y0 ↔ Y1` (the `Y2^6` coefficient 3 breaks full `S_3` on `Y`).

**Lean.** Template with `t = Y0`, `fFiber = t^6-3t^2+7`, `dFiber = 5`. `Witness.lean` 338 lines / 18.7 KB.

**Relations.** Solved deformations: family-trace-zero (§1; triangle → line × conic), fixed-witness-hesse (§4) and
hesse-radial (§5; Hesse cofactor quartics), family-hesse-web (#34; Hesse pencil member `u³+v³+w³-2uvw` in
`(X0²,X1²,Y0²)`), klein-hesse (#22), twoties-crossblock-hesse-cofactor-d4-rational (#20). Open: none.

---

## 7. Relation to the open instances (status as of this checkout: 16 open)

Open ids: cycle-q2, cyclic-triple-residual, dense-rational, diagonal-permutation, eight-plane-configuration,
four-plus-two-residual, klein-coupled-character, off-paired, one-reflection-skew, polarized-k3-remainder,
quadratic-cubic-pullback, quadric-contact-residual, rank-one-biquartic, reflection-coupled-residual,
short-power-presentation, symmetric-sextic-residual.

- **Same variety:** none of the six shares its `F` with an open instance (catalogue `sameVarietyInstances`
  links only fermat-pair ↔ fermat-q2). The pairs named in the task (cycle-three/cycle-q2,
  cyclic-triple-fano/-residual, reflection-coupled-fano/-residual, quadric-contact-broken/-residual,
  one-reflection/-skew, holdout/off-paired) are all pairs *among coupled Fermat deformations*; the Hesse/trace-zero
  varieties of §1, §4–§6 have no open relatives.
- **Deformation relation to the Fermat instances (§2, §3):** every open `F` is the Fermat sextic plus coupling
  monomials (dense-rational: `+ (1/9973)·(236-term Y-even perturbation)`; rank-one-biquartic: `+ (5/7)·(rank-one
  coupling)`; symmetric-sextic-residual, diagonal-permutation: symmetric couplings; eight-plane-configuration,
  short-power-presentation: perturbed diagonal coefficients plus dense terms; klein-coupled-character: Klein
  X-block `X0^5X1+X1^5X2+X2^5X0` with Fermat Y-block). On the undeformed Fermat sextic every open target is served
  by a subset-Euler field (`X0^3Y0^3`, `X0X1X2Y0Y1Y2` and `X0^3Y0^3 - X1^3Y1^3` all with `I = {X0,X1,X2}`,
  where every monomial has `I`-degree 3; cycle-q2's monomial `X0^2X1^2X2^3Y0^2Y1^2Y2` by `I = {X2,Y2}` or by any
  of `{X0,X1}, {Y0,Y1}, {X0,Y0}, {X0,Y1}, {X1,Y0}, {X1,Y1}` — `I`-degree 4 plus `|I| = 2`), so
  fermat-pair/fermat-q2 are the "control" points. **Scan [verified-sympy]: for none of the 16 open instances does any subset `I ⊂ {Z_i}` satisfy
  both `E_I²(F) = 6E_I(F)` and `deg_I A + |I| = 6` for all monomials of `A`.** The couplings destroy
  `I`-homogeneity; this is the precise sense in which "a failed search in one space is not nonexistence": the
  (M1)/(M2)/subset-Euler toolkit of these six does not apply as-is.
- **Specific target coincidences:** fermat-q2's `A` equals holdout-q2's `A` (solved, #124); the open cycle-q2
  is the remaining `q = 2` monomial target on a cyclic deformation (why.md of fermat-q2 frames holdout-q2 and
  block-q2 as the comparisons).
- **Coprimality technique transfer:** the `relPrime_of_specialization ∘ relPrime_of_retract` proof works for any
  future witness whose `D` omits a variable in which `F` is monic; all 16 open `F` contain all six pure sixth
  powers (with constant coefficients), so the technique is available whenever `D` omits a variable.

---

## 8. Independent checks performed in this session

Scripts (outside the repo): `/tmp/hodge/scratch/check6b.py` (all six: `E1`, `E2`, cleared identity, degrees,
homogeneity, parities, legality via the monic-variable argument plus `gcd`, factorization of `D` over ℚ, Smith
normal form of exponent-difference lattices), `mech.py`, `mech2.py` (Piola/Darboux reconstructions),
`fam.py`, `fam2.py` (alternative denominators), `subset_scan.py` (open-instance scan);
instance dump `/tmp/hodge/notes/_instances.txt`. Local sympy 1.12 was used; the `ddg-cas` docker image
(Singular, giac, python-flint; no sympy) was not needed. **No Lean build was run here**; Lean-level facts come
from the accepted records (`acceptance.json`: verifier hashes, Lean 4.33.0, kernel replay, axioms).

Results: all six certificates satisfy `E1 = E2 = 0` exactly (strict), degrees `d + 6q - 5`, required parities,
`gcd(F,D) = 1`. Every stated mechanism reproduces its witness exactly. New computational evidence: hesse-triangle
quartic alternative (§6); fermat-q2 cubic ℚ(i)-denominator with trace equal to the accepted field (§3);
fermat-pair conjugate quartic family (§2); subset-Euler infeasibility for all open instances (§7).

Slips corrected along the way (recorded so nobody repeats them): first reconstruction of §4 used `E` instead of
`E/9` (the record's `D` is `Q/9`); first subset-Euler check in §2 mis-scaled `N` by 6; sympy `cancel` cannot divide
over ℚ(ω) (use `Poly(…, domain=QQ.algebraic_field(sqrt(3)*I))`); `smith_normal_form` needs `domain=ZZ`;
factoring four-variable diagonal sextics over extensions is slow and pointless (they are irreducible).

Conjectures / open: minimal rational denominator degree for fermat-q2's target (≤ 6, cubic exists over ℚ(i));
whether the "complementary sextic" `F - 9SD` of §4 and the "linked research record" quartic of §6 coincide with
the (M1)/(M2) tuples above (plausible, not confirmed).
