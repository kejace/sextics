# Research dossier: `symmetric-sextic-residual` (OPEN)

Prepared 2026-11-01 from the repository at `/home/kejace/devel/sextics` (branch `main`, HEAD `83fbf5c`),
the GitHub discussion #117 (fetched read-only through the GraphQL API on 2026-11-01: **zero comments**;
the saved file `/tmp/hodge/discussions/symmetric-sextic-residual.md` contains only the title line for the
same reason), the book chapters 01, 03, 04, the accepted `solution.md` of the two quotient constructions
that `why.md` points at (`reflection-coupled-fano`, `cyclic-triple-fano`), and the local, unpublished
computations of the attack session of 2026-11-01 (work directory
`/tmp/hodge/results/symmetric-sextic-residual/`, log `log.jsonl`, plus the scratch directory
`/tmp/hodge/scratch-ssr/` created for this dossier). Nothing in the repository was modified; nothing was
posted.

Status labels used below:

* **[P]** proven / exact computation over Q (hence over C when the statement is rational).
* **[C]** computational evidence only (modular rank over p = 2^31 − 1, screened ansatz, bounded search).
* **[H]** heuristic, conjecture, or structural reading not proved.
* **[X]** attempted route that failed or stalled.
* **[N]** new in this dossier (computed on 2026-11-01 during its preparation; scripts in `/tmp/hodge/scratch-ssr/`).

Engine caveat that applies to everything below: the fixed-denominator engine
(`~/devel/sextics-work/hodge/fixed_d.py`) had a sign bug in the `F·N(D)` term until 2026-11-01 ~17:25.
Every number quoted here was produced **after** the fix (`[hodge engine v2]`) or by the independent
symmetry-reduced solver `invsolve.py` (`engine=inv-v1`), which agrees with the corrected engine on every
cross-check and reproduces the Fermat subset-Euler certificate exactly (`invsolve.py` self-test).

---

## 1. Exact target

### 1.1 Statement

* `id: symmetric-sextic-residual`, title "Permutation symmetry · residual sector", `q = 1`,
  variables ordered `X0, X1, X2, Y0, Y1, Y2`.
* `A = X0*X1*X2*Y0*Y1*Y2 = e3(X)·e3(Y)` (one monomial, degree 6 = 6q, Y-degree 3: Y-odd as required).
* `F` has 42 monomials (degree 6, Y-even). **[P]** (`ssr_struct.py`, exact comparison) it is exactly

```
F = p6(X) + p6(Y) + (5/7) · e2(X) · e2(Y) · ( e2(X) + 2·e2(Y) )
  = X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6 + (5/7) e2(X)^2 e2(Y) + (10/7) e2(X) e2(Y)^2 ,
```

  with `e2(X) = X0X1+X0X2+X1X2 = Q_X`, `e2(Y) = Q_Y`, `p6 = sum of sixth powers`. The monomial
  coefficients are 5/7 (type X_i^2X_j^2 Y_kY_l), 10/7 (types X_i^2X_jX_k Y_kY_l and X_iX_j Y_k^2Y_l^2) and
  20/7 (type X_iX_j Y_k^2Y_lY_m); exactly as listed in `instance.rec`.
* Sparse presentation in the weighted quotient coordinates (e1,e2,e3 of each block): with
  `p6 = e1^6 − 6e1^4e2 + 6e1^3e3 + 9e1^2e2^2 − 12e1e2e3 − 2e2^3 + 3e3^2`,

```
F = G(e1X,e2X,e3X,e1Y,e2Y,e3Y),   G = p6(eX) + p6(eY) + (5/7) e2X e2Y (e2X + 2 e2Y),
```

  a weighted sextic in P(1,2,3,1,2,3). The coupling involves **only** e2X and e2Y.
* Degrees for a certificate at q = 1 with `deg D = d`: `deg N_i = d + 1`, `deg Q = 2d − 6`, `deg E = 2d + 6`
  (`E = A D^2 − D Σ N_i ∂_iF + F(D Σ ∂_iN_i − Σ N_i ∂_iD) = F^2 Q`). Parities: with `pD` the Y-parity of D,
  the X-slot numerators have Y-parity `1 − pD`, the Y-slot numerators `pD`, and Q is Y-odd.

### 1.2 Smoothness **[P]**

`checks.json`: `smoothReductionPrime = 5`, reduced partials `X_i^5, Y_j^5`. Indeed 5/7 ≡ 0 mod 5, so
`F ≡ Fermat (mod 5)` and the Jacobian ideal mod 5 is `(X_i^5, Y_j^5)`: smooth reduction, hence smooth over
Q (and C). This is the same mechanism as in the two accepted "coupling 5/7" constructions
(`reflection-coupled-fano`, `cyclic-triple-fano`), and it is the reason the whole family
`F_t = Fermat + t·e2X e2Y(e2X + 2e2Y)` is a natural deformation axis (§5, experiment E6).

`targetNonzeroInReducedJacobianRing = true`: A = Z^(1,1,1,1,1,1) has all exponents ≤ 4, so it is a
nonzero monomial in the Fermat Jacobian ring; by semicontinuity [A] ≠ 0 in R(F)_6, i.e. the class is a genuine
primitive (3,1)-class (`D = 1` sanity row in `log.jsonl`: rank 18 vs augmented 19 — A ∉ J(F)).

### 1.3 Symmetry group and the character of the residue

* Visible linear symmetries: `S3(X) × S3(Y)` (independent permutations of the two blocks) and
  `σ = diag(−1,−1,−1,1,1,1)` (projectively the same as `τ: Y ↦ −Y`). Projective group
  `G = S3 × S3 × ⟨σ⟩`, order 72.
* **[N][P]** The Lie algebra of the linear stabiliser of F (solutions of `Σ M_ij Z_j ∂_iF = c·F`) is
  one-dimensional (the Euler field only): the linear automorphism group is finite (`ssr_struct.py`).
* **[P]** No block exchange: `F(λY, μX) = F(X,Y)` would need `λ^4μ^2 = 2`, `λ^2μ^4 = 1/2`, `λ^6 = μ^6 = 1`,
  impossible (ratio `λ^2/μ^2 = 4` contradicts `λ^6 = μ^6`). The coefficients 5/7 vs 10/7 break the X↔Y symmetry.
* **[P]** Projective diagonal group: `{1, σ}` only (`checks.json`: `projectiveDiagonalGroupOrder = 2`,
  consistent with the exponent-lattice rule of Chapter 1 §"support-coset rule": the mixed monomials
  generate all within-block transfers and `2(e_X − e_Y)`).
* **[P]** No coordinate reflection (`coordinateReflections: []`), no subset-Euler atom
  (`subsetEulerAtoms: []`: F does not split as G(X)+H(Y)), Hessian symmetrizer space scalar
  (`hessianSymmetrizerDimension = 1`: no hidden linear direct-sum decomposition), no junior diagonal element
  (`juniorElementsFixingTarget: []`, but note the screen only looked at the diagonal group of order 2).
* **Residue character [P].** For a signed permutation `g = (π_X, π_Y, signs)`,
  `A(gZ) = (product of signs)·A(Z)` and `vol(gZ) = det(g)·vol = (product of signs)·sgn(π_X)sgn(π_Y)·vol`.
  Hence the residue `Res(A·Ω/F^2)` transforms by

```
χ_res(g) = χ_A(g)·det(g) = sgn(π_X)·sgn(π_Y)        (independent of σ; τ-character +1).
```

  So the target lives in the **(sgn ⊗ sgn, τ = +1)** isotypic component. Every transposition
  (a reflection in a non-coordinate hyperplane `X_i = X_j` or `Y_i = Y_j`) acts by −1 on the residue
  although it fixes A; `why.md` stresses exactly this.
* The quotient map by `S3 × S3` has weights (1,2,3,1,2,3) and Jacobian determinant `±Δ_X Δ_Y`,
  `Δ_X = (X0−X1)(X0−X2)(X1−X2)`. Its density character is (sgn,sgn) as a polynomial, hence
  **trivial** as a residue: the symmetric quotient serves the (triv,triv) sector (§2.1).

### 1.4 Isotypic Hodge numbers **[N][P]**

`mult.py` (session script, rerun here) computes the `S3×S3×⟨τ⟩`-equivariant Hilbert function of the
Jacobian ring in the Fermat monomial model. Because J(F) is generated by a regular sequence transforming
as `V^*`, the equivariant Koszul resolution shows the character of `R_d` is the same for every smooth
G-invariant F, so the table is exact for this F. With `R_6 ≅ H^{3,1}_prim`, `R_12 ≅ H^{2,2}_prim`,
`R_18 ≅ H^{1,3}_prim` and Y-odd numerator ⇔ τ-invariant residue:

| sector (ρ_X, ρ_Y, τ) | h^{3,1} | h^{2,2} | h^{1,3} | note |
|---|---:|---:|---:|---|
| (sgn, sgn, +) **← target** | **1** | **10** | **1** | [A] spans h^{3,1}; K3-type Hodge structure of rank 12 |
| (sgn, sgn, −) | 0 | 6 | 0 | purely (2,2) |
| (triv, triv, +) | 17 | 52 | 17 | served by the symmetric Fano quotient (§2.1) |
| (triv, triv, −) | 26 (+1 for h^{4,0}) | 67 | 26 | contains the (4,0) class |
| (sgn, triv, +) / (triv, sgn, +) | 5 / 5 | 22 / 22 | 5 / 5 | anti-invariant for one block only |
| (std, std, +) | 21 | 92 | 21 | dimension counts multiplicities ×4 |

Total Y-odd R_6 = 208, R_12 = 866 (so H^4_prim(X)^τ has Hodge numbers (208, 1751−885 = 866, 208), consistent
with the full h^{3,1} = 426 = 218 + 208).

Consequences:

* **[P]** The residue class of A is the **unique** (3,1)-class of its character. A certificate for A is
  equivalent (Chapter 1, prescribed-support theorem) to divisor support of the whole transcendental part
  `T` of the rank-12 K3-type Hodge structure `H^4(X)_{(sgn,sgn,+)}`. Conversely, support for *any* class
  of this sector with a nonzero (3,1) component supports [A] (the minimal sub-Hodge structure of the
  (3,1)-line is T, and N^1 is a sub-Hodge structure).
* **[H]** The generalised Hodge conjecture predicts the sector is in N^1 (it has no (4,0) part); the
  challenge is the construction.
* **[P]** Under `S3 × S3` (reflection group generated by the six transposition reflections `s`),
  Chapter 3 §3.3 "Finite groups generated by reflections": `U = Σ_s H^4(X)^s` is divisor-supported and
  `H^4/U` is exactly the (sgn,sgn) isotypic component. The all-target conclusion needs `τ` or `−τ` in the
  reflection group; here τ ∉ S3×S3, and the τ-invariant part of H^4/U is precisely the (1,10,1) sector above.
  **The residual sector of this instance is exactly the sector the reflection-group theorem cannot reach.**

### 1.5 Fermat endpoint **[P]**

At t = 0 (`F_0 = p6(X) + p6(Y)`) the same target is served strictly by the subset-Euler certificate
(Chapter 1; reproduced by `invsolve.py`):

```
D = p6(X) = X0^6+X1^6+X2^6,   N_i = A·X_i/6 (i = 0,1,2),  N_j = 0 (j = 3,4,5),  Q = 0,
```

since `N(F_0) = A·E_X(F_0)/6 = A·p6(X)` and `div N = A − A = 0`. `served_check.py` **[N][C]** confirms
that on Fermat `D = p6(X)` serves 11 of the 19 invariant Y-odd sextic orbits (including A), whereas on the
actual F the same D fails E1 (`inv_scan2.out`: `p6X E1=False`). The Fermat certificate is the certificate
incarnation of the Shioda–Katsura join structure (`p6(X) = 0 = p6(Y)` is the join of two Fermat sextic
curves); the coupling term destroys the split.

---

## 2. Companion and related accepted instances

`provenance.json`: `sameVarietyInstances: []` — no other instance uses this sextic. `why.md` cites
Chapter 3 §§3.3–3.4 (reflection quotients, nonlinear Piola/Fano maps). The relevant accepted applications
and why they do not apply:

### 2.1 The symmetric Fano quotient (Chapter 3 §3.4, no accepted instance uses it) **[P]**

`Φ = (e1X,e2X,e3X,e1Y,e2Y,e3Y)`, degrees (1,2,3,1,2,3), `M = 12`, `δ = det dΦ = ±Δ_XΔ_Y`,
`F = G(Φ)`. The Piola theorem applies (`7 ≤ M ≤ 12`, X smooth): `W = {G = 0} ⊂ P(1,2,3,1,2,3)` is a
normal Fano fourfold (`−K_W = O(M−6) = O(6)`) and every target `A' = δ·b(Φ)` with `wt(b) = 6(q+1) − 12`
is divisor-supported. At q = 1, `wt(b) = 0`: the construction serves **only multiples of Δ_XΔ_Y**
(`why.md`: "that density construction serves a multiple of Δ_X·Δ_Y, not the symmetric monomial chosen
here"). At q = 2 it serves `Δ_XΔ_Y·b(e)`, `wt(b) = 6`: still the (triv,triv) residue sector. By the
averaging-ideal criterion (§3.1) a quotient by K serves ρ iff `ρ^K ≠ 0`; for K = S3×S3 only ρ = triv.

### 2.2 `reflection-coupled-fano` (accepted, `fqlx`, PR #150) — reflection quotient

Same coupling 5/7 (a different sextic). Quotient by `Y2 ↦ −Y2`, `u = Y2^2`, `H(x,u) ⊂ P(2,1,1,1,1,1)`;
char-5 reduction `H_0 = Σ x_i^6 + u^3`, explicit weighted quartic over F_25 with `j^2 = 2`
(`x0 = j(s^4+t^4), x1 = (2+j)(s^3t+st^3), x2 = js^2t^2, x3 = j(s^3t−st^3), x4 = s^4−t^4, u = s^8+s^4t^4+t^8`),
invertible 24×24 deformation matrix, étale lifting, polynomial-family integration, finite trace, legal
clearing. **No explicit D or degree is exported** (existence proof). Here: every transposition is such a
reflection (after the linear change `u = X_i − X_j`), so the reflection theorem applies to each of the six
transpositions — but it serves the targets whose residue is **invariant** under the reflection, i.e. the
sectors with `sgn_X = +1` or `sgn_Y = +1`. Our sector is anti-invariant under all six. **[P]** Not applicable.

### 2.3 `cyclic-triple-fano` (accepted, PR #153, line by `rbajaj5`) — cubic Piola map

`Z = X2^3`, weights (1,1,1,1,1,3), free weighted line `[0:s:0:is:t:it]`, 6×6 incidence Jacobian with
determinant −22500/2401, polynomial-source integration, cubic transport
`D = φ(d), Ñ_j = X2^2 φ(N_j), Ñ_{X2} = φ(N_Z)/3`. It serves the density sector `X2^2·B(…, X2^3, …)`.
Analogue here: the nonlinear polynomial atoms available are the symmetric functions (§2.1) and
products such as `(e2X, e2Y)`-type maps; all have `δ` with a (sgn,·) polynomial character and hence serve
invariant residue sectors. **[P]** Not applicable to (sgn,sgn).

### 2.4 The "residual" siblings (all OPEN)

`reflection-coupled-residual`, `cyclic-triple-residual`, `four-plus-two-residual`,
`quadric-contact-residual` are the same design pattern (a companion construction serves one character; the
residual instance asks for the complementary one). None has an accepted proof; their discussions
(see `/tmp/hodge/notes/dossier-cyclic-triple-residual.md` and the saved threads) establish the same
obstruction vocabulary: quotient densities carry a fixed character; "changing the source B does not change
the character"; the needed transfer is a correspondence that is **not** a quotient. The present instance is
the extreme case: the symmetry group is the largest (order 72) and the residual sector is the
one-dimensional-(3,1) K3-type piece of §1.4.

### 2.5 Fermat-type accepted instances

`fermat-q2`, `fermat-pair-rational-d4-strict`, `holdout*` (Fermat cover machinery of `fqlx`, PR #128) use
the split structure `G(X)+H(Y)` or a Fermat cover. The split is destroyed here; `subsetEulerAtoms: []`.

---

## 3. Prior work

### 3.1 Public record

* Discussion #117 (`https://github.com/hodge-computer/sextics/discussions/117`), opened 2026-10-17 by
  `benediamond` (standard template body). **No comments, no replies** as of 2026-11-01 18:00 (GraphQL
  `comments.totalCount = 0`). There are no linked packets, no posted exclusions, no conjectures.
* `why.md` / `checks.json` (curated screens, see §1.3): no reflection, no subset Euler, no linear split, no
  junior diagonal element; symmetric quotient serves Δ_XΔ_Y; transpositions serve the opposite character;
  "averaging the target into the invariant sector erases it".
* DiamondJM's cross-lane sweep of ~1300 "natural degree-6 supports" at q = 1 (memory note; see the
  dense-rational / off-paired dossiers for the packet links) reported nothing for any open lane; the
  present instance was included only through generic families, not through its symmetry.

### 3.2 Local unpublished work (attack session 2026-11-01, `kejace` agent; not posted)

All entries are in `/tmp/hodge/results/symmetric-sextic-residual/log.jsonl` (131 lines; fields
`D, d, family, e1, consistent, rank, nullity, seconds, note`). "E1" = the necessary condition
`A·D ∈ J(F)` (Chapter 1 normalisation: `A D − q N(F) = F H` is solvable iff `A D ∈ J(F)` since F ∈ J).

**(a) E1 kernel sweeps — all characters** (`e1_kernel.py`, mod p = 2^31−1, `e1_kernel.out`).
`ker_R(d, pD) = dim { D ∈ S_d^{pD} : A·D ∈ J(F) } / J(F)_d`:

| d | pD = 0 (Y-even) | pD = 1 (Y-odd) | kernel elements beyond J |
|---:|---:|---:|---|
| 1–6 | 0 | 0 | — (all E1-passing D of degree ≤ 6 lie in J(F)) |
| 7 | 1 | 1 | `D7 = Δ_X·e2(X^2)`, `D7Y = Δ_Y·e2(Y^2)`, with `e2(X^2) = X0^2X1^2+X0^2X2^2+X1^2X2^2 = e2X^2 − 2e1Xe3X` |
| 8 | 6 | 6 | includes `D7·X_i`, `D7Y·Y_j` (parity 0) resp. `D7·Y_j`, `D7Y·X_i` (parity 1); structure of the rest: see E1 below |

  **Rigour.** A zero kernel mod p implies a zero kernel over Q provided the Jacobian columns have their full
  complete-intersection rank mod p in degree d+6 (rank can only drop mod p, so `dim_p ≥ dim_Q`). The solver
  reports the J-rank; it equals the CI value in every degree used. So **[P]**: every certificate denominator
  not in J(F) has degree ≥ 7, and in degree 7 it is (up to J) one of `D7`, `D7Y`.
  Inconsistency statements below are **[C]** (a rational solution could in principle have p in its
  denominators; with p ≈ 2^31 this is practically excluded but not proved).

**(b) Invariant E1 kernel sweep** (`e1_kernel_inv_sweep.py`, `invsolve.e1_kernel_inv`, d ≤ 16;
`e1inv_d*_p*.pkl`): `ker = dim (J:A)_d^G / J_d^G` for G-invariant D:

| d | 1–7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| pD=0 | 0 | **2** | **1** | 28 | 28 | 50 | 43 | 54 | 41 | 43 |
| pD=1 | 0 | 0 | **1** | 4 | 28 | 26 | 43 | 34 | 41 | 28 |

  Comparison with the representation-forced kernel of `[A]: R_d^{(triv,triv)} → R_{d+6}` (`mult.out`,
  `m_d − m_{d+6}`): forced 0 for d ≤ 9, then 28 (d=10, pD=0), 2 (d=10, pD=1), 28, 28, 50, 26, … So the
  kernels at d = 8, 9 (2, 1, 1) and the excess 2 at (10, pD=1) are **F-specific ("accidental")**; from
  d = 11 on the invariant E1 kernel is entirely the dimension-forced one. Exact rational representatives of
  the accidental classes: `exactK.py` (CRT over 5 primes, e-monomial coordinates) — output pending at the
  time of writing (see §7 addendum if present).

**(c) Fixed-D full solves**, all **[C]** inconsistent (mod p):

| family | D (degree) | solver | rows × cols | rank / nullity | time |
|---|---|---|---|---|---|
| linear, quadrics, cubics, quartics (103 denominators: all symmetric functions, pencils `QX + cQY`, `p3 + c·e3`, `Δ_X`, `Δ_Y`, `Δ_X e1Y`, mixed cubics, K3 cones `X0^4+X1^4+Y0^4+Y1^4`, …) | d ≤ 4 | E1 only | — | **E1 fails for every one** (consistent with (a)) | ≤ 2 s each |
| e1-kernel-deg7 | `D7 = Δ_X e2(X^2)` | full engine v2, fast modular | 25935 × 4497 (N 3861, Q 636) | 3991 / 506, **inconsistent** | 343 s |
| e1-kernel-deg7 | `D7Y = Δ_Y e2(Y^2)` | same | same | 3991 / 506, **inconsistent** | 593 s |
| J6-radial / pencil | `E_X(F)`, `E_Y(F)`, `aE_X(F)+bE_Y(F)` (14 ratios) | inv-v1 | 838 × 161 | 140 / 21, inconsistent | ≤ 1 s |
| **[N]** J-inside, invariant | `ΣF_{X_i}`, `ΣF_{Y_j}` (d=5, = J^G_5) | inv-v1 | 532 × 103 | 87 / 16, inconsistent | 0.2 s |
| **[N]** J-inside, invariant | the full 4-dim `J^G_6` (basis `E_X(F)`, `e1X ΣF_X − E_X(F)`, `E_Y(F)`, `e1Y ΣF_Y − E_Y(F)`; 10 structured + 6 random points) | inv-v1 | 838 × 161 | 140 / 21, inconsistent | ≤ 2 s |
| E1-kernel-inv-d8 | 2 kernel reps, 4 random kernel combinations, 6 random points of `(J:A)_8^G` incl. J-part, 3 random `J_8^G` | inv-v1 | 1858 × 421 | 362 / 59, inconsistent | 4–30 s |
| E1-kernel-inv-d9 | 1 kernel rep, 4 random kernel combinations, 5 random points of `(J:A)_9^G` | inv-v1 | 2648 × 673 | 566 / 107, inconsistent | 10–67 s |
| **[N]** fermat-nbr, disc6, sym6 | `p6X`, `p6Y`, `p6X±p6Y`, `Δ_X^2±Δ_Y^2`, `e3X^2+e3Y^2`, `p3X^2±p3Y^2`, `QX^3±QY^3`, `e3Xe3Y (=A)`, `p3Xp3Y` | inv-v1 | — | **E1 fails** | 0.1 s |

  The nullities (506 at d = 7, 59 and 107 in the invariant systems) are the usual Koszul/Euler gauge
  freedom, not evidence of near-solutions.

**(d) "Sum of supports" span scan** (`served_scan.py`): for 92 invariant denominators of degree 2–6 compute
the space `T(D)^G` of invariant Y-odd sextic targets served by D (targets as unknowns in the invariant
system) and accumulate the span, using `HasCertificate.add` (D1 serves A1, D2 serves A2 ⇒ D1D2 serves A1+A2).
The session run crashed before printing; **[N]** rerun here, and a scaling bug found and fixed
(`served_scan_fixed.py`: a coefficient v on the orbit-representative column m stands for the invariant target
`v·(36/|orbit m|)·orbit_sum(m)`; the session script omitted the orbit-size factor). Corrected result **[C]**:

> every one of the 92 invariant denominators (all symmetric functions of degree ≤ 6, all pencils,
> `Δ_X^2`, `Δ_Y^2`, `E_X(F)`, `E_Y(F)`, `p6X`, `p6Y`, `e3Xe3Y`, …) serves **exactly** `J(F)_6^{G,odd}`
> (dimension 2: `e1Y·ΣF_{X_i}` and `e1X·ΣF_{Y_j}`) and **nothing beyond J** among the 19 invariant
> Y-odd sextic orbits; the cumulative span never leaves J, so [A] (sector dimension 17 = 19 − 2) is not
> reached by any sum.

  Sanity: on Fermat the same code finds `T(p6X)^G` of dimension 11 (including A) and `T(QX)^G = J^G`.

**(e) Representation-theoretic lower bound** (`mult.py`): see §1.4 and (b).

**(f) Structured E1 scans [N]** (`e1_struct8.py`, mod p): within `Δ_X·Sym^4(X)` (dim 4) the E1 kernel is
2-dimensional: `e2^2 − 4e1^2e2 + 2e1^4` and `e1e3 − 2e1^2e2 + e1^4` (in `e_k = e_k(X)`), whose
combination `e2^2 − 2e1e3 = e2(X^2)` is `D7`; the other direction lies in `J_7`. Symmetric in Y.
Degree-8/9 families (`Δ_X·Sym^5(X)`, `Δ_X·Sym^3(X)Sym^2(Y)`, `Δ_XΔ_Y·(quadrics)`, `Δ_X·Sym^6`,
`Δ_X^3`, …) were still running when this dossier was written (background job; results, if finished, are in
`/tmp/claude-1000/.../tasks/bn79dn9nu.output` and should be appended to §7).

**(g) Non-invariant J-inside full solves [N]** (`jtests.py full`: `F_{X0}`, `F_{Y0}`, `F_{X0}+F_{X1}`,
`X0F_{X0}`, `Y0F_{Y0}`, `X0F_{Y0}`, `Y0F_{X0}`, `X1F_{X0}`; full engine v2) — running in background at
writing time (`/tmp/hodge/scratch-ssr/jtests_full.out`); see §7.

### 3.3 New structural results of this dossier

**(h) Age screen over the full 72-element group [N][P]** (`ages.py`). For every `g ∈ G` with
`χ_res(g) = +1` (i.e. [A] ∈ H^4(X)^g), every fixed component `P(V_λ)` of g on P^5 that meets X, and the
age of g on `T_pX` there (ambient eigen-angles `(μ−λ) mod 1` with one zero removed, minus the normal
weight `−6λ`, as in `audit-frontier.py`):

| class of g (ρ = +1) | fixed components meeting X | age on T_pX |
|---|---|---|
| identity | X | 0 |
| σ (= τ) | two plane sextic curves `{X=0}`, `{Y=0}` | 3/2 |
| (3-cycle, 1), (1, 3-cycle) | sextic surface `X ∩ {X0=X1=X2}` (resp. Y) | **1** (isolated eigenpoints `[1:ω:ω^2:0:0:0]` have F = 3 ≠ 0, not on X) |
| (3-cycle, 3-cycle) | three lines `P(V_1), P(V_ω), P(V_{ω^2})` (6 points each) | 2 |
| (transposition, transposition) | sextic surface `X ∩ {X0=X1, Y0=Y1}` and 6 points on `P(V_{−1})` | 1 and 2 |
| σ·(t,t') | curves in `P(V_{±1})` (P^2's) | 3/2 |
| σ·(3-cycle, 1) | curve `X ∩ {X=0}` | 3/2 |
| σ·(3-cycle, 3-cycle) | no fixed point on X | — |

**No element has a fixed point of age < 1.** Hence (Kollár–Larsen, Theorem 2, in both directions as used in
Chapter 4 §4.4): for every subgroup `K ⊂ G` with `[A] ∈ H^4(X)^K` (equivalently `K ⊂ ker χ_res`, an
index-2 subgroup of order 36 containing `A3×A3`, the nine `(t,t')`, and the σ-multiples) the quotient
`X/K` has canonical singularities, numerically trivial canonical class, and is **not uniruled**. The
averaging ideal of §3.1 generated by uniruled quotients of X cannot serve the (sgn,sgn) sector: the
situation is the exact analogue of the cubic-cover table in §4.4 (`V_std` not served by any uniruled
intermediate quotient). Consistency check: the quotient `V = X/ker χ_res` is the double cover of the
symmetric Fano `W` branched along `{disc_X·disc_Y = 0}`, a weighted complete intersection of degrees (6,12)
in P(1,2,3,1,2,3,6) with `K_V = O(6+12−12−6) = O(0)` **[H]** (weighted adjunction, V not quasi-smooth), and
`H^4(V)_− ≅ H^4(X)_{(sgn,sgn)}` is a (1,10,1) ⊕ (0,6,0) piece: a Calabi–Yau-type fourfold whose
anti-invariant middle cohomology is of K3 type.

**(i) Cone-over-plane-curve structure [N][P].** For a ternary denominator `D = D(X0,X1,X2)` the divisor
`E_D = X ∩ V(D)` fibres over the plane curve `C = {D = 0} ⊂ P^2_X` with fibre over x the sextic surface

```
S_x = { p6(x)·s^6 + p6(Y) + (5/7)e2(x)^2 s^4 e2(Y) + (10/7) e2(x) s^2 e2(Y)^2 = 0 } ⊂ P^3_{(s:Y)},
```

which depends on x only through `λ = p6(x)/e2(x)^3` (rescale s): all cone divisors over plane curves are
pulled back from the single `S3(Y)`-symmetric pencil `S_λ = {λ s^6 + p6(Y) + (5/7)s^4 e2(Y) + (10/7)s^2 e2(Y)^2}`
(plus the Fermat-type fibre `p6(x)s^6 + p6(Y)` over `e2(x) = 0`). `p_g(S_λ) = 10`, and as an `S3(Y)`-module
with the residue twist `H^{2,0}(S_λ) ≅ Sym^2(⟨s⟩ ⊕ V_Y) ⊗ sgn_Y` contains **sgn_Y with multiplicity 3**
(from `s^2`, `e1Y^2`, `e2Y`). So the Y-side character is available on every fibre; the X-side character
sgn_X must come from the base curve C and the twist of the Hodge bundle along it. **The unique degree-7
E1-passing denominators are exactly of this type**: `D7 = Δ_X·e2(X^2)` cuts the three reflection planes
and the cone over the 3-nodal (rational) quartic `e2(X^2) = 0` — the sgn_X character appears through Δ_X.
`D7` fails the full system (mod p), so (prescribed-support theorem, `D^q` suffices, q = 1) the class is
**not** supported on `E_{D7}` **[C]**.

---

## 4. What remains unexplored

1. **Degrees 8–9, all characters, exactly.** The full E1 kernel at d = 8 (6 classes per parity) and d = 9
   has never been full-solved except through the invariant subspace (2 + 1 + 1 classes, random points).
   Non-invariant kernel directions (characters (sgn,triv), (triv,sgn), (sgn,sgn), (std,·)) are untested.
   Exact rational representatives (e-monomial form) of the accidental invariant classes are not yet recorded.
2. **Denominators inside J(F) of degree 5–7, non-invariant** (`F_{X0}`, `X_iF_{X_j}`, …): running at
   writing time (§7). Invariant ones (J^G_5, J^G_6) are excluded [C].
3. **Semi-invariant denominators** with characters (sgn,triv), (triv,sgn), (sgn,sgn): `D7·D7Y` (degree 14,
   character (sgn,sgn)), `D7^2`, `D7·Δ_X`, `Δ_XΔ_Y·(invariant)`, `D7·(invariant)` … The invariant solver
   rejects them (`is_invariant` fails), the full solver is too slow beyond degree 8. A **character-twisted
   equivariant solver** (N transforming by `χ_D ⊗ tangent`) has not been written.
4. **"Sum of supports" at degrees 7–12.** The corrected span scan only covers degree ≤ 6 (where nothing
   beyond J is served). The served spaces `T(D)^G` for D in the invariant E1 kernels at d = 8–12 and for
   semi-invariant D were never computed.
5. **Deformation from Fermat.** The family `F_t = Fermat + t·e2X e2Y(e2X+2e2Y)` with the strict
   certificate `(p6X, A X_i/6, 0)` at t = 0 has not been continued: neither the E1 kernel as a function of t
   (where does `p6X` stop passing E1? at generic t the invariant degree-6 kernel beyond J is 0, so the
   degree must jump), nor a perturbative/Q(t) certificate search at degrees 7–12.
6. **Non-rigorous exclusions.** All "inconsistent" results are mod p. No Farkas-type rational certificate of
   inconsistency has been produced for D7, D7Y or the degree-6 J-elements.
7. **Geometric (nonconstructive) routes.** No candidate divisor E with `H^{2,0}(Ẽ)_{(sgn,sgn)} ≠ 0` has
   been identified. The cone-over-curve mechanism (3.3(i)) has not been analysed beyond degree 7: the degrees
   of the sgn_Y-Hodge subbundle of the pencil `S_λ` (rank 3) and the resulting `h^{2,0}` for symmetric plane
   curves C of degree 7–10 with the sgn_X twist are not computed.
8. **Correspondences beyond quotients**: the discriminant-type double cover `V`, the branch divisor
   `{disc_X = 0} ∩ W` (a quotient of a reflection section), relative/twisted cohomology
   `H^4(W∖B, L_{Δ})`. Nothing started.
9. **The full linear automorphism group** is only bounded below (order 72 projectively, finite). A larger
   finite group (e.g. a non-monomial element normalising S3×S3) was not searched for; the Lie algebra is
   trivial, so a finite search over the stabiliser of F's Hessian or of the 42-monomial support would settle it.

---

## 5. Candidate attack plan

Sizes below are for q = 1, `deg D = d`: unknowns `#N = 6·#S_{d+1}^{par} ≈ 3·C(d+6,5)`,
`#Q = #S_{2d−6}^{odd} ≈ C(2d−1,5)/2`, equations `#S_{2d+6}^{odd} ≈ C(2d+11,5)/2`. Invariant reduction
divides rows and columns by ≈ 36 (orbit counts are somewhat larger because of small orbits).

| d | #N | #Q | rows (full) | invariant rows (observed, E1 system) | full-solve cost (engine v2, dense modular rref) |
|---:|---:|---:|---:|---:|---|
| 7 | 3861 | 636 | 25 935 | 268 | 343 s (measured) |
| 8 | 6006 | ≈1 500 | ≈40 400 | 322/376 | ≈ 10–20 min (est.) |
| 9 | 9009 | ≈3 000 | ≈62 000 | 448 | ≈ 1 h (est.) |
| 10 | 13 104 | ≈5 800 | ≈85 000 | 532/604 | hours; use invariant/twisted reduction |
| 12 | 25 704 | ≈16 800 | ≈162 000 | 838/932 | invariant only (seconds) |
| 14 | 46 512 | ≈40 400 | ≈288 000 | 1268/1388 | invariant only |

**E1. Close degrees 8–9 exactly, all characters.** [cheap, 1 day] For each parity, take the 6 (d=8) and the
full (d=9) kernel representatives `(J:A)_d / J_d` (`e1_kernel_d8_p*.pkl`, extend to d=9), decompose into
`G`-isotypic pieces (apply the 36 permutations, project), and full-solve each isotypic class (and random
points of each class) with engine v2. For classes with a nontrivial character, use the twisted solver of E2.
Expected: inconsistent; value: a complete "no denominator of degree ≤ 9" statement, which is what the
dense-rational lane obtained and what a future poster will ask first. Chance of a witness: < 2 %.

**E2. Character-twisted invariant solver.** [1–2 days] Generalise `invsolve.InvSystem` to semi-invariant D
(`g·D = χ_D(g) D`): if a certificate exists for D then averaging `g ↦ χ_D(g)^{-1}·(φ_g^* v)` yields one with
`N_{g(j)} = χ_D(g)·(g·N_j)`, so N is determined by `N_0` (invariant under `Stab(0) = S2 × S3` twisted by
χ_D) and `N_3`, and `Q` transforms by `χ_D^2 = 1`. Then test, in order:
`D7·D7Y` (d = 14, χ = (sgn,sgn), Y-odd), `D7^2` (d = 14, invariant, Y-even), `D7·e1X`, `D7·e1Y`,
`D7·Δ_X` (= Δ_X^2 e2(X^2), d = 10, invariant), `Δ_XΔ_Y·{e2X, e2Y, e1X^2, e1Y^2}` (d = 8),
`Δ_X·Sym^k(X)·Sym^l(Y)` kernels (from `e1_struct8.py`), and `D7·(invariant E1-kernel classes of degree 8, 9)`.
Rationale: products of E1-passing denominators are the natural candidates for `HasCertificate.add`-type
witnesses and are the only degree-14 candidates with a reason to exist. Chance: 5–10 %.

**E3. Served-span scan at degrees 7–12 (sum of supports).** [1 day, after E2] For each class D in the
(twisted-)invariant E1 kernels at d = 7…12 compute `T(D)^{χ}` (targets as unknowns, 19 invariant odd
sextic orbits), accumulate the span over all tested D, and test whether [A] mod J enters the span. Report the
served sector dimensions per degree (currently: 0 beyond J for everything of degree ≤ 6). If [A] enters a
span of served spaces of D1,…,Dk, the witness is `D = ΠD_i` with the added fields (Chapter 1 addition
formula; `HodgeComputer.CertificateOperations`). Chance: 10 % — this is the mechanism most likely to produce
an explicit witness if one of moderate degree exists, because single-D searches miss sums.

**E4. Deformation from Fermat / certificate over Q(t).** [2–3 days] Work with `F_t`. (i) For d = 6…10
compute the invariant E1 kernel beyond J over `F_p(t)` by evaluating at ~10 values of t and interpolating
ranks: find the special t where `dim (J_t : A)_d^G / J` jumps (t = 0 is special at d = 6 since `p6X` is in the
kernel there; t = 5/7 shows the accidental classes at d = 8, 9 — are they special to 5/7 or generic?). If the
d = 8, 9 classes exist for generic t, continue the **certificate** itself: solve the fixed-D system over
`F_p[t]` with D(t) the t-family of kernel representatives, looking for t-values where it becomes consistent
(a consistent t ≠ 0 would give a witness for the corresponding sextic; 5/7 may or may not be among them).
(ii) Perturbative continuation of the Fermat certificate at degree 6k: `D = p6X^k + t D_1 + …` is excluded at
k = 1 by the degree-6 kernel; test k = 2 (d = 12, invariant: rows 838/932) directly at t = 5/7 by solving the
invariant system with D in the span of the invariant E1 kernel at d = 12 (50/26 classes): this is the exact
"does any invariant degree-12 denominator work" question and costs seconds per random point. Chance: 10 %
for a witness at d ≤ 12; high value either way (explains whether 5/7 is generic).

**E5. Rational inconsistency certificates.** [1 day] For D7, D7Y and the 4-dimensional J^G_6 family, compute
the left kernel of the coefficient matrix mod several primes, reconstruct a rational row vector y with
`yM = 0`, `y·b ≠ 0`, and verify over Q. Upgrades the degree-7 exclusion from [C] to [P], and gives a template
for posting rigorous exclusions. No chance of a witness, but it is what a discussion post should contain.

**E6. Cone-over-plane-curve divisors (geometric, then algebraic).** [1 week] Using 3.3(i): compute the
Hodge bundle of the pencil `S_λ` restricted to the `sgn_Y`-part (rank 3; Griffiths residues
`q(s,Y)Ω_3/F_s` with `q ∈ {s^2, e1Y^2, e2Y}` and their λ-dependence), its splitting degrees on P^1_λ, and
for `S3(X)`-semi-invariant plane curves `C = {Δ_X·s(X) = 0}` (s symmetric of degree 4, 5, 6, …) the
dimension of `H^{2,0}(Ẽ_C)_{(sgn,sgn)} = H^0(C, φ^*(Hodge bundle)^{sgn_Y} ⊗ twist)^{sgn_X}`, where
`φ = p6/e2^3 : C → P^1`. The first C with a nonzero answer is the candidate support; then test
`D = Δ_X·s(X)` (and `D·Δ_Y·s'(Y)`) in the certificate system — E1 first (cheap, structured families as in
`e1_struct8.py`), then the twisted solver. Note the E1 computation already singled out `e2(X^2)` at degree 4
for s, which is suggestive but failed the full system. Chance: 15 % for a witness of degree 8–12; this is
the only route in which the target's two sign characters have an identified geometric origin.

**E7. The double cover V and its branch divisor (nonconstructive).** [open-ended] Show that
`T ⊂ H^4(X)_{(sgn,sgn,+)}` is supported on a divisor by a correspondence that is not a quotient:
e.g. the relative cohomology `H^4(V)_− ≅ H^4(W, B)`-type sequence for the double cover `V → W` branched along
`B = B_X ∪ B_Y`, `B_X = (X ∩ {X0 = X1})/Stab` (a quotient of a smooth reflection section, a threefold with
`h^{3,0}` possibly ≠ 0), and a Lefschetz/Gysin argument expressing the (2,0)-part of the sector through
`H^3(B)`. If `H^4(V)_−` is generated by boundary/Gysin classes from B, then the support is the union of the
six reflection sections — but those serve only reflection-invariant classes (Chapter 3 §3.3 argument), so
either the argument fails at the Hodge level (expected: `(1,10,1)` cannot sit in `Σ_s H^4(X)^s`) or it
identifies the missing correspondence. Worth one careful day to settle on paper; chance of a complete proof
this way: < 5 %, but it is the natural "why is this sector special" statement for a discussion post.

**E8. Enlarge the symmetry group / find hidden structure.** [half a day] Compute the finite linear
automorphism group of F exactly (stabiliser of the 42-monomial support under the normaliser of the torus is
monomial, so search monomial matrices with 6th-root entries: 6^6·720 cases, trivial) and check for
non-monomial elements by solving `F(MZ) = F(Z)` with M in the normaliser of `S3×S3` (e.g. `M` acting as a
rotation in the `(e1X, e1Y)` plane). A larger group would change the character analysis and possibly supply
a junior element. Chance of changing the picture: 10 %.

**E9. q = 2 shadow.** [1 day] Although the instance fixes q = 1, Chapter 1 §3.6 ("from a computed support to
every q") shows a support found at q = 2 transfers to q = 1. The symmetric Fano quotient serves
`Δ_XΔ_Y·b(e)` at q = 2 (`wt b = 6`); the twisted solver at q = 2 (degrees `d + 7` numerators) with
semi-invariant D is a different linear system in which Δ-type denominators enter more naturally. Low
priority; chance 5 %.

---

## 6. Feasibility assessment

**Explicit witness.** Everything of degree ≤ 7 is excluded (E1 rigorously for d ≤ 6 outside J, modular for
d = 7 and for the invariant J-elements of degrees 5–6). The invariant sector is barren through degree 9
(random sampling) and all 92 natural invariant denominators of degree ≤ 6 serve nothing beyond J. The
representation theory says the invariant E1 kernel becomes large from degree 10 on (28, 28, 50, …), so a
witness, if of moderate degree, is most likely **invariant or (sgn,sgn)-semi-invariant of degree 10–14**
(products such as `D7·D7Y`, `D7^2`, or an accidental degree-12 class) and found by E2–E4 at a cost of hours
of CPU, not days. The full non-reduced system at d = 14 (288k × 87k) is out of reach for dense rref but the
twisted-invariant reduction (≈ 8k × 2.5k) is trivial. If E2–E4 fail through degree 14, the honest statement
is "no denominator of degree ≤ 14 with a character of the 72-element group", and the explicit route becomes
expensive (the nonlinear search over a 28–50-dimensional E1 kernel at each degree is a variety problem, not a
linear one; random points only see the generic stratum).

**Nonconstructive argument.** The two support mechanisms of the book are dead for this sector by exact
computation: reflection quotients serve the complementary characters (§1.4, [P]), and no quotient by any
subgroup of the 72-element group fixing the residue is uniruled (3.3(h), [P]). A proof therefore needs a new
correspondence — a divisor `E ⊂ X` with `H^{2,0}(Ẽ)_{(sgn,sgn)} ≠ 0` whose Gysin image has nonzero
(3,1)-part. The cone-over-plane-curve family (3.3(i)) is the only identified geometric source in which both
sign characters appear, and E1 independently selected a member of it (D7) at the first admissible degree;
this is the route to develop. The K3-type shape (1,10,1) of the sector suggests looking for a surface with
`p_g = 1` (K3 or abelian) in the `S3×S3`-geometry — e.g. a `(sgn,sgn)`-eigen-surface inside the symmetric
quotient W or the double cover V — but no such surface has been identified.

**What a Lean proof would require.** (a) Explicit witness: `proof.json` with D, N_0…N_5, Q over Q, the
protected `Solution.lean` template (`HodgeComputer.HasCertificate F A q`), a coprimality proof for D and F
over C (for a (sgn,sgn)-semi-invariant D such as `D7·D7Y`, factor by the reflection planes and prove each
factor coprime to the irreducible F — `F` is irreducible over C by the Eisenstein specialisation used in the
accepted 5/7 proofs: set all but `X0` and one `Y` to zero), and the kernel replay of the cleared identity at
degree 2d+6 (for d = 14: ≈ 288k monomials — within the size of accepted dense witnesses only if the
coefficient bit-size is moderate; a (sgn,sgn)-symmetric witness can be verified orbit-wise but the
protected statement needs the full polynomial identity; `leangen.py`/`leancheck.py` pipeline verified on
2026-11-01). (b) Nonconstructive: an existence proof in the style of PR #150/#153 needs a geometric curve or
family on an auxiliary variety with a square incidence seed, a polynomial normal, dominance, finite trace,
and target transport back to F — all of which presuppose the (not yet identified) correspondence; the
shared library (`AveragingProjectors`, `PolynomialCurves`, `ResidueTransport`, `CertificateOperations`)
covers the algebra once the geometry is specified.

---

## 7. Reproduction and addenda

* Session scripts: `/tmp/hodge/results/symmetric-sextic-residual/{h.py, invsolve.py, e1_kernel.py,
  e1_kernel_inv_sweep.py, families1.py, test_d7.py, test_K.py, exactK.py, served_scan.py, inv_scan.py,
  struct8.py, mult.py}`; outputs `*.out`, kernels `*.pkl`, log `log.jsonl`.
* Dossier scripts: `/tmp/hodge/scratch-ssr/{ages.py, ssr_struct.py, served_check.py, served_scan_fixed.py,
  jtests.py, e1_struct8.py}`.
* Engine: `~/devel/sextics-work/hodge/{core.py, fixed_d.py}` (v2, sign fix 2026-11-01 17:25), python-flint,
  primes `2147483647, 2147483629, …`.
* GitHub: discussion #117 fetched read-only (`gh api graphql`); raw JSON at
  `/tmp/hodge/discussions/symmetric-sextic-residual.raw.json`.

Addenda (results that finished after the main text was written) are appended below.
