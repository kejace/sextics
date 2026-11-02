# Wave 2 search report: `symmetric-sextic-residual`

STATUS: NO WITNESS FOUND (wave 2 started 2026-11-01 ~19:10 EST; report regenerated 2026-11-01 20:55 EST from 472 log records).

Log: `/tmp/hodge/results/symmetric-sextic-residual/wave2-search.log.jsonl` (one JSON line per batch; engine tag `w2lib-v1`). Scripts: `/tmp/hodge/results/symmetric-sextic-residual/w2b/` (`w2lib.py`, `w2_varpro.py`, `w2_kernels.py`, `w2_single*.py`, `w2_struct.py`, `scheduler.py`, `report.py`).

## 0. Corrections to the wave-1 record

* `scan_served.py` (wave 1, last 260 log lines, family `quot-curve6X`): the lines "A ENTERS THE SPAN" and
  "A-system consistent=True" are artefacts. The candidate `1EX(F)+1EY(F)` is `6F` (Euler), so its served
  space is everything and the cumulative span is polluted from that point on; and `S.solve()` was called on
  the system that still contained the target columns, so "consistent" is vacuous. Real content of that
  scan: every tested D of degree 3-6 serves nothing beyond `J(F)_6^G`.

## 1. Plan (wave 2)

1. E1-admissible denominator spaces per degree d = 4..14 in the character sectors triv/sgnX/sgnY/sgnsgn
   of S3xS3 (both Y-parities), plus S2xS2 / S3xS2 / S2xS3 types at d = 10 and the unrestricted kernel at d = 9, 10.
2. Single-D exact (mod p, CRT-lifted) solves of every kernel class + random combinations (never done for
   the non-invariant classes).
3. Character-reduced variable-projection (varpro) joint search over the full E1-admissible family
   (kernel + J^chi_d) per (chi, d, pD), thousands of restarts on the remote.
4. Structured families (cones over plane curves with the three G-stable vertex planes, line-vertex cones,
   discriminant/Hessian-block products, products of kernel elements) via restricted E1 kernels + solves.
5. Rational reconstruction + exact re-solve of every candidate with residual < 1e-10.

(Results sections are appended below as batches finish.)

## 2. Methods built for wave 2 (all in the lane dir; engine files untouched)

* `w2lib.CharFamily` — character-reduced *float* E-system for a linear family `D(a) = Σ a_k B_k` of
  χ-semi-invariant denominators (same orbit folding and column meaning as `charsolve.CharSystem`), with the
  variable-projection objective `f(a) = min_x |M(a)x − b(a)|²` on the unit sphere of the admissible subspace
  (E1 kernel inside the family, computed mod p for the dimension and by SVD for the float basis, minus the
  F-multiples). Gradient by the Golub–Pereyra formula; L-BFGS-B restarts; rational reconstruction with
  increasing denominators; exact confirmation by `charsolve.solve_char` (CRT lift) + `core.verify_witness`.
  Reduction factor ≈ 36 (S3×S3): d = 10 invariant system is 2.6k × 0.9k instead of 85k × 19k.
* `w2lib.JacFamily` — the Jacobian-type ("w-space") search of the task, *unrestricted by symmetry*: `D = Σ w_i F_i`
  (+ the E1-kernel classes beyond J with least-squares E1 lifts), explicit lift `N₀ = A·w/q`, unknowns = Koszul
  corrections `c_ij` (degree d−4) and `Q`; equation `T_D(c) − F·Q = −B_D(N₀)` in degree 2d (book ch. 5, full
  correction criterion; complete at q = 1 because the `F·H` freedom in E1 changes `B_D` by `−HD/q`). Sizes:
  d = 6: 3080 × 388 (K = 18), d = 7: ≈ 5.8k × 0.8k (K = 64), d = 8: ≈ 10k × 1.9k (K = 174).
* Validation (logged as `pipeline-validation`): on the Fermat sextic `p6X + p6Y` both engines recover the
  subset-Euler denominator at d = 6 (residual 1e-24) and the exact solvers lift and verify it; gradient
  checks agree to 8 digits; on the real F the invariant `J^G_6` family (dimP = 3) has a single minimum at
  rel 8.3e-3 and the unrestricted w-space at d = 6 has minima at rel ≈ 0.2.

## 3. Structural observation (reduces the semi-invariant sectors to the invariant one)

For the S3 action on three variables every sign-semi-invariant polynomial is `Δ × (symmetric)`, so
`sgnX`-type D = `Δ_X·S`, `sgnY`-type = `Δ_Y·S`, `sgnsgn`-type = `Δ_XΔ_Y·S` with `S` invariant (the orbit
counts confirm it: the sgnsgn family at (d,pD) = (10,1) has `nD = 12` = number of invariant Y-even quartic
orbits). With the sum principle `served(D₁D₂) = served(D₁) + served(D₂)` (coprime factors) and the fact that
the reflection `s = (X0 X1)` fixes the section `X ∩ {X0 = X1}` pointwise — so every class supported on it is
`s`-invariant (geometric action), i.e. `served(X_i − X_j) ⊂ U = Σ_s H⁴(X)^s`; the residue class of an *invariant*
Y-odd sextic `A'` transforms geometrically by `det(s) = −1` (the residue map twists the polynomial action by det), so
`U` meets the sector of `[A]` — residues of the 17-dimensional space of invariant Y-odd sextics mod J, Hodge
numbers (17,52,17); the dossier's label "(sgn,sgn,+) = (1,10,1)" mixed up the two actions — trivially. One gets, for the target sector,

> `[A] ∈ served(Δ_X^a Δ_Y^b S)  ⇔  [A] ∈ served(S)`.

Hence the sign sectors are redundant with the invariant sector at degree d−3 / d−6, and the wave-1 exclusions
of `D7 = Δ_X e2(X²)` and `D7Y` are *implied* by the degree-4 E1 failure of `e2(X²)`, `e2(Y²)`. Compute is
therefore concentrated on `triv` (both parities), the unrestricted w-space at d ≤ 8, and the subgroup
(S2×S2, S3×S2, S2×S3) types, which probe invariant denominators of much higher degree (orbit products).
The sign sectors are still run at low cost as consistency checks of the principle.

## 3a. What the joint search can and cannot establish

* Every family searched is the *complete* E1-admissible space of its type and degree (kernel of `[A]·` on the
  character sector, including `J(F)_d`, minus the F-multiples), so a zero of the varpro functional is equivalent to
  the existence of a certificate with a denominator of that type and degree (up to the numerical/rational
  reconstruction step, which is always followed by an exact re-solve). Conversely, a nonzero global minimum is
  numerical evidence of non-existence for that type and degree, not a proof: multistart L-BFGS can miss narrow
  basins. The number of restarts, the distribution of the minima found, and the fraction of restarts reaching the
  best minimum are therefore recorded per sector (§7), so a reader can judge the coverage.
* Dimension count: a fixed-D system of the invariant sector at d = 10 has ≈ 2.6k equations and ≈ 0.9k unknowns, so
  consistency imposes ≈ 1.7k conditions on the ≈ 80 admissible parameters; a solution can only exist for a geometric
  reason (a divisor with the right (2,0)-forms), never by counting. The search is a necessary check of the low
  degrees, not an expectation.
* By the orbit-product argument (§3), a certificate with an arbitrary denominator of degree d implies one with a
  G-invariant denominator of degree ≤ 36·d; the subgroup-type and unrestricted w-space sectors are probes of
  invariant denominators of degrees far beyond direct reach.
* The sum principle also shows that products of tested denominators need no separate test: if no single admissible
  D of the listed types serves [A], no product of them does either (served(D₁D₂) = served(D₁) + served(D₂)).

## 6. Concurrent agent on this lane (cross-check, not our work)

Another process submits jobs named `search-symmetric-sextic-residual-{dims,dims2,struct,vp-d8triv0,vp-d8sgnX1}` and
`excl-ssr-relax1` (its scripts `w2x_*.py`, `w2_test2.py`, outputs `w2struct.jsonl`, `w2vp_*.jsonl`, `wave2-excl/` live
only on the remote lane dir). From its outputs (pulled 19:50 EST): its own character-reduced varpro ("cvp") on
`triv d=8 pD=0` (6 restarts) reached rel 8.7e-5 and on `sgnX d=8 pD=1` (22 restarts) rel 0.094 — consistent with
our landscape (invariant minima ≈ 1e-4 at d = 7–8, no zero); its structured E1 scan tested 1014 products of
symmetric functions at d = 6–8 (26 pass E1, none full-solved in that file). Our jobs use the prefix
`search-symmetric-sextic-residual-w2b-*` and the `w2b/` directory to avoid name collisions.
The 26 E1-passing candidates of that scan are all `{EXF, EYF, D7, D7Y} × (symmetric factor of degree ≤ 2)`; the
`EXF/EYF` multiples lie in `J(F)_d^G` (d ≤ 8) and the `D7/D7Y` multiples reduce by §3 to invariant denominators of
degree ≤ 5, so they are inside the invariant varpro families `triv d = 6..8` (both parities) searched jointly below.



## 5. E1-admissible denominator spaces (job `kernels`, mod p = 2^31−1)

`dimK / dimJ / beyond` = dim{D of this type : A·D ∈ J(F)} / dim J(F)_d inside the type / kernel classes modulo J.

| sector (χ, pD) | d=4 | d=5 | d=6 | d=7 | d=8 | d=9 | d=10 | d=11 | d=12 | d=13 | d=14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| triv, pD=0 | 0/0/0 | 1/1/0 | 4/4/0 | 8/8/0 | 22/20/2 | 34/33/1 | 92/64/28 | 122/94/28 | 206/156/50 | 254/211/43 | 370/316/54 |
| triv, pD=1 | 0/0/0 | 1/1/0 | 2/2/0 | 8/8/0 | 14/14/0 | 34/33/1 | 55/51/4 | 122/94/28 | 158/132/26 | 254/211/43 | 310/276/34 |
| sgnX, pD=0 | 0/0/0 | 0/0/0 | 1/1/0 | 3/2/1 | 6/6/0 | 19/12/7 | 26/22/4 | 59/40/19 | 74/61/13 | 123/100/23 | 150/136/14 |
| sgnX, pD=1 | 0/0/0 | 0/0/0 | 0/0/0 | 1/1/0 | 4/3/1 | 8/8/0 | 28/18/10 | 37/31/6 | 75/57/18 | 93/82/11 | 150/132/18 |
| sgnY, pD=0 | 0/0/0 | 0/0/0 | 1/1/0 | 1/1/0 | 6/6/0 | 8/8/0 | 26/22/4 | 37/31/6 | 74/61/13 | 93/82/11 | 150/136/14 |
| sgnY, pD=1 | 0/0/0 | 0/0/0 | 0/0/0 | 3/2/1 | 4/3/1 | 19/12/7 | 28/18/10 | 59/40/19 | 75/57/18 | 123/100/23 | 150/132/18 |
| sgnsgn, pD=0 | 0/0/0 | 0/0/0 | 0/0/0 | 0/0/0 | 0/0/0 | 1/1/0 | 2/2/0 | 14/8/6 | 18/13/5 | 37/30/7 | 46/42/4 |
| sgnsgn, pD=1 | 0/0/0 | 0/0/0 | 0/0/0 | 0/0/0 | 0/0/0 | 1/1/0 | 11/4/7 | 14/8/6 | 30/20/10 | 37/30/7 | 64/56/8 |

Subgroup types at d = 10 (`dimK/dimJ/beyond`): S22_++ pD=0: 400/282/118; S22_++ pD=1: 321/260/61; S22_+- pD=0: 252/182/70; S22_+- pD=1: 251/182/69; S22_-+ pD=0: 252/182/70; S22_-+ pD=1: 251/182/69; S22_-- pD=0: 146/104/42; S22_-- pD=1: 191/123/68; S2S3_++ pD=0: 195/138/57; S2S3_++ pD=1: 134/116/18; S2S3_+- pD=0: 47/38/9; S2S3_+- pD=1: 64/38/26; S2S3_-+ pD=0: 129/96/33; S2S3_-+ pD=1: 107/83/24; S2S3_-- pD=0: 23/18/5; S2S3_-- pD=1: 47/24/23; S3S2_++ pD=0: 195/138/57; S3S2_++ pD=1: 134/116/18; S3S2_+- pD=0: 129/96/33; S3S2_+- pD=1: 107/83/24; S3S2_-+ pD=0: 47/38/9; S3S2_-+ pD=1: 64/38/26; S3S2_-- pD=0: 23/18/5; S3S2_-- pD=1: 47/24/23.

* Served-target check at 2 random admissible points per sector (d ≤ 10): 16 sectors, points serving something beyond `J(F)_6^G`: **0**.
* **Fermat deformation** `F_t = p6X + p6Y + t·e2X e2Y (e2X + 2e2Y)`, t ∈ {0, 5/7, 1, 2, 1/3, −1, 7/5, 3/11, −5/7}: beyond-J kernels identical for all t ≠ 0: **True**; they jump only at t = 0. Examples (t=0 vs t≠0): ('triv', 6, 0): 8 vs 0; ('triv', 6, 1): 2 vs 0; ('triv', 7, 0): 9 vs 0; ('triv', 7, 1): 9 vs 0; ('triv', 8, 0): 21 vs 2; ('triv', 8, 1): 10 vs 0; ('triv', 9, 0): 22 vs 1; ('triv', 9, 1): 22 vs 1; ('triv', 10, 0): 38 vs 28; ('triv', 10, 1): 23 vs 4. So 5/7 is not special for E1.

## 7. Joint (varpro) search landscape — all batches so far

Per sector: admissible dimension `dimP` (= dimK − F-multiples), total restarts (all jobs/seeds), best and median relative residual, number of restarts below 1e-6 / 1e-10, and the distinct minima values (rounded, count). A solution has rel ≤ 1e-16.

| sector | nD | dimK (J+ker) | dimP | rows×cols | restarts | best rel | median | <1e-6 | <1e-10 | minima (value: count) | jobs |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---|---|
| jac d=5 pD=0 | 3 | 3 (3+0) | 3 | 1491×105 | 12 | 5.81e-01 | 5.8e-01 | 0 | 0 | 0.581: 9, 0.592: 3 | local1 |
| jac d=5 pD=1 | 3 | 3 (3+0) | 3 | 1491×105 | 12 | 6.03e-01 | 6.0e-01 | 0 | 0 | 0.603: 9, 0.656: 3 | local1 |
| triv d=7 pD=0 | 37 | 8 (8+0) | 7 | 1268×279 | 10 | 4.34e-04 | 4.3e-04 | 0 | 0 | 0.000434: 10 | local2 |
| triv d=7 pD=1 | 37 | 8 (8+0) | 7 | 1268×279 | 12 | 1.54e-04 | 1.5e-04 | 0 | 0 | 0.000154: 12 | local2,local3 |
| triv d=8 pD=0 | 64 | 22 (20+2) | 18 | 1858×421 | 60 | 6.39e-04 | 6.4e-04 | 0 | 0 | 0.000639: 18, 0.000641: 28, 0.0032: 8, 0.00328: 6 | vp1b |
| triv d=8 pD=1 | 46 | 14 (14+0) | 13 | 1858×457 | 60 | 1.00e-04 | 1.0e-04 | 0 | 0 | 0.0001: 60 | vp1b |
| sgnX d=7 pD=0 | 12 | 3 (2+1) | 3 | 1268×279 | 12 | 1.57e-02 | 1.6e-02 | 0 | 0 | 0.0157: 9, 0.357: 3 | local1 |
| sgnY d=7 pD=1 | 12 | 3 (2+1) | 3 | 1268×279 | 12 | 1.82e-02 | 6.1e-02 | 0 | 0 | 0.0182: 5, 0.0615: 7 | local1 |

Candidates with rel < 1e-10 submitted to exact confirmation: 0; verified witnesses: 0.

### 7a. Sectors still running (partial worker output, deduplicated by file)

| sector | job | restarts so far | best rel | median | <1e-6 | <1e-10 | minima (value: count) |
|---|---|---:|---:|---:|---:|---:|---|
| jac d=6 pD=0 | jac1 | 212 | 1.80e-01 | 2.1e-01 | 0 | 0 | 0.18: 12, 0.187: 5, 0.191: 5, 0.192: 21, 0.2: 1 |
| triv d=9 pD=0 | vp1b | 10 | 3.84e-05 | 3.8e-05 | 0 | 0 | 3.84e-05: 7, 5.76e-05: 2, 0.00017: 1 |

## 8. Single-denominator exact solves (character-reduced systems mod p, CRT lift on consistency)

| family | E1 pass & consistent | E1 pass & inconsistent | E1 fails | errors |
|---|---:|---:|---:|---:|
| kernel-class | 0 | 138 | 0 | 0 |
| kernel-random-combo | 0 | 18 | 0 | 0 |
| structured-single | 0 | 8 | 5 | 0 |

Structured candidates: `DX^2*DY^2` d=12 triv: E1 fails; `EXF*EYF` d=12 triv: inconsistent (rank 1729, null 315, 6713×2044); `HXX` d=12 triv: E1 fails; `HXX+HYY` d=12 triv: E1 fails; `HXX-HYY` d=12 triv: E1 fails; `HYY` d=12 triv: E1 fails; `EXF*D7Y` d=13 sgnY: inconsistent (rank 2010, null 899, 8834×2909); `EYF*D7` d=13 sgnX: inconsistent (rank 2010, null 899, 8834×2909); `D7*D7Y` d=14 sgnsgn: inconsistent (rank 2415, null 1627, 11459×4042); `D7*T8` d=15 triv: inconsistent (rank 4434, null 944, 14672×5378); `D7*T8Y` d=15 sgnsgn: inconsistent (rank 3407, null 1971, 14672×5378); `T8*D7Y` d=15 sgnsgn: inconsistent (rank 3407, null 1971, 14672×5378); `T8*T8Y` d=16 sgnsgn: inconsistent (rank 4691, null 2376, 18564×7067)

Kernel classes / combos per sector (tested, consistent, errors): sgnX d=12 pD=0: 16/0/0; sgnX d=12 pD=1: 21/0/0; sgnY d=12 pD=0: 16/0/0; sgnY d=12 pD=1: 21/0/0; triv d=12 pD=0: 53/0/0; triv d=12 pD=1: 29/0/0

## 10. Batches logged per job

`kernels`: 292 records, 0.1 h, `local-test`: 2 records, 0.0 h, `local1`: 4 records, 0.0 h, `local2`: 2 records, 0.0 h, `local3`: 1 records, 0.0 h, `single`: 169 records, 6.7 h, `vp1b`: 2 records, 0.3 h
