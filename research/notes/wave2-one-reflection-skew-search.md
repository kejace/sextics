# Wave 2 — `one-reflection-skew` joint search (report rewritten after every batch)

STATUS: NO WITNESS FOUND (so far). Searches are running on the remote under extreme contention (box load
900–1300 from ~60 other jobs; capacity gate at load 150 never opened, so one 4-process chain was forced
through and a second one followed; two more chains remain queued). Pipeline validated end-to-end on an accepted
witness (family-hesse-web): exact 6→5 reduction, operator (rel 2.6e-15 at the known D), LM variable projection
finds the known D from random starts (≈15–50 % of restarts at K = 32), snap-to-rational + exact confirm
recovers exactly the accepted denominator and verifies the 6-variable lift (`core.verify_witness`).

Summary of results so far (details below): (1) exact E1-admissible dimensions for all four (ρ,σ) sectors at
d = 4..12 (table §1; non-Jacobian directions first appear at d = 7); (2) the Jacobian-type degrees d = 5, 6
(6-variable denominators (w·∇f)² of degree 10/12) are excluded numerically in every sector (global minima of
the variable-projection objective ≥ 6.5e-2 over tiny families, found 70–124 times out of 160 restarts);
(3) degree-6/7 coordinate-cone families are empty or without zeros; (4) kernel-family searches at d = 7–12,
structured candidates at d = 6–12, the half-target A⁺ and the solved-neighbour control are in progress.

## 0. Setting (exact reduction used throughout)

`F = f(X0,X1,X2,Y0,Y1) + Y2^6`, `A = X0^3 Y0^3` (free of Y2), `q = 1`. Averaging a certificate over
`μ6(Y2)` and expanding in powers of `F` (details in §6) gives the **exact equivalence**

> a certificate exists  ⇔  a rational 5-variable field `u = N/D` (deg D = d, deg N_i = d+1) satisfies
> `A D² − D N(f) + (6/5) f (D div N − N(D)) = 0`     (*)

and then `(D², (N_i D)_{i<5}, Y2 (D div N − N(D))/5, Q = 0)` is a strict 6-variable certificate of degree
`2d`. Necessary condition (E1): `A D ∈ J_f`. `D` can be taken in any `(ρ, σ)` character
(`ρ`: sign of (Y0,Y1), `σ`: sign of (X1,Y1)); `N` and the rows then live in determined sectors.
Verified on `family-hesse-web` (accepted degree-6 witness, same target): the reduced field satisfies (*)
exactly, the float operator gives relative residual `3.7e-15` at it, `exact_confirm` recovers and
verifies the 6-variable certificate.

Cross-reference: the exclusion lens of this lane (`wave2-one-reflection-skew-excl.md`, same lane directory,
log `wave2-excl.log.jsonl`) independently proved the *strict* Y2-free form (fringewidth's statement): for Y2-free
`D` a certificate exists iff some Y2-free `N'` of degree `d+1` has `N'(f) = A D` and `D div N' − N'(D) = 0`. The two
forms are equivalent: given `N` solving (*), `g := (N(f) − AD)/f` is a polynomial and `N' = N − (g/6)E'` (`E'` the
5-variable Euler field) satisfies the strict system (`E'(g) = d·g`), and conversely strict solutions solve (*). Their
exact degree-7 count `(J_f : A)_7 = 84 = 75 + 9` over all sectors matches the sector table below (28+19+18+19 = 84,
non-Jacobian 7+1+0+1 = 9). Their rigorous degree-7 Gröbner/msolve exclusion work is complementary to the numerical
joint search here.

Two structural facts used in the search design (both one-line checks from (*)): (i) `(fG, fW)` solves (*) iff
`(G, W)` does, so `f`-multiples are never spurious; (ii) more generally `(qD, qN)` solves (*) for every polynomial
`q` whenever `(D, N)` does, so **the degree-`d` solution set contains `D*·S_{d−d₀}` for every solution `D*` of
lower degree `d₀`** — a degree-`d` kernel search covers all lower degrees (as positive-dimensional families, which
the snap-to-rational step then has to collapse), and "no zero at `d`" is heuristic evidence for all `d' ≤ d`.

Code: `/tmp/hodge/results/one-reflection-skew/w2lib.py` (operator linear in D over monomial coordinates,
variable projection with Tikhonov normal equations, full Golub–Pereyra Jacobian LM, rationalisation +
LLL, exact modular confirm with CRT), `w2_search.py` (multiprocess restarts), `w2_dims.py`,
`w2_struct.py`, `chains/chain{1,2,3,4}.sh`. Log: `wave2-search.log.jsonl` (one line per batch).

## 1. E1-admissible denominator spaces (task 1) — computed (`w2_dims.py`, remote, mod p = 2^31−1 and float SVD agree)

Five-variable problem, `D` in character sector `(ρ,σ)`: `kerD = dim{D : A D ∈ J_f}` / `dim J_f,d` / non-Jacobian
directions (`ker mod J`) / `f·S_{d−6}` inside; `rk[A]` = rank of multiplication by `[A]: (S/J)_d^{(ρ,σ)} → (S/J)_{d+6}`.

| d | nD (0,0) | (0,0) ker/J/extra/fmult | rk[A] | (0,1) ker/J/extra | (1,0) ker/J/extra | (1,1) ker/J/extra |
|---|---|---|---|---|---|---|
| 4 | 22 | 0/0/0/0 | 22 | 0/0/0 | 0/0/0 | 0/0/0 |
| 5 | 36 | 2/2/0/0 | 34 | 1/1/0 | 1/1/0 | 1/1/0 |
| 6 | 60 | 7/7/0/1 | 53 | 6/6/0 | 6/6/0 | 6/6/0 |
| 7 | 90 | **28/21/7**/2 | 62 | 19/18/1 | 18/18/0 | 19/18/1 |
| 8 | 135 | **91/46/45**/6 | 44 | 76/43/33 | 67/43/24 | 76/43/33 |
| 9 | 190 | 161/92/69/11 | 29 | 146/86/60 | 141/86/55 | 146/86/60 |
| 10 | 266 | 250/161/89/22 | 16 | 229/153/76 | 223/153/70 | 229/153/76 |
| 11 | 357 | 349/259/90/36 | 8 | 328/247/81 | 325/247/78 | 328/247/81 |
| 12 | 476 | 473/387/86/60 | 3 | 445/371/74 | 442/371/71 | 445/371/74 |

Jacobian-type degrees (kernel = `J_d`): `d ≤ 6` in every sector and `d = 7` in `(1,0)`. The first non-Jacobian
directions appear at `d = 7` (7 in the main sector `(0,0)`, 1 each in `(0,1)`, `(1,1)`). Consistency: sector sums
reproduce wave-1's ρ-even/ρ-odd totals (47/37 at d=7, 167/143 at d=8, …). Six-variable cross-check with the
engine's `jacobian.mult_kernel` (Y-even D): kernel 0, 3, 18, 71, 247, 572 at d = 4..9 with J_d = 0, 3, 18, 63, 168, 378
(non-Jacobian 0, 0, 0, 8, 79, 194), matching the dossier's totals (135 = 126 + 9 at d = 7, 1144 = 756 + 388 at d = 9 over both parities).
Note: `(S/J_f)_{d+6}` has dimension 255, 185, 121, 70, 35, 15 for d = 7..12, so `E1` becomes almost vacuous
at high degree (at d = 12 only 3 conditions remain in the main sector).

### Exact rational kernel bases (auto-generated; `w2lib.e1_kernel_exact`, cached in `w2cache/`)

<!-- AUTO-KERNELS -->
| d | sector | target | kernel dim | primes used | max height (bits) | seconds |
|---|---|---|---|---|---|---|
| 6 | (0,0) | A | 7 | 3 | 6 | 0.7 |
| 7 | (0,0) | A | 28 | 13 | 181 | 7.1 |
| 7 | (0,0) | Aplus | 21 | 3 | 13 | 5.4 |
| 7 | (0,1) | A | 19 | 3 | 17 | 1.2 |
| 7 | (0,1) | Aplus | 19 | 3 | 17 | 5.3 |
| 7 | (1,0) | A | 18 | 3 | 11 | 1.3 |
| 7 | (1,1) | A | 19 | 8 | 104 | 3.3 |
| 8 | (0,0) | A | 91 | 52 | 674 | 49.5 |
| 8 | (0,0) | Aplus | 66 | 13 | 184 | 52.0 |
| 8 | (0,1) | A | 76 | 52 | 682 | 45.4 |
| 8 | (0,1) | Aplus | 55 | 8 | 103 | 31.6 |
| 8 | (1,0) | A | 67 | 28 | 415 | 25.1 |
| 8 | (1,1) | A | 76 | 52 | 680 | 47.1 |
| 9 | (0,0) | A | 161 | 79 | 1137 | 150.6 |
| 9 | (0,0) | Aplus | 141 | 42 | 525 | 332.8 |
| 9 | (0,1) | A | 146 | 79 | 1132 | 74.7 |
| 9 | (0,1) | Aplus | 126 | 28 | 407 | 128.1 |
| 9 | (1,0) | A | 141 | 79 | 997 | 78.9 |
| 9 | (1,1) | A | 146 | 79 | 1130 | 99.9 |
| 10 | (0,0) | A | 250 | 122 | 1526 | 412.5 |
| 10 | (0,0) | Aplus | 237 | 79 | 1110 | 2274.3 |
<!-- /AUTO-KERNELS -->

Symmetry breaking measured on the E1 kernel (sector (0,0), 5 variables, mod p): the undeformed
`one-reflection` sextic has kernel 20, 46, 98, 164 at d = 6, 7, 8, 9 versus 7, 28, 91, 161 for the skew sextic.
The skew term thus removes 13, 18, 7, 3 admissible directions; these survive numerically as near-kernel
directions (singular values 1e-7 … 1e-15 of the float E1 operator), which is why exact kernel bases are needed
and why a fixed-D test of a *numerically* admissible D can be misleading (E1 must be checked exactly).

## 2. Optimizer calibration (needle = accepted degree-6 witness of family-hesse-web, same target A)

* First attempts (L-BFGS, LM-Kaufman, LM with the full Golub–Pereyra Jacobian — the latter checked
  against finite differences to 3e-9) all converged to recurring local minima (rel 3.056e-2, 4.256e-2,
  2.3e-3, …), 0/26 hits. Basin probe: LM started 5°–45° from the needle still ended at rel 5e-2.
* Cause: the family had the f-multiple subspace `f·S_{d−6}` projected out (as the engine's varpro does
  with `exclude_F_multiples`). The true `D*` lies in the E1 kernel (projection residual 8e-16) but has a
  component along `f`, so the restricted span did not contain it. **Removing the F-multiple subspace is
  wrong**: it removes solutions `D0 + F·G` that are not divisible by F. Recorded in
  `/tmp/hodge/notes/engine-bugs.md`; the lane driver now searches the full kernel (in the strict
  5-variable problem there are no spurious solutions: `(fG, fW)` solves iff `(G, W)` does).
* Sector bookkeeping check (all four `(ρ,σ)` sectors at d = 7, random 12-term `D`): the fixed-D operator reproduces
  the trivially served Jacobian targets `Z_j f_i` exactly (residual ≤ 9e-17 for the explicit solution
  `N = D·Y0·e_X0` of target `Y0 f_X0`), and the served-target analysis returns exactly the 6-dimensional
  Jacobian part with nothing beyond it (`served_mod_J = 0`), as expected for a random `D`.
* With the full kernel (K = 32): **4/8 LM restarts converge to the needle** (rel 2e-13 … 6e-8,
  |cos(u,u*)| = 1.0000); the other restarts end at the usual local minima. So random-restart LM has a
  ~50% per-restart hit rate on a 32-dimensional family with an 8-term solution in that small local sample.
  Remote calibration `calib-hesse-d6` (same family, LM vs L-BFGS alternating, **400 restarts, final**):
  lm: n=200  <1e-08: 24 (0.120 [0.08,0.17])  <1e-06: 30 (0.150 [0.11,0.21])  <0.0001: 30 (0.150 [0.11,0.21])  median s/restart 14.2; lbfgs: n=200  <1e-08: 5 (0.025 [0.01,0.06])  <1e-06: 32 (0.160 [0.12,0.22])  <0.0001: 35 (0.175 [0.13,0.23])  median s/restart 23.6. L-BFGS approaches the zero and stops short of polishing, which the snap step's LM refit then
  finishes, so the effective hit rate of both methods is ≈ 0.15–0.19 per restart at K = 32, with LM ~2× faster.
  Every confirmed hit (22 confirmations, all from τ-sparsified refits) is the accepted 8-term D, verified exactly.
  Working estimate: ≈ 0.15 per restart, so
  (Caveat on `calib-xtype-d6`: that job was configured in sector (0,1) but the accepted xtype-strict-rational-113
  denominator lies in sector (0,0), so its needle is not in the searched family; its result is not a calibration
  datum. A correctly configured rerun `calib-xtype-d6-r0s0` is queued in chain 5.)
  N = 100 restarts miss a comparable needle with probability ≈ 1e-7 … 1e-9.

## 3. Search batches (filled from the log; remote chains `w2-chain1b` (4 proc), `w2-chain2` (4 proc, cones/structured), `w2-chain3` (8 proc, kernels d = 9–12, started 01:54 UTC); box load 900–1700 throughout)

Jacobian-type degrees (task 2; 5-variable `d = 5, 6`, i.e. 6-variable denominators `(w·∇f)²` of degree 10, 12,
never excluded before): full-E2 variable projection over the complete `w`-space, 160 restarts per sector
(LM + L-BFGS alternating; the families have K ≤ 7 so this is exhaustive to numerical precision).

| d | sector | K | restarts | best rel | median rel |
|---|---|---|---|---|---|
| 5 | (0,0) | 2 | 160 | 3.13e-1 | 3.13e-1 |
| 5 | (0,1) | 1 | 160 | 5.36e-1 | — |
| 5 | (1,0) | 1 | 160 | 4.55e-1 | — |
| 5 | (1,1) | 1 | 160 | 6.35e-1 | — |
| 6 | (0,0) | 7 | 160 | 6.50e-2 | 7.6e-2 |
| 6 | (0,1) | 6 | 160 | 3.08e-1 | 3.6e-1 |
| 6 | (1,0), (1,1) | 6 | running | | |

No zero anywhere (a solution shows rel ≤ 1e-10 in the calibration). A "Hessian relaxation" was unnecessary:
the full E2 system over these tiny `w`-spaces is solved directly by the variable projection.

Cost model (one unloaded core, LM with full Golub–Pereyra Jacobian, ≤150 evaluations per restart): d = 6:
~5 s/restart (system 1815×410); d = 7: 19 s (2640×615, 0.22 s per evaluation, K = 28); d = 8: ~1 min (3740×900);
d = 10 ≈ 5 min; d = 12 ≈ 20 min (11594×2975). On the shared remote (load 900–1300 during this run) wall times
were 2–4× longer. Calibration hit rate at K = 32 (family-hesse-web needle): LM ≈ 15–50 % per restart,
L-BFGS ≈ 5–10 % (final numbers in the table).

Early observation at d = 9 (`ker-d9-r0s0`, K = 161, 8 processes): the first 16 LM restarts end at rel
1.9e-3 … 1.7e-2 with a recurring attractor at rel ≈ 2.4e-3 (6/16) — an order of magnitude below the best
minima at d ≤ 8 (6.5e-2 … 0.6), as expected from the larger family, but still seven orders of magnitude
away from a zero; ~40 % of the restarts hit the 150-evaluation cap (≈ 2 min per restart under load).

### Full table of searches (auto-generated from `wave2-search.log.jsonl`)

<!-- AUTO-SEARCH -->
| job | instance | d | sector | family | target | K | system rows×cols | method | restarts | best rel | median rel | #<1e-8 / <1e-6 / <1e-4 | wall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| jac-d5-r0s0 | one-reflection-skew | 5 | (0,0) | jacobian | A | 2 | ?×? | both | 160 | 3.13e-01 | 3.13e-01 | 0 / 0 / 0 | 27.5 s |
| jac-d5-r0s1 | one-reflection-skew | 5 | (0,1) | jacobian | A | 1 | ?×? | both | 160 | 5.36e-01 | 5.36e-01 | 0 / 0 / 0 | 4.5 s |
| jac-d5-r1s0 | one-reflection-skew | 5 | (1,0) | jacobian | A | 1 | ?×? | both | 160 | 4.55e-01 | 4.55e-01 | 0 / 0 / 0 | 4.5 s |
| jac-d5-r1s1 | one-reflection-skew | 5 | (1,1) | jacobian | A | 1 | ?×? | both | 160 | 6.35e-01 | 6.35e-01 | 0 / 0 / 0 | 3.4 s |
| cone-X0X1Y0Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X1,Y0,Y1 | A | 3 | 1815×410 | lm | 120 | 1.88e-01 | 1.88e-01 | 0 / 0 / 0 | 135.1 s |
| cone-X1X2Y0Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X1,X2,Y0,Y1 | A | 1 | ?×? | lm | 120 | 6.31e-01 | 6.31e-01 | 0 / 0 / 0 | 7.8 s |
| jac-d6-r0s0 | one-reflection-skew | 6 | (0,0) | jacobian | A | 7 | ?×? | both | 160 | 6.50e-02 | 7.59e-02 | 0 / 0 / 0 | 163.1 s |
| jac-d6-r0s1 | one-reflection-skew | 6 | (0,1) | jacobian | A | 6 | ?×? | both | 160 | 3.08e-01 | 3.62e-01 | 0 / 0 / 0 | 178.7 s |
| jac-d6-r1s0 | one-reflection-skew | 6 | (1,0) | jacobian | A | 6 | ?×? | both | 160 | 2.93e-01 | 3.85e-01 | 0 / 0 / 0 | 191.7 s |
| jac-d6-r1s1 | one-reflection-skew | 6 | (1,1) | jacobian | A | 6 | ?×? | both | 160 | 3.09e-01 | 3.17e-01 | 0 / 0 / 0 | 212.8 s |
| cone-X0X1Y0Y1-d7 | one-reflection-skew | 7 | (0,0) | cone:X0,X1,Y0,Y1 | A | 6 | ?×? | lm | 120 | 1.04e-01 | 1.49e-01 | 0 / 0 / 0 | 667.6 s |
| cone-X1X2Y0Y1-d7 | one-reflection-skew | 7 | (0,0) | cone:X1,X2,Y0,Y1 | A | 4 | 2640×615 | lm | 120 | 4.45e-01 | 4.45e-01 | 0 / 0 / 0 | 1400.9 s |
| cone-X0X1X2-d8 | one-reflection-skew | 8 | (0,0) | cone:X0,X1,X2 | A | 1 | 3718×890 | lm | 120 | 6.14e-01 | 6.14e-01 | 0 / 0 / 0 | 154.5 s |
| cone-X0X1X2Y0-d8 | one-reflection-skew | 8 | (0,0) | cone:X0,X1,X2,Y0 | A | 11 | 3718×890 | lm | 16 (running) | 2.47e-01 | 2.53e-01 | 0 / 0 / 0 | 605.9 s |
| cone-X0Y0Y1-d8 | one-reflection-skew | 8 | (0,0) | cone:X0,Y0,Y1 | A | 1 | 3718×890 | lm | 120 | 5.96e-01 | 5.96e-01 | 0 / 0 / 0 | 169.9 s |
| ker-d9-r0s0 | one-reflection-skew | 9 | (0,0) | kernel | A | 161 | 5096×1246 | lm | 32 (running) | 1.86e-03 | 2.44e-03 | 0 / 0 / 0 | 663.6 s |
| calib-hesse-d6 | family-hesse-web | 6 | (0,0) | kernel | A | 32 | ?×? | both | 400 | 2.51e-14 | 2.25e-02 | 12 / 33 / 3 | 3454.7 s |
| calib-xtype-d6 | xtype-strict-rational-113 | 6 | (0,1) | kernel | A | 28 | 1815×410 | lm | 32 (running) | 4.30e-03 | 6.63e-03 | 0 / 0 / 0 | 202.9 s |
| cone-X0X1X2-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X1,X2 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X1X2Y0-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X1,X2,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X1X2Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X1,X2,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X1Y0-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X1,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X1Y0-d8 | one-reflection-skew | 8 | (0,0) | cone:X0,X1,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X2Y0-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X2,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X2Y0-d8 | one-reflection-skew | 8 | (0,0) | cone:X0,X2,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X2Y0Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X2,Y0,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X2Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,X2,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0X2Y1-d8 | one-reflection-skew | 8 | (0,0) | cone:X0,X2,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X0Y0Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X0,Y0,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X1X2Y0-d6 | one-reflection-skew | 6 | (0,0) | cone:X1,X2,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X1X2Y0-d8 | one-reflection-skew | 8 | (0,0) | cone:X1,X2,Y0 | A | 0 | | | 0 | — | — | | empty family |
| cone-X1Y0Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X1,Y0,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X1Y0Y1-d8 | one-reflection-skew | 8 | (0,0) | cone:X1,Y0,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X2Y0Y1-d6 | one-reflection-skew | 6 | (0,0) | cone:X2,Y0,Y1 | A | 0 | | | 0 | — | — | | empty family |
| cone-X2Y0Y1-d8 | one-reflection-skew | 8 | (0,0) | cone:X2,Y0,Y1 | A | 0 | | | 0 | — | — | | empty family |

Total restarts on the skew lane so far: 2048.
<!-- /AUTO-SEARCH -->

Kernel families (task 3) at d = 7 … 12 and the half-target/cone/structured runs: queued in the chains; the
E1-kernel bases are now computed **exactly** (mod-p nullspace, canonical rref, CRT over 31-bit primes,
rational reconstruction verified by an extra prime; heights 2^17 … 2^674 bits for d ≤ 8) because the float
SVD kernel is unreliable here: the skew term is a 1/7-perturbation of `one-reflection`, so many E1-kernel
directions of the undeformed sextic survive as near-kernel directions with singular values 1e-7 … 1e-15
(the float basis in the main sector at d = 7 was off by 1.3e-5; in sector (0,1) the 19-dimensional kernel
was not even separable from 4 spurious directions).

### How to read the numbers

`rel = ‖M(D)x − b(D)‖/‖b(D)‖` is the exact least-squares residual of the fixed-D system (SVD) at the end of a
restart; a genuine solution gives `rel ≲ 1e-10` (calibration: 1e-13 … 1e-8 before polishing). Values of
order 1e-2 … 1e-1 are ordinary local minima of the variable-projection landscape (there are many; e.g. the
degree-6 Jacobian family has ≥ 10 distinct minima, the lowest `0.065` being reached by 70/160 restarts). The
searches are heuristic: a zero that exists can be missed, but the calibration bounds the miss probability
per restart (≈ 0.5–0.85 at K = 32), so N independent restarts miss a comparable needle with probability
≈ 0.85^N (N = 100 ⇒ 1e-7). No statement here is a proof of nonexistence; exact exclusions are the business of
the exclusion lane.

## 4. Structured families (task 4)

* Coordinate cones `D ∈ C[subset of variables] ∩ E1-kernel` (vertex a line or a plane; 13 subsets × degrees
  6–12, `w2-chain2`): at `d = 6` every cone family is empty (`K = 0`) except `C[X0,X1,Y0,Y1]` (`K = 3`,
  global minimum rel 0.188 reached by 120/120 restarts) and `C[X1,X2,Y0,Y1]` (`K = 1`, rel 0.631); at `d = 7`
  the non-empty cones are again `C[X0,X1,Y0,Y1]` (`K = 6`, 9 distinct minima, lowest 0.104 reached 23/120)
  and `C[X1,X2,Y0,Y1]` (`K = 4`, lowest 0.445 reached 63/67) — both landscapes fully covered, no zero. The
  3-variable cones (vertex a line) are empty up to d = 7: no cone over a plane curve is E1-admissible there.
  At d = 8 the first admissible 3-variable cone appears: `C[X0,X1,X2] ∩ kernel` is 1-dimensional (a single cone
  over a plane octic with vertex the line {X = 0}); its fixed-D residual is 0.614 — far from a solution.
  Remaining degrees 8–12 pending.
### Structured fixed-D summary (auto-generated)

<!-- AUTO-STRUCT -->
| family | d | candidates | pass E1 | served_mod_J > 0 | serves A |
|---|---|---|---|---|---|
| companion | 6 | 3 | 0 | 0 | 0 |
| cone | 6 | 10 | 0 | 0 | 0 |

No structured candidate serves any target beyond the trivially served Jacobian part.
<!-- /AUTO-STRUCT -->

* Note on provenance: the lane directory is shared with the exclusion lens of this lane; a run of my `w2_struct.py`
  with default arguments (job `w2-struct`, 4386 candidates, d = 6..12) was started on the remote at 01:11 UTC that is
  not one of my chain invocations — its records (same schema, same code) are included in the tables below.
* Fixed-D structured candidates (`w2_struct.py`): companion pieces of the accepted `one-reflection`
  certificate (`f_c`, `h`, `δ`, skew-twisted variants and products), the homogenised binary-cubic covariants
  of `P(U) = G0 U³ + cp U² + skew U + Y6` (degree 12), 3×3 Hessian minors (degree 12), 2×2 minors × invariants,
  block radial derivatives × invariants, products of the skew-adapted invariants at degrees 9–12, norm-type
  cones `f·q + λ g²`, Fermat/twisted 3-variable cones, random generic controls. First results: `f_c`, `f_c + skew`
  and all degree-6 Fermat cones fail E1 (`A D ∉ J_f`). Local exact E1 screen of the 37 degree ≤ 8 companion /
  radial / Hessian-product candidates: the 32 block-radial-derivative products pass (they lie in `J_f`, hence
  inside the Jacobian-type part of the kernel families already searched), while the five irreducible-type
  companion pieces (`h`, `h ± skew·e`, `f_c`, `f_c + skew`) fail E1 exactly — the `one-reflection` companion
  denominators do not even satisfy the necessary condition on the skew sextic.

## 5. Rational reconstruction attempts (task 5)

Finding from the calibration needle: a numerical zero (rel 2e-13) does **not** rationalise directly — the
zero set is positive-dimensional (LM Jacobian at D* has two near-null directions beyond scaling, singular
values 2.7e-10 and 7.9e-6; e.g. direction `X0⁴X1² − X1²Y0⁴`), so a generic numerical zero is an irrational
point of a family of denominators. Remedy implemented (`w2lib.snap_confirm`): drop monomials with
|coefficient| < τ (τ ladder 1e-6 … 1e-1), refit by LM inside the sub-family, rationalise
(`limit_denominator` ladder + LLL), confirm exactly (E1 mod p, fixed-D rref mod p, CRT lift, exact 5-variable
identity, `core.verify_witness` on the 6-variable lift). On the needle: τ = 1e-4 keeps 8 monomials, refit
rel 4e-15, exact solve → the accepted D, verified.

Skew-lane candidate/confirmation events (auto-generated):

<!-- AUTO-CANDS -->
(none: no skew-lane restart reached the 1e-7 confirmation threshold)
<!-- /AUTO-CANDS -->

## 6. Derivation of the reduction (for the record)

Let `v = N/D` with `K_A(v) := A − v(F) + F div v = F² Q/D²`. Average over `g: Y2 ↦ ζ^k Y2` (F, A invariant;
`K_A` is natural): `v̄` is invariant, `K_A(v̄) = F² · (avg of g·(Q/D²))`, denominator `∏ g(D)²` coprime to F.
Invariance: `v̄_i = u_i(x, w)` for i < 5 and `v̄_5 = Y2 a(x, w)`, `w = Y2^6`. With `F = f + w`, regard
everything as rational in `(x, F)`: `K_A = A − ũ(f) − 6(F − f)ã + F(div_x ũ + ∂_F ũ(f) + ã + 6(F − f)ã_F)`.
`F² | K_A` ⇒ coefficients of `F⁰, F¹` vanish: `A − u₀(f) + 6 f a₀ = 0` and `div u₀ − 5 a₀ = 0`
(`u₀ = ũ(x,0)`), hence (*) for `u₀`. Conversely `v = (u, Y2 div u/5)` gives `K_A(v) = A − u(f) + (6/5) f div u = 0`
with denominator `D²` coprime to `F` (F is monic in Y2, D is Y2-free).
