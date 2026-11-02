# Wave 2 — `four-plus-two-residual` — JOINT SEARCH lens (sector-restricted variable projection)

STATUS: **NO WITNESS FOUND** (search in progress, remote jobs running). No candidate with residual < 1e-8 in any family/sector of this instance so far. Best relative residuals of the joint (D, N, Q) variable-projection search: d=6 kernel+J sectors 5.5e-4 / 1.5e-3 (775 restarts each, always the same minimum); d=7 kernel+J sectors 5.0e-4 / 1.4e-3 (75 restarts each so far); all pencils (K=2) scanned exhaustively (best 5e-2 .. 2e-1). **Pipeline validated end-to-end on the positive control** `family-hesse-web` (same target X0^3Y0^3, known d=6 witness): random-restart L-BFGS finds the witness in ~6% of restarts (6/112), LM in ~6% (5/80), the LM polish brings the residual to 3e-30 and the rational reconstruction + exact sector-restricted re-solve recovers a verified witness (5/5 hits) — so a witness of comparable family dimension (dimP ~ 12) would have been found here with high probability in the d=6 sectors (775 restarts each, dimP 21/29); the hit rate for dimP 20-50 is being calibrated (coarser-subgroup controls). Degrees <= 5 are excluded over C by the sibling exclusion lens.

Log: `/tmp/hodge/results/four-plus-two-residual/wave2-search.log.jsonl` (one JSON line per batch; merged from the remote per-job files `wave2-search/*.jsonl` by `wave2-search/fpt_sync.py`). Code: `/tmp/hodge/results/four-plus-two-residual/wave2-search/` (`fptlib.py` library, `fpt_search.py` driver, `fpt_dims.py`, `fpt_control.py`). Engine v2 (post sign-fix). Sibling lens: `wave2-four-plus-two-residual-excl.md` (exact exclusions; deg D <= 5 excluded over C).

## 0. Setup (what is searched and why it is lossless)
* F = H(X0,X1,Y0,Y1) + B(X2,Y2), A = X0^3 Y0^3, q = 1. Projective diagonal group G of F: order 24, generated (weights mod 6) by h' = (5,5,0,5,5,0) (zeta6-scaling of (X2,Y2) up to scalars), r2 = (0,0,3,0,0,0), s = (0,0,0,3,3,0). Sector label of a monomial: (m.h', m.r2, m.s) mod 6 in Z/6 x {0,3} x {0,3}; the Y-parity is determined by the label.
* For D a G-eigenvector (every D supported in one sector), the fixed-D unknowns N_i, Q can be restricted to the characters forced by (A, D) — lossless by averaging over G — which cuts rows/columns by ~24x (e.g. d=6: 392 unknowns instead of 2602; d=7: 472-506 instead of 4497).
* Family per (d, sector): the full E1-admissible space {D in S_d^sector : A D in J(F)} (= kernel of [A] on the Jacobian ring + J(F)_d), computed numerically (SVD, singular-value gap >= 1e-2 vs <= 1e-14) and cross-checked mod two primes; F-multiples (spurious denominators) are projected out of the parameter space.
* Joint search: f(a) = min_{N,Q} |M(a)x - b(a)|^2 on the unit sphere of family coordinates a (variable projection; Golub-Pereyra gradient, verified against finite differences to 1e-10), random-restart L-BFGS-B, Levenberg-Marquardt with the Kaufman Jacobian, deflation (repelling from found minima), sparse-integer initialisation. Candidates with rel < 1e-8 are polished, rationally reconstructed with growing denominators (1..1e6) and re-solved exactly (sector-restricted fixed-D system with CRT lifting, `core.verify_witness`, F-divisibility check).

### 0.1 Koszul (E2-only) formulation — validated, now the production formulation
Every E1 lift is N = N0(D) + sum_{i<j} c_ij (F_j d_i - F_i d_j) (book Ch.5 §5.2; the partials form a regular sequence), so with a linear E1 lift N0(a) (least squares on the sector E1 system, residual ~1e-14) the joint problem becomes the E2 identity D div N - N(D) = F Q in degree 2d with unknowns c_ij (degree d-4) and Q (degree 2d-6): on the control the sector system shrinks from 144x56 to 39x17 (witness residual 5.6e-30, gradient check 2e-10, LM finds the witness), on this instance at d=6/(0,0,0) from 1836x392 to 384x88 — ~10x faster restarts at d=6, ~100x cheaper linear algebra at d<=7. Jobs kmain (kernel+J d=6,7,8; LM; 2000/2000/1000 restarts per sector), kstruct (variable-subset cones, J-only, binary-factor products, Hessian relaxation), kctl (hit rate vs family dimension on the control) replaced the full-formulation jobs at 20:20.

## 1. E1-admissible denominator spaces (dimensions) — DONE (job fptdims; `wave2-search/dims.jsonl`)
Per degree, summed over the 24 character sectors (mod p = 2^31-1 and 2^31-19 agree with the numerical SVD in every sector):

| d | sectors (nonempty) | #monomials | dim E1-admissible | dim J(F)_d | kernel mod J (even,odd) | F-multiples | largest sector K | numerical/mod-p agree |
|---|---|---|---|---|---|---|---|---|
| 4 | 16 | 126 | 0 | 0 | 0, 0 | 0 | 0 | yes |
| 5 | 20 | 252 | 6 | 6 | 0, 0 | 0 | 2 | yes |
| 6 | 23 | 462 | 69 | 36 | 20, 13 | 1 | 30 | yes |
| 7 | 24 | 792 | 252 | 126 | 63, 63 | 6 | 54 | yes |
| 8 | 24 | 1287 | 630 | 336 | 154, 140 | 21 | 93 | yes |
| 9 | 24 | 2002 | 1294 | 756 | 269, 269 | 56 | 128 | yes |
| 10 | 24 | 3003 | 2346 | 1497 | 434, 415 | 126 | 186 | yes |
| 11 | 24 | 4368 | 3828 | 2682 | 573, 573 | 252 | 256 | yes |
| 12 | 24 | 6188 | 5795 | 4437 | 690, 668 | 462 | 382 | yes |

Core kernels ker(a : R(H)_k -> R(H)_{k+6}) mod J(H)_k by (Y0,Y1)-parity (even, odd), two primes agree:
k=5: (0,0) | k=6: (20,13) | k=7: (30,30) | k=8: (41,34) | k=9: (38,38) | k=10: (36,31) | k=11: (26,26) | k=12: (19,16)
Thom-Sebastiani check: the kernel-mod-J totals above equal sum_j K_core(d-j) x dim R(B)_j (binary standard monomials X2^a Y2^b, a,b <= 4) for every d and parity (e.g. d=12: 690, 668; the book's forced lower bound at d=12 is 677).

Per-sector detail for d = 6..8 (nonempty sectors only):

| d | sector (h,r2,s) | Y-par | #mons | dim fam | dim J | kernel mod J | F-mult |
|---|---|---|---|---|---|---|---|
| 6 | (0, 0, 0) | 0 | 48 | 30 | 10 | 20 | 1 |
| 6 | (0, 0, 3) | 1 | 40 | 21 | 8 | 13 | 0 |
| 6 | (0, 3, 0) | 1 | 3 | 2 | 2 | 0 | 0 |
| 6 | (1, 0, 0) | 1 | 28 | 2 | 2 | 0 | 0 |
| 6 | (1, 0, 3) | 0 | 28 | 2 | 2 | 0 | 0 |
| 6 | (1, 3, 0) | 0 | 28 | 2 | 2 | 0 | 0 |
| 6 | (1, 3, 3) | 1 | 28 | 2 | 2 | 0 | 0 |
| 6 | (5, 0, 0) | 1 | 6 | 2 | 2 | 0 | 0 |
| 6 | (5, 0, 3) | 0 | 6 | 2 | 2 | 0 | 0 |
| 6 | (5, 3, 0) | 0 | 6 | 2 | 2 | 0 | 0 |
| 6 | (5, 3, 3) | 1 | 6 | 2 | 2 | 0 | 0 |
| 7 | (0, 0, 0) | 1 | 48 | 31 | 11 | 20 | 1 |
| 7 | (0, 0, 3) | 0 | 40 | 21 | 8 | 13 | 0 |
| 7 | (0, 3, 0) | 0 | 48 | 31 | 11 | 20 | 1 |
| 7 | (0, 3, 3) | 1 | 40 | 21 | 8 | 13 | 0 |
| 7 | (1, 0, 0) | 0 | 56 | 4 | 4 | 0 | 0 |
| 7 | (1, 0, 3) | 1 | 56 | 4 | 4 | 0 | 0 |
| 7 | (1, 3, 0) | 1 | 28 | 2 | 2 | 0 | 0 |
| 7 | (1, 3, 3) | 0 | 28 | 2 | 2 | 0 | 0 |
| 7 | (4, 0, 0) | 1 | 18 | 6 | 6 | 0 | 0 |
| 7 | (4, 0, 3) | 0 | 12 | 4 | 4 | 0 | 0 |
| 7 | (4, 3, 0) | 0 | 18 | 6 | 6 | 0 | 0 |
| 7 | (4, 3, 3) | 1 | 12 | 4 | 4 | 0 | 0 |
| 7 | (5, 0, 0) | 0 | 68 | 54 | 24 | 30 | 2 |
| 7 | (5, 0, 3) | 1 | 68 | 54 | 24 | 30 | 2 |
| 7 | (5, 3, 0) | 1 | 6 | 4 | 4 | 0 | 0 |
| 7 | (5, 3, 3) | 0 | 6 | 4 | 4 | 0 | 0 |
| 8 | (0, 0, 0) | 0 | 93 | 60 | 20 | 40 | 2 |
| 8 | (0, 0, 3) | 1 | 80 | 42 | 16 | 26 | 0 |
| 8 | (0, 3, 0) | 1 | 48 | 32 | 12 | 20 | 1 |
| 8 | (0, 3, 3) | 0 | 40 | 21 | 8 | 13 | 0 |
| 8 | (1, 0, 0) | 1 | 56 | 4 | 4 | 0 | 0 |
| 8 | (1, 0, 3) | 0 | 56 | 4 | 4 | 0 | 0 |
| 8 | (1, 3, 0) | 0 | 56 | 4 | 4 | 0 | 0 |
| 8 | (1, 3, 3) | 1 | 56 | 4 | 4 | 0 | 0 |
| 8 | (3, 0, 0) | 1 | 30 | 10 | 10 | 0 | 0 |
| 8 | (3, 0, 3) | 0 | 30 | 10 | 10 | 0 | 0 |
| 8 | (3, 3, 0) | 0 | 30 | 10 | 10 | 0 | 0 |
| 8 | (3, 3, 3) | 1 | 30 | 10 | 10 | 0 | 0 |
| 8 | (4, 0, 0) | 0 | 109 | 93 | 52 | 41 | 6 |
| 8 | (4, 0, 3) | 1 | 96 | 82 | 48 | 34 | 4 |
| 8 | (4, 3, 0) | 1 | 18 | 12 | 12 | 0 | 0 |
| 8 | (4, 3, 3) | 0 | 12 | 8 | 8 | 0 | 0 |
| 8 | (5, 0, 0) | 1 | 68 | 56 | 26 | 30 | 2 |
| 8 | (5, 0, 3) | 0 | 68 | 56 | 26 | 30 | 2 |
| 8 | (5, 3, 0) | 0 | 68 | 56 | 26 | 30 | 2 |
| 8 | (5, 3, 3) | 1 | 68 | 56 | 26 | 30 | 2 |

## 2. Joint searches run so far (per family / degree / sector)
rel = min_x |M(a)x-b(a)|^2 / |b(a)|^2 at the optimizer's end point (b = -A D^2), F-multiples projected out; "grid" = exhaustive scan for pencils (K<=2); cands = restarts with rel < 1e-8 (then polished / reconstructed / re-solved exactly).

| family | d | sector | K | dimP | rows x cols (restricted) | mode | restarts | best rel | cands<1e-8 | CPU s |
|---|---|---|---|---|---|---|---|---|---|---|
| kernelJ | 5 | (0, 0, 0) | 1 | 1 | 888x64 | grid | 1 | 2.04e-01 | 0 | 0 |
| kernelJ | 5 | (0, 3, 0) | 1 | 1 | 888x64 | grid | 1 | 1.87e-01 | 0 | 0 |
| kernelJ | 5 | (1, 0, 0) | 2 | 2 | 1152x248 | grid | 20 | 7.05e-02 | 0 | 19 |
| kernelJ | 5 | (1, 0, 3) | 2 | 2 | 1152x248 | grid | 20 | 5.29e-02 | 0 | 20 |
| kernelJ | 6 | (0, 0, 0) | 30 | 29 | 1836x392 | lbfgs | 775 | 5.47e-04 | 0 | 11038 |
| kernelJ | 6 | (0, 0, 3) | 21 | 21 | 1836x408 | lbfgs | 775 | 1.49e-03 | 0 | 7932 |
| kernelJ | 6 | (0, 3, 0) | 2 | 2 | 1788x144 | grid | 20 | 1.17e-01 | 0 | 17 |
| kernelJ | 6 | (1, 0, 0) | 2 | 2 | 1788x292 | grid | 20 | 6.93e-02 | 0 | 35 |
| kernelJ | 6 | (1, 0, 3) | 2 | 2 | 1788x292 | grid | 20 | 5.23e-02 | 0 | 33 |
| kernelJ | 6 | (1, 3, 0) | 2 | 2 | 1788x292 | grid | 20 | 6.90e-02 | 0 | 33 |
| kernelJ | 6 | (1, 3, 3) | 2 | 2 | 1788x292 | grid | 20 | 5.21e-02 | 0 | 38 |
| kernelJ | 6 | (5, 0, 0) | 2 | 2 | 1524x146 | grid | 20 | 1.54e-01 | 0 | 13 |
| kernelJ | 6 | (5, 0, 3) | 2 | 2 | 1524x146 | grid | 20 | 1.45e-01 | 0 | 12 |
| kernelJ | 6 | (5, 3, 0) | 2 | 2 | 1524x146 | grid | 20 | 1.41e-01 | 0 | 12 |
| kernelJ | 6 | (5, 3, 3) | 2 | 2 | 1524x146 | grid | 20 | 1.33e-01 | 0 | 16 |
| kernelJ | 7 | (0, 0, 0) | 31 | 30 | 2760x472 | lbfgs | 75 | 5.07e-04 | 0 | 1516 |
| kernelJ | 7 | (0, 0, 3) | 21 | 21 | 2760x493 | lbfgs | 75 | 1.38e-03 | 0 | 1215 |
| kernelJ | 7 | (0, 3, 0) | 31 | 30 | 2760x472 | lbfgs | 75 | 5.01e-04 | 0 | 1831 |
| kernelJ | 7 | (0, 3, 3) | 21 | 21 | 2760x493 | lbfgs | 75 | 1.37e-03 | 0 | 1168 |
| kernelJ | 7 | (1, 0, 0) | 4 | 4 | 2424x506 | lbfgs | 75 | 5.95e-02 | 0 | 396 |
| kernelJ | 7 | (1, 0, 3) | 4 | 4 | 2424x506 | lbfgs | 75 | 3.82e-02 | 0 | 424 |
| kernelJ | 7 | (1, 3, 0) | 2 | 2 | 2424x336 | grid | 20 | 6.86e-02 | 0 | 74 |
| kernelJ | 7 | (1, 3, 3) | 2 | 2 | 2424x336 | grid | 20 | 5.19e-02 | 0 | 71 |
| kernelJ | 7 | (4, 0, 0) | 6 | 6 | 2424x276 | lbfgs | 75 | 3.29e-02 | 0 | 184 |
| kernelJ | 7 | (4, 0, 3) | 4 | 4 | 2424x295 | lbfgs | 75 | 8.48e-02 | 0 | 150 |
| kernelJ | 7 | (4, 3, 0) | 6 | 6 | 2424x276 | lbfgs | 75 | 3.19e-02 | 0 | 133 |
| kernelJ | 7 | (4, 3, 3) | 4 | 4 | 2424x295 | lbfgs | 75 | 7.83e-02 | 0 | 166 |
| kernelJ | 7 | (5, 3, 0) | 4 | 4 | 2824x292 | lbfgs | 75 | 8.78e-02 | 0 | 171 |
| kernelJ | 7 | (5, 3, 3) | 4 | 4 | 2824x292 | lbfgs | 75 | 8.17e-02 | 0 | 197 |

## 3. Positive control (family-hesse-web, A = X0^3Y0^3, d = 6 witness) — pipeline validated
Stabilizer of (A, D_witness): 864 elements; E1-admissible family in the witness sector: K = 13 (mod-p dim 13; J-part 6, F-multiples 1); restricted full system 144 x 56. Witness: projection residual 9e-16, rel residual 3.0e-29 (formulation exact). Witness sits at a rank drop of M(a) (rank 46 vs 48 generic) — but this is not a general signature (2 of 6 accepted witnesses tested show a drop). Basin: perturbations of angular size <= 0.3 polish back to the witness (6/6), 0.5: 2/6 — in dimP = 12 such a basin has measure ~1e-5 of the sphere for *uniform* sampling, yet the optimizers reach it from ~6% of Gaussian starts (the landscape funnels: along the straight line from a random point to the witness the residual decreases monotonically).
* L-BFGS-B (400 it): 6 hits / 112 restarts (rel 4e-15 .. 6e-16 before polish); residual clusters at 7.5e-6 (a dominant local minimum, ~25%) and 1e-3 .. 1e-2.
* LM (Kaufman Jacobian, 200 nfev): 5 hits / ~80 restarts, polished to 3e-30; **5/5 reconstructed and re-solved exactly** (`family-hesse-web` witness recovered up to the F-gauge D -> D + lambda F, which preserves q = 1 certificates with the same N).
* Earlier control runs with 0/148 hits were terminated by a foreign `pkill` matching the generic script names (fixed by unique names); their data are discarded.
* Pending: deflation and sparse-init variants; coarser-subgroup controls (dimP 20-60) to calibrate hit rate vs dimension.

## 4. Interpretation so far
Low residual floors (5e-4 .. 1.5e-3) are reproduced across hundreds of restarts in each d=6,7 sector, so the dominant local minima of these landscapes are known and none is a certificate; but the control shows that random-restart variable projection can miss an existing witness in a 13-dim family, so these searches are evidence, not exclusion. Exclusions over C are the sibling lens' job.
