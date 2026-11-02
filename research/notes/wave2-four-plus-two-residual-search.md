# Wave 2 — `four-plus-two-residual` — JOINT SEARCH lens (sector-restricted variable projection)

STATUS: **NO WITNESS FOUND** (search in progress, remote jobs running). No candidate with residual < 1e-8 in any family/sector searched so far. Best relative residuals of the joint (D, N, Q) variable-projection search: d=6 kernel+J sectors 5.5e-4 / 1.5e-3 (hundreds of restarts each, always the same minimum); d=7 kernel+J sectors 5.0e-4 / 1.4e-3 (first batches). Positive-control finding: on `family-hesse-web` (same target X0^3Y0^3, known d=6 witness) the known witness has residual 3e-29 inside the E1-admissible family but random-restart L-BFGS / LM never reached it in the first ~100 restarts (lands on a 7.5e-6 local minimum) — calibration jobs for four optimizer strategies are running; see §3.

Log: `/tmp/hodge/results/four-plus-two-residual/wave2-search.log.jsonl` (one JSON line per batch; merged from the remote per-job files `wave2-search/*.jsonl` by `wave2-search/fpt_sync.py`). Code: `/tmp/hodge/results/four-plus-two-residual/wave2-search/` (`fptlib.py` library, `fpt_search.py` driver, `fpt_dims.py`, `fpt_control.py`). Engine v2 (post sign-fix). Sibling lens: `wave2-four-plus-two-residual-excl.md` (exact exclusions; deg D <= 5 excluded over C).

## 0. Setup (what is searched and why it is lossless)
* F = H(X0,X1,Y0,Y1) + B(X2,Y2), A = X0^3 Y0^3, q = 1. Projective diagonal group G of F: order 24, generated (weights mod 6) by h' = (5,5,0,5,5,0) (zeta6-scaling of (X2,Y2) up to scalars), r2 = (0,0,3,0,0,0), s = (0,0,0,3,3,0). Sector label of a monomial: (m.h', m.r2, m.s) mod 6 in Z/6 x {0,3} x {0,3}; the Y-parity is determined by the label.
* For D a G-eigenvector (every D supported in one sector), the fixed-D unknowns N_i, Q can be restricted to the characters forced by (A, D) — lossless by averaging over G — which cuts rows/columns by ~24x (e.g. d=6: 392 unknowns instead of 2602; d=7: 472-506 instead of 4497).
* Family per (d, sector): the full E1-admissible space {D in S_d^sector : A D in J(F)} (= kernel of [A] on the Jacobian ring + J(F)_d), computed numerically (SVD, singular-value gap >= 1e-2 vs <= 1e-14) and cross-checked mod two primes; F-multiples (spurious denominators) are projected out of the parameter space.
* Joint search: f(a) = min_{N,Q} |M(a)x - b(a)|^2 on the unit sphere of family coordinates a (variable projection; Golub-Pereyra gradient, verified against finite differences to 1e-10), random-restart L-BFGS-B, Levenberg-Marquardt with the Kaufman Jacobian, deflation (repelling from found minima), sparse-integer initialisation. Candidates with rel < 1e-8 are polished, rationally reconstructed with growing denominators (1..1e6) and re-solved exactly (sector-restricted fixed-D system with CRT lifting, `core.verify_witness`, F-divisibility check).

## 1. E1-admissible denominator spaces (dimensions)
(filled from `wave2-search/dims.jsonl` when the job finishes; the sibling exclusion lens has the exact core kernels: k<=5: 0,0 | 6: 20,13 | 7: 30,30 | 8: 41,34 | 9: 38,38.)
Tensor prediction (Thom-Sebastiani) of dim ker[A] mod J on S_d by Y-parity (even, odd): d=6: 20, 13 | d=7: 63, 63 | ... | d=12: 690, 668 (consistent with the book's forced bound >= 677 at d=12).

## 2. Joint searches run so far (per family / degree / sector)
(see the table in the log summary below; updated after every pull)

## 3. Positive control (family-hesse-web, A = X0^3Y0^3, d = 6 witness)
Stabilizer of (A, D_witness): 864 elements; E1-admissible family in the witness sector: K = 13 (mod-p dim 13; J-part 6, F-multiples 1); restricted system 144 x 56. Witness: projection residual 9e-16, rel residual 3.0e-29 (the formulation is exact). Random restarts: 0/48 hits (L-BFGS, 400 it): residual clusters at 7.5e-6 (25%) and 1e-3..1e-2. Calibration jobs: lbfgs/lm/deflate/sparse x 800-1600 restarts (running).

## 4. Interpretation so far
Low residual floors (5e-4 .. 1.5e-3) are reproduced across hundreds of restarts in each d=6,7 sector, so the dominant local minima of these landscapes are known and none is a certificate; but the control shows that random-restart variable projection can miss an existing witness in a 13-dim family, so these searches are evidence, not exclusion. Exclusions over C are the sibling lens' job.
