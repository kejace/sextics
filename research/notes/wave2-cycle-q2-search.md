# Wave 2 search report: `cycle-q2` (q = 2, A = X0^2 X1^2 X2^3 Y0^2 Y1^2 Y2)

STATUS: no witness found so far. Search infrastructure built and validated (2026-11-01 evening); remote jobs queued (the
128-core box is at load > 200 from other lanes, submissions wait for capacity). This file is rewritten after every batch.

Code (lane dir `/tmp/hodge/results/cycle-q2/`): `w2lib.py` (sector bookkeeping, E1-kernel per sector mod p and numeric,
Koszul-reduced sector-restricted complex variable projection `KoszulVarPro` with L-BFGS and Levenberg-Marquardt drivers,
rational reconstruction, exact confirmation through the wave-1 sector solver `eq.py`), `w2_search.py` (multiprocessing
driver over sectors; families kernel/jac/cone4/cone3/mult4/mult5/sing2/full), `w2_e1dims.py` (step 1), `w2_hess.py`
(Hessian relaxation at d = 5 exactly, J-type relaxation at d = 6). Log: `wave2-search.log.jsonl` (one line per batch).

## Method

E1 (`2 N(F) = A D`) is solved exactly on the family: for D in the sector E1-kernel K_kappa the lift N0(D) is linear in D.
Every lift is N0 + K(c) with the Koszul fields K(c) = sum_{i<j} c_ij (F_j d_i - F_i d_j), deg c_ij = d + 2, and E2 becomes
the affine problem  -B_D(N0) in F^2 S_{2d-6} + im T_D  (docs/05 §5.2).  For D(a) = sum a_k D_k the residual
f(a) = min_{c,Q} |T_{D(a)}(c) - F^2 Q + B_{D(a)}(N0(a))|^2 is minimized over the complex unit sphere (complex solutions are
not missed) with the variable-projection gradient (L-BFGS) or the Kaufman Jacobian (Levenberg-Marquardt).  All unknowns are
restricted to the character sector forced by kappa (lossless by averaging; 72 characters of the order-72 projective
diagonal group, parity included).  System sizes (rows x unknowns per sector): d=6 490x280, d=7 740x435, d=8 1120x670,
d=9 1650x1000, d=10 2360x1450, d=12 4500x2700.  Candidates with rel < 1e-10 are reconstructed over Q (common
denominators up to 2^20, continued fractions) and re-solved exactly (E1, full sector system mod 2^31-1 and a second prime,
CRT lift, `core.verify_witness`, F | D check).

Validation: gradient vs finite differences 5e-10 (cycle-q2 sector 96) and 8e-10 (Fermat); block-q2's witness has
f = 8e-32; fermat-q2's three strict witnesses (subset Euler, I = {X1,X2}, {X1,Y1}, {X2,Y1}) are exact zeros of f in the
6-dimensional pure-power family, and random-restart hit rate is 15-17 % per restart (L-BFGS 15/100, LM 17/100) -- the
landscape has spurious local minima, so dozens of restarts per sector are needed; LM reaches f ~ 1e-64 at true zeros.

## Results so far

### Theorem (exact, this wave): every certificate for cycle-q2 has deg D >= 6, for arbitrary D

E1 forces A D in J(F); the annihilator of [A] in R_d vanishes for d <= 5 (dossier A.3 at p = 10007; wave-1 sector
kernels at p = 2^31-1 give K_5 = J_5), so d <= 4 is impossible and d = 5 forces D = sum_i w_i F_i with constant w.
With the canonical lift N0 = (A/2) w, E2 modulo J(F) reads sum_{i,j} w_i w_j [A F_ij] = 0 in R_16 (im T_D, F^2 Q and
D div N0 all lie in J).  The 21 classes [A F_ij] (i <= j) fall into 21 distinct character sectors of R_16 and each is
nonzero mod 2^31-1 and mod 2^31-19 (`w2_hess.py d5`: rank 21/21); dim R_16 mod p = 951 for both primes, the
characteristic-zero value, so J_16 mod p has full rank and the independence transfers to Q and C.  Hence all w_i w_j = 0,
w = 0, D = 0.  No symmetry or parity assumption on D is needed.  (Closes gap 3 of the dossier; the previous bound was
deg D >= 5 with only semi-invariant degree-5 denominators excluded.)

### Numerical joint search

(pending; remote jobs queued: e1dims d=4..12; kernel d=6 (300 restarts/sector) and d=7 (150); kernel d=8 (60);
structured cone4/mult4/cone3 families d=7..10; J-type relaxation d=6)
