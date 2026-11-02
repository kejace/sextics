# Wave 2 search report: klein-coupled-character (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (search in progress)

Updated 2026-11-01 20:28:25 from 68 log lines in `/tmp/hodge/results/klein-coupled-character/wave2-search.log.jsonl`.


Code (lane dir `/tmp/hodge/results/klein-coupled-character/`): `w2lib.py` (sector-restricted Koszul variable projection,
sector E1 spaces with rational lifts, mod-J relaxation), drivers `w2_kernels.py`, `w2_hess.py`, `w2_varpro.py`,
`w2_struct.py`, control `w2_control.py`, queue `w2_queue.py` (remote jobs `w2-*`), this generator `w2_report.py`.
Conventions: (x,y,z,t,r,s) = (X0..Y2); sigma-weights W = (1,2,4,0,0,0) mod 7; "sector (d, parD, k)" = Y-parity parD
denominators of degree d and sigma-weight k; A has weight 3, Y-odd.  "rel" = squared residual / squared right-hand side.

## Method
* **Sector E1 spaces.** (J(F):A)_d per sector = J(F)_d ⊕ (kernel of [A] on (R/J)_d) computed mod p = 2147483647 with
  canonical bases, lifted to Q by 12-prime CRT + LLL (small-height representatives), exact E1 lifts N_k by CRT.
* **Joint E2 search (varpro, Koszul form).** D(a) = sum a_k D_k, N0(a) = sum a_k N_k, unknowns p_ij (deg d-4) and Q
  (deg 2d-6) character-restricted; Q eliminated by the t-monic remainder mod F.  f(a) = min_p |R(T_{D(a)} p + B_{D(a)}(N0(a)))|^2
  on the unit sphere orthogonal to the F-multiples.  Optimizers: L-BFGS (engine style), LM/TRF with the Kaufman
  Jacobian, joint LM in (a,x), combined (E2 + relaxation) LM.  Candidates with rel < 1e-10 are reconstructed in the
  canonical rem_F chart with denominators up to 1e6 and re-solved exactly (kcc2 Koszul test mod p, CRT lift, verify_witness).
* **Mod-J (Hessian) relaxation.** Necessary condition B_{D(a)}(N0(a)) in J(F)_{2d} (for Jacobian-type D: A (w^T Hess F w) in J).
  Real zeros searched numerically (L-BFGS, hundreds-thousands of restarts; dimension <= 45 for d <= 8 so the global
  minimum is reliably found); the exact quadric system mod 32003 in a rational chart complementary to the F-multiples
  is passed to Singular: `dim = 0` (only the origin) is a **rigorous exclusion** of every sigma-eigen denominator of that
  (d, parity, weight) -- mod-p emptiness of a homogeneous system implies emptiness over Q (closure of the generic fibre).
  `dim > 0` mod p does not imply real/rational zeros (the (6,even,0) component is an octic curve with no real points found).


## (1) E1-admissible sector dimensions (kernel of [A] on (R/J)_d^{parity}[k], mod p; upper bounds over Q)

| d | parity | nD per k=0..6 | J_d per k | kernel mod J per k (nontrivial) | rational reps (max height) | lifts failed |
|---|---|---|---|---|---|---|
| 9 | even | -,156,156,139,-,139,139 | -,63,63,52,-,52,52 | -,3,3,2,-,1,3 | -,-,-,49994,-,31345,272340 | 0 |
| 9 | odd | -,139,-,-,-,136,- | -,49,-,-,-,46,- | -,5,-,-,-,6,- | -,2420,-,-,-,3625,- | 0 |

## (2) Mod-J relaxation per sector (real zeros numerically; Singular exact dimension mod 32003)

| d | parity | k | family | unknowns (dimP) | equations | restarts | best rel | median rel | #rel<1e-10 | Singular dim (vars, eqs, s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | even | 0 | J | 11 | 123 | 100 | 8.21e-03 | 1.53e-02 | 0 | 2 (11, 46, 0.0s) |
| 6 | even | 1 | J | 1 | 123 | 100 | 4.11e-01 | 4.11e-01 | 0 | **0** (1, 2, 0.0s) |
| 6 | even | 2 | J | 1 | 128 | 100 | 8.56e-01 | 8.56e-01 | 0 | **0** (1, 2, 0.0s) |
| 6 | even | 3 | J | 1 | 123 | 100 | 5.79e-01 | 5.79e-01 | 0 | **0** (1, 4, 0.0s) |
| 6 | even | 4 | J | 1 | 123 | 100 | 6.44e-01 | 6.44e-01 | 0 | **0** (1, 2, 0.0s) |
| 6 | even | 5 | J | 1 | 123 | 100 | 4.02e-01 | 4.02e-01 | 0 | **0** (1, 2, 0.0s) |
| 6 | even | 6 | J | 1 | 123 | 100 | 1.25e-01 | 1.25e-01 | 0 | **0** (1, 2, 0.0s) |
| 6 | odd | 0 | J | None | None | None | - | - | None | - |
| 6 | odd | 1 | J | 3 | 123 | 100 | 7.88e-02 | 7.89e-02 | 0 | **0** (3, 8, 0.0s) |
| 6 | odd | 2 | J | 3 | 128 | 100 | 6.76e-02 | 6.76e-02 | 0 | **0** (3, 8, 0.0s) |
| 6 | odd | 3 | J | 3 | 123 | 100 | 5.27e-03 | 3.20e-01 | 0 | **0** (3, 22, 0.0s) |
| 6 | odd | 4 | J | 3 | 123 | 100 | 1.51e-01 | 1.51e-01 | 0 | **0** (3, 12, 0.0s) |
| 6 | odd | 5 | J | 3 | 123 | 100 | 1.27e-02 | 4.59e-01 | 0 | **0** (3, 20, 0.0s) |
| 6 | odd | 6 | J | 3 | 123 | 100 | 3.62e-02 | 3.62e-02 | 0 | **0** (3, 14, 0.0s) |
| 7 | even | 0 | J | 3 | 108 | 200 | 1.05e-02 | 1.05e-02 | 0 | **0** (3, 35, 0.0s) |
| 7 | even | 1 | J | 11 | 108 | 200 | 1.97e-03 | 2.21e-03 | 0 | 2 (11, 46, 0.0s) |
| 7 | even | 2 | J | 11 | 105 | 200 | 1.28e-03 | 1.31e-03 | 0 | 2 (11, 86, 0.0s) |
| 7 | even | 3 | J | 8 | 105 | 200 | 1.30e-03 | 2.55e-03 | 0 | **0** (8, 95, 0.0s) |
| 7 | even | 4 | J | 11 | 105 | 200 | 2.67e-03 | 3.25e-03 | 0 | 2 (11, 91, 0.0s) |
| 7 | even | 5 | J | 8 | 108 | 200 | 2.43e-03 | 2.43e-03 | 0 | **0** (8, 96, 0.0s) |
| 7 | even | 6 | J | 8 | 105 | 200 | 6.66e-03 | 1.18e-02 | 0 | **0** (8, 93, 0.0s) |
| 7 | odd | 0 | J | 24 | 108 | 200 | 7.71e-06 | 1.14e-04 | 0 | timeout (24, 93, 45s) |
| 7 | odd | 1 | J | 6 | 108 | 200 | 3.53e-03 | 3.55e-03 | 0 | **0** (6, 78, 0.0s) |
| 7 | odd | 2 | J | 6 | 105 | 200 | 2.00e-03 | 2.06e-03 | 0 | **0** (6, 73, 0.0s) |
| 7 | odd | 3 | J | 6 | 105 | 200 | 3.55e-03 | 4.55e-03 | 0 | **0** (6, 66, 0.0s) |
| 7 | odd | 4 | J | 6 | 105 | 200 | 8.50e-04 | 2.46e-03 | 0 | **0** (6, 68, 0.0s) |
| 7 | odd | 5 | J | 6 | 108 | 200 | 2.34e-04 | 2.35e-04 | 0 | **0** (6, 67, 0.0s) |
| 7 | odd | 6 | J | 6 | 105 | 200 | 1.71e-03 | 1.72e-03 | 0 | **0** (6, 45, 0.0s) |
| 8 | even | 0 | J | 45 | 72 | 120 | 2.14e-08 | 1.80e-07 | 0 | - |
| 8 | even | 0 | probe-J | 45 | 72 | 60 | 2.20e-08 | 1.56e-07 | None | - |
| 8 | even | 1 | probe-J | 19 | 72 | 60 | 2.45e-05 | 3.66e-05 | None | - |
| 8 | even | 2 | probe-J | 19 | 57 | 60 | 4.13e-06 | 1.90e-05 | None | - |
| 8 | even | 3 | probe-J | 18 | 65 | 60 | 1.94e-05 | 2.54e-05 | None | - |
| 8 | even | 4 | probe-J | 19 | 65 | 60 | 1.58e-05 | 2.72e-05 | None | - |
| 8 | even | 5 | probe-J | 18 | 72 | 60 | 7.54e-06 | 2.89e-05 | None | - |
| 8 | even | 6 | probe-J | 18 | 65 | 60 | 1.11e-05 | 2.87e-05 | None | - |
| 8 | odd | 0 | probe-J | 12 | 72 | 40 | 2.95e-05 | 3.36e-05 | None | - |
| 8 | odd | 1 | probe-J | 27 | 72 | 40 | 1.53e-06 | 9.44e-06 | None | - |
| 8 | odd | 2 | probe-J | 27 | 57 | 40 | 2.58e-06 | 7.34e-06 | None | - |
| 8 | odd | 3 | probe-J | 22 | 65 | 40 | 3.15e-06 | 9.63e-06 | None | - |
| 8 | odd | 4 | probe-J | 27 | 65 | 40 | 6.39e-07 | 4.70e-06 | None | - |
| 8 | odd | 5 | probe-J | 22 | 72 | 40 | 6.02e-06 | 2.81e-05 | None | - |
| 8 | odd | 6 | probe-J | 22 | 65 | 40 | 7.47e-06 | 1.32e-05 | None | - |

Empty sectors (no E1-admissible eigen denominator): (6,odd,0)


## (3) Full-E2 joint searches (sector-restricted Koszul varpro)

| d | parity | k | family | K | dimP | rows x cols (p,Q elim.) | restarts | best rel | q25 | median | #<1e-6 | #<1e-10 | methods | candidates |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6 | even | 0 | J | 12 | 11 | 368x30 | 6 | 3.96e-01 | 3.96e-01 | 1.00e+00 | 0 | 0 | lm | 0 |
| 6 | even | 3 | J | 1 | 1 | 365x70 | 2 | 9.56e-01 | 9.56e-01 | 9.56e-01 | 0 | 0 | lm | 0 |
| 9 | even | 3 | K | 2 | 2 | 1969x240 | 4 | 3.41e-02 | 3.41e-02 | 3.41e-02 | 0 | 0 | lm | 0 |

## Remote job queue (tail)
```
20:02:01 gate closed for w2-control: REMOTE BUSY: load 1093 (max 150), available memory 76 GB (min 60), 89 jobs running; NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:03:34 gate closed for w2-control: REMOTE BUSY: load 1118 (max 150), available memory 71 GB (min 60), 92 jobs running; NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:05:07 gate closed for w2-control: REMOTE BUSY: load 1158 (max 150), available memory 54 GB (min 60), 91 jobs running; NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:20:49 gate closed for w2-control: 00), available memory 237 GB (min 90), 63 jobs running with 124 cpus allocated (max 220); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:22:22 gate closed for w2-control: 00), available memory 235 GB (min 90), 67 jobs running with 152 cpus allocated (max 220); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:23:54 gate closed for w2-control: 00), available memory 229 GB (min 90), 67 jobs running with 148 cpus allocated (max 220); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:25:27 gate closed for w2-control: 00), available memory 190 GB (min 90), 66 jobs running with 146 cpus allocated (max 220); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
20:26:59 gate closed for w2-control: 00), available memory 212 GB (min 90), 67 jobs running with 148 cpus allocated (max 220); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
```

## Notes and interpretation

* Wave-1 (attack1) facts used: E1 kernel is zero for d <= 8 (so every legal eigen denominator of degree <= 8 lies in J(F));
  d=9 kernels 2,3,3,2,10,1,3 (even) / 0,5,0,0,9,6,4 (odd); all tested fixed denominators inconsistent.
* The control test (accepted witness `hesse-triangle`, family {D*} + J_6, 18 parameters): the solver reproduces rel ~1e-30
  at the witness direction and converges to it from perturbations of relative size ~0.3-0.5, but random restarts hit the
  basin only ~1% of the time (0/26 in local tests; remote job `w2-control` quantifies this).  Consequence: for d >= 9 the
  random-restart joint search is a landscape survey, not an exhaustive method; for d <= 8 the decisive tool is the
  low-dimensional mod-J relaxation (necessary condition), which is settled reliably by multistart and rigorously by Singular.
* Degree 6: Singular dim = 0 in 13 of 14 sectors (rigorous exclusion of all sigma-eigen degree-6 denominators of those
  parities/weights); sector (6, even, 0) has a 2-dimensional mod-p cone (one component, an octic in three chart variables)
  but no real zeros numerically (300 restarts, best rel 8.2e-3).
* Degree 7: Singular dim = 0 for even k=0,3,5,6 and odd k=1..6 (rigorous); even k=1,2,4 have 2-dim mod-p cones with no real
  zeros (best rel 1.3e-3 .. 2.7e-3 over 200 restarts); odd k=0 (24 unknowns): global minimum rel 7.709e-6 (stable under LM
  polish from three different minima, full-E2 rel there 0.80) -- a near miss but nonzero; Singular timed out at 45 s locally.
* Exact E2 at degree 6 (`w2_exactE2.py`, F-reduced Koszul bilinear system in (a, p), saturated by the ideal of the
  denominator coordinates, Singular over GF(32003)): sector (6, even, 0) -- the only one not killed by the relaxation --
  has DIM_SAT = -1 (empty): no sigma-eigen Y-even weight-0 denominator of degree 6 admits a certificate mod 32003.
  Together with the 13 relaxation exclusions: **no sigma-eigen denominator of degree 6 (either parity, any weight)**.
  Caveat on rigor: the relaxation statements are homogeneous in the denominator coordinates alone, so mod-p emptiness
  implies emptiness over Q-bar; the saturated E2 statement involves the affine numerator coordinates and a single prime,
  so it is strong evidence (to be repeated at further primes / in characteristic 0), not yet a proof.
* Degree 8 relaxation probes (local, 60-120 restarts, LM-polished minima): (8,even,0) global min rel 2.198e-8 (stable,
  nonzero; full-E2 rel there 0.25), (8,even,1) 2.45e-5, (8,even,4) 1.58e-5.
* (6, even, 0) exact E2 saturation repeated at p = 1000003: again DIM_SAT = -1 (58 s).  Two independent primes.
* Degree 8 relaxation, all 14 sectors probed locally (40-120 restarts each, LM-polished): global minima are nonzero,
  rel between 6.4e-7 (8,odd,4) and 3.0e-5 (8,odd,0); no real zeros.  Remote job `w2-hess8` repeats this with 1500 restarts
  and Singular (dim of the exact quadric system).
* Correction to the wave-1/dossier description: fqlx's degree-9 denominator `Jhat` lies in J(F)_9 (membership verified
  mod 2147483647 and 2147483629), i.e. it is of Jacobian type, not a non-Jacobian E1-kernel element.
* Sector E1 kernels at d = 9: the nontrivial (mod J) kernel directions in sectors (9,even,1) and (9,even,2) have no
  representative of height below ~1e21 (LLL on the mod-m kernel lattice, m = product of 4 primes, finds exactly the 63
  short J-vectors and nothing else), while (9,even,0/3/5/6) have representatives of height 6e3..3e5.  For the large-height
  sectors the joint search uses an orthonormal floating-point kernel basis (`float_kernel_reps`, lifts to 4e-15).

