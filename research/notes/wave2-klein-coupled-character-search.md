# Wave 2 search report: klein-coupled-character (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (search in progress)

Updated 2026-11-01 22:38:37 from 934 log lines in `/tmp/hodge/results/klein-coupled-character/wave2-search.log.jsonl`.


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
| 9 | even | -,156,156,139,156,139,139 | -,63,63,52,63,52,52 | -,3,3,2,10,1,3 | -,-,-,49994,-,31345,272340 | 0 |
| 9 | odd | -,139,-,-,139,136,136 | -,49,-,-,49,46,46 | -,5,-,-,9,6,4 | -,2420,-,-,-,3625,3500 | 0 |

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
| 9 | even | 0 | KJ | 35 | 34 | 1949x240 | 4 | 1.64e-01 | 1.86e-01 | 2.70e-01 | 0 | 0 | lm | 0 |
| 9 | even | 1 | K | 3 | 3 | 1933x272 | 5 | 1.26e+00 | 1.26e+00 | 1.26e+00 | 0 | 0 | lm | 0 |
| 9 | even | 1 | KJ | 66 | 59 | 1933x272 | 2 | 5.38e-02 | 5.38e-02 | 5.38e-02 | 0 | 0 | lm | 0 |
| 9 | even | 2 | K | 3 | 3 | 1918x240 | 5 | 1.42e+00 | 1.42e+00 | 1.42e+00 | 0 | 0 | lm | 0 |
| 9 | even | 2 | KJ | 66 | 59 | 1918x240 | 2 | 4.16e-02 | 5.05e-02 | 5.94e-02 | 0 | 0 | lm | 0 |
| 9 | even | 3 | K | 2 | 2 | 1969x240 | 4 | 3.41e-02 | 3.41e-02 | 3.41e-02 | 0 | 0 | lm | 0 |
| 9 | even | 3 | KJ | 54 | 52 | 1972x240 | 3 | 1.75e-02 | 1.75e-02 | 1.75e-02 | 0 | 0 | lm | 0 |
| 9 | even | 4 | K | 10 | 10 | 1932x354 | 6 | 7.51e-01 | 7.51e-01 | 7.51e-01 | 0 | 0 | lm | 0 |
| 9 | even | 4 | KJ | 73 | 66 | 1932x354 | 2 | 5.61e-01 | 5.61e-01 | 5.61e-01 | 0 | 0 | lm | 0 |
| 9 | even | 5 | KJ | 53 | 51 | 1953x272 | 3 | 1.22e-02 | 4.54e-02 | 7.85e-02 | 0 | 0 | lm | 0 |
| 9 | even | 6 | K | 3 | 3 | 1966x272 | 6 | 1.69e-01 | 8.54e-01 | 8.54e-01 | 0 | 0 | lm | 0 |
| 9 | even | 6 | KJ | 55 | 53 | 1966x272 | 3 | 5.04e-02 | 5.37e-02 | 5.70e-02 | 0 | 0 | lm | 0 |
| 9 | odd | 1 | K | 5 | 5 | 1947x264 | 5 | 8.90e-01 | 8.90e-01 | 8.90e-01 | 0 | 0 | lm | 0 |
| 9 | odd | 1 | KJ | 54 | 51 | 1947x264 | 3 | 7.73e-02 | 1.18e-01 | 1.58e-01 | 0 | 0 | lm | 0 |
| 9 | odd | 4 | K | 9 | 9 | 1973x216 | 6 | 8.66e-01 | 8.66e-01 | 8.66e-01 | 0 | 0 | lm | 0 |
| 9 | odd | 4 | KJ | 58 | 55 | 1973x216 | 3 | 1.28e-01 | 2.51e-01 | 3.75e-01 | 0 | 0 | lm | 0 |
| 9 | odd | 5 | K | 6 | 6 | 1929x264 | 5 | 7.55e-01 | 7.55e-01 | 7.55e-01 | 0 | 0 | lm | 0 |
| 9 | odd | 5 | KJ | 52 | 49 | 1929x264 | 3 | 4.77e-01 | 4.77e-01 | 4.77e-01 | 0 | 0 | lm | 0 |
| 9 | odd | 6 | K | 4 | 4 | 1939x264 | 5 | 1.26e+00 | 1.26e+00 | 1.26e+00 | 0 | 0 | lm | 0 |
| 9 | odd | 6 | KJ | 50 | 47 | 1939x264 | 2 | 5.36e-02 | 5.95e-02 | 6.55e-02 | 0 | 0 | lm | 0 |

## (4) Structured families (E1-admissible dimension per sector; varpro where nonzero)

Cells with zero E1-admissible part (excluded as supports), by family: G-multiples:9, K+G-multiples:9, K-multiples:9, M-multiples:9, Xdeg_0:3, Xdeg_1:9, Xdeg_2:18, Xdeg_3:21, Xdeg_4:21, Xdeg_5:21, Xdeg_6:21, Xdeg_7:14, Xdeg_8:7, Xdeg_le3:42, Ydeg_le3:4, no_r:25, no_s:25, no_t:25, no_x:19, no_y:19, no_z:19, pureX:23, pureY:3, vars_trsx:25, vars_trsy:25, vars_trsz:25, vars_xyzr:27, vars_xyzs:27, vars_xyzt:27

721 (family, d, parity, k) cells tested; 190 with a nonzero E1-admissible part (listed), the rest are excluded as supports (mod-p kernel zero).

| family | d | parity | k | dim V | dim E1 | nontrivial mod J | dimP | restarts | best rel | exact tests |
|---|---|---|---|---|---|---|---|---|---|---|
| EX-multiples | 6 | 0 | 0 | 1 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| EY-multiples | 6 | 0 | 0 | 1 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 0 | 22 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| no_r | 6 | 0 | 0 | 20 | 3 | 0 | 3 | 6 | 3.87e-01 | inconsistent,inconsistent,inconsistent |
| no_s | 6 | 0 | 0 | 20 | 3 | 0 | 3 | 6 | 8.21e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 6 | 0 | 0 | 20 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| no_x | 6 | 0 | 0 | 35 | 2 | 0 | 2 | 6 | 7.76e-01 | inconsistent,inconsistent |
| vars_xyzr | 6 | 0 | 0 | 8 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 6 | 0 | 0 | 8 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzt | 6 | 0 | 0 | 8 | 1 | 0 | 1 | 6 | 5.85e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 1 | 16 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 2 | 16 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 3 | 16 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 4 | 16 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 5 | 16 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 0 | 6 | 16 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 6 | 1 | 3 | 29 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 6 | 1 | 5 | 29 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 6 | 1 | 6 | 29 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 0 | 0 | 24 | 3 | 0 | 3 | 6 | 7.61e-01 | inconsistent,inconsistent,inconsistent |
| EX-multiples | 7 | 0 | 1 | 1 | 1 | 0 | 1 | 6 | 6.82e-01 | inconsistent |
| EY-multiples | 7 | 0 | 1 | 1 | 1 | 0 | 1 | 6 | 8.48e-01 | inconsistent |
| Ydeg_le3 | 7 | 0 | 1 | 23 | 3 | 0 | 3 | 6 | 1.47e-01 | inconsistent,inconsistent,inconsistent |
| no_r | 7 | 0 | 1 | 26 | 3 | 0 | 3 | 6 | 3.71e-01 | inconsistent,inconsistent,inconsistent |
| no_s | 7 | 0 | 1 | 26 | 3 | 0 | 3 | 6 | 2.86e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 7 | 0 | 1 | 26 | 3 | 0 | 3 | 6 | 4.22e-01 | inconsistent,inconsistent,inconsistent |
| no_y | 7 | 0 | 1 | 35 | 2 | 0 | 2 | 6 | 7.70e-01 | inconsistent,inconsistent |
| no_z | 7 | 0 | 1 | 35 | 2 | 0 | 2 | 6 | 8.09e-01 | inconsistent,inconsistent |
| vars_xyzr | 7 | 0 | 1 | 10 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 7 | 0 | 1 | 10 | 1 | 0 | 1 | 6 | 9.85e-01 | inconsistent |
| vars_xyzt | 7 | 0 | 1 | 10 | 1 | 0 | 1 | 6 | 6.71e+00 | inconsistent |
| EX-multiples | 7 | 0 | 2 | 1 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| EY-multiples | 7 | 0 | 2 | 1 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| Ydeg_le3 | 7 | 0 | 2 | 23 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| no_r | 7 | 0 | 2 | 26 | 3 | 0 | 3 | 6 | 3.87e-01 | inconsistent,inconsistent,inconsistent |
| no_s | 7 | 0 | 2 | 26 | 3 | 0 | 3 | 6 | 8.21e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 7 | 0 | 2 | 26 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| no_x | 7 | 0 | 2 | 35 | 2 | 0 | 2 | 6 | 7.76e-01 | inconsistent,inconsistent |
| no_z | 7 | 0 | 2 | 35 | 2 | 0 | 2 | 6 | 8.29e-01 | inconsistent,inconsistent |
| vars_xyzr | 7 | 0 | 2 | 10 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 7 | 0 | 2 | 10 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzt | 7 | 0 | 2 | 10 | 1 | 0 | 1 | 6 | 5.85e+00 | inconsistent |
| Ydeg_le3 | 7 | 0 | 3 | 23 | 2 | 0 | 2 | 6 | 9.98e-01 | inconsistent,inconsistent |
| EX-multiples | 7 | 0 | 4 | 1 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| EY-multiples | 7 | 0 | 4 | 1 | 1 | 0 | 1 | 6 | 1.50e-01 | inconsistent |
| Ydeg_le3 | 7 | 0 | 4 | 23 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| no_r | 7 | 0 | 4 | 26 | 3 | 0 | 3 | 6 | 5.20e-02 | inconsistent,inconsistent,inconsistent |
| no_s | 7 | 0 | 4 | 26 | 3 | 0 | 3 | 6 | 8.15e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 7 | 0 | 4 | 26 | 3 | 0 | 3 | 6 | 1.01e-02 | inconsistent,inconsistent,inconsistent |
| no_x | 7 | 0 | 4 | 35 | 2 | 0 | 2 | 6 | 2.51e-02 | inconsistent,inconsistent |
| no_y | 7 | 0 | 4 | 35 | 2 | 0 | 2 | 6 | 2.06e-02 | inconsistent,inconsistent |
| vars_xyzr | 7 | 0 | 4 | 10 | 1 | 0 | 1 | 6 | 9.52e-01 | inconsistent |
| vars_xyzs | 7 | 0 | 4 | 10 | 1 | 0 | 1 | 6 | 5.61e-02 | inconsistent |
| vars_xyzt | 7 | 0 | 4 | 10 | 1 | 0 | 1 | 6 | 5.85e+00 | inconsistent |
| Ydeg_le3 | 7 | 0 | 5 | 23 | 2 | 0 | 2 | 6 | 8.92e-01 | inconsistent,inconsistent |
| Ydeg_le3 | 7 | 0 | 6 | 23 | 2 | 0 | 2 | 6 | 5.34e-01 | inconsistent,inconsistent |
| EX-multiples | 7 | 1 | 0 | 3 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| EY-multiples | 7 | 1 | 0 | 3 | 3 | 0 | 3 | 6 | 1.55e-02 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 1 | 0 | 42 | 9 | 0 | 9 | 6 | 9.76e-01 | - |
| no_r | 7 | 1 | 0 | 28 | 6 | 0 | 6 | 6 | 3.83e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 7 | 1 | 0 | 28 | 6 | 0 | 6 | 6 | 7.96e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 7 | 1 | 0 | 28 | 6 | 0 | 6 | 6 | 7.04e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_x | 7 | 1 | 0 | 49 | 6 | 0 | 6 | 6 | 4.69e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_y | 7 | 1 | 0 | 49 | 6 | 0 | 6 | 6 | 2.78e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_z | 7 | 1 | 0 | 49 | 6 | 0 | 6 | 6 | 2.37e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| vars_xyzr | 7 | 1 | 0 | 8 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 7 | 1 | 0 | 8 | 1 | 0 | 1 | 6 | 9.99e-01 | inconsistent |
| vars_xyzt | 7 | 1 | 0 | 8 | 1 | 0 | 1 | 6 | 3.16e+01 | inconsistent |
| Ydeg_le3 | 7 | 1 | 1 | 32 | 3 | 0 | 3 | 6 | 9.97e-01 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 1 | 2 | 32 | 3 | 0 | 3 | 6 | 9.39e-01 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 1 | 3 | 32 | 3 | 0 | 3 | 6 | 9.94e-01 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 1 | 4 | 32 | 3 | 0 | 3 | 6 | 1.00e+00 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 1 | 5 | 32 | 3 | 0 | 3 | 6 | 6.61e-01 | inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 7 | 1 | 6 | 32 | 3 | 0 | 3 | 6 | 9.84e-01 | inconsistent,inconsistent,inconsistent |
| EX-multiples | 8 | 0 | 0 | 6 | 6 | 0 | 6 | 6 | 9.45e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| EY-multiples | 8 | 0 | 0 | 6 | 6 | 0 | 6 | 4 | 1.51e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| Ydeg_le3 | 8 | 0 | 0 | 30 | 3 | 0 | 3 | 6 | 5.06e-02 | inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 0 | 42 | 9 | 0 | 9 | 4 | 1.38e-01 | - |
| no_s | 8 | 0 | 0 | 42 | 9 | 0 | 9 | 4 | 5.12e-03 | - |
| no_t | 8 | 0 | 0 | 42 | 9 | 0 | 9 | 4 | 6.05e-03 | - |
| no_x | 8 | 0 | 0 | 67 | 12 | 0 | 12 | 6 | 4.03e-02 | - |
| no_y | 8 | 0 | 0 | 67 | 12 | 0 | 12 | 6 | 1.21e-02 | - |
| no_z | 8 | 0 | 0 | 67 | 12 | 0 | 12 | 6 | 2.89e-03 | - |
| vars_xyzr | 8 | 0 | 0 | 14 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 8 | 0 | 0 | 14 | 1 | 0 | 1 | 6 | 9.99e-01 | inconsistent |
| vars_xyzt | 8 | 0 | 0 | 14 | 1 | 0 | 1 | 6 | 3.33e+01 | inconsistent |
| EX-multiples | 8 | 0 | 1 | 1 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| EY-multiples | 8 | 0 | 1 | 1 | 1 | 0 | 1 | 6 | 1.50e-01 | inconsistent |
| Ydeg_le3 | 8 | 0 | 1 | 31 | 5 | 0 | 5 | 6 | 9.94e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 1 | 36 | 3 | 0 | 3 | 6 | 5.20e-02 | inconsistent,inconsistent,inconsistent |
| no_s | 8 | 0 | 1 | 36 | 4 | 0 | 4 | 6 | 8.15e-01 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 8 | 0 | 1 | 36 | 4 | 0 | 4 | 6 | 1.01e-02 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_x | 8 | 0 | 1 | 50 | 2 | 0 | 2 | 6 | 2.51e-02 | inconsistent,inconsistent |
| no_y | 8 | 0 | 1 | 36 | 2 | 0 | 2 | 6 | 2.06e-02 | inconsistent,inconsistent |
| vars_xyzr | 8 | 0 | 1 | 14 | 1 | 0 | 1 | 6 | 9.52e-01 | inconsistent |
| vars_xyzs | 8 | 0 | 1 | 14 | 1 | 0 | 1 | 6 | 5.61e-02 | inconsistent |
| vars_xyzt | 8 | 0 | 1 | 14 | 1 | 0 | 1 | 6 | 5.85e+00 | inconsistent |
| EX-multiples | 8 | 0 | 2 | 1 | 1 | 0 | 1 | 6 | 1.49e-02 | inconsistent |
| EY-multiples | 8 | 0 | 2 | 1 | 1 | 0 | 1 | 6 | 1.43e-02 | inconsistent |
| Ydeg_le3 | 8 | 0 | 2 | 31 | 5 | 0 | 5 | 6 | 9.86e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 2 | 36 | 4 | 0 | 4 | 6 | 1.25e-02 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 8 | 0 | 2 | 36 | 4 | 0 | 4 | 6 | 6.83e-03 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 8 | 0 | 2 | 36 | 3 | 0 | 3 | 6 | 2.75e-03 | inconsistent,inconsistent,inconsistent |
| no_y | 8 | 0 | 2 | 50 | 2 | 0 | 2 | 6 | 2.15e-02 | inconsistent,inconsistent |
| no_z | 8 | 0 | 2 | 36 | 2 | 0 | 2 | 6 | 2.25e-02 | inconsistent,inconsistent |
| vars_xyzr | 8 | 0 | 2 | 14 | 1 | 0 | 1 | 6 | 3.69e-02 | inconsistent |
| vars_xyzs | 8 | 0 | 2 | 14 | 1 | 0 | 1 | 6 | 3.55e-02 | inconsistent |
| vars_xyzt | 8 | 0 | 2 | 14 | 1 | 0 | 1 | 6 | 8.79e-02 | inconsistent |
| EX-multiples | 8 | 0 | 3 | 1 | 1 | 0 | 1 | 6 | 8.42e-02 | inconsistent |
| EY-multiples | 8 | 0 | 3 | 1 | 1 | 0 | 1 | 6 | 8.11e-02 | inconsistent |
| Ydeg_le3 | 8 | 0 | 3 | 30 | 4 | 0 | 4 | 6 | 6.58e-02 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 3 | 35 | 3 | 0 | 3 | 6 | 7.27e-02 | inconsistent,inconsistent,inconsistent |
| no_s | 8 | 0 | 3 | 35 | 3 | 0 | 3 | 6 | 2.34e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 8 | 0 | 3 | 35 | 3 | 0 | 3 | 6 | 1.86e-01 | inconsistent,inconsistent,inconsistent |
| no_z | 8 | 0 | 3 | 35 | 2 | 0 | 2 | 6 | 4.87e-01 | inconsistent,inconsistent |
| vars_xyzr | 8 | 0 | 3 | 13 | 1 | 0 | 1 | 6 | 5.61e-01 | inconsistent |
| vars_xyzs | 8 | 0 | 3 | 13 | 1 | 0 | 1 | 6 | 4.83e-01 | inconsistent |
| vars_xyzt | 8 | 0 | 3 | 13 | 1 | 0 | 1 | 6 | 3.46e+00 | inconsistent |
| EX-multiples | 8 | 0 | 4 | 1 | 1 | 0 | 1 | 6 | 9.97e-01 | inconsistent |
| EY-multiples | 8 | 0 | 4 | 1 | 1 | 0 | 1 | 6 | 9.97e-01 | inconsistent |
| Ydeg_le3 | 8 | 0 | 4 | 31 | 5 | 0 | 5 | 6 | 9.29e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 4 | 36 | 4 | 0 | 4 | 6 | 3.78e-01 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 8 | 0 | 4 | 36 | 3 | 0 | 3 | 6 | 8.09e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 8 | 0 | 4 | 36 | 4 | 0 | 4 | 6 | 8.82e-01 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_x | 8 | 0 | 4 | 36 | 2 | 0 | 2 | 6 | 7.76e-01 | inconsistent,inconsistent |
| no_z | 8 | 0 | 4 | 50 | 2 | 0 | 2 | 6 | 8.18e-01 | inconsistent,inconsistent |
| vars_xyzr | 8 | 0 | 4 | 14 | 1 | 0 | 1 | 6 | 9.33e-01 | inconsistent |
| vars_xyzs | 8 | 0 | 4 | 14 | 1 | 0 | 1 | 6 | 9.90e-01 | inconsistent |
| vars_xyzt | 8 | 0 | 4 | 14 | 1 | 0 | 1 | 6 | 5.69e+00 | inconsistent |
| EX-multiples | 8 | 0 | 5 | 1 | 1 | 0 | 1 | 6 | 6.81e-01 | inconsistent |
| EY-multiples | 8 | 0 | 5 | 1 | 1 | 0 | 1 | 6 | 1.27e-01 | inconsistent |
| Ydeg_le3 | 8 | 0 | 5 | 30 | 4 | 0 | 4 | 6 | 4.66e-01 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 5 | 35 | 3 | 0 | 3 | 6 | 4.69e-02 | inconsistent,inconsistent,inconsistent |
| no_s | 8 | 0 | 5 | 35 | 3 | 0 | 3 | 6 | 2.86e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 8 | 0 | 5 | 35 | 3 | 0 | 3 | 6 | 7.18e-03 | inconsistent,inconsistent,inconsistent |
| no_y | 8 | 0 | 5 | 35 | 2 | 0 | 2 | 6 | 2.06e-02 | inconsistent,inconsistent |
| vars_xyzr | 8 | 0 | 5 | 13 | 1 | 0 | 1 | 6 | 9.53e-01 | inconsistent |
| vars_xyzs | 8 | 0 | 5 | 13 | 1 | 0 | 1 | 6 | 5.55e-02 | inconsistent |
| vars_xyzt | 8 | 0 | 5 | 13 | 1 | 0 | 1 | 6 | 6.71e+00 | inconsistent |
| EX-multiples | 8 | 0 | 6 | 1 | 1 | 0 | 1 | 6 | 9.82e-01 | inconsistent |
| EY-multiples | 8 | 0 | 6 | 1 | 1 | 0 | 1 | 6 | 1.47e-01 | inconsistent |
| Ydeg_le3 | 8 | 0 | 6 | 30 | 4 | 0 | 4 | 6 | 8.89e-01 | inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 0 | 6 | 35 | 3 | 0 | 3 | 6 | 4.93e-02 | inconsistent,inconsistent,inconsistent |
| no_s | 8 | 0 | 6 | 35 | 3 | 0 | 3 | 6 | 8.15e-01 | inconsistent,inconsistent,inconsistent |
| no_t | 8 | 0 | 6 | 35 | 3 | 0 | 3 | 6 | 9.95e-03 | inconsistent,inconsistent,inconsistent |
| no_x | 8 | 0 | 6 | 35 | 2 | 0 | 2 | 6 | 2.52e-02 | inconsistent,inconsistent |
| vars_xyzr | 8 | 0 | 6 | 13 | 1 | 0 | 1 | 6 | 9.52e-01 | inconsistent |
| vars_xyzs | 8 | 0 | 6 | 13 | 1 | 0 | 1 | 6 | 5.61e-02 | inconsistent |
| vars_xyzt | 8 | 0 | 6 | 13 | 1 | 0 | 1 | 6 | 5.84e+00 | inconsistent |
| Ydeg_le3 | 8 | 1 | 0 | 48 | 9 | 0 | 9 | 6 | 5.49e-01 | - |
| no_r | 8 | 1 | 0 | 30 | 1 | 0 | 1 | 6 | 1.15e+00 | inconsistent |
| no_s | 8 | 1 | 0 | 30 | 1 | 0 | 1 | 6 | 9.35e-01 | inconsistent |
| no_t | 8 | 1 | 0 | 30 | 1 | 0 | 1 | 6 | 8.48e-01 | inconsistent |
| Ydeg_le3 | 8 | 1 | 1 | 45 | 9 | 0 | 9 | 6 | 1.43e-01 | - |
| no_r | 8 | 1 | 1 | 36 | 6 | 0 | 6 | 6 | 3.94e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 8 | 1 | 1 | 36 | 6 | 0 | 6 | 6 | 2.76e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 8 | 1 | 1 | 36 | 6 | 0 | 6 | 6 | 4.58e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_y | 8 | 1 | 1 | 49 | 6 | 0 | 6 | 6 | 2.78e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_z | 8 | 1 | 1 | 49 | 6 | 0 | 6 | 6 | 2.31e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| vars_xyzr | 8 | 1 | 1 | 10 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 8 | 1 | 1 | 10 | 1 | 0 | 1 | 6 | 9.84e-01 | inconsistent |
| vars_xyzt | 8 | 1 | 1 | 10 | 1 | 0 | 1 | 6 | 1.88e+01 | inconsistent |
| Ydeg_le3 | 8 | 1 | 2 | 45 | 9 | 0 | 9 | 6 | 9.26e-01 | - |
| no_r | 8 | 1 | 2 | 36 | 6 | 0 | 6 | 6 | 3.76e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 8 | 1 | 2 | 36 | 6 | 0 | 6 | 6 | 7.96e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 8 | 1 | 2 | 36 | 6 | 0 | 6 | 6 | 7.04e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_x | 8 | 1 | 2 | 49 | 6 | 0 | 6 | 6 | 4.69e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_z | 8 | 1 | 2 | 49 | 6 | 0 | 6 | 6 | 9.41e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| vars_xyzr | 8 | 1 | 2 | 10 | 1 | 0 | 1 | 6 | 1.00e+00 | inconsistent |
| vars_xyzs | 8 | 1 | 2 | 10 | 1 | 0 | 1 | 6 | 9.99e-01 | inconsistent |
| vars_xyzt | 8 | 1 | 2 | 10 | 1 | 0 | 1 | 6 | 3.16e+01 | inconsistent |
| Ydeg_le3 | 8 | 1 | 3 | 45 | 6 | 0 | 6 | 6 | 9.90e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 8 | 1 | 3 | 34 | 1 | 0 | 1 | 6 | 1.59e-01 | inconsistent |
| Ydeg_le3 | 8 | 1 | 4 | 45 | 9 | 0 | 9 | 6 | 9.79e-01 | - |
| no_r | 8 | 1 | 4 | 36 | 6 | 0 | 6 | 6 | 6.51e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_s | 8 | 1 | 4 | 36 | 6 | 0 | 6 | 6 | 7.67e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 8 | 1 | 4 | 36 | 6 | 0 | 6 | 6 | 1.46e-04 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_x | 8 | 1 | 4 | 49 | 6 | 0 | 6 | 6 | 6.27e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_y | 8 | 1 | 4 | 49 | 6 | 0 | 6 | 6 | 6.64e-03 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| vars_xyzr | 8 | 1 | 4 | 10 | 1 | 0 | 1 | 6 | 9.52e-01 | inconsistent |
| vars_xyzs | 8 | 1 | 4 | 10 | 1 | 0 | 1 | 6 | 3.00e-03 | inconsistent |
| vars_xyzt | 8 | 1 | 4 | 10 | 1 | 0 | 1 | 6 | 2.00e+01 | inconsistent |
| Ydeg_le3 | 8 | 1 | 5 | 45 | 6 | 0 | 6 | 6 | 8.28e-02 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_t | 8 | 1 | 5 | 34 | 1 | 0 | 1 | 6 | 9.96e-01 | inconsistent |
| Ydeg_le3 | 8 | 1 | 6 | 45 | 6 | 0 | 6 | 6 | 4.23e-01 | inconsistent,inconsistent,inconsistent,inconsistent,inconsistent,inconsistent |
| no_r | 8 | 1 | 6 | 34 | 1 | 0 | 1 | 6 | 7.45e-02 | inconsistent |
| vars_xyzr | 9 | 0 | 0 | 17 | 1 | 0 | 1 | 6 | 5.36e-01 | inconsistent |
| vars_xyzs | 9 | 0 | 0 | 17 | 1 | 0 | 1 | 6 | 4.55e-01 | inconsistent |
| vars_xyzt | 9 | 0 | 0 | 17 | 1 | 0 | 1 | 6 | 3.17e+00 | inconsistent |
| vars_xyzt | 9 | 0 | 1 | 18 | 2 | 0 | 2 | 6 | 4.50e+00 | inconsistent,inconsistent |

## Remote job queue (tail)
```
22:27:29 gate closed for w2-control: 0), available memory 163 GB (min 110), 61 jobs running with 320 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:29:04 gate closed for w2-control: 0), available memory 134 GB (min 110), 62 jobs running with 328 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:30:38 gate closed for w2-control: 0), available memory 161 GB (min 110), 62 jobs running with 328 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:32:13 gate closed for w2-control: 0), available memory 159 GB (min 110), 61 jobs running with 327 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:33:47 gate closed for w2-control: 0), available memory 168 GB (min 110), 62 jobs running with 329 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:35:21 gate closed for w2-control: 0), available memory 155 GB (min 110), 63 jobs running with 333 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:36:55 gate closed for w2-control: 0), available memory 154 GB (min 110), 62 jobs running with 332 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
22:38:29 gate closed for w2-control: 0), available memory 145 GB (min 110), 62 jobs running with 332 cpus allocated (max 200); NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).
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
* Degree-9 kernel-only families (the non-Jacobian E1 directions, 1-10 dimensional per sector; rational reps where the
  height is small, orthonormal float reps otherwise), full-E2 LM, local runs (3-6 restarts each, all converging to the
  same minima across seeds): best rel (9,e,3) 3.4e-2, (9,e,6) 0.17, (9,e,4) 0.75, (9,o,5) 0.75, (9,o,4) 0.87, (9,o,1) 0.89,
  (9,e,1) 1.26, (9,o,6) 1.26, (9,e,2) 1.42; (9,e,0) 1.03 and (9,e,5) is a single direction (fixed-D, inconsistent in wave 1).
  No zeros: no certificate denominator in the span of the non-Jacobian kernel directions alone.
* Structured families (part 4, local runs at degrees 6-7, continued on the remote): pure-X cones (vertex P^2_Y) and
  pure-Y cones (vertex P^2_X) have no E1-admissible member in any sector at d = 6, 7 (kernel zero mod p => excluded as
  supports); line-vertex cones in (x,y,z,t), (x,y,z,r), (x,y,z,s) have a single E1-admissible direction each at (6,even,0)
  (a J-element), exactly tested: inconsistent; point-vertex cones with x missing: 2 directions at (6,even,0), both
  inconsistent.
* Why random-restart joint search cannot be exhaustive for d >= 9: a certificate (if any) is an isolated point of an
  overdetermined bilinear system (about 10x more equations than unknowns after eliminating Q); on the 18-dimensional
  control family the LM basin covers about 1% of the sphere and for 70-250 dimensional sector families the fraction is
  astronomically small.  The searches at d >= 9 are therefore landscape surveys plus exhaustive searches of the
  low-dimensional kernel-only and structured subfamilies; the decisive statements are the relaxation/exact-E2 exclusions
  at d <= 8 and the E1-kernel dimension table.
* Cross-reference: the parallel exclusion packet in this lane (`wave2-excl/relaxation-verdicts.md`, prime 1000003 and
  exact GB) reaches the same sector verdicts for d <= 7 by a different route: (6,even) = weighted-Euler curve C, all points
  excluded by exact constant-dual quadrics; (7,even,1/2/4) = X_i.C lifted components (excluded via the NZD lemma and the
  sum principle); (7,odd,0) lifted span(Y).C (pending); all other d <= 7 sectors EMPTY.  My relaxation/exact-E2 results
  (mod 32003 and 1000003) are consistent with every entry of that table.
* Structured families, local coverage: degrees 6 and 7 complete for all families (cones with plane/line/point vertex,
  fixed X-degree pieces, X-degree <= 3 / Y-degree <= 3 pieces, multiples of K, G, M, K+G and of the block Euler derivatives
  E_X F, E_Y F) in all 14 sectors; degree 8 complete for plane/line-vertex cones, partial for point-vertex cones.  Every
  nonzero E1-admissible part found so far lies inside J(F)_d (nontrivial mod J = 0) and has no E2 zero (best rel >= 3e-3,
  all exactly tested representatives inconsistent).  The remote job `w2-struct` extends this to degrees 8-12.
* Relaxation summary (best rel over all restarts; S = Singular dim of the exact mod-32003 quadric system, 0 = rigorous
  exclusion): d=6 even k=0 8.2e-3 (S=2, excluded by exact E2 saturation at two primes), k=1..6 0.13..0.86 (S=0); d=6 odd
  k=1..6 5e-3..0.15 (S=0); d=7 even k=0,3,5,6 (S=0), k=1,2,4 1.3e-3..2.7e-3 (S=2 = X_i.C lifts, cf. exclusion packet);
  d=7 odd k=0 7.7e-6 (24 unknowns, stable nonzero minimum), k=1..6 (S=0); d=8 even 2.1e-8 (k=0, 45 unknowns) .. 2.5e-5,
  d=8 odd 6.4e-7 (k=4) .. 3.0e-5 -- no real zeros anywhere at d <= 8.
* Second prime (1000003) for the exact relaxation dimensions at d = 6 and 7: identical verdicts to p = 32003 in all
  26 sectors (dim 0 everywhere except (6,even,0) and (7,even,1/2/4) with dim 2).
* Closest structured near miss: (8, odd, 4), cones not involving Y0 (6-dimensional E1 part inside J_8): 36 LM restarts give
  7 distinct local minima, best rel 1.46e-4 (single restart) / 2.36e-4 (30-restart batch), all exactly tested
  representatives inconsistent -- not a zero.
* Structured families, degree 8 (local): plane-, line- and point-vertex cone families now complete in all 14 sectors;
  every nonzero E1-admissible part lies in J_8 and has no E2 zero (smallest residuals 1.5e-4 .. 3e-3, see table).
* Structured families now complete locally for degrees 6, 7 and 8 (all families, all 14 sectors): no E2 zero anywhere;
  all nonzero E1-admissible parts are Jacobian-type.  Degrees 9-12 of the structured sweep remain for the remote job.
* Degree-9 structured sweep (local, partial): (9,even,0) and (9,even,1) plane/line-vertex cones -- E1 parts of dimension
  1-2 (Jacobian type), best rel >= 0.45; the remaining degree 9-12 cells are left to the remote job `w2-struct`.
* Control hit rate (hesse-triangle, 18-dim family {D*}+J_6, reduced system, varpro-LM): seeds 7/3 gave 0 hits in 26
  restarts; seed 101 found the witness (rel 8e-30) on the 4th restart.  Overall about 1 hit per 30 random restarts in
  18 dimensions; the reconstructed D agrees with rem_F(D*) and the exact re-solve confirms it (pipeline validated end to end).
* Control, seed 202: hit on the 3rd restart (rel 1.5e-24); the rem_F-chart reconstruction returned a rational point of the
  1-dimensional solution set different from D*, and the exact re-solve produced a verified certificate for it
  (confirm_exact -> solved, verify_witness True).  Running tally: 3 hits in 41 random restarts (~7%) on the 18-dim control.

