# Wave 2 — klein-coupled-character — SEARCH lens (joint denominator/numerator search)

STATUS: NO WITNESS. 27 search batches, 14800 restarts total, best relative residual 9.595e-02; 0 exact-confirmation attempts. Updated 2026-11-02 01:00 UTC. Lane dir /tmp/hodge/results/klein-coupled-character (w2fast.py engine, w2_search.py driver, w2_dims.py, w2_validate.py; batch log wave2-search.log.jsonl; per-restart streams wave2-<tag>-d*p*c*-w*.jsonl).

## Method

Target q = 1, A = X0^3 Y0^3 (Y-odd, sigma-character 3; weights W = (1,2,4,0,0,0) mod 7, F character 0). For a sigma-eigen denominator D of degree d, Y-parity parD, character chi the certificate system is restricted losslessly (averaging) to N_i (deg d+1, parity 1+parD+pi_i, char 3+chi+W_i), Koszul corrections p_ij (deg d-4, parity 1+parD+pi_i+pi_j, char 3+chi+W_i+W_j), Q (deg 2d-6, Y-odd, char 3+2chi). Correction-map form (ch. 5.2): with an E1 lift N0 (N0(F) = A D) the condition is B_D(N0) + sum p_ij L_ij = F Q. Joint search = variable projection over a linear family D(a) (E1-admissible sector space, lifts linear in a): f(a) = min_{p,Q} |T(a) p - F Q - b(a)|^2 on the unit sphere, F-multiples projected out, L-BFGS random restarts; rel = f/|b|^2. Inner solve: projected normal equations with eigen-cutoff 1e-7 (checked against LAPACK gelsd to 1e-16). Engine validated on accepted q=1 witnesses (f = 0 at the witness; Koszul-perturbed lifts absorbed to 1e-23 relative), gradient checked by finite differences. Families: full = (J:A)_{d,sector} (float null space of the sector E1 system = J_d + non-Jacobian kernel); w = J_d with lifts A w; hess = Hessian relaxation (N = A w exactly, only Q free). A numerical residual is NOT evidence; only exact re-solves (kcc2 Koszul test mod p + multi-prime lift + core.verify_witness) count.

## E1-admissible denominator spaces per sector (exact mod p = 2147483647; K_float = float null-space dimension cross-check)

| d | parD | chi | dim S_d sector | dim J_d | dim (J:A)_d | non-Jacobian kernel | F-multiples | nN | nP | nQ | rows |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 0 | 18 | 0 | 0 | 0 | 0 | 107 | 3 | 0 | 87 |
| 4 | 0 | 1 | 8 | 0 | 0 | 0 | 0 | 120 | 0 | 0 | 87 |
| 4 | 0 | 2 | 8 | 0 | 0 | 0 | 0 | 107 | 3 | 0 | 69 |
| 4 | 0 | 3 | 8 | 0 | 0 | 0 | 0 | 107 | 3 | 3 | 102 |
| 4 | 0 | 4 | 8 | 0 | 0 | 0 | 0 | 75 | 0 | 3 | 102 |
| 4 | 0 | 5 | 8 | 0 | 0 | 0 | 0 | 120 | 0 | 0 | 87 |
| 4 | 0 | 6 | 8 | 0 | 0 | 0 | 0 | 120 | 0 | 3 | 102 |
| 4 | 1 | 0 | 3 | 0 | 0 | 0 | 0 | 96 | 0 | 0 | 87 |
| 4 | 1 | 1 | 13 | 0 | 0 | 0 | 0 | 102 | 1 | 0 | 87 |
| 4 | 1 | 2 | 13 | 0 | 0 | 0 | 0 | 96 | 0 | 0 | 69 |
| 4 | 1 | 3 | 6 | 0 | 0 | 0 | 0 | 96 | 0 | 3 | 102 |
| 4 | 1 | 4 | 13 | 0 | 0 | 0 | 0 | 162 | 3 | 3 | 102 |
| 4 | 1 | 5 | 6 | 0 | 0 | 0 | 0 | 102 | 1 | 0 | 87 |
| 4 | 1 | 6 | 6 | 0 | 0 | 0 | 0 | 102 | 1 | 3 | 102 |
| 5 | 0 | 0 | 9 | 0 | 0 | 0 | 0 | 181 | 3 | 6 | 209 |
| 5 | 0 | 1 | 24 | 0 | 0 | 0 | 0 | 191 | 6 | 6 | 209 |
| 5 | 0 | 2 | 24 | 0 | 0 | 0 | 0 | 181 | 3 | 3 | 180 |
| 5 | 0 | 3 | 15 | 1 | 1 | 0 | 0 | 181 | 3 | 13 | 228 |
| 5 | 0 | 4 | 24 | 0 | 0 | 0 | 0 | 270 | 18 | 13 | 228 |
| 5 | 0 | 5 | 15 | 1 | 1 | 0 | 0 | 191 | 6 | 6 | 209 |
| 5 | 0 | 6 | 15 | 1 | 1 | 0 | 0 | 191 | 6 | 13 | 228 |
| 5 | 1 | 0 | 30 | 3 | 3 | 0 | 0 | 199 | 11 | 6 | 209 |
| 5 | 1 | 1 | 16 | 0 | 0 | 0 | 0 | 213 | 4 | 6 | 209 |
| 5 | 1 | 2 | 16 | 0 | 0 | 0 | 0 | 199 | 11 | 3 | 180 |
| 5 | 1 | 3 | 16 | 0 | 0 | 0 | 0 | 199 | 11 | 13 | 228 |
| 5 | 1 | 4 | 16 | 0 | 0 | 0 | 0 | 150 | 0 | 13 | 228 |
| 5 | 1 | 5 | 16 | 0 | 0 | 0 | 0 | 213 | 4 | 6 | 209 |
| 5 | 1 | 6 | 16 | 0 | 0 | 0 | 0 | 213 | 4 | 13 | 228 |
| 6 | 0 | 0 | 50 | 12 | 12 | 0 | 1 | 343 | 30 | 29 | 436 |
| 6 | 0 | 1 | 31 | 1 | 1 | 0 | 0 | 357 | 21 | 29 | 436 |
| 6 | 0 | 2 | 31 | 1 | 1 | 0 | 0 | 343 | 30 | 19 | 395 |
| 6 | 0 | 3 | 31 | 1 | 1 | 0 | 0 | 343 | 30 | 40 | 459 |
| 6 | 0 | 4 | 31 | 1 | 1 | 0 | 0 | 276 | 9 | 40 | 459 |
| 6 | 0 | 5 | 31 | 1 | 1 | 0 | 0 | 357 | 21 | 29 | 436 |
| 6 | 0 | 6 | 31 | 1 | 1 | 0 | 0 | 357 | 21 | 40 | 459 |
| 6 | 1 | 0 | 19 | 0 | 0 | 0 | 0 | 317 | 15 | 29 | 436 |
| 6 | 1 | 1 | 40 | 3 | 3 | 0 | 0 | 331 | 20 | 29 | 436 |
| 6 | 1 | 2 | 40 | 3 | 3 | 0 | 0 | 317 | 15 | 19 | 395 |
| 6 | 1 | 3 | 29 | 3 | 3 | 0 | 0 | 317 | 15 | 40 | 459 |
| 6 | 1 | 4 | 40 | 3 | 3 | 0 | 0 | 432 | 48 | 40 | 459 |
| 6 | 1 | 5 | 29 | 3 | 3 | 0 | 0 | 331 | 20 | 29 | 436 |
| 6 | 1 | 6 | 29 | 3 | 3 | 0 | 0 | 331 | 20 | 40 | 459 |
| 7 | 0 | 0 | 39 | 3 | 3 | 0 | 0 | 522 | 48 | 87 | 824 |
| 7 | 0 | 1 | 66 | 12 | 12 | 0 | 1 | 543 | 58 | 87 | 824 |
| 7 | 0 | 2 | 66 | 12 | 12 | 0 | 1 | 522 | 48 | 69 | 771 |
| 7 | 0 | 3 | 53 | 8 | 8 | 0 | 0 | 522 | 48 | 102 | 851 |
| 7 | 0 | 4 | 66 | 12 | 12 | 0 | 1 | 666 | 102 | 102 | 851 |
| 7 | 0 | 5 | 53 | 8 | 8 | 0 | 0 | 543 | 58 | 87 | 824 |
| 7 | 0 | 6 | 53 | 8 | 8 | 0 | 0 | 543 | 58 | 102 | 851 |
| 7 | 1 | 0 | 78 | 27 | 27 | 0 | 3 | 558 | 70 | 87 | 824 |
| 7 | 1 | 1 | 53 | 6 | 6 | 0 | 0 | 571 | 58 | 87 | 824 |
| 7 | 1 | 2 | 53 | 6 | 6 | 0 | 0 | 558 | 70 | 69 | 771 |
| 7 | 1 | 3 | 53 | 6 | 6 | 0 | 0 | 558 | 70 | 102 | 851 |
| 7 | 1 | 4 | 53 | 6 | 6 | 0 | 0 | 474 | 36 | 102 | 851 |
| 7 | 1 | 5 | 53 | 6 | 6 | 0 | 0 | 571 | 58 | 87 | 824 |
| 7 | 1 | 6 | 53 | 6 | 6 | 0 | 0 | 571 | 58 | 102 | 851 |
| 8 | 0 | 0 | 120 | 51 | 51 | 0 | 6 | 868 | 152 | 209 | 1445 |
| 8 | 0 | 1 | 89 | 20 | 20 | 0 | 1 | 879 | 133 | 209 | 1445 |
| 8 | 0 | 2 | 89 | 20 | 20 | 0 | 1 | 868 | 152 | 180 | 1380 |
| 8 | 0 | 3 | 88 | 19 | 19 | 0 | 1 | 868 | 152 | 228 | 1479 |
| 8 | 0 | 4 | 89 | 20 | 20 | 0 | 1 | 765 | 99 | 228 | 1479 |
| 8 | 0 | 5 | 88 | 19 | 19 | 0 | 1 | 879 | 133 | 209 | 1445 |
| 8 | 0 | 6 | 88 | 19 | 19 | 0 | 1 | 879 | 133 | 228 | 1479 |
| 8 | 1 | 0 | 69 | 12 | 12 | 0 | 0 | 819 | 114 | 209 | 1445 |
| 8 | 1 | 1 | 102 | 30 | 30 | 0 | 3 | 851 | 133 | 209 | 1445 |
| 8 | 1 | 2 | 102 | 30 | 30 | 0 | 3 | 819 | 114 | 180 | 1380 |
| 8 | 1 | 3 | 87 | 22 | 22 | 0 | 0 | 819 | 114 | 228 | 1479 |
| 8 | 1 | 4 | 102 | 30 | 30 | 0 | 3 | 996 | 195 | 228 | 1479 |
| 8 | 1 | 5 | 87 | 22 | 22 | 0 | 0 | 851 | 133 | 209 | 1445 |
| 8 | 1 | 6 | 87 | 22 | 22 | 0 | 0 | 851 | 133 | 228 | 1479 |
| 9 | 0 | 0 | 116 | 33 | 35 | 2 | 1 | 1238 | 240 | 436 | 2393 |
| 9 | 0 | 1 | 156 | 63 | 66 | 3 | 7 | 1282 | 272 | 436 | 2393 |
| 9 | 0 | 2 | 156 | 63 | 66 | 3 | 7 | 1238 | 240 | 395 | 2313 |
| 9 | 0 | 3 | 139 | 52 | 54 | 2 | 2 | 1238 | 240 | 459 | 2435 |
| 9 | 0 | 4 | 156 | 63 | 73 | 10 | 7 | 1449 | 354 | 459 | 2435 |
| 9 | 0 | 5 | 139 | 52 | 53 | 1 | 2 | 1282 | 272 | 436 | 2393 |
| 9 | 0 | 6 | 139 | 52 | 55 | 3 | 2 | 1282 | 272 | 459 | 2435 |
| 9 | 1 | 0 | 176 | 93 | 93 | 0 | 10 | 1301 | 294 | 436 | 2393 |
| 9 | 1 | 1 | 139 | 49 | 54 | 5 | 3 | 1310 | 264 | 436 | 2393 |
| 9 | 1 | 2 | 139 | 49 | 49 | 0 | 3 | 1301 | 294 | 395 | 2313 |
| 9 | 1 | 3 | 136 | 46 | 46 | 0 | 3 | 1301 | 294 | 459 | 2435 |
| 9 | 1 | 4 | 139 | 49 | 58 | 9 | 3 | 1176 | 216 | 459 | 2435 |
| 9 | 1 | 5 | 136 | 46 | 52 | 6 | 3 | 1310 | 264 | 436 | 2393 |
| 9 | 1 | 6 | 136 | 46 | 50 | 4 | 3 | 1310 | 264 | 459 | 2435 |
| 10 | 0 | 0 | 255 | 159 | 183 | 24 | 18 | 1890 | 532 | 824 | 3783 |

Totals per (d, parity): d=4 par=0: J=0, (J:A)=0, kernel mod J=0; d=4 par=1: J=0, (J:A)=0, kernel mod J=0; d=5 par=0: J=3, (J:A)=3, kernel mod J=0; d=5 par=1: J=3, (J:A)=3, kernel mod J=0; d=6 par=0: J=18, (J:A)=18, kernel mod J=0; d=6 par=1: J=18, (J:A)=18, kernel mod J=0; d=7 par=0: J=63, (J:A)=63, kernel mod J=0; d=7 par=1: J=63, (J:A)=63, kernel mod J=0; d=8 par=0: J=168, (J:A)=168, kernel mod J=0; d=8 par=1: J=168, (J:A)=168, kernel mod J=0; d=9 par=0: J=378, (J:A)=402, kernel mod J=24; d=9 par=1: J=378, (J:A)=402, kernel mod J=24; d=10 par=0: J=159, (J:A)=183, kernel mod J=24

## Search batches

| family | d | parD | chi | K (dimP) | nP/nQ/rows | restarts (spurious) | best rel | residual histogram | distinct minima (rel, rounded) | rank T at best (nullity) | polish -> rel | s/restart |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hess | 6 | 0 | 0 | 12 (11) | 30/29/436 | 500 (0) | 9.595e-02 | <0.1:148 <0.2:212 <0.5:130 <1:10 | 0.0959, 0.1564, 0.1792, 0.1827, 0.2613, 0.3292 | 0 (30) | 9.59e-02, 9.59e-02, 9.59e-02 | 0 |
| hess | 6 | 0 | 1 | 1 (1) | 21/29/436 | 500 (0) | 9.573e-01 | <1:500 | 0.9573 | 0 (21) | 9.57e-01, 9.57e-01, 9.57e-01 | 0 |
| hess | 6 | 0 | 2 | 1 (1) | 30/19/395 | 500 (0) | 1.000e+00 | <10:500 | 1 | 0 (30) | 1.00e+00, 1.00e+00, 1.00e+00 | 0 |
| full | 6 | 0 | 0 | 12 (11) | 30/29/436 | 600 (0) | 9.595e-02 | <0.1:191 <0.2:240 <0.5:149 <1:20 | 0.0959, 0.1564, 0.1792, 0.1827, 0.2613, 0.3292 | 27 (3) | 9.59e-02, 9.59e-02, 9.59e-02 | 0 |
| hess | 6 | 0 | 3 | 1 (1) | 30/40/459 | 500 (0) | 9.561e-01 | <1:500 | 0.9561 | 0 (30) | 9.56e-01, 9.56e-01, 9.56e-01 | 0 |
| hess | 6 | 0 | 4 | 1 (1) | 9/40/459 | 500 (0) | 9.778e-01 | <1:500 | 0.9778 | 0 (9) | 9.78e-01, 9.78e-01, 9.78e-01 | 0 |
| hess | 6 | 0 | 5 | 1 (1) | 21/29/436 | 500 (0) | 1.000e+00 | <10:500 | 1 | 0 (21) | 1.00e+00, 1.00e+00, 1.00e+00 | 0 |
| hess | 6 | 0 | 6 | 1 (1) | 21/40/459 | 500 (0) | 9.254e-01 | <1:500 | 0.9254 | 0 (21) | 9.25e-01, 9.25e-01, 9.25e-01 | 0 |
| full | 6 | 0 | 1 | 1 (1) | 21/29/436 | 600 (0) | 9.573e-01 | <1:600 | 0.9573 | 18 (3) | 9.57e-01, 9.57e-01, 9.57e-01 | 0 |
| full | 6 | 0 | 2 | 1 (1) | 30/19/395 | 600 (0) | 1.000e+00 | <10:600 | 1 | 27 (3) | 1.00e+00, 1.00e+00, 1.00e+00 | 0 |
| full | 6 | 0 | 3 | 1 (1) | 30/40/459 | 600 (0) | 9.561e-01 | <1:600 | 0.9561 | 27 (3) | 9.56e-01, 9.56e-01, 9.56e-01 | 0 |
| hess | 6 | 1 | 1 | 3 (3) | 20/29/436 | 500 (0) | 9.814e-01 | <1:500 | 0.9814 | 0 (20) | 9.81e-01, 9.81e-01, 9.81e-01 | 0 |
| full | 6 | 0 | 4 | 1 (1) | 9/40/459 | 600 (0) | 9.778e-01 | <1:600 | 0.9778 | 9 (0) | 9.78e-01, 9.78e-01, 9.78e-01 | 0 |
| full | 6 | 0 | 5 | 1 (1) | 21/29/436 | 600 (0) | 1.000e+00 | <10:600 | 1 | 18 (3) | 1.00e+00, 1.00e+00, 1.00e+00 | 0 |
| full | 6 | 0 | 6 | 1 (1) | 21/40/459 | 600 (0) | 9.252e-01 | <1:600 | 0.9252 | 18 (3) | 9.25e-01, 9.25e-01, 9.25e-01 | 0 |
| full | 6 | 1 | 1 | 3 (3) | 20/29/436 | 600 (0) | 9.067e-01 | <1:600 | 0.9067 | 19 (1) | 9.07e-01, 9.07e-01, 9.07e-01 | 0 |
| hess | 6 | 1 | 2 | 3 (3) | 15/19/395 | 500 (0) | 9.956e-01 | <1:500 | 0.9956 | 0 (15) | 9.96e-01, 9.96e-01, 9.96e-01 | 0 |
| hess | 6 | 1 | 3 | 3 (3) | 15/40/459 | 500 (0) | 9.546e-01 | <1:331 <10:169 | 0.9546, 0.9771, 1 | 0 (15) | 9.55e-01, 9.55e-01, 9.55e-01 | 0 |
| full | 6 | 1 | 2 | 3 (3) | 15/19/395 | 600 (0) | 9.956e-01 | <1:600 | 0.9956 | 14 (1) | 9.96e-01, 9.96e-01, 9.96e-01 | 0 |
| full | 6 | 1 | 3 | 3 (3) | 15/40/459 | 600 (0) | 9.541e-01 | <1:415 <10:185 | 0.9541, 0.9692, 1 | 14 (1) | 9.54e-01, 9.54e-01, 9.54e-01 | 0 |
| hess | 6 | 1 | 4 | 3 (3) | 48/40/459 | 500 (0) | 9.825e-01 | <1:500 | 0.9825 | 0 (48) | 9.83e-01, 9.83e-01, 9.83e-01 | 0 |
| hess | 6 | 1 | 5 | 3 (3) | 20/29/436 | 500 (0) | 9.628e-01 | <1:351 <10:149 | 0.9628, 0.9964, 1 | 0 (20) | 9.63e-01, 9.63e-01, 9.63e-01 | 0 |
| hess | 6 | 1 | 6 | 3 (3) | 20/40/459 | 500 (0) | 9.509e-01 | <1:500 | 0.9509, 0.9572, 0.9573 | 0 (20) | 9.51e-01, 9.51e-01, 9.51e-01 | 0 |
| full | 6 | 1 | 4 | 3 (3) | 48/40/459 | 600 (0) | 9.749e-01 | <1:600 | 0.9749, 0.9775 | 36 (12) | 9.75e-01, 9.75e-01, 9.75e-01 | 1 |
| full | 6 | 1 | 5 | 3 (3) | 20/29/436 | 600 (0) | 9.627e-01 | <1:404 <10:196 | 0.9627, 0.9798, 1 | 19 (1) | 9.63e-01, 9.63e-01, 9.63e-01 | 1 |
| hess | 7 | 0 | 0 | 3 (3) | 48/87/824 | 500 (0) | 9.343e-01 | <1:500 | 0.9343 | 0 (48) | 9.34e-01, 9.34e-01, 9.34e-01 | 0 |
| full | 6 | 1 | 6 | 3 (3) | 20/40/459 | 600 (0) | 9.508e-01 | <1:600 | 0.9508, 0.9572, 0.9573 | 19 (1) | 9.51e-01, 9.51e-01, 9.51e-01 | 0 |

Empty / trivial sectors: d=6 par=1 chi=0 (no E1-admissible D in sector); d=6 par=1 chi=0 (no E1-admissible D in sector)


## Other batches (structured families / sweeps / notes)

- {"lens": "structured", "d": 6, "parD": 0, "chi": 0, "K": 12, "n_candidates": 12, "e1_fail": 0, "F_multiples": 0, "n_in_family": 12, "best_rel": 0.6956689911877479, "median_rel": 0.975660744590755, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.2}
- {"lens": "structured", "d": 6, "parD": 0, "chi": 1, "K": 1, "n_candidates": 1, "e1_fail": 0, "F_multiples": 0, "n_in_family": 1, "best_rel": 0.9573112110577245, "median_rel": 0.9573112110577245, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.1}
- {"lens": "structured", "d": 6, "parD": 0, "chi": 2, "K": 1, "n_candidates": 1, "e1_fail": 0, "F_multiples": 0, "n_in_family": 1, "best_rel": 1.0, "median_rel": 1.0, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 12.7}
- {"lens": "structured", "d": 6, "parD": 0, "chi": 3, "K": 1, "n_candidates": 1, "e1_fail": 0, "F_multiples": 0, "n_in_family": 1, "best_rel": 0.95614541696718, "median_rel": 0.95614541696718, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 6.0}
- {"lens": "structured", "d": 6, "parD": 0, "chi": 4, "K": 1, "n_candidates": 1, "e1_fail": 0, "F_multiples": 0, "n_in_family": 1, "best_rel": 0.9777507455658452, "median_rel": 0.9777507455658452, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.5}
- {"lens": "structured", "d": 6, "parD": 0, "chi": 5, "K": 1, "n_candidates": 1, "e1_fail": 0, "F_multiples": 0, "n_in_family": 1, "best_rel": 1.0, "median_rel": 1.0, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.1}
- {"lens": "structured", "d": 6, "parD": 0, "chi": 6, "K": 1, "n_candidates": 1, "e1_fail": 0, "F_multiples": 0, "n_in_family": 1, "best_rel": 0.9251782211080002, "median_rel": 0.9251782211080002, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.5}
- {"lens": "structured", "d": 6, "parD": 1, "chi": 1, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.9067219028513659, "median_rel": 1.0, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.3}
- {"lens": "structured", "d": 6, "parD": 1, "chi": 2, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.9955669193391246, "median_rel": 1.0, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.3}
- {"lens": "structured", "d": 6, "parD": 1, "chi": 3, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.9540540147835398, "median_rel": 0.9691852825426057, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.2}
- {"lens": "structured", "d": 6, "parD": 1, "chi": 4, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.9749365077377793, "median_rel": 1.0, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.4}
- {"lens": "structured", "d": 6, "parD": 1, "chi": 5, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.9627367760499989, "median_rel": 0.9797512855702651, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.2}
- {"lens": "structured", "d": 6, "parD": 1, "chi": 6, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.9507700292485158, "median_rel": 0.9571691202201009, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 1.3}
- {"lens": "structured", "d": 7, "parD": 0, "chi": 0, "K": 3, "n_candidates": 3, "e1_fail": 0, "F_multiples": 0, "n_in_family": 3, "best_rel": 0.7047010672589801, "median_rel": 0.9131922500944042, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 20.0}
- {"lens": "structured", "d": 7, "parD": 0, "chi": 1, "K": 12, "n_candidates": 12, "e1_fail": 0, "F_multiples": 0, "n_in_family": 12, "best_rel": 0.682379205674801, "median_rel": 0.8994011128253099, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 50.9}
- {"lens": "structured", "d": 7, "parD": 0, "chi": 2, "K": 12, "n_candidates": 12, "e1_fail": 0, "F_multiples": 0, "n_in_family": 12, "best_rel": 0.6918877961045277, "median_rel": 0.9721497149425821, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 0.5}
- {"lens": "structured", "d": 7, "parD": 0, "chi": 3, "K": 8, "n_candidates": 8, "e1_fail": 0, "F_multiples": 0, "n_in_family": 8, "best_rel": 0.29762892976746713, "median_rel": 0.929548041096925, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 13.0}
- {"lens": "structured", "d": 7, "parD": 0, "chi": 4, "K": 12, "n_candidates": 12, "e1_fail": 0, "F_multiples": 0, "n_in_family": 12, "best_rel": 0.6912431447210485, "median_rel": 0.922863401852494, "exact_tests": [], "families": ["polar-x-monomial"], "seconds": 18.2}

## Interpretation / caveats

- A local minimum with rel > 1e-6 says nothing rigorous; a rel ~ 1e-1 .. 2e-1 plateau reached from many restarts suggests no solution in that sector at this degree but is not an exclusion.
- Y-even and Y-odd D are both covered; non-eigen D (mixing characters) are not covered by sector-restricted searches.
- Exclusions of degree <= 8 outside J(F) and the Jhat / polar supports are rigorous prior results (dossier, attack1, verify lens); this lens adds numerical evidence only unless a witness is confirmed exactly.
