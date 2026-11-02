# Wave 2 — quadric-contact-residual — SEARCH lens (joint (D,N,Q) search)

STATUS: NO WITNESS. Updated 2026-11-01 20:47 local. Batches logged: 115 (19 varpro sector batches, 1 Hessian-relaxation batches, 6 E1 profiles, 84 structured candidates); random restarts so far: 2363; best full-E2 relative residual: 3.20e-04.

Target: q = 1, A = X0^3 Y0^3 on the quadric-contact sextic F_b. Residue character trivial on G_aff (order 48), odd under s0 (X0<->Y0) and r1 (X1<->Y1). Prior rigorous exclusions (fqlx #125, replicated in wave 1): deg D <= 6; deg 7 in every pure parity/character sector; D8* alone at degree 8; any certificate keeping (Y1/X1, Y2/X2) passive.

Logs: /tmp/hodge/results/quadric-contact-residual/wave2-search.log.jsonl (one JSON line per batch; remote lines merged by w2_report.py). Scripts (lane dir): w2lib.py, w2_varpro.py, w2_e1profile.py, w2_hess.py, w2_struct.py, w2_report.py. Remote job logs: /tmp/hodge/remote/_jobs/search-quadric-contact-residual-*.log

## Method (what is new in wave 2)

* Exact E1 families per G_aff character sector: (J:A)_d[chi] = {D in S_d[chi] : A D in J_{d+6}} over Q (python-flint integer nullspace), split into the Jacobian part J_d[chi] (sparse generators m*F_i) and the non-Jacobian extras (exact integer polynomials).
* Koszul-form variable projection: D(a) = sum a_k D_k with E1-lifts N_k (exact A m e_i for generators, least-squares lifts for extras); every E1 solution is N = sum a_k N_k + K with K = sum_{r<s} Q_rs (F_s e_r - F_r e_s), and the certificate reduces (q = 1) to B_D(N) = F Q. f(a) = min over (Q_rs, Q) of ||B_{D(a)}(sum a_k N_k + K) - F Q||^2 on the unit sphere (L-BFGS, exact Golub-Pereyra gradient, Q-block eliminated by orthogonal projection, lstsq gelsy with cutoff 1e-9). Spurious directions (D(a) = 0, D(a) in F*S_{d-6}) are projected out (monomial-orthogonality to F-multiples). rel = f/|b|^2 with b = -sum a_k a_l B_{D_k}(N_l).
* Hessian relaxation for character-MIXED families: T(a,a) = NF_{J_{2d}}(sum a_k a_l B_{D_k}(N_l)) = 0 (necessary; fqlx's Hessian congruence), minimized numerically with the same restart machinery; near-zeros are polished with the full Koszul system.
* Every minimizer with rel < 1e-10 is rationally reconstructed in monomial coordinates (denominators up to 1e6) and re-solved EXACTLY (engine symmetry.solve_fixed_D_symmetric / fixed_d solvers, core.verify_witness, coprimality by univariate specialization). Only an exact re-solve counts.

## Jobs (remote hodge-sextics, 128 cores shared with ~50 agents)

* `excl-reverify` (submitted 2026-11-01T19:28:24): `cd wave2-excl && HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-reverify.log.jsonl python3 reverify.py 1 2 3 4 5`
* `search-quadric-contact-residual-pipeB` (submitted 2026-11-01T19:50:19): `bash w2_pipeline_B.sh`
* `search-quadric-contact-residual-pipeA2` (submitted 2026-11-01T20:29:50): `bash w2_pipeline_A.sh`
  - search-quadric-contact-residual-e1prof.log: last line: `{"task": "e1-profile", "d": 10, "n_sectors": 24, "tot_dimS": 3003, "tot_dimJ": 1497, "tot_kernel": 2053, "tot_extra": 556, "tot_gauge": 126, "per_sector_extra":`
  - search-quadric-contact-residual-pipeA.log: last line: `[koszul-varpro] d=9 chi=0011 K=41 rows=1356 cols=427 (K 187, Q 240) Mstack nnz=275283 Vstack nnz=12423 dimP=38 excluded=3 built 1.2s`
  - search-quadric-contact-residual-pipeA2.log: last line: `[build_family] {'d': 8, 'chi': '1050', 'parity': None, 'dimS': 57, 'rows': 498, 'cols': 597, 'nullity': 118, 'dim_kernel': 18, 'dimJ': 18, 'n_extra': 0, 'dim_ga`
  - search-quadric-contact-residual-pipeB.log: last line: `=== 2026-11-02 00:47:52 w2_hess.py --d 7 --parity 0 --restarts 1500 --procs 8 --seed 11 --tag hd7e`
  - search-quadric-contact-residual-pipeB2.log: last line: `=== 2026-11-02 01:29:56 w2_hess.py --d 7 --parity 0 --restarts 1500 --procs 8 --seed 11 --tag hd7e`
  - search-quadric-contact-residual-pipeC.log: last line: `E1 FAIL d=8 chi=(0, 0, 0, 0) X0*X2**5*Y0*Y2 + X0*X2*Y0*Y2**5 {'rank': `
  - search-quadric-contact-residual-vp-d8.log: last line: `[koszul-varpro] d=8 chi=0000 K=33 rows=1002 cols=303 (K 130, Q 173) Mstack nnz=155331 Vstack nnz=43948 dimP=29 excluded=4 built 0.7s`

## (1) E1-admissible denominator spaces per degree and character sector (exact over Q)

| d | sectors done | dim S_d | dim J_d | dim (J:A)_d | non-Jacobian | gauge F*S_{d-6} | sectors with extras (chi:count) |
|---|---|---|---|---|---|---|---|
| 4 | 24 | 126 | 0 | 0 | 0 | 0 |  |
| 5 | 24 | 252 | 6 | 6 | 0 | 0 |  |
| 6 | 24 | 462 | 36 | 36 | 0 | 1 |  |
| 7 | 24 | 792 | 126 | 126 | 0 | 6 |  |
| 8 | 24 | 1287 | 336 | 337 | 1 | 21 | 0000:1 |
| 9 | 24 | 2002 | 756 | 760 | 4 | 56 | 0010:1, 0051:1, 0100:1, 0101:1 |
| 10 | 24 | 3003 | 1497 | 2053 | 556 | 126 | 0000:20, 0001:16, 0020:35, 0021:21, 0040:35, 0041:21, 0110:19, 0111:19, 0130:30, 0131:30, 0150:19, 0151:19, 1010:19, 1011:19, 1030:30, 1031:30, 1050:19, 1051:19, 1100:16, 1101:16, 1120:26, 1121:26, 1140:26, 1141:26 |

Cross-checks: d <= 9 totals agree with fqlx (#125) and wave 1 (dim (J:A)_d = 0,0,0,0,0,6,36,126,337,760,2053,...); the per-sector split is new.
Per-sector detail (d: chi -> dimS/J/kernel/extra) for degrees with non-Jacobian directions:
* d=8: 0000: 87/32/33/1
* d=9: 0010: 95/41/42/1; 0051: 95/41/42/1; 0100: 117/58/59/1; 0101: 117/58/59/1
* d=10: 0000: 180/108/128/20; 0001: 173/103/119/16; 0020: 120/48/83/35; 0021: 113/49/70/21; 0040: 120/48/83/35; 0041: 113/49/70/21; 0110: 131/73/92/19; 0111: 131/73/92/19; 0130: 102/36/66/30; 0131: 102/36/66/30; 0150: 131/73/92/19; 0151: 131/73/92/19; 1010: 131/73/92/19; 1011: 131/73/92/19; 1030: 102/36/66/30; 1031: 102/36/66/30; 1050: 131/73/92/19; 1051: 131/73/92/19; 1100: 160/98/114/16; 1101: 160/98/114/16; 1120: 102/42/68/26; 1121: 102/42/68/26; 1140: 102/42/68/26; 1141: 102/42/68/26

### Cone / coordinate-subspace families (D supported on a coordinate variable subset, per sector)

* d=6: 0 (subset, sector) pairs with an E1-kernel beyond the gauge
* d=7: 0 (subset, sector) pairs with an E1-kernel beyond the gauge
* d=8: 4 (subset, sector) pairs with an E1-kernel beyond the gauge: X0X1X2Y0Y1:0000 ker=4 J=4 extra=0 gauge=3, X0X1X2Y0Y2:0000 ker=5 J=5 extra=0 gauge=4, X0X1Y0Y1Y2:0000 ker=4 J=4 extra=0 gauge=3, X0X2Y0Y1Y2:0000 ker=5 J=5 extra=0 gauge=4

## (2)-(3) Full-E2 variable-projection searches (Koszul form), per degree and sector

| tag | d | chi | vars | members K | dimP | rows | nK+nQ | t_eval s | restarts | spurious | best rel | median rel | n<1e-6 | n<1e-10 | mean iters | seconds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| w2b | 8 | 0000 | - | 33 | 29 | 1002 | 130+173 | 0.3123 | 9 | 0 | 0.001466 | 0.05118 | 0 | 0 | 119.3 | 1858.5 |
| w2b | 8 | 0000 | - | 33 | 29 | 1002 | 130+173 | 0.0284 | 290 | 0 | 0.001436 | 0.001657 | 0 | 0 | 106.8 | 671.8 |
| w2b | 8 | 0001 | - | 30 | 27 | 1000 | 124+173 | 0.0211 | 270 | 0 | 0.006748 | 0.006748 | 0 | 0 | 97.05 | 341.6 |
| w2b | 8 | 0020 | - | 8 | 7 | 812 | 55+113 | 0.0075 | 70 | 0 | 0.1252 | 0.268 | 0 | 0 | 23.06 | 10.7 |
| w2b | 8 | 0021 | - | 9 | 9 | 813 | 49+113 | 0.0047 | 90 | 0 | 0.03661 | 0.03681 | 0 | 0 | 31.94 | 14.3 |
| w2b | 8 | 0040 | - | 8 | 7 | 812 | 55+113 | 0.007 | 70 | 0 | 0.1252 | 0.268 | 0 | 0 | 23.67 | 10.2 |
| w2b | 8 | 0041 | - | 9 | 9 | 813 | 49+113 | 0.0053 | 90 | 0 | 0.03661 | 0.03681 | 0 | 0 | 30.89 | 13.4 |
| w2b | 8 | 0110 | - | 18 | 17 | 812 | 100+113 | 0.0118 | 170 | 0 | 0.03388 | 0.03388 | 0 | 0 | 45.75 | 68.0 |
| w2b | 8 | 0111 | - | 18 | 17 | 812 | 100+113 | 0.0092 | 170 | 0 | 0.03382 | 0.03382 | 0 | 0 | 43.94 | 53.5 |
| w2b | 8 | 0130 | - | 4 | 4 | 1002 | 38+173 | 0.0043 | 40 | 0 | 0.3298 | 0.3315 | 0 | 0 | 14.3 | 4.1 |
| w2b | 8 | 0131 | - | 4 | 4 | 1002 | 38+173 | 0.0048 | 40 | 0 | 0.3298 | 0.3315 | 0 | 0 | 14.53 | 4.2 |
| w2b | 8 | 0150 | - | 18 | 17 | 812 | 100+113 | 0.0112 | 170 | 0 | 0.03382 | 0.03382 | 0 | 0 | 44.58 | 53.8 |
| w2b | 8 | 0151 | - | 18 | 17 | 812 | 100+113 | 0.0118 | 170 | 0 | 0.03388 | 0.03388 | 0 | 0 | 45.68 | 61.8 |
| w2b | 8 | 1010 | - | 18 | 17 | 812 | 100+113 | 0.0097 | 170 | 0 | 0.01338 | 0.06012 | 0 | 0 | 59.15 | 66.4 |
| w2b | 8 | 1011 | - | 18 | 17 | 812 | 100+113 | 0.0085 | 170 | 0 | 0.01338 | 0.06012 | 0 | 0 | 59.85 | 66.7 |
| w2b | 8 | 1030 | - | 4 | 4 | 1002 | 38+173 | 0.0049 | 40 | 0 | 0.168 | 0.168 | 0 | 0 | 20.93 | 5.1 |
| w2b | 8 | 1031 | - | 4 | 4 | 1002 | 38+173 | 0.0046 | 40 | 0 | 0.168 | 0.168 | 0 | 0 | 21.77 | 4.7 |
| w2b | 8 | 1050 | - | 18 | 17 | 812 | 100+113 | 0.011 | 170 | 0 | 0.01338 | 0.06012 | 0 | 0 | 58.69 | 73.1 |
| w2a | 9 | 0010 | - | 42 | 38 | 1356 | 190+240 | 0.0374 | 120 | 0 | 3.20e-04 | 4.00e-04 | 0 | 0 | 155.2 | 1135.3 |

Reconstruction attempts (rel < 1e-10): 0

## Hessian relaxation (character-mixed families, necessary condition mod J_{2d})

| tag | d | mode | members | dimP | dim R_{2d} | restarts | best |T|^2/scale | median | n<1e-8 | n<1e-12 | n<1e-16 | polished (full rel) | seconds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| smoke | 5 | par0 | 3 | 3 | 696 | 4 | 0.5164 | 0.5164 | 0 | 0 | 0 | - | 4.8 |

## (4) Structured candidates (geometry-suggested denominators; E1 screen + exact sector solve mod p)

| d | candidate | chi | E1 | consistent | rank | nullity | rows x cols | s |
|---|---|---|---|---|---|---|---|---|
| 8 | D8r (r1-contour) | None | False | None | 5052 | 954 | NonexNone | 164.5 |
| 8 | pol_r1+ * pol_r1- | None | False | None | 5052 | 954 | NonexNone | 153.3 |
| 8 | minor(X0,Y0) | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | minor(X1,Y1) | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | minor(X2,Y2) | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | minor(X0,X1) | 0040 | False | None | 395 | 55 | NonexNone | 0.3 |
| 8 | minor(X0,Y1) | 0020 | False | None | 395 | 55 | NonexNone | 0.4 |
| 8 | minor(Y0,X1) | 0040 | False | None | 395 | 55 | NonexNone | 0.4 |
| 8 | minor(X0,X2) | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | minor(X1,X2) | 0040 | False | None | 395 | 55 | NonexNone | 0.4 |
| 8 | rad0*S0 | 0000 | True | False | 718 | 109 | 3703x827 | 4.6 |
| 8 | rad0*delta | 0000 | True | False | 718 | 109 | 3703x827 | 4.5 |
| 8 | rad0*X0^2 | 0000 | True | False | 718 | 109 | 3703x827 | 4.1 |
| 8 | rad1*S1 | None | True | False | 6693 | 804 | 40310x7497 | 421.4 |
| 8 | rad1*d1 | None | True | False | 6693 | 804 | 40310x7497 | 426.1 |
| 8 | rad2*S2 | 0000 | True | False | 718 | 109 | 3716x827 | 4.5 |
| 8 | rad1*S0 | 0000 | True | False | 717 | 110 | 3703x827 | 4.4 |
| 8 | rad0*S1 | None | True | False | 6703 | 794 | 40262x7497 | 409.6 |
| 8 | rad0*S2 | 0000 | True | False | 718 | 109 | 3704x827 | 4.0 |
| 8 | rad0*P1 | 0001 | True | False | 716 | 115 | 3704x831 | 3.8 |
| 8 | GS*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | GS*delta^2 | 0000 | True | False | 713 | 114 | 3703x827 | 5.1 |
| 8 | GS*P0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | GP*P0*S0 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | GS*Q1 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | D8* | 0000 | True | False | 718 | 109 | 3720x827 | 3.5 |
| 8 | GP*GS | 0001 | False | None | 534 | 124 | NonexNone | 0.3 |
| 8 | GP*P0P1 | 0001 | False | None | 534 | 124 | NonexNone | 0.3 |
| 8 | GP*P0P2 | 0001 | False | None | 534 | 124 | NonexNone | 0.3 |
| 8 | GP*P0^2 | 0001 | False | None | 534 | 124 | NonexNone | 0.3 |
| 8 | GP*P1^2 | 0001 | False | None | 534 | 124 | NonexNone | 0.3 |
| 8 | GP*P2^2 | 0001 | False | None | 534 | 124 | NonexNone | 0.4 |
| 8 | GP*Q1 | 0001 | False | None | 534 | 124 | NonexNone | 0.4 |
| 8 | GP*Q2 | 0000 | False | None | 524 | 130 | NonexNone | 0.4 |
| 8 | GP*S0S1 | None | False | None | 5070 | 936 | NonexNone | 100.1 |
| 8 | GP*S0S2 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | GP*S0^2 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | GP*S1S2 | None | False | None | 5070 | 936 | NonexNone | 156.2 |
| 8 | GP*X0^4+Y0^4 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | GP*X1^4+Y1^4 | None | False | None | 5070 | 936 | NonexNone | 155.3 |
| 8 | GP*X2^4+Y2^4 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | GP*delta^2 | 0001 | True | False | 708 | 123 | 3703x831 | 4.1 |
| 8 | GS*P0P1 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | GS*P0P2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | GS*P1^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | GS*P2^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | GS*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | GS*S0S1 | None | False | None | 5052 | 954 | NonexNone | 106.7 |
| 8 | GS*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.4 |
| 8 | GS*S1S2 | None | False | None | 5052 | 954 | NonexNone | 110.1 |
| 8 | GS*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.5 |
| 8 | GS*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 124.4 |
| 8 | GS*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P1*P0P2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P1*P0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P1*P1^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P1*P2^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P1*Q1 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P1*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.9 |
| 8 | P0P1*S0S1 | None | False | None | 5052 | 954 | NonexNone | 190.6 |
| 8 | P0P1*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P1*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P1*S1S2 | None | False | None | 5052 | 954 | NonexNone | 177.0 |
| 8 | P0P1*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P1*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 178.1 |
| 8 | P0P1*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P1*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P2*P0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P2*P1^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P2*P2^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P2*Q1 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | P0P2*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.8 |
| 8 | P0P2*S0S1 | None | False | None | 5052 | 954 | NonexNone | 131.7 |
| 8 | P0P2*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P2*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P2*S1S2 | None | False | None | 5052 | 954 | NonexNone | 167.6 |
| 8 | P0P2*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0P2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 142.9 |
| 8 | P0P2*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0P2*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0^2*P1^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0^2*P2^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0^2*Q1 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0^2*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.7 |

## Best residual landscape (full-E2 rel, log10 histogram of local minima over all restarts)

* d=8: 1e-3: 452, 1e-2: 1164, 1e-1: 623
* d=9: 1e-4: 90, 1e-2: 30

Reading: a genuine witness in a family shows up as rel ~ 1e-14..1e-12 (float64 floor); the smallest values seen so far are listed above; all minima are separated from the floor by many orders of magnitude unless stated otherwise.

## Structure notes (from the exact E1 tables)

* (J:A)/J has exactly ONE minimal generator beyond J in degrees <= 9: D8* = 4G_S^2 - G_P^2 (= const * G_u G_v, the product of the two reflection-contour quartics); the four degree-9 directions are D8*.X1, D8*.Y1, D8*.X2, D8*.Y2 (D8*.X0, D8*.Y0 lie in J_9). At degree 10 only 10 of the 556 non-Jacobian directions come from D8*.S_2; the other 546 are new minimal generators of the colon ideal (computed mod p this session). From d = 10 on E1 is essentially the Gorenstein-dual constraint codim (J:A)_d = dim R_{d+6} (minus 1 at d = 10): 950 = 951 - 1, 666, 426 for d = 10, 11, 12.
* Gorenstein duality (socle degree 24) pairs ker(.A: R_d -> R_{d+6}) with coker(.A: R_{18-d} -> R_{24-d}); the single kernel direction D8* at d = 8 is dual to the single missing direction of A.R_10 in R_16, consistent with the 556 = 1506 - 950 kernel at d = 10.
* Sector bookkeeping: G_aff characters chi = (e0+e3, e2+e5 mod 2; e1-e4 mod 6; Y-parity); the target residue is G_aff-invariant and the identity for a chi-eigen D lives in sector chi_A + 2 chi, so restricting all unknowns to sectors is lossless for eigen-denominators (24 sectors per degree parity).

## Errors / notes

* Engine: no engine files were modified; lane-local library w2lib.py. Observations about the shared engine are in /tmp/hodge/notes/engine-bugs.md.
