# Wave 2 — quadric-contact-residual — SEARCH lens (joint (D,N,Q) search)

STATUS: NO WITNESS. Updated 2026-11-01 22:39 local. Batches logged: 213 (38 varpro sector batches, 2 Hessian-relaxation batches, 6 E1 profiles, 161 structured candidates); random restarts so far: 6145; best full-E2 relative residual: 6.34e-07.

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
* `excl-misc` (submitted 2026-11-01T21:24:07): `cd wave2-excl && bash misc.sh`
* `excl-lattice7` (submitted 2026-11-01T21:26:18): `cd wave2-excl && bash lattice7.sh`
  - search-quadric-contact-residual-e1prof.log: last line: `{"task": "e1-profile", "d": 11, "n_sectors": 24, "tot_dimS": 4368, "tot_dimJ": 2682, "tot_kernel": 3702, "tot_extra": 1020, "tot_gauge": 252, "per_sector_extra"`
  - search-quadric-contact-residual-pipeA.log: last line: `[koszul-varpro] d=9 chi=0011 K=41 rows=1356 cols=427 (K 187, Q 240) Mstack nnz=275283 Vstack nnz=12423 dimP=38 excluded=3 built 1.2s`
  - search-quadric-contact-residual-pipeA2.log: last line: `{"task": "varpro", "`
  - search-quadric-contact-residual-pipeB.log: last line: `=== 2026-11-02 00:47:52 w2_hess.py --d 7 --parity 0 --restarts 1500 --procs 8 --seed 11 --tag hd7e`
  - search-quadric-contact-residual-pipeB2.log: last line: `=== 2026-11-02 03:00:57 w2_hess.py --d 7 --parity 1 --restarts 1500 --procs 8 --seed 12 --tag hd7o`
  - search-quadric-contact-residual-pipeC.log: last line: `=== 2026-11-02 03:19:45 w2_struct.py 10`
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
| 11 | 24 | 4368 | 2682 | 3702 | 1020 | 252 | 0010:36, 0011:38, 0030:56, 0031:56, 0050:38, 0051:36, 0100:32, 0101:32, 0120:49, 0121:49, 0140:49, 0141:49, 1000:32, 1001:32, 1020:49, 1021:49, 1040:49, 1041:49, 1110:34, 1111:34, 1130:52, 1131:52, 1150:34, 1151:34 |

Cross-checks: d <= 9 totals agree with fqlx (#125) and wave 1 (dim (J:A)_d = 0,0,0,0,0,6,36,126,337,760,2053,...); the per-sector split is new.
Per-sector detail (d: chi -> dimS/J/kernel/extra) for degrees with non-Jacobian directions:
* d=8: 0000: 87/32/33/1
* d=9: 0010: 95/41/42/1; 0051: 95/41/42/1; 0100: 117/58/59/1; 0101: 117/58/59/1
* d=10: 0000: 180/108/128/20; 0001: 173/103/119/16; 0020: 120/48/83/35; 0021: 113/49/70/21; 0040: 120/48/83/35; 0041: 113/49/70/21; 0110: 131/73/92/19; 0111: 131/73/92/19; 0130: 102/36/66/30; 0131: 102/36/66/30; 0150: 131/73/92/19; 0151: 131/73/92/19; 1010: 131/73/92/19; 1011: 131/73/92/19; 1030: 102/36/66/30; 1031: 102/36/66/30; 1050: 131/73/92/19; 1051: 131/73/92/19; 1100: 160/98/114/16; 1101: 160/98/114/16; 1120: 102/42/68/26; 1121: 102/42/68/26; 1140: 102/42/68/26; 1141: 102/42/68/26
* d=11: 0010: 200/132/168/36; 0011: 198/130/168/38; 0030: 162/81/137/56; 0031: 162/81/137/56; 0050: 198/130/168/38; 0051: 200/132/168/36; 0100: 236/165/197/32; 0101: 236/165/197/32; 0120: 162/89/138/49; 0121: 162/89/138/49; 0140: 162/89/138/49; 0141: 162/89/138/49; 1000: 236/165/197/32; 1001: 236/165/197/32; 1020: 162/89/138/49; 1021: 162/89/138/49; 1040: 162/89/138/49; 1041: 162/89/138/49; 1110: 180/120/154/34; 1111: 180/120/154/34; 1130: 144/72/124/52; 1131: 144/72/124/52; 1150: 180/120/154/34; 1151: 180/120/154/34

### Cone / coordinate-subspace families (D supported on a coordinate variable subset, per sector)

* d=6: 0 (subset, sector) pairs with an E1-kernel beyond the gauge
* d=7: 0 (subset, sector) pairs with an E1-kernel beyond the gauge
* d=8: 4 (subset, sector) pairs with an E1-kernel beyond the gauge: X0X1X2Y0Y1:0000 ker=4 J=4 extra=0 gauge=3, X0X1X2Y0Y2:0000 ker=5 J=5 extra=0 gauge=4, X0X1Y0Y1Y2:0000 ker=4 J=4 extra=0 gauge=3, X0X2Y0Y1Y2:0000 ker=5 J=5 extra=0 gauge=4
* d=10: 81 (subset, sector) pairs with an E1-kernel beyond the gauge: X0X2Y0Y2:0000 ker=24 J=12 extra=12 gauge=11, X0X2Y0Y2:0001 ker=16 J=10 extra=6 gauge=8, X0X2Y0Y2:1100 ker=24 J=10 extra=14 gauge=8, X0X2Y0Y2:1101 ker=24 J=10 extra=14 gauge=8, X1X2Y1Y2:0001 ker=5 J=5 extra=0 gauge=4, X0X1X2Y0Y1:0000 ker=14 J=12 extra=2 gauge=8, X0X1X2Y0Y1:0001 ker=10 J=10 extra=0 gauge=6, X0X1X2Y0Y1:0020 ker=12 J=5 extra=7 gauge=4, X0X1X2Y0Y1:0021 ker=3 J=2 extra=1 gauge=2, X0X1X2Y0Y1:0040 ker=12 J=5 extra=7 gauge=4, X0X1X2Y0Y1:0041 ker=3 J=2 extra=1 gauge=2, X0X1X2Y0Y1:0110 ker=4 J=4 extra=0 gauge=3, X0X1X2Y0Y1:0151 ker=4 J=4 extra=0 gauge=3, X0X1X2Y0Y1:1010 ker=10 J=6 extra=4 gauge=4, X0X1X2Y0Y1:1011 ker=10 J=6 extra=4 gauge=4, X0X1X2Y0Y1:1030 ker=6 J=2 extra=4 gauge=2, X0X1X2Y0Y1:1031 ker=6 J=2 extra=4 gauge=2, X0X1X2Y0Y1:1050 ker=10 J=6 extra=4 gauge=4, X0X1X2Y0Y1:1051 ker=10 J=6 extra=4 gauge=4, X0X1X2Y0Y1:1100 ker=6 J=6 extra=0 gauge=4, X0X1X2Y0Y1:1101 ker=6 J=6 extra=0 gauge=4, X0X1X2Y0Y2:0000 ker=35 J=23 extra=12 gauge=11, X0X1X2Y0Y2:0001 ker=24 J=18 extra=6 gauge=8, X0X1X2Y0Y2:0020 ker=12 J=5 extra=7 gauge=4, X0X1X2Y0Y2:0021 ker=3 J=2 extra=1 gauge=2, X0X1X2Y0Y2:0110 ker=21 J=8 extra=13 gauge=5, X0X1X2Y0Y2:0111 ker=21 J=8 extra=13 gauge=5, X0X1X2Y0Y2:1010 ker=21 J=7 extra=14 gauge=5, X0X1X2Y0Y2:1011 ker=21 J=7 extra=14 gauge=5, X0X1X2Y0Y2:1100 ker=32 J=18 extra=14 gauge=8, X0X1X2Y0Y2:1101 ker=32 J=18 extra=14 gauge=8, X0X1X2Y0Y2:1120 ker=8 J=2 extra=6 gauge=2, X0X1X2Y0Y2:1121 ker=8 J=2 extra=6 gauge=2, X0X1X2Y1Y2:0000 ker=11 J=8 extra=3 gauge=8, X0X1X2Y1Y2:0001 ker=7 J=7 extra=0 gauge=6, X0X1X2Y1Y2:0020 ker=12 J=4 extra=8 gauge=4, X0X1X2Y1Y2:0040 ker=12 J=4 extra=8 gauge=4, X0X1X2Y1Y2:0110 ker=10 J=4 extra=6 gauge=4, X0X1X2Y1Y2:0111 ker=10 J=4 extra=6 gauge=4, X0X1X2Y1Y2:0130 ker=6 J=2 extra=4 gauge=2

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
| w2b | 8 | 1051 | - | 18 | 17 | 812 | 100+113 | 0.0096 | 170 | 0 | 0.01338 | 0.06012 | 0 | 0 | 59.85 | 85.4 |
| w2b | 8 | 1100 | - | 28 | 26 | 1000 | 130+173 | 0.0157 | 260 | 0 | 0.01311 | 0.05618 | 0 | 0 | 90.05 | 287.7 |
| w2b | 8 | 1101 | - | 28 | 26 | 1000 | 130+173 | 0.015 | 260 | 0 | 0.01311 | 0.05765 | 0 | 0 | 89.12 | 311.2 |
| w2b | 8 | 1120 | - | 6 | 6 | 813 | 54+113 | 0.0055 | 60 | 0 | 0.2902 | 0.3161 | 0 | 0 | 18.68 | 6.2 |
| w2b | 8 | 1121 | - | 6 | 6 | 813 | 54+113 | 0.0067 | 60 | 0 | 0.2902 | 0.3161 | 0 | 0 | 17.62 | 4.5 |
| w2b | 8 | 1140 | - | 6 | 6 | 813 | 54+113 | 0.0054 | 60 | 0 | 0.2902 | 0.3161 | 0 | 0 | 18.27 | 5.4 |
| w2b | 8 | 1141 | - | 6 | 6 | 813 | 54+113 | 0.0055 | 60 | 0 | 0.2902 | 0.2905 | 0 | 0 | 17.93 | 4.5 |
| w2a | 9 | 0010 | - | 42 | 38 | 1356 | 190+240 | 0.0374 | 120 | 0 | 3.20e-04 | 4.00e-04 | 0 | 0 | 155.2 | 1135.3 |
| w2a | 9 | 0010 | - | 42 | 38 | 1356 | 190+240 | 0.0452 | 120 | 0 | 3.20e-04 | 4.00e-04 | 0 | 0 | 155.2 | 455.0 |
| w2a | 9 | 0011 | - | 41 | 38 | 1356 | 187+240 | 0.0471 | 120 | 0 | 0.001484 | 0.001487 | 0 | 0 | 97.88 | 392.4 |
| w2a | 9 | 0030 | - | 16 | 15 | 1611 | 93+332 | 0.0402 | 90 | 0 | 0.0238 | 0.0238 | 0 | 0 | 46.14 | 110.9 |
| w2a | 9 | 0031 | - | 16 | 15 | 1611 | 93+332 | 0.0242 | 90 | 0 | 0.0238 | 0.0238 | 0 | 0 | 44.94 | 120.8 |
| w2a | 9 | 0050 | - | 41 | 38 | 1356 | 187+240 | 0.0482 | 120 | 0 | 0.001003 | 0.001487 | 0 | 0 | 96.23 | 431.6 |
| w2a | 9 | 0051 | - | 42 | 38 | 1356 | 190+240 | 0.063 | 120 | 0 | 1.33e-04 | 1.36e-04 | 0 | 0 | 164.7 | 541.4 |
| w2a | 9 | 0100 | - | 59 | 53 | 1611 | 234+332 | 0.113 | 120 | 0 | 6.34e-07 | 6.37e-07 | 120 | 0 | 192.8 | 1058.5 |
| w2a | 9 | 0101 | - | 59 | 52 | 1611 | 234+332 | 0.0905 | 116 | 0 | 5.58e-04 | 0.008087 | 0 | 0 | 227.8 | 1503.8 |
| w2a | 9 | 0120 | - | 20 | 19 | 1356 | 118+240 | 0.0305 | 114 | 0 | 0.0387 | 0.0387 | 0 | 0 | 44.75 | 150.7 |
| w2a | 9 | 0121 | - | 20 | 19 | 1356 | 118+240 | 0.0346 | 114 | 0 | 0.03835 | 0.03835 | 0 | 0 | 45.7 | 124.1 |
| w2a | 9 | 0140 | - | 20 | 19 | 1356 | 118+240 | 0.0256 | 114 | 0 | 0.0387 | 0.0387 | 0 | 0 | 43.71 | 121.7 |
| w2a | 9 | 0141 | - | 20 | 19 | 1356 | 118+240 | 0.0266 | 114 | 0 | 0.03835 | 0.03835 | 0 | 0 | 45.04 | 114.8 |

Reconstruction attempts (rel < 1e-10): 0

## Hessian relaxation (character-mixed families, necessary condition mod J_{2d})

| tag | d | mode | members | dimP | dim R_{2d} | restarts | best |T|^2/scale | median | n<1e-8 | n<1e-12 | n<1e-16 | polished (full rel) | seconds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| smoke | 5 | par0 | 3 | 3 | 696 | 4 | 0.5164 | 0.5164 | 0 | 0 | 0 | - | 4.8 |
| hd7e | 7 | par0 | 63 | 60 | 740 | 1500 | 4.64e-04 | 4.64e-04 | 0 | 0 | 0 | - | 5460.1 |

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
| 8 | P0^2*S0S1 | None | False | None | 5052 | 954 | NonexNone | 145.8 |
| 8 | P0^2*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0^2*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | P0^2*S1S2 | None | False | None | 5052 | 954 | NonexNone | 144.7 |
| 8 | P0^2*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0^2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 148.4 |
| 8 | P0^2*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P0^2*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*P2^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*Q1 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.7 |
| 8 | P1^2*S0S1 | None | False | None | 5052 | 954 | NonexNone | 151.1 |
| 8 | P1^2*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*S1S2 | None | False | None | 5052 | 954 | NonexNone | 156.5 |
| 8 | P1^2*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 156.0 |
| 8 | P1^2*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P1^2*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P2^2*Q1 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P2^2*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.7 |
| 8 | P2^2*S0S1 | None | False | None | 5052 | 954 | NonexNone | 156.7 |
| 8 | P2^2*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P2^2*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P2^2*S1S2 | None | False | None | 5052 | 954 | NonexNone | 152.8 |
| 8 | P2^2*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P2^2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 151.3 |
| 8 | P2^2*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | P2^2*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | Q1*Q2 | 0001 | False | None | 534 | 124 | NonexNone | 0.7 |
| 8 | Q1*S0S1 | None | False | None | 5052 | 954 | NonexNone | 149.8 |
| 8 | Q1*S0S2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | Q1*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | Q1*S1S2 | None | False | None | 5052 | 954 | NonexNone | 150.5 |
| 8 | Q1*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | Q1*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 151.3 |
| 8 | Q1*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | Q1*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.4 |
| 8 | Q2*S0S1 | None | False | None | 5070 | 936 | NonexNone | 151.5 |
| 8 | Q2*S0S2 | 0001 | False | None | 534 | 124 | NonexNone | 0.7 |
| 8 | Q2*S0^2 | 0001 | False | None | 534 | 124 | NonexNone | 0.7 |
| 8 | Q2*S1S2 | None | False | None | 5070 | 936 | NonexNone | 173.6 |
| 8 | Q2*X0^4+Y0^4 | 0001 | False | None | 534 | 124 | NonexNone | 0.8 |
| 8 | Q2*X1^4+Y1^4 | None | False | None | 5070 | 936 | NonexNone | 177.3 |
| 8 | Q2*X2^4+Y2^4 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | Q2*delta^2 | 0001 | False | None | 534 | 124 | NonexNone | 0.6 |
| 8 | S0S1*S0S2 | None | False | None | 5052 | 954 | NonexNone | 156.9 |
| 8 | S0S1*S0^2 | None | False | None | 5052 | 954 | NonexNone | 166.7 |
| 8 | S0S1*S1S2 | None | False | None | 5052 | 954 | NonexNone | 161.5 |
| 8 | S0S1*X0^4+Y0^4 | None | False | None | 5052 | 954 | NonexNone | 151.7 |
| 8 | S0S1*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 155.5 |
| 8 | S0S1*X2^4+Y2^4 | None | False | None | 5052 | 954 | NonexNone | 162.4 |
| 8 | S0S1*delta^2 | None | False | None | 5052 | 954 | NonexNone | 147.1 |
| 8 | S0S2*S0^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | S0S2*S1S2 | None | False | None | 5052 | 954 | NonexNone | 144.8 |
| 8 | S0S2*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.6 |
| 8 | S0S2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 155.6 |
| 8 | S0S2*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | S0S2*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | S0^2*S1S2 | None | False | None | 5052 | 954 | NonexNone | 176.0 |
| 8 | S0^2*X0^4+Y0^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | S0^2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 184.2 |
| 8 | S0^2*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | S0^2*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | S1S2*X0^4+Y0^4 | None | False | None | 5052 | 954 | NonexNone | 178.3 |
| 8 | S1S2*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 181.7 |
| 8 | S1S2*X2^4+Y2^4 | None | False | None | 5052 | 954 | NonexNone | 175.4 |
| 8 | S1S2*delta^2 | None | False | None | 5052 | 954 | NonexNone | 188.4 |
| 8 | X0^4+Y0^4*X1^4+Y1^4 | None | False | None | 5052 | 954 | NonexNone | 185.4 |
| 8 | X0^4+Y0^4*X2^4+Y2^4 | 0000 | False | None | 524 | 130 | NonexNone | 0.7 |
| 8 | X0^4+Y0^4*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | X1^4+Y1^4*X2^4+Y2^4 | None | False | None | 5052 | 954 | NonexNone | 185.4 |
| 8 | X1^4+Y1^4*delta^2 | None | False | None | 5052 | 954 | NonexNone | 179.0 |
| 8 | X2^4+Y2^4*delta^2 | 0000 | False | None | 524 | 130 | NonexNone | 0.8 |
| 8 | pol_r1+*pol_r1- | None | False | None | 5052 | 954 | NonexNone | 177.1 |
| 8 | pol_s0+*pol_s0- | 0000 | True | False | 718 | 109 | 3720x827 | 5.3 |
| 10 | D8* * S0 | 0000 | True | False | 1666 | 268 | 7639x1934 | 34.1 |

## Best residual landscape (full-E2 rel, log10 histogram of local minima over all restarts)

* d=8: 1e-3: 452, 1e-2: 1854, 1e-1: 863
* d=9: 1e-7: 120, 1e-4: 316, 1e-3: 284, 1e-2: 751, 1e-1: 1

Reading: a genuine witness in a family shows up as rel ~ 1e-14..1e-12 (float64 floor); the smallest values seen so far are listed above; all minima are separated from the floor by many orders of magnitude unless stated otherwise.

## Notable minima (checked by hand)

* d=9, sector 0100 (family = J_9[0100] (58) + D8*.X2, dimP 53): every one of 120 restarts converges to rel = 6.34e-7 (median 6.37e-7), three orders of magnitude below every other sector. Re-evaluated locally with lstsq cutoffs 1e-9 .. 1e-13: identical value (the 54 structural zero singular values of the Koszul block are cleanly separated from the rest), and polishing from the best point with 600 more L-BFGS iterations changes nothing: it is a genuine flat local-minimum valley at 6.3e-7, not a numerically masked zero (a witness would sit at ~1e-14). The mirror sectors 0101 / 0010 / 0051 (D8*.Y2, D8*.X1, D8*.Y1) sit at 1.3e-4 .. 3.2e-4.

## Structure notes (from the exact E1 tables)

* (J:A)/J has exactly ONE minimal generator beyond J in degrees <= 9: D8* = 4G_S^2 - G_P^2 (= const * G_u G_v, the product of the two reflection-contour quartics); the four degree-9 directions are D8*.X1, D8*.Y1, D8*.X2, D8*.Y2 (D8*.X0, D8*.Y0 lie in J_9). At degree 10 only 10 of the 556 non-Jacobian directions come from D8*.S_2; the other 546 are new minimal generators of the colon ideal (computed mod p this session). From d = 10 on E1 is essentially the Gorenstein-dual constraint codim (J:A)_d = dim R_{d+6} (minus 1 at d = 10): 950 = 951 - 1, 666, 426 for d = 10, 11, 12.
* Gorenstein duality (socle degree 24) pairs ker(.A: R_d -> R_{d+6}) with coker(.A: R_{18-d} -> R_{24-d}); the single kernel direction D8* at d = 8 is dual to the single missing direction of A.R_10 in R_16, consistent with the 556 = 1506 - 950 kernel at d = 10.
* Sector bookkeeping: G_aff characters chi = (e0+e3, e2+e5 mod 2; e1-e4 mod 6; Y-parity); the target residue is G_aff-invariant and the identity for a chi-eigen D lives in sector chi_A + 2 chi, so restricting all unknowns to sectors is lossless for eigen-denominators (24 sectors per degree parity).

## Errors / notes

* Engine: no engine files were modified; lane-local library w2lib.py. Observations about the shared engine are in /tmp/hodge/notes/engine-bugs.md.
