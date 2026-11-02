# Wave 2 — one-reflection-skew — SEARCH lens

STATUS: STARTED 2026-11-02 (VM). No witness. Reading done (dossier, wave-1 logs: 318 fixed-D 5-var tests d<=11, all inconsistent;
E1-kernel tables). Building the sector-restricted joint-search library (w2skew.py). Nothing searched yet.

Lane dir: /tmp/hodge/results/one-reflection-skew/ (wave-2 scripts w2skew.py, w2_*.py; jobs log to /tmp/hodge/jobs/search-one-reflection-skew-*.log).
Batch log: /tmp/hodge/results/one-reflection-skew/wave2-search.log.jsonl (one JSON line per batch).

## Plan
1. E1-admissible denominator spaces per degree d = 4..12: (a) five-variable reduction (D free of Y2; sectors (rho,sigma) in {+,-}^2),
   (b) six-variable per diagonal character (e5 mod 6, rho, sigma), with Jacobian part and non-Jacobian quotient dimensions.
2. Jacobian-type degrees: Hessian relaxation (A w^T Hess(F) w in J) then full-E2 varpro over the w-space (6-var, sector-restricted).
3. First non-Jacobian degrees (d = 7, 8, 9, ...): varpro over kernel + J_d (5-var, sector-restricted), thousands of restarts.
4. Structured families (cones over 3/4-variable subsets, bihomogeneous, invariant-ring products, binary-cubic covariants
   h_hom/k_hom, companion h/delta transported) via 5-var fixed-D tests and sub-family varpro.
5. Rational reconstruction + exact re-solve (5-var exact lift -> 6-var strict witness with D^2) for every candidate with relres < 1e-10.

## Batches
| ts | batch | family | d | sector | restarts | best rel | min/median iters | secs | note |
|---|---|---|---|---|---|---|---|---|---|
| 11-02T00:21 | start | - | None | None | 0 | None |  | 0 | wave-2 search started on VM; library construction |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 4 | rho0sig0 | 0 | None |  | 0.1 | nD=22 kerD=0 dimJ_d=0 nonJ=0 matrix 244x178 sv_gap=None |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 4 | rho0sig1 | 0 | None |  | 0.1 | nD=16 kerD=0 dimJ_d=0 nonJ=0 matrix 242x172 sv_gap=None |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 4 | rho1sig0 | 0 | None |  | 0.1 | nD=16 kerD=0 dimJ_d=0 nonJ=0 matrix 248x178 sv_gap=None |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 4 | rho1sig1 | 0 | None |  | 0.1 | nD=16 kerD=0 dimJ_d=0 nonJ=0 matrix 242x172 sv_gap=None |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 5 | rho0sig0 | 0 | None |  | 1.4 | nD=36 kerD=2 dimJ_d=2 nonJ=0 matrix 335x296 sv_gap=(0.00041212169351212635, 3.5527136788005025e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 5 | rho0sig1 | 0 | None |  | 0.5 | nD=30 kerD=1 dimJ_d=1 nonJ=0 matrix 334x290 sv_gap=(0.0003096390396150395, 4.100581074357064e-16) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 5 | rho1sig0 | 0 | None |  | 0.1 | nD=30 kerD=1 dimJ_d=1 nonJ=0 matrix 352x300 sv_gap=(0.000385110643052802, 3.552713678800501e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 5 | rho1sig1 | 0 | None |  | 0.1 | nD=30 kerD=1 dimJ_d=1 nonJ=0 matrix 335x290 sv_gap=(0.00015029153628498586, 8.698585526996496e-16) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 6 | rho0sig0 | 0 | None |  | 0.3 | nD=60 kerD=7 dimJ_d=7 nonJ=0 matrix 448x470 sv_gap=(7.2359173463574204e-06, 1.1909405287546993e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 6 | rho0sig1 | 0 | None |  | 0.9 | nD=50 kerD=6 dimJ_d=6 nonJ=0 matrix 448x460 sv_gap=(0.00014892553783067315, 1.2026186550753806e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 6 | rho1sig0 | 0 | None |  | 4.7 | nD=50 kerD=6 dimJ_d=6 nonJ=0 matrix 473x470 sv_gap=(7.864439401534822e-05, 1.6771845985915685e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 6 | rho1sig1 | 0 | None |  | 0.3 | nD=50 kerD=6 dimJ_d=6 nonJ=0 matrix 448x460 sv_gap=(6.24367336191745e-05, 3.507815069467059e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 7 | rho0sig0 | 0 | None |  | 3.0 | nD=90 kerD=28 dimJ_d=21 nonJ=7 matrix 588x705 sv_gap=(7.053834259805499e-10, 3.5210606109193014e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 7 | rho0sig1 | 0 | None |  | 4.1 | nD=80 kerD=19 dimJ_d=18 nonJ=1 matrix 588x695 sv_gap=(6.718892297484834e-15, 1.4515108722020692e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 7 | rho1sig0 | 0 | None |  | 3.1 | nD=80 kerD=18 dimJ_d=18 nonJ=0 matrix 616x710 sv_gap=(6.803086746021771e-10, 2.8655951429004052e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 7 | rho1sig1 | 0 | None |  | 0.6 | nD=80 kerD=19 dimJ_d=18 nonJ=1 matrix 588x695 sv_gap=(5.5500041248361314e-15, 3.5072803366595137e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 8 | rho0sig0 | 0 | None |  | 2.1 | nD=135 kerD=91 dimJ_d=46 nonJ=45 matrix 756x1025 sv_gap=(7.797812810775675e-05, 7.000059955173007e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 8 | rho0sig1 | 0 | None |  | 7.3 | nD=120 kerD=76 dimJ_d=43 nonJ=33 matrix 756x1010 sv_gap=(7.340408164801743e-05, 7.002261854364943e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 8 | rho1sig0 | 0 | None |  | 4.9 | nD=120 kerD=67 dimJ_d=43 nonJ=24 matrix 792x1025 sv_gap=(8.005053411721793e-06, 7.004823164754639e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 8 | rho1sig1 | 0 | None |  | 1.9 | nD=120 kerD=76 dimJ_d=43 nonJ=33 matrix 756x1010 sv_gap=(0.00014304945713989234, 6.9654687646712834e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 9 | rho0sig0 | 0 | None |  | 14.7 | nD=190 kerD=161 dimJ_d=92 nonJ=69 matrix 960x1436 sv_gap=(0.0003864068195329069, 6.963684843870582e-15) |
| 11-02T00:25 | e1-kernel5 | E1 kernel (5-var) | 9 | rho0sig1 | 0 | None |  | 5.2 | nD=175 kerD=146 dimJ_d=86 nonJ=60 matrix 960x1421 sv_gap=(0.00016550511232228112, 6.967952098221713e-15) |
| 11-02T00:26 | e1-kernel5 | E1 kernel (5-var) | 9 | rho1sig0 | 0 | None |  | 5.0 | nD=175 kerD=141 dimJ_d=86 nonJ=55 matrix 996x1442 sv_gap=(0.00043428162614897784, 7.000607495074884e-15) |
| 11-02T00:26 | e1-kernel5 | E1 kernel (5-var) | 9 | rho1sig1 | 0 | None |  | 3.9 | nD=175 kerD=146 dimJ_d=86 nonJ=60 matrix 960x1421 sv_gap=(0.0003247213419745124, 6.949252157606642e-15) |
| 11-02T00:26 | e1-kernel5 | E1 kernel (5-var) | 10 | rho0sig0 | 0 | None |  | 6.5 | nD=266 kerD=250 dimJ_d=161 nonJ=89 matrix 1200x1967 sv_gap=(0.006765765151149101, 6.8307429371494335e-15) |
| 11-02T00:26 | e1-kernel5 | E1 kernel (5-var) | 10 | rho0sig1 | 0 | None |  | 18.0 | nD=245 kerD=229 dimJ_d=153 nonJ=76 matrix 1200x1946 sv_gap=(0.018695602008891362, 6.7911138305558276e-15) |
| 11-02T00:26 | e1-kernel5 | E1 kernel (5-var) | 10 | rho1sig0 | 0 | None |  | 13.3 | nD=245 kerD=223 dimJ_d=153 nonJ=70 matrix 1245x1967 sv_gap=(0.00043503413489084176, 6.964437748219748e-15) |
| 11-02T00:27 | e1-kernel5 | E1 kernel (5-var) | 10 | rho1sig1 | 0 | None |  | 18.6 | nD=245 kerD=229 dimJ_d=153 nonJ=76 matrix 1200x1946 sv_gap=(0.0188506902799893, 6.88568017368748e-15) |
