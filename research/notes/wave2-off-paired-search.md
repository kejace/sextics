# Wave 2 — off-paired — SEARCH lens

STATUS: IN PROGRESS (2026-11-01 19:35 local). No witness. Step 1 (E1 dimension table d=4..12, all 36 characters) DONE.
Sector VarPro library validated on the Fermat positive control (exact f6 recovered to 7e-14, confirmed exactly).
Remote jobs queued (remote at load ~250/128, waiting for capacity): Fermat success-rate calibration, d=8 (0,0) colon search.
NOTE: a second agent (lens `excl`) works on this lane in wave2-excl/ and runs the Delta_t full solve and the two-prime
receipts for wave-1/dossier block products and 2x2 Hessian minors; those are NOT repeated here. My files live in
/tmp/hodge/results/off-paired/wave2-search/ (job names search-off-paired-s2-*).

Lane dir: /tmp/hodge/results/off-paired/wave2-search/ (w2lib.py = sector VarPro + exact bases; w2_dims.py; w2_varpro.py = LM
random-restart driver; w2_batch_sectors.py; w2_hess.py = Hessian relaxation; w2_homotopy.py; w2_structured.py; w2_control.py).
Batch log: /tmp/hodge/results/off-paired/wave2-search.log.jsonl (one JSON line per batch; remote jobs write
wave2-<job>.jsonl + wave2-search.remote.log.jsonl, merged here after each pull).

## Plan
1. E1-admissible denominator spaces (J:A)_{d,chi} for d = 4..12, all 36 diagonal characters (exact rational bases,
   gauge F*R_{d-6} removed), Jacobian part J_{d,chi} and the non-Jacobian quotient.
2. Jacobian-type (w-space) search: Hessian relaxation (N_i = A w_i, Q free) then full E2 varpro on J_{d,chi}.
3. Full colon-family varpro per (d, chi), d = 8..12 (d = 8 only the six open (a,0) sectors), thousands of restarts.
4. Structured families (sub-linear-families and fixed candidates) via varpro / sweep; pending wave-1 Delta_t solve.
5. Rational reconstruction + exact re-solve for every candidate with relres < 1e-10.

## Step 1 — E1-admissible denominator spaces (mod p = 2^31-1; verified against the dossier's exact colon dims 1,11,54,177,457,1014,2173,3753,5780)

| d | total colon | J part | non-J (ker mult-[A]) | gauge | eff (colon-gauge) | (0,0): colon/J/nonJ/gauge/eff | largest sectors (eff/nonJ) |
|---|---|---|---|---|---|---|---|
| 4 | 1 | 0 | 1 | 0 | 1 | 1/0/1/0/1 | (0,0):1/1 |
| 5 | 11 | 6 | 5 | 0 | 11 | 1/1/0/0/1 | (3,0):2/1 (0,0):1/0 (0,1):1/1 (0,5):1/0 (1,0):1/1 (2,0):1/1 |
| 6 | 54 | 36 | 18 | 1 | 53 | 7/6/1/1/6 | (0,0):6/1 (3,0):5/3 (4,0):5/3 (2,0):4/2 (5,0):4/1 (0,5):3/1 |
| 7 | 177 | 126 | 51 | 6 | 171 | 15/10/5/1/14 | (0,0):14/5 (3,0):13/4 (1,0):12/3 (2,0):11/2 (5,0):11/5 (4,0):10/4 |
| 8 | 457 | 336 | 121 | 21 | 436 | 30/20/10/2/28 | (0,0):28/10 (3,0):28/7 (2,0):25/7 (4,0):24/6 (1,0):23/7 (5,0):22/6 |
| 9 | 1014 | 756 | 258 | 56 | 958 | 55/42/13/5/50 | (0,0):50/13 (3,0):50/13 (4,0):44/9 (5,0):44/9 (1,0):43/11 (2,0):43/11 |
| 10 | 2173 | 1497 | 676 | 126 | 2047 | 95/76/19/10/85 | (0,0):85/19 (3,0):82/21 (2,0):77/19 (1,0):76/17 (4,0):76/18 (5,0):75/16 |
| 11 | 3753 | 2682 | 1071 | 252 | 3501 | 142/116/26/15/127 | (3,0):128/27 (0,0):127/26 (2,0):120/24 (1,0):119/23 (4,0):119/26 (5,0):118/25 |
| 12 | 5780 | 4437 | 1343 | 462 | 5318 | 211/177/34/26/185 | (0,0):185/34 (3,0):183/31 (4,0):175/31 (2,0):174/30 (5,0):173/31 (1,0):172/30 |

All 36 characters have eff > 0 from d = 8 on (at d = 8 only the six h-invariant sectors (a,0) are still open after deanodotdotdot's exclusions).
The non-Jacobian part (kernel of multiplication by [A] on R/J) is nonzero in EVERY degree from d = 4 (D0) on, so this lane has no
'Jacobian-type' degree: the colon family (J:A)_{d,chi} = J_{d,chi} + (non-J lifts) is the right search space everywhere.
Full table: /tmp/hodge/results/off-paired/wave2-search/dims.md (and dims.json).

## Optimizer calibration (local, Fermat F0 = f6+g6, d = 6, sector (0,0), K = 19, gauge 1, dimP = 18)

* objective at the true solution f6 (gauge-fixed): relres 3e-16; analytic gradient = finite differences.
* L-BFGS-B from 8 random starts: 0 successes (local minima at relres 0.0155 (x4), 0.125, 0.129, 0.185).
* Levenberg-Marquardt (Kaufman VarPro Jacobian) from 12 random starts: 1 success (seed 6: relres 7e-14 in 22 iterations, exact f6 reconstructed and re-solved),
  11 stalls at the same local minima. LM is the production optimizer; the landscape needs many restarts per family.
