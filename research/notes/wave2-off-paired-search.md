# Wave 2 — off-paired — SEARCH lens

STATUS: IN PROGRESS (2026-11-01 20:20 local). No witness; no nontrivially served class found so far.
Remote (128 cores) is at load 1000-1300 from ~50 agents (10x oversubscribed), so throughput is ~1/10 of nominal; jobs are
force-submitted through a local queue runner (op_s2_queue.py, <= 4 running, prefix search-off-paired-s2-).
DONE: step 1 (E1 dimension table); Fermat calibration d=6 (16/96 exact hits), d=7 (10/96); exact served-space scan of 3922
structured denominators (all 3/4/5-variable-subset forms of degree 4-6 in every G-character, Fermat-type cones, named blocks,
random colon elements at d=9..12): 0 nontrivially served (3,0)-classes; d=8 (0,0) colon search, first 52 restarts (exact-QR
solver): best relres 1.3173e-2 (a hyperplane-arrangement attractor, see landscape note), 25 distinct minima; d=8 (1,0): 25 restarts,
best 7.3e-2. Solver sped up 3x (ridge normal equations + SVD polish, identical minima); d=8 jobs relaunched (d8c00b, d8h0b).
RUNNING/QUEUED: d8c00b, d8h0b, d9h0, d10h0, hess, struct-hess3, homotopy(+served tracking, multi-start), d11, ft-d9h0, d8hmix,
d12, d9h0more, struct-cones, d8gmix, d10rest, ft-d8h0, struct-hess4, d7tau, d10h0more; servedscan-all (all 18 Y-odd target sectors).
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


## Methods (implemented in /tmp/hodge/results/off-paired/wave2-search/w2lib.py)

* **Exact bases.** (J:A)_{d,chi} and J_{d,chi} are computed as canonical rref bases over Q (modular rref with 2^31-sized primes,
  CRT, rational reconstruction, verification at a further prime); the Euler gauge F*R_{d-6} is removed by an orthogonal slice.
* **Sector VarPro.** For D(a) in the family, the fixed-D system (all six numerators of degree d+1 and Q of degree 2d-6, restricted
  to the characters forced by chi_D, which is lossless) is assembled matrix-free from a precomputed sparse tensor; the Q-block is
  projected out once (QR); f(a) = min_x ||M(a)x - b(a)||^2 with b = -A D(a)^2, on the unit sphere of the gauge complement
  (so D is never an F-multiple). Optimizer: Levenberg-Marquardt with the Kaufman variable-projection Jacobian (checked
  against finite differences); L-BFGS was 0/8 on the control, LM 1/12.
* **Free-target variant.** Parameters (D, A') with A' in the 24-dim target sector modulo its 2-dim Jacobian part (22 free
  directions): a zero means D serves SOME nontrivial G-invariant (3,1)-class. By the sum principle served(D1 D2) =
  served(D1) + served(D2), several such D could be multiplied to serve A. Exact confirmation of a hit: round D, compute
  served(D_r) exactly (sector target_image mod p), round-and-project if trivial. On the Fermat control the hits lie on a
  curved positive-dimensional variety and rounding leaves it (served = J only); a hit on off-paired would need structural
  analysis before exact confirmation.
* **Hessian relaxation (w2_hess.py).** D = sum w_i F_i, N_i = A w_i (exact E1), residual = distance of
  D*div(A w) - A*w.grad(D) to F*R_{2d-6}; LM over the w-sphere modulo Koszul/gauge directions. Sufficient only.
* **Homotopy (w2_homotopy.py).** F_t = Fermat + t*Qx tracked from the Fermat solution f6*X0^{d-6} in the full sector family.
* **Structured candidates (w2_structured.py).** 3x3 / 4x4 Hessian minors (degrees 12, 16; not in wave 1 or the excl lens
  which does 2x2), mixed Fermat-type cones times blocks (degrees 8-12). Full sector-reduced fixed-D solve after E1.

## Job queue (remote, prefix search-off-paired-s2-; local runner w2_queue.py keeps <= 4 running)
control (Fermat calibration d=6,7,8), d8c00 (800 restarts), d8h0 (sectors (a,0), a=1..5, 600 each), d9all (36 sectors x 300),
ft-d9 (free-target, 36 sectors x 200), hess (d=8..12, six h-invariant sectors, 600 each), struct-hess3, homotopy (d=8,9,10),
struct-cones, struct-hess4; then d10all / d11h0 / d12h0 sized from the calibration.

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

## Landscape note (d = 8, sector (0,0), colon family, K = 30, gauge 2, dimP = 28)

The best local minimum found so far (relres 1.3173e-2, reached by several independent restarts; reproduced locally from seed 0)
is essentially a binary octic in (X0, Y1), even in both variables:
D ~ -0.654 X0^4 Y1^4 + 0.489 X0^2 Y1^6 + 0.381 X0^6 Y1^2 - 0.259 Y1^8 - 0.169 X0^8 + small corrections
(sigma2-even, 74 % of its norm in J_8). Geometrically this is a union of eight hyperplanes X0 = c_i Y1 through the P^3
{X0 = Y1 = 0}; such unions cannot support (3,1)-classes (book ch. 3), so the dominant attractor of the landscape is a known
non-solution type that merely minimizes the defect (20 % of the residual sits on five degree-22 monomials
X0^2 X1 Y0^4 Y1^9 * {Y0^6, Y1^4, Y2^6, X2^6, X1^6}). Consequence: genuine structured solutions, if any, have small basins;
restart counts must be large, and the free-target variant (any served class) is the more sensitive probe.
