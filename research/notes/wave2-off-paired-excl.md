# Wave 2 — lane `off-paired` — lens EXCLUSION / RIGOR

STATUS: no witness; last update 2026-11-01 21:33; 1803 logged experiments.

Packet: `/tmp/hodge/results/off-paired/wave2-excl/` (README.md there explains every receipt and how to replay). Log: `/tmp/hodge/results/off-paired/wave2-excl.log.jsonl`. Wave-2 primes 2147483629 / 2147483587 (fixed-D), 2147483647 (msolve).

## Rigor conventions
* modular rank <= rational rank; compressed rank <= modular rank.
* Fixed-D exclusion is PROVEN when rank_p([M|b]) = rank_p(M)+1 and rank_Q(M) <= rank_p(M) is certified by the explicit rational kernel (Euler gauge L*Z and F-multiples (F W, B_D(W)), dimension e computed exactly) or by exact verification of a rationally reconstructed basis of the whole modular kernel.
* E1 failure is PROVEN when rank_p(J_{d+6}[chi]) equals the exact Koszul dimension of J(F)_{d+6}[chi] and the augmented rank is one more.
* Bilinear sector results (msolve, charts) are modular SCREENS unless a rational certificate is recorded.

## Fixed-denominator receipts by family (counts of (label, prime) records)

| family | exp | results |
|---|---|---|
| block-product | e1 | PROVEN-E1-fail: 132, e1-pass-modp: 281 |
| block-product | fixedD | PROVEN-inconsistent: 554, inconsistent-modp(kernel-lift-incomplete): 8 |
| colon-quartic | e1 | e1-pass-modp: 2 |
| colon-quartic | fixedD | PROVEN-inconsistent: 4 |
| fermat-cone | e1 | PROVEN-E1-fail: 4 |
| partial*monomial | e1 | e1-pass-modp: 2 |
| partial*monomial | fixedD | PROVEN-inconsistent: 2, inconsistent-modp(nullity>explicit): 2 |

Distinct denominators tested: 350; PROVEN excluded (rational rank sandwich): 346; excluded modulo p only: 4; solvable mod p: []

## Bilinear sector systems (D free in a character sector; msolve on charts)

| d | char | mode | tsub | nt | ns | neq | charts run | no-solution | other | max s |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | [0, 0] | sector | all | 6 | 16 | 127 | 6 | 6 |  | 1 |
| 6 | [0, 1] | sector | all | 2 | 14 | 81 | 2 | 2 |  | 0 |
| 7 | [0, 0] | sector | all | 14 | 45 | 373 | 28 | 28 |  | 0 |
| 7 | [0, 1] | sector | all | 8 | 30 | 240 | 16 | 16 |  | 0 |
| 7 | [0, 2] | sector | all | 2 | 15 | 123 | 2 | 2 |  | 0 |
| 7 | [0, 3] | sector | all | 2 | 9 | 118 | 2 | 2 |  | 0 |
| 7 | [0, 4] | sector | all | 2 | 15 | 119 | 2 | 2 |  | 0 |
| 7 | [0, 5] | sector | all | 7 | 37 | 299 | 7 | 7 |  | 0 |
| 7 | [1, 0] | sector | all | 12 | 38 | 325 | 12 | 12 |  | 0 |
| 7 | [1, 1] | sector | all | 5 | 28 | 260 | 5 | 5 |  | 0 |
| 7 | [1, 2] | sector | all | 1 | 15 | 55 | 1 | 1 |  | 0 |
| 7 | [1, 3] | sector | all | 0 | 7 | None | 1 | 0 | None:empty-t-space | 0 |
| 7 | [1, 4] | sector | all | 1 | 15 | 78 | 1 | 1 |  | 0 |
| 7 | [1, 5] | sector | all | 6 | 33 | 241 | 6 | 6 |  | 0 |
| 7 | [2, 0] | sector | all | 11 | 38 | 327 | 11 | 11 |  | 0 |
| 7 | [2, 1] | sector | all | 5 | 33 | 259 | 5 | 5 |  | 0 |
| 7 | [2, 2] | sector | all | 1 | 15 | 53 | 1 | 1 |  | 0 |
| 7 | [2, 3] | sector | all | 0 | 7 | None | 1 | 0 | None:empty-t-space | 0 |
| 7 | [2, 4] | sector | all | 1 | 15 | 78 | 1 | 1 |  | 0 |
| 7 | [2, 5] | sector | all | 6 | 28 | 229 | 6 | 6 |  | 0 |
| 7 | [3, 0] | sector | all | 13 | 45 | 371 | 13 | 13 |  | 0 |
| 7 | [3, 1] | sector | all | 7 | 37 | 305 | 7 | 7 |  | 0 |
| 7 | [3, 2] | sector | all | 3 | 15 | 131 | 3 | 3 |  | 0 |
| 7 | [3, 3] | sector | all | 2 | 9 | 190 | 2 | 2 |  | 0 |
| 7 | [3, 4] | sector | all | 3 | 15 | 150 | 3 | 3 |  | 0 |
| 7 | [3, 5] | sector | all | 8 | 30 | 234 | 8 | 8 |  | 0 |
| 7 | [4, 0] | sector | all | 10 | 45 | 343 | 10 | 10 |  | 0 |
| 7 | [4, 1] | sector | all | 6 | 33 | 247 | 6 | 6 |  | 0 |
| 7 | [4, 2] | sector | all | 2 | 10 | 135 | 2 | 2 |  | 0 |
| 7 | [4, 3] | sector | all | 0 | 7 | None | 1 | 0 | None:empty-t-space | 0 |
| 7 | [4, 4] | sector | all | 2 | 10 | 144 | 2 | 2 |  | 0 |
| 7 | [4, 5] | sector | all | 5 | 28 | 249 | 5 | 5 |  | 0 |
| 7 | [5, 0] | sector | all | 11 | 45 | 338 | 11 | 11 |  | 0 |
| 7 | [5, 1] | sector | all | 6 | 28 | 238 | 6 | 6 |  | 0 |
| 7 | [5, 2] | sector | all | 2 | 10 | 102 | 2 | 2 |  | 0 |
| 7 | [5, 3] | sector | all | 0 | 7 | None | 1 | 0 | None:empty-t-space | 0 |
| 7 | [5, 4] | sector | all | 2 | 10 | 102 | 2 | 2 |  | 0 |
| 7 | [5, 5] | sector | all | 5 | 33 | 253 | 5 | 5 |  | 0 |
| 8 | [0, 0] | sector | all | 28 | 94 | 587 | 28 | 28 |  | 2 |
| 8 | [0, 1] | sector | all | 18 | 70 | 495 | 18 | 18 |  | 1 |
| 8 | [0, 2] | sector | all | 9 | 37 | 388 | 9 | 9 |  | 0 |
| 8 | [0, 3] | sector | all | 5 | 28 | 498 | 5 | 5 |  | 0 |
| 8 | [0, 4] | sector | all | 9 | 37 | 377 | 9 | 9 |  | 0 |
| 8 | [0, 5] | sector | all | 17 | 70 | 494 | 17 | 17 |  | 1 |
| 8 | [1, 0] | sector | all | 23 | 85 | 552 | 23 | 23 |  | 1 |
| 8 | [1, 1] | sector | all | 15 | 62 | 455 | 15 | 15 |  | 1 |
| 8 | [1, 2] | sector | all | 5 | 32 | 361 | 5 | 5 |  | 0 |
| 8 | [1, 3] | sector | all | 3 | 24 | 418 | 3 | 3 |  | 0 |
| 8 | [1, 4] | sector | all | 5 | 32 | 402 | 5 | 5 |  | 0 |
| 8 | [1, 5] | sector | all | 13 | 68 | 457 | 13 | 13 |  | 0 |
| 8 | [2, 0] | sector | all | 25 | 83 | 552 | 25 | 25 |  | 2 |
| 8 | [2, 1] | sector | all | 14 | 62 | 462 | 14 | 14 |  | 0 |
| 8 | [2, 2] | sector | all | 6 | 36 | 380 | 6 | 6 |  | 0 |
| 8 | [2, 3] | sector | all | 3 | 24 | 401 | 3 | 3 |  | 0 |
| 8 | [2, 4] | sector | all | 6 | 36 | 446 | 6 | 6 |  | 0 |
| 8 | [2, 5] | sector | all | 14 | 68 | 462 | 14 | 14 |  | 0 |
| 8 | [3, 0] | sector | all | 28 | 84 | 574 | 28 | 28 |  | 3 |
| 8 | [3, 1] | sector | all | 17 | 70 | 494 | 17 | 17 |  | 1 |
| 8 | [3, 2] | sector | all | 8 | 39 | 418 | 8 | 8 |  | 0 |
| 8 | [3, 3] | sector | all | 5 | 28 | 445 | 5 | 5 |  | 0 |
| 8 | [3, 4] | sector | all | 8 | 39 | 475 | 8 | 8 |  | 0 |
| 8 | [3, 5] | sector | all | 18 | 70 | 495 | 18 | 18 |  | 0 |
| 8 | [4, 0] | sector | all | 24 | 83 | 554 | 24 | 24 |  | 2 |
| 8 | [4, 1] | sector | all | 13 | 68 | 466 | 13 | 13 |  | 0 |
| 8 | [4, 2] | sector | all | 7 | 36 | 402 | 7 | 7 |  | 0 |
| 8 | [4, 3] | sector | all | 3 | 24 | 397 | 3 | 3 |  | 0 |
| 8 | [4, 4] | sector | all | 7 | 36 | 451 | 7 | 7 |  | 0 |
| 8 | [4, 5] | sector | all | 15 | 62 | 454 | 15 | 15 |  | 0 |
| 8 | [5, 0] | sector | all | 22 | 85 | 554 | 22 | 22 |  | 2 |
| 8 | [5, 1] | sector | all | 14 | 68 | 467 | 14 | 14 |  | 0 |
| 8 | [5, 2] | sector | all | 6 | 32 | 356 | 6 | 6 |  | 0 |
| 8 | [5, 3] | sector | all | 3 | 24 | 382 | 3 | 3 |  | 0 |
| 8 | [5, 4] | sector | all | 6 | 32 | 373 | 6 | 6 |  | 0 |
| 8 | [5, 5] | sector | all | 14 | 62 | 458 | 14 | 14 |  | 1 |
| 9 | [0, 0] | sector | all | 50 | 159 | 881 | 29 | 27 | 22:timeout 21:timeout | 60 |

## Colon / non-Jacobian E1-kernel dimensions (colon mod p = upper bound for Q; dim J exact by Koszul)

| d | colon | dim J_d | non-Jacobian kernel | gauge | prime | nonjac per sector (nonzero) |
|---|---|---|---|---|---|---|
| 4 | 1 | 0 | 1 | 0 | 2147483629 | 00:1 |
| 5 | 11 | 6 | 5 | 0 | 2147483629 | 01:1 10:1 20:1 30:1 35:1 |
| 6 | 54 | 36 | 18 | 1 | 2147483629 | 00:1 02:1 04:1 05:1 11:1 20:2 21:1 30:3 31:1 40:3 45:1 50:1 55:1 |
| 7 | 177 | 126 | 51 | 6 | 2147483629 | 00:5 01:1 03:1 05:3 10:3 12:1 14:1 15:3 20:2 21:2 22:1 24:1 25:1 30:4 31:3 32:1 33:1 34:1 35:1 40:4 41:3 50:5 51:1 55:2 |
| 8 | 457 | 336 | 121 | 21 | 2147483629 | 00:10 01:5 02:2 03:1 04:2 05:4 10:7 11:3 13:1 15:4 20:7 21:2 22:3 23:1 24:3 25:5 30:7 31:4 32:3 33:1 34:3 35:5 40:6 41:4 42:4 43:1 44:4 45:3 50:6 51:5 52:1 53:1 54:1 55:2 |
| 9 | 1014 | 756 | 258 | 56 | 2147483629 | 00:13 01:10 02:7 03:5 04:7 05:7 10:11 11:7 12:5 13:4 14:5 15:7 20:11 21:7 22:5 23:4 24:5 25:7 30:13 31:7 32:7 33:5 34:7 35:10 40:9 41:7 42:7 43:4 44:7 45:7 50:9 51:7 52:7 53:4 54:7 55:7 |
| 10 | 2173 | 1497 | 676 | 126 | 2147483629 | 00:19 01:18 02:22 03:21 04:22 05:17 10:17 11:17 12:18 13:21 14:18 15:15 20:19 21:16 22:21 23:21 24:21 25:16 30:21 31:17 32:19 33:21 34:19 35:18 40:18 41:15 42:22 43:21 44:22 45:17 50:16 51:16 52:19 53:21 54:19 55:16 |
| 11 | 3753 | 2682 | 1071 | 252 | 2147483629 | 00:26 01:25 02:34 03:37 04:34 05:27 10:23 11:25 12:33 13:36 14:33 15:25 20:24 21:26 22:33 23:36 24:33 25:24 30:27 31:27 32:33 33:37 34:33 35:25 40:26 41:25 42:33 43:36 44:33 45:25 50:25 51:24 52:33 53:36 54:33 55:26 |

## Notable single receipts

* {"label": "D0", "D": "D0", "d": 4, "char": [0, 0], "exp": "e1", "prime": 2147483629, "rows": 102, "cols": 72, "rank_aug": 69, "rank_J": 69, "dimJ_exact": 69, "result": "e1-pass-modp", "seconds": 0.0}
* {"label": "D0", "D": "D0", "d": 4, "char": [0, 0], "exp": "fixedD", "prime": 2147483629, "rows": 344, "cols": 75, "rank_M": 67, "rank_aug": 68, "explicit_kernel": 8, "result": "PROVEN-inconsistent", "seconds": 0.0}
* {"label": "D0", "D": "D0", "d": 4, "char": [0, 0], "exp": "fixedD", "prime": 2147483587, "rows": 344, "cols": 75, "rank_M": 67, "rank_aug": 68, "explicit_kernel": 8, "result": "PROVEN-inconsistent", "seconds": 0.1}
* {"label": "f6", "D": "f6", "d": 6, "char": [0, 0], "exp": "e1", "prime": 2147483629, "rows": 216, "cols": 192, "rank_aug": 177, "rank_J": 176, "dimJ_exact": 176, "result": "PROVEN-E1-fail", "seconds": 0.2}
* {"label": "f6", "D": "f6", "d": 6, "char": [0, 0], "exp": "e1", "prime": 2147483587, "rows": 216, "cols": 192, "rank_aug": 177, "rank_J": 176, "dimJ_exact": 176, "result": "PROVEN-E1-fail", "seconds": 0.1}
* {"label": "D0", "D": "D0", "d": 4, "char": [0, 0], "exp": "e1", "prime": 2147483629, "rows": 102, "cols": 72, "rank_aug": 69, "rank_J": 69, "dimJ_exact": 69, "result": "e1-pass-modp", "seconds": 0.0}
* {"label": "D0", "D": "D0", "d": 4, "char": [0, 0], "exp": "fixedD", "prime": 2147483629, "rows": 344, "cols": 75, "rank_M": 67, "rank_aug": 68, "explicit_kernel": 8, "result": "PROVEN-inconsistent", "rigor": "explicit-kernel", "seconds": 0.0}
* {"label": "D0", "D": "D0", "d": 4, "char": [0, 0], "exp": "fixedD", "prime": 2147483587, "rows": 344, "cols": 75, "rank_M": 67, "rank_aug": 68, "explicit_kernel": 8, "result": "PROVEN-inconsistent", "rigor": "explicit-kernel", "seconds": 0.0}
* {"label": "f6", "D": "f6", "d": 6, "char": [0, 0], "exp": "e1", "prime": 2147483629, "rows": 216, "cols": 192, "rank_aug": 177, "rank_J": 176, "dimJ_exact": 176, "result": "PROVEN-E1-fail", "seconds": 0.1}
* {"label": "f6", "D": "f6", "d": 6, "char": [0, 0], "exp": "e1", "prime": 2147483587, "rows": 216, "cols": 192, "rank_aug": 177, "rank_J": 176, "dimJ_exact": 176, "result": "PROVEN-E1-fail", "seconds": 0.0}

## Narrative (exclusion / rigor lens, wave 2)

### Best rigorous bound and what is new
* Unrestricted: deg D >= 7 for every legal denominator (fqlx, exact; re-verified here at degrees 4-6: D0 rank sandwich PROVEN, degree-6
  sectors chart-complete mod p).
* Single diagonal character (D a (g,h)-eigenvector): **degree 7 excluded for all 36 characters** — deanodotdotdot's theorem re-derived
  (i) chart-complete mod 2147483647 (171 charts + 4 empty colon sectors) and (ii) **exactly over Q by the relaxation pipeline for 34 of 36
  characters** ((1,0), (3,0): Singular runs pending on the remote); **degree 8 excluded for all 36 characters modulo p** (436 charts,
  chart-complete, including the six h-invariant characters (a,0) that the thread left open) and **exactly over Q for 21 characters** so far
  (9 attempted sectors have a positive-dimensional relaxation residual whose exact decomposition exceeded the local cap; the six (a,0)
  sectors and these 9 are queued in the remote chain `excl-off-paired-chain` with 2-hour Singular caps).  So, at the modular level, every
  single-character denominator has degree >= 9; at the theorem level the new content is the exact closure of 21 + 34 sectors by an
  independent method, and the method itself (below).
* Degree 9 single character: 27 of 50 charts of (0,0) closed locally (the rest need more than 60 s per chart); the full degree-9 sweep
  is in the queued remote chain.
* Structured fixed denominators (two new primes 2147483629 / 2147483587, explicit-kernel rank sandwich): see the receipts table — every
  tested block product of degree 7-14 processed so far is PROVEN excluded (either E1 fails against the exact Koszul dimension, or the full
  system is inconsistent with the rational rank certified); 8 records (3 denominators of degree 12-13) await the multi-prime
  re-certification (`w2_recert.py`, queued).  The wave-1 structured families (621) and Hessian minors (213) are queued after the block
  products in the same remote job (`excl-off-paired-batch1`, running).

### Method notes (why the exact sector exclusions are theorems)
* Fixed-D rank sandwich: rank_p(aug) = rank_p(M)+1 and rank_Q(M) <= rank_p(M) certified by the explicit kernel {(L Z,0)} + {(F W, B_D(W))}
  (exact dimension) or by exact verification of a reconstructed kernel basis.
* E1: modular rank of the Jacobian generators = exact Koszul dimension dim J(F)_{d+6}[chi] (F smooth) and augmented rank +1.
* Sector exclusion over Q: rational colon basis D_a modulo the Euler gauge and rational Koszul kernel S_k (verified exactly); exact bilinear
  equations; the degree-2 relaxation J2 (constant-coefficient combinations whose Koszul terms cancel; exact integer kernel, LLL-reduced) is a
  homogeneous ideal in Q[t]; if its reduction mod p is m-primary (msolve GB has a pure power of every variable) it is m-primary over Q
  (Hilbert function over Q <= Hilbert function mod p), so t = 0 and D is in the gauge.  Otherwise Singular over Q decomposes V(J2) and each
  component is excluded on its linear span by exact degree-3/4 localized relaxations (m-primary over Q).  Fallback: certified span of V(J2)
  (exact nilpotency certificates of the defining linear forms) + localized relaxation.  Degree-7 (0,0) and (2,0) were closed this way.
* Chart-complete msolve results are modular screens; a GB = {1} on a chart does not transfer to Q (the px-1 phenomenon), which is why the
  exact pipeline works with HOMOGENEOUS ideals in t only.

### Pending / how to pick up
* Remote jobs: `excl-off-paired-batch1` (running at session end; 346/680 block products done, then the 621 wave-1 denominators and 213
  Hessian minors follow automatically).  The combined chain (deltat full solve, exact rigor d=7/8 all sectors with 2-hour Singular caps,
  degree-9 chart sweep, multi-prime re-certification, kernel dims d=12) could NOT be submitted: the remote stayed above the runner's load
  guard (load > 150, < 60 GB free) for the whole session.  To launch it when capacity frees:
  `cd /home/kejace/devel/sextics-work && python3 -m hodge.rjob submit --lane off-paired --name excl-off-paired-chain --cpus 56 --wait 7200 -- 'bash /tmp/hodge/results/off-paired/wave2-excl/w2_chain.sh'`
  (run `python3 wave2-excl/w2_sync.py` first so the local log is a superset of the remote one; the chain writes per-process logs).  After they run: `python3 w2_sync.py && python3 w2_report.py
  report_notes.md` regenerates this report; receipts land in `/tmp/hodge/remote/off-paired/wave2-excl/receipts_*.jsonl`.
* Open at the end of this session: degree 7 unrestricted (mixed characters) — not attacked (dense 171-variable relaxation; the
  Hessian/Jacobian relaxation alone is not m-primary at degree 7); degree 9+ single character (modular, in progress); Delta_t (degree 18,
  E1 passes, queued).

