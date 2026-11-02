# Wave 2 — lane `off-paired` — lens EXCLUSION / RIGOR

STATUS: no witness; last update 2026-11-01 19:23; 24 logged experiments.

Packet: `/tmp/hodge/results/off-paired/wave2-excl/` (README.md there explains every receipt and how to replay). Log: `/tmp/hodge/results/off-paired/wave2-excl.log.jsonl`. Wave-2 primes 2147483629 / 2147483587 (fixed-D), 2147483647 (msolve).

## Rigor conventions
* modular rank <= rational rank; compressed rank <= modular rank.
* Fixed-D exclusion is PROVEN when rank_p([M|b]) = rank_p(M)+1 and rank_Q(M) <= rank_p(M) is certified by the explicit rational kernel (Euler gauge L*Z and F-multiples (F W, B_D(W)), dimension e computed exactly) or by exact verification of a rationally reconstructed basis of the whole modular kernel.
* E1 failure is PROVEN when rank_p(J_{d+6}[chi]) equals the exact Koszul dimension of J(F)_{d+6}[chi] and the augmented rank is one more.
* Bilinear sector results (msolve, charts) are modular SCREENS unless a rational certificate is recorded.

## Fixed-denominator receipts by family (counts of (label, prime) records)

| family | exp | results |
|---|---|---|
| colon-quartic | e1 | e1-pass-modp: 2 |
| colon-quartic | fixedD | PROVEN-inconsistent: 4 |
| fermat-cone | e1 | PROVEN-E1-fail: 4 |
| partial*monomial | e1 | e1-pass-modp: 2 |
| partial*monomial | fixedD | PROVEN-inconsistent: 2, inconsistent-modp(nullity>explicit): 2 |

Distinct denominators tested: 3; PROVEN excluded (rational rank sandwich): 3; excluded modulo p only: 0; solvable mod p: []

## Bilinear sector systems (D free in a character sector; msolve on charts)

| d | char | mode | tsub | nt | ns | neq | charts run | no-solution | other | max s |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | [0, 0] | sector | all | 6 | 16 | 127 | 6 | 6 |  | 1 |
| 6 | [0, 1] | sector | all | 2 | 14 | 81 | 2 | 2 |  | 0 |

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
