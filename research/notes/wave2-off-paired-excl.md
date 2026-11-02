# Wave 2 — lane `off-paired` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness (nothing consistent found so far).

Packet: `/tmp/hodge/results/off-paired/wave2-excl/` (scripts `w2lib.py`, `w2_cands.py`, `w2_batch.py`, `w2_deltat.py`; receipts `receipts_<job>.jsonl`).
Log: `/tmp/hodge/results/off-paired/wave2-excl.log.jsonl`. Wave-2 primes: 2147483629, 2147483587 (wave 1 used 1000003 / 2147483647; msolve 1073741827).

## Rigor conventions (used in every receipt)
* modular rank <= rational rank; compressed rank <= modular rank.
* Fixed-D exclusion is PROVEN when rank_p([M|b]) = rank_p(M) + 1 and rank_Q(M) <= rank_p(M) is certified either by the
  explicit rational kernel (Euler gauge L·Z, F-multiples (FW, B_D(W)); dimension e computed exactly) or by exact
  verification of a rationally reconstructed basis of the whole modular kernel ("kernel-certified").
* E1 failure is PROVEN when rank_p(J_{d+6}[chi]) equals the exact Koszul dimension of J(F)_{d+6}[chi] (F smooth) and the
  augmented rank is one more.

## Results so far
| experiment | result |
|---|---|
| D0 = F_X0/(6X0) (degree 4) full system, primes 2147483629/2147483587 | rank 67 / aug 68, nullity 8 = explicit kernel 8: **PROVEN inconsistent** (reproduces fqlx 67/68 and wave 1) |
| f6 E1 (degree 6) | rank_J 176 = exact dim J_12[(3,0)], aug 177: **PROVEN E1 failure** (two primes) |
| F_X0·X1·Y1 (degree 7) | rank 287 / aug 288, nullity 47 (explicit 43 + 4 extra, all 47 kernel vectors reconstructed with height <= 30 and verified exactly): **PROVEN inconsistent** |

## Running
* remote job `excl-off-paired-batch1`: rigorous two-prime receipts for the 680 dossier block products (degrees 7-14), the 621 wave-1 structured denominators and 213 Hessian 2x2 minors (E1).
* remote job `excl-off-paired-deltat`: full solve of the cubic-cover discriminant Delta_t (degree 18, G-invariant), pending since wave 1.
