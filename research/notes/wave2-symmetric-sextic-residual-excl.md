# Wave 2 — `symmetric-sextic-residual` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; started 2026-11-01 ~19:00). Log: `/tmp/hodge/results/symmetric-sextic-residual/wave2-excl.log.jsonl`. Packet: `/tmp/hodge/results/symmetric-sextic-residual/wave2-excl/`.

## 0. Re-verification findings so far (reading wave-1 outputs)

* `scan_served.out` / `log.jsonl` (family `J6-radial-pencil`, entry `1EX(F)+1EY(F)`): this denominator is `E_X(F)+E_Y(F) = 6F` (Euler), i.e. **divisible by F (illegal)**; it trivially "serves" the whole 17-dim sector, and from that line on every `cum_beyondJ=17 / A_in_span=True / "!!!!! A ENTERS THE SPAN"` flag in the wave-1 scan is contamination by this one entry. The per-denominator results (`beyondJ=0` for every legal D) are unaffected. The attack-1 report did not claim a span hit, so nothing published is wrong, but the raw output must not be read naively.
* `log.jsonl` rows with `"consistent": true` in families `sgn3, sgn4, sgn5, sgn6, sgnsgn*, J6-radial-pencil, quot-curve6X, …` (rows ≥ 366) come from `CharSystem(..., targets=True)`: that system has the target A' as an unknown and is always consistent (A' = −A). They are **not** certificates. Genuine fixed-D results for the same D are `beyondJ=0`.

## 1. Plan
(a) re-verify rigorous claims with a second prime / exact arithmetic; (b) rigorous degree bound via the Hessian relaxation (all D for d ≤ 6 incl. J-part; pure-parity d = 7; sectors d = 8, 9) with msolve/Singular on the remote; (c) sector families d = 9–12 screened with the twisted E1+E2 system at two primes; (d) packet; (e) draft post.
