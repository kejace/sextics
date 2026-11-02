# symmetric-sextic-residual — wave-2 exclusion / rigor packet

Instance: `symmetric-sextic-residual` (q = 1, A = X0 X1 X2 Y0 Y1 Y2,
F = p6(X) + p6(Y) + (5/7) e2(X) e2(Y) (e2(X) + 2 e2(Y))). Variables X0,X1,X2,Y0,Y1,Y2.
Engine: `~/devel/sextics-work/hodge` (v2). Primes: 2147483647 (wave 1), 2147483629, 1000000007 (this packet).
All scripts take absolute paths; replay from this directory with `python3 <script> ...`.

## Contents
* `xlib.py` — helpers: vectorised Jacobian matrices mod p, complete-intersection (CI) dimension targets, logging.
* `a_verify.py` → `a_*.json` — adversarial re-verification of the rigorous claims (second primes, engine code path,
  exact rational membership A*D7 ∈ J(F) by CRT lift + exact re-multiplication, age table, automorphism group).
* `b_hessian.py d par [--kernel D7|D7Y|none] [--fermat]` → `b_hessian_*.json`, `b_quadrics_*.ms/.sing`, `b_groebner_*.json` —
  Jacobian-type reduction + Hessian relaxation (see docstring for the exact quadratic condition and the rigor argument).
* `c_sectors.py dmin dmax [types]` → log lines `c:*` — character-sector full E1+E2 solves of E1-kernel elements.
* Log of every experiment: `../wave2-excl.log.jsonl` (one JSON line each; fields include prime, ranks, seconds).

## Rigor conventions
* "rank mod p = CI value" for J(F)_k means the Jacobian columns have their characteristic-0 rank mod p; since
  rank_p ≤ rank_Q and the Q-rank of a smooth sextic's Jacobian ideal is the CI value, equality is checked, not assumed.
* dim of an E1 kernel mod p is an UPPER bound for the Q-dimension once rank J_k mod p = CI value in degrees k = d, d+6
  (rank argument on [J | A·S_d]); exact rational membership is the LOWER bound. Both are recorded.
* "inconsistent (mod p)" fixed-D results are modular evidence only, unless the packet gives a rational argument.
* Hessian-relaxation exclusions are rigorous over C: full rank of the quadric coefficient matrix mod p, or a mod-p
  Hilbert function vanishing in some degree, implies the same over Q (rank_p ≤ rank_Q on integer matrices) and
  hence over C.

(Statements and tables are filled in as receipts land; see the note /tmp/hodge/notes/wave2-symmetric-sextic-residual-excl.md.)
