# Wave 2 — `four-plus-two-residual` — EXCLUSION / RIGOR lens

STATUS: no witness (in progress). Jobs running on the remote: excl-kernels (E1-kernel re-verification, exact reps), excl-d5 (degree-5 Hessian relaxation, exact).

Log: `/tmp/hodge/results/four-plus-two-residual/wave2-excl.log.jsonl` (merged from remote per-job files `wave2-excl/remote/*.jsonl`).
Packet: `/tmp/hodge/results/four-plus-two-residual/wave2-excl/` (w2lib.py, w2_*.py, remote/ receipts, README.md).

## Rigor conventions
* [EXACT] computed over Q with fraction-free integer linear algebra (python-flint fmpz_mat); ranks are field-independent, so statements hold over C.
* [RANK] modular-rank argument for a fixed rational D: rank_p(M) <= rank_Q(M); the exact nullity of M is computed from the explicit parametrisation of the homogeneous solutions (N = m z/6 + K Koszul, D div K - K(D) = F Q); if it equals the modular nullity and the modular augmented rank exceeds it, the system is inconsistent over C.
* [MODP] a modular computation only (evidence, not proof).

## Results so far
(none yet — see log)
