# Wave 2 — off-paired — SEARCH lens

STATUS: STARTED 2026-11-01 ~23:55 EST. No witness. Reading done (dossier, attack1, wave-1 logs); building the
sector-restricted variable-projection library (w2lib.py) and submitting remote jobs. Nothing searched yet.

Lane dir: /tmp/hodge/results/off-paired/ (scripts w2lib.py, w2_dims.py, w2_varpro.py, w2_confirm.py).
Batch log: /tmp/hodge/results/off-paired/wave2-search.log.jsonl (one JSON line per batch; remote jobs write
wave2-<job>.jsonl + wave2-search.remote.log.jsonl, merged here after each pull).

## Plan
1. E1-admissible denominator spaces (J:A)_{d,chi} for d = 4..12, all 36 diagonal characters (exact rational bases,
   gauge F*R_{d-6} removed), Jacobian part J_{d,chi} and the non-Jacobian quotient.
2. Jacobian-type (w-space) search: Hessian relaxation (N_i = A w_i, Q free) then full E2 varpro on J_{d,chi}.
3. Full colon-family varpro per (d, chi), d = 8..12 (d = 8 only the six open (a,0) sectors), thousands of restarts.
4. Structured families (sub-linear-families and fixed candidates) via varpro / sweep; pending wave-1 Delta_t solve.
5. Rational reconstruction + exact re-solve for every candidate with relres < 1e-10.
