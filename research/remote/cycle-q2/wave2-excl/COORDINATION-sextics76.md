# Coordination note (agent session `sextics-76`, machine brunos, lens EXCLUSION/RIGOR, lane cycle-q2)

Another agent (engine `xlib-v1`, files `xlib.py`, `a_*.py`, `b_*.py`, `commentary.md`, `farkas/`, `receipts/`) works in this
same packet directory with the same lens.  To avoid duplicated remote compute, this is what `sextics-76` runs (files prefixed
as listed; all log lines carry `"lens": "excl", "engine": 2` and `"job": "excl-cycle-q2-*"`):

* `w2lib.py` — engine-v2 sector solver wrapper + exact Farkas certificates (`farkas()`), certificates in `certs/`.
* `reverify.py` (remote job `excl-cycle-q2-reverify`, 32 workers): 190 candidates = dossier A.5 (F_i, F_i^2, B6, B6^2, G_rad,i, G_rad,i^2,
  all 36 Z_j F_i, 2 monomials, K10, all 21 F_i F_j) + wave-1 structural families (pure-power, coupling-sum, blocks, pair-Fermat, edges,
  single couplings, all 15 Fermat 4-cones, all 20 Fermat 3-cones (cubic, quartic) and their squares).  Results: `reverify.results.jsonl`.
* `hess5.py` (done, local): 21 classes [H_ij] independent in R_16 at p = 10^9+7 and 998244353 => deg D = 5 excluded. `hess5.receipt.json`.
* `hess6.py --deg 6` + `gb_check.py` (remote job `excl-cycle-q2-hess6`, msolve 32 threads): full 62-dim E1 kernel at d = 6, 426 quadrics
  (Hessian relaxation, all sectors coupled), Groebner basis => decides deg D = 6 for ALL D (not only semi-invariant).
* `ann.sing` (remote job `excl-cycle-q2-ann`): Singular E1 annihilator dims d = 0..9 at p = 2^31-1 and 10^9+7, vdim(J) smoothness.
* planned: `kernels.py` / `sectortest.py` — sector-wise non-Jacobian E1-kernel elements at d = 7..9 (basis elements + Farkas), per-sector
  Hessian relaxation varieties; structured irreducible families at d <= 12.

Report (brunos): `/tmp/hodge/notes/wave2-cycle-q2-excl.md`.  Please do not overwrite; append your own section or use a different file name.

**WARNING (2026-11-01 ~20:50): log-file race.** `rjob submit`/`up` rsyncs the *local* lane dir over the remote one (no `--update`), so a local copy of
`wave2-excl.log.jsonl` overwrites the remote copy and silently drops lines that remote jobs appended since the last pull. This bit my resumed
`sectors.py --deg 6` run (records recovered from `sectors6.out` and `certs/`). My remote jobs now log to `wave2-excl.remote.log.jsonl`
(remote-only; merged into the canonical brunos log by `merge.py`). If your jobs append to the remote `wave2-excl.log.jsonl`, consider the same
split, or pull before every submit.
