# Wave 2 — cyclic-triple-residual — SEARCH lens

STATUS: STARTED 2026-11-02 00:20 (VM local). No witness. Reading done (dossier, attack1, wave-1 logs, wave-1 outputs
of kdims/d8limit/planes). Building the sector-reduced, E1-eliminated variable-projection library
(/tmp/hodge/results/cyclic-triple-residual/wave2/w2lib.py); no search batch run yet.

Lane dir: /tmp/hodge/results/cyclic-triple-residual/ (wave-1 scripts: sect.py sector solver, kernel8/9.py, planes.py ...)
Wave-2 dir: /tmp/hodge/results/cyclic-triple-residual/wave2/
Batch log: /tmp/hodge/results/cyclic-triple-residual/wave2-search.log.jsonl (one JSON line per batch)

## Plan
1. E1-admissible spaces (J(F)_d^chi + kernel reps) for d = 4..12, all 18 characters of the order-18 diagonal group
   (both Y-parities), dimensions logged.
2. Jacobian-type degrees (d = 5, 6, 7; kernel mod J = 0): Hessian relaxation (N = A w, Q free) then full E2
   varpro (Koszul-correction unknowns z) over the w-space, many restarts in parallel jobs.
3. Non-Jacobian degrees d = 8 (D8 + J_8), 9 (120 + J_9), 10, 11, 12: sector-restricted varpro over kernel + J_d.
4. Structured families (cones, discriminants, block Hessians, block radial derivatives, companion-transported
   denominators) at degrees 6-12 via sweep / sector solver.
5. Rational reconstruction + exact re-solve for every candidate with relative residual < 1e-10.

## Results
(none yet)
