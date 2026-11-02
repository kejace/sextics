# Wave 2 search report — short-power-presentation (q=1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (wave 2 started 2026-11-02 00:10; in progress)

Log: /tmp/hodge/results/short-power-presentation/wave2-search.log.jsonl (one JSON line per batch).
Scripts: /tmp/hodge/results/short-power-presentation/w2_*.py (remote jobs run there via hodge.rjob).

## Plan
1. E1-admissible spaces K_d^± (mult-kernel mod J) for d = 4..12, both Y-parities (the diagonal group of F is {1, tau}: Y-parity is the only character).
2. d <= 8: exact classification — E1 kernel + quadratic (Hessian) relaxation of E2 over the full family K_d + J(F)_d via Groebner bases mod p, then exact linear algebra on every relaxation solution.
3. Reduced variable-projection (Koszul-correction unknowns only) numeric search over K_d + J(F)_d, thousands of restarts, d = 7, 8 (and 9 if feasible).
4. Structured candidates at d = 7..12 screened exactly by E1 (normal form against the kernel basis), full solve for survivors.
5. Rational reconstruction + exact re-solve for every numeric candidate with small residual.

## Results so far
(none yet)
