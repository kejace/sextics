# Wave 2 — dense-rational — SEARCH lens

STATUS: STARTED 2026-11-02 00:15 — infrastructure being written; no witness. (Rewritten after every batch.)

Target: q=1, A = X0^3 Y0^3, F = Fermat + (1..17)/9973 on all 230 non-pure even sextic monomials. Only symmetry: tau (Y -> -Y);
the only "character sectors" are the tau-parities of D (Y-even / Y-odd / mixed).

Log: /tmp/hodge/results/dense-rational/wave2-search.log.jsonl (one JSON line per batch). Scripts: w2s_*.py in the lane dir.

## Plan
1. E1-admissible spaces per degree 4..12, both parities (numerical SVD + exact cross-checks); kernel structure.
2. Jacobian-type degrees 7, 8, 9: Hessian relaxation (thousands of restarts, real + complex, pure + mixed sectors), then
   full normalized-E2 variable projection over the gauge-fixed w-space (inner unknowns: Koszul corrections K, rows reduced mod F).
3. Degree 10 (first non-Jacobian degree): same E2 varpro over K_full(10) = ker[A] + J_10 (F-multiples removed) and the kernel-only part.
4. Structured candidates (Hessian covariants of the cubic cover, 3x3 Hessian minors, gradient-norm forms, cones) at d=10,12: E1 then E2.
5. Rational reconstruction + exact re-solve for every candidate with residual < 1e-10.

## Results
(none yet)
