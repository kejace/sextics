# Wave 2 search report: `rank-one-biquartic` (q=1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (wave 2 search started 2026-11-02 ~01:00 UTC on hodge-sextics; report rewritten after every batch).

Log: `/tmp/hodge/results/rank-one-biquartic/wave2-search.log.jsonl` (one JSON line per batch).
Code: `/tmp/hodge/results/rank-one-biquartic/w2/` (w2lib.py = reduced-E2 variable projection; jobs write `w2/out-<job>.jsonl`).

## Plan
1. E1-admissible denominator spaces (J:A)_d per parity, d = 4..12 (exact mod p dims; float bases for the search).
2. Jacobian-type degrees (d = 5, 6; w-space D = dF(G), N0 = A G): Hessian relaxation (no Koszul), then full E2 varpro.
3. First non-Jacobian degrees (d = 7, 8, 9): varpro over kernel + J_d (full E1-admissible space, F-multiples removed), thousands of restarts.
4. Structured families (E1-filtered) via sweep / fixed-D solver.
5. Rational reconstruction + exact re-solve for every candidate with relative residual < 1e-10.

(nothing searched yet)
