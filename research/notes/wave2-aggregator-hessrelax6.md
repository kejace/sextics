# Aggregator cross-check: mixed-parity degree-6 Hessian relaxation (receipts)

STATUS: 2026-11-02 03:45 UTC — eight-plane-configuration mixed-parity degree 6 CLOSED (m-primary); diagonal-permutation reproduced (control); quadratic-cubic-pullback running; short-power-presentation, rank-one-biquartic, dense-rational (control) systems being built.

Why: the verify lenses (eight-plane row 20, qcp C28) point out that the degree-6 Hessian-relaxation exclusions of several lanes were run per Y-parity sector only, so a *mixed-parity* linear vector field `w = w₊ + w₋` (35 parameters after the Euler gauge) was not covered, and "`deg D ≥ 7`" held only for pure-parity denominators. This note closes that gap independently of the lanes.

Statement checked, per lane: if the E1 kernel mod J vanishes in degree 6 for both parities (so every E1-admissible degree-6 D is `D = w(F)` with w a linear vector field, 36 coefficients `c_ik`), then a certificate forces `A·wᵀHess(F)w ∈ J_12`, hence (injectivity of `[A]` on `R_6^±`) `wᵀHess(F)w ∈ J_6` — both parity parts separately. The Euler field `E` is the illegal gauge (`D = 6F`; the relaxation polynomial changes by `10h·w(F) + 30h²F ∈ J` under `w ↦ w + hE`), removed by the slice `c_55 = 0`. The ideal of the 426 quadrics (coordinates of `NF_J(wᵀHw)` in the standard-monomial basis of `R_6`) on the 35 slice variables has only the origin iff its leading ideal (msolve `-g 1`, grevlex) contains a pure power of every variable; mod p this transfers to ℂ when `dim J_6 mod p` equals the complete-intersection value (218/208 ring dimensions — checked in every run).

| lane | E1 kernel d = 6 (even, odd) | Hilbert check | #quadrics | leading ideal | verdict |
|---|---|---|---|---|---|
| diagonal-permutation (control; the excl lens reported m-primary for its own mixed 35/404 system) | 0, 0 | ✓ | 424 | 421 leading monomials, pure powers of all 35 variables (max exponent 3) | **m-primary — agrees with the lens** |
| eight-plane-configuration | 0, 0 | ✓ | 426 | 405 leading monomials, pure powers of all 35 variables (max exponent 2) | **m-primary — mixed-parity degree 6 excluded; `deg D ≥ 7` for arbitrary D** |
| quadratic-cubic-pullback | 0, 0 | ✓ | 426 | 405 leading monomials, pure powers of all 35 variables (max exponent 2) | **m-primary — mixed-parity degree 6 excluded; `deg D ≥ 7` for arbitrary D** (verify C28 covered pure parities only) |
| short-power-presentation | 0, 0 | ✓ | 426 (only ~132 independent: the Hessian of a sum of 8 sixth powers has rank ≤ 8) | running (msolve; the lens reports > 100 s and used P^k slices) | |
| rank-one-biquartic | 0, 0 | ✓ | 419 | 707 leading monomials, pure powers of all 35 variables (max exponent 3) | **m-primary — mixed-parity degree 6 excluded; with the excl lens's exact degree-5 exclusion (any parity): `deg D ≥ 7` for arbitrary D** (new: no lens had run the degree-6 relaxation) |
| dense-rational (control; lens: m-primary at p = 1073741789) | 0, 0 | ✓ | 426 | 405 leading monomials, pure powers of all 35 variables (max exponent 2) | **m-primary — agrees with the lens** |

Files: `/tmp/hodge/results/_aggregator/hessrelax6.py` (builder + output checker), `/tmp/hodge/results/_aggregator/hess6/hessrelax6-<lane>.{ms,json,g1.out}` (msolve inputs written with `^2`, never `c*c`; JSON summary with the kernel/Hilbert checks; leading-ideal output), remote job logs `/tmp/hodge/jobs/agg-hess6-*.log`. Prime 1000000007; msolve `-t 4 -g 1` (1 s per system on the remote).

Replay: `python3 hessrelax6.py <lane>` (builds the .ms locally, ~1–3 min), `msolve -g 1 -f hessrelax6-<lane>.ms -o out`, `python3 hessrelax6.py --check out 35`.
