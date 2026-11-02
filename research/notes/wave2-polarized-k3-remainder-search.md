# Wave 2 — `polarized-k3-remainder` joint search (remote 128-core box, 2026-11-01)

**STATUS: NO WITNESS FOUND.** (Report regenerated 2026-11-01 19:41 from 12 log records; search jobs listed at the end.)


Lane: q = 1, A = X0X1X2Y0Y1Y2, F = ½(f(U)+f(V)), U = X+Y, V = X−Y. Diagonal symmetry group = {1, τ} only, so the
relevant "character sectors" are the two Y-parities of D (Y-even = swap-symmetric in (U,V)). Wave-1 facts used:
[A] has maximal rank on (R/J)_d in every degree (injective d ≤ 9 ⇒ every E1-admissible D of degree ≤ 9 is
Jacobian-type D = Σ a_i F_i; surjective d ≥ 10), degrees ≤ 7 excluded (thread #115), sum principle ⇒ only D with an
irreducible factor of degree ≥ 8 can serve the target.

**Method (new in wave 2, `/tmp/hodge/results/polarized-k3-remainder/w2lib.py`).** Reduced Koszul E2 system of book
§5.2: for an E1-admissible D with exact lift N0 (dF(N0) = A·D; for Jacobian-type D = Σ a_i F_i take N0_i = A a_i),
E2 ⟺ B_D(N0) + T_D(c) ∈ F·S_{2d−6} with c antisymmetric of degree d−4 (the slack G in dF(N) = AD + FG drops out of
E2, so WLOG G = 0). The Q-block is eliminated exactly by the X0-remainder map (F is monic in X0), the c-block is a
dense least-squares problem, and the denominator coordinates a are optimised by Levenberg–Marquardt on the unit
sphere with the Kaufman (variable-projection) Jacobian; F-multiples are projected out of the a-space (lossless by the
prescribed-support theorem). Validated exactly at d = 6 (identity E = F·B_D(N) and the assembled system, mod p).
*Hessian relaxation* = the same system modulo J(F): A·Hess_F(a,a) ∈ J(F)_{2d} (necessary), solved numerically with
deflation against the known spurious family D = g·h (g = f(U)−f(V) ∈ J(F)_6, h a quadric), and exactly (Singular,
mod p) as a system of quadrics with saturation by the spurious component.


## 1. E1-admissible denominator spaces (engine `mult_kernel`, p = 2^31−1), both parities

| d | parity | dim R_d | dim J(F)_d | dim (R/J)_d | admissible dim {D: AD ∈ J} | kernel of [A] mod J | seconds |
|---|---|---|---|---|---|---|---|
| 4 | even | 66 | 0 | 66 | 0 | 0 | 0.7 |
| 4 | odd | 60 | 0 | 60 | 0 | 0 | 0.6 |
| 5 | even | 126 | 3 | 123 | 3 | 0 | 2.2 |
| 5 | odd | 126 | 3 | 123 | 3 | 0 | 2.2 |
| 6 | even | 236 | 18 | 218 | 18 | 0 | 7.2 |
| 6 | odd | 226 | 18 | 208 | 18 | 0 | 7.2 |
| 7 | even | 396 | 63 | 333 | 63 | 0 | 21.9 |
| 7 | odd | 396 | 63 | 333 | 63 | 0 | 22.1 |
| 8 | even | 651 | 168 | 483 | 168 | 0 | 94.5 |
| 8 | odd | 636 | 168 | 468 | 168 | 0 | 107.4 |
| 9 | even | 1001 | 378 | 623 | 378 | 0 | 248.2 |
| 9 | odd | 1001 | 378 | 623 | 378 | 0 | 264.2 |

For d ≤ 9 the admissible space is exactly J(F)_d (Jacobian type); the extra kernel appears at d = 10.

## Jobs / reproduction

Remote job logs are mirrored in `/tmp/hodge/remote/_jobs/`; scripts in `/tmp/hodge/results/polarized-k3-remainder/w2_*.py`; queued/launched commands in `w2_jobs.txt`. Each job appends to its own `w2-*.jsonl` (merged by `w2_merge.py`).
