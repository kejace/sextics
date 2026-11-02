# Wave 2 search report — `eight-plane-configuration` (q = 1, A = X0^3 Y0^3)

STATUS: STARTED (2026-11-01, late evening). No witness. Building the E2-reduced joint-search toolkit; no remote jobs submitted yet.

## Inherited facts (wave 1, engine v2, mod p)
- E1-kernel mod J(F) is 0 for every degree d <= 9 and both Y-parities: every E1-admissible denominator of degree <= 9 is Jacobian-type, D = w(F) = sum_i w_i dF/dz_i.
- d = 10: kernel mod J = 294 (even) / 261 (odd) = exactly the forced minimum (mult. by [A] has maximal rank); bases on the remote (p = 2147483629).
- Fixed-D tests (all inconsistent): six partials (d=5), block radial E_i(F), X-Euler, X0*dF/dX0 (d=6), random J_5..J_8 elements; target_image(dF/dX0) = target_image(block radial) = 18 = dim J_6^- (these serve only the Jacobian sector).

## Plan (this wave)
1. Dimensions of the E1-admissible spaces for d = 4..12 (kernel + J_d), both parities.
2. Jacobian degrees d = 5..9: Hessian relaxation (necessary condition A*Hess_F(w,w) in J, exact mod p / msolve) and full-E2 variable projection over the w-space with Koszul corrections (reduced system in degree 2d instead of the cleared identity in degree 2d+6), many restarts on the remote.
3. d = 10: varpro over kernel + J_10 (float kernel basis, particular E1 solutions), parity-restricted.
4. Structured families at d = 10..12 intersected with the admissible space, then sweep.
5. Rational reconstruction + exact re-solve of any candidate with residual < 1e-10.
