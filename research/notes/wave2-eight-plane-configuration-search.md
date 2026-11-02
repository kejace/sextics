# Wave 2 search report — `eight-plane-configuration` (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS. Updated 2026-11-01T20:07:51. 0 search restarts collected over 0 batches; 0 Hessian-relaxation systems, 0 new kernel dimensions, 0 family intersections, 0 reconstruction attempts. 

## Inherited facts (wave 1, engine v2, mod p)
- E1-kernel mod J(F) is 0 for every degree d <= 9 and both Y-parities: every E1-admissible denominator of degree <= 9 is Jacobian-type, D = w(F) = sum_i w_i dF/dz_i (rigorous over Q by the rank argument).
- d = 10: kernel mod J = 294 (even) / 261 (odd) = exactly the forced minimum (multiplication by [A] has maximal rank); bases on the remote (p = 2147483629).
- Fixed-D tests (all inconsistent): six partials (d = 5), block radial E_i(F), X-Euler, X0*dF/dX0 (d = 6), random J_5..J_8 elements (both parities); target_image(dF/dX0) = target_image(block radial) = 18 = dim J_6^-: these denominators serve only the Jacobian sector of targets.
- Wave-1 E1-only sweeps of ~60 geometric candidates of degree <= 6 (Q-products, K/H quartics, Fermat parts, cones, odd low-degree forms): all fail E1, as the kernel theorem predicts.

## Method of this wave (lane dir `w2lib.py`, `w2_*.py`)
Normalize N(F) = A D exactly (always possible by subtracting (G/6)*Euler; Q unchanged). Since F is smooth the syzygies of (F_0..F_5) are Koszul, so E1 is solved analytically: N = N_0(a) + K_c with K_c = sum_{i<j} c_ij (F_j d_i - F_i d_j), deg c = d-4. E2 becomes
    D div N - N(D) = F Q      (degree 2d, odd parity)
instead of the cleared identity in degree 2d+6 (e.g. 3080 x 384 at d = 6 instead of 16797 x 2602). Variable projection: f(a) = min_{c,Q} ||Omega(D(a))(N_0(a) + K c) - F Q||^2 on the sphere, Q eliminated by an orthogonal projector, c by SVD, Kaufman Jacobian + Levenberg-Marquardt, random restarts in parallel processes. The spurious subspace {a : D(a) in F*S_{d-6}} is projected out — lossless, because (D + F g, N + (A g/6) E, Q') is again a certificate (checked symbolically), so only D mod F matters.
Necessary condition ("Hessian relaxation"): reducing E2 modulo J(F) kills the Koszul and F Q terms, leaving [D div N_0 - N_0(D)] = 0 in R_{2d}; for Jacobian type this is A * Hess_F(w,w) in J(F)_{2d}, i.e. [Hess_F(w,w)] in Ann_R(A)_{2d-6}: a system of quadrics in w. Number of independent quadrics = rank([A]: R_{2d-6}^+ -> R_{2d}^-) = 218 (d=6, 18 vars), 483 (d=7, 63 vars), 468 (d=8, 168 vars), 208 (d=9, 378 vars) — overdetermined for d <= 8 (generically empty, decidable exactly), underdetermined at d = 9 (variety of dimension >= 169).

## New rigorous consequence: degree-5 denominators are impossible
For d = 5, D = d_a F (a a constant vector, any parity) and the Hessian condition is A * d_a^2 F in J_10, i.e. [d_a^2 F] in Ann_R(A)_4 = 0 (wave-1 kernel computation, rigorous over Q), so d_a^2 F = 0 identically, which makes [a] a singular point of X. Hence no certificate with deg D = 5 exists, extending the wave-1 result (six coordinate partials) to the whole family. Degrees <= 4 were already excluded (no E1-admissible D). Together with the wave-1 kernel theorem, every certificate has deg D >= 6, and for 6 <= deg D <= 9 it is Jacobian-type.

## Self-tests (local, d = 6)
Assembled residual equals the exact Fraction computation (max diff 2e-11 at scale 4e4); spurious subspace at d = 6 even is exactly the Euler direction; R-part projector (dim 866 = dim R_12^-) annihilates Koszul and F Q terms; Kaufman Jacobian matches a fixed-x finite difference to 5e-10; LM converges in 50-70 evaluations (local minima at rel ~0.3 for random starts).

## Exact results obtained locally (wave 2, mod p = 2147483629; Singular for Hilbert data)
- d = 6, odd parity: 18 w-variables, 153 independent quadrics; Hilbert function 1, 18, 18, 2, 0, ... — the ideal contains all quartics, the projective variety is EMPTY. Rigorous over Q (Macaulay rank can only grow in characteristic 0): **no odd degree-6 denominator exists**.
- d = 6, even parity: 135 independent quadrics; projective dimension 0, degree 38; primary decomposition: a single component, the Euler point w = Euler (D = 6F, spurious), with multiplicity 38. Hence **no genuine even degree-6 Jacobian denominator** (mod p; every characteristic-0 point would have to reduce to the Euler point).
- d = 6, mixed parity (36 variables): 405 independent quadrics; projective dimension 0, degree 923; the quadrics have no linear term in the odd directions at the Euler point (H(E, v) = 5 v(F) in J), and the radical of the (dehomogenized) ideal has vector-space dimension 1: the scheme is supported at the single point w = Euler. Hence **no mixed-parity degree-6 Jacobian denominator** either (mod p).
- **Degree-6 verdict: no certificate with deg D = 6** — odd parity rigorous over Q; even and mixed parity modulo the standard mod-p caveat (confirmed with a second prime for the even sector: rank 135, degree 38).
- d = 5: excluded rigorously (Hessian condition forces a singular point), see above.
- d = 7 (both parities, local, mod p): 63 w-variables, exactly 483 independent quadrics each (= rank of [A]: R_8^+ -> R_14^-, as predicted), dim R_14^- = 744. The spurious plane {g * Euler : g linear} lies in the variety, so emptiness of the complement cannot be certified by Macaulay ranks (degree-4 matrix ~ 10^6 x 6*10^5, dense); decided numerically instead (stage-1 searches with the spurious subspace projected out). Structured probe: 344 natural quadratic fields (linear form x {Euler, block scalings, block rotations, coordinate fields}, Q_j x constant fields, gradients of all cubic monomials and of z_a Q_j) evaluated on the exact quadrics: only the six spurious fields z_a * Euler are Hessian-admissible.


## 1. E1-admissible denominator spaces per degree (mod p = 2147483629)

| d | parity | dim S_d | dim J_d | dim R_d | kernel mod J (E1-admissible non-Jacobian) | admissible total | source |
|---:|:-:|---:|---:|---:|---:|---:|---|
| 4 | 0/1 | 66/60 | 0 | 66/60 | 0 | 0 | wave 1 |
| 5 | 0/1 | 126/126 | 3/3 | 123/123 | 0 | 3/3 | wave 1 |
| 6 | 0/1 | 236/226 | 18/18 | 218/208 | 0 | 18/18 | wave 1 |
| 7 | 0/1 | 396/396 | 63/63 | 333/333 | 0 | 63/63 | wave 1 |
| 8 | 0/1 | 651/636 | 168/168 | 483/468 | 0 | 168/168 | wave 1 |
| 9 | 0/1 | 1001/1001 | 378/378 | 623/623 | 0 | 378/378 | wave 1 |
| 10 | 0/1 | 1512/1491 | 750/747 | 762/744 | 294/261 | 1044/1008 | wave 1 (remote) |

## 2. Hessian relaxation (necessary condition A*Hess_F(w,w) in J(F)_{2d}), exact mod p

(pending)

## 3. Full-E2 variable-projection searches (residual landscape)

rel = ||E2 residual|| / ||Omega(D) N_0|| after optimal Koszul corrections and Q; r_part = component outside J(F)_{2d}.

| batch | family | d | par | params (free) | rows | restarts | best rel | 10% | median | worst | best r_part | <1e-6 | <1e-10 | s/restart | nfev | recon |
|---|---|---:|:-:|---|---:|---:|---|---|---|---|---|---:|---:|---:|---:|---:|

## 4. Structured families at d = 10 intersected with the admissible space (mod p)

(pending)

## 5. Float admissible bases (d = 10)


## 6. Rational reconstruction attempts

None yet (no restart reached the reconstruction threshold rel < 1e-9).
