# Wave 2 search report — `eight-plane-configuration` (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS. Updated 2026-11-01T20:21:35. 29 search restarts collected over 2 batches; 0 Hessian-relaxation systems, 0 new kernel dimensions, 12 family intersections, 0 reconstruction attempts. 

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

- d = 10 structured families vs the E1-admissible space (kernel + J_10, exact mod p, local): G_k(Q) (par 0): dim 21, ∩admissible 0, ∩J 0; C[Q,K,H,diag] (par 0): dim 144, ∩admissible 12, ∩J 12; (Q)^4*S_2 (par 0): dim 156, ∩admissible 0, ∩J 0; (Q)^3*S_4 (par 0): dim 486, ∩admissible 18, ∩J 18; I_Pi^5 cap I_Pibar^5 (par 0): dim 231, ∩admissible 0, ∩J 0; F_i*F_j (par 0): dim 12, ∩admissible 12, ∩J 12; X-only (cone over plane curve, vertex Y-plane) (par 0): dim 66, ∩admissible 0, ∩J 0; Y-only (cone, vertex X-plane) (par 0): dim 66, ∩admissible 0, ∩J 0; oddquad*G_k(Q) (par 1): dim 135, ∩admissible 0, ∩J 0; X0Y0*C[Q,K,H,diag] (par 1): dim 72, ∩admissible 3, ∩J 3; I_Pi^5 cap I_Pibar^5 (par 1): dim 210, ∩admissible 0, ∩J 0; F_i*F_j (par 1): dim 9, ∩admissible 9, ∩J 9. **No structured family contributes an admissible direction outside J(F)**; the J-intersections are the spurious products F*(quartic) and F_i*F_j (sum principle).

## Structural observations on the Hessian map (used in the searches)
- Phi(w) := [A * Hess_F(w,w)] in R_{2d}^- satisfies Phi(w + g*Euler) = Phi(w) for every g of degree d-6, because Hess_F(Euler, v) = 5 v(F) in J and Hess_F(Euler, Euler) = 30 F in J. So the Hessian variety V_d is a cylinder over the spurious subspace Z_d = {g*Euler} (dim 0/1/3/12/28 for d = 5/6/7/8/9 per parity) and the search only has to run on the complement, which is exactly the projection used by the variable-projection code (lossless, since D and D + F g serve the same classes).
- The same identity explains the fat Euler point at d = 6 (multiplicity 38 even, 923 mixed): the quadrics have no linear terms in the odd directions at Euler.
- Dimension counts (per parity): #independent quadrics = rank([A]: R_{2d-6}^+ -> R_{2d}^-) = 218, 483, 468, 208 for d = 6, 7, 8, 9 (computed exactly for d = 6, 7; the d = 8, 9 values follow from the kernel theorem and Gorenstein duality R_16^- = R_8^-, R_18^- = R_6^-), against complement dimensions 17, 60, 156, 350. Generic expectation: empty for d <= 8 (confirmed exactly at d = 6), a variety of dimension >= 141 at d = 9, so d = 9 is the first Jacobian degree where the full E2 (Koszul corrections + Q) is the real test; its linear part has ~4970 unknowns against ~16.8k equations in degree 18, i.e. generically no solution unless (F, A) is special.
- Mixed parity (w = w_+ + w_-) is a separate system; it was excluded at d = 6 (single point) and is searched numerically at d = 7, 8 (queued).
- d = 10, ideal of the eight planes (exact mod p, local): (Q)*S_8 has dim 1311/1512 (even) and meets the admissible space in 843 dimensions (573 inside J); (Q)^2*S_6 (vanishing to order 2 along all eight planes): dim 924, meets the admissible space in 456 (237 in J); odd parity: 807/570 and 420/237. These are exactly the generic counts (1044 - codim), so the admissible space is not specially aligned with the plane arrangement; the order-2 family (456 / 420 parameters) is kept as the geometrically motivated restricted search space for d = 10. (Q)^3*S_4 meets the admissible space only inside J (18 = the F*(Q)^2-type spurious elements).

## Predicted E1-admissible dimensions for d = 10..12 (parity-refined Koszul complex + maximal rank of [A], the latter verified mod p through d = 10)
| d | parity | dim S_d | dim J_d | dim R_d | dim R_{d+6} (target parity) | kernel mod J | admissible total |
|---:|:-:|---:|---:|---:|---:|---:|---:|
| 10 | 0 | 1512 | 750 | 762 | 468 | 294 | 1044 |
| 10 | 1 | 1491 | 747 | 744 | 483 | 261 | 1008 |
| 11 | 0 | 2184 | 1341 | 843 | 333 | 510 | 1851 |
| 11 | 1 | 2184 | 1341 | 843 | 333 | 510 | 1851 |
| 12 | 0 | 3108 | 2223 | 885 | 208 | 677 | 2900 |
| 12 | 1 | 3080 | 2214 | 866 | 218 | 648 | 2862 |
Exact confirmation for d = 11, 12 is queued on the remote (`w2_kernel.py`); the d = 10 line reproduces the computed 294 / 261.


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

| batch | family | d | par | params (free) | rows | c+Q unknowns (Koszul-block rank) | restarts | best rel | 10% | median | worst | best r_part | <1e-6 | <1e-10 | s/restart | nfev | recon |
|---|---|---:|:-:|---|---:|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|
| local-search-d6-p0 (running) | jacobian-d6-p0 | 6 | 0 | 18 (17) | 3080 | 162+226 (144-144) | 13 | 0.17 | 0.222 | 0.592 | 0.813 | - | 0 | 0 | 25.2 | 78.8 | 0 |
| local-stage1-d7-p0 (stage 1 only) (running) | jacobian-d7-p0 | 7 | 0 | 63 (60) | 5796 | 420+636 (None-None) | 16 | 0.00441 | 0.00503 | 0.0148 | 0.0185 | 0.00441 | 0 | 0 | 17.1 | 98.9 | 0 |

## 4. Structured families at d = 10 intersected with the admissible space (mod p)

| d | par | family | #gens | dim fam | fam ∩ admissible | fam ∩ J | new directions | lifted |
|---:|:-:|---|---:|---:|---:|---:|---:|---:|
| 10 | 0 | (Q)^3*S_4 | 660 | 486 | 18 | 18 | 0 | 153 |
| 10 | 0 | (Q)^4*S_2 | 180 | 156 | 0 | 0 | 0 | - |
| 10 | 0 | C[Q,K,H,diag] | 246 | 144 | 12 | 12 | 0 | 114 |
| 10 | 0 | F_i*F_j | 12 | 12 | 12 | 12 | 0 | 12 |
| 10 | 0 | G_k(Q) | 21 | 21 | 0 | 0 | 0 | - |
| 10 | 0 | I_Pi^5 cap I_Pibar^5 | 231 | 231 | 0 | 0 | 0 | - |
| 10 | 0 | X-only (cone over plane curve, vertex Y-plane) | 66 | 66 | 0 | 0 | 0 | - |
| 10 | 0 | Y-only (cone, vertex X-plane) | 66 | 66 | 0 | 0 | 0 | - |
| 10 | 1 | F_i*F_j | 9 | 9 | 9 | 9 | 0 | 9 |
| 10 | 1 | I_Pi^5 cap I_Pibar^5 | 210 | 210 | 0 | 0 | 0 | - |
| 10 | 1 | X0Y0*C[Q,K,H,diag] | 114 | 72 | 3 | 3 | 0 | 45 |
| 10 | 1 | oddquad*G_k(Q) | 135 | 135 | 0 | 0 | 0 | - |

## 5. Float admissible bases (d = 10)


## 6. Rational reconstruction attempts

None yet (no restart reached the reconstruction threshold rel < 1e-9).
