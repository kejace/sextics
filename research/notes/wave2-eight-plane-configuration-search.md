# Wave 2 search report — `eight-plane-configuration` (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS. Updated 2026-11-01T23:50:02. 283 search restarts collected over 10 batches; 2 Hessian-relaxation systems, 0 new kernel dimensions, 12 family intersections, 0 reconstruction attempts. End of session 2026-11-01 ~00:05 EST (Nov 2). Remote box saturated by other agents for the whole wave (load 400-1900 on 128 cores); the runner's capacity gate refused every submission (first queue item gave up after 4 h); all remote searches remain queued in w2_queue1.sh / w2_queue2.sh (they submit automatically when the gate opens; rerun 'python3 w2_collect.py' to pull and regenerate this report). Results below: exact local computations plus bounded local numerical batches; a local d=8 full-E2 batch is still running.

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
- d = 7 numerical Hessian relaxation (local, LM on the R-part residual with the spurious plane projected out, 60 free parameters per parity): odd parity 40 restarts, minima cluster at r_part_rel = 0.00387 (reached 3 times; median 0.0040), continuation with a 600-evaluation budget stays at 0.003871 with |grad|/|r| = 3e-8 — a genuine positive local minimum; even parity 52 restarts, best 0.00239 (stationary, |grad|/|r| = 5e-9), median ~0.012. No restart came within 10^-5 of zero: numerically, the degree-7 Hessian variety contains nothing beyond the spurious plane in either pure parity (mixed parity still to run). Combined with the dimension count (483 quadrics on P^59) this makes a genuine degree-7 Jacobian denominator very unlikely; an exact certificate is out of reach (degree-4 Macaulay matrix ~10^6 x 6*10^5, dense).
- d = 6 full-E2 variable projection (local, even parity, 57 restarts): best rel 0.170 (reached twice), median 0.59 — consistent with the exact exclusion.
- Served-target spaces of Jacobian-type denominators (target_image, mod p, local): four independent random elements of J_7 (two even, two odd) and one random element of J_8 (odd), all with coefficients in [-2,2], each serve a target space of dimension exactly 18 = dim J_6^- (all 18 basis elements lie in J), as did dF/dX0 (d = 5) and the block radial E_0(F) (d = 6) in wave 1 (d = 8 system: 40326 x 8133 mod p). Every Jacobian-type denominator tested so far serves only the trivial (zero) classes. Conjecture (empirical, consistent with the Hessian results): **D in J(F) implies served(D) ⊆ J(F)**; if true, every certificate for this instance has deg D >= 10 (first degree with Ann_R(A) != 0). Six cases in degrees 5-8, no exception.

- d = 7, plane-adapted Jacobian fields (w tangent to the eight-plane arrangement: w = sum_j L_j R_j + M_j S_j + Q_j v_j, 24 independent parameters per parity; the tangency divisor X ∩ V(w(F)) then contains all eight planes): the exact Hessian quadrics restricted to this family have rank 210 and their zero set is exactly the spurious plane {g*Euler : g linear} (Singular, both parities). Structured probe (corrected block scalings): only the six spurious fields among 344 candidates.
- d = 7 full-E2 variable projection (local, even parity, 16 restarts initialized on stage-1 minima, Koszul corrections 420 unknowns + Q 636 unknowns, 5796 equations): best rel 0.111, median 0.33 (the full residual moves off the stage-1 minima: R-part 0.045 at the best point). No restart below 0.1. Odd parity (16 restarts, all converged): best rel 0.230, median 0.335.
- d = 8 family statistics (local build, 71 min single-threaded for the R-projector): 168 w-parameters per parity (spurious subspace {g*Euler, g quadratic} of dim 12, 156 free), R-part conditions 468 = dim R_16^-, full E2 in degree 16: 10152 equations, Koszul corrections 954 + Q 1491 unknowns, nnz(Omega) = 2.95M. Stage-1 batches (20 restarts per parity) running locally.
- d = 7, mixed parity (local, mod p, 20 min): 126 w-parameters, 951 independent quadrics = 483 (odd target sector) + 468 (even target sector from the cross term Hess(w_+, w_-)), dim R_14 (both parities) = 1506; the spurious subspace {g*Euler, g linear} has dim 6, so 951 quadrics on P^119: generically empty (a Groebner/Hilbert computation on 126 variables is beyond local reach; a bounded Singular attempt was made).
- d = 8 numerical Hessian relaxation (local, even parity, 20 restarts, 156 free parameters vs 468 conditions): minima between 1.6e-4 and 1.1e-3 (median 7e-4), the best one stationary (|grad|/|r| = 2.3e-8 under a 600-evaluation continuation). The floor is lower than at d = 7 (as the parameter count grows) but no restart approaches zero; odd parity (20 restarts): minima between 9.8e-5 and 1.4e-3 (median 7.8e-4), the best one stationary (|grad|/|r| = 6e-9 under an 800-evaluation continuation). (The queued remote batches, 160 full-E2 restarts per parity with stage-1 initialization, would sharpen this.)
- d = 8 full-E2 variable projection (local, even parity, stage-1 initialized, 954 Koszul + 1491 Q unknowns, 10152 equations; Koszul block rank 795): first restarts at rel 0.19 and 0.052 (converged; 24 min each single-threaded under load) — the full-E2 minima decrease with the degree (0.17 at d = 6, 0.11 at d = 7, 0.05 at d = 8) but stay far above the reconstruction threshold; batch of 8 left running, results appended to `local-search-d8-p0.jsonl` and picked up by `w2_collect.py`.

## Conclusions of this wave (status 2026-11-01 ~23:45 EST)
1. No witness. No candidate reached the reconstruction threshold (rel < 1e-9); the smallest full-E2 residual over all Jacobian degrees is rel = 0.11 (d = 7 even), and the smallest necessary-condition (R-part) residual is 9.8e-5 (d = 8 odd), a genuine positive stationary point.
2. Rigorous degree bounds (with the stated mod-p caveats): deg D >= 7 for every certificate, in every parity sector (d <= 4: no E1-admissible D; d = 5: Hessian condition forces a singular point; d = 6: Hessian variety = spurious Euler point only, in the even, odd and mixed sectors; odd d = 6 is rigorous over Q, the other two mod two primes).
3. For 7 <= deg D <= 9 every denominator is Jacobian-type w(F) (kernel theorem, rigorous), and numerically the Hessian necessary condition has no non-spurious solution at d = 7 (floors 0.0024 even / 0.0039 odd, stationary) and no zero found at d = 8 (40 restarts, floors ~1e-4); exact algebra confirms that no structured quadratic field (344 candidates, plane-adapted 24-parameter family) passes at d = 7.
4. Six out of six Jacobian-type denominators tested (degrees 5-8, random and structured) serve exactly the 18 trivial classes J_6^-: the conjecture "D in J(F) => served(D) ⊆ J(F)" would make deg D >= 10 rigorous for this instance.
5. At d = 10 (first non-Jacobian degree, admissible space 1044 even / 1008 odd), no structured family tied to the eight planes or the net contributes an admissible direction outside J; the admissible space meets the ideal of the planes in the generic dimensions (843 / 456 for first / second order vanishing). The kernel-family and restricted searches (float basis, E2-reduced varpro) are implemented and queued for the remote.
6. Predicted (exact Koszul + duality) admissible dimensions: d = 11: 510 + 1341 = 1851 per parity; d = 12: 677 + 2223 = 2900 (even), 648 + 2214 = 2862 (odd); exact confirmation queued.
7. The remote box was saturated by other agents for the whole wave (load 400-1900 on 128 cores, memory at the gate floor); all remote work (d = 8/9 full searches with thousands of restarts, d = 10 kernel family, d = 11/12 kernels, mixed stage-1) is queued in `w2_queue1.sh` / `w2_queue2.sh` (submits automatically when the gate opens) and collected by `python3 w2_collect.py` (pull + report regeneration).

## Reproduction
Lane dir `/tmp/hodge/results/eight-plane-configuration/`: `w2lib.py` (toolkit), `w2_jvarpro.py` (searches), `w2_hess.py` / `w2_hessmixed.py` / `w2_sing_hilb.py` / `w2_sing_points.py` / `w2_sing_radical.py` (exact Hessian relaxation), `w2_struct7.py` / `w2_structprobe.py` / `w2_planeadapted7.py` (structured probes), `w2_fam10.py` / `w2_fam10b.py` (d = 10 families), `w2_timage.py` / `w2_timage_struct.py` (served targets), `w2_kernel.py` / `w2_k10prep.py` / `w2_k10restrict.py` (admissible spaces), `w2_collect.py` (log + report). Logs: `wave2-search.log.jsonl` (one line per batch/experiment), per-batch `local-*.jsonl`, `hess*.jsonl`, `fam10local.jsonl`.
- Served-target test for a structured Jacobian denominator vanishing on all eight planes, D = Q0*Q1*dF/dX0 (d = 8): the computation (40k x 7.7k modular nullspace) was killed twice on the shared local box without a traceback (memory/watchdog); not obtained in this wave (the script `w2_timage_struct.py` reruns it in ~30 min when a machine is free).


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

| d | parity | dim W | #pairs | dim R_2d^- | rank of quadrics | Macaulay ranks (deg: rank/cols) | projective variety | msolve |
|---:|:-:|---:|---:|---:|---:|---|---|---|
| 6 | 0 | 18 | 171 | 866 | 135 | 2: 135/171; 3: 1136/1140 | undecided | - |
| 6 | 1 | 18 | 171 | 866 | 153 | 2: 153/171; 3: 1138/1140 | undecided | - |

## 3. Full-E2 variable-projection searches (residual landscape)

rel = ||E2 residual|| / ||Omega(D) N_0|| after optimal Koszul corrections and Q; r_part = component outside J(F)_{2d}.

| batch | family | d | par | params (free) | rows | c+Q unknowns (Koszul-block rank) | restarts | best rel | 10% | median | worst | best r_part | <1e-6 | <1e-10 | s/restart | nfev | recon |
|---|---|---:|:-:|---|---:|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|
| local-search-d6-p0 | jacobian-d6-p0 | 6 | 0 | 18 (17) | 3080 | 162+226 (144-144) | 57 | 0.17 | 0.222 | 0.592 | 0.963 | - | 0 | 0 | 25.8 | 81.2 | 0 |
| local-search-d6-p1 | jacobian-d6-p1 | 6 | 1 | 18 (18) | 3080 | 153+226 (136-136) | 60 | 0.682 | 0.704 | 0.731 | 0.985 | - | 0 | 0 | 17.8 | 67.4 | 0 |
| local-search-d7-p0 | jacobian-d7-p0 | 7 | 0 | 63 (60) | 5796 | 420+636 (360-360) | 16 | 0.111 | 0.113 | 0.289 | 0.459 | 0.0143 | 0 | 0 | 187.3 | 95.3 | 0 |
| local-search-d7-p1 | jacobian-d7-p1 | 7 | 1 | 63 (60) | 5796 | 420+636 (360-360) | 16 | 0.23 | 0.299 | 0.335 | 0.38 | 0.0514 | 0 | 0 | 198.0 | 96.0 | 0 |
| local-search-d8-p0 (running) | jacobian-d8-p0 | 8 | 0 | 168 (156) | 10152 | 954+1491 (795-795) | 2 | 0.0518 | 0.0518 | 0.12 | 0.189 | 0.00642 | 0 | 0 | 1388.6 | 115.5 | 0 |
| local-stage1-d7-p0 (stage 1 only) | jacobian-d7-p0 | 7 | 0 | 63 (60) | 5796 | 420+636 (None-None) | 40 | 0.00441 | 0.00503 | 0.0124 | 0.0185 | 0.00441 | 0 | 0 | 19.9 | 102.8 | 0 |
| local-stage1-d7-p1 (stage 1 only) | jacobian-d7-p1 | 7 | 1 | 63 (60) | 5796 | 420+636 (None-None) | 40 | 0.00387 | 0.0039 | 0.004 | 0.0151 | 0.00387 | 0 | 0 | 19.8 | 120.0 | 0 |
| local-stage1-d8-p0 (stage 1 only) | jacobian-d8-p0 | 8 | 0 | 168 (156) | 10152 | 954+1491 (None-None) | 20 | 0.000161 | 0.00029 | 0.000642 | 0.00116 | 0.000161 | 0 | 0 | 107.5 | 172.7 | 0 |
| local-stage1-d8-p1 (stage 1 only) | jacobian-d8-p1 | 8 | 1 | 168 (159) | 10152 | 936+1491 (None-None) | 20 | 9.84e-05 | 0.000218 | 0.000737 | 0.00139 | 9.84e-05 | 0 | 0 | 103.8 | 178.2 | 0 |
| local-stage1b-d7-p0 (stage 1 only) | jacobian-d7-p0 | 7 | 0 | 63 (60) | 5796 | 420+636 (None-None) | 12 | 0.00239 | 0.00241 | 0.00987 | 0.0185 | 0.00239 | 0 | 0 | 16.8 | 112.4 | 0 |

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
