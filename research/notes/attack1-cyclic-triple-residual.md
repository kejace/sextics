# Attack 1: `cyclic-triple-residual` (q = 1, A = X0^3 Y0^3) — resumed run

**No witness found.** The only "witness-candidate" file from the crashed first run
(`/tmp/hodge/results/cyclic-triple-residual/witness-candidate-1793572767.json`) is the *illegal control*
D = 6F (D must be coprime to F); it was produced on purpose to validate the solver and is not a certificate.

Work directory: `/tmp/hodge/results/cyclic-triple-residual/` (scripts, `log.jsonl` with one JSON line per tested
denominator, kernel elements `kernel_d*.json`). All linear algebra is exact over Q (python-flint `fmpz_mat`/`fmpq_mat`)
unless marked "mod p" (p = 2147483647, second prime 2147483629). The engine's sign bug in the F·N(D) term was fixed
before any result quoted here (all entries carry `"engine": 2`).

Headline results (details in sections 2-4):

1. **E1 is a sharp low-degree filter.** With R = Q[Z]/J(F) the Jacobian ring, the annihilator of [A] in R is
   zero in degrees d <= 7 and exactly one-dimensional in degree 8 (rigorous over Q: modular kernel dimensions are
   upper bounds for the kernel over Q, and dim J_d is the same over Q and F_p). Hence **every denominator of degree <= 7
   lies in J(F)** (in particular there is no denominator of degree <= 4 at all, since J_d = 0 there), and every
   denominator of degree 8 is congruent mod J(F)_8 to a multiple of one explicit polynomial `D8` (25 terms, sector
   (0,4) of the order-18 group). The degree-8 kernel is **not** forced by dimensions (sector dimensions 79 -> 97),
   so `D8` is a sporadic annihilator tied to the coupling structure; at degree 9 the kernel mod J has dimension 120 of
   which 116 are dimension-forced.
2. `D8` itself is excluded exactly over Q (sector-reduced system 2357 x 1188, rank 1155, augmented 1156), and random
   points of the coset D8 + J(F)_8 are excluded mod p (see table). The degree-9 sporadic kernel elements are excluded
   as listed below.
3. **Fermat calibration identifies the mechanism this class uses at t = 0**: on the Fermat sextic the target
   X0^3 Y0^3 is served by the *join divisor* D0 = X0^6 + X1^6 + X2^6 (equivalently Y0^6+Y1^6+Y2^6 = F - D0) with the
   strict block-Euler certificate N = A·(X0, X1, X2, 0, 0, 0)/6, Q = 0; it is **not** served by the (X0,Y0)-cone
   hyperplanes X0^6 + Y0^6 (the classical Shioda cones), nor by the reflection hyperplane X0^k, nor by any other
   Fermat-type sum tested. Geometrically: X ∩ V(D0) is the join of the plane sextic curves C_X, C_Y, whose resolution is
   a P^1-bundle over C_X x C_Y with H^{2,0} = H^{1,0}(C_X) ⊗ H^{1,0}(C_Y) (dimension 100), and the class of
   X0^3Y0^3 has Fermat character (4,1,1 | 4,1,1) of "join type" (both triples sum to 0 mod 6).
   The natural deformation of D0 to the coupled F, D = E_X(F) = Σ X_i ∂F/∂X_i (the ramification divisor of the
   degree-6 map X -> P^2 x P^2), fails (reduced Koszul system rank 370 of 388 columns, augmented 371, exact), and so do
   all twelve 0/±1 weighted Euler derivatives Σ m_i Z_i F_i (previous run) — the coupling terms of X-degree 2 reduce the
   multiplicity of T along the fixed curves from 6 (cone over a sextic, irrational) to 2 (transversal A_k, rational).

