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

## Results so far (updated 2026-11-01 19:44)

### Fixed facts (wave 1 + wave 2 so far)
* E1 kernel K_d^pm (mod J): d<=6: 0 (both parities, rigorous mod p); d=7: 1 (even), 1 (odd); d=8 even: 6. [wave 1, two primes]
* **D7** (the even d=7 kernel element, rational, 168 terms): irreducible; depends on Y only through M = Y0+3Y1+5Y2; M-degrees 0,2,4; min X-degree 3.
  So V(D7) is a cone with vertex the line {X0=X1=X2=M=0}. Full fixed-D solve (remote, engine v2, p=2^31-1): **inconsistent** (rank 4035, nullity 462 of 26532 x 4497);
  served space = J_6^- only (dim 18): D7 alone serves no primitive target.
* Structural fact (wave 2): 2F_{X0}-F_{X1}, 3F_{X0}-F_{X2}, 3F_{Y0}-F_{Y1}, 5F_{Y0}-F_{Y2} cancel the coupling, so mod J(F):
  X1^5 = 2 X0^5, X2^5 = 3 X0^5, Y1^5 = 3 Y0^5, Y2^5 = 5 Y0^5, and X0^5 = -H(L,M)/6, Y0^5 = -G(L,M)/6 with H = d_L C, G = d_M C,
  C = c((L+M)^6+(L-M)^6).  Consequently many "binary/ternary" forms lie in J_d (e.g. X0^7 - 2 X0^5 X1^2 in J_7, verified).
* Gauge invariance (proved in wave 2, used to shrink all searches): if (D, N, Q) is a certificate with N(F) = A D then
  (D + hF, N + hA E/6, Q + h div N - N(h) + h^2 A) is a certificate for every h in S_{d-6}.  So denominators are searched modulo F S_{d-6}.
* E2-mod-J relaxation (exact, correct form A*Phi in J_{2d}; gauge-fixed; msolve chart + Singular, p = 2^31-1):
  d=6 even (17 vars, 599 quadrics): no projective solutions (msolve [-1]; Singular dim 0, vdim 1216).
  d=6 odd (18 vars, 601 quadrics): no projective solutions (msolve [-1]; Singular dim 0, vdim 540).
  => no certificate with deg D = 6 of pure parity (reconfirms wave 1 with the corrected relaxation).
* Structured-family screen at d=7 even (152 subring families = monomials in 2..5 of {X0,X1,X2,Y0,Y1,Y2,L,M}; exact mod two primes):
  the only families meeting K_7 + J_7 outside J_7 are the ones spanning <X0,X1,X2,M> (they contain D7); all other families meet only J_7.


### d = 7 structure (wave 2, exact)
* K_7^+ = <D7> mod J_7, D7 in C[X0,X1,X2,M] = P7(X) + P5(X) M^2 + P3(X) M^4 (cone, vertex line {X=0, M=0}); rational, 4-prime numerators.
* K_7^- = <D7'> mod J_7, D7' in C[Y0,Y1,Y2,L] = Q7(Y) + Q5(Y) L^2 + Q3(Y) L^4 + Q1(Y) L^6 (cone, vertex line {Y=0, L=0}); rational (3 primes, 93-bit modulus),
  canonical cone representative in D7odd_cone.json (42 coefficients, denominators ~1e8).  Note the asymmetry (M-degree <= 4 vs L-degree <= 6).
* The exact subring-family screen shows no other cone/plane/point-vertex family meets K_7 + J_7 outside J_7 (both parities).
* Kernel numerators: even one reconstructed with 4 primes; odd NF representative's numerators have heights > 124 bits (not needed:
  the relaxation computes numerators mod p and the varpro uses float least-squares numerators, which is exact-equivalent since
  Koszul corrections span the particular-solution ambiguity).

### Tooling validated (wave 2)
* w2lib.JRef / e1_kernel: reproduces wave-1 ranks (J_10^- = 747, kernel 0 at d=4).
* w2_varpro (reduced variable projection: unknowns = Koszul corrections c_ij (deg d-4) + Q; rows = degree-2d monomials):
  quadratic term, Koszul block and Q block checked against direct float polynomial arithmetic (errors 1e-13..1e-17);
  gradient checked by finite differences (ratio 1.000000); Fermat control (D = X0^6+X1^6+X2^6 serves X0^3Y0^3) gives f = 1e-34.
* Landscape note: residual f/scale0 at random points is ~1; local minima at d=6 have f/scale0 ~ 0.95 (no basin of a zero; consistent with the exact exclusion).

### Remote status
* Job search-spp-j1 (done): odd d=7 kernel vector found mod 2 primes but rational reconstruction failed (height); D7 tests as above; e1vec 8 0 running/finished.
* Jobs search-spp-kA (kernels d=7,8 both parities, 3 primes, numerators) and search-spp-kBC (d=9, d=10) are queued with --wait: the remote has been
  refusing submissions (load 300-550 on 128 cores, 40-74 jobs from other agents) since 00:20.
* Local fallback (brunos, 1-2 threads): d=7 even relaxation over D7 + J_7 running; d=7 odd kernel (2 primes) running.
