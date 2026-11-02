# Draft discussion comment (reply under fqlx's comment in #125) — quadric-contact-residual, exclusion/rigor update

**Status: no certificate; rigorous lower bound still deg D >= 7; new exclusions at degrees 7 and 8; replication of fqlx's
exact results at a second prime.** Everything below is a *necessary-condition* computation (Hessian relaxation) or an exact
fixed-denominator solve; nothing here claims nonexistence of a certificate.

### Replication (new prime p = 2147483629, independent code)
* dim (J:A)_d for d = 0..9 = 0,0,0,0,0,6,36,126,337,760 (agrees with fqlx); the non-Jacobian E1 directions below degree 10 are
  exactly D8* = 4G_S^2 - G_P^2 (d = 8, invariant sector) and X1 D8*, Y1 D8*, X2 D8*, Y2 D8* (d = 9).
* Degree 5 (21 Hessian columns, rank 21 in R_10) and degree 6 (35 gauge-fixed unknowns, 630 pair columns of rank 334 in R_12,
  28 square-column kills exactly as listed by fqlx, then Groebner dim 0): reproduced. Degree <= 6 is excluded over C.
* The degree-8 direction D8* and X1 D8*, ..., Y2 D8* (degree 9): full sector-restricted fixed-D systems are inconsistent mod p
  (wave-1 result, re-run at the second prime) — i.e. no numerators exist for these denominators.

### New exclusions (Hessian relaxation: for D = sum c_k D_k in the E1 kernel with E1 lifts N_k, every certificate forces
sum c_k c_l B_{D_k}(N_l) in J_{2d}; the resulting homogeneous quadric ideal in c is shown to be m-primary mod p — msolve -g 1
leading ideal, or Singular dim 0 — which transfers to C by properness of Proj over Spec Z since J mod p keeps its Hilbert function)
* **Degree 7.** The 120-unknown family (D = C(F), C in S_2^6, Euler gauge removed) has exactly 950 independent quadrics
  (= dim R_8 - 1, via C^T Hess(F) C in (J:A)_8 = J_8 + <D8*>). Beyond fqlx's 8 pure parity classes, every subfamily in which
  two of the three involution characters (sigma0, sigma2, tau_Y) are fixed and everything else (incl. all six h1 characters)
  mixes is excluded: 12 subfamilies of 24-32 unknowns. So a degree-7 denominator would have to mix both eigenvalues of at
  least two of sigma0, sigma2, tau_Y simultaneously. [Single-involution subfamilies (60 unknowns): RUNNING.]
  The fully mixed system is out of Groebner reach: the Froeberg count for 950 quadrics in 120 unknowns gives a Hilbert series
  1, 120, 6310, 181240, ~2.9e6, <0, so a Nullstellensatz certificate needs degree >= 5 (Macaulay matrices ~3e8 x 2e8).
* **Degree 8.** (i) The Gamma-invariant family (G_aff plus the exchanges s0, r1; 20 parameters = 3 gauge + 16 Jacobian + D8*)
  is excluded (17 unknowns, 52 quadrics of rank 35, m-primary leading ideal with max pure power 5). (ii) 20 of the 24 pure
  G_aff character sectors are excluded (full gauge-fixed E1-admissible family of the sector: 4-17 unknowns, 37-52 quadrics).
  [The four large sectors (0,0,0,0), (0,0,0,1), (1,1,0,0), (1,1,0,1), 26-29 unknowns: RUNNING.]

### Structured supports (sum principle)
served(D) := space of Y-odd degree-6 targets admitting a certificate with denominator D; served(D1 D2) = served(D1) + served(D2)
(Gysin images are additive over components; served(D^k) = served(D)). For all of: the six coordinate hyperplanes, the quadrics
S0, delta = X0^2 - Y0^2, S1, S2, S0 +- S2, S0 + S1, X0^2 + X2^2, X0^2 - X1^2, X1^2 - Y1^2, X2^2 - Y2^2, the quartics G_S, G_P,
Q1, Q2, X0^4 + Y0^4, S0^2 + 2P0^2, X1^4 + Y1^4, the six polars dF/dZ_i, [f0 = S0^3 + 2 S0 P0^2, rad_i, X1^6+Y1^6, X2^6+Y2^6: pending]
served(D) = J_6^odd (dimension 18, the trivially served targets). Hence no product of any of these blocks, of any degree, serves
any nontrivial target, in particular not X0^3 Y0^3.

### Caveats / tooling
* Modular ranks are lower bounds for rational ranks; all "excluded" statements above are of the type "a modular rank is
  maximal" or "the mod-p ideal is m-primary", which transfer to characteristic 0.
* msolve misparses a product written `c1*c1` (it is NOT read as a square); squares must be written `c1^2`. Earlier msolve runs
  in this lane that used `c*c` syntax must be discarded (Singular runs are unaffected).

Reproduction: /tmp/hodge/results/quadric-contact-residual/wave2-excl/ (README.md, scripts, JSON receipts with primes and ranks).
