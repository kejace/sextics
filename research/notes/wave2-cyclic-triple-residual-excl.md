# Wave 2 — cyclic-triple-residual — lens EXCLUSION / RIGOR

STATUS: in progress (started 2026-11-01 ~19:05 local); no witness. RIGOROUS BOUND SO FAR: deg D >= 7 (no certificate with deg D <= 6, any parity/character). Degree-7 denominators lie in J(F); D8 is the unique degree-8 annihilator (exact). Gauge-fixed relaxation runs for d = 7, 8 in progress on the remote.

Log: /tmp/hodge/results/cyclic-triple-residual/wave2-excl.log.jsonl
Packet: /tmp/hodge/results/cyclic-triple-residual/wave2-excl/
Target: F = Fermat sextic + (5/7)(six couplings), A = X0^3 Y0^3, q = 1.

## Plan
(a) re-verify wave-1 rigorous claims with other primes / methods on the remote;
(b) rigorous lower bound on deg D via E1 (Jacobian-type reduction) + Hessian relaxation (Groebner, msolve);
(c) exhaustive sector-restricted tests of E1-kernel elements and structured families, d = 9..12;
(d) packet; (e) draft post.

## Results so far

### (a) Re-verification of wave-1 claims (job excl-ctr-kerverify, remote; receipts in wave2-excl.log.jsonl)
* E1-kernel dimensions dim{D in S_d : A D in J}/J_d recomputed at p = 1000003 and p = 32003 (wave 1: p = 2147483647):
  d = 5,6,7: 0;  d = 8: 1 (sector (0,4));  d = 9: 120;  d = 10: 556;  d = 11: 1020 — identical at all three primes,
  sector-by-sector identical. dim J_d mod p equals the complete-intersection Hilbert function in every degree 5..11 at both
  primes (6, 36, 126, 336, 756, 1497, 2682), which is what makes the modular kernel dimension an upper bound for the kernel
  over Q (rank_p J = rank_Q J and rank_p [A|J] <= rank_Q [A|J]).  AGREES with wave 1.
* Exact over Q (fmpz/fmpq): D8 is not in J(F)_8 (rank 39 -> 40 with D8 adjoined); A*D8 is in J(F)_14 with an explicit integer
  solution P (sum_i P_i d_i(7F) = A*D8, verified exactly; saved wave2-excl/e1_d8_exact.json).  Hence the degree-8 kernel is
  exactly one-dimensional over Q (<= 1 modular, >= 1 exact).  AGREES with wave 1 and sharpens it (wave 1 had only the modular bound).
* Consequence (rigorous): no denominator of degree <= 4 exists; every denominator of degree 5, 6, 7 lies in J(F)
  (D = G(F) for a vector field G of degree d-5); every denominator of degree 8 is c*D8 + G(F).

* Dimension-forced lower bounds (exact, from the G-graded Hilbert function of R = Fermat monomial count, local run,
  wave2-excl/forced_dims.json): forced kernel dimensions 0, 116, 555, 1020, 1325 at d = 8, 9, 10, 11, 12 versus modular
  kernel dimensions 1, 120, 556, 1020 (d = 8..11).  Sporadic (non-forced) annihilators: 1 at d = 8 (D8), 4 at d = 9
  (sectors (0,0), (2,5), (4,1), (4,4)), 1 at d = 10 (sector (2,1)), none at d = 11.  AGREES with wave 1 (116 + 4 at d = 9).
* Structure of the sporadic degree-9 annihilators (mod 1000003, local): X0*D8 and Y2*D8 lie in J(F); X2*D8 is the sporadic
  element of sector (2,5); Y0*D8, Y1*D8 account for the sporadic direction of sector (0,0); X1*D8 is a non-forced-looking
  but dimension-forced element of sector (0,3).  The sporadic elements of sectors (4,1) and (4,4) are NOT multiples of D8
  (wrong characters): two genuinely new degree-9 annihilators (exact representatives in kernel_exact_d9.json, heights
  ~1e14 and ~1e125).

### (b) Hessian relaxation (necessary condition A * G^T Hess(F) G in J(F)_{2d}, see wave2-excl/w2lib.py docstring)
* d = 5: 6 unknowns, 50 quadrics (mod 1000003); reduced Groebner basis = all 21 quadratic monomials, i.e. the ideal is
  m-primary mod p. Since the ideal is homogeneous and the mod-p system is the reduction of the rational one (dim J_10 mod p =
  c.i. value), properness of Proj over Spec Z gives V(I_Q) = {0}: NO DENOMINATOR OF DEGREE 5 (rigorous).
* Gauge fixing (makes the m-primary test applicable): the relaxation polynomial changes by an element of J under
  G -> G + h E (E^T Hess F = 5 grad F, E^T Hess F E = 30 F), so V(I) is a union of translates of the illegal locus L = {G = hE}
  (and {c = 0} at d = 8).  On the slice {theta_(Y2-component, Y2 | m) = 0}, which is complementary to L, V(I) subset L
  <=> the restricted ideal is m-primary.  m-primary mod p transfers to Q-bar: a nonzero projective zero over Q-bar, scaled to
  a primitive p-integral vector, reduces to a nonzero zero mod a prime above p.  (Requires the mod-p quadrics to be reductions
  of the rational ones: dim J_{2d} mod p = c.i. value, checked.)
* d = 6 (ALL parities and characters, 36 unknowns, 675 quadrics mod 1000003): the full relaxation ideal has V(I_p) = L
  (per-form Rabinowitsch tests, 35/35 empty) and the gauge-fixed ideal (35 unknowns) is m-primary (pure powers of all 35
  variables in the leading ideal; msolve -g 1, 4 s).  => NO DENOMINATOR OF DEGREE 6 (rigorous).  Receipts:
  wave2-excl/relax_d6_p1000003{,_g}*.json, *_gauge_lead.txt.
* Control (W2_FERMAT=1, log entries with job 'excl-ctr-relax6_fermat'): on the Fermat sextic F0 = sum Z_i^6, where the
  join-type sums X0^6+X1^6+X2^6 etc. are legal degree-6 denominators serving A, the same d = 6 pipeline gives only 52
  quadrics and pure powers for 12/35 variables, i.e. NOT m-primary — the test does not trivially report exclusion.
* d = 7, 8: running (job excl-ctr-relax78; d = 7: 126 unknowns, 1441 quadrics, msolve GB in progress since 00:36 UTC).
* Diagnostic (wave2-excl/wspan.py, p = 1000003, random family members): the degree-0 Koszul relaxation (functionals killing
  F*S and every Koszul image D div(S.gradF) - (S.gradF)(D)) has 2543 functionals at d = 6 versus 1751 for the Hessian
  relaxation, but only <= 1615 versus 1506 at d = 7 (Koszul span rank 8726 vs pi(J) rank 8835 in the 10341-dim quotient).
  So for d >= 7 the Hessian relaxation is essentially the whole degree-0 relaxation; the remaining strength of E2 sits in
  the theta-dependent part (full bilinear system).  Consequence for d = 9: the Hessian relaxation has only dim R_18 = 426
  quadrics for 876 unknowns, so degree-9 exclusions need the full parametric system (sector by sector) — see (c).

### (c) Degree 9+: exact kernel elements and sector families
* CRT/rational reconstruction of the canonical (rref) kernel basis fails at d = 9 for most sectors even with 6 primes
  (heights too large; job excl-ctr-kernel9, killed).  Replaced by exact integer nullspace + LLL (w2lib.kernel_exact_lll):
  small-height exact representatives with exact E1 lifts; reproduces D8 at d = 8 (max |coeff| 740880, verified).
  Job excl-ctr-kexact9: all 120 degree-9 kernel elements exactly, each tested with the sector system (modular ranks +
  exact null vectors, exact fmpz rank fallback).
* Note: a parallel excl-lens session running directly on the remote box (files lib.py, b2_hessian.py, receipts/ in the same
  packet dir) has the Y-EVEN-restricted Hessian relaxation at d = 5, 6 (m-primary, p = 1073741789) and is running d = 7
  Y-even; the runs here are for ALL parities (a general D need not be Y-even: the parity reduction D -> (D tau D)_red doubles
  the degree), so the two sets of receipts are complementary.
