* Wave-1 (attack1) facts used: E1 kernel is zero for d <= 8 (so every legal eigen denominator of degree <= 8 lies in J(F));
  d=9 kernels 2,3,3,2,10,1,3 (even) / 0,5,0,0,9,6,4 (odd); all tested fixed denominators inconsistent.
* The control test (accepted witness `hesse-triangle`, family {D*} + J_6, 18 parameters): the solver reproduces rel ~1e-30
  at the witness direction and converges to it from perturbations of relative size ~0.3-0.5, but random restarts hit the
  basin only ~1% of the time (0/26 in local tests; remote job `w2-control` quantifies this).  Consequence: for d >= 9 the
  random-restart joint search is a landscape survey, not an exhaustive method; for d <= 8 the decisive tool is the
  low-dimensional mod-J relaxation (necessary condition), which is settled reliably by multistart and rigorously by Singular.
* Degree 6: Singular dim = 0 in 13 of 14 sectors (rigorous exclusion of all sigma-eigen degree-6 denominators of those
  parities/weights); sector (6, even, 0) has a 2-dimensional mod-p cone (one component, an octic in three chart variables)
  but no real zeros numerically (300 restarts, best rel 8.2e-3).
* Degree 7: Singular dim = 0 for even k=0,3,5,6 and odd k=1..6 (rigorous); even k=1,2,4 have 2-dim mod-p cones with no real
  zeros (best rel 1.3e-3 .. 2.7e-3 over 200 restarts); odd k=0 (24 unknowns): global minimum rel 7.709e-6 (stable under LM
  polish from three different minima, full-E2 rel there 0.80) -- a near miss but nonzero; Singular timed out at 45 s locally.
* Exact E2 at degree 6 (`w2_exactE2.py`, F-reduced Koszul bilinear system in (a, p), saturated by the ideal of the
  denominator coordinates, Singular over GF(32003)): sector (6, even, 0) -- the only one not killed by the relaxation --
  has DIM_SAT = -1 (empty): no sigma-eigen Y-even weight-0 denominator of degree 6 admits a certificate mod 32003.
  Together with the 13 relaxation exclusions: **no sigma-eigen denominator of degree 6 (either parity, any weight)**.
  Caveat on rigor: the relaxation statements are homogeneous in the denominator coordinates alone, so mod-p emptiness
  implies emptiness over Q-bar; the saturated E2 statement involves the affine numerator coordinates and a single prime,
  so it is strong evidence (to be repeated at further primes / in characteristic 0), not yet a proof.
* Degree 8 relaxation probes (local, 60-120 restarts, LM-polished minima): (8,even,0) global min rel 2.198e-8 (stable,
  nonzero; full-E2 rel there 0.25), (8,even,1) 2.45e-5, (8,even,4) 1.58e-5.
* (6, even, 0) exact E2 saturation repeated at p = 1000003: again DIM_SAT = -1 (58 s).  Two independent primes.
* Degree 8 relaxation, all 14 sectors probed locally (40-120 restarts each, LM-polished): global minima are nonzero,
  rel between 6.4e-7 (8,odd,4) and 3.0e-5 (8,odd,0); no real zeros.  Remote job `w2-hess8` repeats this with 1500 restarts
  and Singular (dim of the exact quadric system).
* Correction to the wave-1/dossier description: fqlx's degree-9 denominator `Jhat` lies in J(F)_9 (membership verified
  mod 2147483647 and 2147483629), i.e. it is of Jacobian type, not a non-Jacobian E1-kernel element.
* Sector E1 kernels at d = 9: the nontrivial (mod J) kernel directions in sectors (9,even,1) and (9,even,2) have no
  representative of height below ~1e21 (LLL on the mod-m kernel lattice, m = product of 4 primes, finds exactly the 63
  short J-vectors and nothing else), while (9,even,0/3/5/6) have representatives of height 6e3..3e5.  For the large-height
  sectors the joint search uses an orthonormal floating-point kernel basis (`float_kernel_reps`, lifts to 4e-15).
* Degree-9 kernel-only families (the non-Jacobian E1 directions, 1-10 dimensional per sector; rational reps where the
  height is small, orthonormal float reps otherwise), full-E2 LM, local runs (3-6 restarts each, all converging to the
  same minima across seeds): best rel (9,e,3) 3.4e-2, (9,e,6) 0.17, (9,e,4) 0.75, (9,o,5) 0.75, (9,o,4) 0.87, (9,o,1) 0.89,
  (9,e,1) 1.26, (9,o,6) 1.26, (9,e,2) 1.42; (9,e,0) 1.03 and (9,e,5) is a single direction (fixed-D, inconsistent in wave 1).
  No zeros: no certificate denominator in the span of the non-Jacobian kernel directions alone.
* Structured families (part 4, local runs at degrees 6-7, continued on the remote): pure-X cones (vertex P^2_Y) and
  pure-Y cones (vertex P^2_X) have no E1-admissible member in any sector at d = 6, 7 (kernel zero mod p => excluded as
  supports); line-vertex cones in (x,y,z,t), (x,y,z,r), (x,y,z,s) have a single E1-admissible direction each at (6,even,0)
  (a J-element), exactly tested: inconsistent; point-vertex cones with x missing: 2 directions at (6,even,0), both
  inconsistent.
* Why random-restart joint search cannot be exhaustive for d >= 9: a certificate (if any) is an isolated point of an
  overdetermined bilinear system (about 10x more equations than unknowns after eliminating Q); on the 18-dimensional
  control family the LM basin covers about 1% of the sphere and for 70-250 dimensional sector families the fraction is
  astronomically small.  The searches at d >= 9 are therefore landscape surveys plus exhaustive searches of the
  low-dimensional kernel-only and structured subfamilies; the decisive statements are the relaxation/exact-E2 exclusions
  at d <= 8 and the E1-kernel dimension table.
* Cross-reference: the parallel exclusion packet in this lane (`wave2-excl/relaxation-verdicts.md`, prime 1000003 and
  exact GB) reaches the same sector verdicts for d <= 7 by a different route: (6,even) = weighted-Euler curve C, all points
  excluded by exact constant-dual quadrics; (7,even,1/2/4) = X_i.C lifted components (excluded via the NZD lemma and the
  sum principle); (7,odd,0) lifted span(Y).C (pending); all other d <= 7 sectors EMPTY.  My relaxation/exact-E2 results
  (mod 32003 and 1000003) are consistent with every entry of that table.
* Structured families, local coverage: degrees 6 and 7 complete for all families (cones with plane/line/point vertex,
  fixed X-degree pieces, X-degree <= 3 / Y-degree <= 3 pieces, multiples of K, G, M, K+G and of the block Euler derivatives
  E_X F, E_Y F) in all 14 sectors; degree 8 complete for plane/line-vertex cones, partial for point-vertex cones.  Every
  nonzero E1-admissible part found so far lies inside J(F)_d (nontrivial mod J = 0) and has no E2 zero (best rel >= 3e-3,
  all exactly tested representatives inconsistent).  The remote job `w2-struct` extends this to degrees 8-12.
* Relaxation summary (best rel over all restarts; S = Singular dim of the exact mod-32003 quadric system, 0 = rigorous
  exclusion): d=6 even k=0 8.2e-3 (S=2, excluded by exact E2 saturation at two primes), k=1..6 0.13..0.86 (S=0); d=6 odd
  k=1..6 5e-3..0.15 (S=0); d=7 even k=0,3,5,6 (S=0), k=1,2,4 1.3e-3..2.7e-3 (S=2 = X_i.C lifts, cf. exclusion packet);
  d=7 odd k=0 7.7e-6 (24 unknowns, stable nonzero minimum), k=1..6 (S=0); d=8 even 2.1e-8 (k=0, 45 unknowns) .. 2.5e-5,
  d=8 odd 6.4e-7 (k=4) .. 3.0e-5 -- no real zeros anywhere at d <= 8.
* Second prime (1000003) for the exact relaxation dimensions at d = 6 and 7: identical verdicts to p = 32003 in all
  26 sectors (dim 0 everywhere except (6,even,0) and (7,even,1/2/4) with dim 2).
* Closest structured near miss: (8, odd, 4), cones not involving Y0 (6-dimensional E1 part inside J_8): 36 LM restarts give
  7 distinct local minima, best rel 1.46e-4 (single restart) / 2.36e-4 (30-restart batch), all exactly tested
  representatives inconsistent -- not a zero.
* Structured families, degree 8 (local): plane-, line- and point-vertex cone families now complete in all 14 sectors;
  every nonzero E1-admissible part lies in J_8 and has no E2 zero (smallest residuals 1.5e-4 .. 3e-3, see table).
* Structured families now complete locally for degrees 6, 7 and 8 (all families, all 14 sectors): no E2 zero anywhere;
  all nonzero E1-admissible parts are Jacobian-type.  Degrees 9-12 of the structured sweep remain for the remote job.
* Degree-9 structured sweep (local, partial): (9,even,0) and (9,even,1) plane/line-vertex cones -- E1 parts of dimension
  1-2 (Jacobian type), best rel >= 0.45; the remaining degree 9-12 cells are left to the remote job `w2-struct`.
* Control hit rate (hesse-triangle, 18-dim family {D*}+J_6, reduced system, varpro-LM): seeds 7/3 gave 0 hits in 26
  restarts; seed 101 found the witness (rel 8e-30) on the 4th restart.  Overall about 1 hit per 30 random restarts in
  18 dimensions; the reconstructed D agrees with rem_F(D*) and the exact re-solve confirms it (pipeline validated end to end).
* Control, seed 202: hit on the 3rd restart (rel 1.5e-24); the rem_F-chart reconstruction returned a rational point of the
  1-dimensional solution set different from D*, and the exact re-solve produced a verified certificate for it
  (confirm_exact -> solved, verify_witness True).  Running tally: 3 hits in 41 random restarts (~7%) on the 18-dim control.
* Degree-9 full sector families (kernel + J, or J alone where the kernel is zero), local LM samples of 2-4 restarts per
  sector (43-55 s per restart for K = 50-93): best rel (9,e,0) 0.16, (9,e,1) 0.054, (9,e,2) 0.042, (9,e,3) 0.018,
  (9,e,4) 0.56, (9,e,5) 0.012, (9,e,6) 0.050, (9,o,0) 6.3e-3 (J only, 93-dim), (9,o,1) 0.077, (9,o,2) 0.25, (9,o,3) 0.38,
  (9,o,4) 0.13, (9,o,5) 0.48, (9,o,6) 0.054 -- a landscape survey only (a few restarts per 50-90-dim family); the
  thousands-of-restarts runs are queued as remote jobs `w2-vp9`, `w2-vp9sparse`.
