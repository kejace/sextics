# Wave 2 — quadric-contact-residual — lens EXCLUSION / RIGOR

STATUS: COMPLETE for this session (remote follow-ups running). No witness. Best rigorous bound: every certificate denominator has deg D >= 7 (fqlx; d <= 6 independently replayed here at
two new primes). NEW rigorous exclusions inside degrees 7, 8, 9 (items 3-4, 7-8, 12 below: at degree 7 every pair of the eight
(sigma0, sigma2, tau_Y)-parity classes is excluded, so a degree-7 denominator must spread over >= 3 parity classes; at degree 8 the
Gamma-invariant family and 20 of 24 pure G_aff sectors; at degree 9 the 4 closable sectors) and a structured-support theorem (item 9).
Full mixing of all G_aff sectors at degree 7 is provably out of reach of the mod-J (Hessian) relaxation by Groebner methods (item 5).
Remote jobs still running at the time of writing (remote load 1000-2000, i.e. 8-15x oversubscribed, so they are slow):
excl-sub7 (degree-7 half-cubes = single-involution subfamilies, 60 unknowns, msolve -g 1 then block-Macaulay k=4 fallback),
excl-d8large (the 4 large degree-8 sectors), excl-planes7 (the 8 diagonal affine planes of the parity cube, 60 unknowns each).
excl-lattice7 (bilinear tests) was stopped once the zero-pair hypothesis was refuted (item 18).  excl-misc finished its replays and served-space parts; its degree-10 extras screen was stopped to free the
slot (16 of the 556 directions were solved locally/remotely, all inconsistent).  Results append to
/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-{sub7,d8large,planes7,misc}.log.jsonl on the remote (merge with
wave2-excl/merge_remote_logs.py).  Queued but not submitted: parity8.sh (8 parity classes at degree 8, regularity 5).

Lane: q = 1, A = X0^3 Y0^3. Primes: P2 = 2147483629 (main), 1000000007 (confirmation); wave 1 used 2147483647, the dossier
1000000009, Singular 32003.  Code + receipts: /tmp/hodge/results/quadric-contact-residual/wave2-excl/ (README.md there).
Log: /tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl (one JSON line per experiment, prime recorded).

## 10-line summary (best rigorous bound first)
1. Best rigorous bound unchanged: deg D >= 7 for every certificate denominator (fqlx), now replayed independently at two primes
   (E1 colon dims d <= 9, degree-5/6 closures) with Singular and msolve verdicts agreeing.
2. NEW degree 7: all 28 pairs of (sigma0,sigma2,tau_Y)-parity classes excluded (24-32 unknowns each, two primes) — a degree-7
   denominator must spread over >= 3 parity classes.  (A difference-multiplicity lemma would have reduced larger supports to the
   14 affine planes, but its hypothesis fails: zero pairs of the cross forms exist — item 18.)  The six half-cubes and the eight
   diagonal planes (60-unknown problems) are running on the remote.
3. Full mixing at degree 7 is provably beyond the mod-J relaxation + Groebner (950 quadrics / 120 unknowns: regularity 5).
4. NEW degree 8: the Gamma-invariant family (incl. D8*) and 20 of 24 pure G_aff sectors excluded (two primes); 4 large sectors running.
5. NEW degree 9: the 4 closable sectors excluded; the other 20 are positive-dimensional under the relaxation (fewer quadrics than unknowns).
6. Structured supports: served(D) = J_6^odd (trivial) for 6 hyperplanes, 11 quadrics, 7 quartics, 6 polars, f0, rad_i, Fermat blocks
   and D8* itself (V(D8*) ∩ X supports no nontrivial class) => no product of these blocks, at any degree, serves A (sum principle).
7. Replays: D8*, Z_i D8* fixed-D systems inconsistent at P2 with identical ranks (sector solver and engine); fqlx's K3 Gauss-Manin
   socle test replayed exactly over Q; dossier's epsilon-sector dims replayed; dossier's "(J:A)_13 = S_13" corrected (8322 = 8568-246).
8. Tooling: msolve misparses `c*c` squares and FGLM fails for p >= 2^31; wave-1 msolve verdicts discarded; use `c^2` + `msolve -g 1`.
9. Screens: 7 of 20 trivial-sector and 9 of 16 (0,0,0,1)-sector non-Jacobian directions at d = 10 solved exactly: all inconsistent;
   the full 556-direction screen is running on the remote (excl-misc).
10. No witness found; nothing here claims nonexistence of a certificate.  Packet: /tmp/hodge/results/quadric-contact-residual/wave2-excl/.

## Rigor conventions
* A modular rank is a lower bound for the rational rank; "independent columns force zero" transfers to Q and C.
* "Only the origin" for a homogeneous quadric ideal computed mod p (msolve -g 1 leading ideal containing a pure power of every
  variable, or Singular dim 0) transfers to Q-bar by properness of Proj(Z[c]/I) -> Spec Z (empty special fibre => empty generic
  fibre), PROVIDED the mod-p system is the reduction of the rational one.  For the Hessian relaxation this holds because dim J_m
  mod p equals the rational Hilbert function of the Jacobian ring (checked in every run), so the annihilator lattice of J_m reduces
  correctly and the mod-p quadrics are reductions of integral quadrics.
* Hessian relaxation (necessary condition): for D = sum c_k D_k in the E1 kernel with E1 lifts N_k, every certificate has
  N = sum c_k N_k + Koszul and B_D(N) = F Q in J_{2d}, hence sum_{k,l} c_k c_l B_{D_k}(N_l) in J_{2d}; Koszul terms and F*Q drop
  mod J; gauge directions s*F (with N = A s Z/6) drop as well, so their coefficients are set to 0.
* For D = C(F) in J_d (C a vector field of degree d-5) the relaxation reads A*C^T Hess(F) C in J_{2d}, i.e.
  C^T Hess(F) C in (J:A)_{2d-6}; at d = 7 this is C^T Hess(F) C in (J:A)_8 = J_8 + <D8*>, which gives exactly 950 = 951 - 1
  independent quadrics (dim R_8 = 951).

## Results (all at P2 unless stated; "two primes" = also at 1000000007)
1. [replayed, agrees] E1 colon dims dim (J:A)_d = 0,0,0,0,0,6,36,126,337,760 for d = 0..9; the non-Jacobian E1 directions below
   degree 10 are D8* = 4G_S^2 - G_P^2 (d = 8, trivial sector) and X1 D8*, Y1 D8*, X2 D8*, Y2 D8* (d = 9).  [exp E1-colon]
2. [replayed, agrees] d = 5: 21 Hessian columns rank 21 in R_10, Singular dim 0.  d = 6: 35 gauge-fixed unknowns, 630 pair
   columns of rank 334 in R_12, the 28 off-diagonal m_ij except (0,3),(3,0) killed by square columns, Singular dim 0 and msolve
   m-primary.  => deg D <= 6 EXCLUDED over C.  [exp closure-d5, closure-d6]
3. [NEW, two primes] d = 7: the 12 subfamilies in which two of the three involution characters (b0 = sigma0, b2 = sigma2,
   y = tau_Y) are fixed and everything else (incl. all six h1-characters) mixes: 24-32 unknowns, 227-258 quadrics, msolve -g 1
   m-primary (max pure power 3-4), Singular dim 0 on d7_b01b21 (vdim 209) => EXCLUDED over C.  Consequence: a degree-7
   denominator must mix both eigenvalues of at least two of sigma0, sigma2, tau_Y.  [exp sub-closure-msolve, tags d7_*, d7p3_*]
   Pending (remote excl-sub7): the 6 subfamilies with ONE involution character fixed (60 unknowns, 477 quadrics, Froeberg
   regularity 4; msolve -g 1 with block-Macaulay k=4 fallback).
4. [NEW, two primes] d = 7, the parity-class cube: the 8 classes (b0,b2,y) in (Z/2)^3 are the vertices of a cube; item 3 covers
   the 12 edges (pairs differing in one coordinate; fqlx's single vertices are contained).  The 12 face diagonals and the 4 space
   diagonals (pairs differing in two or three coordinates; 28-32 unknowns, 227-246 quadrics) are also EXCLUDED (msolve -g 1,
   max pure power 3-4).  Hence: a degree-7 denominator must have G_aff-character support on at least THREE of the eight parity
   classes.  [exp sub-closure-msolve, tags d7_fd_*, d7_sd_*]  Larger supports: see item 18 (bilinear lemma); the brute-force
   rungs (56 triples ~45 unknowns, 70 quadruples ~60, 56 quintuples ~75) are Froeberg regularity 4, sextuples (90) regularity 5.
5. [computed, two primes] the full d = 7 system (120 unknowns): 950 independent quadrics, mean 58 terms (max 298), no square
   c_a^2 in I_2, only 8 product monomials c_a c_b in I_2 (the same 8 at both primes); branch-and-kill on I_2 does not cascade
   (one kill per node).  Froeberg count for 950 quadrics in 120 unknowns: HS = 1, 120, 6310, 181240, 2632405, <0 => a
   Nullstellensatz certificate needs degree >= 5, i.e. Macaulay matrices of size ~2.8e8 x 2.3e8: INFEASIBLE.  The same count
   for the 60-unknown single-involution subfamilies gives regularity 4 (feasible), for 32 unknowns regularity 3.
6. [tooling] msolve misparses `c*c` (not read as a square) and crashes in FGLM for primes >= ~2^31; wave-1's msolve verdicts in
   this lane were on misparsed systems and are discarded; Singular runs are unaffected.  Fix: `c^2` syntax and `msolve -g 1`
   (leading ideal; m-primary test).  Recorded in /tmp/hodge/notes/engine-bugs.md.  target_image validated on collision-h3 and
   hesse-triangle (A in served(D_witness)).
7. [NEW, two primes] d = 8, Gamma-invariant family (G_aff and the exchanges s0, r1; contains D8* and the gauge F*S_2^Gamma):
   20 parameters = 3 gauge + 16 Jacobian + 1 extra (D8*), 17 unknowns, 52 quadrics of rank 35 (wave-1 numbers reproduced);
   msolve -g 1 m-primary (max pure power 5) => EXCLUDED over C.  (Closes wave-1's pending item.)  [exp closure-Gamma-invariant]
8. [NEW, two primes] d = 8, pure G_aff character sectors (gauge-fixed full E1-admissible family of the sector): 20 of 24 sectors
   EXCLUDED — all sectors with h1-character != 0, i.e. (b0,b2,h,y) with h in {1,...,5}: unknowns 4-17, quadrics 37-52.
   Pending (remote excl-d8large): (0,0,0,0) [29 unknowns], (0,0,0,1) [27], (1,1,0,0) [26], (1,1,0,1) [26] (52 quadrics each;
   Froeberg regularity 6-7).  [exp closure-pure-sector, d=8]
9. [NEW] served-target spaces (sum principle: served(D1 D2) = served(D1) + served(D2), served(D^k) = served(D)): for the 6
   coordinate hyperplanes, 11 quadrics (S0, delta = X0^2-Y0^2, S1, S2, S0+-S2, S0+S1, X0^2+X2^2, X0^2-X1^2, X1^2-Y1^2, X2^2-Y2^2),
   7 quartics (G_S, G_P, Q1, Q2, X0^4+Y0^4, S0^2+2P0^2, X1^4+Y1^4), the six polars dF/dZ_i, f0 = S0^3 + 2 S0 P0^2, the three
   Euler-type sextics rad_i = X_i F_Xi + Y_i F_Yi and the Fermat blocks X1^6+Y1^6, X2^6+Y2^6 (remote, P2):
   AND the degree-8 contour divisor D8* itself (remote, 40326 x 7723 system, 1113 s): served(D) = J_6^odd (dim 18, the trivially
   served targets) in every case — in particular V(D8*) ∩ X supports NO nontrivial residue class.  The total span is J_6^odd and
   A is not in it => NO product of these blocks (any degree, any multiplicities) serves A or any other nontrivial target.
   Rigor: served_p contains J_6^odd (dim 18) and dim served_p >= dim served_Q (rank can only drop mod p), so served_Q = J_6^odd exactly.
   [exp served, served-span; served_local.json]
10. [replayed at P2, agrees] sector-restricted full fixed-D solves (own implementation) of D8* (rank 718, nullity 109) and
    X1 D8*, Y1 D8* (rank 926), X2 D8*, Y2 D8*, X0 D8*, Y0 D8* (rank 1098): all inconsistent, identical ranks to wave 1.  With the
    prescribed-support theorem and item 9 this excludes D8*^k times any product of the item-9 blocks at every degree.
11. [computed] d = 9 pure sectors, (unknowns, quadrics) after gauge: (38,14)x4, (15,30)x2, (53,30)x2, (19,14)x8, (52,30)x2,
    (34,14)x4, (12,30)x2 — in 20 of 24 sectors the relaxation has fewer quadrics than unknowns (dim R_18 = 426 in total), so its
    solution set is positive-dimensional and it cannot close them; from degree 9 on exhaustive exclusion needs the full bilinear
    E1+E2 system.  Froeberg regularities for degree 8 beyond pure sectors: parity classes (~40 unknowns, ~120 quadrics): 5;
    two involutions fixed (~80, ~240): 8; full (316, 951): 14 — only pure sectors and (with effort) parity classes are in reach.
12. [NEW, two primes] d = 9: the four sectors with enough quadrics, (0,0,3,0), (0,0,3,1) [15 unknowns, 30 quadrics] and
    (1,1,3,0), (1,1,3,1) [12, 30], are EXCLUDED (msolve -g 1 m-primary).  [exp closure-pure-sector, d=9]
13. [screen] d = 10, trivial sector: 20 non-Jacobian E1 directions (kernel 128 = 108 + 20); the first 7 basis elements solved
    exactly (sector-restricted full E1+E2 system, 1934 unknowns, rank 1667): inconsistent.  Not a family exclusion; the
    remaining 549 extras at d = 10 are being screened by excl-misc.
14. [second method] Singular dim 0 on the msolve-closed systems d7_b01b21 (vdim 209), d7_fd_00y0_11y0 (288), d7_sd_000_111 (288),
    d9_chi1130 (65), d8_chi0130 (5); Singular times out (80 s) on the 17-unknown Gamma-invariant system where msolve's F4 needs 3 s.

15. [replayed exactly over Q] fqlx's weighted-K3 Gauss-Manin obstruction: with Singular in characteristic 0, weights (1,1,2,2),
    t = (2,3): Hilbert function of R(F_sub) = 1,2,5,8,12,14,16,14,12,8,5,2,1,0 and A*dF_sub/dt1, A*dF_sub/dt2 are NOT in
    J_12(F_sub) (the socle pairings are nonzero) — agrees with fqlx §4.  [/tmp/mstest/k3.sing]
16. [replayed, agrees] dossier's epsilon-sector dimensions (G_aff character of A, s0- and r1-invariant): degree 6: 25/3/22,
    degree 12: 191/140/51 (S/J/R).  [exp epsilon-sector-dims]
17. [replayed / corrected] dim (J:A)_10 = 2053 (556 non-Jacobian directions; dim R_16 = 951, i.e. xA: R_10 -> R_16 has corank 1)
    agrees with the dossier.  The dossier's remark "(J:A)_d = S_d for d >= 13 since R_19 = 0" is wrong as stated: R_19 has
    dimension 246, and by Gorenstein duality xA: R_13 -> R_19 is dual to xA: R_5 -> R_11, which is injective (item 1), so
    dim (J:A)_13 = 8568 - 246 = 8322 — CONFIRMED numerically at P2 (1440 non-Jacobian directions; exp E1-colon d=13, 471 s).
    Correct statement: (J:A)_d = S_d only for d >= 19.

18. [lemma; hypothesis REFUTED by computation] Difference-multiplicity lemma for the degree-7 cube.  Write c = sum_s c_s over
    the parity classes s in (Z/2)^3 and let S = {s : c_s != 0}.  The relaxation splits by (Z/2)^3-character: for u != 0,
    E_u : sum_{ {a,b} subset S, a+b = u } H(c_a, c_b) = 0, and E_0 : sum_s H(c_s, c_s) = 0.  If a pair {a,b} subset S is the
    only pair in S with difference a+b, then H(c_a, c_b) = 0 on its own (all pairs when |S| = 3; all six pairs of the 56
    non-planar 4-sets; the four spoke pairs of every 5-set).  This would reduce the possible supports to the 14 affine planes
    and the 6-,7-,8-sets IF no class pair (s,t) admitted a "zero pair" (c_s != 0, c_t != 0, H(c_s,c_t) = 0).
    COMPUTED: zero pairs DO exist.  For (s,t) = ((0,0,0),(0,0,1)) the bidegree-(2,2) Macaulay test does not close (bilinear
    rank 103 of 256 cross monomials; 637 s), and the h1-micro-block analysis explains why: the block V_s[h=1] x V_t[h=3]
    (7 x 1) has rank 6 < 7 with a rank-one kernel tensor, likewise (h=3,h=5); for ((0,1,0),(1,0,0)) the (0,0) block is 12 x 12
    of rank 45.  So the cross terms alone do not exclude larger supports; the lemma only says that on a 3-support the three
    cross terms vanish separately, i.e. (c_a,c_b) must lie on the (non-empty) zero-pair varieties Z_ab — a strong but not
    conclusive constraint.  The brute-force rungs remain: 56 triples (~45 unknowns), 14 affine planes (60; 6 half-cubes in
    excl-sub7, 8 diagonal planes in excl-planes7), quadruples/quintuples (regularity 4), sextuples and beyond (regularity 5).
    [exp bilinear-pair; bilinear_pairs2.py; micro-block analysis in the session log]

## What remains open / recommended next steps
* Degree 7: the six single-involution subfamilies (running); the fully mixed 120-unknown problem needs a new idea — more
  quadratic consequences of E2 than the 950 mod-J ones (the Koszul image B_D(Kos) + F S_8 has dimension <= 2127 << 7062 = dim J_14,
  so the mod-J relaxation discards most of E2, but the extra conditions depend on D), or a geometric argument on
  Z = X ∩ V(C(F)) (non-normal / non-rational singularities needed for a (3,1) support).
* Degree 8: the four large pure sectors (running); the 8 parity classes (regularity 5, ~1.1e6 monomials of degree 5) are the
  next feasible target; the gauge-fixed full family (316 unknowns) is hopeless for Groebner (regularity 14).
* Degree >= 9: only screens and structured families; the served-space method (item 9) is the right tool for structured
  supports and should be extended to the non-eigenvector blocks (needs a mixed-parity fixed-D solver) and to the remaining
  degree-6/8 blocks (running).
