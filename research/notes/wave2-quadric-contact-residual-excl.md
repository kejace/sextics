# Wave 2 — quadric-contact-residual — lens EXCLUSION / RIGOR

STATUS: in progress — no witness; rigorous bound: deg D >= 7 (fqlx; d<=6 independently replayed at P2). NEW at d=7: every subfamily with two of the three involution characters (sigma0, sigma2, tau_Y) fixed is excluded (12 subfamilies, msolve -g 1, m-primary). Full mixing at d=7 is beyond Groebner reach (Froeberg regularity >= 5 for 950 quadrics in 120 unknowns).

Lane: q = 1, A = X0^3 Y0^3, F = quadric-contact sextic (see dossier).  New prime P2 = 2147483629 (wave 1: 2147483647; dossier: 1000000009; Singular: 32003).
Code: /tmp/hodge/results/quadric-contact-residual/wave2-excl/{w2lib.py, family.py, reverify.py, fixedD_batch.py, d7mixed.py, msverdict.py}
Log: /tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl (+ per-job logs wave2-excl/excl-*.log.jsonl merged after pull)

## Rigor conventions
* A modular rank is a lower bound for the rational rank; "independent columns force zero" transfers to Q and C.
* "Only the origin" for a homogeneous quadric ideal computed mod p (Singular dim 0 / msolve leading ideal m-primary) transfers
  to Q-bar by properness of Proj(Z[c]/I) -> Spec Z (empty special fibre => empty generic fibre), PROVIDED the mod-p system is the
  reduction of the rational one.  For the Hessian relaxation this holds because dim J_m mod p equals the rational Hilbert
  function (checked in every run: the Jacobian ideal of F mod p is still a complete intersection with the same Hilbert series),
  so the annihilator lattice of J_m reduces correctly.
* Hessian relaxation (necessary condition): for D = sum c_k D_k in the E1 kernel with E1 lifts N_k, every certificate has
  N = sum c_k N_k + Koszul, and B_D(N) = F Q in J_{2d}, hence sum c_k c_l B_{D_k}(N_l) in J_{2d}; gauge directions s*F drop.

## Results so far (all at P2 unless stated)
1. [replayed, agrees] d = 5: 6 unknowns, 21 Hessian columns rank 21 in R_10, all squares forced, Singular dim 0 -> EXCLUDED.
2. [replayed, agrees] d = 6: 35 unknowns (gauge m_55 = 0), rank 334 among 630 pair monomials, 28 square-column kills
   (exactly the off-diagonal m_ij except (0,3),(3,0)), Singular dim 0 -> EXCLUDED over C.  msolve -g 1: m-primary leading ideal.
3. [tooling] msolve misparses `c*c` squares (see /tmp/hodge/notes/engine-bugs.md); wave-1 msolve verdicts in this lane
   (closure_quads_d8*.ms, dossier d7_system.ms) were computed on misparsed systems and must be discarded; Singular runs unaffected.

4. [NEW, rigorous] d = 7, the 12 subfamilies in which two of (b0, b2, y) are fixed (all six h1-characters and the third
   involution character mixed): 24-32 unknowns, 227-258 quadrics each, msolve -g 1 leading ideal m-primary (max pure power 3-4)
   => only c = 0 => EXCLUDED over C (log: exp = sub-closure-msolve, tags d7_b00b20 ... d7_b21y1).  Consequence: a degree-7
   denominator must mix both values of at least two of the three characters sigma0, sigma2, tau_Y.
5. [computed] full d = 7 system (120 unknowns): exactly 950 independent quadrics (= dim R_8 - 1 = 951 - 1, via the degree-8
   formulation C^T Hess F C in (J:A)_8 = J_8 + <D8*>), no square c_a^2 in I_2, only 8 product monomials c_a c_b in I_2;
   branch-and-kill on I_2 does not cascade.  Froeberg count for 950 quadrics / 120 vars: HS(S/I) = 1, 120, 6310, 181240,
   2.9e6, <0  => a Nullstellensatz certificate needs degree >= 5 (Macaulay matrices ~2.8e8 x 2.3e8): infeasible.
6. [tooling check] target_image (served-target space, sum principle) validated on collision-h3 and hesse-triangle (A in served(D_witness)).

7. [NEW, rigorous] d = 8, Gamma-invariant family (G_aff and s0, r1 invariant; contains D8* and the F*S_2^Gamma gauge):
   20 parameters = 3 gauge + 16 Jacobian + 1 extra (D8*), 17 unknowns, 52 quadrics of rank 35 (wave-1 numbers reproduced at P2);
   msolve -g 1: m-primary leading ideal (max pure power 5) => EXCLUDED over C.  This closes wave-1's pending Groebner item
   (wave-1's own msolve runs were on misparsed input).
8. [NEW, rigorous] d = 8, pure G_aff character sectors: 20 of the 24 sectors EXCLUDED (gauge-fixed full E1-admissible family,
   4-17 unknowns, 37-52 quadrics, msolve -g 1 m-primary; log exp = closure-pure-sector).  Remaining: (0,0,0,0) [29 unknowns],
   (0,0,0,1) [27], (1,1,0,0) [26], (1,1,0,1) [26] (52 quadrics each; Froeberg regularity ~5-6; queued for the remote).

9. [NEW, rigorous] served-target spaces (sum principle): for the 6 coordinate hyperplanes, 11 quadrics (S0, delta, S1, S2,
   S0+-S2, S0+S1, X0^2+X2^2, X0^2-X1^2, X1^2-Y1^2, X2^2-Y2^2), 7 quartics (GS, GP, Q1, Q2, T0, S0^2+2P0^2, X1^4+Y1^4), the six
   polars F_Zi and f0 = S0^3 + 2 S0 P0^2: served(D) = J_6^odd (dim 18, the trivially served targets).  Span of all of them = J_6^odd,
   A not in it => no product of these blocks (any degree) serves A.  [served_local.py; log exp = served, served-span]
   Pending on the remote: rad_i, X1^6+Y1^6, X2^6+Y2^6, D8* (degree 6 and 8 blocks).

10. [computed] Degree-9 pure sectors: (unknowns, quadrics) per sector after gauge = (38,14)x4, (15,30)x2, (53,30)x2, (19,14)x8,
    (52,30)x2, (34,14)x4, (12,30)x2 — in 20 of 24 sectors there are fewer quadrics than unknowns, so the Hessian relaxation
    (quadrics in R_18, dim 426 total) cannot close them (solution sets are positive-dimensional); only the four (.,.,3,.) sectors
    can be closed by it.  From degree 9 on, exhaustive exclusion needs the full bilinear E1+E2 system.
11. [computed] Froeberg regularity (first degree where the generic Hilbert series of m quadrics in n unknowns turns <= 0):
    d7 full (120, 950): 5; d7 tau-even (60, 477): 4; d7 two-involutions-fixed (32, 241): 3; d8 Gamma-invariant (17, 52): 4;
    d8 sectors (29,52): 7, (27,52): 6, (26,52): 6.  These bound the Groebner effort from below (generic case).

## Running (remote)
* excl-reverify: E1 colon dims d<=9, d=7 parity sectors, all 24 pure G_aff sectors at d=8 and d=9 (Singular verdicts).
* excl-fixedD: sector re-solves of D8*, Z_i*D8* at P2; engine cross-check; served-target spaces of structured blocks (sum principle).
* excl-d7mixed: 120-unknown degree-7 mixed-sector system (degree-8 Hessian formulation), I_2/I_3 monomial analysis, msolve -g 1.
