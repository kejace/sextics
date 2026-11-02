# Wave 2 — quadric-contact-residual — lens EXCLUSION / RIGOR

STATUS: in progress — no witness; rigorous bound so far: deg D >= 7 (fqlx, now independently replayed at d <= 6); degree-7 full closure running.

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

## Running (remote)
* excl-reverify: E1 colon dims d<=9, d=7 parity sectors, all 24 pure G_aff sectors at d=8 and d=9 (Singular verdicts).
* excl-fixedD: sector re-solves of D8*, Z_i*D8* at P2; engine cross-check; served-target spaces of structured blocks (sum principle).
* excl-d7mixed: 120-unknown degree-7 mixed-sector system (degree-8 Hessian formulation), I_2/I_3 monomial analysis, msolve -g 1.
