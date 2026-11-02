# Wave 2 — lane `quadratic-cubic-pullback` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness.  Rigorous so far: K_d = J(F)_d re-verified at prime 2147483629 for d in {3,5,6(odd)} (local smoke), NO Y-even polar denominator of degree 5 exists (exact relaxation, pi_W(ker)=0 over Z), all wave-1 structural claims reconfirmed at a second prime.  Remote jobs (kernels d<=12, 58 structured polars, relaxation d=5/6) queued (remote saturated at 19:30).

Started 2026-11-01 (wave 2).  Log: `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl.log.jsonl`.
Packet: `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/`.

## Results so far
- Local check: F = sum Q_i^3 exactly; the fixed-D kernel of a polar denominator D = dF/dX0 (d=5) has the explicit
  rational basis (G z, 0), G in S_5^odd (126) and (g F, D div g - g.grad D), g constant on the Y-slots (3): 129 elements,
  equal to the modular nullity (engine reports 128 = rank of augmented matrix convention).  Consequence: the three
  wave-1 degree-5 polar inconsistencies are rigorous over Q (modular rank of M = rational rank, sandwich argument).
- Engine note written to /tmp/hodge/notes/engine-bugs.md (reported nullity off by one for inconsistent systems).
- Dossier correction: E_I F = sum_{i in I} z_i dF/dz_i lies in J(F) for EVERY index set I (not only I = all), so the
  block-Euler denominators EI([0,1,2]) etc. listed in wave-1 `listC` were never tested; they are natural
  degree-6 polar candidates (the Fermat join X0^6+X1^6+X2^6 is exactly E_X F_0 / 6) and are queued now.

## Re-verification of wave-1 claims (goal a) — all at prime p2 = 2147483629 (wave 1 used 2147483647)
| claim | wave 1 | wave 2 (p2) | status |
|---|---|---|---|
| h^{3,1}(X)^+ = 226 - 18 = 208; A, delta not in J_6, independent mod J_6 | yes | J6_odd rank 18 = Hilbert, A,delta independent | rigorous (agrees) |
| projective diagonal symmetry group = Z/2 | Smith form (sympy) | engine `symmetry.diagonal_group`: order 2 | agrees |
| F not a join: rank of 21 second partials = 21 | p1 | p2: 21 | rigorous (agrees) |
| cubics in Q meet J_6^even in exactly C.F | p1 | p2: intersection dim 1 | agrees (mod-p dim is an upper bound, F gives the lower bound) |
| two-plane (4,4)-cancellation system has nullity 0 | sympy exact | fmpq exact: rank 12 of 12 | rigorous (agrees) |
| delta, join012, join345, F60, F06, fermatX, Q0^3+Q3^3, Q0^2, F-Q0^3 not in J_6 (so excluded at every pole order by K_6=J_6) | p1 | p2: all not in J, J-rank certified | rigorous (agrees) |
| E1 kernel K_d = J_d (even d<=8, odd d<=6) | p1 | p2 (smoke: d=3,5 both; 6 odd) — remote job covers d<=12 | agrees so far |
| three degree-5 polar denominators serve nothing | p1, modular | structural kernel 129 = modular nullity => rigorous over Q | upgraded to rigorous |
| E_X F etc. "in J only for I = all" (dossier §3.3) | — | E_X F is in J_6 (certified) | dossier statement corrected |

## New rigorous results
* **No Y-even degree-5 denominator** (local, exact): polar family D = sum_{i<3} w_i F_i, relaxation matrix 1491 x 201
  (6 W + 135 V + 60 Q unknowns), exact integer nullspace has dimension 12 and projects to ZERO on the W-coordinates,
  so no w gives a certificate.  Self-test of the E2 reduction E = F (Phi(w) + L(c,w)) passed exactly.
* gradX0sq = 2 X0 F_0 (d = 6, reducible control): serves nothing, rank-certified (nullity 244 = structural).
