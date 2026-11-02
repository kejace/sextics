# Wave 2 search report: `cyclic-triple-residual` (q = 1, A = X0^3 Y0^3)

STATUS (2026-11-01 20:02 local): **No witness found.** Rigorous this wave: no certificate with deg D <= 5 (d <= 4: no E1-admissible D; d = 5: relaxation quadrics span all quadratic forms in both parities); no Y-odd denominator of degree 6 (HF_4 = 0); for Y-even degree 6 the relaxation variety modulo the F-gauge is at most ONE projective point (HF_4 = 1 on the complement of the Euler line; E_X(F) is not it), identification (HF_5 + numerical search + exact tests) was running locally when the session was closed. E1-admissible spaces d = 4..10 computed exactly per G-sector (bases saved as `w2-e1space-d*-c*.json` in the remote mirror). Structured sweep: 70/102 candidates finished, none consistent (every degree-6 structured candidate fails E1; the Jacobian-type degree-8 pieces pass E1 and fail E2 exactly). Remote joint searches (Fermat calibration, Jacobian-type d = 5..7, sector searches d = 8, 9, served-space sum-principle test) were written, tested, and queued but did not start: the shared 128-core box ran at load 700-940 from other agents and the runner refuses submissions above load 150; the follow-up submission scripts (`w2_submit3.sh`, `w2_submit_served.sh`, `w2_submit_next.sh`) are in the lane dir and resume the plan unchanged.

Code (all in `/tmp/hodge/results/cyclic-triple-residual/`): `w2core.py` (exact E1-admissible spaces per G-sector via fmpz nullspace; Koszul-reduced, sector-restricted variable-projection search; mod-J relaxation; reconstruction + exact confirmation through the wave-1 sector solver `sect.py`), `w2run.py` (multiprocessing driver), jobs `w2_e1space.py`, `w2_jac.py` (+ `w2_jacjob.sh`), `w2_sector.py`, `w2_sweep.py`, `w2_served.py` (served target spaces / sum principle), `w2_macblocks.py` (rigorous relaxation exclusion by G-graded block Macaulay ranks), `w2_macaulay.py` + `w2_hilb.py` (msolve cross-check), `w2_report.py`/`w2_mkreport.py` (this report). Canonical log: `/tmp/hodge/results/cyclic-triple-residual/wave2-search.log.jsonl` (merged from per-process `w2log-*.jsonl` files by `w2_pull.sh`; remote mirror `/tmp/hodge/remote/cyclic-triple-residual/`).

## Method
Split form E1: q N(F) = A D, E2: D div N - N(D) = F Q. Every E1 lift is N0 + K(c), K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j), deg c_ij = d - 4 (book ch. 5.2).
For a linear family D(a) = sum a_k D_k with exact lifts N0_k the E2 residual r = Bil(D(a), N0(a) + K(c)) - F Q is linear in (c, Q):
f(a) = min_{c,Q} |r|^2 on the unit sphere in a (variable projection, analytic Golub-Pereyra gradient, L-BFGS-B random restarts; the subspace of F-multiples is projected out, which is lossless by the gauge (D, N) -> (D + lF, N + A l Z/6); rel = f/|Bil(D, N0)|^2).
Mod-J ("Hessian") relaxation: f_J(a) = |P_perp^T Bil(D(a), N0(a))|^2 with P_perp an orthonormal basis of the complement of J(F)_{2d} in the E2 row space -- a necessary condition (im T_D and F S lie in J).
Character sectors of the order-18 diagonal group G (generators mu3 = (0,0,2,0,0,0), s = (1,5,1,2,2,0); labels (m.mu3, m.s) mod 6; chi_A = (0,3); chi(N_i) = chi_A + c + chi(e_i), chi(Q) = chi_A + 2c).
Rigorous exclusion at Jacobian-type degrees (d <= 7, where Ann([A])_d = 0 so every E1-admissible D is w.grad F): the relaxation is a system of dim R_{2d}^- quadrics q_s(w) in the w-space; the degree-k Macaulay matrix is block-diagonal by the G-character of w; if every block has full column rank mod p (p good for J_{2d}: the mod-p rank of J_{2d}^- equals the Hilbert-series value), the ideal contains all degree-k forms over F_p, hence over Q, so w = 0 is the only solution. The spurious locus L = {l Z/6} (D = lF) consists of solutions; V(q) is invariant under translation by L (gauge), so restricting to a complement of L (dropping the coordinates (0, m) with X0 | m) is lossless. Macaulay's bound HF(k+1) <= HF(k)^<k> turns HF_k = 1 into "at most one projective point".

Jacobian ring dimensions (Hilbert series (1+t+..+t^4)^6, Y-odd parts): R_10^- = 744, R_12^- = 866, R_14^- = 744, R_16^- = 468, R_18^- = 208, R_20^- = 60, R_22^- = 9, R_24^- = 0.

## Self-test (local)
Vectorized residual vs direct evaluation: 1.3e-16. Gradient FD check: 8e-10 / 1.7e-8. Fermat control (F0 = Fermat sextic, d=6, sector (0,0), family dim 31): the known witness direction (X0^6+X1^6+X2^6, gauge-shifted by F0/2) has residual 0 in the family; 6 random restarts missed it (best rel 1e-2) -> hit-rate calibration with 400 restarts is the first batch of the jac job.

## Results
### Rigorous consequences (this wave; all mod-p statements are of the "rank over F_p <= rank over Q" type, so they are theorems over Q/C)

1. deg D <= 4: no E1-admissible denominator at all (Ann([A])_d = 0 and J_d = 0), confirmed by the exact E1-space computation.
2. deg D = 5 (both Y-parities): every E1-admissible D is a polar combination a.grad F; the relaxation quadrics already span all quadratic forms in the 3 unknowns (HF_2 = 0 in each parity), so no certificate with deg D = 5 exists.
3. deg D = 6, D Y-odd: E1-admissible D = w.grad F with w a Y-odd/even linear field (18 unknowns); the degree-4 Macaulay matrix of the 866 relaxation quadrics has full column rank mod p (HF_4 = 0): no Y-odd denominator of degree 6.
4. deg D = 6, D Y-even: the solutions of the relaxation form, modulo the F-gauge line (D = cF), a scheme with Hilbert function HF_3 = 11, HF_4 = 1, hence (Macaulay bound) at most one projective point; E_X(F) = sum_{i<3} X_i dF/dX_i (gauge-equivalent to E_Y(F)) is not that point (Bil(E_X F, A X) is a 5-term polynomial -- A times the five mixed couplings -- not in J_12, exactly over Q). See the identification entry below for the status of this single candidate point.

### E1-admissible spaces by degree (G-sectors; dims mod p = rigorous upper bounds, bases exact over Q)

| d | dim S_d | dim J_d | kerA in S_d | ker mod J | dimension-forced | sporadic | sectors with ker mod J > 0 | time |
|---|---|---|---|---|---|---|---|---|
| 4 | 126 | 0 | 0 | 0 | 0 | 0 | - | 1.2s |
| 5 | 252 | 6 | 6 | 0 | 0 | 0 | - | 2.4s |
| 6 | 462 | 36 | 36 | 0 | 0 | 0 | - | 3.8s |
| 7 | 792 | 126 | 126 | 0 | 0 | 0 | - | 8.1s |
| 8 | 1287 | 336 | 337 | 1 | 0 | 1 | 04:1 | 20.1s |
| 9 | 2002 | 756 | 876 | 120 | 116 | 4 | 00:19, 01:20, 02:17, 03:19, 04:16, 05:19, 25:1, 40:3, 41:1, 42:3, 44:2 | 83.0s |
| 10 | 3003 | 1497 | 2053 | 556 | 555 | 1 | 00:52, 01:50, 02:53, 03:51, 04:53, 05:52, 20:19, 21:19, 22:20, 23:17, 24:19, 25:16, 40:24, 41:18, 42:27, 43:21, 44:27, 45:18 | 210.9s |

### Rigorous relaxation exclusions (block Macaulay ranks mod p; HF_k = 0 is rigorous over Q)

| d | Y-parity | restricted to complement of F-gauge locus | #w vars | quadric span rank | k | HF_k | rigorous exclusion | time |
|---|---|---|---|---|---|---|---|---|
| 5 | 0 | False | 3 | - | 3 | 0 | YES | 2.4s |
| 5 | 1 | False | 3 | - | 3 | 0 | YES | 2.5s |
| 6 | 0 | False | 18 | - | 3 | 80 | no | 23.8s |
| 6 | 1 | False | 18 | - | 3 | 2 | no | 23.1s |
| 6 | 0 | False | 18 | 102 | 4 | 81 | no | 49.8s |
| 6 | 1 | False | 18 | 134 | 4 | 0 | YES | 53.3s |
| 6 | 0 | True | 17 | 102 | 3 | 11 | no | 18.2s |
| 6 | 0 | True | 17 | 102 | 4 | 1 | no | 49.9s |
| 6 | 1 | False | 18 | 134 | 3 | 2 | no | 29.0s |

- relaxation-point-d6: {"d": 6, "parity": 0, "Bil_EX_in_J12_modp": false, "Bil_EX_in_J12_exactQ": false, "Bil_EX_terms": 5, "Bil_EX_degree": 12, "exact_ranks": "J_12 sector generators rank 365, with Bil(E_X F, A X) 366 (fmpz exact)", "relax_rel_at_EX": 0.10121255560322683, "fullE2_rel_at_EX_numerical": 0.20922440273498893, "note": "d=6 Y-even: relaxation variety modulo the F-gauge line is at most one projective point (HF_4 = 1, Macaulay bound); E_X(F) = sum X_i dF/dX_i (gauge-equivalent to E_Y(F)) is NOT that point: Bil(E_X F, A X) = A * (five mixed couplings) is not in J_12 (exactly over Q). Values recorded from the printed output of the local check (the logging call was cut by the 110 s local timeout)."}

### Rational reconstructions tested exactly

| label | d | descr | recon err | num rel | E1 | consistent | rank/aug | witness |
|---|---|---|---|---|---|---|---|---|
| fermat-control | 6 | coef-limit_den(1) | 0.4858 | 0.01087 | True | False | 44/45 | no |
| fermat-control | 6 | coef-limit_den(6) | 0.04378 | 0.01087 | True | False | 44/45 | no |
| fermat-control | 6 | coef-limit_den(12) | 0.01347 | 0.01087 | True | False | 44/45 | no |
| fermat-control | 6 | coef-limit_den(60) | 0.001967 | 0.01087 | True | False | 44/45 | no |
| fermat-control | 6 | coef-limit_den(420) | 7.51e-05 | 0.01087 | True | False | 44/45 | no |
| fermat-control | 6 | coef-limit_den(2520) | 2.44e-07 | 0.01087 | True | False | 44/45 | no |

### Structured candidates (sweep; exact sector solver for d <= 9, two primes above)

| family | candidate | d | terms | F|D | E1 | consistent | rank/aug (cols) | exact | time |
|---|---|---|---|---|---|---|---|---|---|
| polar-pencil | polar F_a, a = e_Y1 | 5 | 6 | False | True | False | 256/257 (257) | True | 7.2s |
| polar-pencil | polar F_a, a on line {X0=Y1=X1=0} | 5 | 8 | False | True | False | 713/714 (715) | True | 9.3s |
| polar-pencil | polar F_a, a=e_X0+e_X1 | 5 | 10 | False | True | False | 723/724 (726) | True | 9.9s |
| F-pieces | F - F|_{X2=0} | 6 | 4 | False | False | False | None/None (None) | None | 23.1s |
| F-pieces | Fermat part sum Z_i^6 | 6 | 6 | False | False | False | None/None (None) | None | 22.6s |
| F-pieces | coupling part (F - Fermat) | 6 | 6 | False | False | False | None/None (None) | None | 20.5s |
| F-restriction-1 | F|_{X1=0} | 6 | 8 | False | False | False | None/None (None) | None | 21.4s |
| F-restriction-1 | F|_{Y1=0} | 6 | 6 | False | False | False | None/None (None) | None | 21.8s |
| F-restriction-1 | F|_{X2=0} | 6 | 8 | False | False | False | None/None (None) | None | 23.7s |
| F-restriction-1 | F|_{Y0=0} | 6 | 8 | False | False | False | None/None (None) | None | 25.3s |
| F-restriction-1 | F|_{X0=0} | 6 | 6 | False | False | False | None/None (None) | None | 28.2s |
| F-restriction-1 | F|_{Y2=0} | 6 | 8 | False | False | False | None/None (None) | None | 23.4s |
| F-restriction-2 | F|_{X0,X2=0} | 6 | 4 | False | False | False | None/None (None) | None | 21.9s |
| F-restriction-2 | F|_{X0,X1=0} | 6 | 4 | False | False | False | None/None (None) | None | 26.5s |
| F-restriction-2 | F|_{X0,Y0=0} | 6 | 4 | False | False | False | None/None (None) | None | 20.2s |
| F-restriction-2 | F|_{X0,Y1=0} | 6 | 4 | False | False | False | None/None (None) | None | 22.4s |
| F-restriction-2 | F|_{X0,Y2=0} | 6 | 5 | False | False | False | None/None (None) | None | 21.5s |
| F-restriction-2 | F|_{X1,X2=0} | 6 | 5 | False | False | False | None/None (None) | None | 21.9s |
| F-restriction-2 | F|_{X1,Y0=0} | 6 | 6 | False | False | False | None/None (None) | None | 23.9s |
| F-restriction-2 | F|_{X1,Y2=0} | 6 | 5 | False | False | False | None/None (None) | None | 21.0s |
| F-restriction-2 | F|_{X1,Y1=0} | 6 | 5 | False | False | False | None/None (None) | None | 23.7s |
| F-restriction-2 | F|_{X2,Y0=0} | 6 | 5 | False | False | False | None/None (None) | None | 24.6s |
| F-restriction-2 | F|_{X2,Y1=0} | 6 | 4 | False | False | False | None/None (None) | None | 21.6s |
| F-restriction-2 | F|_{Y0,Y1=0} | 6 | 5 | False | False | False | None/None (None) | None | 20.1s |
| F-restriction-2 | F|_{Y0,Y2=0} | 6 | 5 | False | False | False | None/None (None) | None | 21.1s |
| F-restriction-2 | F|_{X2,Y2=0} | 6 | 5 | False | False | False | None/None (None) | None | 27.0s |
| F-restriction-2 | F|_{Y1,Y2=0} | 6 | 5 | False | False | False | None/None (None) | None | 25.1s |
| F-restriction-3 | F|_{X0,X1,X2=0} | 6 | 3 | False | False | False | None/None (None) | None | 21.9s |
| F-restriction-3 | F|_{X0,X1,Y0=0} | 6 | 3 | False | False | False | None/None (None) | None | 21.3s |
| F-restriction-3 | F|_{X0,X1,Y1=0} | 6 | 3 | False | False | False | None/None (None) | None | 24.4s |
| F-restriction-3 | F|_{X0,X1,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 22.1s |
| F-restriction-3 | F|_{X0,X2,Y0=0} | 6 | 3 | False | False | False | None/None (None) | None | 21.1s |
| F-restriction-3 | F|_{X0,X2,Y1=0} | 6 | 3 | False | False | False | None/None (None) | None | 22.7s |
| F-restriction-3 | F|_{X0,Y0,Y1=0} | 6 | 3 | False | False | False | None/None (None) | None | 22.4s |
| F-restriction-3 | F|_{X0,X2,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 26.6s |
| F-restriction-3 | F|_{X0,Y0,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 23.1s |
| F-restriction-3 | F|_{X0,Y1,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 22.8s |
| F-restriction-3 | F|_{X1,X2,Y0=0} | 6 | 3 | False | False | False | None/None (None) | None | 22.4s |
| F-restriction-3 | F|_{X1,X2,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 21.0s |
| F-restriction-3 | F|_{X1,X2,Y1=0} | 6 | 3 | False | False | False | None/None (None) | None | 24.7s |
| F-restriction-3 | F|_{X1,Y0,Y1=0} | 6 | 4 | False | False | False | None/None (None) | None | 21.7s |
| F-restriction-3 | F|_{X1,Y0,Y2=0} | 6 | 4 | False | False | False | None/None (None) | None | 21.8s |
| F-restriction-3 | F|_{X2,Y0,Y1=0} | 6 | 3 | False | False | False | None/None (None) | None | 22.5s |
| F-restriction-3 | F|_{X2,Y0,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 21.4s |
| F-restriction-3 | F|_{X1,Y1,Y2=0} | 6 | 4 | False | False | False | None/None (None) | None | 24.9s |
| F-restriction-3 | F|_{X2,Y1,Y2=0} | 6 | 3 | False | False | False | None/None (None) | None | 25.0s |
| F-restriction-3 | F|_{Y0,Y1,Y2=0} | 6 | 4 | False | False | False | None/None (None) | None | 25.7s |
| FX-pencil | X0^6+X1^6+X2^6 + 0 X0^3X2^3 | 6 | 3 | False | False | False | None/None (None) | None | 20.7s |
| FX-pencil | X0^6+X1^6+X2^6 + 1 X0^3X2^3 | 6 | 4 | False | False | False | None/None (None) | None | 23.8s |
| FX-pencil | X0^6+X1^6+X2^6 + -5/7 X0^3X2^3 | 6 | 4 | False | False | False | None/None (None) | None | 21.7s |
| FX-pencil | X0^6+X1^6+X2^6 + 10/7 X0^3X2^3 | 6 | 4 | False | False | False | None/None (None) | None | 20.1s |
| FX-pencil | X0^6+X1^6+X2^6 + 1/7 X0^3X2^3 | 6 | 4 | False | False | False | None/None (None) | None | 23.0s |
| join-pencil | F_X + -1 F_Y | 6 | 7 | False | False | False | None/None (None) | None | 21.6s |
| join-pencil | F_X + 1 F_Y | 6 | 7 | False | False | False | None/None (None) | None | 21.9s |
| join-pencil | F_X + 1/2 F_Y | 6 | 7 | False | False | False | None/None (None) | None | 22.1s |
| join-pencil | F_X + 2 F_Y | 6 | 7 | False | False | False | None/None (None) | None | 23.1s |
| pieces-8 | X2^2*H' | 8 | 12 | False | True | False | 685/686 (698) | True | 159.2s |
| pieces-8 | X2^2*C*Z' | 8 | 9 | False | True | False | 685/686 (698) | True | 162.4s |
| pieces-8 | X2^2*H | 8 | 8 | False | False | False | None/None (None) | None | 181.6s |
| pieces-8 | X2^2*C^2 | 8 | 6 | False | False | False | None/None (None) | None | 180.9s |
| pieces-8 | X2^2*Z'^2 | 8 | 10 | False | True | False | 685/686 (698) | True | 188.0s |
| D8-products | X2*D8 | 9 | None | None | True | False | 1155/1156 (None) | None | 565.6s |
| pieces-9 | X2^3*H | 9 | None | None | False | False | None/None (None) | None | 538.7s |
| pieces-9 | X2^3*H' | 9 | None | None | True | False | 1738/1739 (None) | None | 576.7s |
| pieces-9 | C*H' | 9 | None | None | False | False | None/None (None) | None | 558.8s |
| pieces-9 | Z'*H' | 9 | None | None | False | False | None/None (None) | None | 553.1s |
| pieces-9 | C^3 | 9 | None | None | False | False | None/None (None) | None | 573.3s |
| pieces-9 | Z'^3 | 9 | None | None | False | False | None/None (None) | None | 566.9s |
| pieces-9 | C^2*Z' | 9 | None | None | False | False | None/None (None) | None | 584.6s |
| pieces-9 | C*Z'^2 | 9 | None | None | False | False | None/None (None) | None | 582.0s |

### Other log entries

- selftest w2core self-test: vectorized Koszul-varpro residual vs direct E2 evaluation rel err 1.3e-16; gradient FD check 8e-10 (full) / 1.7e-8 (relax); lane F d=6 sector (0,0): S=49 J=8 kerA=8 ker_mod_J=0, 3 restarts best rel 0.19 (full), 8.2e-4 (relax, dim R_12^(0,3)=119); Fermat control F0 d=6 sector (0,0): S=49 J=8 kerA=31 ker_mod_J=23; target direction X0^6+X1^6+X2^6-F0/2 has residual exactly 0 in the f
- macaulay-prep 
- macaulay HF(k)=0 for some k means the relaxation quadrics span all degree-k forms mod p, hence over Q: no E1-admissible D of this degree/parity satisfies E2
- served 
- served 
- served 
- served-span 
- macaulay-blocks degree-3 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-3 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 5, Y-parity 0
- macaulay-blocks degree-3 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-3 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 5, Y-parity 1
- macaulay-blocks degree-3 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-3 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 0
- macaulay-blocks degree-3 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-3 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 1
- macaulay-blocks degree-4 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-4 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 0
- macaulay-blocks degree-4 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-4 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 1
- macaulay-blocks degree-3 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-3 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 0
- macaulay-blocks degree-4 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-4 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 0
- note remote capacity: load 700-770 on 128 cores, 130-140 jobs from other agents; the runner refuses submissions above load 150; jobs jac (rigorous d=7 exclusion + Fermat calibration + Jacobian-type searches), sec8, sec9, served are queued with --wait; local work restricted to single niced processes (rigorous relaxation exclusions at d=5,6 done locally).
- macaulay-blocks degree-3 Macaulay matrix of the relaxation quadrics, block-diagonal by G-character of w; HF_k = 0 means all degree-3 forms in w lie in the ideal mod p, hence over Q: no certificate with deg D = 6, Y-parity 1
- relaxation-point-d6 d=6 Y-even: relaxation variety modulo the F-gauge line is at most one projective point (HF_4 = 1, Macaulay bound); E_X(F) = sum X_i dF/dX_i (gauge-equivalent to E_Y(F)) is NOT that point: Bil(E_X F, A X) = A * (five mixed couplings) is not in J_12 (exactly over Q). Values recorded from the printed output of the local check (the logging call was cut by the 110 s local timeout).
