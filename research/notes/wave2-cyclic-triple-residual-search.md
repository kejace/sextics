# Wave 2 search report: `cyclic-triple-residual` (q = 1, A = X0^3 Y0^3)

STATUS (2026-11-01 19:30 local): **no witness found so far**. Remote job `search-cyclic-triple-residual-aux` running (E1-admissible spaces d=4..10 done with exact bases; structured sweep in progress). Jobs `jac` (rigorous mod-J exclusion for d=5,6,7 via msolve + Fermat calibration + Jacobian-type joint searches), `sec8`, `sec9` (sector joint searches at d=8, 9) are queued behind the remote load cap (load ~300 on 128 cores from other agents; the runner now refuses submissions above load 150 and `--wait` polls for capacity).

Code: `/tmp/hodge/results/cyclic-triple-residual/w2core.py` (exact E1-admissible spaces per G-sector via fmpz nullspace; Koszul-reduced, sector-restricted variable-projection search; mod-J relaxation; reconstruction + exact confirmation through the wave-1 sector solver `sect.py`), `w2run.py` (multiprocessing driver), job scripts `w2_e1space.py`, `w2_jac.py`, `w2_sector.py`, `w2_sweep.py`, `w2_macaulay.py` + `w2_hilb.py` (rigorous relaxation exclusion), `w2_report.py` (renders the tables below from the log). Canonical log: `/tmp/hodge/results/cyclic-triple-residual/wave2-search.log.jsonl` (merged from per-process `w2log-*.jsonl` files by `w2_pull.sh`).

## Method
Split form E1: q N(F) = A D, E2: D div N - N(D) = F Q. Every E1 lift is N0 + K(c), K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j), deg c_ij = d - 4 (book ch. 5.2).
For a linear family D(a) = sum a_k D_k with exact lifts N0_k the E2 residual r = Bil(D(a), N0(a) + K(c)) - F Q is linear in (c, Q):
f(a) = min_{c,Q} |r|^2 on the unit sphere in a (variable projection with the analytic Golub-Pereyra gradient; the subspace of F-multiples is projected out; L-BFGS-B random restarts; rel = f/|Bil(D, N0)|^2).
Mod-J ("Hessian") relaxation: f_J(a) = |P_perp^T Bil(D(a), N0(a))|^2 with P_perp an orthonormal basis of the complement of J(F)_{2d} in the E2 row space -- a necessary condition (im T_D and F S lie in J).
Character sectors of the order-18 diagonal group G (generators mu3 = (0,0,2,0,0,0), s = (1,5,1,2,2,0); labels (m.mu3, m.s) mod 6; chi_A = (0,3); chi(N_i) = chi_A + c + chi(e_i), chi(Q) = chi_A + 2c).
Rigorous exclusion at Jacobian-type degrees: for d <= 7 every E1-admissible D is w.grad F; the relaxation is a system of dim R_{2d}^- quadrics in the w-space; if msolve's leading ideal mod p shows the Hilbert function vanishing in some degree, the Macaulay matrix has full column rank mod p, hence over Q, so w = 0 is the only solution (no certificate of that degree/parity). Verified on d=5 (Y-even): the 6 quadratic monomials are all leading terms.

Jacobian ring dimensions (Hilbert series (1+t+..+t^4)^6, Y-odd parts): R_10^- = 744, R_12^- = 866, R_14^- = 744, R_16^- = 468, R_18^- = 208, R_20^- = 60, R_22^- = 9, R_24^- = 0.

## Self-test (local)
Vectorized residual vs direct evaluation: 1.3e-16. Gradient FD check: 8e-10 / 1.7e-8. Fermat control (F0 = Fermat sextic, d=6, sector (0,0), family dim 31): the known witness direction (X0^6+X1^6+X2^6, gauge-shifted by F0/2) has residual 0 in the family; 6 random restarts missed it (best rel 1e-2) -> hit-rate calibration with 400 restarts is part of the jac job.

## Results so far
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

### Other log entries

- selftest w2core self-test: vectorized Koszul-varpro residual vs direct E2 evaluation rel err 1.3e-16; gradient FD check 8e-10 (full) / 1.7e-8 (relax); lane F d=6 sector (0,0): S=49 J=8 kerA=8 ker_mod_J=0, 3 restarts best rel 0.19 (full), 8.2e-4 (relax, dim R_12^(0,3)=119); Fermat control F0 d=6 sector (0,0): S=49 J=8 kerA=31 ker_mod_J=23; target direction X0^6+X1^6+X2^6-F0/2 has residual exactly 0 in the f
- macaulay-prep 
- macaulay HF(k)=0 for some k means the relaxation quadrics span all degree-k forms mod p, hence over Q: no E1-admissible D of this degree/parity satisfies E2
