# Wave 2 — `rank-one-biquartic` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; re-verification and rigorous degree bound in progress)

Log: `/tmp/hodge/results/rank-one-biquartic/wave2-excl.log.jsonl` (one line per experiment).
Packet: `/tmp/hodge/results/rank-one-biquartic/wave2-excl/`.

## Results so far

### R1. Exact identity re-checks (local, Fractions)
* `F = f + g + (5/7) L^2 M` re-verified; `B_{E_X(F)}(A X) = 8 t A L^2 M` (dossier 1.6/4.5) re-verified exactly;
  `B_F(A eps/6) = A F` (the spurious `D = F` solution of the normalized system) verified; Koszul fields tangent. AGREE.

### R2. Degree 5 is excluded over C (PROVED, all parities including mixed-parity D)
Method (`wave2-excl/w2lib.py`, `fam_run.py 5 1000000007 --cert`): Euler gauge + Koszul parametrization; the E1-kernel
`V_5 = {D in S_5 : A D in J}` equals `J_5 = <F_0,...,F_5>` (both parities; wave-1 modular ranks with the J-block at the
regular-sequence rank, re-verification job V1 pending). Family `D = sum a_k F_k`, `N = sum a_k A e_k + K(c)` (90 Koszul
unknowns, linear forms), `Q in S_4` (126). The weight-2 Macaulay matrix of E2 (2924 rows x 687 cols) has rank 507 mod
1000000007, and ALL 21 monomials `a_k a_l` lie in its row space; exact rational certificates `lambda^T E = a_k^2`
(k = 0..5; supports 51/248 rows, heights <= 197 digits) were verified over Q. Hence `a = 0` on the solution variety:
**no certificate with deg D = 5 exists** (no spurious locus in degree 5). Receipts: `wave2-excl/fam_d5_p1000000007_all.{json,certs.json,log.jsonl}`.

### R3. Degree 6: Hessian relaxation (in progress; modular part done)
Family `V_6 = J_6` (36 directions, both parities), spurious line `D = F` (`G = Euler/6`). Observation used throughout:
every Koszul correction `B_D(K(c))` and every `F Q` lies in `J_{2d}` (`div K(c) = sum (d_i c_ij) F_j` by antisymmetry
of `c` and symmetry of the Hessian), so any functional on `(S/J)_{2d}` applied to `B_D(N_0(D))` gives a quadratic form in
the family coordinates `a` that must vanish for a certificate ("Hessian relaxation"); for a polar `D = dF(G)` the class is
`[B_D(A G)] = -[A Hess_F(G,G)] = -[A G(D)]` in `(S/J)_{2d}` (exact identity, re-verified mod p on random pairs).
* `hess_run.py 6 1000000007 --blocks 1 --cert`: block of Y-odd equations: rank 2368 = 2214 (`J_12^-`) + 154; 25 of the 35
  coordinate squares `u^2` of `V_6/<F>` lie in the weight-2 ideal mod p and all 25 were certified exactly over Q
  (`hess_d6_p1000000007_all.blk1.certs.json`). The 10 failing coordinates: the five `Y_j F_{3+k}` directions
  `(3,Y2),(3,Y1),(4,Y0),(5,Y1),(5,Y0)` and the five differences inside the spurious group.
* `quadforms_dump.py 6` + Singular mod p: the 304 quadratic forms (154 same-parity + 150 cross-parity) generate an ideal whose
  zero set is exactly the spurious line: `dim V = 1` (affine) and every coordinate `u` of `V_6/<F>` satisfies
  `dim(std(I + (u - 1))) = -1`. So **degree 6 is excluded modulo p** at the Hessian-relaxation level; the characteristic-zero
  version (exact forms via CRT lift + exact certificates, then Singular over Q) is running (`hesslift_d6`).

## Re-verification ledger (goal (a)) — updated as jobs finish
| claim (source) | wave-1 / dossier value | wave-2 check | status |
|---|---|---|---|
| `F = f + g + (5/7) L^2 M`, `M = L(Y)` (dossier 1.1) | exact | exact Fractions (R1) | AGREE |
| `B_{E_X(F)}(A X) = 8 t A L^2 M` (dossier 1.6, 4.5) | exact | exact Fractions (R1) | AGREE |
| E1-kernel `V_d = J_d` for `d <= 6` both parities; `dim K_7 = 3 (even), 0 (odd)`; `K_8 = 9+9`; `K_9 = 27+27` (wave-1 log, dossier 4.2) | p = 2^31-1, 10007, 65521 | `verify_e1kernel.py` at p = 10^9+7 (`hodge.jacobian.mult_kernel`, different code path) + J-block rank check | pending (local queue) |
| `(J:A)` minimal generators 3 (deg 7), 3 (deg 8), 456 (deg 10); HF of `S/(J:A)`; `Disc in J` (dossier 4.6) | Singular p = 10007 | `quotJA_p65521.sing` | pending (local queue) |
| fixed-D inconsistency of `E_X(F), E_Y(F), X0 F_0, F_Xi, ...` (wave-1 log, engine v2 compressed solver) | p = 2^31-1 | `verify_fixedD.py` dense rref at p = 2147483423 | pending (remote saturated) |
| `K_7 subset S_7(X) + M S_5(X)` (wave-1 kstruct3, 5 primes) | mod p | `k7exact.py 7`: exact Q-computation of `W = {b in B : A b in J}` (dim 21 = 18 + 3 expected) | pending (local queue) |
