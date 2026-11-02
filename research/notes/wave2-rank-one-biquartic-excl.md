# Wave 2 — `rank-one-biquartic` — lens EXCLUSION / RIGOR

STATUS: no witness. PROVED over Q: no certificate with deg D <= 5 (all parities; d<=4 has no E1-admissible D at all, d=5 by exact weight-2 certificates). Degree 6: excluded mod 10^9+7 at the Hessian-relaxation level (V(forms) = spurious line D ~ F); exact char-0 version: exact step still running (local queue; outputs: hpipe_d6.out, hpipe_d6.*.json, hpipe_d6.D.out). Degree 7+: local queue running (k7exact 7 -> hpipe_d7 -> Singular (J:A) p=65521 -> verify_e1kernel p=10^9+7); outputs land in wave2-excl/*.out and *.log.jsonl (NOT yet merged into the main log if the session ended before they finished).

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
| `D in J` / `A D in J` membership table: `L, M, f, g, L^2M, f-g, A L, A M, A L^2, A M^2, A L M` not in `J`; `E_X(F), E_Y(F), f-2g, F` in `J`; degree-12 products `fg, fL^2M, gL^2M, (L^2M)^2, f^2, g^2` not in `J` (dossier 4.2/4.6) | p = 10007 (Singular) | `hodge.jacobian.in_jacobian` at p = 10^9+7 (J blocks at regular-sequence rank, so non-membership is rigorous over Q) | AGREE |
| `(J:A)` minimal generators 3 (deg 7), 3 (deg 8), 456 (deg 10); HF of `S/(J:A)`; `Disc in J` (dossier 4.6) | Singular p = 10007 | `quotJA_p65521.sing` | pending (local queue) |
| fixed-D inconsistency of `E_X(F), E_Y(F), X0 F_0, F_Xi, ...` (wave-1 log, engine v2 compressed solver) | p = 2^31-1 | `verify_fixedD.py` dense rref at p = 2147483423 | pending (remote saturated) |
| `K_7 subset S_7(X) + M S_5(X)` (wave-1 kstruct3, 5 primes) | mod p | `k7exact.py 7`: exact Q-computation of `W = {b in B : A b in J}` (dim 21 = 18 + 3 expected) | pending (local queue) |

## Method (precise statement of what the exact computations prove)
Let `S = Q[X0..Y2]`, `J = (F_0..F_5)` (regular sequence since `X` is smooth), `q = 1`, `A = X0^3 Y0^3`.
Suppose `(D, N, Q)` is a certificate with `d = deg D`, `D` coprime to `F` (any Y-parity, mixed allowed).
1. (Euler gauge) Reducing the identity mod `F` gives `D (A D - N(F)) == 0 (mod F)`, so `N(F) = A D + F h`, `deg h = d`.
   The field `N' = N - (h/6) E` (`E` = Euler field) satisfies `N'(F) = A D` and the cleared identity is unchanged
   (`div(phi E) = (deg phi + 6) phi`, `(phi E)(D) = d phi D`, so the `F`-terms cancel). Hence WLOG `N(F) = A D` exactly,
   and in particular `A D in J`: `D in V_d := {D in S_d : A D in J_{d+6}}`.
2. (Koszul) If `(D_k, N_k)` is a Q-basis of `V_d` with `N_k(F) = A D_k`, then `N = sum a_k N_k + K(c)` with
   `K(c)_i = sum_j c_ij F_j`, `c_ij = -c_ji in S_{d-4}` (syzygies of a regular sequence are Koszul), and the identity
   becomes `E2(a, c, Q): D div N - N(D) - F Q = 0` in `S_{2d}`: a system of polynomial equations of weight 2
   (`a:1, c:1, Q:2`) in the unknowns; let `I` be the ideal they generate in `Q[a, c, Q]`.
3. (Relaxation) `div K(c) = sum_ij (d_i c_ij) F_j in J` (antisymmetry of `c`, symmetry of the Hessian) and `K(c)(D) in J`,
   so `B_D(K(c)) in J_{2d}` and `F Q in J_{2d}`. Hence for every linear functional `lambda` on `S_{2d}` vanishing on
   `J_{2d}`, `phi_lambda(a) := lambda(B_{D(a)}(sum a_k N_k))` is a quadratic form in `a` lying in `I`
   (`phi_lambda = lambda(E2)` because `lambda` kills the `c`- and `Q`-terms). For polar directions `D_k = dF(G_k)`,
   `N_k = A G_k` the exact identity `B_{D_k}(N_l) = D_k (G_l(A) + A div G_l) - A sum_ij G_{l,i} (d_i G_{k,j}) F_j
   - A Hess_F(G_l, G_k)` shows that modulo `J_{2d}` the quadratic column of the pair `(k,l)` is `-A [Hess_F(G_k,G_l) +
   Hess_F(G_l,G_k)]` (used as the sparse column representative; a random-pair check mod p is in the log).
4. (Certificates) A vector `lambda in Q^{S_{2d}}` with `lambda(J_{2d}) = 0` and `lambda(E2-columns) = u^2` (for a
   coordinate `u` of `V_d / (F S_{d-6})`) proves `u^2 in I`, i.e. `u = 0` on `V(I)`. Such `lambda` are found from the
   pivot structure mod p (an invertible square submatrix mod p is invertible over Q), solved exactly with `fmpq_mat`,
   and verified exactly against every column. Quadratic forms `phi in I` are likewise certified exactly (CRT lift of the
   mod-p rref rows, then the same exact certificate). Finally Singular over Q (char 0) shows that the ideal generated by
   the certified `u`'s and `phi`'s contains every coordinate of `V_d/(F S_{d-6})` in its radical: `V(I)` is contained in
   the spurious locus `{D in F S_{d-6}}`, contradicting coprimality. Hence no certificate of degree `d` exists.
5. (E1-kernel over Q) With `p` not dividing 7: `rank_p J_{d+6} = rank_Q J_{d+6}` whenever the modular rank equals the
   regular-sequence value (checked in every run), and then `dim_Q V_d <= dim_p V_d`; `J_d subset V_d` always. So
   `V_d = J_d` over Q for `d <= 6` (modular kernel = `dim J_d`). Non-Jacobian generators (`d = 7, 8`) are computed
   exactly over Q inside the base-type subspace (`k7exact.py`), with exact E1 lifts; completeness of the basis uses the
   modular upper bound `dim_Q V_d <= dim_p V_d` again.

### Last state of the degree-6 pipeline (hpipe_d6.out, 20:21:52)
```
[A blk0 p=1000000007] rank 2373 rank_J 2223 forms 150 (26s)
[A blk0 p=999999937] rank 2373 rank_J 2223 forms 150 (48s)
[A blk0 p=999999929] rank 2373 rank_J 2223 forms 150 (68s)
[A blk1 p=1000000007] rank 2368 rank_J 2214 forms 154 (90s)
[A blk1 p=999999937] rank 2368 rank_J 2214 forms 154 (110s)
[A blk1 p=999999929] rank 2368 rank_J 2214 forms 154 (132s)
[A] coordinate squares in row space mod p: 25 of 35; failing: ['a18', 'a19', 'a26', 'a31', 'a32', 'a5-a10', 'a5-a15', 'a5-a20', 'a5-a25', 'a5-a30']
[B] residual vars 11 (10 failing + 1 centers); Singular mod p: dim V(I): 1  #forms: 39  time: 0; not in radical: 0 (132s)
[B] blk0: restricted span rank 0, selected 0 full rows
[B] blk1: restricted span rank 39, selected 39 full rows
[C blk0] lifted 0/0 selected forms with 3 primes; certifying 0 targets (132s)
[C blk1] lifted 36/39 selected forms with 3 primes; certifying 61 targets (132s)
```

Local queue log (`wave2-excl/localq.log`):
```
[20:15:31] START timeout 1800 python3 hess_pipeline.py 6 3 --sparse --tag hpipe_d6
```

## How to resume (for the next agent)
* The local queue `wave2-excl/localq.sh` runs: `hess_pipeline.py 6 3 --sparse` -> `k7exact.py 7` -> `hess_pipeline.py 7 3 --sparse --nonj nonj_d7.json`
  -> `Singular -q quotJA_p65521.sing` -> `verify_e1kernel.py 1000000007`.  Check `localq.log`, then merge `wave2-excl/*.log.jsonl` into
  `wave2-excl.log.jsonl` and update this report.  Remote jobs to submit when the 128-core box is free: `wave2-excl/remote_jobs.sh`
  (degree-8 exact generators + pipeline, degree-9 modular screening, dense fixed-D re-verification, Singular quotient at 32003).
* Three remote submissions (`excl-rank-one-biquartic-fam6`, `-K7lift`, `-V1e1kernel`) were queued with `--wait 7200` while the remote
  was saturated (load > 1000); if they never started, resubmit (the local queue covers V1 and the K7 lift exactly).
