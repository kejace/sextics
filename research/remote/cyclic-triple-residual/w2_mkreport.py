"""Compose the wave-2 report: STATUS header (from w2_status.txt) + method + generated tables."""
import subprocess, time, os
WD = '/tmp/hodge/results/cyclic-triple-residual'
status = open(f'{WD}/w2_status.txt').read().strip() if os.path.exists(f'{WD}/w2_status.txt') else 'no witness found so far.'
rep = subprocess.run(['python3', f'{WD}/w2_report.py'], capture_output=True, text=True, cwd=WD).stdout
head = f"""# Wave 2 search report: `cyclic-triple-residual` (q = 1, A = X0^3 Y0^3)

STATUS ({time.strftime('%Y-%m-%d %H:%M')} local): {status}

Code (all in `{WD}/`): `w2core.py` (exact E1-admissible spaces per G-sector via fmpz nullspace; Koszul-reduced, sector-restricted variable-projection search; mod-J relaxation; reconstruction + exact confirmation through the wave-1 sector solver `sect.py`), `w2run.py` (multiprocessing driver), jobs `w2_e1space.py`, `w2_jac.py` (+ `w2_jacjob.sh`), `w2_sector.py`, `w2_sweep.py`, `w2_served.py` (served target spaces / sum principle), `w2_macblocks.py` (rigorous relaxation exclusion by G-graded block Macaulay ranks), `w2_macaulay.py` + `w2_hilb.py` (msolve cross-check), `w2_report.py`/`w2_mkreport.py` (this report). Canonical log: `{WD}/wave2-search.log.jsonl` (merged from per-process `w2log-*.jsonl` files by `w2_pull.sh`; remote mirror `/tmp/hodge/remote/cyclic-triple-residual/`).

## Method
Split form E1: q N(F) = A D, E2: D div N - N(D) = F Q. Every E1 lift is N0 + K(c), K(c) = sum_{{i<j}} c_ij (F_j e_i - F_i e_j), deg c_ij = d - 4 (book ch. 5.2).
For a linear family D(a) = sum a_k D_k with exact lifts N0_k the E2 residual r = Bil(D(a), N0(a) + K(c)) - F Q is linear in (c, Q):
f(a) = min_{{c,Q}} |r|^2 on the unit sphere in a (variable projection, analytic Golub-Pereyra gradient, L-BFGS-B random restarts; the subspace of F-multiples is projected out, which is lossless by the gauge (D, N) -> (D + lF, N + A l Z/6); rel = f/|Bil(D, N0)|^2).
Mod-J ("Hessian") relaxation: f_J(a) = |P_perp^T Bil(D(a), N0(a))|^2 with P_perp an orthonormal basis of the complement of J(F)_{{2d}} in the E2 row space -- a necessary condition (im T_D and F S lie in J).
Character sectors of the order-18 diagonal group G (generators mu3 = (0,0,2,0,0,0), s = (1,5,1,2,2,0); labels (m.mu3, m.s) mod 6; chi_A = (0,3); chi(N_i) = chi_A + c + chi(e_i), chi(Q) = chi_A + 2c).
Rigorous exclusion at Jacobian-type degrees (d <= 7, where Ann([A])_d = 0 so every E1-admissible D is w.grad F): the relaxation is a system of dim R_{{2d}}^- quadrics q_s(w) in the w-space; the degree-k Macaulay matrix is block-diagonal by the G-character of w; if every block has full column rank mod p (p good for J_{{2d}}: the mod-p rank of J_{{2d}}^- equals the Hilbert-series value), the ideal contains all degree-k forms over F_p, hence over Q, so w = 0 is the only solution. The spurious locus L = {{l Z/6}} (D = lF) consists of solutions; V(q) is invariant under translation by L (gauge), so restricting to a complement of L (dropping the coordinates (0, m) with X0 | m) is lossless. Macaulay's bound HF(k+1) <= HF(k)^<k> turns HF_k = 1 into "at most one projective point".

Jacobian ring dimensions (Hilbert series (1+t+..+t^4)^6, Y-odd parts): R_10^- = 744, R_12^- = 866, R_14^- = 744, R_16^- = 468, R_18^- = 208, R_20^- = 60, R_22^- = 9, R_24^- = 0.

## Self-test (local)
Vectorized residual vs direct evaluation: 1.3e-16. Gradient FD check: 8e-10 / 1.7e-8. Fermat control (F0 = Fermat sextic, d=6, sector (0,0), family dim 31): the known witness direction (X0^6+X1^6+X2^6, gauge-shifted by F0/2) has residual 0 in the family; 6 random restarts missed it (best rel 1e-2) -> hit-rate calibration with 400 restarts is the first batch of the jac job.

## Results
"""
rig = open(f'{WD}/w2_rigor.md').read() if os.path.exists(f'{WD}/w2_rigor.md') else ''
open('/tmp/hodge/notes/wave2-cyclic-triple-residual-search.md', 'w').write(head + rig + '\n' + rep)
print('report written:', len(head + rep), 'chars')
