# Wave 2 search report: `cyclic-triple-residual` (q = 1, A = X0^3 Y0^3)

STATUS: no witness found so far. Four remote jobs submitted at 19:35 local (aux: E1 spaces d=4..10 + structured sweep; jac: Fermat calibration + Jacobian-type d=5,6,7; sector joint searches at d=8 and d=9). Results pending.

Code: `/tmp/hodge/results/cyclic-triple-residual/w2core.py` (exact E1-admissible spaces per G-sector via fmpz nullspace; Koszul-reduced, sector-restricted variable-projection search; mod-J relaxation; reconstruction + exact confirmation through the wave-1 sector solver `sect.py`), `w2run.py` (multiprocessing driver), job scripts `w2_e1space.py`, `w2_jac.py`, `w2_sector.py`, `w2_sweep.py`. Log: `/tmp/hodge/results/cyclic-triple-residual/wave2-search.log.jsonl`.

## Method
Split form E1: q N(F) = A D, E2: D div N - N(D) = F Q. Every E1 lift is N0 + K(c), K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j), deg c_ij = d - 4.
For a linear family D(a) = sum a_k D_k with exact lifts N0_k the E2 residual r = Bil(D(a), N0(a) + K(c)) - F Q is linear in (c, Q):
f(a) = min_{c,Q} |r|^2 on the unit sphere (variable projection, analytic gradient), F-multiples projected out.
Mod-J relaxation: f_J(a) = |P_perp^T Bil(D(a), N0(a))|^2 (orthogonal complement of J(F)_{2d} in the E2 row space) -- a necessary condition.
Character sectors of the order-18 diagonal group G (generators mu3 = (0,0,2,0,0,0), s = (1,5,1,2,2,0); labels (m.mu3, m.s) mod 6; chi_A = (0,3)).

## Self-test (local)
Vectorized residual vs direct evaluation: 1.3e-16. Gradient FD check: 8e-10 / 1.7e-8. Fermat control (F0, d=6, sector (0,0), family dim 31): the known witness direction has residual 0; 6 random restarts missed it (best rel 1e-2) -> hit-rate calibration with 400 restarts running on the remote.

## Results
(pending)
