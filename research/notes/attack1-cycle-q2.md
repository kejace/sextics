# Attack 1 — `cycle-q2` (resumed session 2026-11-01)

**No witness found.** (Status line updated at the end of the session; see §3.)

Instance: `F = Σ_i (X_i^6 + Y_i^6) + ½(X0^3 X1 Y0 Y1 + X1^3 X2 Y1 Y2 + X0 X2^3 Y0 Y2)`, `q = 2`,
`A = X0^2 X1^2 X2^3 Y0^2 Y1^2 Y2` (degree 12). Same hypersurface as the accepted `cycle-three`
(`A = X0^3 Y0^3`, `q = 1`, nonconstructive Lean proof: no explicit `D`).

Work directory: `/tmp/hodge/results/cycle-q2/` (`log.jsonl` = every tested denominator; `eq.py` =
character-sector fixed-D solver; `modlu.py` = blocked float64 modular LU; `exp_*.py` = experiments).

## 0. Attack plan (written before running)

1. Re-validate the fast modular path (`modlu`) against flint on an identical matrix.
2. Exact rational E1 kernel `K_d = {D : A·D ∈ J(F)}` per character sector, reduced mod `J_d`, for
   `d = 6` (26 nontrivial dims mod p from the previous run), then full fixed-D tests of every
   complement representative and small pencils, with `D` and `D²`.
3. Cone families suggested by the (3,1)-support constraint (book ch.3 §3.1): sextic Fermat cones in 4
   variables (vertex a line `X ∩ L` = 6 points, exceptional surface = Fermat sextic surface, p_g = 10),
   cubic/quartic/quintic/sextic Fermat cones in 3 variables (vertex a plane, exceptional surface
   `C × B` with `C` the Fermat sextic curve, p_g = 10·g(B)).
4. Structural pieces of `F` (edges, couplings, F minus couplings) at `D²` (the previous run covered
   only 5 of them), Jacobian family (`F_i`, `x_j F_i`, radial combos, `F_i F_j`).
5. Mathematics: where the (3,1) part can come from in the (3,3)-surface fibration picture; what a
   certificate must look like.
