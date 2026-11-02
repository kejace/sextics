# Wave 2 search report: klein-coupled-character (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (search in progress; this file is rewritten after every batch).

Log: `/tmp/hodge/results/klein-coupled-character/wave2-search.log.jsonl` (one JSON line per batch; engine tag
`w2-koszul-varpro`).  Code: `w2lib.py` (sector-restricted Koszul variable projection, sector E1 spaces with rational
lifts, Hessian relaxation), drivers `w2_kernels.py`, `w2_hess.py`, `w2_varpro.py`, `w2_struct.py`, control test
`w2_control.py`.  Remote jobs are named `w2-*`.

## Method (what is new relative to wave 1)
* Joint search over D in a *linear family inside one sigma-sector* (Y-parity, weight k mod 7), in the Koszul
  correction-map formulation: D(a) = sum a_k D_k with exact E1 lifts N_k, unknowns p_ij (degree d-4) and Q
  (degree 2d-6), both character-restricted; Q is eliminated by the t-monic remainder modulo F (rows = monomials of
  degree 2d not divisible by Y0^6).  Systems are 6-25x smaller than the engine's unrestricted VarProFamily.
* Optimizers: L-BFGS on the sphere (engine style), Levenberg-Marquardt/TRF with the Kaufman variable-projection
  Jacobian, joint LM over (a, x).  Candidates with rel < 1e-10 are reconstructed in the canonical rem_F chart and
  re-solved exactly (kcc2 Koszul test mod p, 12-prime CRT lift, core.verify_witness).
* Hessian (mod-J) relaxation per sector: |proj_(R/J) B_D(N0)|^2 minimized numerically (real), and the exact quadric
  system mod 32003 in a rational chart complementary to the F-multiples fed to Singular (dim(std(I)) = 0 means
  no sigma-eigen denominator of that degree/parity/weight exists, over Q by semicontinuity).

## Results so far
(filled in as batches finish)
