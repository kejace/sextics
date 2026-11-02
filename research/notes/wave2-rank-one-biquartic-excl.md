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
