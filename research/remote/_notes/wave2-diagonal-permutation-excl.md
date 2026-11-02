# Wave 2 — diagonal-permutation — EXCLUSION / RIGOR lens

STATUS: no witness (exclusion lens) -- in progress.  Best rigorous bound so far: deg D >= 7 for EVERY homogeneous D (any Y-parity, mixed included); wave 1 had only d<=4 rigorous + d<=7 for Jacobian-type D tested individually.

Started 2026-11-01 (wave 2). Packet: /tmp/hodge/results/diagonal-permutation/wave2-excl/ ; log: /tmp/hodge/results/diagonal-permutation/wave2-excl.log.jsonl

## Plan
(a) re-verify wave-1 rigorous claims at new primes / new method (E1 kernels d<=8, D8, served(partials)=0).
(b) rigorous lower bound on deg D via E1 + "E2 mod J" relaxation (Lemma below) + msolve Groebner runs, per sector
    (Y-parity, mixed parity, with/without the D8 line).
(c) non-Jacobian E1-kernel elements (D8; d=9 kernel) and sector kernels d=10..12; structured families; full E1+E2 tests.
(d) packet + README; (e) draft post.

## Lemma (E2 modulo J)  [to be proven in README]
If (D,N,Q) is a certificate at q=1 then N(F)=A D + F H, N = u + K + (H/6)E with u(F)=A D, K Koszul; and
Psi(D) := D div u - u(D)  lies in J(F) (independent of the choice of u modulo J).  For D = w(F): Psi = -A * (w^T Hess(F) w) mod J.

## Results so far (2026-11-01, wave 2)

### R1. Rigorous lower bound deg D >= 7 (all parities)
* d <= 4: E1 kernel = 0 mod P1=1000000007 (both parities; recomputed locally, e1k job repeats at P2) and J_d = 0  => no D.
* d = 5: all admissible D are w(F), w constant; relaxation Psi = -A w^T Hess(F) w in J_10 = requires w^T Hess(F) w = 0 as a quartic
  (Ann(A)_4 = 0); the 21 relaxation equations in the 6 unknowns span ALL quadratic monomials (receipt relax_d5_mixed.ms),
  so w = 0.  (Hand proof: d^2_w F = 0 forces [w] to be a point of multiplicity >= 5 on X, impossible for smooth X.)
* d = 6: all admissible D are w(F), w linear (E1 kernel 0 in both parities).  Relaxation systems (quadratic, Euler family removed):
  even sector 17 unknowns / 134 indep. eqs, odd 18 / 153, MIXED 35 / 404; msolve reduced Groebner bases are zero-dimensional
  (pure power of every variable among leading terms; receipts relax_d6{e,o,m}.{ms,gb,json}, p = 1000000007, good reduction in
  degree 12 checked).  Hence V = {0} over Q as well (Hilbert-function argument) and the only degree-6 w(F) passing E2 mod J
  are multiples of F.  => no certificate with deg D = 6.
  Cross-check: Singular 4.3.2 (std, option(redSB)) on the same three .ms systems gives reduced bases of the same sizes
  (139 / 153 / 437) with pure powers of all 17 / 18 / 35 variables among the leading terms (singular_d6.out).
* Note on wave 1: its status line "no certificate with deg D <= 7" was not proven there (only individual vector fields at d=5,6
  were tested); the gap is closed here for d <= 6 and d = 7 is running (relax7F job).

### R1b. Re-verification of the wave-1 E1-kernel table at the fresh prime P1 = 1000000007 (job e1kF, in progress)
* d = 1..7, both parities: kernel 0 at P1 = 1000000007 AND at P2 = 998244353 (agrees with wave 1 at 2147483647); good
  reduction in degrees d and d+6 confirmed at both primes (dim R_n mod p = Koszul value).  d = 8: kernel 1 (even), 0 (odd)
  at both P1 and P2 (agrees); d = 9 even at P1: kernel 21 (agrees with the sibling value at 2147483647).
  A canonical representative of the d=8 line mod P1 is in e1k_e1kF.json.

### R1c. Representation theory re-verified (w2_reptheory.py, independent code, P1)
S3-characters on Jacobian-ring pieces (trace of the transposition / 3-cycle on standard-monomial bases mod P1):
R_6^- = 208 = 39 triv + 31 sign + 69 x std  (dossier table confirmed); (01)-eigenspaces 108 / 100.
R_12^- = 866 = 154 + 136 + 288 x std; R_6^+ = 218 = 46 + 28 + 72 x std (confirmed).  New: R_8^+ = 94/67/161,
R_8^- = 84/72/156, R_9^+ = 115/94/207, R_10^+ = 144/110/254, R_10^- = 132/116/248; (01)-dims in reptheory.json.
Signed-permutation symmetries of F: 12 projective elements (S3 x <tau>) [re-verified]; diagonal group order as in the engine.

### Claims checklist (wave-1 report + dossier)  [status]
| claim | status |
|---|---|
| F smooth; J(F) has Koszul Hilbert function | re-verified: dim R_n mod P1 = Koszul value in every degree used (receipts) |
| symmetry group S3 x <tau>, no reflections, diagonal group {1,tau} | re-verified (signed permutations enumerated; engine diagonal_group) |
| R_6^- = 39 triv + 31 sign + 69 std | re-verified at P1 |
| E1 kernel 0 for d <= 7 (both parities) at 2147483647 | re-verified at P1 (and P2 pending): agrees |
| E1 kernel at d = 8: 1 (even), 0 (odd) | re-verified at P1 (even); odd pending |
| D8: A*D8 in J over Q | pending exact (d8F) |
| served(F_Xi) etc. = zero class (served_dim 18) | re-verified by a different method: Ann_6([H(w,w)]) = 0 for w = d_i, Euler_X, pair-radial, ... (relaxation: served(w(F))/J_6 inside Ann_6([H(w,w)]) ); also superseded by the rigorous d=5 result |
| "no certificate with deg D <= 7" (wave-1 status line) | NOT established by wave 1 (only individual D tested); now proven for d <= 6, d = 7 pending |
| duality k_10^+ = 294, k_10^- = 262, k_11 = 510, k_12 = 677/648 | follows from k_8 and k_d (d<=7) by Gorenstein duality (argument checked); direct sector computation pending |

### R1d. Sum-principle catalogue (w2_ann.py, local, P1; in progress)
For every structured non-Jacobian factor tested so far (quadrics, cubics, quartics, Vandermonde squares, ...: 36 so far),
Ann_6(D) = {A' in R_6^- : A' D in J} = 0 (rank 208 = dim R_6^-).  Since served(D)/J_6 is inside Ann_6(D) and
served(D1 D2) = served(D1) + served(D2), every product of catalogued factors serves only the zero class, in every degree.

### R1e. Equivariant Hilbert series and the structure of the degree-9 kernel (exact, by hand + script)
Graded characters chi_g(t) = det(1 - t^5 g)/det(1 - t g) for g in {1, tau, (01), (01)tau} give dim R_n^{parity,(01)-sign}
(table in the log / w2-excl notes).  Multiplication by A (Y-odd, (01)-anti-invariant) maps the sector (par, sigma) to
(1-par, -sigma).  Forced kernels (source - target when positive):
  d <= 8: none in any sector;  d = 9: (even,+): 322 -> 301 forces >= 21; (even,-): 301 -> 322 forces 0; same for odd.
  d = 10: (even,+) >= 170, (even,-) >= 124 [sum 294 = k_10^+], (odd,+) >= 152, (odd,-) >= 109 [sum 261 = k_10^- - 1].
  d = 11: >= 275 / 235 per parity; d = 12: (even) >= 361 + 316, (odd) >= 342 + 306.
Consequence: the 21-dimensional degree-9 kernel (both parities) is exactly the forced one and lies entirely in the
(01)-INVARIANT sector provided the anti-invariant sector kernel is 0 (e1kS run: d=9 (odd,+): kernel 21 with rank 301 = full rank onto the target, i.e. exactly the forced kernel; other sectors pending).  Every non-Jacobian
E1-admissible denominator of degree 9 is then (01)-invariant modulo J, so the (01)-equivariant reduced E2' system applies.

### Job log (wave 2, this lens)
* e1kF: full-parity kernels done (d<=9 at P1, d<=8 at P2); crashed in the (01)-sector part (bug: products with representative
  monomials instead of orbit sums; fixed in w2lib, verified locally: d=8 sectors (even,+)=1, (even,-)=0, (odd,+)=0 with sector
  dims 255/228/240 matching the equivariant Hilbert series).  Resubmitted as e1kS (sectors d=8..12).
* d8F: crashed -- rational reconstruction of u8 failed after 6 primes (heights too large).  Resubmitted as d8G with a mod-p
  fallback; the mod-p relaxation verdict is rigorous for D8 by the regular-sequence argument (README §1.3, "D8 without an exact lift").
* relax7F: d=7 even system built (60 unknowns, 480 independent quadrics in 1830 pair monomials, good reduction in degree 14);
  msolve running (4 threads; remote load average 1300 on 128 cores at 01:10 UTC, so a multi-hour run is expected; a sibling
  lane's comparable d=7 msolve has >190 CPU-minutes).  Odd and mixed d=7 sectors queued behind it in the same job.

### R2. Sibling-lens data on this lane seen on the remote (search lens, not mine; to be re-verified)
* d = 9 E1 kernel (Y-even) primitive dimension 21 (p = 2147483647), hence 21 odd as well by duality (e1kernel_d9.out).
* D8 full fixed-D system inconsistent at p = 2147483647 and 2147483629 (rank 6703, nullity 794) (w2-d8d9.jsonl).
