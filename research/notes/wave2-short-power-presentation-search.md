# Wave 2 search report — short-power-presentation (q = 1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (wave 2, 2026-11-01 19:10 – ongoing; this file is rewritten after every batch)

Log: /tmp/hodge/results/short-power-presentation/wave2-search.log.jsonl (one JSON line per batch/experiment).
Scripts (all in the lane dir): w2lib.py (exact mod-p machinery), w2_kernel.py (E1 kernels + multi-prime reconstruction),
w2_relax.py (E2-mod-J quadratic relaxation over K_d + J_d, gauge-fixed, msolve/Singular), w2_relax_sub.py / w2_subfam.py
(relaxation restricted to structured subfamilies), w2_varpro.py (reduced variable projection with Koszul-correction unknowns,
deflation), w2_families.py (exact subring-family screen), w2_cands.py (named candidates + kernel structure), w2_polar.py,
w2_annl.py, w2_annpow.py (annihilator surveys), w2_fulltest_modp.py + msparse.py (mod-p full test of relaxation points).
Remote jobs (hodge.rjob): search-spp-j1 (done); search-spp-kA / search-spp-kBC queued with --wait since 00:20 UTC — the remote
has been refusing all submissions (load 300–1150 on 128 cores, 40–93 jobs of other agents).  Everything below after job j1 was
therefore computed locally on brunos (1–2 threads per process, load 60–100).

## 1. Exact results (mod p, two to five primes agree; p in {536870909, 536870879, 536870869, 2^31-1, 2147483629})

### 1.1 E1-admissible spaces K_d = Ann_R(A)_d (E1 kernel modulo J(F)_d)
The E1 kernel in all degrees is the annihilator ideal Ann_R(A) of the Gorenstein Jacobian ring R (socle degree 24; A odd), so
dim Ann_d^+ = dim R_d^+ - (dim R_{18-d}^- - dim Ann_{18-d}^-), dim Ann_d^- = dim R_d^- - (dim R_{18-d}^+ - dim Ann_{18-d}^+).

| d | dim R_d (+,-) | dim J_d (+,-) | dim K_d mod J (+,-) | how |
|---|---|---|---|---|
| 1..6 | ..,(218,208) | ..,(18,18) | (0,0) | direct, wave 1 (two primes) |
| 7 | (333,333) | (63,63) | (1,1) | direct; D7 ∈ C[X0,X1,X2,M], D7' ∈ C[Y0,Y1,Y2,L] |
| 8 | (483,468) | (168,168) | (6,6) | direct, 5 primes; odd = S_1·K_7; even = 4 + 2 new (sums of two cones) |
| 9 | (623,623) | (378,378) | (18,18) | direct, 4 primes each; K_9 = S_1·K_8 exactly (no new elements) |
| 10 | (762,744) | (750,747) | (300,267) | duality from d=8 |
| 11 | (843,843) | (1341,1341) | (511,511) | duality from d=7 |
| 12 | (885,866) | (2223,2214) | (677,648) | duality from d=6 |

So E1 is a strong filter only for d <= 9; at d >= 10 most of R_d is annihilated by A.

### 1.2 Structure of the kernel elements (all exact)
* D7 (168 terms, rational, heights 27/21 bits): irreducible, depends on Y only through M = Y0+3Y1+5Y2, M-degrees {0,2,4}:
  a cone with vertex the line {X=0, M=0} (the intersection of Λ with the Fermat coordinate plane z_4 z_5 z_6).
  Full fixed-D solve (remote, engine v2): inconsistent (rank 4035, nullity 462); served space = J_6^- only.
* D7' (odd): in C[Y0,Y1,Y2,L], L-degrees {0,2,4,6}: cone with vertex the line {Y=0, L=0}; canonical cone representative in
  D7odd_cone.json (42 coefficients); rational NF representative in kernel_d7_p1.json (3 primes).
* Jacobian structure: 2F_X0 - F_X1, 3F_X0 - F_X2, 3F_Y0 - F_Y1, 5F_Y0 - F_Y2 cancel the coupling, so in R:
  X1^5 = 2X0^5, X2^5 = 3X0^5, Y1^5 = 3Y0^5, Y2^5 = 5Y0^5, X0^5 = -H(L,M)/6, Y0^5 = -G(L,M)/6 (H = d_L C, G = d_M C,
  C = c((L+M)^6+(L-M)^6)); equivalently the eight fifth powers l_k^5 span the 2-dim space <alpha, beta> of R_5.
* Polar residuals: L·D7 ∈ J_8 and M·D7' ∈ J_8 (these are exactly the two relations among the six products X_i D7, Y_j D7' mod J_8),
  i.e. D7 = w(F)/L, D7' = w'(F)/M are residual components of reducible polar divisors.  Likewise Ann_R(L)_6^+ = <E6> with E6 a
  degree-6 cone over {X=M=0} (E6.json) and Ann_R(L)_7^+ = span(X_i E6 [2-dim], D7).
* d = 8: K_8^- = S_1 K_7 exactly.  K_8^+ ⊂ C[X,M]_8 + C[Y,L]_8: four directions from S_1 K_7 and two NEW elements which are genuine
  sums P(X,M) + Q(Y,L) (rational, 5 primes, K8even_twocone_reps.json).  Trivial but useful: Ann_R(X0) + Ann_R(Y0) ⊂ Ann_R(A);
  dim Ann_R(X0)_8^+ = dim Ann_R(Y0)_8^+ = 1 (the Lefschetz maps by the coordinate forms fail to be injective at d = 8).
* Annihilator surveys (mod p): Ann(L)_d^+ = 1,3,11,22,48,71 (d=6..11), Ann(M)_d^+ = 1,3,11,25,48,77; Ann(X0)_d^+ = 0,0,1,1,7,7;
  Ann(Y0)_d^+ = 0,0,1,3,7,13; odd: Ann(X0)_9^- = 3, Ann(Y0)_9^- = 1, Ann(L)_8^- = 9, Ann(L)_10^- = 42, Ann(X0)_10^- = 3, Ann(X0)_11^- = 13.
* Monomial annihilators (w2_annpow.py, mod p): D7 ∈ Ann_R(X0^3) (X0^3 D7 ∈ J_10) and D7' ∈ Ann_R(Y0^3) — the degree-7 cones are
  killed by a single cube.  At d = 8 (even): Ann(X0^k)_8 = k, Ann(Y0^k)_8 = k (k=1,2,3), Ann(X0^a Y0^b)_8 = a+b-? (table in the log:
  X0Y0:2, X0Y0^2:3, X0Y0^3:4, X0^2Y0^2:4, X0^2Y0^3:5, X0^3Y0^2:5, X0^3Y0^3 = K_8^+:6); odd d = 8: Ann(X0^a Y0^b)_8^- = 3 for every
  multiplier with a = 3 or b = 3.  So K_8^+ = Ann(X0^3)_8 + Ann(Y0^3)_8 (4-dim) + two elements needing both cubes.
* Named candidates: all 2x2 minors of Hess(F) (99 even, 72 odd at d=8) fail E1; E_X(F), E_Y(F) times forms lie in J.

### 1.3 Gauge invariance (proved, used everywhere)
If (D, N, Q) is a certificate normalized by N(F) = A D (always possible), then for every h ∈ S_{d-6}
(D + hF, N + hA E/6, Q + h div N - N(h) + h^2 A) is a certificate.  Hence denominators are searched modulo F S_{d-6}
("gauge-fixed": the coordinates (Y2-component, h·Y2) of the Jacobian field are removed).

### 1.4 E2-mod-J relaxation (exact necessary condition): for D = sum a_k G_k with particular numerators P_k (P_k(F) = A G_k),
B(a) = sum_{k,l} a_k a_l (G_k div P_l - P_l(G_k)) must lie in J_{2d}; Koszul corrections and Q drop out mod J.
* d = 6 (pure Jacobian family, 17/18 gauge-fixed variables, ~600 quadrics): no projective solutions (msolve [-1]; Singular dim 0,
  vdim 1216 / 540) for both parities  =>  no certificate with deg D = 6 of pure parity (and none of mixed parity <= 3).
* d = 7 (64 generators = D7 + J_7, 61 variables after gauge, 734 even / 725 odd quadrics): the full Gröbner bases
  (msolve F4 on a random affine chart; Singular std; msolve leading ideal) have been running for > 1 h CPU each with 7–11 GB:
  the degree-3 Macaulay step did not certify emptiness, degree 4 (635k columns) is in progress.  Partial exact results:
  - D7 alone fails (696/734 equations nonzero); the pencils span(D7, w(F)) for the polar-residual fields w(F) = L E6 (even) and
    w(F) = M E6M (odd) contain no solutions;
  - Fermat subfamily Λ = {w_j = z_j q_j} + <D7> (16 variables; the family that solves the relaxation identically at c = 0):
    only the origin (Singular dim 0, vdim 194 / 191; msolve chart [-1]) for both parities;
  - subring subfamilies ∩ (K_7+J_7) (w2_subfam.py; 2..5 forms, two-cone family): 80 (even) / 81 (odd) families done, every one
    with msolve chart "[-1]" (no projective solutions), including the full cone families <X0,X1,X2,M> (22 vars, contains D7) and
    <Y0,Y1,Y2,L> (22 vars, contains D7') and the 24–27-variable 5-form families; the remaining larger families are running
    (restarted with a 20 s Singular cap since msolve's answer is decisive).
  - the full 61-variable GB: the local msolve runs exceed the memory available on brunos at the F4 degree-4 step (the even run died
    with the job timeout/OOM; the odd run is at ~70 CPU-min, 7 GB); a 16-thread remote job (search-spp-gb7) is queued.
* d = 8 (174 generators, 162 variables, 468/483 equations): the full GB is hopeless (13203 quadratic monomials vs ~470 equations);
  the equations are being built (both parities; reducers at degree 16 done) and a watcher then runs the subfamily restrictions.
* d = 7 mixed parity (128 generators, 122 variables, 1475 equations): built; subfamily restrictions running (w2_subfam2.py).

## 2. Structured families screened exactly against E1 (w2_families.py, w2_cands.py)
152 (even) / 134 (odd) subring families at d = 7, 45+ at d = 8 (deduplicated by span), plus the two-cone family and the named
candidates.  At d = 7 the only families meeting K_7 + J_7 outside J_7 are <X0,X1,X2,M> (even, D7) and <Y0,Y1,Y2,L> (odd, D7').
At d = 8: <X0,X1,X2,M> and <Y0,Y1,Y2,L> contain only the products X_i D7 / Y_i D7'; the composite family C[X,M]+C[Y,L] contains all of
K_8^+.  Many families meet J_d nontrivially (binomial relations above) — those are Jacobian-type and covered by the w-space.

## 3. Numerical joint search (reduced variable projection, w2_varpro.py)
Formulation: unknowns for fixed a are only the Koszul coefficients c_ij (degree d-4) and Q; rows = degree-2d monomials; each
evaluation is a 5796 x 1056 least-squares at d = 7 (vs 26532 x 4497 in the engine's VarProFamily).  Validated against direct
polynomial arithmetic (1e-13), finite-difference gradient (ratio 1.000000) and the Fermat control certificate.
* d = 7 even, 6 restarts (2 procs, maxiter 300): every restart converges to the SAME basin f/scale0 ≈ 3.2e-6 (f ≈ 1.1e-3), ~100
  iterations, 10–15 min per restart at load 60–100.  The basin is the Fermat certificate X_i E_Y(F) (direct evaluation: f/scale0 =
  2.0e-6..2.6e-6 at w = X_i E_Y, 0.9e-5..1.2e-5 at X_i E_X; D7 alone: 4e-4).  Explanation: at c = 0 the family w_j = z_j q_j solves
  the relaxation identically and contains genuine Fermat certificates; the coupling c = 5/7^7 breaks them only at relative size ~c.
  Hence residuals ~1e-6 are expected artefacts; no candidate reached the confirmation threshold (1e-9) and none was rationally
  reconstructed.  Deflation was added for the remote batches (queued).

## 4. What is excluded, what is open
* Rigorous (mod p, several primes): no certificate with deg D <= 6 of pure Y-parity; deg D <= 3 mixed.  deg D = 7: E1 forces
  D ∈ <D7> + J_7 (even) or <D7'> + J_7 (odd); D7, D7' alone and the Fermat/polar-residual subfamilies are excluded; the full
  61-variable relaxation is undecided (GB running).
* Open: the full d = 7 families, d = 8 (162-dim families), mixed parity at d = 7 (122 variables; build running), d >= 9.
