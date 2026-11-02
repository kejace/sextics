# Wave 2 search report: `cycle-q2` (q = 2, A = X0^2 X1^2 X2^3 Y0^2 Y1^2 Y2)

STATUS: no witness found so far. New exact result: deg D >= 6 for every certificate (any D). Search infrastructure built
and validated; four remote job chains running since 19:43-19:50 local (box at load 500-750 from ~60 agent jobs, so
everything runs slowly). This file is rewritten after every batch. A parallel EXCLUSION agent (notes
`wave2-cycle-q2-excl.md`) and a second SEARCH agent running on the VM itself (`w2s_*.py`) work on the same lane; the
degree-5 theorem below was obtained independently by both of us; the d = 6 Groebner computation of the relaxation ideal
is theirs (job `excl-cycle-q2-hess6`); my `w2_relax6.py` reproduces its input (62-dim kernel, 426 quadrics, Psi(F) = 0).

Code (lane dir `/tmp/hodge/results/cycle-q2/`): `w2lib.py` (sector bookkeeping, E1-kernel per sector mod p and numeric,
Koszul-reduced sector-restricted complex variable projection `KoszulVarPro` with L-BFGS and Levenberg-Marquardt drivers,
rational reconstruction, exact confirmation through the wave-1 sector solver `eq.py`), `w2_search.py` (multiprocessing
driver over sectors; families kernel/jac/cone4/cone3/mult4/mult5/sing2/full), `w2_e1dims.py` (step 1), `w2_hess.py`
(Hessian relaxation at d = 5 exactly, J-type relaxation at d = 6). Log: `wave2-search.log.jsonl` (one line per batch).

## Method

E1 (`2 N(F) = A D`) is solved exactly on the family: for D in the sector E1-kernel K_kappa the lift N0(D) is linear in D.
Every lift is N0 + K(c) with the Koszul fields K(c) = sum_{i<j} c_ij (F_j d_i - F_i d_j), deg c_ij = d + 2, and E2 becomes
the affine problem  -B_D(N0) in F^2 S_{2d-6} + im T_D  (docs/05 §5.2).  For D(a) = sum a_k D_k the residual
f(a) = min_{c,Q} |T_{D(a)}(c) - F^2 Q + B_{D(a)}(N0(a))|^2 is minimized over the complex unit sphere (complex solutions are
not missed) with the variable-projection gradient (L-BFGS) or the Kaufman Jacobian (Levenberg-Marquardt).  All unknowns are
restricted to the character sector forced by kappa (lossless by averaging; 72 characters of the order-72 projective
diagonal group, parity included).  System sizes (rows x unknowns per sector): d=6 490x280, d=7 740x435, d=8 1120x670,
d=9 1650x1000, d=10 2360x1450, d=12 4500x2700.  Candidates with rel < 1e-10 are reconstructed over Q (common
denominators up to 2^20, continued fractions) and re-solved exactly (E1, full sector system mod 2^31-1 and a second prime,
CRT lift, `core.verify_witness`, F | D check).

Validation: gradient vs finite differences 5e-10 (cycle-q2 sector 96) and 8e-10 (Fermat); block-q2's witness has
f = 8e-32; fermat-q2's three strict witnesses (subset Euler, I = {X1,X2}, {X1,Y1}, {X2,Y1}) are exact zeros of f in the
6-dimensional pure-power family, and random-restart hit rate is 15-17 % per restart (L-BFGS 15/100, LM 17/100) -- the
landscape has spurious local minima, so dozens of restarts per sector are needed; LM reaches f ~ 1e-64 at true zeros.

## Results so far

### Theorem (exact, this wave): every certificate for cycle-q2 has deg D >= 6, for arbitrary D

E1 forces A D in J(F); the annihilator of [A] in R_d vanishes for d <= 5 (dossier A.3 at p = 10007; wave-1 sector
kernels at p = 2^31-1 give K_5 = J_5), so d <= 4 is impossible and d = 5 forces D = sum_i w_i F_i with constant w.
With the canonical lift N0 = (A/2) w, E2 modulo J(F) reads sum_{i,j} w_i w_j [A F_ij] = 0 in R_16 (im T_D, F^2 Q and
D div N0 all lie in J).  The 21 classes [A F_ij] (i <= j) fall into 21 distinct character sectors of R_16 and each is
nonzero mod 2^31-1 and mod 2^31-19 (`w2_hess.py d5`: rank 21/21); dim R_16 mod p = 951 for both primes, the
characteristic-zero value, so J_16 mod p has full rank and the independence transfers to Q and C.  Hence all w_i w_j = 0,
w = 0, D = 0.  No symmetry or parity assumption on D is needed.  (Closes gap 3 of the dossier; the previous bound was
deg D >= 5 with only semi-invariant degree-5 denominators excluded.)

### Step 1: E1-admissible denominator spaces per degree (mod 2^31-1, all 72 character sectors; `w2_e1dims.py`)

K_d = {D in S_d : A D in J(F)} decomposes over the 72 characters of the projective diagonal group (parity included).

| d | dim S_d | dim K_d | dim J_d | non-Jacobian (K mod J) | Y-even: S / K / J | sectors with K > J | max sector dim K |
|---|---:|---:|---:|---:|---|---:|---:|
| 4 | 126 | 0 | 0 | 0 | 66 / 0 / 0 | 0 / 69 | 0 |
| 5 | 252 | 6 | 6 | 0 | 126 / 3 / 3 | 0 / 72 | 1 |
| 6 | 462 | 62 | 36 | 26 | 236 / 36 / 18 | 19 / 72 | 6 (trivial sector = J) |
| 7 | 792 | 546 | 126 | 420 | 396 / 273 / 63 | 72 / 72 | 12 |
| 8 | 1287 | 1161 | 336 | 825 | 651 / 591 / 168 | 72 / 72 | 22 |
| 9 | 2002 | 1946 | 756 | 1190 | 1001 / 973 / 378 | 72 / 72 | 32 |
| 10 | 3003 | 2982 | 1497 | 1485 | 1512 / 1503 / 750 | 72 / 72 | 48 |
| 11 | 4368 | 4362 | 2682 | 1680 | 2184 / 2181 / 1341 | 72 / 72 | 66 |
| 12 | (running) | | | | | | |

The non-Jacobian dimensions reproduce the annihilator dimensions of the dossier (26, 420, 825, 1190, 1485, 1680 =
dim R_d - dim R_{d+12}, i.e. multiplication by [A] is surjective onto R_{d+12} for d >= 8).  From d = 8 on, E1 is almost
vacuous (codimension 126, 56, 21, 6 in S_d), so the per-sector search spaces are 16-66 dimensional.

### Rank / obstruction statistics of the reduced E2 systems (generic point of the sector kernel family; `rank-stats` log lines)

| d | sector key | family dim K | rows | unknowns (c_ij + Q) | rank | Koszul gauge nullity | linear obstructions (rows - rank) |
|---|---|---:|---:|---:|---:|---:|---:|
| 6 | 0 (trivial, = J part) | 6 | 480 | 271 + 7 | 245 | 32 | 235 |
| 6 | 96 | 2 | 492 | 271 + 9 | 247 | 32 | 245 |
| 6 | 597 | 3 | 471 | 265 + 6 | 237 | 34 | 234 |
| 7 | 6 | 10 | 750 | 409 + 18 | 380 | 46 | 370 |

So at d = 6 the right-hand side -B_D(N0) must satisfy ~240 independent linear conditions while the family has <= 6
parameters: a solution needs a genuine geometric reason (this quantifies why every tested D fails with rank/aug = +1).
The Koszul unknowns carry a 32-34 dimensional gauge (second syzygies), handled by the rank-revealing least squares.

### Steps 2-5: joint numerical search (remote jobs)

Mid-course corrections (20:00-20:40): the box reached load 1000-1300 (70 agent sessions on 128 cores), so the first chains
(300/150/60 restarts per sector at d = 6/7/8, 200 LM evaluations per restart, 8-16 workers, 50 structured sub-family
invocations per degree) made almost no progress in an hour and were killed and resubmitted twice with 32 worker processes,
at most 60 LM evaluations per restart (LM converges in 20-35 at true zeros, see the Fermat calibration), fewer restarts,
and a `multi` mode that builds each sector system once and searches all structured sub-families on it through linear
sub-family views (`KoszulVarPro.subfamily`).  Second, decisive correction (21:30): the workers were running 32 BLAS threads each (the runner exports
`OPENBLAS_NUM_THREADS=<cpus>` and `w2_search.py` imported numpy before the thread cap), i.e. ~1000 runnable threads for my
four chains alone (load 1800, 11 % CPU per worker, no finished batch in an hour).  All scripts now force single-threaded
BLAS before importing numpy (see `/tmp/hodge/notes/engine-bugs.md`), and the chains were resubmitted as A4/B4/D3/E3.
Third correction (22:45): with ~80 agent sessions the box stays at load ~600 and each of my workers gets ~6 % of a
core (a 10-restart batch at d = 7 takes 52 minutes of wall time), so chain A was reduced once more (A5: d = 7 with 10
restarts per sector, skipping the 33 sectors already sampled; d = 8 with 5; d = 9 with 3; 48 workers).  Coverage at
d >= 7 is therefore a sample (10-15 % hit probability per restart for an isolated zero in the Fermat calibration), not
an exhaustive search; the per-sector best residuals are still informative about the landscape.
Final job set:

* chain A3 `search-cycle-q2-chainA3`: kernel d = 6 (40 LM restarts per sector, 72 sectors) -> trivial sector d = 6 (the
  6-dim family D = sum a_i Z_i F_i whose Jacobian relaxation is the degree-8 surface V0 of the excl notes, 400 extra
  restarts) -> kernel d = 7 (30) -> d = 8 (15) -> d = 9 (8).
* chain B3 `search-cycle-q2-chainB3`: structured sub-families inside every sector kernel at d = 7, 8, 9, 10: cone4 (15
  line vertices), mult4 (15 coordinate lines: multiplicity >= 4 at the six points of X ∩ line, a linear condition over Q
  because the six points form one G-orbit), cone3 (20 plane vertices); 15 restarts each.  Then kernel d = 10 (12
  restarts), the product/square relaxation prod6 at d = 12 (span of products of degree-6 kernel elements of complementary
  characters), kernel d = 11 (8) and d = 12 (6).
* chain D2 `search-cycle-q2-chainD2`: modular sector tests (two primes) of 420 structured semi-invariant candidates:
  117 Hessian 2x2 minors (d = 8), 210 Hessian 3x3 minors (d = 12), 63 subset-Euler forms E_I(F) (d = 6) and their squares
  (d = 12), block-q2 cubics in each pair, pair-Fermat products, F_i F_j; then an exact scan of the nine 2-dimensional
  E1-kernel sectors at d = 6 (all F_p-rational points of each pencil through a maximal-minor polynomial, `w2_pencil.py`).
* chain E `search-cycle-q2-chainE`: q = 1 sibling (F, B = X0 X1^2 Y0 Y1^2) kernel search at d = 8..12 (systems
  ~300 x 75 at d = 8).  Its first results (d = 9, 21:00) were 512 "zeros" with rel ~ 1e-30 in 8 sectors, all of which turned
  out to be the trivial q = 1 solutions D = F G (N = 0, Q = B G^2; F S_3 lies in every sector kernel): the numerical
  certificate is exact (full identity residual 6e-14) but D is divisible by F to 1e-7 and every rational reconstruction is
  inconsistent.  At q = 1 this linear subspace of spurious solutions attracts every restart, so the sibling search was
  restarted (chain E2) with a deflated residual r / dist(a, L)^2, L = {a : F | D(a)}, which removes the attractor without
  discarding any legal denominator (projecting the family onto a complement of L would).  At q = 2 F-multiples are not
  automatically consistent (D = F G needs [A][G] = 0 in R, i.e. E1 for G), and no such attractor was observed.

### Results of the joint search so far (22:00; updated as jobs finish)

Structured modular sweep (chain D3, two primes, sector solver `eq.test_D_fast`): 110/117 Hessian 2x2 minors tested at
d = 8 (104 fail E1, 6 pass E1 and are inconsistent, ranks 909-1062), 201/210 Hessian 3x3 minors at d = 12 (199 pass E1,
all inconsistent, ranks 2614-3510), 24/63 subset-Euler forms E_I(F) at d = 6 (all pass E1, all inconsistent, ranks
595-629) and 16 of their squares at d = 12 (inconsistent, ranks 2818-3003); block cubics, F_i F_j and the pencil scans
follow.

| q | target | d | family | sectors | restarts | system (rows x unknowns) | best rel over all sectors | median of sector-best rel | candidates (rel < 1e-10) |
|---|---|---:|---|---:|---:|---|---:|---:|---:|
| 2 | A | 6 | kernel (all 48 nonempty sectors, K = 1..6) | 48 | 1920 | <= 492 x 304 | 7.5e-3 (key 777, K = 1: a fixed D, its least-squares defect) | 0.76 | 0 |
| 2 | A | 6 | trivial sector D = sum a_i Z_i F_i (K = 6) | 1 | 40 + 2704 | 480 x 278 | 0.18 (median over 32 independent parts 0.24) | - | 0 |

| 1 | B (sibling) | 8 | kernel, deflated (chain E3, finished) | 72 | 10800 | <= 297 x 79 | 4.5e-7 (best non-spurious) | 0.73 | 0 |
| 1 | B (sibling) | 9 | kernel, deflated (chain E3, running) | 16 of 72 | 700 | <= 492 x 150 | 9.2e-7 (non-spurious) | 4.4e-6 | 0 |
| 2 | A | 7 | kernel (chain A4 sample; A5 completes the other 40 sectors) | 32 of 72 | 320 | <= 766 x 458 | 1.2e-3 | 0.040 | 0 |
| 2 | A | 7 | mult4 sub-families (mult >= 4 at X ∩ coordinate line; chain B5) | 65 | 885 | <= 766 x 458 | 9.9e-4 | 0.097 | 0 |
| 2 | A | 7 | cone4 sub-families (line-vertex cones; chain B5) | 42 | 531 | <= 766 x 454 | 0.26 | 0.86 | 0 |
| 2 | A | 8 | cone4 sub-families (chain B5, running) | 19 | 285 | <= 1152 x 690 | 0.17 | 0.51 | 0 |

Sibling note: with the deflated residual r/dist(a,L)^2, 20-40 % of the restarts in sectors containing F-multiples still
end within 5 % of L (flagged spurious; e.g. sector 609: 51 of 150 restarts below 1e-10 are all spurious), and the best
non-spurious values 4e-7..6e-7 (sectors 12, 609, 528) sit just outside the 5 % exclusion radius, consistent with the
f ~ dist^4 funnel of the trivial solutions rather than with isolated certificates; no non-spurious restart reached 1e-10.

Interpretation: rel = min_{c,Q} |E2 residual|^2 / |B_D(N0)|^2 on the unit sphere of the family.  Values of order 0.2-0.8
mean that the best Koszul correction removes only 20-80 % of the obstruction norm; a true certificate gives rel ~ 1e-25
(Fermat calibration).  The two-dimensional sectors (keys 486, 876, 1062, 96, 576, 936, 516, 822; K = 2) have minima
0.18-0.51, the three-dimensional sector 597 has 0.54, the trivial six-dimensional sector 0.23: no degree-6
semi-invariant denominator comes close, in agreement with the exact exclusions (Farkas certificates of the EXCLUSION agent,
and the pencil scans of chain D3).

