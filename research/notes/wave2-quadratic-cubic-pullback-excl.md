# Wave 2 — lane `quadratic-cubic-pullback` — lens EXCLUSION / RIGOR

STATUS: no witness.  BEST RIGOROUS BOUND (exact, Hilbert-certified): every certificate for (F, 1, X0^3 Y0^3) has deg D >= 7 (K_d = J_d for d <= 6 in both parities; the complete polar families of degree 5 and 6 — all parity mixtures — are excluded exactly: degree 5 by pi_W(ker) = 0 / K_H = 0, degree 6 by the exact Hessian relaxation whose candidate cone is exactly the Euler line D ~ F).  17 structured degree-7 Jacobian-type candidates excluded rigorously.  Remote jobs never ran (box saturated, load 300-1600).

Started 2026-11-01 (wave 2).  Log: `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl.log.jsonl`.
Packet: `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/`.

## Results so far
- Local check: F = sum Q_i^3 exactly; the fixed-D kernel of a polar denominator D = dF/dX0 (d=5) has the explicit
  rational basis (G z, 0), G in S_5^odd (126) and (g F, D div g - g.grad D), g constant on the Y-slots (3): 129 elements,
  equal to the modular nullity (engine reports 128 = rank of augmented matrix convention).  Consequence: the three
  wave-1 degree-5 polar inconsistencies are rigorous over Q (modular rank of M = rational rank, sandwich argument).
- Engine note written to /tmp/hodge/notes/engine-bugs.md (reported nullity off by one for inconsistent systems).
- Dossier correction: E_I F = sum_{i in I} z_i dF/dz_i lies in J(F) for EVERY index set I (not only I = all), so the
  block-Euler denominators EI([0,1,2]) etc. listed in wave-1 `listC` were never tested; they are natural
  degree-6 polar candidates (the Fermat join X0^6+X1^6+X2^6 is exactly E_X F_0 / 6) and are queued now.

## Re-verification of wave-1 claims (goal a) — all at prime p2 = 2147483629 (wave 1 used 2147483647)
| claim | wave 1 | wave 2 (p2) | status |
|---|---|---|---|
| h^{3,1}(X)^+ = 226 - 18 = 208; A, delta not in J_6, independent mod J_6 | yes | J6_odd rank 18 = Hilbert, A,delta independent | rigorous (agrees) |
| projective diagonal symmetry group = Z/2 | Smith form (sympy) | engine `symmetry.diagonal_group`: order 2 | agrees |
| F not a join: rank of 21 second partials = 21 | p1 | p2: 21 | rigorous (agrees) |
| cubics in Q meet J_6^even in exactly C.F | p1 | p2: intersection dim 1 | agrees (mod-p dim is an upper bound, F gives the lower bound) |
| two-plane (4,4)-cancellation system has nullity 0 | sympy exact | fmpq exact: rank 12 of 12 | rigorous (agrees) |
| delta, join012, join345, F60, F06, fermatX, Q0^3+Q3^3, Q0^2, F-Q0^3 not in J_6 (so excluded at every pole order by K_6=J_6) | p1 | p2: all not in J, J-rank certified | rigorous (agrees) |
| E1 kernel K_d = J_d (even d<=8, odd d<=6) | p1 | p2 (smoke: d=3,5 both; 6 odd) — remote job covers d<=12 | agrees so far |
| three degree-5 polar denominators serve nothing | p1, modular | structural kernel 129 = modular nullity => rigorous over Q | upgraded to rigorous |
| E_X F etc. "in J only for I = all" (dossier §3.3) | — | E_X F is in J_6 (certified) | dossier statement corrected |

## New rigorous results
* **No Y-even degree-5 denominator** (local, exact): polar family D = sum_{i<3} w_i F_i, relaxation matrix 1491 x 201
  (6 W + 135 V + 60 Q unknowns), exact integer nullspace has dimension 12 and projects to ZERO on the W-coordinates,
  so no w gives a certificate.  Self-test of the E2 reduction E = F (Phi(w) + L(c,w)) passed exactly.
* gradX0sq = 2 X0 F_0 (d = 6, reducible control): serves nothing, rank-certified (nullity 244 = structural).

## Tools and statements developed in this lens (mathematics)

### T1. Certification lemma for modular fixed-D results
Let M be the (N,Q)-coefficient matrix of the fixed-D system over Q and M_p its reduction.  rank M_p <= rank M.
If an explicit rational kernel of M of dimension n_str = ncols - rank M_p is known, then rank M = rank M_p and
(i) "M_p x = b_p inconsistent" implies inconsistency over Q; (ii) the modular served dimension (= #targets -
(rank[M|B]_p - rank M_p)) is an UPPER bound for the rational served dimension.  Explicit kernel for every D
(deg d, N_i of parity par(N_i)):  (G z, 0) with G in S_d of parity par(A)-par(D)  [E(Gz,0) = 0 identically],
(g F, D div g - g.grad D) with g_i in S_{d-5} of parity par(N_i), overlap {G = hF, g = h z}.
n_str(d) = dim S_d^{parG} + sum_i dim S_{d-5}^{par N_i} - dim S_{d-6}^{parG}: 129 (d=5), 244 (d=6).
A fixed-D verdict is logged "rigorous" only when the modular nullity equals n_str.

### T2. E1 rigor
"A D not in J(F)" mod p is rigorous when the modular rank of the degree-(d+6) Jacobian block equals the Hilbert
value dim J_{d+6}^{par} = dim S - dim R (parity-split Hilbert series (1-t^5)^3(1-ut^5)^3/((1-t)^3(1-ut)^3), u^2=1;
totals 1,6,21,56,126,246,426,666,951,1246,1506,...; splits R_6 = 218+208, R_10 = 762+744, R_16 = 483+468).
Kernel statement: dim K_d(p) >= dim K_d(Q) >= dim J_d, so K_d(p) = J_d (with certified J-rank) proves K_d(Q) = J_d.

### T3. Jacobian-type (polar) reduction and the exact relaxation
If K_d = J_d then every admissible D of degree d is D = w.grad F.  Write N = A w + H z/6 + c.grad F (complete:
(N - A w).grad F in (F), Koszul syzygies).  Then E = F.(Phi(w) + L(c,w)) with
   Phi(w) = D (grad A.w + A div w) - A w.grad D,      L(c,w) = sum_{i,j} F_j [(d_i c_ij) D - c_ij d_i D],
all H-terms cancel, and the certificate condition is Phi(w) + L(c,w) = F Q.  (Checked exactly on random data.)
Relaxation: W = w w^T, V = c (x) w independent => linear system in (W, V, Q).  Admissible w satisfy
w w^T in pi_W(ker).  pi_W(ker) = 0 (exact integer nullspace) proves that NO polar denominator of that degree/sector
admits a certificate for A.  Result at d = 5: see below.

### Result R1 (rigorous, exact): no denominator of degree 5
Sector even (w = (w0,w1,w2,0,0,0)): relaxation 1491 x 201, exact nullity 12, pi_W(ker) = 0.
Sector odd (w = (0,0,0,w3,w4,w5)): 1491 x 201, exact nullity 12, pi_W(ker) = 0.
Unrestricted w in C^6: 3003 x 687, exact nullity 156, pi_W(ker) = 0.   (log entries relaxL5o, relaxL5a, smokerelax)
Together with K_5 = J_5 (both parities; mod-p kernel dimension is an upper bound, J-rank Hilbert-certified) this
proves: no certificate with deg D = 5 exists.  With K_d = 0 for d <= 4: deg D >= 6 for every certificate.

### Hilbert bookkeeping (exact; Koszul series of the regular sequence dF, parity-graded)
| d | R_d^even, R_d^odd | R_{d+6}^odd, R_{d+6}^even | Gorenstein lower bound on dim K_d/J_d (even, odd) | dim J_d (even, odd) | structural fixed-D nullity n_str (even, odd D) |
|---|---|---|---|---|---|
| 5 | 123, 123 | 843, 843 | 0, 0 | 3, 3 | 129, 129 |
| 6 | 218, 208 | 866, 885 | 0, 0 | 18, 18 | 244, 253 |
| 7 | 333, 333 | 843, 843 | 0, 0 | 63, 63 | 456, 456 |
| 8 | 483, 468 | 744, 762 | 0, 0 | 168, 168 | 795, 807 |
| 9 | 623, 623 | 623, 623 | 0, 0 | 378, 378 | 1351, 1351 |
| 10 | 762, 744 | 468, 483 | 294, 261 | 750, 747 | 2187, 2202 |
| 11 | 843, 843 | 333, 333 | 510, 510 | 1341, 1341 | 3444, 3444 |
| 12 | 885, 866 | 208, 218 | 677, 648 | 2223, 2214 | 5230, 5248 |
(dim K_d/J_d computed mod p is an upper bound; when it equals the Gorenstein lower bound the rational value is pinned.)

### T4. Hessian relaxation (rigorous, Hilbert-certified) — the certifiable form of T3
Since L(c,w) and F Q lie in J(F), E2 forces Phi(w) in J(F)_{2d}, and Phi(w) = -A w.grad D (mod J)
= -A sum_{k,l} w_k w_l F_kl (mod J).  So every admissible polar w satisfies  **A (w^T Hess(F) w) in J(F)_{2d}**.
Relaxing W_kl = w_k w_l gives a linear condition whose kernel K_H is computed mod p with the Jacobian block at its
Hilbert rank, hence dim K_H(Q) <= dim K_H(p) rigorously (same logic as the E1 sweeps).  Script `wave2_excl_hess.py`.
Degree 6 (p = 268435399, rank J_12^odd = 2214 = Hilbert, J_12^even = 2223 = Hilbert):
| sector | w-coords | W-unknowns | dim K_H (upper bound) | structure |
|---|---|---|---|---|
| Y-even D | 18 | 171 | 36 | contains the 18 exact "fake" directions e_a (.) u (u = z/6, D = F) |
| Y-odd D | 18 | 171 | 18 | only antisymmetric off-diagonal 2x2 minors: no diagonal W_aa => w = 0 |
| cross (mixed parity) | 36 | 324 | 126 | coefficients +-1 |
Full relaxation (T3, mod 4 primes): pi_W(ker) has dim 19 (even) / 18 (odd); the odd one coincides with K_H.
msolve note: msolve silently miscomputes for primes near 2^31 ("Prime ... too large" on small inputs); use p < 2^28.

### Result R2 (rigorous, exact): no Y-pure denominator of degree 6
Every degree-6 admissible D is a polar D = w.grad F (K_6 = J_6, both parities, primes 2147483629 and 2147483587,
Jacobian blocks at Hilbert rank).  Hessian relaxation (T4) at p = 268435399 (J_12^odd rank 2214 = Hilbert):
* Y-odd sector: dim K_H = 18 (upper bound); the 18 symmetric-lifted basis vectors (entries in {0,+-1}) are verified
  EXACTLY to satisfy A sum W_ab m_a m_b F_ab in J_12 (exact rational solve on a 2214x2214 pivot system + full check),
  so K_H(Q) is exactly their span; the exact annihilator (153 quadrics) has, by Singular over Q (std + minAssGTZ),
  the single zero-dimensional component {w = 0}.  => no Y-odd polar denominator of degree 6.
* Y-even sector: dim K_H = 36 (upper bound, same at p = 2147483629); all 36 lifted vectors (entries in {0,+-1,+-2})
  verified exactly; the exact annihilator (135 quadrics) defines, by Singular over Q, exactly ONE component of
  dimension 1: the Euler line w2 = w4 = w6 = w11 = w13 = w15 (the six coordinates z_i e_i), all other coordinates 0,
  i.e. w = lambda z, D = 6 lambda F, which is not coprime to F.  => no Y-even polar denominator of degree 6 coprime to F.
Scripts: wave2_excl_hess.py (modular, certified upper bound), wave2_excl_hess_exact.py (exact lift + Singular);
receipts: wave2-excl.hessL6.jsonl, wave2-excl.hessX6.jsonl, KH files wave2-excl-hess-hessL6-d6-par*-KH.json,
Singular scripts wave2-excl-hess-exact-hessX6-d6-par*.sing.
Caveat found and fixed: msolve parses `w*w` as `w` (and rejects primes near 2^31); every msolve-based count in this
lens was discarded and redone with Singular over Q.

### T5. Clean form of the Hessian condition, and the exact rigor logic
For a polar denominator D = d_w F (d_w = w.grad the derivation along the polynomial vector field w), E2 implies
      A . d_w D = A . d_w^2 F  in  J(F)   (degree 2 deg D).
(Proof: Phi(w) = D(grad A.w + A div w) - A w.grad D with D in J, and w.grad D = sum_k (w.grad w_k) F_k + d_w^2 F.)
Rigor of the modular computation: a rational w with A d_w^2 F in J(F) reduces to a solution mod p whenever the Jacobian
block has its Hilbert rank mod p (p-integral Cramer solution), so dim K_H(Q) <= dim K_H(p); the exact lift of the
mod-p basis (verified by an exact solve) then pins K_H(Q) and the candidate cone {w : w w^T in K_H(Q)} exactly; Singular
over Q decides the cone.  A Singular computation over F_p alone (no exact lift) is screening: V(F_p-bar) = empty does
not prove V(Q-bar) = empty (finitely many bad primes).
Trivial/fake directions: u = z/6 gives D = F (not coprime), and u (.) w2 is in K_H for every w2 (polarization identity,
verified exactly); the Y-odd d=6 kernel consists of the 18 antisymmetric "minor" directions e_(i,Y_a) (.) e_(j,Y_b)
- e_(i,Y_b) (.) e_(j,Y_a) (and X-analogues), which have no diagonal part, forcing w = 0.

### Degree 5 via the Hessian relaxation alone (rigorous, 3 s): K_H = 0 in the even, odd and unrestricted sectors
(p = 268435399, J_10 blocks at Hilbert rank 747/747/1497).  So "no degree-5 denominator" needs only a certified
rank computation, independently of the exact full relaxation of R1.

### Degree 6, mixed parity (joint 36-coordinate relaxation): dim K_H = 261 (certified upper bound; = 36 + 18 + 126
+ 81 coupling directions between the two Y-odd blocks), coefficients in {0,+-1,+-2}; exact lift + Singular over Q
running (receipt wave2-excl.hessX6.jsonl when done).

### Degree 7 (Y-even polars): dim K_H = 1533 of 2016 (certified upper bound), basis coefficients large (~1e8 mod p), so
the exact lift is not cheap; the cone {w : A d_w^2 F in J} is cut by 483 dense quadrics in 63 variables; Singular std (stopped after ~30 min) and
msolve (two charts) over F_p did not finish locally.  Status: inconclusive (neither excluded nor screened); the exact
lift is out of reach locally (basis coefficients ~1e8 mod p).

## Summary of the lens (to be finalized)
1. Best rigorous bound: deg D >= 6 for every certificate (E1 kernels + exact relaxation / Hessian relaxation at d = 5);
   no Y-pure degree-6 denominator coprime to F (exact); mixed-parity degree 6: [[pending exact run]] -> deg D >= 7.
2. All wave-1 structural claims re-verified at a second prime / exactly; three wave-1 degree-5 polar results upgraded
   to rigorous; one dossier error corrected (E_I F in J for every I).
3. New tools: certification lemma for fixed-D verdicts (structural kernel), E2 reduction for polars (H cancels),
   exact relaxation (d=5) and Hilbert-certified Hessian relaxation A d_w^2 F in J (d = 5, 6, 7 computed).
4. Not executed (remote saturated all session, load 300-1060): E1 kernels d = 7..9 odd / d = 9, 10, 11, 12 both parities
   at new primes; 74 structured polar candidates (certified served-space tests); K_10 structured intersections and
   candidates (incl. K3-surface double windows); all scripts are in the packet and queued (submit_queue.sh).
5. msolve pitfalls documented (c*c parsed as c; primes near 2^31 silently wrong); Singular over Q used for all
   candidate-locus statements.

## Ten-line summary (draft; final numbers filled in after the last local runs)
1. Best rigorous bound: every certificate for (F, 1, X0^3Y0^3) has deg D >= 7 (exact; all-parity degree-6 cone = Euler line only);
   deg D >= 6 is fully settled (exact), and no Y-pure degree-6 denominator coprime to F exists (exact).
2. Degree-6 proof: K_6 = J_6 (both parities, two primes, Hilbert-certified) + Hessian relaxation A d_w^2 F in J_12
   (kernels 18 / 36 / 261 for odd / even / mixed, exact lifts verified) + Singular over Q (cone = {0} / Euler line).
3. Degree-5 proof: exact integer nullspace of the full relaxation (pi_W(ker) = 0), and independently K_H = 0.
4. All wave-1 rigorous claims re-verified at a second prime or exactly; three wave-1 polar verdicts upgraded to rigorous
   via the explicit fixed-D kernel; dossier error corrected (E_I F in J for all I; E_X F tested: excluded).
5. Degree-7 Y-even polars: K_H has dim 1533 (certified); cone of 483 quadrics in 63 variables undecided locally;
   all 17 structured degree-7 polars excluded rigorously (A d_w D not in J_14, Hilbert-certified); degree 8: see hesstest8.
6. Not executed (remote saturated all session, load 300-1060): new-prime E1 kernels d = 7(odd)..12, 74 certified fixed-D
   tests of structured polars (d = 6..10), K_10/J_10 bases and structured intersections (quintics in Q, delta windows,
   K3-curve/K3-surface double windows, plane windows); scripts in the packet, submission queue still waiting.
7. Tool pitfalls documented: engine nullity off-by-one for inconsistent systems; msolve parses c*c as c and silently
   miscomputes for primes near 2^31; Singular minAssGTZ hilb overflow on the 36-variable cone (chart method instead).
8. Packet: /tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/ (README with exact statements, scripts, receipts).
9. Draft post: /tmp/hodge/notes/wave2-quadratic-cubic-pullback-draftpost.md.
10. No witness; no claim of nonexistence beyond the stated degree bounds.

### Result R3 (rigorous): 17 structured degree-7 Jacobian-type candidates excluded by the Hessian condition
Script wave2_excl_hesstest.py, p = 2147483629, degree-14 Y-odd Jacobian block rank 5052 = Hilbert (certified):
A d_w(D) not in J_14 for every w in {grad(X0 Q0), grad(X0 Q3), grad(X1 Q0), grad(Y0 Q0) [odd], grad(X0^3), grad(X0 Y0 Y1),
(Q0,Q1,Q2,0,0,0), (0,0,0,Q3,Q4,Q5) [odd], (Q1,Q2,Q0,0,0,0), (Q3,Q4,Q5,0,0,0), (0,0,0,Q0,Q1,Q2) [odd], Q_k e_k (k=0..5)}.
Hence none of D = w.grad F admits a certificate (rigorous; a failure of the necessary condition T5).  At degree 6 the
same test excludes E_X F and grad Q_k . grad F (k = 0..5) independently of R2.  Degree 8 running.

### Result R2' (exact): mixed-parity degree 6 — the whole 36-coordinate polar family
Joint Hessian relaxation (all parities): dim K_H = 261 (certified upper bound, p = 268435399; J_12 blocks 2223/2214 =
Hilbert); all 261 lifted vectors (entries in {0,+-1,+-2}) verified exactly in both parity blocks (exact solves,
224 s); exact annihilator: 405 quadrics in 36 variables; Singular over Q: dim(std(I)) = 1 (a finite union of lines
through 0; minAssGTZ crashed with an internal hilb overflow).  Chart L1 = sum r_a w_a = 1 (random r in [1,50]):
zero-dimensional of length 923 whose RADICAL is a single rational point, w5 = w10 = w15 = w20 = w25 = w30 = 1/85,
all other 30 coordinates 0 — exactly the Euler point z/6 (scaled).  HYPERPLANE1: dim(I + (L1)) = 0 (Singular over Q), so no line of the cone lies in {L1 = 0}: the Euler line is the ONLY line.
(Confirmation with a second random chart L2: again length 923 with radical a single point.)
**Conclusion: no degree-6 denominator of any parity is coprime to F.  Together with R1 and K_d = J_d (d <= 6): deg D >= 7 for every certificate.**
Chart argument (exact): the cone C = {w : w w^T in K_H(Q)} is a closed cone of affine dimension 1, hence a finite
union of lines through 0.  A line meets the affine chart {L1 = 1} unless it lies in the hyperplane {L1 = 0}; the chart
contains exactly one point (radical of length-1), the Euler point; and dim(I + (L1)) = 0 would show C ∩ {L1 = 0} = {0},
i.e. no line is hidden.  Together: C = the Euler line = {w = lambda z}, D = 6 lambda F.
