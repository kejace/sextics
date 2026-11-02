# wave-2 exclusion packet — klein-coupled-character (q = 1, A = X0^3 Y0^3)

Instance: `data/klein-coupled-character` (F = Klein sextic block K(X0,X1,X2) + three quartic-times-quadric couplings
+ Fermat block in Y; sigma = diag(zeta, zeta^2, zeta^4, 1, 1, 1), zeta^7 = 1; A has sigma-character 3).
Variables X0,X1,X2,Y0,Y1,Y2; "character" = sigma-weight mod 7 with weights (1,2,4,0,0,0); "parity" = Y-parity.

## Exact statements proved here (q = 1)

Let D be a nonzero homogeneous polynomial coprime to F for which a certificate (N, Q) exists, i.e.
A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q with deg N_i = deg D + 1.

T1 (E1-admissibility, re-verified).  After Euler normalization N(F) = A D, so A D in J(F).  The colon space
(J : A)_d equals J(F)_d for every d <= 8 in both Y-parities and all seven sigma-characters.  Proof: modular
computation at p = 1000003 (and, in wave 1, at 32003 and 2^31 - 1) of dim (J:A)_d mod p = dim J_d mod p, plus
the flatness check dim J_d mod p = complete-intersection value in every (parity, character) sector; since
(J:A)_Q reduces into (J:A)_p and dim J_d is the same over Q and mod p, dim_Q (J:A)_d <= dim_p (J:A)_d = dim J_d.
Receipts: jobA_p1000003/colon_sweep.json.

T2 (Hessian relaxation).  If D = sum_i w_i F_i (deg w_i = d - 5) then every certificate numerator has the form
N = A w + K(p) with p antisymmetric (Koszul), and E2 modulo J reads  A w^T Hess(F) w in J, i.e.
w^T Hess(F) w in (J : A)_{2d-6}.  The gauge w -> w + l z (D -> D + 6 l F, N -> N + l A z/6) preserves
solvability, so w may be taken in a complement of {l z}.  The resulting homogeneous quadratic system in the
w-coefficients has only the trivial solution in the cases listed in `relaxation-verdicts.md`; projective
emptiness mod p implies projective emptiness over Qbar (reduce a Qbar-point modulo a prime above p).

T3 (degree bound, rigorous). No certificate denominator of degree <= 6 exists, in any Y-parity (including mixed):
  d <= 4: (J:A)_d = 0.  d = 5: relaxation empty (Macaulay degree 2 full rank, local_d6/ and local_nzd/nzd.out).
  d = 6, Y-odd: relaxation empty (local_d6/d6_odd.*).  d = 6, Y-even and mixed: the exact relaxation ideal contains
  u_j^2 for every coordinate outside the weighted-Euler family D_lambda = sum_i lambda_i Z_i F_i (local_d6/rad_*.sing),
  and on that family the 8 exact constant-dual quadrics (local_d6/euler_quadrics_exact.json, duals in
  euler_duals_exact.json) cut out exactly the gauge line lambda = c(1,...,1) (D = 6F): local_d6/exact_duals.out.
T3' (linear factors). Every nonzero l in span(X0,X1,X2) or span(Y0,Y1,Y2) is a non-zero-divisor on R_6
  (local_nzd/*.sing: annihilator of the quadratic-syzygy module of the restricted partials is (1) in every chart);
  hence served(l D) = served(D) mod J (sum principle).  The relaxation locus C at d=6 is the complete intersection
  of the 3 quadrics in local_d6/C_quadrics_exact.json (weighted-Euler P^4).
T3'' (degree 7, sigma-eigen, pure parity). Relaxation empty for Y-even characters 0,3,5,6 and Y-odd characters 1..6
  (local_d7/*.lead, mod p => Qbar).  For Y-even characters 1,2,4 the relaxation locus equals the lifted curve X_j.C
  EXACTLY over Q (local_radical_exact/: msolve char 0 GB + Singular radical membership), hence excluded by T3+T3'.
  Y-odd character 0: lifted family span(Y).C (exact replay: see relaxation-verdicts.md).

T4 (specific denominators, exact rational duals).  For each D listed in `duals/`, a rational vector y with
y^T M(D) = 0 and y^T b(D) = 1 (verified over Q) proves that no certificate with denominator D exists; by the
prescribed-support theorem no certificate with poles inside V(D) exists.

## Replay
All scripts take the prime as first argument; run from the lane directory with the engine on PYTHONPATH:
`python3 wave2-excl/w2_jobA.py 1000003`, `python3 wave2-excl/w2_jobC.py 1000003 24`, `python3 wave2-excl/w2_jobD.py 1000003 16`,
`bash wave2-excl/w2_jobB.sh 48 36000 wave2-excl/jobA_msolve/d7_even_all.ms ...`.  Local (light) pieces: `local_d6/`.
Every experiment appends one JSON line to `../wave2-excl.log.jsonl` (fields exp, D/family, d, par, k, prime, result,
rank, seconds, job, note).

## Proof chain for a sigma-eigen denominator of degree 7 or 8 (used sector by sector)
Let D be a certificate denominator, sigma-eigen of character k and pure Y-parity, deg D = d in {7, 8}.
1. E1 => A D in J => D in (J:A)_d = J_d (T1) => D = sum_i w_i F_i with w in (S_{d-5})^6, sector-compatible.
2. Gauge: (D, N) -> (D + l'F, N + l'Az/6) preserves E1 and E2 mod F for any l' in S_{d-6} of the parity/character of D
   (ch. 1 denominator gauge, verified by direct expansion), so w may be taken in the gauge-fixed slice W'.
3. Relaxation (T2): w lies in V = V(I_relax) over Qbar.  Either V is empty (mod-p emptiness transfers to Qbar), or the
   EXACT (characteristic 0) radical-membership computation shows V ⊆ L = { m w_lambda : m in span(sector quadrics or
   linear forms), lambda in C }, i.e. D = m D_lambda with lambda in C ⊂ {lambda_0 = 0} and m of degree d-6.
4. Sum principle (ch. 3, Gysin description of N^1, (3,1)-part): served(m D_lambda) = served(m) + served(D_lambda) mod J.
   NZD lemmas: served(m) ∩ (character-3 classes) = 0 for every linear m in span(X) ∪ span(Y) (T3') and every Y-quadric m
   (local_nzd/nzd_yquadrics.out); monomial quadrics X_i X_j, X_i Y_j are products of two such linear forms.
5. Hence [A] is supported on V(D_lambda), i.e. a certificate with denominator D_lambda exists (prescribed-support
   theorem + Verify.add for the J-part) — contradicting T3 (every lambda with lambda_0 = 0, lambda != 0 is excluded by
   the exact constant-dual quadrics).  Therefore no such D exists.

## Hilbert-function sandwich (rigorous containment from a mod-p Groebner basis)
For the lifted family L (exact ideal I(L) computed over Q by elimination from its parametrization) and the relaxation
ideal I_Q (exact) with mod-p reduction I_p:  HF(S/I(L)) <= HF(S/rad I_Q) <= HF(S/I_Q) <= HF(S/I_p)  (first: L ⊆ V(I_Q);
last: a mod-p Macaulay matrix is the reduction of an integer one, so rank_Q >= rank_p).  If the Hilbert series of
S/in(I_p) (from the msolve leading ideal) and of S/I(L) differ by a polynomial — checked by divisibility of the difference
of first-Hilbert-series numerators by (1-t)^n — then rad(I_Q) = I(L) in all large degrees, hence V(I_Q) = L over Qbar.
Script: w2_hfsandwich.py; receipts local_hfsandwich/<tag>.out (log exp = hilbert-sandwich).  Cross-checked against the
characteristic-0 Groebner computation on the d=7 even sectors k = 1, 2, 4.
