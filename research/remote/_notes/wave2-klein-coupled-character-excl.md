# wave 2 / EXCLUSION lens — klein-coupled-character (q=1, A = X0^3 Y0^3)

STATUS: no witness; in progress (2026-11-01 ~20:20 local). RIGOROUS so far: deg D >= 7 for every certificate denominator (any Y-parity). Remote jobs A/C queued (remote load > 400, 60+ jobs);
the decisive low-degree work is being done locally (small systems).

Log: /tmp/hodge/results/klein-coupled-character/wave2-excl.log.jsonl
Packet: /tmp/hodge/results/klein-coupled-character/wave2-excl/ (w2lib.py, w2_job{A,B,C,D,E}.py, local_d6/, local_d7/)

## Rigorous results so far (q = 1)
Method ("Hessian relaxation", T2 of the packet README): E1 forces A D in J(F); for d <= 8, (J:A)_d = J_d (mod p
+ flatness => over Q), so D = sum_i w_i F_i; E2 modulo J then reads w^T Hess(F) w in (J:A)_{2d-6}, a homogeneous
quadratic system in the coefficients of w (gauge w -> w + l z removed). Projective emptiness mod p transfers to Qbar.
* d = 5 (both parities): relaxation variety EMPTY (Macaulay degree 2 full rank). => no degree-5 denominator.
* d = 6, Y-odd: relaxation variety EMPTY (msolve leading ideal irrelevant; Singular dim 0, vdim 57). => no Y-odd degree-6 denominator.
* d = 6, Y-even: relaxation variety = ONE irreducible curve C in P^16 (Singular minAssGTZ over GF(p) and over Q);
  C lies in the "weighted Euler" family D_lambda = sum_i lambda_i Z_i F_i (12 off-diagonal coordinates vanish on C),
  birational to a plane octic in (lambda_t : lambda_r : lambda_s). For this sigma-invariant family the Koszul-form
  system (59 unknowns, 294 rows in the character-0 sector) has 150 CONSTANT left-kernel vectors y (y^T M_i = 0 for all
  six blocks), giving 8 independent necessary quadrics y^T b(lambda) = 0; with the gauge lambda_x = 0 the quadrics have
  no common projective zero (Singular dim 0). Since C lies in lambda_x = 0, NO point of C carries a certificate.
  => no Y-even degree-6 denominator. EXACT over Q (local_d6/exact_duals.out): the 8 quadrics cut out exactly the
  gauge line {l0=..=l5} (D=6F) and nothing with l_x = 0; radical membership (u_j^2 in I) over Q for all 12 off-diagonal
  even coordinates and all 30 off-diagonal/odd coordinates of the mixed system => mixed-parity degree 6 excluded too.
* d = 5 mixed: EMPTY (Macaulay degree 2). d <= 4: (J:A)_d = J_d = 0.
* THEOREM: deg D >= 7 for every certificate denominator (any Y-parity, incl. mixed).
* Structural fact: relaxation solutions lift, V_{d+1} ⊇ S_1·V_d, so degree-7/8 relaxation varieties contain the
  lifted families l·D_lambda (l linear, lambda in C). Those are handled by the sum principle (ch. 3) once
  (J : l)_6 = J_6 for the linear forms involved.
* LEMMA (rigorous, local_nzd/): every nonzero l in span(X0,X1,X2) or span(Y0,Y1,Y2) is a non-zero-divisor on R_6
  (quadratic syzygies of the restricted partials: module annihilator = (1) in every chart, mod p => Qbar). Hence
  served(l D) = served(D) mod J for all D (sum principle): linear factors from these spans can be stripped.
* d = 7, sigma-eigen sectors (mod p, msolve/Singular; local_d7/): relaxation EMPTY (rigorous) for Y-even characters
  0,3,5,6 and Y-odd characters 1..6; NONEMPTY exactly where lifts exist: Y-even chars 1,2,4 (curves X0.C, X1.C, X2.C)
  and Y-odd char 0 (the 3-dimensional family span(Y).C). Mod-p radical membership: V(sector) subset of the lifted
  family (even k=1,2,4 done). Rigor caveat: containment mod p does not transfer to Q-bar (emptiness does), so the
  containments are replayed EXACTLY over Q (msolve char 0 + Singular; local_radical_exact/): even k=1,2,4 DONE over Q
  (exact GB size 94 each, every generator of I_L has a power in the exact ideal); odd k=0 (24 vars) running.
  Once exact: lifted points are l.D_lambda with lambda in C, excluded by the NZD lemma + sum principle + the d=6 result
  => no sigma-eigen denominator of degree 7 of pure Y-parity.
* METHOD (Hilbert-function sandwich, w2_hfsandwich.py): containment V(I_Q) = L (lifted family) follows rigorously
  from a mod-p Groebner basis when the Hilbert series of S/in(I_p) and of the exact ideal I(L) differ by a polynomial;
  cross-checked against the characteristic-0 Groebner bases on the d=7 even sectors k=1,2,4 (all agree: dim 2, degree 8).
  This replaces the (slower) characteristic-0 replays for the remaining sectors.
* LEMMA (rigorous, local_nzd/nzd_yquadrics.out): every Y-quadric q(Y) satisfies (J:q)_6 ∩ (Y-odd, char 3) = J_6^-[3]
  (needed for the degree-8 lifts q(Y).D_lambda in the even character-0 sector).
* Validation: the Koszul-form solver underlying the constant-dual proof reproduces 7 accepted witnesses (incl. two of
  degree 6 with 388 unknowns) and rejects perturbed denominators (local_d6/koszul_validation3.out).
* d = 7 full (non-eigen, 60 unknowns): msolve GB mod p running locally (local_d7full/); d = 8 eigen sectors running.

## Re-verification of wave 1 (so far)
* Hilbert function of S/J in every (parity, character) sector equals the complete-intersection value mod 1000003
  for d = 5..10 (flatness); kernel of [A] zero in all sectors for d <= 8; at d = 9: even [2,3,3,2,10,1,3],
  odd [0,5,0,0,9,6,4]; d = 10 even [24,46,39,39,54,46,46] — all agree with wave 1 (32003 and 2^31-1); NEW: d = 10 odd
  kernel [36,39,39,36,39,36,36] (total 261). Full fixed-D tests of the d = 9 kernel representatives and random
  combinations at p = 1000003 are running (jobC_local, 132 systems).
* The wave-1 "served-targets" entries are mislabeled (character 4, trivial elements X2 F_j); to be redone (job D).
* The wave-1 "scan-1param consistent" hit is the degenerate point 6F.

## Pending
jobs A (colon sweep d<=12, eigen d=7,8 relaxations), B (msolve full d=7/8), C (d=9,10 sectors + exact duals),
D (structured families, served targets, exact kernel reps), E (d=11,12 probes) — queued on the remote.
Local: exact Q proof of the d=6 even exclusion; eigen d=7 sectors; (J:Z_i)_6 checks.
