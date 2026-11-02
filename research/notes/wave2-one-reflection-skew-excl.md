# wave2 / one-reflection-skew / lens EXCLUSION-RIGOR

STATUS: no witness; in progress (local job5_degree.py 7 [Y2-free strict pipeline] and remote jobs excl-ors-d7sub / excl-ors-families-f3c / excl-ors-extradirs still running at hand-off; their receipts land in wave2-excl/ and remote-*.log.jsonl)
Best rigorous bound: deg D >= 7 (degrees <= 6 re-verified exactly over Q by an independent implementation); at degree 7 all Y2-odd mixtures excluded rigorously (independent of the thread); Y2-free/mixed even cases pending.

Started 2026-11-01 ~19:00 local. Log: /tmp/hodge/results/one-reflection-skew/wave2-excl.log.jsonl
Packet: /tmp/hodge/results/one-reflection-skew/wave2-excl/

### (b) degree 7 (in progress)
* Stage A (exact, local): gauge-fixed E1 kernel K_7/F.S_1 = 129 = 120 (Jacobian) + 9 (non-Jacobian), per fine sector
  (24 characters of the diagonal group Z/6 x Z/2 x Z/2); agrees with the thread (135 = 126 + 9 before the 6 F-gauges).
  Non-Jacobian directions live in sectors (0,0,1):1, (3,0,4):7, (3,3,1):1 (all Y2-even; 2-torsion classes 0010, 1110, 1000).
  Sector table (fine char -> e5 mod 6, 2-torsion class, dim K): wave2-excl/d7_sector_table.json, d7_stageB.json.
* Stage B (mod-J Hessian relaxation, exact duals): 1019 quadrics, rank 748 on 5754 nonzero pair-monomials; NO square in
  the span (so no coordinate is killed at this stage, unlike d=5,6); 56 off-diagonal pair-monomials in the span
  (all between different fine sectors; involve 35 Jacobian + all 9 extra directions).
* Stage C (relaxed mod-F with all-character Koszul corrections on the full 129-dim family): running locally.
* Remote job excl-ors-d7sub: relaxed quadrics + msolve zero-dimensionality for the sub-families odd (all Y2-odd mixtures),
  tor:0010, tor:0100, tor:1000, tor:1110 (pure 2-torsion classes), even (all Y2-even mixtures), full.

### Y2-free reduction (proved here; used for the t-free sectors)
Lemma. Let D != 0 be Y2-free homogeneous of degree d and A Y2-free. A certificate (D,N,Q) for (F = f + Y2^6, A, q=1)
exists iff there is a Y2-free N' in Q[x]^5 of degree d+1 with  N'(f) = A D  and  D div N' - N'(D) = 0.
Proof sketch. (<=) N = (N',0), Q = 0 satisfies the cleared identity. (=>) average the certificate over mu6 acting on
Y2 (F, A, D invariant): N_i (i<5) become polynomials in x and s = Y2^6, N_5 = Y2 a(x,s), Q = Q(x,s). Euler-normalize
(N(F) = A D exactly; this keeps the shape). Both E1 and E2 (B_D(N) = (f+s) Q) are then polynomial identities in s.
Specialize s -> -f: with N~_i = N_i(x,-f), a~ = a(x,-f) one gets N~(f) = A D + 6 f a~ and (using d/ds of E1 at s=-f)
B_D(N~) = 5 a~ D.  Finally N' = N~ - a~ E' (E' = five-variable Euler field) gives N'(f) = A D and B_D(N') = 0.
Consequence: Y2-free denominators are tested by a strict 5-variable system with NO quotient Q and NO F-multiples;
the corresponding relaxation (5-variable Koszul corrections only) is much tighter than the 6-variable one.
(The strict form is the one fringewidth stated; the dossier could not reproduce it - the derivation above does.)

## Plan
(a) re-verify rigorous claims of dossier + wave-1 (different primes / methods, remote)
(b) push rigorous lower bound on deg D (currently deg D >= 7 from the thread; degree 7 open in even mixed sectors)
(c) exhaustive test of non-Jacobian E1-kernel elements and structured families, symmetry-restricted, full E1+E2
(d) packet, (e) draft post

## Results so far
### (a) re-verification of the thread's rigorous claims (independent exact implementation, python-flint, F scaled to 7F)
* deg <= 4: E1 kernel {D : A D in J} is 0 (exact nullspace, sector by sector). deg 5: kernel = J_5 (dim 6), exact.
  Hessian relaxation mod J: rank 15 on 15 nonzero pair-monomials (21 total), all six squares in the Q-span (certified
  exactly) => a = 0 => no degree-5 denominator.  Agrees with rbajaj5 C1 (rank 15; six square witnesses).
* deg 6: gauge-fixed kernel = J_6/F (35 params). mod-J stage kills 28 then 2 squares (exact certificates), 5 diagonal
  survivors (sector (0,0,0)), weak Hessian rank 4; relaxed mod-F stage with all 315 Koszul corrections (exact duals):
  9 independent quadrics; Macaulay matrix mod p: deg 3 rank 35/35, deg 4 70/70 => V_C = {0} => degree 6 excluded.
  Agrees with rbajaj5 C2 (nine quadrics, 35/35, 70/70).  Status: deg D >= 7 RIGOROUS (re-verified).
* wave-1 log: 318 experiments, all 5-variable fixed-D screens mod one prime (P1 = 2^31-1), engine 2; no rigorous
  claims beyond E1-kernel dimensions (5-var: d=7 par0 47=39+8, par1 37=36+1 -> total 84 = 75 + 9, consistent with the
  thread's (J_f:A)_7 = 84).  Nothing to retract; those screens are evidence, not proofs (single prime, no dual certificate).
* E1-kernel dimension table of the dossier (mod 32003 there) re-computed EXACTLY over Q (fmpz nullspace, sector by
  sector; remote job excl-ors-e1dims; receipts wave2-excl/e1dims_d<d>.json, bases e1basis_d<d>.pkl):
  d=5..14: dim K_d = 6, 36, 135, 480, 1144, 2196, 3711, 5762, 8322, 11502; dim J_d = 6, 36, 126, 336, 756, 1497, 2682,
  4437, 6882, 10122; non-Jacobian = 0, 0, 9, 144, 388, 699, 1029, 1325, 1440, 1380.  All agree with the dossier.
* Degree-7 E1 form of fqlx's addendum (D = h + t l.grad f + t^2 c.grad f + 6 t^5 R2, 84+25+5+15 = 129 after the six
  F-gauges) confirmed exactly by the Y2-degree decomposition of my gauge-fixed kernel.
* Degree-7 Y2-odd mixtures (all 40 odd directions at once): mod-J relaxation kills 23 then 12 squares (certified),
  5 survivors; strict relaxed stage: 9 quadrics, Macaulay 35/35 mod p => EXCLUDED.  Independent confirmation of the
  addendum's 'every degree-7 Y2-odd denominator is excluded, mixtures included' (remote job excl-ors-d7sub, case odd).
  The five odd survivors are exactly Y2 * (the five diagonal degree-6 survivors), i.e. the rbajaj5 sector-0001 reduction
  D = Y2 * sum a_i Z_i dF_i reappears, and the same nine quadrics / 35-35 Macaulay argument closes it.
* Degree-7 pure even 2-torsion classes 0010/0100/1000 (the addendum's Groebner claims, not replayed by the dossier):
  my 6-variable relaxation leaves exactly one coordinate axis per class (Y1*Y2^6, Y0*Y2^6, X1*Y2^6 ~ the f*Z_k directions),
  with multiplicity 495/67/419; each axis is EXCLUDED-exact by the full fixed-D test.  The 6-variable relaxation cannot
  close these classes alone (see (b)); the Y2-free reduction below does the job on the t-free part.
