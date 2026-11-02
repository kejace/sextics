# wave2 / one-reflection-skew / lens EXCLUSION-RIGOR

STATUS: in progress (no witness; this lens does not search for witnesses, it certifies exclusions)

Started 2026-11-01 ~19:00 local. Log: /tmp/hodge/results/one-reflection-skew/wave2-excl.log.jsonl
Packet: /tmp/hodge/results/one-reflection-skew/wave2-excl/

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
