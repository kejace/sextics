# wave2 / one-reflection-skew / lens EXCLUSION-RIGOR

STATUS: in progress (no witness; this lens does not search for witnesses, it certifies exclusions)

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
