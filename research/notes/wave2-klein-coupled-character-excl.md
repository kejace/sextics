# wave 2 / EXCLUSION lens — klein-coupled-character (q=1, A = X0^3 Y0^3)

STATUS: in progress (no witness; wave-2 exclusion program started 2026-11-01 19:25 local)

Log: /tmp/hodge/results/klein-coupled-character/wave2-excl.log.jsonl
Packet: /tmp/hodge/results/klein-coupled-character/wave2-excl/ (w2lib.py, w2_jobA.py, w2_jobB.sh, w2_jobC.py, receipts jobA_*/jobC_*)

## Plan
(a) re-verify wave-1 rigorous claims with prime 1000003 and new code (colon ideal (J:A)_d per parity/character,
    Hilbert-function flatness check, Jhat and polar exclusions as exact rational duals);
(b) rigorous lower bound on deg D via the E2-mod-J (Hessian) relaxation for D in J_d (d <= 8), msolve/Macaulay;
(c) d = 9, 10 sector tests (all kernel reps, random combos) + structured families, new prime;
(d) packet; (e) draft post.

## Results so far
(none yet — jobs A and C submitted to the remote)

## Findings from re-reading wave 1 (before any new computation)
* The wave-1 "served-targets" log entries (family served-targets) are mislabeled: the listed served polynomials
  have sigma-character 4 (not 3) and are exactly X2*F_3, X2*F_4, X2*F_5 (trivial elements of J), i.e. the
  computation ran in the wrong sector (target_weight defaulted to the weight of the first odd monomial). To be redone.
* The wave-1 "scan-1param consistent=1" entry is the degenerate point radX + radY = 6F (F | D), not a witness.
