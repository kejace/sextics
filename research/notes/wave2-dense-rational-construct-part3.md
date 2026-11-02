# wave2 dense-rational CONSTRUCT — part 3: computational results (updated as jobs finish)

All lines are in /tmp/hodge/results/dense-rational/wave2-construct.log.jsonl; inputs/outputs in
/tmp/hodge/results/dense-rational/wave2-construct/ (msolve .ms inputs, .g1.out leading ideals, .ltcheck.txt receipts).

## 3.1 Projective loci of the Hessian and the reflection condition (G = 9973 F, p = 1073741789, msolve -g 1, DRL)
| locus | equations | degree | leading ideal | pure powers | wall time / RSS | verdict |
| V1 = {rank Hess G <= 1} | 120 (2x2 minors) | 8 | 1923 monomials, max deg 11 | all six (deg 8..11) | 4.9 s / 55 MB | EMPTY over C [rigorous] |
| V2 = {rank Hess G <= 2} | 210 (3x3 minors) | 12 | 10053 monomials, max deg 16 | all six (deg 12..16) | 374 s / 1.2 GB | EMPTY over C [rigorous] |
| V3 = {D^3G(e) vanishes on ker dG(e)} | 56 (T'_{ijk}) | 15 | running (job construct-V3) | | | pending |
Sanity checks on the Fermat sextic with the same generators: V1 gives the 15 binomials e_i^4 e_j^4 (six coordinate
points, no pure power), V3 is non-empty (525 leading monomials, no pure power) — as required, since Fermat has all
exchange covers and is even in every coordinate.  Log lines sanity-fermat-V1, sanity-fermat-V3.

Meaning (part 2): V1 = V2 = ∅ => (i) no holdout-type exchange cover in any linear coordinates, any block split (r = 1..4),
any base form g; (ii) no quasi-reflection of order 3, 5 or 6 in any coordinates (order 4 impossible on a smooth sextic).
V3 = ∅ would add: no linear reflection (F even in no linear form), closing the "hidden reflection" question of the dossier
(§1.2 "not done") for the one case that matters (reflections are the only involutions with age < 1).

## 3.2 Ternary cones D = g(X), g(Y): E1 test  (A g ∈ J_{k+6}, parity-split Hilbert caps, p = 2147483647)
| k | sector | m = k+6 / parity | cols | gens J | rank J | cap | #cones | E1-admissible cones | rigorous | secs |
| 7 | X | 13 / odd | 4284 | 3861 | 3441 | 3441 | 36 | 0 | yes | 30+26 |
| 7 | Y | 13 / even | 4284 | 3861 | 3441 | 3441 | 36 | 0 | yes | 36+49 |
| 8 | X | 14 / odd | 5796 | 6006 | 5052 | 5052 | 45 | 0 | yes | 139+130 |
| 8 | Y | 14 / odd | 5796 | 6006 | 5052 | 5052 | 45 | 0 | yes | (local) |
| 9 | X | 15 / odd | 7752 | 9009 | 7129 | 7129 | 55 | 0 | yes | (local) |
| 9 | Y | 15 / even | 7752 | 9009 | 7129 | 7129 | 55 | 0 | yes | (local) |
| 10 | X | 16 / odd | 10152 | 13104 | 9684 | 9684 | 66 | 0 | yes | (local) |
| 10 | Y | 16 / odd | 10152 | 13104 | 9684 | 9684 | 66 | 0 | yes | (local) |
| 11, 12 | X, Y | | | | | | | remote job construct-dense-rational-k1012-univ (queued; remote at capacity) | | pending |
(k <= 6: cones are in J_k ∩ C[X]_k = 0 by (J:A)_k = J_k and the plane-sextic Jacobian ring; F60, F06 ∉ J_6 is referee claim 23.)

## 3.3 Universal degree-12 denominators (R' and the products of the pieces)  — remote job construct-dense-rational-univ12
pending.
