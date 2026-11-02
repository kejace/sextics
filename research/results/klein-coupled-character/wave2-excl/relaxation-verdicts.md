# Relaxation verdicts (Hessian relaxation w^T Hess(F) w in (J:A)_{2d-6}, gauge-fixed), prime 1000003 unless "exact"

Legend: EMPTY = projective variety empty mod p (=> empty over Qbar; rigorous exclusion of the whole class);
LIFTED = equals the lifted family S_{d-6}.C (C = degree-6 weighted-Euler curve), proof status as indicated;
"exact" = characteristic-0 Groebner basis (msolve char 0 + Singular) for containment statements.

| d | parity | characters | unknowns | eqs | verdict | status |
|---|---|---|---|---|---|---|
| <=4 | any | all | - | - | (J:A)_d = 0 | rigorous (mod p + flatness) |
| 5 | even / odd / mixed | all | 3 / 3 / 6 | 12 / 6 / 30 | EMPTY (Macaulay deg 2) | rigorous |
| 6 | odd | all | 18 | 153 | EMPTY (GB; dim 0, vdim 57) | rigorous |
| 6 | even | all | 17 | 105 | curve C in weighted-Euler family; all points excluded by 8 exact constant-dual quadrics | rigorous (exact) |
| 6 | mixed | all | 35 | 313 | contained in the even curve C (exact radical membership) | rigorous (exact) |
| 7 | even | 0 | 3 | 15 | EMPTY | rigorous |
| 7 | even | 1, 2, 4 | 11 | 27/27/39 | LIFTED (X0.C, X1.C, X2.C) | rigorous (exact GB, size 94) + NZD lemma + sum principle |
| 7 | even | 3, 5, 6 | 8 | 53/51/51 | EMPTY | rigorous |
| 7 | odd | 0 | 24 | 59 | LIFTED (span(Y).C expected) | mod p / exact: running |
| 7 | odd | 1..6 | 6 | 24-31 | EMPTY | rigorous |
| 7 | even / odd | all (non-eigen) | 60 / 60 | 447 / 420 | GB running | pending |
| 7 | mixed | all | 120 | - | not attempted (infeasible GB) | open |
| 8 | even / odd | eigen sectors | 18-45 | - | running | pending |
| 8 | any | non-eigen | 156+ | 468 | infeasible (generic vanishing degree >= 12) | open |

Lemmas used for lifted components: (NZD-lin) every nonzero l in span(X) or span(Y) is a NZD on R_6;
(NZD-Yquad) every Y-quadric q satisfies (J:q)_6 ∩ (Y-odd, char 3) = J_6^-[3]; sum principle served(D1 D2) = served(D1)+served(D2) mod J.
