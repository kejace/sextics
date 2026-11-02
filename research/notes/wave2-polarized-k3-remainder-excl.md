# Wave 2 — `polarized-k3-remainder` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; this lens only produces exclusions / re-verifications)

Started 2026-11-01 ~23:59 EST. Lane dir `/tmp/hodge/results/polarized-k3-remainder/`, packet dir `wave2-excl/`,
log `wave2-excl.log.jsonl` (one JSON line per experiment).

## Plan
(a) adversarial re-verification of wave-1 / thread claims (second prime, independent implementation);
(b) rigorous lower bound on deg D: replicate d<=7 exclusion for A; NEW: same for the K3-sector target
    A_K3 = (U0 U1 U2)^2 - (V0 V1 V2)^2 (needed to make the sum-principle reduction to irreducible supports valid);
(c) non-Jacobian E1-kernel of [A_K3] at d=8,9 + structured irreducible families (fixed-D, remote, mod p);
(d) packet wave2-excl/ with receipts; (e) draft post.

## Results so far
(none yet)

### R1. Degree-6 exclusion for A independently REPLICATED (exact over Q)  [DONE 00:40]
Own implementation `wave2-excl/relax.py` ((U,V) frame, G = 7f(U)+7f(V), T = prod(U_i^2-V_i^2), L = (G)+J(G)^2+(R_ijk),
Euler gauge drops the z5-divisible coefficient of a_5). Prime 65521 and exact rationals. Results identical to the
thread's DEGREE_SIX.md (p = 101): block L-ranks 603/558/396/260/396/558 (mod p = over Q), augmented Q-ranks
713/715/617/470/617/715, constraint rank 437; exact square elimination kills the 30 off-diagonal entries, then the 5
diagonal ones (rank 15) => D = 0 or D = c*G. Degree-5 smoke test: all 6 parameters killed (rank 21). Receipts
`wave2-excl/relax-full-d6-p65521-local.json`, `-exactQ.json`, `smoke-d5-p101.json`.
