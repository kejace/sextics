# Wave-2 exclusion packet — `four-plus-two-residual` (q = 1, A = X0^3 Y0^3)

Instance: `data/four-plus-two-residual/instance.rec`, F = H(X0,X1,Y0,Y1) + B(X2,Y2) (variables X0,X1,X2,Y0,Y1,Y2).
Certificate (cleared form, q = 1): A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q, deg N_i = deg D + 1, deg Q = 2 deg D - 6.

## Contents
* `w2lib.py` — exact/modular linear algebra: Jacobian-ring Hilbert functions, multiplication kernels (exact over Q and mod p),
  small-height verified kernel representatives, the fixed-D system builder (optionally parity-free or character-sector
  restricted), the exact Koszul-gauge nullity, served-target sectors.
* `w2_kernels.py` — E1-kernel structure: regular-sequence checks, core kernels k = 0..16 (3 primes + exact), direct 6-variable kernels d <= 8 vs Thom–Sebastiani.
* `w2_d5.py` — degree 5: the Hessian relaxation (exact over Q) + rank-one test (Singular), partials and random J(F)_5 members.
* `w2_d6.py` — degree 6: full unrestricted system, 2 primes, exact [RANK] nullity, 33 kernel reps, structural D, random family members (pure and mixed parity), pencils.
* `w2_symsweep.py` — core-only kernel reps k = 7..12, full 6-variable system restricted to the character sector of the diagonal group (lossless), exact ranks over Q.
* `w2_sum.py` — served-target sectors and the q = 1 sum principle (products of core and binary factors).
* `w2_generic.py` — random members of the whole E1-admissible family at degrees 7..10.
* `w2_pull.py` — pulls the remote mirror and merges the per-job JSONL receipts into `../wave2-excl.log.jsonl`.
* `receipts/` — the per-job JSONL logs and JSON receipts copied from the remote (primes, ranks, dimensions, timings).

## Replay
    cd /tmp/hodge/results/four-plus-two-residual && python3 wave2-excl/w2_kernels.py excl-kernels
    python3 wave2-excl/w2_d5.py excl-d5
    python3 wave2-excl/w2_d6.py excl-d6 16            # 16 worker processes
    python3 wave2-excl/w2_symsweep.py excl-symsweep 32 7 12
    python3 wave2-excl/w2_sum.py excl-sum 24
    python3 wave2-excl/w2_generic.py excl-generic 8 7 10 2
(engine on PYTHONPATH: /home/kejace/devel/sextics-work; python-flint >= 0.9, sympy, numpy, scipy; Singular for the rank-one test.)
Each script appends one JSON line per experiment to `remote/<job>.jsonl`.

## Statements and their rigor labels
Labels: [EXACT] exact over Q (hence over C); [RANK] modular rank + exact nullity (rigorous over C for the specific
rational D); [MODP] modulo one prime (evidence).  Full report: `/tmp/hodge/notes/wave2-four-plus-two-residual-excl.md`;
draft discussion comment: `/tmp/hodge/notes/wave2-four-plus-two-residual-draftpost.md`.

T1 [EXACT] (receipts `excl-kernels*.jsonl`, `excl-kernels-corereps.json`).  H and B have regular Jacobian sequences over Q and
    mod 2^31-1, 10^9+7, 998244353.  ker(a. on R(H)_k) modulo J(H)_k, (Y0,Y1)-even/odd: k<=5: 0/0; 6: 20/13; 7: 30/30; 8: 41/34;
    9: 38/38; 10: 36/31; 11: 26/26; 12: 19/16 (= all of R(H)_k for k >= 11).  Consequently {D : A D in J(F)}_d =
    sum_k ker_{d-k}(H) (x) C[X2,Y2]_k + J(F)_d; no E1-admissible D of degree <= 4; degree 5: D in J(F)_5 only.
T2 [EXACT] (receipts `excl-d5-local*`, `excl-d5-control-local.jsonl`).  No certificate with deg D = 5 exists over C
    (Hessian relaxation of the reduced E2 for D = v.grad F; W of dim 8 with no rank-one point).  With T1: deg D >= 6.
T3 [EXACT] (receipts `excl-d6jac-local*`, `excl-d6jac2-local*`).  No G-eigenvector Jacobian-type denominator of degree 6
    (13 character classes, 10+8+11x2 parameters) admits a certificate; the only rank-one points are the Euler line D ~ F.
T4 [RANK]/[EXACT] (receipts `excl-symsweep-local-*.jsonl`).  Individual core-only kernel representatives (exact rational,
    verified) of degrees 6..12 in the sector-restricted full 6-variable system -- see the log for the per-D verdicts.
T5 [MODP] (receipts `excl-sum-local.jsonl`, `-served.json`).  Served-target sectors and sum-principle membership tests.
T6 [RANK]/[MODP] (receipts `excl-d6.jsonl`, `excl-generic.jsonl`, `excl-struct-local.jsonl`, when available).
