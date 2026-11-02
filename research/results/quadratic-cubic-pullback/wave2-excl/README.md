# Wave-2 EXCLUSION / RIGOR packet — lane `quadratic-cubic-pullback` (F = sum_{i=0}^5 Q_i^3, q = 1, A = X0^3 Y0^3)

Everything here is reproducible with the engine at `/home/kejace/devel/sextics-work` (python-flint 0.9) and the lane
helper `qcp.py` (quadrics Q_0..Q_5, delta = det(dQ/dZ)).  Variables X0,X1,X2,Y0,Y1,Y2; "Y-parity" = parity of the
total Y-degree; A is Y-odd, F is Y-even.

## Files
* `w2x.py` — helpers: parity-split Hilbert function of R = S/J(F); structural fixed-D kernel; rank-certified E1 test;
  rank-certified served-space computation (`served_certified`); candidate constructors (`polar`, `gradpolar`, `euler`).
* `wave2_excl_reverify.py` — re-verification of wave-1 structural claims at the prime 2147483629 and exactly.
* `wave2_excl_kernel.py` — E1 kernels K_d = {D : A D in J(F)} per parity and prime, Hilbert-certified; saves bases at d >= 10.
* `wave2_excl_relax.py` — exact (integer) relaxation of the polar family D = w.grad F at degree d (see statement S2).
* `wave2_excl_polar.py` — certified served-space tests of structured Jacobian-type denominators (degrees 6-10).
* `wave2_excl_k10.py` — degree-10 kernel: intersections with structured subspaces, candidate extraction, tests.
* `receipts/` — JSON lines copied from `../wave2-excl.log.jsonl` (every experiment: prime, ranks, Hilbert values,
  structural nullities, timings) and the kernel basis files.

## Rigor conventions
* R1 (modular rank).  A rank computed mod p is a lower bound for the rational rank.  A kernel/served dimension computed
  mod p is therefore an UPPER bound for the rational one *provided* the rank of the block that is subtracted is the same
  mod p and over Q.  We certify that block rank by (a) a Hilbert value (Jacobian blocks: dF is a regular sequence, so
  dim J_k^{par} = dim S_k^{par} - dim R_k^{par} with R's series (1-t^5)^3(1-ut^5)^3/((1-t)^3(1-ut)^3), u^2 = 1), or
  (b) an explicit rational kernel of the right dimension (fixed-D block: the elements (G z, 0) and (g F, D div g - g.grad D)).
* R2 (E1).  A D not in J(F) mod p, with the Jacobian block of degree d+6 at its Hilbert rank, implies A D not in J(F) over Q.
  At q = 1 this excludes V(D) ∩ X as a support at every pole order (prescribed-support theorem).
* R3 (exact).  The relaxation nullspaces are computed over Z (flint fmpz_mat.nullspace), no primes involved.

## Statements (see ../wave2-excl.log.jsonl and the per-job receipts in receipts/)
S1 (rigorous). K_d = J(F)_d for Y-even d <= 8 and Y-odd d <= 6 (prime 2147483647, wave 1) and for d = 3, 5, 6(odd)
    also at 2147483629 (wave 2); Jacobian blocks at Hilbert rank in every case.  [d = 7, 8 odd and d = 9: queued remotely]
    Consequence: every denominator of degree <= 6 is a polar D = w.grad F; every D of degree <= 6 not in J(F) is excluded
    at every pole order.
S2 (rigorous, exact). No certificate has deg D = 5: for D = w.grad F, w in C^6 (no parity assumption), the exact integer
    nullspace of the full relaxation (W = w w^T, V = c (x) w, Q) projects to zero on W.  Hence deg D >= 6.
S3 (rigorous, exact). No Y-pure denominator of degree 6: Hessian relaxation A (w^T Hess F w) in J_12, kernel dims 18 (odd)
    / 36 (even) as Hilbert-certified upper bounds, bases lifted and verified exactly, Singular over Q: the locus
    {w : w w^T in K_H} is {0} (odd) and exactly the Euler line w = lambda z, D = 6 lambda F (even), which is not coprime.
S3' (in progress). Mixed-parity degree 6: joint 36-coordinate relaxation (666 W-unknowns), then exact lift + Singular.
S4 (screened). Degree-7 Y-even polars: dim K_H = 1533 of 2016 (certified upper bound); candidate locus mod p pending.
S5 (queued remotely). Structured polars of degree 6..10 (74 candidates, certified served-space tests), K_10 bases and
    their intersections with the structured subspaces (quintics in Q, delta x quartics, K3-curve and K3-surface windows,
    plane windows), d = 11, 12 kernels.
Structural identities used (verified exactly): E = F (Phi(w) + L(c,w)) for N = A w + c.grad F;  H-terms cancel;
Phi(u, w2) = (1/2) F (grad A . w2 + A div w2) for u = z/6 (so u (.) w2 are kernel directions for every w2: the trivial
solution D = F and its polarizations); the fixed-D kernel elements (G z, 0) and (g F, D div g - g.grad D).

## How to replay
    cd /tmp/hodge/results/quadratic-cubic-pullback
    python3 wave2_excl_reverify.py
    python3 wave2_excl_kernel.py replay 4 "9:0:1,9:1:1,10:0:1:b,10:1:1:b"
    python3 wave2_excl_relax.py replay5 5 all        # exact, ~15 s
    python3 wave2_excl_relax.py replay6 6 0          # exact, minutes
    python3 wave2_excl_polar.py replayP 8 1 9        # structured polars up to degree 9
    python3 wave2_excl_k10.py replayK 8 0 1 10 6     # needs the d=10 basis file from the kernel run
