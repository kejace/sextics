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

## Statements (filled in as the jobs finish; see ../wave2-excl.log.jsonl)
S1. K_d = J(F)_d for d <= [[9]] in both parity sectors (primes 2147483629 and 2147483587; wave 1: 2147483647).
    Hence every certificate has deg D >= 10 or D = w.grad F.
S2. No certificate has deg D = 5: for the full family D = w.grad F, w in C^6, the exact relaxation has pi_W(ker) = 0.
    With K_d = 0 for d <= 4: deg D >= 6 for every certificate.
S3. Degree 6 polar sectors: [[relaxation result]].
S4. Structured polars of degree 6..10: [[table]] (served dimension mod J, certified or modular).
S5. Degree 10: dim K_10/J_10 = [[294 even / 261 odd]]; intersections with structured subspaces: [[...]].

## How to replay
    cd /tmp/hodge/results/quadratic-cubic-pullback
    python3 wave2_excl_reverify.py
    python3 wave2_excl_kernel.py replay 4 "9:0:1,9:1:1,10:0:1:b,10:1:1:b"
    python3 wave2_excl_relax.py replay5 5 all        # exact, ~15 s
    python3 wave2_excl_relax.py replay6 6 0          # exact, minutes
    python3 wave2_excl_polar.py replayP 8 1 9        # structured polars up to degree 9
    python3 wave2_excl_k10.py replayK 8 0 1 10 6     # needs the d=10 basis file from the kernel run
