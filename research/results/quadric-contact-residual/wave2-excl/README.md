# quadric-contact-residual — wave 2 EXCLUSION packet

Instance: `quadric-contact-residual` (q = 1, A = X0^3 Y0^3). Prime used throughout: P2 = 2147483629 (= hodge.core.PRIMES[1]);
wave 1 used 2147483647, the dossier 1000000009, Singular runs 32003. All polynomial arithmetic is exact mod P2 (python-flint).
Variables X0, X1, X2, Y0, Y1, Y2; G_aff character of a monomial chi(e) = (e0+e3 mod 2, e2+e5 mod 2, e1-e4 mod 6, e3+e4+e5 mod 2).

## Exact statements proved (all transfer from F_P2 to C by the arguments in `wave2-quadric-contact-residual-excl.md`)
* R1 (replay): dim (J:A)_d = 0,0,0,0,0,6,36,126,337,760 for d = 0..9; the only non-Jacobian E1 directions below degree 10
  are D8* (d = 8, trivial sector) and Z_i*D8* (i = X1,Y1,X2,Y2) at d = 9.  [reverify.py part 1; log exp=E1-colon]
* R2 (replay): no certificate with deg D = 5 or 6.  [reverify.py parts 2-3; log exp=closure-d5, closure-d6]
* R3 (new): no certificate with deg D = 7 whose denominator's G_aff-character support fixes two of the three involution
  characters (sigma0-parity b0, sigma2-parity b2, tau_Y-parity y) — 12 subfamilies, msolve -g 1 leading ideal m-primary.
  [sub_closure.py on d7mixed.quads.local.json; log exp=sub-closure-msolve, tags d7_b00b20 ... d7_b21y1]
* R4 (new): no certificate with deg D = 8 and D Gamma-invariant (G_aff + s0 + r1), in particular none of the form
  c*D8* + (Gamma-invariant Jacobian part) + F*(invariant quadric).  [inv_closure.py 8; log exp=closure-Gamma-invariant]
* R5 (new): no certificate with deg D = 8 and D a G_aff eigenvector of any of the 20 characters listed in the log
  (exp=closure-pure-sector, d=8, result EXCLUDED); the four characters (0,0,0,0), (0,0,0,1), (1,1,0,0), (1,1,0,1) are pending.
* R6 (new): served-target spaces.  For every block D in served_local.json (linear forms, 11 quadrics, 7 quartics, polars F_Xi)
  the space of Y-odd degree-6 targets served by D equals J_6^odd (dim 18): no nontrivial class is supported on any of these
  divisors, and by the sum principle no product of them serves A.  [served_local.py; log exp=served, served-span]

* R7 (replays of the geometric inputs): fqlx's K3 socle test replayed over Q (Singular, /tmp/mstest/k3.sing copied here as k3.sing);
  the dossier's epsilon-sector dimensions (25/3/22, 191/140/51) replayed; the dossier's "(J:A)_d = S_d for d >= 13" corrected
  (dim R_19 = 246; dim (J:A)_13 = 8322 by duality).
* Second prime: d7mixed.quads.p1000000007.json and the tags d7p3_* / the local_d8 runs with "prime": 1000000007 in the log
  reproduce R3, R4, R5 at p = 1000000007.

* R8 (new, degree 7 cube lattice): all 28 pairs of (Z/2)^3 parity classes excluded (12 edges = item R3, 12 face diagonals,
  4 space diagonals; tags d7_fd_*, d7_sd_*, two primes).  Difference-multiplicity lemma (report item 18) + bilinear zero-pair
  tests (bilinear_pairs2.py, bidegree (2,2) Macaulay rank per h1-block; log exp=bilinear-pair) reduce the possible supports
  of a degree-7 solution to the 14 affine planes and the 6-, 7-, 8-sets.

## How to replay
```
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
python3 reverify.py 1 2 3                 # E1 colon dims, d=5, d=6 closures (Singular verdicts)
python3 d7root_local.py                   # builds d7mixed.quads.local.json (950 quadrics, 120 unknowns)
python3 sub_closure.py d7mixed.quads.local.json d7_b00b20 'b0=0&b2=0' 2     # one of the 12 subfamilies (msolve -g 1)
python3 inv_closure.py 8 2 600            # Gamma-invariant degree-8 family
python3 d8sectors_local.py 0020,0021 60   # pure sectors at d=8
python3 served_local.py 5                 # served-target spaces and the span test
python3 macaulay.py <quads.json> 4 8      # block Macaulay closure (I_4 = S_4 test) for a quads JSON
python3 sub_closure.py d7mixed.quads.local.json d7_sd_000_111 'b0=0&b2=0&y0|b0=1&b2=1&y1' 2   # a space diagonal (| = union of classes)
python3 bilinear_pairs2.py d7mixed.quads.local.json 2 2 all   # bilinear zero-pair tests for all 28 class pairs
```
Verdict rule for msolve: `-g 1` prints the grevlex leading ideal; the homogeneous quadric ideal has V(I) = {0} iff every
variable has a pure power among the leading monomials (msverdict.py).  NOTE: msolve misparses `c*c`; squares are written `c^2`.
Singular verdict: `dim(std(I)) == 0`.

## Files
w2lib.py (sectors, normal forms, sector fixed-D solver), family.py (E1-admissible families with explicit gauge, Hessian
relaxation), reverify.py, fixedD_batch.py, d7mixed.py / d7root_local.py, sub_closure.py, macaulay.py, branchkill.py,
inv_closure.py, d8sectors_local.py, parity_closure.py, extras_screen.py, served_local.py, msverdict.py.
Receipts: ../wave2-excl.log.jsonl (one JSON line per experiment, prime recorded), local_d8/*.ms (+ .lead), sub_*.ms (+ .lead),
inv/inv_d8.*, reverify/*, served_local.json.
