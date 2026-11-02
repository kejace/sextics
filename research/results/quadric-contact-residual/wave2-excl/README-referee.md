# quadric-contact-residual — wave-2 exclusion packet: referee replay guide

Written by the referee pass (Claude, Anthropic, for `kejace`; no human review). It tells a stranger how to replay every
rigorous claim of this packet with pinned tools and says which claims the referee replayed, with which method, and what
transfers to characteristic 0.  The claim-by-claim verdict table is `referee-report.md` (same directory).

## Pinned engine
* python 3 with **python-flint 0.9.0** (`import flint; flint.__version__`), **msolve 0.6.5** (`msolve -h`), **Singular 4.3.2**
  (`Singular --version`), `pyyaml`.  The polynomial dict arithmetic and the instance loader come from the project engine
  `hodge.core` (`/home/kejace/devel/sextics-work`, mirrored at `research/engine/` of the fork branch); `core.load_instance` reads
  `data/quadric-contact-residual/instance.rec` of the main repository (`F`, `A = X0^3 Y0^3`, `q = 1`).
* Primes.  Lens runs: P2 = 2147483629 (main) and 1000000007 (confirmation); wave 1: 2147483647; dossier: 1000000009;
  aggregator cross-checks: 998244353; **referee: 1000000021** (unused before).  Exact rational parts use `fmpq_mat.rref`.
* Variables `X0, X1, X2, Y0, Y1, Y2`; G_aff character of a monomial `chi(e) = (e0+e3 mod 2, e2+e5 mod 2, e1-e4 mod 6, e3+e4+e5 mod 2)`
  (sigma0, sigma2, h1, tau_Y of fqlx); parity class = (chi_1, chi_2, chi_4) in (Z/2)^3; `A` has character (0,0,0,1).

## Referee code (independent re-implementation): `referee/ref.py`
```
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl/referee      # or research/results/quadric-contact-residual/wave2-excl/referee
python3 ref.py hilbert            # dim J_m mod p = CI Hilbert function, m = 5..18, at 1000000021, 2147483629, 1000000007, 998244353
python3 ref.py e1                 # dim (J:A)_d, d = 0..9, per G_aff sector (expect 0,0,0,0,0,6,36,126,337,760)
python3 ref.py d5exact            # EXACT over Q: the 21 columns A*d_i d_j F have rank 21 modulo J_10   (=> no degree-5 denominator)
python3 ref.py d8star             # exact identities N*(F) = A D8*, D8* (X0^2-Y0^2) = F_X0^2 - F_Y0^2; D8* not in J_8
python3 ref.py d6                 # degree-6 Hessian relaxation (35 unknowns): msolve -g 1 m-primary AND Singular dim 0
python3 ref.py d7pair 000 001     # one parity-class pair at degree 7 (writes d7_pair_000_001.{ms,sing,lead}); any of the 28 pairs
python3 ref.py d8sector 0110,0021 # pure G_aff sectors at degree 8 (gauge-fixed E1-admissible family, msolve -g 1)
python3 ref.py d8inv              # Gamma-invariant degree-8 family: 3 gauge + 16 Jacobian + D8*, 17 unknowns, 52 quadrics
python3 ref.py fixedD-exact 'D8*' # EXACT over Q: sector-restricted fixed-D system for A with D = D8* (rank M vs rank [M|b])
python3 ref.py served-exact S0 delta GS F_X0 ...   # EXACT over Q: served(D) for all 24 Y-odd target characters (eigenvector D)
python3 ref.py served-exact-full S1                # EXACT over Q, full (unsplit) system, for blocks that are not eigenvectors
```
Every command appends one JSON line per result to `referee/referee.log.jsonl` (prime, sizes, ranks, verdicts, seconds);
msolve/Singular inputs and the msolve leading ideals (`*.lead`) are kept next to it.  Set `REF_P=<prime>` to change the prime.

## How each kind of claim is checked and why it transfers to C
1. **Hessian relaxation closures (degrees 5-9).**  For `D = sum c_k D_k` in the E1 kernel with E1 lifts `N_k`, every certificate
   forces `B_D(N) in J_{2d}` (fqlx #125 §1), a homogeneous quadric system `I` in `c`.  `msolve -g 1` prints the grevlex leading
   ideal; a pure power of every variable there means `F_p[c]/I` is finite, i.e. `V(I) = {0}` over the algebraic closure of
   `F_p` (Singular `dim(std(I)) == 0` is the second method).  Transfer: the mod-`p` quadrics are the reductions of the rational
   ones because (i) the family and its lifts are rational (`m F_j`, `A m e_j`, `D8*`, `N*`) or any mod-`p` E1 lift differs from
   a rational one by a mod-`p` syzygy, which is a reduction of a rational (Koszul) syzygy when `dim J_{d+6} mod p` has its
   rational value, and (ii) `dim J_{2d} mod p` has its rational value, so `J_{2d} mod p` is the reduction of the saturated
   lattice `J_{2d} ∩ Z[Z]`.  Both are the `hilbert` check (all primes used, `m = 5..18`).  Then `Proj(Z[c]/I_Z) -> Spec Z`
   is projective with empty fibre at `p`, hence empty generic fibre: only `c = 0` over `C`.  Gauge: `D -> D + sF`,
   `N -> N + A s Z/6` acts trivially on the quadrics (uses `A D in J`), so a complement of the gauge subspace may be fixed;
   the dropped coordinates are chosen so that the gauge basis projects isomorphically onto them.
2. **Modular rank = lower bound** for the rational rank: "columns independent mod p" (the degree-5 columns, the E1 kernel
   dimensions `dim (J:A)_d mod p >= dim over Q`, independence of `Z_i D8*` mod `J_9`) transfers; for the E1 kernels the
   exact rational lower bounds `J_d`, `J_8 + <D8*>`, `J_9 + <Z_i D8*>` meet the modular upper bounds, so the values
   `0,0,0,0,0,6,36,126,337,760` are exact.
3. **Fixed-D inconsistency / served spaces mod p do NOT transfer by rank monotonicity alone** (`M = [p], b = [1]` is consistent
   over Q and inconsistent mod p).  They transfer when `rank_p(M) = rank_Q(M)`, which a modular run cannot certify.  The
   referee therefore recomputed these **exactly over Q** (sector-restricted; restriction to G_aff isotypic components is
   lossless by averaging over the group of order 48): `fixedD-exact` and `served-exact*`.
4. **Sum principle and powers** are theorems of the project book (`docs/01-certificates.md`, "Prescribed support and pole
   multiplicity" and "Combining certificates"): at `q = 1` a certificate with denominator `D` exists iff `Res(A Omega/F^2)`
   vanishes on `X \ V(D)`, so served(D) depends only on `V(D)` (served(D^m) = served(D)); for `V(D1 D2) = Z1 ∪ Z2` the kernel
   of restriction is the Gysin image of a resolution of `Z1 ⊔ Z2` (weight argument), so served(D1 D2) = served(D1) + served(D2).
   Hence "served(D_i) = J_6^odd for every block" implies that no product of the blocks serves any nontrivial target.
