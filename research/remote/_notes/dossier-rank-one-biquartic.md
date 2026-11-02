# Research dossier: `rank-one-biquartic` (OPEN)

Prepared 2026-11-01 from the repository at `/home/kejace/devel/sextics` (branch `main`, HEAD `83fbf5c`), the saved
discussion thread `/tmp/hodge/discussions/rank-one-biquartic.md` (GitHub discussion #118), the reproduction packet
linked from that thread (gist `fqlx/93cce6ff81f26815d56d1c2865a0b9a9`, archive SHA-256
`b5ac61677d33ffb9174fc55b99ef842ab0e0906cf3403c89ed113220036b29a4`, 97 files, unpacked read-only under
`/tmp/hodge/gist-rank-one/unpacked/`), the companion instances, and the book chapters 01, 02, 03, 03a, 04, 05.
New computations made for this dossier live in `/tmp/hodge/scratch/` (scripts `rank_one_compute.py`, `ann7.py`,
`ann7_e2.py`, `annE2.py`, `tstruct.py`, `polar_sweep2.py`, `sing_r1/quot_p*.sing`) with their logs next to them.

Revision note (2026-11-01, second session). The first writing session crashed while several long computations were
still running in the background; this revision replaces every "pending" entry in section 4 with what the logs
actually contain (a killed run is reported as *not completed*, never as a result) and records the two computations
that were re-run afterwards (Singular `(J:A)` mod 10007; degree-8 annihilator `E_2` test). Nothing else changed.

Conventions. Variables are ordered `X0,X1,X2,Y0,Y1,Y2`; `S = C[X,Y]`; `J = J(F) = (F_0,...,F_5)` is the Jacobian
ideal; `tau` is the block-sign involution `Y -> -Y`; "Y-even/Y-odd" is parity of the total `Y`-degree;
`t = 5/7`; `E_X = sum_i X_i d/dX_i`, `E_Y = sum_j Y_j d/dY_j` are the block Euler fields.
Status words used below: **proved** (exact argument or exact computation whose lift to characteristic zero is
justified), **computational evidence** (exact modular computation whose lift is not justified, or a bounded search),
**written classical argument** (packet proofs not formalized and depending on listed inputs), **conjecture**.

Headline (for readers in a hurry).

* Open; no certificate, no nonconstructive proof. Discussion #118 has one post (fqlx, 2026-10-22) and no replies.
* The instance is the one-parameter deformation `F_t = f(X) + g(Y) + t L(X)^2 M(Y)` at `t = 5/7` of a split
  sextic for which the target `A = X0^3 Y0^3` has the *trivial* strict degree-6 certificate
  `D = f`, `N = A (X0,X1,X2,0,0,0)/6` (subset Euler, book 2.2). Everything about this instance is "what replaces
  that certificate when the product structure is broken".
* fqlx's packet proves (modulo written classical steps) that no certificate can have poles only on divisors pulled
  back from the "base" of the natural cubic cover `X/tau -> P^2_X x P^2_Y` or from the X-shape/scale threefold,
  and exhibits the exact obstruction left by the Euler-type partial current. It does not exclude denominators that
  genuinely correlate the X- and Y-curve coordinates, and it constructs nothing.
* New here: the exact Jacobian-ring annihilator profile of `[A]`. For Y-even denominators of degree `d <= 6`
  the necessary condition `A D in J` forces `D in J` (proved); the first non-Jacobian admissible denominators
  appear in degree 7 (a 3-dimensional space modulo `J`), then 9 in degree 8, 27 in degree 9 (modular upper bounds,
  stable under two primes at degree 7). The full annihilator ideal computed by Singular modulo 10007 is
  `(J : A) = J + (3 generators of degree 7, 3 of degree 8, 456 of degree 10)`, with `(J:A)/J` of dimension
  `3, 18, 54, 573, 1023, 1325` in degrees `7..12` (both parities) and multiplication by `A` surjective onto
  `(S/J)_{d+6}` for `d >= 12` (the latter is a proved statement over `Q`). Every member of the degree-7 family and
  every tested member of the degree-8 family fails the complete `E_2` system modulo 10007. All natural low-degree denominators (`L, M, L^2, LM, M^2, F_Xi, f, g, f+g,
  E_X(F), ...`) fail the complete fixed-denominator system modulo 10007. The exact deformation residual of the
  Euler field is `8 t A L^2 M / E_X(F)^2`, which is precisely the source fqlx's "elliptic19 passive current"
  fails to integrate.

---

## 1. Exact target

### 1.1 The record

`data/rank-one-biquartic/instance.rec` (SHA-256 of the record as quoted by the packet:
`f02e6ff8d1ed3b07c9c6715d26038aeff1eb940d3c067038ea6186ec24621420`):

* `q = 1`
* `A = X0**3*Y0**3`
* `F` = 54 monomials (quoted in full in the record). Verified exactly (sympy/flint, `rank_one_compute.py structure`):

```
F = f(X) + g(Y) + (5/7) * L(X)^2 * M(Y)

f(X) = X0^6 + X1^6 + X2^6 + (5/7) * ( X0^5 X1 + X0^2 X1^3 X2 + 3 X0 X2^5 + 2 X1^5 X2 )
g(Y) = Y0^6 + Y1^6 + Y2^6 + (5/7) * ( 2 Y0^5 Y1 + 2 Y0^3 Y1 Y2^2 + Y0 Y2^5 + 3 Y1^5 Y2 )
L(X) = X0^2 + X0 X1 + 2 X1 X2 + 3 X2^2            (a smooth conic; disc of its matrix is -14, a 101-unit)
M(Y) = Y0^2 + Y0 Y1 + 2 Y1 Y2 + 3 Y2^2 = L(Y)      (same quadratic form in the Y block)
```

Monomial support by block bidegree `(deg_X, deg_Y)`: 7 monomials of bidegree (6,0) (= `f`), 7 of bidegree (0,6)
(= `g`), 40 of bidegree (4,2) (= the 10 monomials of `L^2` times the 4 monomials of `M`). There is **no (2,4)
component** and no (3,3), (5,1), (1,5) component. The mixed part has tensor rank one as an element of
`S^4(X) (x) S^2(Y)`. All non-Fermat coefficients lie in `(5/7) Z`, so `F == sum Z_i^6 (mod 5)`.

Hidden relations (the second was found in fqlx's packet, both re-verified here):

```
M(Y) = L(Y)
g(Y) = f(Y2, Y0, Y1) + (5/7) * Y0^3 Y1 Y2^2
```

so with `Z = (Y2,Y0,Y1)` (cyclic relabelling) one can also write
`F = f(X) + f(Z) + (5/7) ( L(X)^2 N(Z) + Z0^2 Z1^3 Z2 )`, `N(Z) = M(Z1,Z2,Z0) = 3 Z0^2 + 2 Z0 Z2 + Z1^2 + Z1 Z2`,
`A = X0^3 Z1^3`. This is *not* a symmetry: the coupling `L(X)^2 N(Z)` is not symmetric under `X <-> Z` and the
extra monomial `Z0^2 Z1^3 Z2` breaks `f(X)+f(Z)` symmetry (see 1.4).

Pencils attached to the blocks (packet terminology): `u = f/L^3` on `P^2_X`, `k = g/M^3` on `P^2_Y`. Both are
Lefschetz pencils of plane sextics with 12 transverse base points (`f = L = 0`, resp. `g = M = 0`), exactly 31
finite nonzero Morse critical values each (critical-value polynomials of degree 31, irreducible over `Q`:
`X` mod 29, `Y` mod 17; packet `two_pencil_monodromy`, `CM_graph/critical_lifting.md`,
`elliptic_D18/D18_pencils_101.log`). The conic `L = 0` is parametrized by
`(X0,X1,X2) = (u(u+2v), -u^2-3v^2, v(u+2v))` (packet `exact_geometry.json`; verified).

### 1.2 Degrees, parity, characters

* `deg F = 6`, `deg A = 6 = 6q`. `F(X,-Y) = F`, `A(X,-Y) = -A`.
* Certificate degrees at `q = 1`: `deg D = d`, `deg N_i = d + 1`, `deg Q = 2d - 6` (so `Q = 0` is forced for
  `d <= 2`, and for `d = 3` by parity: `Q` would be a Y-odd constant), `deg E = 2d + 6`.
* The denominator may be taken Y-even without loss (book ch. 1 "One denominator"): if `v` serves `A` then
  `tau_* v` serves `-A`, so `(v - tau_* v)/2` serves `A` with denominator dividing `D * D(X,-Y)`, of degree `2d`.
  With `D` Y-even the component parities are: `N_0,N_1,N_2` Y-odd, `N_3,N_4,N_5` Y-even, `Q` Y-odd.
* Character of the residue. `tau` acts on `A` by `-1` and on the projective volume `Omega` (three `Y`'s) by `-1`;
  the residue of `A Omega/F^2` is therefore **tau-invariant**: the class lies in `H^4(X)^tau`, the sector the
  split theorem (3A.1) supports when `t = 0`. In Jacobian-ring terms the target is a Y-odd element of degree 6,
  i.e. a class with nonzero `H^{3,1}` component (Griffiths: degree `6 <-> H^{3,1}`), so it is not algebraic.
  `A` is nonzero in `S/J`: modulo 5, `J = (Z_i^5)` and `X0^3 Y0^3` survives; the packet also gives a `37 x 37`
  unit minor over `Z_(5)` (`two_pencil_base_support/source_bridge.md`). The `(4,0)` part is absent by parity.
* `checks.json` (regenerated by `scripts/audit-frontier.py`): `smoothReductionPrime = 5`, reduced partials
  `Z_i^5`, `targetNonzeroInReducedJacobianRing = true`, `coordinateReflections = []`,
  `subsetEulerAtoms = []`, `hessianSymmetrizerDimension = 1`, `projectiveDiagonalGroupOrder = 2`,
  `juniorElementsFixingTarget = []`.

### 1.3 Smoothness

Modulo 5 every non-pure coefficient vanishes, so `F == sum Z_i^6`, whose partials `Z_i^5` have no common
projective zero (`5` does not divide 6). Smoothness of the reduction implies smoothness over `Q` and `C`
(properness of the singular locus over `Z_(5)`). Hence `F` is irreducible, `S/J` is Gorenstein with socle degree
24 and Hilbert series `((1-t^5)/(1-t))^6`, and the parity-split Hilbert function equals the Fermat one:

| degree k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dim (S/J)_k^+ (Y-even) | 1 | 3 | 12 | 28 | 66 | 123 | 218 | 333 | 483 | 623 | 762 | 843 | 885 |
| dim (S/J)_k^- (Y-odd)  | 0 | 3 | 9 | 28 | 60 | 123 | 208 | 333 | 468 | 623 | 744 | 843 | 866 |

(symmetric about 12). The tau-invariant primitive Hodge numbers are `(208, 866, 208)` (Y-odd pieces in degrees
6, 12, 18), the anti-invariant `(1+218, 885, 218+1)` including `H^{4,0}`.

### 1.4 Symmetry (all recomputed independently; `rank_one_compute.py structure`)

* Projective diagonal symmetry group: since all six `Z_i^6` occur, diagonal symmetries have `lambda_i^6` equal,
  i.e. weights in `mu_6`; the brute-force screen over `mu_6^5` (last weight normalized) finds only the identity
  and `(3,3,3,0,0,0)`, i.e. `tau`. Order 2. Residue character of `A` under `tau`: polynomial character `3`,
  determinant character `3`, total `0 (mod 6)`: invariant (as above).
* Monomial symmetries (permutation times diagonal): the only permutation of the six variables preserving the
  monomial support of `F` is the identity; the full monomial automorphism group of `F` (up to scalars) is `{1, tau}`.
  In particular there is **no coordinate reflection, no block exchange, no permutation symmetry**, hence no
  character decomposition finer than tau-parity (book ch. 1: for a full-support Y-even sextic all Y-odd targets of
  one degree lie in a single character coset).
* Hessian symmetrizer space (matrices `T` with `Hess(F) T` symmetric): dimension 1 over `Q` (scalars). Hence no
  nontrivial linear direct-sum decomposition exists, even over `C`: the `3+3` split cannot be restored by a linear
  change of coordinates (why.md; confirmed).
* No "junior" diagonal element fixes the target (no McKay-type quotient route).
* Nobody in the thread or packet reports any non-diagonal symmetry; the packet explicitly checks that the
  paired-sign symmetries `(X_i,Y_i) -> (-X_i,-Y_i)` of the cubic-Jacobian method fail (e.g. the paired-0 difference
  has coefficient `-10/7` on `X0^5 X1`).

### 1.5 Block / quotient structure

Write a point of `P^5` off the fixed curves as `[X : rho Y]` with `X,Y` shapes and `rho` the block scale. Then

```
F(X, rho Y) = f(X) + t L(X)^2 M(Y) rho^2 + g(Y) rho^6 .
```

With `s = rho^{-2}` (the "inverse squared block scale" of why.md) and homogenizing variable `u`:

```
f(X) s^3 + t L(X)^2 M(Y) s^2 u + g(Y) u^3 = 0 ,        disc = -g * ( 4 t^3 L^6 M^3 + 27 f^2 g ) .
```

Thus `Q := X/tau -> B := P^2_X x P^2_Y` is a **triple cover**, branched along `g = 0` (where the double root
`s = 0` is the fixed curve `X = 0`) and along the genuine discriminant `Delta' := 4 t^3 L^6 M^3 + 27 f^2 g`
(bidegree (12,6); `4t^3 = 500/343`). At `t = 0` the cover is the *cyclic* cover `(s/u)^3 = -g/f`, i.e.
`(T_f x T_g)/mu_3` with triple planes `T_f: u^3 = f`, `T_g: v^3 = g` (3A.3). For `t != 0` it is **non-Galois**
with monodromy `S_3` (book 4.7, "natural cubic cover"): its rank-3 permutation local system splits as
`Q (+) V_std`, and the primitive problem lives in the standard summand; neither `Q` (= `C_2`-quotient of the Galois
closure) nor the discriminant double cover (`A_3`-quotient) nor `B` can serve `V_std` through the averaging ideal
of 3.1, because the first two are not uniruled (Kollar-Larsen: `Q` has torsion canonical class and terminal
`1/2(1,1,1)` points) and the last two have no `V_std`-invariants.

Two-pencil model (packet). On `L M != 0` the quotient is also the family, over the `(t_s,k)`-plane or over the
incidence surface `u + c t_s + k t_s^3 = 0` (`c = 5/7`, `t_s = M(Y)/L(X) * rho^2`, `u = f/L^3`, `k = g/M^3`),
of products of pencil curves `C_X(u) = {f - u L^3 = 0}` and `C_Y(k) = {g - k M^3 = 0}` (plane sextics,
genus 10, `H^1` of rank 20, punctured at the 12 base points rank 31). Equivalently `Q` is fibered over the
`k`-line with fibers `Z_k x C_k`, where

```
Z_k : f(X) + c L(X)^2 S + k S^3 = 0   in  P(1,1,1,2)_{X,S}      (S = M(Y) rho^2)
```

is a weighted sextic surface, `K = O(1)`, primitive `H^2` of rank 42 with Hodge numbers `(3,36,3)` (Jacobian
degrees 1, 7, 13), geometric generic Picard rank 1 (packet, see 3.2). The "rank-840 compact equation" is the
search for a rational primitive in `H^2_prim(Z_k) (x) H^1(C_k)` over `C(k)`; the "rank-441" version uses the
punctured curves; the "degree-12 coefficient field" is `C(x_0,x_1)` over `C(L,f)` (trace-zero rank 11, times
20 = 220).

### 1.6 The trivial certificate at `t = 0` and its exact deformation defect (new, elementary)

Because `deg_X A = 3 = |X-block|`, the subset Euler atom of book 2.2 applies to `F_0 = f + g`:
`v_0 = A E_X / E_X(F_0) = A E_X/(6 f)`, i.e. `D = f`, `N = A (X0,X1,X2,0,0,0)/6`, `Q = 0`, strict. (Also
`D = g` with `E_Y`.) So `rank-one-biquartic` is a deformation of an instance with a degree-6 strict witness.

For `t != 0` the same recipe keeps `E_1`: with `D_X := E_X(F) = 6 f + 4 t L^2 M = 6f + (20/7) L^2 M` and
`N = A X_i` (`i < 3`), `sum_i N_i F_i = A D_X` exactly, legal and of the right parity. Its normalized `E_2`
defect is exact and small:

```
B_{D_X}(N) = D_X div N - N(D_X) = A ( 6 D_X - E_X(D_X) ) = 8 t A L^2 M ,
div( A E_X / E_X(F) ) = 8 t A L^2 M / E_X(F)^2 .
```

The same residual `8 t A L^2 M / E_Y(F)^2` arises from `E_Y`, `E_Y(F) = 6 g + (10/7) L^2 M`. The certificate
problem is therefore equivalent to: find a rational vector field `w` tangent to `X` modulo `F`
(`dF(w) = F h`) with `div w - h == -8 t A L^2 M / E_X(F)^2 (mod F)`. This is exactly the source that fqlx's
"elliptic19 passive current" transports to a degree-6 elliptic cover (the packet's `H = V(F) = 6F6 + 4F4` is the
pull-back of `E_X(F)` and its retained residual `r = 8 a F4/(D H^2)` is this `8 t A L^2 M/E_X(F)^2`); see 3.2(10).

---

## 2. Companion and related solved instances

| instance | status | mechanism | denominator | why it does not transfer |
|---|---|---|---|---|
| `split-ternary-q2` (same `f,g`, `t=0`, `q=2`, `A=X0^4X1^3Y0^3Y1^2`) | accepted (fqlx, PR #159, 2026-10-20), existential, strict `Q=0` | 3A.2-3A.3: triple planes `T_f,T_g`; Fermat elliptic companions `W_f: f+S^3+T^3=0` in `P(1,1,1,2,2)` with moving-conic seeds mod 5 (13 coefficient equations, Jacobian det 1 in `F_5`, free/normal block det 4); second normal residue `R_2` for the third-order pole (`HigherResidues`); radial descent with weights (1,1,1,2); elliptic contraction with the Fermat curve `a^3+b^3+1=0`; normalized trace; legal clearing | not exported (`deg N = deg D + 7`) | the whole chain rests on `X/tau ~ (T_f x T_g)/mu_3`; for `t != 0` the fiber cubic is non-Galois, there is no product, no second triple plane, no fixed elliptic contraction (why.md, 4.7) |
| the `t=0` member for *this* target (`F_0=f+g`, `A=X0^3Y0^3`, not catalogued) | trivially solved | subset Euler, book 2.2 | `D = f` (deg 6), strict | `E_1` persists with `D = E_X(F)`; `E_2` defect `8tAL^2M` (1.6) |
| `crossed-ternary-asymmetric` (`F=P(X0,X1,Y0)+P(Y1,Y2,X2)`, `q=1`, `A=X0^4X1Y1`) | accepted (1337isnot1337, PR #149) | involution-stable ternary decomposition (3A.3): companion cone `P+S^3+T^3`, moving conic seed mod 5, elliptic kernels `K_0,K_1`, residue functional with infinity correction, translated elliptic contraction, trace; libraries `EllipticCompanions`, `PolynomialCurves` | not exported; stronger normalization `sum N_i F_i = A D`, `B_D(N) = F Q` with `Q != 0` allowed | needs a *linear* splitting into two ternary blocks (Hessian symmetrizer dim 2 there); here it is 1 |
| `quadric-contact-broken` (deformation `t=5/7` of `two-reflection-cubic-surface-d2`) | accepted, strict | cubic reduction `S=2(x^2+y^2), T=4xy, B=a^2, C=c^2`, `G=16F`; a line on the cubic over a finite coefficient algebra `E = O[r]/(P)`, `deg P = 27`; adjugate clearing | `deg D = 2n+58` (!), old `D = X0^2+Y0^2` excluded by `E_1` (rank 336/337) | model of "repairing a broken contact": the *old* denominator dies at first order and the fix is an entirely new finite-algebra geometry, not continuation |
| `holdout-ternary-deformation` (holdout + ternary terms, `A=X0^3Y0^3`) | accepted (hodgepodgeconjecture, PR #151) | Fermat cover `T^6 = Y0^6+Y2^6`, exchange `y <-> T`, weighted quotient `G` in `P(1,1,1,2,2)`, replacement conic seed mod 5 (13x13 unit Jacobian), lifting, polynomial integration, trace, deck character `T -> zeta T` | not exported | the Fermat-cover exchange requires the block `Y0^6+Y2^6` and a chain coupling; nothing of the sort survives here (`g` is dense in the `Y` block) |
| `cyclic-triple-fano`, `reflection-coupled-fano` | accepted, existential | weighted Fano quotients (`Z=X2^3`, resp. `u=Y2^2`) + free curve (rbajaj5's line) + polynomial-family integration (3.4, 3.6) | not exported | require a finite Piola map / reflection: none here (1.4) |
| `cycle-three`, `mixed-biweight`, `holdout-moving-branch` | accepted | paired signs + diagonal character; cubic reduction `u_i = X_i^2`; moving Weierstrass map, `TransitiveCubic`, `ResidueTransport` | not exported | paired signs absent (packet check, 1.4); the squared chart is not available |
| `four-plus-two-residual` (shares `f`) | open | - | - | benediamond (disc. #120): binary factor is a 3-dim Prym `(2,1)`; the visible elliptic quotient has the wrong character |
| `polarized-k3-remainder` (shares `f`) | open | - | - | 1337isnot1337 (disc. #115): `End_Hdg(T)=Q`, Picard `1<=rho<=3`; **all legal denominators of degree <= 7 excluded over C** by Euler gauge + Koszul parametrization + saturation; degree 8 undecided. This methodology is directly reusable here (section 6, item 2) |

Precisely why the 3A machinery fails (why.md + 04 section 4.7, restated): the split theorem needs `H^4(X)^tau` to
be assembled from `[V_f (x) V_g]^{mu_3}` with each factor supported by a Fano *threefold* companion
`(T x E)/mu_3` of dimension 3; the deformation destroys the product and replaces the cyclic cover by an `S_3`
cover whose standard summand is carried by non-uniruled quotients. "Producing an actual support correspondence
across this deformation is the unresolved step" (why.md). The quotient `Q` is not uniruled, so nothing like the
Fano curve theorem applies directly to it.

---

## 3. Prior work, exhaustively

### 3.1 The discussion thread (#118)

One top-level comment, no replies.

**fqlx, 2026-10-22T00:19:35Z, "Research update: checked reductions, but no certificate"** (AI-assisted, Codex and
ChatGPT). Content: three parallel leads returned without a certificate or an unrestricted existence proof.
Checked progress: (a) geometric generic Picard rank one of `Z_k` after finite base change (three holomorphic
residue seeds `Res(X_i Omega_w/H)` each generate all 42 primitive coordinates; cyclic determinants `95, 9, 73
mod 101` at `k=1`; integral Macaulay pencil units `38, 43, 64`; 74,088 connection Taylor entries audited by a
second jet method); (b) a scoped obstruction to all rational "compact-Y base currents" over the X/S threefold,
including the faithful degree-12 coefficient field (rank-220 trace-zero presentation) with arbitrary coefficient
poles; the simpler compact rank-840 equation also excluded; explicitly *not* excluding a divisor correlating the
X and Y curve coordinates; (c) an explicit elliptic partial current cancelling both compact moments but leaving a
nonzero 17-dimensional trace-zero residue (`r_0 = 51 mod 101`), not removable by cover trace; logarithmic shape
repair reduces to the known angular loop; bounded corrections and genus-ten matching tests fail within their
families. Fifteen focused replay commands pass. Three later "Pro" reports (mixed-support, conductor, residue,
degree-8/degree-9/global-conductor claims) were *not* replayed and are not counted. The post restates the exact
remaining contract (homogeneous `D != 0` coprime to `F`, `deg N_i = deg D + 1`, polynomial `Q`, no denominator
bound). Link: https://gist.github.com/fqlx/93cce6ff81f26815d56d1c2865a0b9a9 .

### 3.2 The packet, component by component (all under `research/rank-one-biquartic/` in the archive)

Replay harness: `verify.py` runs 15 commands (`REPLAY.json`), Python 3.12, SymPy 1.14.0, python-flint 0.9.0.
Every component binds the target by SHA-256 before computing (`relative_contraction/target.py`).

1. **`weighted_X_connection`** (build.py, witness.json 404 KB). The primitive rank-42 Gauss-Manin connection of
   `Z_k` over `Q(k)` by integral Macaulay pencils of `7H = 7f + 5L^2 S + 7k S^3` in weighted degrees 7, 13, 19
   (sizes 70, 308, 825; remainder columns 36, 3, 0; determinants at `k=1`: 38, 43, 64 mod 101). Block-column
   formula with factors 7, 2, 7/2 from the residue reductions (Batyrev-Cox Thm 10.13/Rem 10.14). The compact
   first variation `v_X = -(a h/6k) Res(Omega_w/H^2)`, `a = (5/7)L^2`, `h = X0^3`, is cyclic of rank 42
   (determinant 41 mod 101). Exact two-step Hodge identity `C_X B_X = 2 * reversal` over `Q(k)`.
   *Status:* exact modular arithmetic with a specified characteristic-zero lift (all entries in
   `Z_(101)[k,k^-1,d_7^-1,d_13^-1,d_19^-1]`). Explicitly not: irreducibility of the differential module, a primitive.
2. **`weighted_X_audit`**: independent reconstruction of everything in 1 (74,088 entries), quasi-smoothness of
   generic `Z_k` (the `X_i^19` argument plus the vertex `[0:0:0:1]` excluded for `k != 0`), Hilbert series
   `(1-t^5)^3(1-t^4)/((1-t)^3(1-t^2))`. *Proved* (exact).
3. **`holomorphic_picard`** + **`holomorphic_picard_audit`**: the three *holomorphic* seeds (degree-one monomials,
   pole 1) are cyclic (determinants 95, 9, 73 mod 101; two independent jet methods `d/dk + A` and `T' = T A`).
   Written classical argument: a divisor class over any finite extension of `C(k)` spreads to a flat `(1,1)`
   class, pairs to zero with a holomorphic seed and all its derivatives, hence its primitive part is zero:
   `NS(Z_{bar C(k)}) (x) Q = Q h`, `h^2 = 3`. *Status:* exact arithmetic + written argument; says nothing about
   special fibers or the integral Picard group.
4. **`two_pencil_monodromy`**: both pencils have 31 Morse critical points with distinct nonzero critical values
   (input from `pro_results/lead1/.../exact_geometry.json`, Singular lattices; resultant degree 43, base degree
   12); degree-31 critical polynomials irreducible over `Q` (Rabin checks mod 29 and mod 17 with Bezout
   identities); arithmetic transitivity on the 31 nodes; invariant-cycle theorem on a rational resolution gives zero
   fixed part; intersection graph connected (number of components divides 31 and is <= 10); Deligne Weil I Lemma
   5.11 gives geometric monodromy `Sp_20` for each compact `H^1`. *Status:* exact arithmetic + written argument.
5. **`CM_graph/critical_lifting.md`**: properness argument lifting the finite-field Lefschetz-pencil checks on the
   triple plane `T_P` (93 = 3 x 31 critical values of `u^3 = P` pencil) to characteristic zero; Katz SGA7 XVIII
   6.7 for finite-index absolute irreducibility. **`CM_graph/incidence_gysin.md`**: for any tau-invariant divisor
   whose quotient is dominated by a resolution of `S_X x_C S_Y` with mixed local system without invariants, the
   tau-invariant `H^{2,0}` Gysin image in `H^{3,1}(X)` vanishes (incidence class
   `c_2(E)(2H + k) = (H^2 + hH + h^2)(2H + k)` is a product of algebraic classes). Applied to earlier supports
   "D68" and "D18" (whose packets `D68_support_exclusion.md`, `D68_normalization.py`, `zero_fixed_part`,
   `moving_jacobian_route`, `x_hamiltonian_correction` are *referenced but not included*).
6. **`compact_vertical_support`**: the compact rank-840 equation `nabla z = -v_X (x) e dk` in
   `H^2(Z_k)_prim (x) H^1(C_k)` over `C(k)` has no rational solution. Proof sketch: graph residue identifies the
   physical double-pole density with the Leray class; a solution would kill the restriction of `eta` to the open
   `X2 Y2 M != 0` minus finitely many `D_k = {g - k M^3 = 0}`; every removed vertical support has zero
   tau-invariant `H^{2,0}` Gysin image (nodal `Z_k`, `k = 0` rational, `k = infinity` is `M = 0` whose quotient is
   covered by `T_f x C_B` with nontrivial `mu_3` character on all of `H^{2,0}`, coordinate hyperplanes by weak
   Lefschetz); contradiction with `eta != 0`. *Status:* written classical argument on stated inputs.
7. **`two_pencil_base_support`** + **`source_bridge.md`**: extends 6 to the full rank-441 equation with arbitrary
   rational coefficient poles on *curves of the two-dimensional base*: for a dominant base curve `Gamma`, the
   support is birational to `S_X x_Gamma S_Y`; Sp20 irreducibility kills cross `H^{2,0}` (a `(2,0)` tensor has rank
   <= 10 < 20); every holomorphic 2-form splits into surface pullbacks (multisection trace argument); each surface
   summand has zero physical Gysin image (classes `(H^2 + hH + h^2)(2H + T)`); boundaries `L = 0`, `M = 0`
   handled by cyclic covers with nontrivial character; affine Leray in total degree 4; regular relative de Rham
   source bridge without deleting the coordinate polars `J_X, J_Y`; the 37-column Jacobian witness for `A notin J`.
   *Status:* written classical argument; explicitly does not cover supports depending on both curve coordinates.
8. **`rotated_shape_correspondence`** + **`compactY_support.md`**:
   (a) the relation `g(Y) = f(Y2,Y0,Y1) + c Y0^3 Y1 Y2^2`;
   (b) an explicit degree-12 rational field solving `E_1` exactly: in the three variables `(X0, X1, Z2)` let
   `H = grad(U) x grad(q)` with `U = f(X) + f(Z)`, `q = L(X)^2 N(Z)`, `J = f_{X0} L_{X1} - f_{X1} L_{X0}`,
   `h = X0^3`; then `V = [h/(2c Z0^2 L N J)] H` satisfies `V(F) = A` with legal data
   `D = 2c Z0^2 L N J` (degree 12), `N_i = h H_i`, other `N_i = 0`. Its divergence is nonzero (value
   `501093/259210` at `X=(1,0,1), Z=(1,1,1)`) and "retained rather than declared divisible by `F`": **the complete
   `E_2` system for this `D` (with all Koszul corrections and `Q`) was not run** (see section 5). The field differs
   from the old angular current `V_old = A H_L/J` by `[h/(cZ0^2)] T_z` with `T_z` a `q`-preserving lift of
   `d/dZ2`;
   (c) the mixed support `D_rot = L(X)^2 M(Y) + Y0^3 Y1 Y2^2` (on which `F` restricts to the symmetric split
   `f(X) + f(RY) = 0`): its parity quotient is birational to `E x_{P^1_u} S_X` with `E` the rational double plane
   `w^2 = -Y0 Y1 M(Y)`, connected generic fibers (Stein factorization degree argument with `H`-degrees 12 and 18),
   Sp20 kills cross `H^{2,0}`, so `H^{2,0} = 0` and `D_rot` cannot support the class at any multiplicity;
   (d) `compactY_support.md`: no rational base 2-form current on the X/S threefold `T: f + cL^2 S + kS^3 = 0`
   with values in the compact rank-20 Y-module whose covariant differential is the physical source, including the
   degree-12 field model with arbitrary coefficient poles. Uses inputs 3, 4, 6, 7, Picard control of an individual
   Stein component, and the multisection/trace splitting. *Status:* written classical argument.
9. **`two_pencil_shape_current`**: exact degree-12 field split `rho = D_v Hbar + rho_log`,
   `rho_log = 56 z^5/P_z = -8/J` (fixed-`ell` obstruction: residue 8 at `z = 0`), repaired by the closed
   trace-zero form `lambda = (1/2) dlog((x0+2)/(x0+x1-2))`; the repaired current equals the known angular
   current (`q_0 = -X0^3/12 - L X0/7`, `U = L G^{-1} grad q_0 - 3 q_0 X`) up to an exact current and an explicit
   higher-jet current `chi wedge dk e_2`; no independent cancellation. Cyclic determinant 8 mod 101 for the full
   rank-20 Y source. *Status:* exact identities, checked.
10. **`elliptic19_transfer`** / **`elliptic19_passive`** / **`elliptic19_parent_bridge`**: a degree-6 cover
    `X_i = v z_i`, `Y = (q w - c_y/2, c_y, 1)`, `q^2 = beta(c_y) = 2c_y + 3 - c_y^2/4`, `w^2 = v^3 - 1`
    (fixed elliptic curve `E_0`, `omega_0 = dv/(2w)`), with `F_E` of degree 18 in `v` (18 simple zeros on `E_0`),
    exact density `a omega_0/F_E^2`, `a = -3q v^8 z0^3 (qw - c_y/2)^3`; the 18 residues `rho` reproduced by the
    trace kernel `kappa(P,Q) = (w_P + w_Q)/(v_P - v_Q) omega_P`; compact Hermite remainder `(c0,c1) = (74,13)
    mod 101` nonzero; the explicit passive current `t_i = -z_i a/(D H)` (`H = V(F) = 6F6 + 4F4`, `V` the
    X-Euler field: this is the pull-back of the Euler-polar field of 1.6) cancels both compact moments exactly and
    leaves the residual `r = 8 a F4/(D H^2)` (= the pull-back of `8 t A L^2 M/E_X(F)^2`), trace-zero projection
    `r_0 = 51 mod 101`, invariant under the deck groups `mu_3, mu_2` (so cover trace does not remove it).
    Two bounded searches (rigorous within their families via nonzero augmented minors):
    (i) `t_i = phi_i a/(D F_i) + phi_E t_i^{Euler}`, `phi` polynomial in `z0,z1,z2,c_y,q` of degree `<= d`, `d = 0..3`
    (5, 30, 100, 250 unknowns, 864 equations, augmented rank one more each time);
    (ii) `t_i = sum_{r<=R} phi_{ir} e_{ir}(v,w) a/(D H^r)`, `R = 1,2,3`, `phi` in the 6 affine conic functions,
    `e` in `L(5 infinity)` (120, 240, 360 unknowns, augmented minor determinants 18, 20, 24 mod 101).
    `parent_bridge`: the residual density `6T^2 B_rem/J_Y` takes different values (38 vs 22 mod 101) on two
    physical branches over the same point of the "positive Jacobian parent" (`xi, eta, U, V` isogeny formulas),
    so it is not a scalar source on that parent; also the literal squared-coordinate chart of the accepted
    transitive-cubic proof needs paired signs that `F` lacks.
11. **`quadric_section_matching`**: constant coordinate permutations with `k`-independent quadric sections cannot
    match the `g`-pencil curve on the rotated surface (rank-2 obstruction modulo `M`); first-order deformation of
    the graph `X' = Y, S = -M` in the direction `c(T + N^2 M)` around the actual `g` seed: normalized `28 x 15`
    matrix of rank 15, augmented 16, dual functional `Lambda(T + N^2 M) = 11883702676272 != 0`, sole kernel the
    weighted scaling; genus-10 sections of `|2H|` on `Z_k` can only be matched by degree-one maps (Riemann-Hurwitz).
12. **`pro_results`**: `lead1/.../exact_geometry.json` (Singular lattices for the two pencils, conic parametrization,
    Bezout witnesses mod 5) and `lead2/lead2-new/broader_matrix/witness.json` (degree-12 field: primitive
    polynomial with ascending coefficients `7(ell-7)^6, 178(ell-7)^5, ...`, trace-free rank 11, joint Y rank 220,
    `certificate_found: false`).

### 3.3 Classification of the prior results

* **Proved / exact and lifted:** decomposition identities; pencil critical structure (31 Morse points, distinct
  values, irreducibility of the degree-31 polynomials); connection data and cyclic ranks over `Q(k)`;
  `A notin J`; smoothness; the bounded-ansatz exclusions in 10(i),(ii) (full augmented column rank with nonzero
  minors is a characteristic-zero statement *for those finite families*); the first-order non-deformation in 11.
* **Written classical arguments (not Lean, conditional on the listed inputs, independently unreviewed):** generic
  Picard rank 1 of `Z_k`; `Sp_20` monodromy; the three support obstructions (compact rank-840; rank-441 with base
  curve poles; compact-Y base currents over `T` with the degree-12 field); exclusion of the support `D_rot`.
  Their common scope: denominators whose zero loci are pulled back from the base of the cubic cover or from the
  X-shape/scale threefold `T`, or are the boundaries `L = 0`, `M = 0`, coordinate hyperplanes, or `D_rot`.
  **Not covered:** any divisor depending jointly on the X-curve and Y-curve coordinates (e.g. any non-bihomogeneous
  `D(X,Y)` not of the form polynomial in `X, M(Y), g(Y)`).
* **Screened-method failures (no nonexistence content beyond the family):** the two elliptic19 ansatz tables; the
  logarithmic repair loop; the quadric-section matching; the parent-bridge non-descent.
* **Not replayed / not claimed:** three later Pro reports (degree-8, degree-9, global-conductor claims).

### 3.4 Links

* Thread: https://github.com/hodge-computer/sextics/discussions/118
* Packet: https://gist.github.com/fqlx/93cce6ff81f26815d56d1c2865a0b9a9 (files `README.md`, `REPLAY.json`,
  `unpack.py`, `rank-one-checked-research.zip.b64`; raw zip SHA-256 `b5ac6167...29a4` verified here)
* References cited by the packet: Batyrev-Cox https://arxiv.org/pdf/alg-geom/9306011#page=34 ; Deligne Weil I
  https://www.numdam.org/item/PMIHES_1974__43__273_0.pdf#page=22 ; Katz SGA7 XVIII
  https://web.math.princeton.edu/~nmk/old/etudecoho.pdf#page=73 ; Stacks 0BTX.
* Reusable methodology from a neighbour: 1337isnot1337's degree-by-degree exclusion packet for
  `polarized-k3-remainder`: https://gist.github.com/1337isnot1337/d504e828c4f50d8dade0687792ea9ec0 .
* Accepted constructions referenced above: PR #159 (`split-ternary-q2`), #149 (`crossed-ternary-asymmetric`),
  #151 (`holdout-ternary-deformation`), #128 (`holdout`), #165 (`cycle-three`), #170 (`mixed-biweight`), #153
  (`cyclic-triple-fano`).

---

## 4. New computations for this dossier

All scripts in `/tmp/hodge/scratch/`; host tools python-flint 0.9.0, sympy 1.12, Singular 4.4.1 (also in the
`ddg-cas` docker image). Nothing was written inside the repository.

### 4.1 Structure (proved, exact rational arithmetic)

`F = f + g + (5/7)L^2 M`, `M = L(Y)`, `g(Y) = f(Y2,Y0,Y1) + (5/7)Y0^3Y1Y2^2`; coordinate reflections none;
diagonal group `{1, tau}`; monomial automorphism group `{1, tau}` x scalars; Hessian symmetrizer dimension 1;
parity-split Hilbert function as in 1.3.

### 4.2 The `E_1` (annihilator) profile of the target (`rank_one_compute.py e1`, `p = 10007`)

For a certificate with Y-even `D` of degree `d`, `E_1` (book 5.2) forces `A D in J`, i.e. `[D] in Ann_{S/J}([A])`.
The table gives, for each `d`, the kernel of `S_d^+ -> (S/J)_{d+6}^-`, `D -> [A D]`, computed modulo 10007, its
trivial part `J_d^+`, and the difference `dim Ann_{S/J}([A])_d^+`.

| d | rows (Y-odd deg d+6) | J columns | rank_p(J) = expected dim J_{d+6}^- | dim S_d^+ | kernel | dim J_d^+ | **dim Ann mod J** |
|---|---|---|---|---|---|---|---|
| 0 | 226 | 18 | 18 | 1 | 0 | 0 | 0 |
| 1 | 396 | 63 | 63 | 3 | 0 | 0 | 0 |
| 2 | 636 | 168 | 168 | 12 | 0 | 0 | 0 |
| 3 | 1001 | 378 | 378 | 28 | 0 | 0 | 0 |
| 4 | 1491 | 756 | 747 | 66 | 0 | 0 | 0 |
| 5 | 2184 | 1386 | 1341 | 126 | 3 | 3 | 0 |
| 6 | 3080 | 2376 | 2214 | 236 | 18 | 18 | 0 |
| 7 | 4284 | 3861 | 3441 | 396 | 66 | 63 | **3** |
| 8 | 5796 | 6006 | 5052 | 651 | 177 | 168 | **9** |
| 9 | 7752 | 9009 | 7129 | 1001 | 405 | 378 | **27** |
| 10 | not computed (the `d = 10` run, `e1_profile_d10.log`, was killed at start). Dimension bound, exact over any field: kernel `>= dim S_10^+ - dim (S/J)_16^- = 1512 - 468 = 1044` (`(S/J)_16^- ~ (S/J)_8^-` by Gorenstein duality, socle Y-even), `dim J_10^+ = 1512 - 762 = 750` | | | 1512 | >= 1044 | 750 | **>= 294** (bound) |

Both parities together (Singular, `(J:A)` mod 10007, 4.6): `dim ((J:A)/J)_d = HF(S/J)_d - HF(S/(J:A))_d` is
`0, 0, 0, 0, 0, 0, 0, 3, 18, 54, 573, 1023, 1325, 1260, 1260, 1120, 895, 645, 420, 245, 126, 56, 21, 6, 1` for
`d = 0..24`. Comparison with the Y-even column: degree 7 is entirely Y-even (3 = 3); degree 8 splits `9 + 9`; degree 9
splits `27 + 27`. From degree 12 on the annihilator has exactly the dimension forced by `dim (S/J)_d - dim (S/J)_{d+6}`
(multiplication by `A` is surjective in those degrees); in degrees 10 and 11 it exceeds the forced value by 18 and 3.

Rigor. In every row `rank_p(J block) = dim_Q J_{d+6}^-` (the regular-sequence value, which is the same over `Q` by
the connectedness/semicontinuity argument), so `dim_Q ker <= dim_p ker`; since `J_d^+ subset ker` over any field,
rows with `kernel = dim J_d^+` are **proved over C**: for `d <= 6`, every Y-even `D` with `A D in J` lies in `J`.
Consequences (proved): no certificate whose Y-even denominator has degree `<= 4`; no certificate with any
denominator of degree `<= 2`; a certificate with Y-even `D` of degree 5 or 6 must have `D = dF(G)` for a block-linear
polynomial field `G` (`G_i in <X0,X1,X2>` for `i < 3`, `G_j in <Y0,Y1,Y2>` for `j >= 3`), 3 resp. 18 parameters.
For `d = 7, 8, 9` the nonzero values are modular *upper* bounds; the degree-7 value 3 was recomputed at `p = 65521`
with identical dimension (`ann7_p65521.log`). The sequence 3, 9, 27 is *not* generated in degree 7 alone: Singular (4.6) finds minimal generators of `(J:A)`
beyond `J` in degrees 7 (three), 8 (three) and 10 (456), none in degree 9, 11 or higher.

Structure of the degree-7 family (`ann7.py`, `ann7_e2.py`): the three representatives modulo `J_7^+` are
irreducible (flint factorization mod p), have 106-109 terms, leading monomial `X0^7`, and are **not** of the form
`l(X) P` or `l(Y) Q`: for every linear form `l` the space `{P in S_6^(+/-) : l P passes E_1}` is exactly `J_6`
(dimension 18). Natural products `X_i f, X_i g, X_i L^2 M, X_i L^3, X_i L M^2, X_i M^3` all fail `E_1`.
Base-type test (`tstruct.py`, `p = 10007`, `tstruct_d7.log`): the span of the 73 polynomials
`P(X) M(Y)^a g(Y)^b` of degree 7 (depending on `Y` only through `M, g`, i.e. pulled back from the packet's
X-shape/scale threefold `T`) meets the 66-dimensional `E_1` kernel in a 27-dimensional space, of which 24 dimensions
are `J_7^+ cap span`; so **all 3 annihilator directions of degree 7 are `T`-base-type** (the pattern
`X0^5 * (1,1,2,3)(Y0^2,Y0Y1,Y1Y2,Y2^2) = X0^5 M` visible in the representatives is the general phenomenon).
Consequence: every `E_1`-admissible Y-even denominator of degree `<= 7` is either in `J` (a polar `dF(G)`) or a
`T`-base-type polynomial; both kinds lie inside the scope of the packet's compact-Y obstruction (3.2(8d)), which
therefore predicts that `E_2` fails on the whole degree-7 family. The direct `E_2` tests agree (4.3).
Status of the follow-ups: the second-prime confirmation of `dim Ann mod J = 3` at degree 7 is done
(`ann7_p65521.log`, identical dimensions and term counts 106/109/106); the degree-7 base-type test at `p = 65521`
(`tstruct_d7_p65521.log`) and the degree-8 and degree-9 base-type tests (`tstruct_d8.log`, `tstruct_d9.log`) were
killed by the session crash before producing output and were **not completed**. So "base-type" is established
(modulo the single prime 10007) only in degree 7; in degrees 8 and 9 it is an open question whether the 9- and
27-dimensional annihilator families contain genuinely mixed denominators.

### 4.3 Complete fixed-denominator systems (`rank_one_compute.py fixed`, `p = 10007`)

For a fixed `D` the full cleared system `E = F^2 Q` is linear in `(N, Q)` with the component parities of 1.2
(all monomials included). `rank(M)` / `rank([M|b])`:

| D | deg | rows | N unknowns | Q unknowns | rank M | rank [M b] | result |
|---|---|---|---|---|---|---|---|
| L, M, L+M, L-M | 2 | 1491 | 168 | 0 | 159 | 160 | inconsistent |
| L^2, LM, M^2 | 4 | 5796 | 756 | 9 | 705 | 706 | inconsistent |
| F_X0, F_X1, F_X2, F_X0+F_X1+F_X2 | 5 | 10152 | 1386 | 60 | 1316 | 1317 | inconsistent |
| f, g | 6 | 16797 | 2376 | 226 | 2337 | 2338 | inconsistent |
| f+g (= F - tL^2M) | 6 | 16797 | 2376 | 226 | 2352 | 2353 | inconsistent |
| D_X = E_X(F) = 6f + (20/7)L^2M | 6 | 16797 | 2376 | 226 | 2352 | 2353 | inconsistent |
| D_Y = E_Y(F) = 6g + (10/7)L^2M | 6 | 16797 | 2376 | 226 | 2352 | 2353 | inconsistent |
| f-g (not in `J`, see note) | 6 | 16797 | 2376 | 226 | 2352 | 2353 | inconsistent (`fixed_D.log`, last completed row) |
| L^3, L^2M, LM^2, M^3 | 6 | | | | | | not run (killed); **excluded by proof**: none is in `J` (`ann7_p10007.log`: `X_i * D notin J` for each), so `A D notin J` by 4.2 |
| F-g = f + tL^2M, F-f = g + tL^2M | 6 | | | | | | not run (killed); **excluded by proof**: `J cap span{f, g, L^2M} = span{E_X(F), E_Y(F)}` (it contains these two; if 3-dimensional then `f in J`, contradicting 4.2), and `f + tL^2M`, `g + tL^2M` are not in that span (`4t/6 != t`, `2t/6 != t`) |
| D_X - D_Y = 6(f-g) + (10/7)L^2M, two random block-linear polars dF(G) | 6 | | | | | | **not completed** (killed by the crash); these are in `J`, so `E_1` holds and only the `E_2` test is informative |

Note on `f - g`: `E_X(F) - 2E_Y(F) = 6(f - 2g)`, so `f - 2g in J`, but `f - g = (f - 2g) + g` is **not** in `J` since
`g notin J`; its row is therefore also covered by the 4.2 proof, and the modular rank pattern `2352/2353` (the same as
for the `J`-members `f+g, D_X, D_Y`) is a coincidence of ranks, not of mechanism.

The modular inconsistency is only evidence in each row (the coefficient matrix has a kernel, e.g. 9-dimensional for
`d = 2`), but for `L, M, L^2, LM, M^2, f, g, f+g` the exclusion is in fact **proved** by 4.2 (they are not in `J`,
so `A D notin J`). For `F_Xi`, `D_X`, `D_Y` (which are in `J`, so `E_1` holds trivially) the `E_2` failure is modular
evidence; a characteristic-zero proof would need the rank over `Q` of the 16797 x 2602 system (feasible with
fraction-free elimination, not done).

Normalized-system version for the degree-7 annihilator family (`ann7_e2.py`; `E_1` lift read off the nullspace,
`E_2` in degree 14: 5796 rows, 420 Koszul unknowns `c_ij` cubic with parity, 636 unknowns `C in S_8^-`):
all three representatives and a random member of the 66-dimensional affine `E_1` family are inconsistent
(`rank 996 / 997` in every case; 260-809 s each). The run was killed by the session crash after random member #0;
further random members and the pure polar members of `J_7^+` (`dF(G)`, `G` quadratic block-parity) were **not
tested** (`ann7_e2.log`). Also recorded there: for each variable `l in {X_i, Y_j}`, `{P in S_6 of the right parity :
l P passes E_1}` has dimension exactly 18 = `dim J_6`, i.e. no degree-7 admissible denominator is `l * P` with
`P notin J`.

### 4.4 Exhaustive sweep of the degree-5 polar family (`polar_sweep2.py`)

`D = l0 F_X0 + l1 F_X1 + l2 F_X2`, `N_0 = (l0 A, l1 A, l2 A, 0,0,0)`, normalized `E_2` with Koszul gauge
(`45` linear `c_ij` with parity, `60` unknowns `C in S_4^-`, 1491 rows), swept over all points of `P^2(F_p)`,
`p = 101, 103, 107` (10303 / 10713 / 11557 points). Result: **not completed.** The run was killed by the session
crash after the precomputation and the first point; the only datum is the point `(1,0,0)` at `p = 101` with ranks
`101 / 102` (inconsistent), consistent with the `D = F_X0` row of the table above (`polar_sweep2.log`). The sweep
remains to be done (estimated 2 s per 100 points after precomputation, i.e. minutes per prime). A consistent point
would be a candidate to lift; absence at three primes would be evidence (not proof) that the whole 3-parameter family
fails over `C`. Note that the family is tiny: by 4.2 it is the *complete* set of `E_1`-admissible Y-even
denominators of degree 5 (the kernel is exactly `J_5^+`, dimension 3).

### 4.5 The Euler deformation residual (proved, 1.6)

`B_{E_X(F)}(A X) = 8 t A L^2 M`, `div(A E_X/E_X(F)) = 8 t A L^2 M / E_X(F)^2`; identical form for `E_Y`.

### 4.6 Singular ideal quotient `(J : A)` (`sing_r1/quot_p10007.sing`, log `quot_p10007_resume.log`)

Re-run in the second session (Singular 4.4.1, `Z/10007[X0..Y2]`, degrevlex). `std(J)`: 1912 elements, 34 s;
`quotient(J, A)`: 604 s. Results, all modulo 10007:

* **Minimal generators of `(J:A)`: 468.** Six of degree 5 (the partials `F_i`, i.e. `J` itself), **3 of degree 7**
  (208 terms each, leading monomials `X1^5 X2^2`, `X0 X1^4 X2^2`, `X0^3 X1^2 X2^2`, all Y-even, each irreducible
  over `F_10007` and not in `J`), **3 of degree 8** (351 terms each, leading monomials `X1^3 X2^5`, `X0 X1^2 X2^5`,
  `X0^2 X1 X2^5`, Y-even, irreducible, not in `J`), none of degree 9, and **456 of degree 10** (459-475 terms each,
  not in `J`); none of degree `>= 11`. Hence `Ann_{S/J}([A])` is generated by 3 + 3 + 456 classes in degrees 7, 8, 10.
* **Hilbert function of `S/(J:A)`**, degrees 0..18:
  `1, 6, 21, 56, 126, 246, 426, 663, 933, 1192, 933, 663, 426, 246, 126, 56, 21, 6, 1`, zero from degree 19 on.
  It is symmetric about 9 with socle degree 18, as it must be: `S/(J:A) ~ A (S/J)`, a principal ideal in the
  Gorenstein ring `S/J` (socle degree 24), which is Gorenstein with socle degree `24 - 6 = 18`. The value
  `HF(S/(J:A))_k` is the rank of multiplication by `A` from `(S/J)_k` to `(S/J)_{k+6}`; the modular rank is a lower
  bound for the rank over `Q`, so the modular annihilator dimensions `dim (S/J)_k - HF_k` are upper bounds over `Q`,
  and where they equal the trivial lower bound `dim (S/J)_k - dim (S/J)_{k+6}` (degrees `>= 12`) they are exact.
  Hilbert function of `S/J` (degrees 0..24): `1, 6, 21, 56, 126, 246, 426, 666, 951, 1246, 1506, 1686, 1751, 1686,
  ..., 1`, the Fermat value, as expected.
* **Membership tests `A D in J`** (1 = yes): `L, M, f, g, L^2M`: 0 (confirms 4.2); `D_X, D_Y`: 1; `D_X in J`: 1;
  `f - 2g in J`: 1; the degree-12 base-type products `f g, f L^2M, g L^2M, (L^2M)^2, f^2, g^2`: all 0, so none of the
  natural degree-12 polynomials in `(f, g, L^2 M)` is `E_1`-admissible (the admissible degree-12 space is the
  1325-dimensional one above, but it is not generated by these products). Further: `A L^3, A M^3, A LM, A L^2, A M^2`: 0. **`Delta' = 4t^3 L^6 M^3 + 27 f^2 g` (degree 18) lies in `J`
  modulo 10007** (`Disc in J: 1`; hence trivially `A Delta'`, `A g Delta' in J`). This is modular evidence that the
  genuine branch divisor of the cubic cover is a polar `dF(G)` with `G` of degree 13, i.e. an `E_1`-admissible
  denominator *for the trivial reason*; its `E_2` test (degree 36, ~`1.1e5` rows) was not run. A characteristic-zero
  proof of `Delta' in J` would be a 33649 x 33649-ish linear solve or a Groebner reduction over `Q` (not done).
  The `det Hess(F)` test (degree 24) had not printed when this revision was written (see the end of the log; `A det Hess`
  is trivially in `J` by degree).
* The characteristic-zero run (`quot_p0.sing`) was killed in `std` by the session crash and was **not** re-run.

Interpretation. The annihilator is tiny through degree 9 (3, 18, 54 total, half of each Y-even for `d = 8, 9`) and
becomes large exactly at degree 10 (573, with 456 fresh generators). Every `E_1`-admissible denominator of degree
`<= 9` is therefore a combination `sum_i l_i G_i^{(7)} + sum_j m_j G_j^{(8)} + (element of J)` with the six low
generators, and the first place where "new" denominators not built from these six appear is degree 10, matching
the search-space statement in sections 5-6 (the degree-10 Y-even family is at least 294-dimensional modulo `J`).

### 4.7 Degree-8 annihilator family, complete `E_2` test (`annE2.py 10007 8`, log `annE2_d8.log`)

Re-run in the second session with a 1700 s limit. `E_1` nullspace at `d = 8`: nullity 1131 (108 s), kernel of
`D -> [AD]` of dimension 177 with `dim J_8^+ = 168`, so `dim Ann mod J = 9` as in 4.2. Normalized `E_2` in degree
16: 10152 rows, 954 Koszul unknowns, 1491 unknowns `C in S_10^-`. Results (each ~200-280 s):

| member of the degree-8 family (mod `J_8^+`) | D in J | rank M | rank [M b] | result |
|---|---|---|---|---|
| rep0 | no | 2286 | 2287 | inconsistent |
| rep1 | no | 2286 | 2287 | inconsistent |
| rep2 | no | 2286 | 2287 | inconsistent |
| rep3 | no | 2286 | 2287 | inconsistent |
| rep4 .. rep8, two random members of the 177-dim affine family | | | | **not reached** before the 1700 s limit |

Members not reached before the time limit were **not tested**. The pattern (`rank [M|b] = rank M + 1` for every
tested member, in degrees 7 and 8 alike) is the signature expected if the obstruction is a single linear functional
(the Gysin/residue functional of the packet's written argument) that does not vanish on the whole family.

---

## 5. What remains unexplored (concrete gaps)

1. **No complete denominator exclusion by degree.** The packet never ran the fixed-denominator linear algebra of
   book 5.1 for *any* `D`, nor the `E_1` profile. (4.2-4.4 start this; the rigorous part stops at degree 6 for `E_1`
   and at nothing for `E_2`.) The polarized-K3 packet shows the full methodology (Euler gauge, Koszul
   parametrization `N = A a + K(c)`, enlarged correction ideal `L = (F) + J^2 + (R_ijk)`, saturation, projective
   special-fibre argument, msolve/Singular) up to degree 7; nothing comparable exists here.
2. **The degree-7/8/9 annihilator families were unknown.** Their generators, their geometric meaning (base-type
   or genuinely mixed), and the `E_2` consistency locus inside the 66-, 177-, 405-dimensional affine `E_1`-families
   have never been examined; the family of *all* `E_1`-admissible denominators of degree 10 (dimension `>= 294`
   modulo `J`, plus `J_10^+`) even less.
3. **The 18-parameter block-linear polar family** `D = dF(G)` (the complete degree-6 continuation family of the
   `t = 0` witness) has only been sampled at a few points (`D_X`, `D_Y`, ...). The joint quadratic system
   (18 + 162 + 226 unknowns, 3080 equations) was never solved or excluded.
4. **The packet's own degree-12 Nambu field** `V = h H/(2c Z0^2 L N J)` was never subjected to the complete `E_2`
   test (all Koszul corrections of degree 8 and `C in S_18^-`); only its raw divergence was evaluated. Cost:
   normalized `E_2` in degree 24: ~59k rows x (9675 + 16797) unknowns; modulo `F` first, ~59k x 9.7k; feasible
   with a sparse modular solver, hours with dense flint.
5. **Denominators correlating X and Y curve coordinates.** Everything excluded so far is base-type or boundary
   (3.3). No candidate family of genuinely mixed denominators has been proposed except `D_rot` (excluded).
   Candidates never tested: non-bihomogeneous combinations such as `alpha f + beta g + gamma L^2 M + (mixed
   bidegree pieces)`, the degree-7/8/9 annihilators (if they turn out mixed), products with the fixed-curve
   equations, polars `dF(G)` for non-linear `G` (degree 7+), `det Hess(F)` (degree 24, passes `E_1`).
6. **Galois closure of the cubic cover.** Nobody has written the `S_3` Galois closure `Z -> Q` (the double cover of
   `Q` branched along `Delta'`) with its `mu_3` deck action over the discriminant double cover `B'`, nor tried a
   single Fermat elliptic companion `(Z x E)/mu_3` (a fivefold) or a curve/surface fibration of `Z` serving `V_std`.
   Book 4.7 only shows that the *intermediate quotients* of `Z` cannot serve `V_std`.
7. **Deformation in `t` as a function-field problem.** A certificate over `Q(t)` for the family `F_t` (then
   specialized at `t = 5/7`) has not been attempted; neither has a formal power-series continuation of the
   `t = 0` witness with growing denominator degree (which is what a `Q(t)`-certificate with denominators in `t`
   amounts to).
8. **q is fixed at 1**, so the "raise q" tricks (book 4.5) and `D^2` presentations (prescribed-support theorem,
   `q = 2`) are not available; conversely nobody has tried to find a certificate for a *different* Y-odd target on
   the same `F` (e.g. one in the degree-6 image of a known support) and to relate targets through the Jacobian
   ring: since `A` is a *single* monomial and the Jacobian ring has no finer characters, any support that serves
   one nonzero tau-invariant `H^{3,1}` class might serve `A` only through a full-rank argument (book 4.4
   "one generic nonzero target suffices" needs a dominating family, not available for a fixed `F`).
9. **Symmetry sectors:** there are none beyond `tau`; the near-symmetry `g = f(RY) + cB` has been used only for
   the rotated-support attempt.
10. **Nonconstructive routes** have not been pursued beyond fqlx's negative results: no attempt to prove
    `H^4(X)^tau subset N^1 H^4` by a geometric support construction specific to this non-Galois triple cover.

---

## 6. Candidate attack plan

Sizes below use the exact parity counts (q = 1, Y-even `D` of degree `d`): `N` unknowns `3 dim S_{d+1}^- +
3 dim S_{d+1}^+`, `Q` unknowns `dim S_{2d-6}^-`, cleared rows `dim S_{2d+6}^-`; normalized form: `E_1` in degree
`d + 6`, `E_2` in degree `2d` with Koszul unknowns `c_ij` of degree `d - 4` (parity) and `C in S_{2d-6}^-`.

| d | dim S_d^+ | N unknowns | Q unknowns | cleared rows S_{2d+6}^- | Koszul unknowns | E_2 rows S_{2d}^- |
|---|---|---|---|---|---|---|
| 5 | 126 | 1386 | 60 | 10152 | 45 | 1491 |
| 6 | 236 | 2376 | 226 | 16797 | 162 | 3080 |
| 7 | 396 | 3861 | 636 | 26532 | 420 | 5796 |
| 8 | 651 | 6006 | 1491 | 40326 | 954 | 10152 |
| 9 | 1001 | 9009 | 3080 | 59332 | 1890 | 16797 |
| 10 | 1512 | 13104 | 5796 | 84903 | 3480 | 26532 |
| 12 | 3108 | 25704 | 16797 | 162248 | 9675 | 59332 |

1. **Finish the exact low-degree map (degrees 5-9), rigorously.** Idea: combine the `E_1` profile with the
   correction criterion. For `d = 5, 6` the admissible denominators are `dF(G)`, `G` block-linear; the certificate
   equation becomes the quadratic system `B_{dF(G)}(A G) + T_{dF(G)}(c) = F C` in `(G, c, C)` (3/18 + 45/162 +
   60/226 unknowns, 1491/3080 equations). Solve it exactly: Groebner basis over `F_p` and over `Q` (Singular,
   msolve as in the K3 packet), or eliminate `(c, C)` and compute the rank-drop locus in `G`. For `d = 7, 8, 9`
   do the same on the affine `E_1`-families (66, 177, 405 parameters; the `E_1` lift `N_0` comes with the nullspace).
   Chance of finding a witness: low (every sampled point fails); chance of a *proved* exclusion up to degree 9: high
   with a few CPU-days, and it would be the first rigorous low-degree statement for this instance. Required for
   any honest "small denominator" discussion.
2. **Identify the annihilator ideal `(J : A)` exactly and interpret it.** Singular `quotient` over `Q` (running);
   factor the degree-7 generators, test base-type membership (`tstruct.py`), and compute the Hilbert function of
   `(J:A)/J` to predict where `E_1`-admissible denominators become plentiful (degree 10: `>= 294`). If the degree-7
   generators are `T`-base-type, the packet's obstruction says they cannot work and the first genuinely mixed
   admissible denominators must be searched in degree `>= 8`. Cheap, high information value.
3. **Joint nonlinear search in the forced-kernel degrees 10-12 restricted to mixed denominators.** At `d = 10`
   the admissible space modulo `J` is `>= 294`-dimensional; restrict to the complement of base-type polynomials
   (those in the span of `P(X) M^a g^b` and `Q(Y) L^a f^b`) and to candidates vanishing to order along the fixed
   curves... then run variable projection (least squares in `(c, C)` for each `D`, 26532 x (3480 + 5796)) with
   Levenberg-Marquardt over the `D`-parameters, many random starts, exact reconstruction of any numerical zero
   (book 5.4). Chance: unknown but this is the first search space not excluded by anything; cost: hours per run.
4. **Complete `E_2` test of the packet's degree-12 field** (gap 4). Compute `B_D(N_0)` for
   `D = 2c Z0^2 L N J`, `N_0 = h H`, and decide `-B_D(N_0) in F S_18^- + im T_D` by sparse elimination modulo 10007
   (reduce modulo `F` first). Prediction from the packet's geometry: inconsistent (all components base-type or
   boundary). Value: an independent arithmetic check of the written obstruction; low chance of a witness.
5. **Integrate the Euler residual directly on the cubic-cover geometry.** The problem is reduced (1.6) to finding
   a tangent field `w` with `div w - h == -8 t A L^2 M/E_X(F)^2 (mod F)`. On the two-pencil model the source is
   `8t (X0^3 L^2)(Y0^3 M)/E_X(F)^2`, a product of an `X`-form and a `Y`-form divided by a `T`-base-type square.
   fqlx integrated it along a *fixed* elliptic curve in the `Y` direction and got a rank-17 obstruction. The
   untried variant: a *moving* genus-one or genus-10 fibration, i.e. the `Y`-pencil curve `C_k` itself with its
   Jacobian, coupled to the `X`-surface `Z_k` (a current with values in `H^1(C_k)` whose coefficients depend on
   the `X`-curve coordinate as well — precisely the "correlated" class the obstructions leave open). Concretely:
   choose the finite cover where `S = M(Y) rho^2` and `T^3 = g - M^3` (the packet's `elliptic19_parent_bridge`
   coordinates), write the source in `(X, S, T, Y0)`, and look for primitives whose coefficients are polynomial in
   the *X-curve* coordinate on `C_X(u)` and rational in the `Y`-coordinates. Chance: moderate as a research
   direction, low as a quick experiment.
6. **Galois closure companion.** Build `Z = {(s, s') : two roots of f s^3 + tL^2M s^2 + g = 0}` (a double cover of
   `Q` branched along `Delta'`) with its cyclic `mu_3` action over `B'`, then the single companion
   `W = (Z x E)/mu_3` with the Fermat curve `E` and `(u, a, b) -> (zeta u, zeta^-1 a, zeta^-1 b)`. At `t = 0` this
   degenerates to the product of the two 3A companions. Questions to settle before any computation: is `W` (a
   fivefold) uniruled / Fano in a weighted model, and does the `V_std (x) H^1(E)` sector have an algebraic
   contraction analogous to 3A.3? Only if both answers are yes does the Chapter 3 curve machinery (which needs
   surface fibers with algebraic `H^2`, 3.2) have a chance; dimension 5 is beyond everything used so far. Chance of
   a complete proof: low-to-moderate, long.
7. **A two-pencil correspondence of "diagonal" type.** The target class lives in `H^1(P^1_k, R^3)` with
   `R^3 = H^2(Z_k)_prim (x) H^1(C_k)` (rank 840) at the compact level. `H^2(Z_k)_prim` has rank 42 and
   `H^1(C_k)` rank 20; both pencils are Lefschetz with `Sp_20` monodromy on the `C` side; the `X`-pencil
   `C_X(u)` also has genus 10 and `Sp_20`. A correspondence of the kind that supports split-type classes is a
   *translated diagonal* in a product of two curve families over a common base after a base change making the two
   pencil parameters dependent (`u = -ct_s - kt_s^3`). The packet's `D68`/`D18` supports were of this fibre-product
   type and were excluded by the Gysin computation (`incidence_gysin.md`). A support **not** of fibre-product type
   over a base curve — e.g. the graph of a correspondence between `C_X(u)` and `C_Y(k)` varying with *both*
   `(u, k)` — is the open case; the natural candidates are Hecke-type/`l`-isogeny correspondences between the two
   genus-10 Jacobians along a surface in the `(u,k)`-plane where they are isogenous (a Noether-Lefschetz-type
   locus). Testing isogeny loci numerically (periods) is feasible; building a certificate from them is not.
   Chance: speculative.
8. **Certificate over `Q(t)`.** Set up the normalized system for `F_t` with unknown polynomial coefficients in `t`
   and denominator degrees growing with the `t`-order (power-series continuation of `D = f`, `N = AX/6`). The first
   order already shows the obstruction `8 t A L^2 M` must be absorbed by `E_2`-corrections with `D_1` of degree 6
   in `J(F_0)`; by 4.2 applied to `F_0` (same Hilbert data) the degree-6 admissible space is again
   `dF_0(G)`, so a bounded-degree series cannot exist and any solution must have denominators that are rational
   functions of `t`. Worth a one-day experiment to see whether the order-by-order systems are solvable with
   denominator degree growing linearly in the order (which would point at a pole in `t` and a genuinely new
   support). Chance: low for a witness, moderate for insight.
9. **Use the sister instances.** `four-plus-two-residual` and `polarized-k3-remainder` share `f`; their packets
   contain `f`-specific geometry (Picard/Frobenius data for `w^2 = f`, Prym computations). Nothing of this is
   usable until a support for the coupled sextic is proposed; but the K3 packet's exclusion code is the cheapest
   way to upgrade items 1-2 to proofs.

Overall assessment. Low-degree explicit witnesses (degree `<= 9`) look very unlikely: all admissible denominators
through degree 6 are polar (in `J`) and fail `E_2` at every sampled point, and the small annihilator families in
degrees 7-9 fail on all tested members. The first unexplored and unexcluded space is "genuinely mixed admissible
denominators of degree `>= 10`", a `>= 294`-parameter nonlinear problem. The geometric obstruction program of fqlx is
converging on the statement that the support must correlate the two pencil curves; no construction of such a
correlation exists in the project's toolkit (every accepted proof uses a product, a reflection, a Fano quotient, a
Fermat cover exchange, or paired signs, none of which is present).

---

## 7. Feasibility: explicit witness vs. nonconstructive argument; what a Lean proof needs

* **Explicit witness.** Would require a denominator of degree at least 10 (if the modular bounds hold and the
  degree 7-9 families are fully excluded), with numerators of degree `>= 11` and `Q` of degree `>= 14`; the cleared
  identity has `>= 85k` monomials. Verification in Lean of such a witness is routine via the `Verify` interface
  (the accepted `klein-rational` witness has `deg D = 30`), but *finding* it is the problem; the only
  search spaces with a chance are the nonlinear ones of items 3, 5, 8.
* **Existential (finite-algebra) proof.** Every accepted "nonconstructive" proof in the repository is in fact a
  finite algebraic construction: a dominant polynomial curve family on an auxiliary quasismooth hypersurface with
  a polynomial normal (`PolynomialCurves`), integration in the curve parameter, normalized trace, Piola transport,
  legal clearing (`HoldoutCover`, `EllipticCompanions`, `HigherResidues`, `ResidueTransport`, `TransitiveCubic`).
  All of them need an auxiliary variety with (i) a polynomial source matching `A` after density correction and
  (ii) a curve family whose existence is seeded by a nonsingular incidence modulo a prime. For this instance the
  missing object is the auxiliary variety itself (item 6 or 5 above). If a companion `W` with a `mu_3`-type
  character selection and a dominant polynomial curve family were found, the Lean side would reuse
  `HoldoutGeneralization.certificate_of_source_pair`-style entry points plus a new descent (the deck group of the
  non-Galois cover is `S_3`, so the trace/character bookkeeping is new) — several thousand lines, comparable to
  PR #159.
* **Hodge-theoretic nonconstructive proof** (support via cycles without an explicit curve family) is not
  admissible unless formalized; no accepted proof has gone that way and the required inputs (resolution,
  Gysin maps, Hodge theory) are far outside mathlib. Not feasible in the near term.
* **Nonexistence** is not a goal of the challenge and is not claimed by anyone; all negative results are scoped.

---

## Appendix A. Reproduction of the new computations

```
python3 /tmp/hodge/scratch/rank_one_compute.py structure          # 2 s
python3 /tmp/hodge/scratch/rank_one_compute.py e1 9              # ~12 min; e1 10 10 for degree 10
python3 /tmp/hodge/scratch/rank_one_compute.py fixed [names...]  # minutes per degree-6 denominator
python3 /tmp/hodge/scratch/ann7.py 10007 7 ; python3 ann7.py 65521 7
python3 /tmp/hodge/scratch/ann7_e2.py 10007
python3 /tmp/hodge/scratch/annE2.py 10007 8 2   # ~30 min; degree-9: annE2.py 10007 9 (not run)
python3 /tmp/hodge/scratch/tstruct.py 10007 7
python3 /tmp/hodge/scratch/polar_sweep2.py 101
Singular -q /tmp/hodge/scratch/sing_r1/quot_p10007.sing   # ~11 min to (J:A) and generators; output quot_p10007_resume.log
```

All linear algebra is exact modulo the stated prime (python-flint `nmod_mat.rank/nullspace`), polynomial arithmetic
with `nmod_mpoly`/`fmpq_mpoly`; the rigor of each lift to characteristic zero is stated in 4.2-4.3.

## Appendix B. Key formulas collected

```
F_t = f(X) + g(Y) + t L(X)^2 M(Y),      t = 5/7
E_X(F) = 6 f + 4 t L^2 M,   E_Y(F) = 6 g + 2 t L^2 M,   E_X(F) + E_Y(F) = 6 F
f - 2 g = (E_X(F) - 2 E_Y(F))/6  in J
t = 0 witness:  D = f,  N = A (X0, X1, X2, 0, 0, 0)/6,  Q = 0
t != 0:  D = E_X(F), N = A (X0,X1,X2,0,0,0)  solves E_1;  B_D(N) = 8 t A L^2 M
cubic cover of P^2 x P^2:  f s^3 + t L^2 M s^2 u + g u^3 = 0,  disc = -g (4 t^3 L^6 M^3 + 27 f^2 g)
Z_k : f + t L^2 S + k S^3 = 0  in P(1,1,1,2);  C_k : g - k M^3 = 0;  C_X(u) : f - u L^3 = 0
packet degree-12 E_1 field:  V = X0^3 (grad U x grad q)/(2 t Z0^2 L N J), U = f(X)+f(Z), q = L(X)^2 N(Z),
    J = f_{X0} L_{X1} - f_{X1} L_{X0},  Z = (Y2, Y0, Y1)
```
