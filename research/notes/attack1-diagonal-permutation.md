# Attack 1: `diagonal-permutation` (q = 1, A = X0^3 Y0^3 - X1^3 Y1^3)

Status line (updated as the session proceeds): **no witness found yet**.

Work directory: `/tmp/hodge/results/diagonal-permutation/` (log of every tested denominator in `log.jsonl`,
scripts `hlib.py`, `e1kernel*.py`, `jideal.py`, `families.py`, `batch*.py`).

## 1. The obstruction and the geometric mechanism a certificate needs

(see the final version at the end of the session)

## 2. Experiments

(tables filled in below as they complete)

### 1.1 Setup

* `F = P + 5 G` with `P = sum X_i^6 + sum Y_i^6` (Fermat) and
  `G = (1/7) sum_{i != j} X_i^3 X_j Y_i Y_j + (1/11) sum_{(i,j,k)} X_i^3 X_j Y_k^2 + (1/13) sum_{(i,j,k)} X_i^2 Y_j^3 Y_k`
  (sums over ordered pairs `(i,j)`, `k` the third index).  X-degree parts: `F_6 = sum X^6`, `F_4` (Y-degree 2),
  `F_2` (Y-degree 4), `F_0 = sum Y^6`.
* Symmetry group of `F` (projective, linear): `S_3 x <tau>`, where `S_3` permutes the pairs `(X_i,Y_i)` simultaneously and
  `tau: Y -> -Y`.  Computed: the only projective diagonal symmetries are `tau` and the identity (lattice `L(F)` has
  Smith form `(1,1,1,1,2,0)`), so there are **no reflections, no paired signs, no cyclic diagonal symmetries** beyond `tau`.
* `A = T_0 - T_1`, `T_i = X_i^3 Y_i^3`.  Under `S_3` the span of `{T_0 - T_1, T_1 - T_2}` is the 2-dimensional standard
  representation; `det = +1` for simultaneous permutations, so the residue class itself transforms in `std` (no volume sign).
  `A` is `tau`-odd, `(01)`-anti-invariant; `B_2 = T_0 + T_1 - 2 T_2` spans the `(01)`-invariant line of the same `std` copy.
* `S_3`-isotypic decomposition of the primitive target space `R_6^- = (S/J(F))_6^-` (dim 208): `39 triv + 31 sign + 69 std`
  (std counted with multiplicity, i.e. 138 dimensions).  `A` is one direction inside the 69-dimensional multiplicity space
  of `std` — there is **no small sector** that the symmetry alone isolates.

### 1.2 Why the quotient constructions of the book are all unavailable (and why every volume-preserving quotient is)

* Reflection theorem (§3.3) needs a reflection; there is none.  Age criterion (§3.6): every nontrivial fixed stratum of
  `S_3` and of `tau`, `tau*sigma` has age `>= 1` (checks.json; e.g. a transposition has codimension-2 fixed locus with
  eigenvalues `(-1,-1)`, age 1).
* `S_3` and all its subgroups preserve the holomorphic 4-form `Omega/F` (det = +1).  Hence any (crepant resolution of a)
  quotient `X/K`, `K <= S_3`, still has `p_g = 1`.  Consequently **no surface fibration with `p_g = q = 0` fibres and no
  uniruled model can exist for any such quotient** (both would force `H^{4,0} = 0`).  This explains structurally why
  `benediamond`'s transposition-quotient fibration (discussion #177) has fibres with `p_g = 3`: it is forced, not an accident
  of the chosen fibration.  Only quotients involving `tau` kill the 4-form, and `X/tau`, `X/<tau,sigma>` have Kodaira
  dimension 0 (book §4.7), not uniruled.
* Averaging ideal (§3.1): `std^{K} = 0` for `K = S_3`; `std^{(01)}` is 1-dimensional (spanned by `B_2`), which is the
  content of `benediamond`'s reduction: it suffices to support `B_2` on `X/(01)`, but that quotient is a `p_g = 1` fourfold.

### 1.3 What a certificate must look like geometrically

By the Gysin description (§3.1), a divisor `Y = X ∩ V(D)` supporting the nonzero `(3,1)` class `Res(A Omega/F^2)` must have a
component whose resolution has `h^{2,0} > 0`; since `H^2(Y, O_Y) = 0` for every effective divisor on `X`, this forces `Y`
to be **non-normal or to have non-rational singularities** (cone points of multiplicity `>= 4` with a surface base of
`p_g > 0`, curves of singularities with non-rational transverse type such as simple-elliptic `x^3+y^3+z^3`, or a double
surface giving a non-trivial double cover).

The decisive structural fact found in this session (§2.1 below) is that **for `deg D <= 8` (both parities) every
`D` with `A D ∈ J(F)` lies in the Jacobian ideal `J(F)`, with exactly one exception at degree 8.**  So all the "cone" families
of the experiment menu (quartic cones with a line vertex, cubic cones with a plane vertex, products of linear forms,
Fermat-type sums of powers, Hessian minors, partial sums of `F`) fail the necessary condition `E1` outright
up to degree 8, and the only denominators of degree `<= 8` that can work are `D = G(F) = sum_i G_i dF/dZ_i`
(the tangency locus of the polynomial vector field `G` with `X`), plus the single exceptional `D_8`.

## 2. Experiments

Engine note: the shared engine `/home/kejace/devel/sextics-work/hodge` had a sign bug in the `F*N(D)` term that was
fixed at ~17:25 (engine v2).  All fixed-D consistency results below were (re)computed with engine v2; E1 tests and the
Jacobian-ring linear algebra never used the buggy term.  Prime for all mod-p computations: `p = 2147483647`
(E1 kernel reps verified with 2 further primes).  Every tested denominator is in
`/tmp/hodge/results/diagonal-permutation/log.jsonl` (field `engine`).

### 2.1 E1 kernel sweep: `K_d = {D in S_d^{+/-} : A*D in J(F)}` modulo `J(F)_d`  (scripts `e1kernel.py`, `e1kernel_odd.py`)

| d | dim S_d^+ | dim K_d (Y-even) | dim J_d^+ | primitive kernel (even) | dim S_d^- | primitive kernel (odd) | time |
|---|---|---|---|---|---|---|---|
| 1 | 3 | 0 | 0 | **0** | 3 | **0** | <1 s |
| 2 | 12 | 0 | 0 | **0** | 9 | **0** | <1 s |
| 3 | 28 | 0 | 0 | **0** | 28 | **0** | <1 s |
| 4 | 66 | 0 | 0 | **0** | 60 | **0** | 1 s |
| 5 | 126 | 3 | 3 | **0** | 126 | **0** (dim K = dim J = 3) | 3 s |
| 6 | 236 | 18 | 18 | **0** | 226 | **0** (18 = 18) | 10 s |
| 7 | 396 | 63 | 63 | **0** | 396 | **0** (63 = 63) | 30 s |
| 8 | 651 | 169 | 168 | **1** (= `D_8`, below) | 636 | **0** (168 = 168) | 80 s / 240 s |
| 9 | 1001 | (running) | 378 | ? | | | |

Consequences (rigorous mod p, hence over Q since ranks can only drop mod p):
* **No certificate with `deg D <= 4` exists** (no D of either parity passes E1).
* For `5 <= deg D <= 8`, every admissible `D` lies in the Jacobian ideal `J(F)` (i.e. `D = sum_i G_i dF/dZ_i`,
  `deg G_i = d-5`, parities forced), **except** for one 1-dimensional exception at `d = 8` (Y-even).
* By dimension counting (`R_d^+ -> R_{d+6}^-`: 483 -> 744 at d=8, 623 -> 623 at d=9, 762 -> 468 at d=10) a nonzero
  primitive E1 kernel is forced only from `d = 10` on (dimension `>= 294`); the degree-8 element is a genuine
  special relation `A * D_8 = 0` in the Jacobian ring, not a dimension-count artefact.
* `S_3` character of the relevant Jacobian-ring pieces (`triv, sign, std` multiplicities):
  `R_6^- = (39, 31, 69)`, `R_6^+ = (46, 28, 72)`, `R_8^+ = (94, 67, 161)`, `R_14^- = (132, 116, 248)`.

### 2.2 Low-degree explicit candidates (batch1, all fail E1 as the sweep predicts)

Quadrics (23 of them: coordinate sums of squares, `X_i^2 +/- Y_i^2`, `X0X1 +/- Y0Y1`, transposition eigenforms
`(X0-X1)^2 +/- (Y0-Y1)^2`, `(X0+X1)^2+(Y0+Y1)^2`, `sum (X_i-X_j)^2`, `X0Y1 +/- X1Y0`, ...) and Y-even cubics (16: `X0^3+X1^3+X2^3`,
`X0X1X2`, Vandermonde `(X0-X1)(X1-X2)(X2-X0)`, `X0Y0^2+X1Y1^2+X2Y2^2`, `X0Y1Y2+...`, `X0Y0^2-X1Y1^2`, Hesse cubic, ...):
**all fail E1** (`A*D` not in `J(F)`), each in `< 0.5 s`.  Together with the sweep this covers *every* `D` of degree `<= 4`.

### 2.3 Jacobian-ideal sub-families with prescribed singularities (`jideal.py`, `families.py`, `dims_scan.py`)

Since any `D` of degree `5..8` must be `G(F)`, the geometric mechanism has to come from the vector field `G`.
`Y_D = X ∩ V(G(F))` is the locus where `G` is tangent to `X`.  Imposing non-rational singularities is a *linear*
condition on `G`: `G(F) ∈ I^k` for a linear ideal `I`.  Dimensions of `{D ∈ J_d^+ : D ∈ (l_1..l_m)^k}` modulo
`S_{d-6}·F` (adapted coordinates: `p,r = X0±X1`, `q,s = Y0±Y1`, `a = X2`, `b = Y2`; `e = X0+X1+X2`, `r1 = X0-X1`,
`r2 = X1-X2`, similarly `s1, s2`):

| ideal (singular locus of `Y_D`, transverse type) | d=7 | d=8 | d=9 |
|---|---|---|---|
| `(X0,X1,X2)^3`: triple along the genus-10 curve `Z_Y = X∩{X=0}`, cubic cone (simple elliptic) | **21** | 39 | 254 |
| `(Y0,Y1,Y2)^4`: quadruple along `Z_X = X∩{Y=0}`, quartic cone (genus-3) | 0 | **31** | 138 |
| `(X0,X1,X2)^4` | 0 | 39 | 63 |
| `(X0,Y0,X1)^3`, `(X0,Y0,X2)^3` (Fermat-type cone planes of `T_0`) | 0 | 37 | 181 |
| `(X0,Y0,Y1)^3`, `(X0,Y0,Y2)^3` | 0 | 32 | 170 |
| `(X0,Y0,X1)^4`, `(X0,Y0,Y1)^4` | 0 | 0 | 61 / 51 |
| `(X0,Y0)^2`, `(X2,Y2)^2`, `(X0,X1)^2`, `(Y0,Y1)^2`: non-normal (double) along a sextic surface `X∩P^3` | 9 | 46 / 46 / 47 / 45 | 145 / 145 / 148 / 144 |
| `(X0,Y0)^3`, `(X2,Y2)^3` | 0 | 12 | 61 |
| `(X0,X1,X2,Y0)^4`: quadruple points at the 6 points of `X ∩ {X=0,Y0=0}`... (vertex line) | 17 | 84 | 244 |
| `(X0,X1,Y0,Y1)^4`: quadruple points at `X ∩ (X2,Y2)-line` | 10 | 68 | 245 |
| `(r,s)^2`: double along the `(01)`-fixed surface `S = X∩{X0=X1,Y0=Y1}` | 17 | 66 | 185 |
| `(r,s)^3`: triple along `S` | **3** | 24 | 91 |
| `(p,q,a,b)^4`: quadruple points at the 6 points of the `(01)`-anti-invariant line | 0 | 52 | 245 |
| `(p,q,a,b)^3` | 20 | 112 | 305 |
| `(r,s,a)^3`, `(r,s,b)^3`: elliptic cone along `X∩{X0=X1,Y0=Y1,X2=0}` resp. `Y2=0` | **6 / 6** | 48 / 48 | 178 / 176 |
| `(p,q,a)^3`, `(p,q,b)^3` | 0 / 2 | 30 / 38 | 174 / 168 |
| `(r,a,s,b)^4` | 0 | 51 | 245 |
| `(r1,r2,s1,s2)^4`: quadruple points at the 6 `S_3`-fixed points `X ∩ {X0=X1=X2, Y0=Y1=Y2}` | **2** | 53 | 245 |
| `(r1,r2,s1,s2)^3` | 22 | 112 | 305 |
| `(r1,r2)^2`, `(s1,s2)^2`: double along `X∩{X0=X1=X2}`, `X∩{Y0=Y1=Y2}` | 9 / 9 | 45 / 45 | 146 / 144 |

Explicit description of the 21-dimensional `(X)^3 ∩ J_7` family (derived by hand, verified): with `{a,b} = {0,1,2} \ {j}`,
`D = sum_j [alpha_j (Y_a^2 - Y_b^2) + Q_j(X)] F_{X_j} + (5/39) sum_j alpha_j X_j (Y_a F_{Y_b} - Y_b F_{Y_a})`,
`alpha ∈ Q^3` (the permutation representation of `S_3`: `triv + std`), `Q_j` arbitrary X-quadratics (18 parameters; this
part is trivially in `(X)^3` because every `F_{X_j}` vanishes on `{X = 0}`, a consequence of Y-evenness).  The
`alpha`-part uses the specific shape `F_2 = (5/13) sum_i X_i^2 (Y_a^3 Y_b + Y_a Y_b^3)` of the X-degree-2 part of `F`;
for a general Y-even sextic the analogous space is zero.

### 2.4 Attack plan for the resumed session (18:05; the first session crashed with batch3/5/2b/4 and the d=9 kernel unfinished)

Priority order (two solver processes at a time, memory on the machine is shared and tight, ~15 GB free):
1. **`D8`** — the unique E1-admissible non-Jacobian denominator of degree <= 8 — full fixed-D test (engine v2) incl. served
   target dimension; its `(01)`-symmetrisations; a structural probe of the coset `D8 + J_8^+` (sparse representatives,
   Fermat-annihilator ancestry, vanishing loci).
2. Engine-v2 re-run of the Jacobian-ideal members of degree 5 and 6 (`batch2b`: partials, block Euler fields, pair-radial,
   swap/rotation/cyclic fields, all `Z_i F_{Z_j}`), with served-target dimensions.
3. `batch3` (d=7 triple-along-`Z_Y` family, 21-dim, incl. the explicit `D_alpha`), `batch5` (d=7 adapted families:
   double/triple along the `(01)`-fixed surface, elliptic cones, quadruple points), `batch4` (d=8 `(Y)^4`, `(X)^4` families).
4. E1 kernel at d=9 (both parities if time permits); if the primitive kernel is small, extract and test every element.
5. Mathematics: duality interpretation of `Ann([A])`, the vector-field reformulation of the full system for `D = G(F)`,
   and what geometric mechanism remains at degree 9-10.
