# Attack 1: `diagonal-permutation` (q = 1, A = X0^3 Y0^3 - X1^3 Y1^3)

**Status: NO WITNESS FOUND.**  Rigorous output of this session: *no certificate exists with `deg D <= 7` (either Y-parity),
and the only possible denominator of degree 8 modulo `J(F)` is a single explicit polynomial `D8`, which is tested below.*
(Session 1 crashed at 17:40; resumed 18:00 on a machine at load ~300 shared with ~15 sibling agents, so throughput was low.)

Work directory: `/tmp/hodge/results/diagonal-permutation/` — every tested denominator is in `log.jsonl` (field `engine`: 1 = buggy
first engine, results unreliable and all re-done; 2 = fixed engine).  Scripts: `hlib.py` (test harness), `e1kernel*.py`,
`e1k8_extract.py` (exact D8), `d8_probe.py`, `family_kernel.py`, `jideal.py`/`families.py`/`dims_scan.py` (Jacobian-ideal
sub-families), `batch*.py`, `run_D8.py`.  Prime for all mod-p work: `p = 2147483647` (D8 verified with 4 further primes).

## 1. The obstruction and the geometric mechanism a certificate needs

### 1.1 Setup
* `F = P + 5 G`, `P = ΣX_i^6 + ΣY_i^6` (Fermat), `G = (1/7) Σ_{i≠j} X_i^3 X_j Y_i Y_j + (1/11) Σ X_i^3 X_j Y_k^2 + (1/13) Σ X_i^2 Y_j^3 Y_k`
  (`(i,j,k)` a permutation of `(0,1,2)`).  Equivalently `F = Fermat + (5/7) O1 + (5/11) O2 + (5/13) O3` with the three `S3×τ`
  orbit sums `O1, O2, O3` (verified to reproduce `F` exactly).
* Symmetry group of `F`: `S_3 × <τ>` (`S_3` permutes the pairs `(X_i,Y_i)` simultaneously, `τ: Y -> -Y`).  The only projective diagonal
  symmetries are `τ` and the identity (lattice of diagonal symmetries has Smith form `(1,1,1,1,2,0)`): **no reflections, no paired
  signs, no diagonal cyclic symmetries**.  Every nontrivial element of `S_3` has age `>= 1` on its fixed loci.
* `A = T_0 - T_1`, `T_i = X_i^3 Y_i^3`; `{T_0-T_1, T_1-T_2}` spans the standard representation `std` of `S_3` (det = +1 for simultaneous
  permutations, so the residue class itself transforms in `std`).  `A` is `τ`-odd and `(01)`-anti-invariant.
* `S_3`-isotypic decomposition of the primitive Y-odd sextics `R_6^- = (S/J(F))_6^-` (dim 208): `39 triv + 31 sign + 69 std`; `A` is one
  direction inside the 69-dimensional multiplicity space of `std` — the symmetry isolates no small sector.

### 1.2 Hodge-theoretic meaning and why all book mechanisms are unavailable
`ω_A = Res(A Ω/F^2)` is a `τ`-odd primitive `(3,1)`-class.  `H^4(X)^-` has no `(4,0)` part (the 4-form is `τ`-even), so it has level 2
and the generalized Hodge conjecture predicts a divisor `Y` with `ω_A ∈ Gysin(H^2(Ỹ)(-1))`; a certificate `(D,N,Q)` realises this with
`Y = X ∩ V(D)`.  Since `H^{2,0}(Ỹ) = 0` for every effective divisor with rational singularities, **`Y_D` must have non-rational
singularities**: cone points of multiplicity `>= 4` with surface base of `p_g > 0`, curves of singularities of non-rational transverse
type (simple-elliptic `x^3+y^3+z^3`, genus-3 cone `x^4+y^4+z^4`, ...), or non-normal loci giving nontrivial covers.
* Reflection theorem: no reflection.  Age criterion: every fixed stratum has age `>= 1`, so no quotient `X/K` is uniruled.
* `S_3` preserves `Ω/F` (det = +1): every quotient `X/K`, `K <= S_3`, has `p_g = 1`, hence **no surface fibration with `p_g = q = 0` fibres and no
  uniruled model exists for any such quotient**; this is why `benediamond`'s transposition-quotient fibration (discussion #177) has fibres with
  `p_g = 3` — forced, not accidental.  Only quotients involving `τ` kill the 4-form, and `X/τ`, `X/<τ,σ>` have Kodaira dimension 0.
* Averaging ideal: `std^{S_3} = 0`; `std^{(01)}` is spanned by `B_2 = T_0 + T_1 - 2T_2`, which is `benediamond`'s reduction (it suffices
  to support `B_2` on `X/(01)`, a `p_g = 1` fourfold).  At the certificate level: a `(01)`-invariant `D` serving `A` can be taken with
  `(01)`-equivariant `(N,Q)` (average the two certificates), so the `(01)`-anti-invariant sector of the fixed-`D` system is all that matters.

### 1.3 Algebraic structure: Gorenstein duality for the E1 kernel
`R = S/J(F)` is Gorenstein with Y-even socle in degree 24, so `R_d^ε × R_{24-d}^ε -> R_24` is perfect and
`rank(A : R_d^ε -> R_{d+6}^{1-ε}) = rank(A : R_{18-d}^{1-ε} -> R_{24-d}^ε)`.  With `k_d^ε = dim Ann(A)_d^ε` (= dimension of the
primitive E1 kernel) this gives `k_d^ε - k_{18-d}^{1-ε} = dim R_d^ε - dim R_{18-d}^{1-ε}`.  Parity-split Hilbert function of `R`
(Koszul, `((1-t^5)/(1-t))^3 · ((1-(ty)^5)/(1-ty))^3`):

| d | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `R_d^+` | 1 | 3 | 12 | 28 | 66 | 123 | 218 | 333 | 483 | 623 | 762 | 843 | 885 | 843 | 762 | 623 | 483 | 333 | 218 |
| `R_d^-` | 0 | 3 | 9 | 28 | 60 | 123 | 208 | 333 | 468 | 623 | 744 | 843 | 866 | 843 | 744 | 623 | 468 | 333 | 208 |

Dimension counting forces a nonzero kernel only from `d = 10` (`k_10^+ >= 294`, `k_10^- >= 261`, `k_11^± >= 510`, ...).  Combining
with the computed values `k_d = 0` (`d <= 7`), `k_8^+ = 1`, `k_8^- = 0` (§2.1) and duality: **`k_10^+ = 294`, `k_10^- = 262`, `k_11^± = 510`,
`k_12^+ = 677`, `k_12^- = 648`, ..., `k_9^+ = k_9^-` (unknown; the only degree not fixed by the `d <= 8` computation).**
So the E1-admissible denominators modulo `J(F)` are completely known except at `d = 9`: none for `d <= 7`, the single line `D8` at `d = 8`.

## 2. Experiments

### 2.1 E1 kernel sweep `K_d = {D ∈ S_d^± : A·D ∈ J(F)}` modulo `J(F)_d` (`e1kernel.py`, `e1kernel_odd.py`; exact mod p, hence valid over Q)

| d | dim S_d^+ | dim K_d^+ | dim J_d^+ | primitive kernel (Y-even) | dim S_d^- | primitive kernel (Y-odd) | time (unloaded machine) |
|---|---|---|---|---|---|---|---|
| 1 | 3 | 0 | 0 | **0** | 3 | **0** | <1 s |
| 2 | 12 | 0 | 0 | **0** | 9 | **0** | <1 s |
| 3 | 28 | 0 | 0 | **0** | 28 | **0** | <1 s |
| 4 | 66 | 0 | 0 | **0** | 60 | **0** | 1 s |
| 5 | 126 | 3 | 3 | **0** | 126 | **0** | 3 s |
| 6 | 236 | 18 | 18 | **0** | 226 | **0** | 10 s |
| 7 | 396 | 63 | 63 | **0** | 396 | **0** | 30 s |
| 8 | 651 | 169 | 168 | **1** (`D8`) | 636 | **0** | 80 s / 240 s |
| 9 | 1001 | | 378 | see §2.7 | 1001 | (= Y-even by duality) | |

Rigorous consequences: **no certificate with `deg D <= 4`**; for `5 <= deg D <= 8` every admissible `D` is `Σ G_i ∂_i F ∈ J(F)` (a vector-field
tangency locus), except the single Y-even line `D8` at `d = 8`.  A nonzero primitive kernel is forced by dimensions only from `d = 10`
(§1.3), so `D8` is a genuine special relation `A·D8 = 0` in the Jacobian ring.  `D8` was lifted to `Q` by CRT over 3 primes and rational
reconstruction (`e1k8_extract.py`, canonical rref representative modulo `J_8^+`: 414 terms, `D8_e1kernel.json`) and `A·D8 ∈ J` re-verified
at 2 further primes.

### 2.2 Low-degree explicit candidates (batch1; all fail E1, as the sweep predicts)
23 quadrics (coordinate sums of squares, `X_i^2 ± Y_i^2`, `X0X1 ± Y0Y1`, `(X0-X1)^2 ± (Y0-Y1)^2`, `(X0+X1)^2+(Y0+Y1)^2`, `Σ(X_i-X_j)^2`,
`X0Y1 ± X1Y0`, ...) and 16 Y-even cubics (`ΣX_i^3`, `X0X1X2`, Vandermonde, `ΣX_iY_i^2`, `X0Y1Y2+...`, Hesse cubic, ...): **all fail E1** in `< 0.5 s`.
Together with §2.1 this covers *every* `D` of degree `<= 4`, including all companion denominators (`X0^4 - X1^2Y0^2`, `X0^2+X1^2-2Y0^2`,
`X2^2 - Y2^2`, ...) and all their coordinate images.  The degree-6 Fermat/Hesse-web denominators (`family-hesse-web`, `family-block-cube`)
are not in `J(F)` and hence fail E1 too (covered by the `d = 6` row).

### 2.3 Jacobian-ideal sub-families with prescribed singularities (`jideal.py`, `families.py`, `dims_scan.py`)
For `5 <= d <= 8` any admissible `D ≠ D8` is `G(F)`; `Y_D` is the tangency locus of the vector field `G` with `X`.  Imposing non-rational
singularities is linear in `G`: dimensions of `{D ∈ J_d^+ : D ∈ I^k}` modulo `S_{d-6}·F` (coordinates `p,r = X0±X1`, `q,s = Y0±Y1`,
`a = X2`, `b = Y2`; `r1 = X0-X1`, `r2 = X1-X2`, `s1, s2` likewise):

| ideal (singular locus of `Y_D`, transverse type) | d=7 | d=8 | d=9 |
|---|---|---|---|
| `(X0,X1,X2)^3`: triple along the genus-10 curve `Z_Y = X∩{X=0}`, elliptic cone | **21** | 39 | 254 |
| `(Y0,Y1,Y2)^4`: quadruple along `Z_X = X∩{Y=0}`, genus-3 cone | 0 | **31** | 138 |
| `(X0,X1,X2)^4` | 0 | 39 | 63 |
| `(X0,Y0,X1)^3`, `(X0,Y0,X2)^3` / `(X0,Y0,Y1)^3`, `(X0,Y0,Y2)^3` (elliptic cone along a plane sextic curve) | 0 | 37 / 32 | 181 / 170 |
| `(X0,Y0,X1)^4`, `(X0,Y0,Y1)^4` | 0 | 0 | 61 / 51 |
| `(X0,Y0)^2`, `(X2,Y2)^2`, `(X0,X1)^2`, `(Y0,Y1)^2`: double along a sextic surface `X∩P^3` | 9 | 46/46/47/45 | 145/145/148/144 |
| `(X0,Y0)^3`, `(X2,Y2)^3` | 0 | 12 | 61 |
| `(X0,X1,X2,Y0)^4`: quadruple points at the 6 points of `X∩{X=0,Y0=0}` | 17 | 84 | 244 |
| `(X0,X1,Y0,Y1)^4`: quadruple points on the `(X2,Y2)`-line | 10 | 68 | 245 |
| `(r,s)^2` / `(r,s)^3`: double / triple along the `(01)`-fixed surface `X∩{X0=X1,Y0=Y1}` | 17 / **3** | 66 / 24 | 185 / 91 |
| `(p,q,a,b)^4` / `(p,q,a,b)^3` | 0 / 20 | 52 / 112 | 245 / 305 |
| `(r,s,a)^3`, `(r,s,b)^3`: elliptic cone along `X∩{X0=X1,Y0=Y1,X2=0}` resp. `Y2=0` | **6 / 6** | 48 / 48 | 178 / 176 |
| `(p,q,a)^3`, `(p,q,b)^3` | 0 / 2 | 30 / 38 | 174 / 168 |
| `(r1,r2,s1,s2)^4`: quadruple points at the 6 `S_3`-fixed points | **2** | 53 | 245 |
| `(r1,r2,s1,s2)^3` | 22 | 112 | 305 |
| `(r1,r2)^2`, `(s1,s2)^2`: double along `X∩{X0=X1=X2}`, `X∩{Y0=Y1=Y2}` | 9 / 9 | 45 / 45 | 146 / 144 |

Explicit description of the 21-dimensional `(X)^3 ∩ J_7` family (derived by hand, verified): with `{a,b} = {0,1,2} \ {j}`,
`D = Σ_j [α_j (Y_a^2 - Y_b^2) + Q_j(X)] F_{X_j} + (5/39) Σ_j α_j X_j (Y_a F_{Y_b} - Y_b F_{Y_a})`, `α ∈ Q^3`, `Q_j` arbitrary X-quadratics.
The `α`-part uses the specific shape `F_2 = (5/13) Σ_i X_i^2 (Y_a^3 Y_b + Y_a Y_b^3)` of the X-degree-2 part of `F`.

### 2.4 Full fixed-D tests with the fixed engine (v2): Jacobian-ideal denominators of degree 5 and 6 (`batch2b.py`)
Method: E1, then the compressed target-image test (`hlib.target_image_fast`): the full homogeneous system in `(N, Q, A')` with `A'` ranging over
all Y-odd sextics, random row compression to `cols + 40` rows, nullspace, projection to the `A'`-coordinates; `served_dim` is the dimension of
the space of targets served by `D` and `consistent` is whether `A` lies in it.  Baseline `served_dim = 18 = dim J_6^-` (the zero class).
System sizes: `d = 5`: 1386 + 60 unknowns, ~9900 equations; `d = 6`: 2376 + 226 unknowns, ~16800 equations.

| D (degree) | E1 | consistent | served_dim | time (loaded machine) |
|---|---|---|---|---|
| `∂F/∂X_i`, `∂F/∂Y_i` (6 partials, d=5) | yes | **no** | 18 | 33–79 s each |
| `F_X0+F_X1+F_X2`, `F_X0-F_X1`, `F_Y0-F_Y1`, `F_Y0+F_Y1+F_Y2` (d=5) | yes | **no** | 18 | 31–66 s |
| block Euler `Σ X_i F_{X_i}` (d=6) | yes | **no** | 18 | 105 s |
| pair-radial `X0F_X0 + Y0F_Y0` (d=6) | yes | **no** | 18 | 459 s |
| `g_A(F) = (X0∂_0+Y0∂_3-X1∂_1-Y1∂_4)F` (d=6, the "std" radial field) | yes | **no** | 18 | 425 s |
| swap field `(X0∂_1+X1∂_0+Y0∂_4+Y1∂_3)F` (d=6) | yes | **no** | 18 | 476 s |
| rotation field, cyclic-shift field, `Z_i F_{Z_j}` basis (d=6) | yes | see `log.jsonl` | | |

Every one of these vector-field denominators serves **only the zero class** (`served_dim = 18`), i.e. `Y_{G(F)}` supports no odd `(3,1)`-class
at all — consistent with the geometry: polar-type sections have isolated (rational) singularities.  (Session-1 entries with `engine: 1` for the
same `D` reported `rank 1318 / nullity 128` (d=5) and `2359 / 243` (d=6), inconsistent; they are superseded by these v2 results.)

### 2.5 The exceptional degree-8 denominator `D8` (`d8_probe.py`, `run_D8.py`)
Structure of the coset `D8 + J_8^+` (mod p; `dim J_8^+ = 168` inside the 651 Y-even octics):
* `(01)·D8 ≡ +D8 mod J` (the kernel is `(01)`-stable and 1-dimensional, so `D8` is a `(01)`-eigenvector; the sign is `+`, so a `(01)`-invariant
  representative exists); `cyc·D8 ∉ span(J, D8)` (the three-cycle moves `A`).
* None of the natural Fermat-annihilator monomials/forms lies in `span(J, D8)`: `X0^2X1^2Y0^2Y1^2`, its `S_3`-orbit sum, `X0^4Y0^4 ± X1^4Y1^4`,
  `X2^4Y2^4`, `X2^8`, `Y2^8`, `A·(T_0+T_1)`, `A·T_2`, `A·B_2`, `T_0^2 - T_1^2` — all are non-admissible (`∉ J`, `∉ span(J,D8)`).
* `D8 ∉ I^k + J_8` for every singular-locus ideal of §2.3 (`(X)^3, (X)^4, (Y)^3, (Y)^4, (X0,Y0)^2, (X0,Y0)^3, (r,s)^2, (r,s)^3, (X0,X1,X2,Y0)^4,
  (X0,X1,Y0,Y1)^4, (X0,Y0,X1)^3, (r,s,a)^3, (p,q,a,b)^4, (r1,r2,s1,s2)^{3,4}, ...`; the only `True` is `(r,s,a,b)^4 + J = S_8^+`, vacuous).
  So **no representative of `D8` has any of the natural high-order vanishing patterns**: `D8` is not a visibly geometric denominator.
* Full fixed-D test of `D8` (system: `6006 + 1491` unknowns, 40326 equations, compressed): **[RESULT-D8]**

### 2.6 Is the degree-8 exception special to the instance? (`family_kernel.py`, family `Fermat + a·O1 + b·O2 + c·O3`, `d = 8`, Y-even)

| point `(a,b,c)` | dim K_8 | dim J_8^+ | primitive kernel |
|---|---|---|---|
| instance `(5/7, 5/11, 5/13)` | 169 | 168 | **1** |
| Fermat `(0,0,0)` | 370 | 168 | 202 |
| [FAMILY-ROWS] |

### 2.7 Pending / killed runs
[PENDING]

## 3. Witness
None found.  No `witness-candidate.json` was written.

## 4. Mathematical ideas developed

1. **Complete E1 picture via duality (§1.3).**  Only `d = 9` is undetermined.  Any certificate has `deg D >= 8`, and a degree-8 certificate must
   use `D8` (up to `J_8^+`).  If the `D8` test is inconsistent (mod p, confirmed at a second prime) then `deg D >= 9`, and by the prescribed-support
   theorem at `q = 1` all powers `D8^m` are excluded as well.
2. **Vector-field reformulation for Jacobian denominators.**  At `q = 1` the identity forces `N(F) = A·D + F·H` (`deg H = d`) and the certificate
   condition becomes `F | D·div(N) - N(D) - D·H`.  For `D = G(F)` the particular solution `N^0 = A·G`, `H^0 = 0` reduces the problem to
   `F | A·G(F)·div(G) + G(F)·G(A) - A·G^2(F) + (Koszul and `F·M` corrections)`.  All tested fields `G` (partials, block/pair Euler, the `std`
   radial field `g_A`, swap field) serve only the zero class; the served-target computation makes this a statement about the *divisor*
   `Y_{G(F)}`, not just about `A`.
3. **Which Jacobian denominators could work at all.**  Since served targets are a property of `Y_D`, the only Jacobian denominators worth
   testing are those with non-rational singularities (§2.3 table): the 21-dim `(X)^3 ∩ J_7` family (triple along the genus-10 Fermat curve
   `Z_Y`, elliptic-cone transverse type — the exceptional divisor of the resolution is an elliptic fibration over a genus-10 curve, which does carry
   `(2,0)`-forms), the 31-dim `(Y)^4 ∩ J_8` family (genus-3 cone along `Z_X`), the `(r,s)^3 ∩ J_7` family (triple along the `(01)`-fixed sextic
   surface) and the elliptic-cone families `(r,s,a)^3`, `(r,s,b)^3`.  These are exactly the `batch3/batch5/batch4` tests (§2.7).
4. **Fermat ancestry of `D8`.**  At Fermat `Ann(A)_8^+` is 202-dimensional (spanned by monomials in `(X0^2,Y0^2) ∩ (X1^2,Y1^2)`); along the
   deformation it collapses to one line, and that line contains no natural monomial representative (§2.5).  `D8` is therefore a "generic survivor",
   not the deformation of a visible Fermat annihilator such as `X0^2X1^2Y0^2Y1^2`.
5. **The `(01)`-sector reduction at the certificate level** (§1.2) halves the fixed-`D` systems for `(01)`-invariant `D` and is the right setting
   for the degree 9–10 joint search: `D` in the `(01)`-invariant part of `Ann(A)_d` (dimension `k_d^{(01)}`), `N` `(01)`-equivariant.
6. **What would be needed to finish.**  A geometric source of a `(2,0)`-form on a divisor of `X` whose Gysin image is `std`-isotypic.  Natural
   candidates not yet realisable as low-degree `D`: (a) the elliptic fibration over `Z_Y` (needs a `D ∈ (X)^3`, available in `J_7`: 21-dim, being tested),
   (b) the `(01)`-fixed surface `S = X ∩ {X0=X1, Y0=Y1}` (a smooth sextic surface in `P^3`, `p_g = 10`) made non-normal in `Y_D` (`(r,s)^2`-families),
   (c) a degree-10 `D` from the 294-dimensional forced kernel with imposed singularities (the linear conditions of §2.3 can be imposed on `Ann(A)_10`
   directly once a basis is extracted: `294 - codim` survivors).

## 5. Ranked next steps (with estimated compute on an unloaded 40-core machine)
1. **Finish the `D8` fixed-D test at a second prime and over `Q`** (if consistent: exact lift and `verify_witness`; if inconsistent: rigorous
   "`deg D >= 9`" statement). ~5 min.
2. **`d = 9` E1 kernel, both parities** (one computation by duality): `e1kernel.py 9 9`; if the primitive kernel is `<= 20`, extract every element
   over `Q` and run the fixed-D test on each (`d = 9`: `9009 + 2376` unknowns, ~60k equations; ~3–5 min and ~3 GB each). ~1 h total.
3. **Geometric Jacobian families at `d = 7, 8`** (`batch3`, `batch5`, `batch4`): random and symmetry-adapted members of `(X)^3∩J_7`, `(r,s)^3∩J_7`,
   `(r,s,a)^3∩J_7`, `(Y)^4∩J_8`, `(X)^4∩J_8`, with served dimensions.  ~2 min per d=7 test, ~5 min per d=8 test; ~1 h.
4. **`d = 10` E1 kernel (294-dim, Y-even) with imposed singular loci**: compute `Ann(A)_10^+ ∩ I^k` for the ideals of §2.3 (linear algebra on
   `13104`-row systems; ~10 min each) and test survivors (`d = 10`: `13104 + 5796` unknowns, 85k equations, ~15 min and ~10 GB each).  Half the cost with
   the `(01)`-sector reduction.  ~half a day.
5. **Joint bilinear search at `d = 10–11`** inside the `(01)`-invariant E1 kernel (variable projection / Gauss–Newton on the linear subproblem),
   as in the dossier's E3. Days.
6. **Theory**: a sector theorem for the 50-dimensional `S_3×τ`-invariant family (the IVHS algebra commutes with `S_3`, so the generic-target
   argument must be redone isotypic block by block); and the matrix-valued descent needed to use a cover with a reflection for a 2-dimensional
   representation.
