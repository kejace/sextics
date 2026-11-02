# Wave 2 search report: `diagonal-permutation` (q = 1, A = X0^3 Y0^3 - X1^3 Y1^3)

STATUS: NO WITNESS FOUND (search in progress; this file is rewritten after every batch).

Lane dir: `/tmp/hodge/results/diagonal-permutation/`; log: `wave2-search.log.jsonl` (one line per experiment/batch);
remote outputs mirrored to `/tmp/hodge/remote/diagonal-permutation/` (`w2-*.jsonl`, `w2-*.json`).
Scripts: `w2pt.py` (pointwise variable projection, see §2), `w2pt_search.py` (multi-process driver), `w2_d8_d9.py`
(D8 fixed-D test + d=9 kernel), `w2_kernels.py` (E1 kernels d=10..12).

## 1. Results so far

| experiment | result |
|---|---|
| D8 (unique non-Jacobian E1-kernel element at d=8), full fixed-D system (40326 eq, 7497 unknowns), compressed, mod p=2147483647 and p=2147483629 | **inconsistent at both primes** (rank 6703, nullity 794 both) ⇒ the only non-Jacobian denominator of degree ≤ 8 is excluded; every degree-8 denominator is Jacobian-type `w(F)` (w cubic) |
| D8 served-target space (compressed homogeneous system with 226 target columns) | served dim = 18 = `J_6^-` (only the zero class); A not served |
| E1 kernel modulo `J(F)` at d=9, Y-even (`mult_kernel`, rows 7752, cols 10010) | **dim 21** (kernel 399 = 378 + 21): the first non-Jacobian denominators appear at **d = 9**, not 10 (duality ⇒ Y-odd kernel also 21) |
| E1 kernel modulo `J(F)` at d=9, Y-odd | **dim 21** (computed; matches duality) |
| Float basis of the degree-9 admissible space (SVD of the 7752×10010 E1 matrix over R, `w2_kerfam.py`), Y-even | admissible dim 399 = 378 + 21; singular-value gap 3.2e-3 / 1.5e-14; E1 residual of the kernel lifts ~3e-15 |
| same, Y-odd | admissible dim 399 = 378 + 21; gap 4.0e-3 / 1.7e-14 |
| **Structure of the degree-9 Y-even kernel** (`w2_k9_exact.py`, 4 primes) | dim 21 at all 4 primes; the kernel is **(01)-stable and entirely (01)-invariant** (invariant part 21, anti-invariant part 0); 11 of the 21 classes have representatives in `(X0,X1,X2)^3` (triple along the Fermat curve `Z_Y = X∩{X=0}`: elliptic-cone transverse type), 4 in `(Y0,Y1,Y2)^4` (quadruple along `Z_X`, genus-3 cone), (these two are *significant*: the ambient spaces `((X)^3+J)/J`, `((Y)^4+J)/J` have dims 602, 592 < 623, so no intersection is forced); the counts 17 in `(X0−X1,Y0−Y1,X2)^3`, 12 in `(X0−X1,Y0−Y1,Y2)^3`, 13 in `(X0,X1,X2,Y0)^4` are within 1 of the dimension-forced lower bounds (ambient dims 618, 614, 615) and carry no information; none in `(X)^4`, `(X_i,Y_i)^2`, `(r,s)^2`, `(p,q)^2`; rational reconstruction of the canonical rref basis needs more primes (heights ≥ 25 digits; 1/21 reconstructed) |
| d=7 Jacobian family, pointwise varpro, 553 unrestricted restarts (j7a) | best relative residual 9.9e-2, median 0.21: no candidate; stopped to free compute (the excl agent runs the rigorous d=7 Hessian relaxation by msolve) |

## 2. Method for the joint search (reduced E2, pointwise)

For `D` in a linear family with E1 lifts `N^k` (`N^k(F) = A D_k`), every E1 lift is `N = N^0(a) + (H/6)E + K_c`
(Euler field `E`, Koszul field `K_c` with skew `c_ij` of degree `d-4`); the `H`-term cancels in E2 identically, so E2 is
`B_D(N^0) + T_D(c) ≡ 0 (mod F)`, `B_D(N) = D div N − N(D)`, `T_D(c) = Σ_{i<j} D[(∂_i c_ij)F_j − (∂_j c_ij)F_i] − c_ij(D_iF_j − D_jF_i)`.
`mod F` is tested by evaluation at `r` random complex points of `V(F)` (F irreducible), with `r > dim a + dim c`, so the
compressed system has no false zeros; the reduced functional `f(a) = min_c |ρ|²` is minimized on the unit sphere
(orthogonal to the spurious directions `D ∈ F·S_{d−6}`, `D = 0`) by L-BFGS with the exact variable-projection gradient.
For Jacobian-type `D = w(F)` the lift is `N^0 = A w`.  Reduction modulo `J(F)` of the same identity gives the Hessian
relaxation `A·wᵀHess(F)w ∈ J(F)` (necessary; 208 quadrics in 63 unknowns at d=7, 60 in 168 at d=8, 9 in 378 at d=9).

## 3. Batches (filled as they finish)

| batch | family | d | sector | restarts | best rel. residual | median | status |
|---|---|---|---|---|---|---|---|
| j7a | Jacobian `w(F)`, w quadratic (63 params, 3 gauge removed; 420 Koszul unknowns; r = 543 points) | 7 | all | 498 / 3200 | 9.9e-2 | 0.21 | running |
| j8a | Jacobian `w(F)` (w cubic, 168) ⊕ D8 (169 params, 12 gauge; 954 Koszul; r = 1183) | 8 | all | 0 (killed: too slow under remote load ~900; no restart finished in 25 min) | – | – | replaced by j8b |
| j8b | same family, r = 662 points (2r > dim a + dim c), maxiter 300 | 8 | all | started 01:10 UTC | – | – | running |
| j9ep2 | **full admissible family at d=9, Y-even**: J_9^+ (w quartic, 378) ⊕ 21 kernel members (float E1 lifts) = 399 params, restricted to (01)-invariant D (201 params after gauge) with (01)-anti-equivariant Koszul corrections (950 of 1890); r = 684 points | 9 | (01)-inv | started 01:17 UTC (8×150) | – | – | running (first attempt j9ep crashed on a float-member bug, fixed) |
| j9op2 | same for **Y-odd D** (J_9^- ⊕ 21 odd kernel members), (01)-invariant sector | 9 | (01)-inv | started 01:21 UTC (8×150) | – | – | running |
| (planned) j9y4 | admissible ∩ `(Y0,Y1,Y2)^4` at d=9 (the 4 kernel classes with genus-3-cone singularities along `Z_X` plus `(Y)^4 ∩ J_9`), (01)-invariant sector | 9 | (01)-inv, vanish Y4 | queued | | | |
| (planned) j9x3 | admissible ∩ `(X0,X1,X2)^3` at d=9 (the 11 kernel classes with elliptic-cone singularities along `Z_Y` plus `(X)^3 ∩ J_9`), (01)-invariant sector | 9 | (01)-inv, vanish X3 | queued | | | |

Relative residual = `min_c |ρ|² / |B_D(N^0)|²` at the normalized point; values ~0.1–0.3 are ordinary local minima of an
infeasible system.  **Control** (`w2_control.py`, accepted witness `fermat-pair-rational-d4-strict`, d=4): with the exact E1
lift the pointwise residual is `f = 3e-30` (here `B_D(u)` itself vanishes on `V(F)`, so the *relative* value is meaningless
and an absolute criterion `f / median|B_0|²` was added to the driver); perturbing the lift by `(1/50)·monomial` or `D` by
`(1/20)·monomial` gives relative residuals 0.90 / 0.92.

## 3b. E1-admissible denominator spaces per degree (rank/nullity statistics)

`k_d^±` = dim of the E1 kernel modulo `J(F)_d` (Y-even / Y-odd D); `Adm_d^± = k_d^± + dim J_d^±` is the dimension of the
E1-admissible denominator space inside `S_d^±`.  Computed mod p = 2147483647 (zero kernels are rigorous over Q; nonzero
dimensions are confirmed at a second prime where stated); `d ≥ 10` by Gorenstein duality `k_d^ε − k_{18−d}^{1−ε} = R_d^ε − R_{18−d}^{1−ε}`
(cross-checked against the direct computation of the `kern10/11/12` jobs when they finish).

| d | dim S_d^+ | dim J_d^+ | R_d^+ | k_d^+ | Adm_d^+ | dim S_d^- | dim J_d^- | R_d^- | k_d^- | Adm_d^- |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 66 | 0 | 66 | 0 | 0 | 60 | 0 | 60 | 0 | 0 |
| 5 | 126 | 3 | 123 | 0 | 3 | 126 | 3 | 123 | 0 | 3 |
| 6 | 236 | 18 | 218 | 0 | 18 | 226 | 18 | 208 | 0 | 18 |
| 7 | 396 | 63 | 333 | 0 | 63 | 396 | 63 | 333 | 0 | 63 |
| 8 | 651 | 168 | 483 | **1** (D8, excluded by the full fixed-D test) | 169 | 636 | 168 | 468 | 0 | 168 |
| 9 | 1001 | 378 | 623 | **21** (2 primes; all (01)-invariant) | 399 | 1001 | 378 | 623 | **21** | 399 |
| 10 | 1512 | 750 | 762 | 294 | 1044 | 1491 | 747 | 744 | 262 | 1009 |
| 11 | 2184 | 1341 | 843 | 510 | 1851 | 2184 | 1341 | 843 | 510 | 1851 |
| 12 | 3108 | 2223 | 885 | 677 | 2900 | 3080 | 2214 | 866 | 648 | 2862 |

**Why 21, and why (01)-invariant.**  With the (01)-equivariant Koszul series (σ has eigenvalues (1,1,1,1,−1,−1) on C^6,
στ has (1,−1,1,−1,1,−1)) the Jacobian ring splits as `R_d^{±,σ±}` with dimensions (Y-even σ+, Y-even σ−, Y-odd σ+, Y-odd σ−) =
(322, 301, 322, 301) at d = 9 and (398, 364, 380, 364) at d = 10, (255, 228, 240, 228) at d = 16, (174, 159, 174, 159) at d = 17.
Since A is Y-odd and σ-anti-invariant, `[A]` maps `R_9^{+,σ+}` (322) into `R_15^{−,σ−}` (301): a kernel of dimension ≥ 21 is
**forced by the transposition symmetry alone**, and the computed kernel is exactly this forced part (21, all σ-invariant); the
σ-anti-invariant part (301 → 322) is injective.  At d = 10 the same count gives 170 + 124 = 294 (Y-even) and 152 + 109 (+1 from D8)
= 262 (Y-odd), reproducing the duality values.  So the degree-9 kernel is the first place where the symmetry forces E1-admissible
non-Jacobian denominators; it is not a coincidence special to the coefficients 5/7, 5/11, 5/13.

Jacobian-type degrees (every admissible D is `w(F)`): d ≤ 7 (both parities) and d = 8 (Y-odd); d = 8 Y-even has the single extra
line D8 (now excluded).  First genuinely non-Jacobian degree: **d = 9** (21 + 21 classes), one lower than the dimension count forces.
Full fixed-D system of D8: 40326 equations, 7497 unknowns (N 6006, Q 1491), rank 6703, nullity 794, inconsistent at two primes.

## 3c. Mathematics of the method (for the record)

1. *Reduced E2.*  With `N(F) = A D + F H` and `F_0..F_5` a regular sequence, `N = N^0 + (H/6)E + K_c` for any fixed E1 lift
   `N^0` (`N^0(F) = A D`); substituting into the cleared identity gives `F·(B_D(N^0) + T_D(c) − F Q) = 0`, the H-terms cancel
   identically.  Unknowns: skew `c_ij` of degree `d − 4` and `Q`; equations in degree `2d + 6`.
2. *Points instead of rows.*  `P ≡ 0 (mod F)` iff `P` vanishes on `V(F)` (F irreducible).  Evaluating at `r` random points
   `z_t ∈ V(F)` (five coordinates on the unit circle, the sixth a root of the univariate sextic, |F(z)| ~ 1e-14) turns the
   identity into `r` complex (2r real) equations; the dimension count `2r > dim a + dim c` makes false zeros impossible for
   generic points, and every quantity entering the system (D, ∂_iD, N^0, div N^0, F_j, c, ∂c) is a pointwise value, so no
   large matrix is ever assembled: one dense least squares `(2r) × nK` per evaluation, exact variable-projection gradient.
3. *(01)-sector.*  σ = (01) fixes F, negates A.  If σD = D and (N, Q) certifies A then `(N − σ_*N)/2` does too, so one may take
   `N` σ-anti-equivariant: `w` σ-equivariant and `c` in the −1 eigenspace of `(σ·c)_{ij}(z) = c_{σi σj}(σz)` (identity
   `σ_*(K_c) = K_{σ·c}` verified exactly).  This halves both the parameter and the unknown count.  Caveat: the σ-invariant
   sector is a genuine sub-search.  Every non-Jacobian admissible degree-9 *class* is σ-invariant modulo J, but a denominator
   itself is `D = D_inv + w_−(F)` with `D_inv` σ-invariant admissible and `w_−(F)` a σ-anti-invariant Jacobian element; the
   sector search fixes `w_− = 0`.  The unrestricted degree-9 searches (`j9e`, `j9o`, ~4× the cost per evaluation) cover the
   full admissible family.
4. *Hessian relaxation* (reduction of the same identity modulo J instead of F): `A·wᵀHess(F)w ∈ J(F)` for Jacobian-type D; used
   rigorously (msolve) by the exclusion agent.

## 4. Coordination (other agents on this lane, same remote lane dir)
* `wave2-excl/` (exclusion/rigor agent): Hessian relaxation by msolve (d=5,6 all sectors: `V(I)={0}` rigorous mod p with
  good reduction), E1 kernels at two primes (`e1kF`: d=8 odd kernel 0, d=9 even kernel 21 confirmed at p=1000000007),
  D8 exact/rank-certified, structured families (`w2_families.py`), individual degree-9 kernel tests (`w2_d9.py`).
* `w2search/` (second search agent): `fastvp.py` variable projection on Jacobian families from d=5 upward, structured
  candidate E1 filter (`cands-structured.jsonl`).
* `kern10/kern11/kern12` jobs: `w2_kernels.py` (E1 kernel bases d=10,11,12 mod p) launched by the search agent at 00:20 UTC;
  my duplicate `k10` was killed.  My structured sweep (`w2_struct.py`) and Hessian numerics (`w2_hess.py`) are therefore
  NOT run (covered by the other agents); this agent concentrates on the pointwise joint searches.
