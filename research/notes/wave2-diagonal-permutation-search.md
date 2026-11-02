# Wave 2 search report: `diagonal-permutation` (q = 1, A = X0^3 Y0^3 − X1^3 Y1^3)

STATUS: NO WITNESS FOUND.  No candidate with residual below 1e-3 (relative) or 5e-5 (absolute-scaled) appeared in any
family; every low-residual value is reproduced by a random control target.  Rigorous by-products: D8 excluded (two primes),
the degree-9 E1 kernel (21 + 21 classes) computed, lifted to R and explained by the transposition symmetry.

Lane dir `/tmp/hodge/results/diagonal-permutation/` (scripts `w2pt.py`, `w2pt_search.py`, `w2sym.py`, `w2_kerfam.py`,
`w2_refine.py`, `w2_k9_exact.py`, `w2_d8_d9.py`, `w2_kernels.py`, `w2_forced_kernels.py`, `w2_control.py`, `w2_merge.py`,
`w2_report_stats.py`); log `wave2-search.log.jsonl`; remote outputs mirrored to `/tmp/hodge/remote/diagonal-permutation/`
(`w2-*.jsonl`, `w2-kernel-*.json`, `w2-refine-*`, per-restart files `w2-pt-d<d>-<tag>.jsonl`; batch statistics in
`w2-batch-stats.json`).  Jobs `search-diagonal-permutation-{j9eps,j9ops,j9x3s,calAs}` were still running on the remote when
this report was written (04:40 UTC); `python3 -m hodge.rjob pull --lane diagonal-permutation` + `w2_report_stats.py` refresh the numbers.

## 1. Summary of results

1. **D8** (the unique non-Jacobian E1-kernel element of degree 8, wave 1): full fixed-D system (40326 equations, 7497 unknowns)
   **inconsistent at p = 2147483647 and p = 2147483629** (rank 6703, nullity 794 both); served-target space = `J_6^-` only
   (dimension 18).  Hence every degree-≤8 denominator is Jacobian type `w(F)`.
2. **Degree 9 is the first non-Jacobian degree**: `ker([A]: R_9^± → R_15^∓)` has dimension **21 for both Y-parities**
   (mult_kernel at 4 primes; duality cross-check), one degree below the dimension-count bound (d = 10).  Explanation: with the
   (01)-equivariant Koszul series the σ-invariant Y-even part of `R_9` has dimension 322 and its target (σ-anti, Y-odd part of `R_15`)
   301, so a 21-dimensional kernel is forced by the transposition symmetry; the computed kernel is exactly this forced part
   (entirely σ-invariant; the σ-anti part 301 → 322 is injective).  The same count gives 170 + 124 = 294 at d = 10 (Y-even) and
   152 + 109 + 1 = 262 (Y-odd), matching duality.  Structure (probes mod p): 11 of the 21 Y-even classes have representatives in
   `(X0,X1,X2)^3` (triple along the Fermat curve `X ∩ {X = 0}`: elliptic-cone transverse type), 4 in `(Y0,Y1,Y2)^4` (quadruple along
   `X ∩ {Y = 0}`: genus-3 cone); none in `(X)^4`, `(X_i,Y_i)^2`, `(r,s)^k`, `(p,q)^2`; the other probed ideals meet the kernel only in
   the dimension-forced amount.  A float basis with E1 lifts was built by SVD (gap 3e-3 / 1.5e-14; E1 residual ~3e-15).
3. **Joint searches (pointwise variable projection, §3)**: d = 7 Jacobian family 558 restarts (floor 0.10–0.26); d = 8 Jacobian ⊕ D8
   family: no restart finished under remote load > 1000; d = 9: 84 restarts over the full admissible families (Y-even, Y-odd) and
   the `(X)^3`, `(Y)^4` sub-families in the (01)-invariant sector, at two equation surpluses, plus two long refinements.  Best values
   for the real target: relative 1.9e-3 / absolute-scaled 1.4e-5 (small surplus) and 5.6e-3 / 5.5e-5 (doubled surplus).
   **Controls**: a random σ-anti-invariant target A′ gives 3.3e-3 / 1.2e-5 (small surplus, 6 restarts) and 3.1e-2 / 6.1e-5 (doubled
   surplus, 3 restarts); tripling the point count on the best start raised its converged minimum from f ≈ 17 (still decreasing) to
   f = 156.1 (gradient 0.02).  So the degree-9 residual floor is a property of the compressed system, not of A, and no basin of a true
   zero (f ~ 1e-15) was entered.  Degree 10 (admissible dimension 1044) was not searched jointly: each evaluation costs seconds
   under the load, and the (01)-sector version needs ~20 GB per worker; the kernel bases at d = 10, 11, 12 were computed by the
   sibling `kern10/11/12` jobs.

## 2. E1-admissible denominator spaces per degree (rank/nullity statistics)

`k_d^±` = dim ker([A]) modulo `J(F)_d`; `Adm_d^± = k_d^± + dim J_d^±`.  Mod p = 2147483647; zero kernels are rigorous over Q,
nonzero ones confirmed at a second prime (d = 8, 9) or by duality (d ≥ 10).

| d | dim S_d^+ | dim J_d^+ | R_d^+ | k_d^+ | Adm_d^+ | dim S_d^- | dim J_d^- | R_d^- | k_d^- | Adm_d^- |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 66 | 0 | 66 | 0 | 0 | 60 | 0 | 60 | 0 | 0 |
| 5 | 126 | 3 | 123 | 0 | 3 | 126 | 3 | 123 | 0 | 3 |
| 6 | 236 | 18 | 218 | 0 | 18 | 226 | 18 | 208 | 0 | 18 |
| 7 | 396 | 63 | 333 | 0 | 63 | 396 | 63 | 333 | 0 | 63 |
| 8 | 651 | 168 | 483 | **1** (D8, excluded) | 169 | 636 | 168 | 468 | 0 | 168 |
| 9 | 1001 | 378 | 623 | **21** (all (01)-inv.) | 399 | 1001 | 378 | 623 | **21** | 399 |
| 10 | 1512 | 750 | 762 | 294 | 1044 | 1491 | 747 | 744 | 262 | 1009 |
| 11 | 2184 | 1341 | 843 | 510 | 1851 | 2184 | 1341 | 843 | 510 | 1851 |
| 12 | 3108 | 2223 | 885 | 677 | 2900 | 3080 | 2214 | 866 | 648 | 2862 |

(01)-sector dimensions of `R_d` (Y-even σ+, Y-even σ−, Y-odd σ+, Y-odd σ−): d=9: 322/301/322/301; d=10: 398/364/380/364;
d=15: 322/301/322/301; d=16: 255/228/240/228 (`w2_forced_kernels.py`).  S3 × τ isotypic multiplicities (triv/sign/std) of `R_d^+`:
d=9: 115/94/207, d=10: 144/110/254, d=11: 153/128/281, d=12: 167/130/294; of `R_d^-`: d=6: 39/31/69, d=15: 115/94/207, d=16: 84/72/156,
d=17: 63/48/111, d=18: 39/31/69.  For S3-invariant denominators the E1 kernel is forced only from d = 11 (≥ 42) and d = 12 (≥ 98).

## 3. Method (pointwise variable projection on V(F))

1. *Reduced E2.*  With `N(F) = A D + F H` and `F_0..F_5` a regular sequence, every E1 lift is `N = N^0 + (H/6)E + K_c` (`N^0(F) = A D`,
   `E` Euler field, `K_c` Koszul field with skew `c_ij` of degree d − 4); the H-term cancels in the cleared identity, so
   E2 ⇔ `B_D(N^0) + T_D(c) ≡ 0 (mod F)`, `B_D(N) = D div N − N(D)`, `T_D(c) = Σ_{i<j} D[(∂_i c_ij)F_j − (∂_j c_ij)F_i] − c_ij(D_iF_j − D_jF_i)`.
   For Jacobian members `D = m F_i`: `N^0 = A m ∂_i`; for kernel members the float lift from the E1 SVD.
2. *Points instead of rows.*  `P ≡ 0 (mod F)` iff P vanishes on V(F) (F irreducible).  r random complex points of V(F) (five unit-modulus
   coordinates, sixth a root of the univariate sextic, |F(z)| ~ 1e-14) give 2r real equations; with `2r > dim a + dim c` false zeros are
   impossible for generic points.  Everything is pointwise (D, ∂D, N^0, div N^0, F_j, c, ∂c), no large matrix is stored; one dense least
   squares `(2r) × nK` per evaluation (gelsy), exact variable-projection gradient, L-BFGS on the unit sphere of the family coordinates
   orthogonal to the spurious directions (`D ∈ F·S_{d−6}`, `D = 0`).  Controls: the accepted witness of `fermat-pair-rational-d4-strict`
   gives f = 3e-30; perturbing its lift or D by 1/50 … 1/20 of a monomial gives relative residuals 0.90 / 0.92.
3. *(01)-sector.*  σ = (01) fixes F and negates A; for σ-invariant D one may take N σ-anti-equivariant, i.e. w equivariant and c in the
   −1 eigenspace of `(σ·c)_{ij}(z) = c_{σi σj}(σz)` (identity `σ_*(K_c) = K_{σ·c}` verified exactly); this halves parameters and
   unknowns.  It is a sub-search: a general degree-9 denominator is `D_inv + w_−(F)` with a σ-anti-invariant Jacobian part.
4. *Equation surplus matters.*  With `2r` only ~20 % above `dim a + dim c` the compressed landscape is a nearly flat valley whose floor
   (relative ~1e-3) carries no information (shown by the control target and by the point-count test); with `2r ≈ 2(dim a + dim c)` the
   floor is O(1e-2 … 0.2) and converges in 75–400 iterations, as at degree 7.  The Hessian relaxation (same identity mod J:
   `A·wᵀHess(F)w ∈ J`) was pursued rigorously by the exclusion agent (msolve).

## 4. Batches

| batch | family | d | sector | points r | restarts | best rel. / abs.-scaled residual | median rel. | notes |
|---|---|---|---|---|---|---|---|---|
| j7a | Jacobian `w(F)`, w quadratic (63 params, 3 gauge; 420 Koszul) | 7 | all | 543 | 558 | 9.9e-2 / – | 0.215 | converged in 52–214 iterations; stopped to free compute |
| j8a/j8b | Jacobian (168) ⊕ D8, 954 Koszul | 8 | all | 1183 / 662 | 0 | – | – | no restart finished under load 900–1900; killed |
| j9ep2 | J_9^+ ⊕ 21 kernel members (399 → 201 σ-inv. params; 950 Koszul) | 9 | (01)-inv | 684 | 8 | 3.7e-3 / 2.6e-5 | 1.2e-2 | 300-iteration cap; floor |
| j9op2 | Y-odd version | 9 | (01)-inv | 684 | 15 | **1.9e-3 / 1.4e-5** | 3.4e-3 | 300-iteration cap; floor |
| j9x3 | admissible ∩ (X)^3 (141 params) | 9 | (01)-inv | 684 | 10 | 3.9e-3 / 2.3e-5 | 7.3e-3 | 300-iteration cap |
| calA | **control** A′, Y-odd, as j9op2 | 9 | (01)-inv | 684 | 6 | 3.3e-3 / 1.2e-5 | 6.5e-3 | same floor as A |
| ref9op | refinement of j9op2/w4/r0 (varpro, 1500 it.) | 9 | (01)-inv | 684 | 1 | f: 1.8e5 → 153 → 64 → 27.6 → 12.5 at 26/101/251/501 it. | – | slow valley, gradient 3.7; stopped |
| cal9 | same start, r = 984 | 9 | (01)-inv | 984 | 1 | converged f = 156.1 (gradient 0.02) | – | floor rises 6× with 300 more points |
| j9eps | full Y-even family, doubled surplus, 400 it. | 9 | (01)-inv | 1260 | 9 | 5.6e-3 / 7.1e-5 | 1.1e-2 | ~370 iterations, ~60 min/restart; running |
| j9ops | full Y-odd family, doubled surplus | 9 | (01)-inv | 1260 | 8 | 6.7e-3 / 5.5e-5 | 1.0e-2 | running |
| j9x3s | admissible ∩ (X)^3, doubled surplus | 9 | (01)-inv | 1260 | 15 | 1.6e-2 / 1.3e-4 | 1.8e-2 | converged (~215 it.); running |
| j9y4s | admissible ∩ (Y)^4, doubled surplus | 9 | (01)-inv | 1260 | 19 | 4.6e-3 / 6.0e-3 | 0.175 | 18/19 at the 0.17 floor; stopped |
| calAs | **control** A′, Y-odd, doubled surplus | 9 | (01)-inv | 1260 | 3 | 3.1e-2 / 6.1e-5 | 3.2e-2 | same absolute-scaled floor as A; running |

Relative residual = `min_c|ρ|² / |B_D(N^0)|²` at the final point; absolute-scaled = `min_c|ρ|²` divided by the median of `|B_D(N^0)|²`
over random directions (a true zero gives ~1e-15 for both).  Histograms per batch are in `w2-batch-stats.json`.

## 5. Coordination (other agents on the same lane dir)
* `wave2-excl/` (exclusion/rigor agent): Hessian relaxation by msolve (d = 5, 6: `V(I) = {0}`, d = 7 running), E1 kernels at two primes
  (d = 9 even kernel 21 confirmed at p = 1000000007), D8 exact/rank certification, structured families, degree-9 kernel line tests.
* `w2search/` (second search agent): `fastvp.py` variable projection from d = 5 upward with structured E1-filtered candidates.
* `kern10/kern11/kern12` (E1 kernel bases at d = 10–12) launched by the search agent from `w2_kernels.py`; my duplicate was killed.
  My structured sweep (`w2_struct.py`) and Hessian numerics (`w2_hess.py`) were therefore not run (covered elsewhere).
* Coordination note left in `w2-NOTE-pointwise-search-agent.txt`.

## 6. What remains / suggested next steps
1. Unrestricted degree-9 searches (σ-anti-invariant Jacobian part included; ~4× the cost per evaluation) and the σ-anti sector of the
   Jacobian family (`--sector -1`), with doubled surplus, once the remote is less loaded (the code supports them: `w2pt_search.py --d 9
   --extra kerfam [--sector -1]`).
2. Degree 10 in the (01)-sector for the equivalent target B = T0 + T1 − 2T2 (B invariant, N equivariant) would need the float kernel
   family (`w2_kerfam.py 10 0`, SVD 10152 × 14616) and ~20 GB per worker; feasible on an idle machine, not under this load.
3. Exact rational basis of the degree-9 kernel (CRT with ≥ 8 primes; heights ≥ 25 digits) and the geometry of the `(X)^3` classes
   (transverse cubics along `Z_Y`): a candidate geometric source of the (2,0)-forms the certificate needs.
4. The `mod J` relaxation for non-Jacobian D is one scalar cubic condition `res(A·N^0(D)(D)) = 0` per `A' ∈ Ann_6(D)`; worth adding as a
   cheap filter in a future exact search.
