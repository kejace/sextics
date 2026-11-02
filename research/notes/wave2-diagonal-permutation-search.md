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
| Float basis of the degree-9 admissible space (SVD of the 7752×10010 E1 matrix over R, `w2_kerfam.py`) | admissible dim 399 = 378 + 21; singular-value gap 3.2e-3 / 1.5e-14; E1 residual of the kernel lifts ~3e-15 |

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
| j7a | Jacobian `w(F)`, w quadratic (63 params, 3 gauge removed; 420 Koszul unknowns; r = 543 points) | 7 | all | 368 / 3200 | 9.9e-2 | 0.22 | running |
| j8a | Jacobian `w(F)` (w cubic, 168) ⊕ D8 (169 params, 12 gauge; 954 Koszul; r = 1183) | 8 | all | started | – | – | running |

Relative residual = `min_c |ρ|² / |B_D(N^0)|²` at the normalized point; values ~0.1–0.3 are ordinary local minima of an
infeasible system.  **Control** (`w2_control.py`, accepted witness `fermat-pair-rational-d4-strict`, d=4): with the exact E1
lift the pointwise residual is `f = 3e-30` (here `B_D(u)` itself vanishes on `V(F)`, so the *relative* value is meaningless
and an absolute criterion `f / median|B_0|²` was added to the driver); perturbing the lift by `(1/50)·monomial` or `D` by
`(1/20)·monomial` gives relative residuals 0.90 / 0.92.

## 4. Coordination (other agents on this lane, same remote lane dir)
* `wave2-excl/` (exclusion/rigor agent): Hessian relaxation by msolve (d=5,6 all sectors: `V(I)={0}` rigorous mod p with
  good reduction), E1 kernels at two primes (`e1kF`: d=8 odd kernel 0, d=9 even kernel 21 confirmed at p=1000000007),
  D8 exact/rank-certified, structured families (`w2_families.py`), individual degree-9 kernel tests (`w2_d9.py`).
* `w2search/` (second search agent): `fastvp.py` variable projection on Jacobian families from d=5 upward, structured
  candidate E1 filter (`cands-structured.jsonl`).
* `kern10/kern11/kern12` jobs: `w2_kernels.py` (E1 kernel bases d=10,11,12 mod p) launched by the search agent at 00:20 UTC;
  my duplicate `k10` was killed.  My structured sweep (`w2_struct.py`) and Hessian numerics (`w2_hess.py`) are therefore
  NOT run (covered by the other agents); this agent concentrates on the pointwise joint searches.
