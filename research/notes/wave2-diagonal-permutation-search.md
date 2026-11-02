# Wave 2 search report: `diagonal-permutation` (q = 1, A = X0^3 Y0^3 - X1^3 Y1^3)

STATUS: NO WITNESS FOUND (search in progress; this file is rewritten after every batch).

Lane dir: `/tmp/hodge/results/diagonal-permutation/`; log: `wave2-search.log.jsonl` (one line per experiment/batch);
remote outputs mirrored to `/tmp/hodge/remote/diagonal-permutation/` (`w2-*.jsonl`, `w2-*.json`).
Scripts: `w2pt.py` (pointwise variable projection, see §2), `w2pt_search.py` (multi-process driver), `w2_d8_d9.py`
(D8 fixed-D test + d=9 kernel), `w2_kernels.py` (E1 kernels d=10..12).

## 1. Results so far

| experiment | result |
|---|---|
| D8 (unique non-Jacobian E1-kernel element at d=8), full fixed-D system mod p=2147483647, compressed (40326 eq, 7497 unknowns) | **inconsistent** (rank 6703, nullity 794) — second prime pending |

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

(none yet)
