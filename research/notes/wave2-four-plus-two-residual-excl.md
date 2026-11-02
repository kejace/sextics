# Wave 2 — `four-plus-two-residual` — EXCLUSION / RIGOR lens

STATUS: no witness (in progress). BEST RIGOROUS BOUND: **no certificate with deg D <= 5 exists over C** (E1 kernel exact + degree-5 Hessian relaxation exact, §1–2); at degree 6 the 10-parameter eigenvector Jacobian family is excluded exactly (§3.1). Remote box saturated (load 1100+/128): excl-kernels running there; excl-d5/excl-d6 submits waiting; sector-restricted sweeps run locally in a single niced process (chain: symsweep k=6..8 → sum → struct → symsweep k=9..12).

Log: `/tmp/hodge/results/four-plus-two-residual/wave2-excl.log.jsonl` (merged from the remote per-job files `wave2-excl/remote/*.jsonl` by `wave2-excl/w2_pull.py`).
Packet: `/tmp/hodge/results/four-plus-two-residual/wave2-excl/` (w2lib.py, w2_*.py, receipts/, README.md).

## Rigor conventions
* [EXACT] computed over Q with fraction-free integer linear algebra (python-flint fmpz_mat); ranks/kernels are field-independent, so the statements hold over C.
* [RANK] modular-rank argument for a fixed rational D: rank_p(M) <= rank_Q(M); the exact nullity of M is computed from the explicit parametrisation of the homogeneous solutions (N = m z/6 + K, K Koszul, D div K - K(D) = F Q; bijective for deg D <= 8); if it equals the modular nullity and the modular augmented rank exceeds the modular rank, the system is inconsistent over Q, hence over C.
* [MODP] a modular computation only (evidence, not proof).

## 1. E1-kernel structure, re-verified [EXACT] (job excl-kernels)
* H(X0,X1,Y0,Y1) and B(X2,Y2) have Artinian Jacobian rings with the complete-intersection Hilbert functions (1,4,10,20,35,52,68,80,85,80,68,52,35,20,10,4,1,0) and (1,2,3,4,5,4,3,2,1,0): mod 2^31-1, 10^9+7, 998244353 AND exactly over Q. So the partials are regular sequences over Q and modulo these primes (smooth reductions), and R(F) = R(H) (x) R(B) (Thom–Sebastiani, J(F) = J(H) + J(B) in disjoint variables).
* ker(a· : R(H)_k -> R(H)_{k+6}), a = X0^3Y0^3, by (Y0,Y1)-parity (even, odd), dims modulo J(H)_k — identical for the three primes and EXACT over Q:
  k<=4: 0,0 | k=5: 0,0 (kernel = J(H)_5) | k=6: 20,13 | k=7: 30,30 | k=8: 41,34 | k=9: 38,38 | (k=10..16 pending in the log).
* Consequence (exact algebra): {D : A·D in J(F)} in degree d equals sum_k ker_{d-k}(H) (x) C[X2,Y2]_k + J(F)_d. Hence **no E1-admissible denominator of degree <= 4 exists (any Y-parity, any coefficient field), and every E1-admissible denominator of degree 5 lies in J(F)_5 = {F_v = v·grad F : v in C^6}.** (Wave 1 had this modulo one prime for Y-even D only.)
* Exact rational representatives of the kernels modulo J(H)_k were produced two ways: fmpz nullspace (huge heights, up to 10^300 at k=8) and a modular-lattice reduction (`w2lib.small_kernel_reps`: CRT of three 61-bit-prime kernels, LLL, exact verification) giving verified reps of height ~1e11 (k=6 even), ~4e14 (k=6 odd), ~3e8 (k=7); all later fixed-D tests use these exact rational D's.

## 2. Degree 5 — EXCLUDED over C [EXACT] (local run `excl-d5-local`, remote replay `excl-d5` queued)
For D = F_v = v·grad F (the only E1-admissible shape at degree 5, §1) the general solution of E1 is N = A v + (m/6) z + K with K = sum_{i<j} l_ij (F_j d_i − F_i d_j) (Koszul; l_ij linear forms), and E2 reduces to
    F_v A_v − A F_vv + F_v div K − K(F_v) = F Q        (m drops out; A_v = v(A), F_vv = v^T Hess(F) v),
an identity verified exactly on random data (`d5-reduction-identity-check`).
**Hessian relaxation:** replace v⊗v by a free S ∈ Sym²(C^6) (21 unknowns) and v⊗l by a free T ∈ C^6⊗C^{15×6} (540 unknowns), Q ∈ R_4 (126): a homogeneous LINEAR system (2078 equations) solved exactly over Q: rank 435, nullity 252, and its projection W to Sym² has dimension 8. Every certificate of degree 5 gives v⊗v ∈ W with v ≠ 0. The annihilator of W is spanned by 13 quadrics whose Groebner basis is the monomial ideal (v0², v1², v2², v3², v4², v5², v0v1, v0v3, v1v3, v0v4, v3v4, v1v4, v2v5) — Singular dim 0 and an independent sympy Groebner check agree — so W contains no nonzero rank-one tensor (W = span{e_a e_b : a ∈ core, b ∈ {2,5}} ∪ …, only mixed core×binary pairs). **Hence no certificate with deg D = 5 exists over C.**
Positive control: for the Jacobian target Y0·∂F/∂X0 (served by every D) the same code gives W = Sym² (dim 21) and the cone is all of C^6, as it must.
Together with §1: **deg D ≥ 6 for every certificate of this instance** (any parity, any coefficient field ⊂ C). Receipts: `wave2-excl/receipts/excl-d5-local*.{jsonl,json,sing}`.

## 3. Degree 6
### 3.1 Eigenvector Jacobian family — EXCLUDED over C [EXACT] (local run `excl-d6jac-local`)
E1-admissible denominators of degree 6 are D = F_w + c with w ∈ (R_1)^6 and c ∈ ker_6 (33 reps). The sub-family of G-eigenvectors of Jacobian type is the 10-parameter family D = Σ_a w_a F_a, w_a = Σ_{c ∈ class(a)} w_ac z_c, classes {X0,X1}, {Y0,Y1}, {X2}, {Y2} (it contains all weighted Euler polars Σ λ_a z_a F_a). E2 reduces to F_w w(A) + A(F_w div w − w(F_w)) + F_w div K − K(F_w) = F Q with K Koszul (quadratic coefficients, sector-restricted: 48 unknowns) and Q in its sector (40). The Hessian relaxation (S ∈ Sym²(C^10): 55, T = w⊗k: 480, Q: 40; 384 equations) solved exactly over Q has a 26-dimensional projection W; the rank-one locus {w : w⊗w ∈ W} is cut out by 29 quadrics whose minimal associated prime (Singular `minAssGTZ`) is the single line w_00 = w_11 = w_33 = w_44 = w_22 = w_55, w_01 = w_10 = w_34 = w_43 = 0, i.e. D = λ Σ_a z_a F_a = 6λF — not coprime to F. **Hence no denominator of this 10-parameter family has a certificate.** Receipts: `receipts/excl-d6jac-local*.{jsonl,json,sing}`.
### 3.2 Individual degree-6 denominators (full system, two primes) — pending (job excl-d6, queued on the saturated remote); the 33 core-kernel representatives are being excluded exactly by the sector-restricted sweep (§4, local chain).

## 4. Core-only kernel elements, degrees 7–12, sector-restricted full system [EXACT] — pending (job excl-symsweep)
## 5. Sum principle / served sectors — pending (job excl-sum)
## 6. Generic E1-admissible members, degrees 7–10 — pending (job excl-generic)
