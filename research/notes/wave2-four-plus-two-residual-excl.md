# Wave 2 — `four-plus-two-residual` — EXCLUSION / RIGOR lens

STATUS: no witness. BEST RIGOROUS BOUND: **no certificate with deg D <= 5 exists over C** [EXACT] (E1 kernel exact + degree-5 Hessian relaxation); at degree 6 every G-eigenvector Jacobian-type denominator (13 character classes) is excluded [EXACT] (§3.1) and every tested core-kernel representative is excluded [RANK] (§4; the local sweep keeps running in the background via `wave2-excl/w2_local_chain.sh`, merge its outputs with `python3 wave2-excl/w2_pull.py --nopull`). The remote box was saturated (load ~1050/128) throughout this session: `excl-kernels` ran there (receipts pulled); `excl-d5` (replay), `excl-d6` (full unrestricted system), `excl-generic`, `excl-struct` could NOT be submitted (capacity guard) — scripts are in the packet, replay commands in README.md.

Log: `/tmp/hodge/results/four-plus-two-residual/wave2-excl.log.jsonl` (merged from the remote per-job files `wave2-excl/remote/*.jsonl` by `wave2-excl/w2_pull.py`).
Packet: `/tmp/hodge/results/four-plus-two-residual/wave2-excl/` (w2lib.py, w2_*.py, receipts/, README.md).

## Rigor conventions
* [EXACT] computed over Q with fraction-free integer linear algebra (python-flint fmpz_mat); ranks/kernels are field-independent, so the statements hold over C.
* [RANK] modular-rank argument for a fixed rational D: rank_p(M) <= rank_Q(M); the exact nullity of M is computed from the explicit parametrisation of the homogeneous solutions (N = m z/6 + K, K Koszul, D div K - K(D) = F Q; bijective for deg D <= 8); if it equals the modular nullity and the modular augmented rank exceeds the modular rank, the system is inconsistent over Q, hence over C.
* [MODP] a modular computation only (evidence, not proof).

## 1. E1-kernel structure, re-verified [EXACT] (job excl-kernels)
* H(X0,X1,Y0,Y1) and B(X2,Y2) have Artinian Jacobian rings with the complete-intersection Hilbert functions (1,4,10,20,35,52,68,80,85,80,68,52,35,20,10,4,1,0) and (1,2,3,4,5,4,3,2,1,0): mod 2^31-1, 10^9+7, 998244353 AND exactly over Q. So the partials are regular sequences over Q and modulo these primes (smooth reductions), and R(F) = R(H) (x) R(B) (Thom–Sebastiani, J(F) = J(H) + J(B) in disjoint variables).
* ker(a· : R(H)_k -> R(H)_{k+6}), a = X0^3Y0^3, by (Y0,Y1)-parity (even, odd), dims modulo J(H)_k — identical for the three primes and EXACT over Q:
  k<=4: 0,0 | k=5: 0,0 (kernel = J(H)_5) | k=6: 20,13 | k=7: 30,30 | k=8: 41,34 | k=9: 38,38 | k=10: 36,31 | k=11: 26,26 | k=12: 19,16 (for k ≥ 11 the kernel is all of R(H)_k since R(H)_{k+6} = 0 for k+6 > 16; k=12 odd and the three-prime table for k ≤ 16 in the log).
* Consequence (exact algebra): {D : A·D in J(F)} in degree d equals sum_k ker_{d-k}(H) (x) C[X2,Y2]_k + J(F)_d. Hence **no E1-admissible denominator of degree <= 4 exists (any Y-parity, any coefficient field), and every E1-admissible denominator of degree 5 lies in J(F)_5 = {F_v = v·grad F : v in C^6}.** (Wave 1 had this modulo one prime for Y-even D only.)
* Dimensions of the full E1-admissible family {D : A·D ∈ J(F)} by Y-parity (even/odd), mod 10^9+7 (J-part) + exact kernel part: d=6: 38/31 (J 18 + ker 20 | J 18 + ker 13); d=7: 126/126 (63+63); d=8: 322/308 (168+154 | 168+140); within 236/226, 396/396, 651/636 monomials.
* Exact rational representatives of the kernels modulo J(H)_k were produced two ways: fmpz nullspace (huge heights, up to 10^300 at k=8) and a modular-lattice reduction (`w2lib.small_kernel_reps`: CRT of three 61-bit-prime kernels, LLL, exact verification) giving verified reps of height ~1e11 (k=6 even), ~4e14 (k=6 odd), ~3e8 (k=7); all later fixed-D tests use these exact rational D's.

## 2. Degree 5 — EXCLUDED over C [EXACT] (local run `excl-d5-local`, remote replay `excl-d5` queued)
For D = F_v = v·grad F (the only E1-admissible shape at degree 5, §1) the general solution of E1 is N = A v + (m/6) z + K with K = sum_{i<j} l_ij (F_j d_i − F_i d_j) (Koszul; l_ij linear forms), and E2 reduces to
    F_v A_v − A F_vv + F_v div K − K(F_v) = F Q        (m drops out; A_v = v(A), F_vv = v^T Hess(F) v),
an identity verified exactly on random data (`d5-reduction-identity-check`).
**Hessian relaxation:** replace v⊗v by a free S ∈ Sym²(C^6) (21 unknowns) and v⊗l by a free T ∈ C^6⊗C^{15×6} (540 unknowns), Q ∈ R_4 (126): a homogeneous LINEAR system (2078 equations) solved exactly over Q: rank 435, nullity 252, and its projection W to Sym² has dimension 8. Every certificate of degree 5 gives v⊗v ∈ W with v ≠ 0. The annihilator of W is spanned by 13 quadrics whose Groebner basis is the monomial ideal (v0², v1², v2², v3², v4², v5², v0v1, v0v3, v1v3, v0v4, v3v4, v1v4, v2v5) — Singular dim 0 and an independent sympy Groebner check agree — so W contains no nonzero rank-one tensor (W = span{e_a e_b : a ∈ core, b ∈ {2,5}} ∪ …, only mixed core×binary pairs). **Hence no certificate with deg D = 5 exists over C.**
Positive control: for the Jacobian target Y0·∂F/∂X0 (served by every D) the same code gives W = Sym² (dim 21) and the cone is all of C^6, as it must.
Together with §1: **deg D ≥ 6 for every certificate of this instance** (any parity, any coefficient field ⊂ C). Receipts: `wave2-excl/receipts/excl-d5-local*.{jsonl,json,sing}`.

## 3. Degree 6
### 3.1 ALL G-eigenvector Jacobian-type denominators of degree 6 — EXCLUDED over C [EXACT] (local runs `excl-d6jac-local`, `excl-d6jac2-local`)
E1-admissible denominators of degree 6 are D = F_w + c with w ∈ (R_1)^6 and c ∈ ker_6 (33 reps). A Jacobian-type D = Σ_a w_a F_a = Σ w_ac z_c F_a is a G-eigenvector iff the character difference wt(z_c) − wt(z_a) is the same for all pairs (a,c) used; the pairs fall into 13 difference classes (variable classes {X0,X1}, {Y0,Y1}, {X2}, {Y2}): the 10-parameter class δ=0 (w_a ∈ span of its own class; contains all weighted Euler polars Σ λ_a z_a F_a), the 8-parameter Y-odd class {z_c F_a : a∈{X0,X1}, c∈{Y0,Y1}} ∪ {a∈{Y0,Y1}, c∈{X0,X1}}, the 2-parameter classes X2·∇_{X0,X1}F, Y2·∇_{X0,X1}F, X2·∇_{Y0,Y1}F, Y2·∇_{Y0,Y1}F, (X0,X1)·F_X2, (X0,X1)·F_Y2, (Y0,Y1)·F_X2, (Y0,Y1)·F_Y2, {X2 F_Y2, Y2 F_X2}.
For each class, E2 reduces to F_w w(A) + A(F_w div w − w(F_w)) + F_w div K − K(F_w) = F Q with K Koszul (quadratic coefficients restricted to the forced sector) and Q in its sector; the Hessian relaxation (S ∈ Sym²(params), T = params⊗k, Q) is solved exactly over Q and the rank-one locus of its Sym²-projection W is computed by Groebner bases (sympy; Singular `minAssGTZ` for the 10-parameter class):
* δ=0, 10 params: W has dim 26; the rank-one locus is exactly the Euler line w_00 = w_11 = w_33 = w_44 = w_22 = w_55, off-diagonal 0, i.e. D = λ Σ z_a F_a = 6λF, not coprime to F;
* 8-param Y-odd class: dim W = 3, 33 quadrics, only v = 0;
* the two 4-parameter halves separately: dim W = 1, 9 quadrics, only v = 0;
* all 2-parameter classes: W = 0.
**Hence no G-eigenvector Jacobian-type denominator of degree 6 admits a certificate.** (Non-eigenvector Jacobian members, 36 parameters, and the kernel part c are only sampled, §3.2–3.3.) Receipts: `receipts/excl-d6jac*-local*`.

### 3.2 Individual degree-6 denominators (full system, two primes) — pending (job excl-d6, queued on the saturated remote); the 33 core-kernel representatives are being excluded exactly by the sector-restricted sweep (§4, local chain).

## 4. Core-only kernel elements, degrees 6–12, sector-restricted full 6-variable system [RANK] — running (local chain `excl-symsweep-local-k6-8`, then `excl-symsweep-local-k9-12`)
Every exact rational kernel representative (verified by `w2lib.small_kernel_reps`: modular-lattice reduction + exact membership) is tested in the full 6-variable system restricted to the forced character sector (lossless for these G-eigenvectors): modular ranks at 10^9+7 plus the exact Koszul nullity ([RANK]) for d ≤ 8, exact rank of M for d ≥ 9. Results so far (186 tests planned for k = 6..8, ~1.5 min each on the shared box): all 'inconsistent over C [RANK]' with exact nullity equal to the modular nullity (e.g. ker8 rep #0 even: rows 4204, cols 996, rank 856/857, nullity 140 = 96 + 44). Per-D records: `wave2-excl.log.jsonl` (exp = sym-fixedD), raw: `wave2-excl/remote/excl-symsweep-local-*.jsonl`.
## 5. Sum principle / served sectors — pending (job excl-sum)
## 6. Generic E1-admissible members, degrees 7–10 — pending (job excl-generic)
