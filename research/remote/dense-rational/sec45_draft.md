### 4.5 Mixed τ-parity denominators reduce to τ-pure ones of double degree
Legality of `D` depends only on the reduced support `V(D)` (prescribed-support theorem), and the set of legal
denominators is τ-stable (`(D,N,Q) ↦ (τD, −τ'N, −τQ)` with `τ'` the parity-twisted pullback of a vector field serves the
same `A`). If `D = D_+ + D_-` is legal with `V(D) ≠ τV(D)`, then `V(D) ∪ τV(D)` is a τ-stable support of a legal
denominator and its squarefree equation `G` satisfies `τG = ±G`, `deg G ≤ 2 deg D`. Hence
"no τ-pure legal denominator of degree ≤ 2d" ⇒ "no legal denominator of degree ≤ d" (the τ-pure case of degree ≤ d is
included). For d = 7 the τ-pure sectors are decided by §2.4; the genuinely mixed case would need τ-pure degree ≤ 14,
which is far beyond E1-bookkeeping (kernel 702 at d = 14). The honest statement is therefore:
**no degree-7 denominator with τ-stable support** (rigorous) + numerical evidence in the mixed sector.

### 4.6 Why the Hessian relaxation is a Gröbner problem only at degree 7
Hilbert-series bookkeeping for `m` generic quadrics in `n` variables: emptiness is certified at the first Macaulay
degree `k` where `Σ_j (−1)^j C(m,j) C(n−1+k−2j, k−2j) < 0`. d = 7: `(n,m) = (60,483)`: negative at `k = 4`
(`C(63,4) − 483·C(61,2) + C(483,2) = −171,822`) → matrix `883,890 × 595,665`, feasible (this is what msolve did).
d = 8: `(156, 468)` / `(159, 468)` per sector: the alternating count stays positive through `k = 8` at least (`k = 4: +2.0·10^7`, `k = 5: +5.4·10^8`, `k = 6: +1.1·10^10`, `k = 8: +2.7·10^12`), i.e. generic emptiness is not certified below Macaulay degree 9 (`C(164,8) ≈ 10^13` columns already at degree 8). So at degree 8 the relaxation is not decidable by F4 and must be attacked numerically or by a smarter reduction (e.g. the *bilinear* structure in `(w, K)` of the full E2 system, dossier §5.4).
Mixed d = 7 sector `(120, 951)`: `k = 4`: `+2.6·10^6`; first negative at `k = 5` (`−1.4·10^6`, barely) with `C(124,5) = 2.3·10^8` columns — infeasible on this machine without a reduction.

### 4.7 Geometric meaning of the degree-7 result
A τ-pure degree-7 denominator is a polar divisor `Y_w = {z ∈ X : w(z) ∈ T_z X}` for a quadratic vector field `w`
(§4.1). The failed condition `w(w(F)) ∈ J` says: there is no quadratic field `w` on `P^5` such that the second
"tangency derivative" `w²(F)` vanishes on `X` to first order along `X` modulo the Jacobian ideal. Geometrically, the
(2,0)-budget (§1.3–1.4) requires `Y_w` to be non-normal along a surface `Σ` with `h^1(I_Σ(7)) ≥ 208` — a surface in `P^5`
of degree ≤ 42 (`Y_w ∩ (Y_w)^{sing}`) failing 7-normality by 208. The msolve result excludes this *without* needing
the Hodge-theoretic simplicity assumption: it is a pure E1/E2 statement.
