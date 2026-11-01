# Attack 1 — `four-plus-two-residual` (q = 1, A = X0^3 Y0^3, F = H(X0,X1,Y0,Y1) + B(X2,Y2))

**STATUS: NO WITNESS FOUND.** Every tested denominator is inconsistent (§2). Main results of this
session: (i) a *rigorous* degree bound from the E1 kernel: **no denominator of degree ≤ 4 exists, and at
degree 5 only Jacobian elements pass E1 (all tested, inconsistent)**; (ii) the E1-survivors at degree 6
are exactly the 20 (Y-even) + 13 (Y-odd) core-only sextics `D(X0,X1,Y0,Y1)` with `A·D ∈ J(H)`, forced by
Macaulay duality — nothing geometric singles out a candidate; (iii) a **reduction to a threefold**: for a
core-only `D`, `V(D) ∩ X` is fibred over `P^1 = P(X2,Y2)` with fibres `V(D) ∩ Y_H`, `Y_H : U^6 = −H` the
cyclic sextic threefold, and the certificate problem reduces to supporting the class `[aU] ∈ H^{2,1}(Y_H)_{ζ6^2}`
on the surface `V(D) ∩ Y_H`; the algebraic form of this reduction is proved in one direction (§4.2) and
validated numerically on two solved instances (§2.3); the reduced linear system is ~20× smaller, so the
whole core E1-kernel was swept through degree 8 (plus degrees 9–10 partially): all inconsistent; (iv) a
sector computation (§1) locating the target in the Shioda decomposition: it lies in the *same* μ6-sector
as the holomorphic 4-form (the "triple-cover" sector `(ζ6^2, ζ6^4)`), separated from it only by the
Y-sign `ι_1 : (Y0,Y1) ↦ (−Y0,−Y1)`; the τ-invariant sub-Hodge structure containing the target is
`M' = H^3(T_H)^{ι_1 = −1}_χ ⊗_K H^1(P)_{χ^{-1}} ⊕ conj`, of rank 288 with `(h^{3,1},h^{2,2},h^{1,3}) = (64,160,64)`
(K = Q(ζ3), T_H the CY threefold `v^3 = H` in P(1,1,1,1,2), P the 3-dimensional Prym of benediamond's
genus-4 curve `C_B : w^3 = B`). This is level 2, so the generalized Hodge conjecture *predicts* a divisor
support — but the only known mechanism (Fano companion + elliptic contraction) fails for a dimension
reason (§1.4), not for a character reason.

Work directory `/tmp/hodge/results/four-plus-two-residual/`: `log.jsonl` (one JSON per tested denominator;
`"engine": 2` = 6-variable engine after the 2026-11-01 sign fix; `"engine": "threefold-v1"` = reduced
model of §4.2), `h.py` (6-var helpers), `e1kernel.py` + `e1kernel*.json(l)` (6-var E1-kernel sweep),
`corekernel.py/.json` (4-var kernel by core degree), `threefold.py` (reduced solver + core kernel +
rational lift), `threefold_struct.py`, `threefold_cones.py`, `split_neighbour.py`, `kernel6.py`,
`batch1/3.py`, `prym.py` (point counts), `*.out` (raw logs).

Conventions. Variables `X0,X1,X2,Y0,Y1,Y2`. `F = H + B` with

    H = f(X0,X1) + g(Y0,Y1) + K,   f = X0^6 + 5/7 X0^5X1 + 10/7 X0X1^5 + X1^6,
    g = Y0^6 + 5/7 Y0^5Y1 + 15/7 Y0Y1^5 + Y1^6,   K = 5/7 X0^3X1Y0Y1 + 10/7 X0X1Y0^3Y1,
    B = X2^6 + 15/7 X2^4Y2^2 + 10/7 X2^2Y2^4 + Y2^6 = P2(X2^2, Y2^2).

`a := X0^3Y0^3` (the target, a core polynomial). `J(F) = J(H) + J(B)` and `R(F) = R(H) ⊗ R(B)` (Thom–
Sebastiani); Hilbert series `R(H) ~ (1+t+..+t^4)^4 = (1,4,10,20,35,52,68,80,85,80,…)`, `R(B) ~ (1,2,3,4,5,4,3,2,1)`.
Symmetry group of F (mod scalars): `(Z/2)^3 = ⟨ι_1, ι_{X2}, ι_{Y2}⟩`, `ι_1 = τ∘ι_{Y2}` (no diagonal
`μ_k` beyond this: the two coupling monomials of K have bidegrees (4,2) and (2,4), so `(ζX, ηY)` preserves H
only for `η = ±ζ`). Residue characters: target `(ι_1, ι_{X2}, ι_{Y2}) = (−,−,−)`, τ-invariant; the 4-form
`ω_X`: `(+,−,−)`, τ-odd.

---

## 1. The obstruction and the geometric mechanism a certificate would need

### 1.1 Shioda sectors (verified bookkeeping)

For `F = F_1(x) + F_2(y)` with `n_1, n_2` variables, `R(F_1)_i ⊗ R(F_2)_j ⊂ R(F)_{i+j}` sits in the sector
`(c_1, c_2)`, `c_k ≡ −(deg + n_k) (mod 6)`, of `H^{n_1−1}(Y_1)_{c_1} ⊗ H^{n_2−1}(Y_2)_{c_2}`, `Y_k : u^6 = F_k`
(with the piece `R(F_k)_i · u^{c_k − 1}` of the cover's Jacobian ring). Here:

| class | core piece | binary piece | sectors | cover numerators |
|---|---|---|---|---|
| target `a ⊗ 1` | `a ∈ R(H)_6` | `1 ∈ R(B)_0` | `(2, 4)` | `aU ∈ H^{2,1}(Y_H)`, `V^3 ↔ dt/w ∈ H^{1,0}(C_B)` |
| 4-form `1 ⊗ 1` | `1 ∈ R(H)_0` | `1 ∈ R(B)_0` | `(2, 4)` | `U ∈ H^{3,0}(Y_H)`, `V^3` |

(This reproduces benediamond's characters `ζ6^2`, `ζ6^4` in the discussion.) Descending by `u = U^2`,
`w = V^2`: `T_H : v^3 = H ⊂ P(1,1,1,1,2)` is a **smooth CY threefold** (`h^{2,1} = 103`, the weighted
Jacobian ring is `R(H) ⊕ R(H)·v`), its μ3-sector χ (the one containing `ω_{T_H}`) has Hodge numbers
`(h^{3,0},h^{2,1},h^{1,2},h^{0,3}) = (1, 68, 35, 0)` with `H^{2,1}_χ ≅ R(H)_6`, `H^{1,2}_χ ≅ R(H)_{12}·v`;
`C_B : w^3 = B` has genus 4 with `H^{1,0}_{χ^{-1}} = ⟨dt/w⟩` (1-dim), `H^{1,0}_χ` 3-dim; `H^1(C_B) = H^1(E) ⊕ H^1(P)`
with `E` the j = 0 elliptic quotient and `P` the Prym, `H^{1,0}(P) = ⟨dt/w, dt/w^2, t^2dt/w^2⟩`.

So **the target and the holomorphic 4-form lie in the same μ6-sector** `M = H^3(T_H)_χ ⊗_K H^1(C_B)_{χ^{-1}} ⊕ conj`,
and both have binary piece `dt/w ∈ H^1(P)` (both are `ι_{X2}`-, `ι_{Y2}`-anti-invariant). What separates
them is the Y-sign `ι_1` on the core: `ω_{T_H}` is `ι_1`-invariant, `[a]` is `ι_1`-anti-invariant
(`a` is Y-odd, the core volume form is Y-even). Hence the τ-invariant sub-Hodge structure containing the
target is

    M' = H^3(T_H)^{ι_1 = −1}_χ ⊗_K H^1(P)_{χ^{-1}} ⊕ conj.

**Hodge numbers (verified via parity-graded Hilbert functions).** `H^3(T_H)^−_χ`: `(0, 32, 16, 0)`
(`R(H)_6^{odd} = 32`, `R(H)_{12}^{odd} = 16`; the parity-graded Hilbert function of `R(H)` equals that of the
split `R(f) ⊗ R(g)` because J(H) is generated by a regular sequence of parity-homogeneous quintics).
`H^1(P)_{χ^{-1}}`: `(1,0)`: 1, `(0,1)`: 2. Therefore `M'_χ` has `(h^{3,1}, h^{2,2}, h^{1,3}) = (32, 80, 32)` and
`M'` (over Q) `= (64, 160, 64)`, rank 288, **level 2 (no (4,0))**. In Jacobian-ring terms the 64
τ-invariant (3,1)-directions of this sector are `R(H)_6^{odd} ⊗ 1` (32, containing `a`) and
`R(H)_4^{odd} ⊗ ⟨X2^2, Y2^2⟩` (32). None of them is odd in `X2` or `Y2`, so the book's reflection theorem
(Ch. 3 §3.3) serves nothing in this sector, as why.md says; the two-reflection theorem of Ch. 2 (needs
`A` odd in `x` and `y`) likewise.

### 1.2 What GHC predicts

`M'` has no (4,0) part, so a divisor support is *predicted* by the generalized Hodge conjecture (coniveau
1 for the sub-HS generated by `[A]`, which is probably all of `M'` — `H^3(T_H)^−` is a non-rigid level-1
structure (48-dim over K, "Weil-type" signature (32,16)) and `P` is a non-rigid Prym of signature (1,2)).
Conversely, by the book's §3.1 support theorem a certificate forces exactly such a support, so the
instance is a genuine (conjecturally true, currently unproved) GHC case — same logical status as
`polarized-k3-remainder` (Λ²T there), but with an abelian-type piece `H^3(T_H)^− ⊗_K H^1(P)` rather than
an exterior square of a K3 lattice.

### 1.3 Why the split neighbour is solvable and this one is open

Setting `K = 0`: `F_0 = f + g + B` is *exactly* in the Chapter 2a family `P0(X0,X1) + P1(Y0,Y1) + P2(X2^2,Y2^2)`
(all Y-odd targets served, unbounded denominator degree). Then `T_{f+g}` is a join of two genus-4 curves
and `H^3(T_{f+g})^−_χ` splits into Shioda pieces `H^1(C_f)_{c_f} ⊗ H^1(C_g)_{c_g} ⊗ (u)` with `c_f + c_g ≡ 2`,
`c_g` odd: `(1,1)` (16-dim, CM⊗CM, all of type (2,1): contains `[a] = X0^3 ⊗ Y0^3`), `(3,5)`, `(5,3)` (16 each,
hyperelliptic ⊗ CM). The coupling K merges these three pieces into one 48-dimensional varying structure
(the only surviving core symmetry is `ι_1`). The 2a certificates for `X0^3Y0^3` use the two-block contact
packet `D = r_f^2 + r_g^2`, `A = J_f J_g` (`f = C_f^3 + r_f^2` over a 40-sheet cover) and separate the
target over sheets, so `deg D_0` is a multiple of 6 times a sheet count — not a low-degree object, and not
transferable (E1 for F forces `D` into the small kernels of §2.1). "Deformation from the solved neighbour"
(experiment 4) therefore has no low-degree seed; the §2.6 control shows the split neighbour itself has
no core-only certificate of degree ≤ 8 either.

### 1.4 Why the 3A companion mechanism fails here (dimension count, my analysis)

The 3+3 proof realizes `H^1(C_f) ⊗ H^1(C_g)` classes as `Γ_*(α ⊗ β)` with `α, β` supported on divisors
`D_f ⊂ W_f`, `D_g ⊂ W_g` of the two Fano companions `W_f : f + s^3 + t^3 = 0 ⊂ P(1,1,1,2,2)` (support from
rational connectedness) through the contraction correspondence `Γ ⊂ W_f × W_g × X` of dimension 5; the
image of `Γ ∩ (D_f × D_g × X)` has dimension `2 + 2 − 1 = 3`: a divisor of X. For 4+2 the companions are
`W_H : H + s^3 + t^3 ⊂ P(1,1,1,1,2,2)` (Fano fourfold, −K = O(2): supported) and `W_B : B + s^3 + t^3 ⊂ P(1,1,2,2)`
(K3, Picard rank 16; its transcendental lattice is exactly `T(W_B) = H^1(P) ⊗_K H^1(E_0)`, rank 6, and
the target's binary piece `dt/w` projects nontrivially to it, so β is *never* supported on a curve). With
only `D_H` available, `Γ ∩ (D_H × W_B × X)` has dimension `3 + 2 − 1 = 4` and its image is all of X
(`D_H` dominates `P^3` because `W_H → P^3` has 1-dimensional fibres). So the obstruction is not a character
mismatch (the characters match) but the absence of a *second* supported factor — exactly benediamond's
"missing support contribution", made quantitative.

### 1.5 The fibration / threefold reduction (new; §4.2 for the algebra)

Projecting from the plane `{X2 = Y2 = 0}` (base locus the sextic surface `S_H = {H = 0} ⊂ P^3`) gives
`X̃ = Bl_{S_H} X → P^1 = P(X2,Y2)` with fibre over `[t:1]` the sextic threefold `H(x) + B(t,1)Y2^6 = 0 ≅ Y_H`
(twisted by `B(t,1)^{1/6}`; degenerate fibres over the six roots of B are cones over `S_H`). The Shioda
sector `M` is `H^1(P^1, R^3π_*Q)` with the local system splitting by μ6-characters; the class
`[A] = [aU] ⊗ [dt/w]` is "`H^1` of the base curve `C_B` with values in the sector-χ fibre class `[aU]`".
A **core-only** denominator `D(X0,X1,Y0,Y1)` cuts the horizontal divisor `Z_D = ∪_t (V(D) ∩ Y_t)`; `Z_D`
supports `[A]` iff the fibre class `[aU] ∈ H^{2,1}(Y_H)_{ζ6^2}` is supported on the surface
`S_D := V(D) ∩ Y_H` (equivariantly for μ6, automatic). `S_D` is the μ6-cover of the surface `V(D) ⊂ P^3`
branched along the curve `V(D, H)`; for `[aU]` (a (2,1)-class) to come from `H^1(S̃_D)(−1)` one needs
irregularity `q(S̃_D) > 0` in the `ζ6^2`-eigenspace — impossible for a smooth (or rational-singular,
normal) `V(D)` by Lefschetz/vanishing, so `V(D)` must be singular in a way that creates irregularity
(cones with vertex on `S_H`, unions of planes through curves of `S_H`, non-normal surfaces, …). This is the
3-variable-smaller analogue of the book's (3,1)-support constraint, and it is what the cone families of
§2.5 probe. The reduced certificate problem (target `aU`, hypersurface `U^6 + H` in 5 variables, core-only
`D`) is a 5-variable linear system; a solution of its strict form gives a strict 6-variable certificate
(§4.2, proved), and the known witnesses of two solved instances with core-only D are consistent in the
reduced model (§2.3), so the reduced model is a faithful, ~20× cheaper proxy.

---

## 2. Experiments

All linear algebra mod `p = 2^31 − 1` (dense flint rref) unless marked exact. "rank/nullity" refer to the
coefficient matrix (unknowns N, Q) of the fixed-D system; "consistent" = augmented rank equals rank.
Timings are wall-clock on a heavily shared 40-core box (load ≈ 41, jobs niced), so they overstate cost by
3–5×.

### 2.1 E1 kernel sweeps (6 variables: `{D : A·D ∈ J(F)}` mod `J(F)_d`; 4 variables: `{D : a·D ∈ J(H)}` mod `J(H)_k`)

6-variable, Y-even D (file `e1kernel.jsonl`):

| d | #monomials | dim E1-kernel | dim J_d | kernel mod J | rows × cols | s |
|---|---|---|---|---|---|---|
| 1 | 3 | 0 | 0 | **0** | 172×66 | 0.1 |
| 2 | 12 | 0 | 0 | **0** | 402×180 | 0.1 |
| 3 | 28 | 0 | 0 | **0** | 792×406 | 2 |
| 4 | 66 | 0 | 0 | **0** | 1357×822 | 7 |
| 5 | 126 | 3 | 3 | **0** | 2111×1512 | 58 |
| 6 | 236 | 38 | 18 | **20** | 3054×2612 | 96 |
| 7 | 396 | 126 | 63 | **63** | 4279×4257 | 452 |

4-variable core kernel `ker(a· : R(H)_k → R(H)_{k+6})` by Y-parity (file `corekernel.json`, two primes agree):

| k | even: dim mod J(H) | odd: dim mod J(H) | note |
|---|---|---|---|
| 0–4 | 0 | 0 | |
| 5 | 0 | 0 | kernel = J(H)_5 only |
| 6 | **20** | **13** | 20+13 = 33 = 68 − 35: `a·: R(H)_6 → R(H)_12` is surjective (its dual `a·: R(H)_4 → R(H)_10` is injective) |
| 7 | 30 | 30 | |
| 8 | 41 | 34 | |
| 9 | 38 (even) | (see sweep log) | |

Consequences (rigorous, mod-p ranks reproduced with two primes at k ≤ 8): **no denominator of degree ≤ 4
exists for this instance**; at degree 5 the only E1-survivors are `D ∈ J(F)_5 = span(F_i)`; at degree 6
the survivors are `J(F)_6 ⊕ (core kernel_6 ⊗ 1)`; at degree 7 they are `J(F)_7 ⊕ (core kernel_7 ⊗ 1) ⊕
(core kernel_6^{even} ⊗ X2) ⊕ (core kernel_6^{odd} ⊗ Y2)` (63 = 30 + 20 + 13), i.e. the only non-core-only
degree-7 candidates are reducible `c(x)·X2`, `c(x)·Y2` (support `V(c) ∪ {fibre}`). Rational reconstruction
of the rref kernel bases fails even with 6 primes (≈186 bits): the bases have large heights; all tests
below on kernel elements are therefore mod p only (a hit would be lifted with more primes).

### 2.2 Full 6-variable fixed-D tests (engine 2)

| D | d | family | E1 | consistent | rank | nullity | rows×cols | s |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | sanity | no | — | 18 | 0 | | 0 |
| X0²+2Y0², X0²−2Y0², X0²+Y0², X0X1+Y0Y1 | 2 | structural | no | — | 168 | 0 | | 1 |
| X0X1Y0Y1 | 4 | structural | no | — | 747 | 9 | | 11 |
| K (coupling) | 6 | structural | no | — | 2214 | 162 | | 79 |
| f(X0,X1), g(Y0,Y1) | 6 | structural | no | — | 2214 | 162 | | 150–195 |
| B (binary block) | 6 | structural (∈ J) | yes | **no** | 2299 | 303 | 15088×2602 | 221 |
| H (core block) | 6 | structural (∈ J) | yes | **no** | 2299 | 303 | 16370×2602 | 260 |
| ∂F/∂X0, ∂X1, ∂Y0, ∂Y1 | 5 | partials (∈ J) | yes | **no** (exact solve) | 1316 | 130 | | ~100 |
| ∂F/∂X2, ∂F/∂Y2 | 5 | partials (∈ J) | yes | **no** | 1302 | 144 | | ~100 |
| block radial derivatives, generic J_5 combos | 5 | radial / J5-combo | yes | see log (batch3) | | | | |
| random combos #0,#1,(#2) of the 20 d=6 kernel elements | 6 | e1kernel-d6-random | yes | **no** | 2353 | 249 | 16797×2602 | 480–610 |
| the 20 d=6 kernel rref basis elements | 6 | e1kernel-d6-basis | yes | see log / §2.3 (reduced model: all no) | | | | |
| Y0³+Y1³+l²Y0, l = X0 − tX1, f(t,1)=0 | 3 | cone3-vertex-on-S_H | — | **no** | 351 | 27 | 2681×378 | 4 |
| ~110 cubics in all 3-variable subsets (Fermat, cone, random) | 3 | cone3-cubic | no (all) | — | 378 | 0 | | 1–5 |

(The cubic batch is redundant with the E1 sweep — kernel 0 at d ≤ 4 — and was stopped.)

### 2.3 Reduced (threefold) model: validation and core-kernel sweeps

Positive controls (solved instances whose witness D is core-only in a 3-variable core, reduced target `a·U^2`):

| instance | core | witness D | reduced model | negative controls |
|---|---|---|---|---|
| family-hesse-web | (X0,X1,Y0) | X0⁶−2X0⁴X1²+3X0⁴Y0²−4X0²X1²Y0²+3X0²Y0⁴+4X1⁶−2X1²Y0⁴+Y0⁶ | **consistent** (rank 137, nullity 38) | X0⁴, X0⁶+X1⁶+Y0⁶: inconsistent |
| f3-block-x0-quartic-strict | (X0,X1,Y0) | X0⁴ − X1²Y0² | **consistent** (rank 69, nullity 15) | same: inconsistent |

Our instance (core `H`, reduced target `aU`, `F' = H + U^6`); sanity: `D = 1, H_{X0}, X0^6, H` all inconsistent.

RESULTS_TABLE_SWEEP

### 2.4 Reduced model: structural families

RESULTS_TABLE_STRUCT

### 2.5 Reduced model: cones with vertex on `S_H`, hyperplane unions, `a·L^6 ∈ J(H)` tests

RESULTS_TABLE_CONES

### 2.6 Control: the split neighbour `H_0 = f + g` (certificate exists by the 2a theorem)

RESULTS_TABLE_SPLIT

### 2.7 Prym splitting of `C_B` (point counts)

RESULTS_PRYM

---

## 3. Witness

None. `witness-candidate.json` was not created.

---

## 4. Mathematical ideas developed

### 4.1 Where the target lives (summary of §1)

`[A] = [aU] ⊗ [dt/w] ∈ M' = H^3(T_H)^{ι_1=−1}_χ ⊗_K H^1(P)_{χ^{-1}} ⊕ conj`, rank 288, types (64,160,64).
Both factors are non-rigid, abelian-type Hodge structures with `Q(ζ3)`-action: the core factor is the
Y-anti-invariant part of the intermediate-Jacobian-type structure of the CY threefold `T_H` (an abelian
48-fold of Weil type, signature (32,16)), the binary factor is the Prym `P` (signature (1,2)). The
target's sub-HS is of "`H^1(A_1) ⊗_K H^1(A_2)`" type for abelian varieties `A_1, A_2` with K-action; GHC
predicts coniveau 1 in `H^4(X)`, i.e. a threefold `Z ⊂ X` whose resolution has `H^2 ⊇` (this piece).

### 4.2 The threefold reduction (proved direction)

Let `D` be core-only of degree d. Suppose a vector field `ñ = (ñ_0, ñ_1, ñ_3, ñ_4)` on the core
(`deg ñ_i = d + 1`) satisfies

    (E1')  ñ(H) = a·D,          (E2')  ñ(D) = D · div ñ   (ñ is "divergence-logarithmic" along D).

Then with `m := 0` … more precisely put `N := ñ + m·E` where `E = Σ_{all six} z_i ∂_i` and `m` is any core
polynomial of degree `d` (`N_{X2} = mX2`, `N_{Y2} = mY2`): `N(F) = aD + 6mF`, `N(D) = ñ(D) + d·mD`,
`div N = div ñ + (d+6)m`, hence `D·div N − N(D) − D·(6m) = 0`: **a strict certificate (Q = 0) for X**.
Conversely, expanding the 5-variable identity for `(U^6 + H, aU, D)` in powers of `U` (the target is in
the `U ↦ ζU` character-1 sector, so `N'_x ∈ U·C[x] ⊕ U^7·C[x] ⊕ …`, `N'_U ∈ U^2·C[x] ⊕ …`, `Q' ∈ U·C[x] ⊕ …`)
shows that for `d ≤ 6` the only extra freedom beyond (E1'),(E2') is a constant vector `n_7 ∈ C^4` (term
`U^7 n_7`) with `q_0 = −n_7(D)`; for larger d analogous `U^{6k+1}`-terms appear. The reduced solver
(`threefold.py: reduced_test`) solves the full 5-variable system in the character sector (unknowns
`N'_core` with `U`-degree ≡ 1, `N'_U` ≡ 2, `Q'` ≡ 1 mod 6; `deg N' = d + 2`, `deg Q' = 2d − 5`), so it is
at least as permissive as (E1')+(E2'); the positive controls of §2.3 show it detects actual witnesses.
Geometrically (§1.5) it is the statement "`[aU]` is supported on `V(D) ∩ Y_H`", which by the book's
prescribed-support theorem (applied to the threefold `Y_H`) is equivalent to the reduced certificate.

What this buys: the fourfold problem for core-only D (the *only* E1-survivors at d = 6, and all survivors
at every degree up to the reducible `c·X2, c·Y2` products) is equivalent to a **coniveau question on the CY
threefold** `T_H`: is the (2,1)-class `[a] ∈ H^{2,1}(T_H)^{ι_1=−1}_χ` supported on a surface of `T_H`? Since
`H^3(T_H)^−` is of level 1, GHC predicts yes (its intermediate Jacobian should be an abelian variety
dominated by Jacobians/Pryms of curves on surfaces in `T_H`). The natural surfaces are the fixed loci of
`ι_1`: the two genus-4 curves `C_f = {X0=X1=0}`, `C_g = {Y0=Y1=0}` on `T_H` (both of the form `v^3 = binary
sextic`), and surfaces swept out by curves meeting them. In the split limit the support is built from
contact presentations of `f` and `g` over a 40-sheet cover (Chapter 2a), i.e. from the 3-torsion geometry
of the genus-2 Jacobians `w^2 = f`, `w^2 = g` — this is what has to be *continued* in the coupling
parameters `(5/7, 10/7)`.

### 4.3 Reducible denominators and the joining mechanism

Degree-7 survivors `c(x)·X2` have support `Z_c ∪ Y_0` with `Y_0 = {X2 = 0} ∩ X ≅ Y_H` a fibre (which alone
carries nothing: `H^2(Y_H) = Q(−1)`). The union can still support more than its parts via Mayer–Vietoris:
`H_4(Z_c ∪ Y_0) ⊃ ker(H_3(Z_c ∩ Y_0) → H_3(Z_c) ⊕ H_3(Y_0))`, i.e. classes from `H^1` of the surface
`S_c = V(c) ∩ Y_0` glued across the two components (the book's "joining divisor" in 3A.4–3A.5 is of this
kind). These products were not tested in this session (each 6-variable d=7 run costs ~20 min on the shared
machine); see next steps.

### 4.4 Negative structural facts (publishable)

1. `deg D ≤ 4`: impossible (E1 kernel zero, two primes). `deg D = 5`: only `D ∈ J(F)_5`; the six partials
   and the block radial derivatives are inconsistent (exact rational solve), as are generic combinations.
2. `deg D = 6`: all E1-survivors are `J(F)_6 ⊕ (core kernel ⊗ 1)`; `J(F)_6`-elements `B`, `H`, generic
   `Σ ℓ_i H_i` inconsistent; random elements and every rref basis element of the 20-dim even core kernel and
   the 13-dim odd one inconsistent (reduced model; two random combos also in the full 6-variable system).
3. `deg D = 7, 8` core-only: all rref basis elements and random combinations inconsistent (reduced model);
   `deg 9, 10`: random combinations and the first basis elements inconsistent.
4. Cones over plane curves (degrees 2,3,4,6; Fermat, `g`, random; vertex at the F_p-rational points of
   `S_H` on the coordinate lines and at a general point of `S_H`) and six-plane unions `H + λL^6`, `L_1^6 ± L_2^6`
   are inconsistent (most fail E1 already: `a·L^6 ∉ J(H)` for every tested `L`, including the special ones).
5. The split neighbour `f + g` has no core-only certificate of degree ≤ 8 either (§2.6), consistent with
   the unbounded-degree nature of the 2a construction.

---

## 5. Ranked next steps (with estimated compute)

1. **Reducible d=7 products `c·X2`, `c·Y2`** (33 candidates, 6-variable engine, compressed solver): random
   combinations first (3 runs, ~20 min each on the loaded machine, ~5 min when idle; 1–2 GB each).
   This is the only degree-≤7 family not covered by the reduced model.
2. **Finish the reduced-model sweep at degrees 9–10 and extend to 11–12** (core kernels ≈ 40–60 per
   parity; 4–10 min per test at d = 9–10, ~20–40 min at d = 12): ~1 day of a single core; memory < 4 GB.
   If anything is consistent, lift D with 8–12 primes and run the full 6-variable exact solve.
3. **Contact presentations mod p for `f` and `g`** (`f = C^3 + r^2`: 7 equations, 7 unknowns over F_p; some
   of the 40 sheets are F_p-rational for suitable p) and test 2a-type denominators `r_f^2 + r_g^2`,
   `r_f^2 + g`, `f + r_g^2`, their products, and small perturbations by `K`, on **both** `H_0` (must eventually
   work by theorem, possibly only after sheet-products/trace) and `H` (reduced model, seconds each). This
   is the first concrete probe of "continuation in the coupling": compare the ranks of the two systems
   along `H_s = f + g + sK`, `s ∈ F_p`, to see whether the consistency locus is `s = 0` only.
4. **Point counts for the core** (expensive but decisive for the GHC framing): `#T_H(F_{p^k})` is a
   3-fold count; instead count the genus-4 curves `C_f, C_g` and the surface `S_H` to test whether
   `H^3(T_H)^−` has extra endomorphisms/CM (which would make an explicit support more plausible); the
   L-polynomial of the Prym `P` (running: `prym.py`) answers the binary side.
5. **Mixed denominators of degree 8** `c_4(x)·b_4(X2,Y2)` etc. (core kernel at k=4 is zero, so the first
   mixed survivors are `kernel_6 ⊗ R(B)_2` at d = 8: 20·3 + 13·3 candidates; the (3,1)-sector shows the
   `R(H)_4^{odd} ⊗ ⟨X2², Y2²⟩` directions are part of the target's Hodge structure, so denominators mixing
   `X2², Y2²` with core quartics are natural): 6-variable d=8 runs, ~GB and ~30–60 min each when idle.
6. Theory: identify the intermediate Jacobian piece `J(T_H)^−` with a Prym of a curve on `T_H` (candidates:
   the double cover of `T_H/ι_1` structure, the curves `C_f`, `C_g`, and conics/lines meeting them); a
   correspondence of that kind would supply the missing core-side support and, combined with the
   fibration of §1.5, a core-only D. This is the analogue of the 3A "elliptic contraction" with the
   curve-side support coming from the core instead of the binary block.
