# Kronecker-substitution checker

`HodgeComputer.Reflection.Kronecker` certifies large exact identities between explicit polynomials in `MvPolynomial (Fin 6) ℂ` by kernel computation on big integers. It replaces the `MvPolynomial.funext; norm_num <;> ring` template, whose cost grows superlinearly in the number of monomial products and which fails near 10⁵ products within the verifier's 14 GiB budget, by a check whose kernel cost is linear in the number of witness terms. It was developed in a capacity study on the accepted `m3-coupled` and `contact46006` witnesses and synthetic certificate identities; the numbers below come from that study.

## Mechanism

An explicit polynomial with integer coefficients is a balanced binary tree of monomials, `PT`, whose leaves are sorted by their Kronecker position `pos m` (base `m`, with `X 5` dehomogenized). `toC t` is its denotation in `Poly`, written in the same `C c * (X 0 ^ a * … * X 5 ^ f)` shape as certificate modules. Expressions `E` combine trees by `mul`, `add`, `sub`, `smulN` (natural scalar), `pow`, and `pd i t` (the partial derivative of a tree, computed on the leaves by `pderivPT`); `denC e` is their denotation.

For an identity `L = R` the kernel decides three Boolean or arithmetic facts:

- `homE L`, `homE R`, and `degCompat (degE L) (degE R)`: both sides are homogeneous of one degree `n`, and `n < m`;
- `l1E L + l1E R < 2 ^ W`: the ℓ¹ norms bound every coefficient of `L - R` below `2 ^ W`;
- `psiEq m W L R`: the Kronecker values `psiEN m W L` and `psiEN m W R` agree, where the value of a polynomial is its evaluation at `X i ↦ 2 ^ (W * m ^ i)` for `i < 5` and `X 5 ↦ 1`.

The soundness theorem is

```lean
theorem soundN {m W : ℕ} (L R : E) (hL : homE L = true) (hR : homE R = true)
    (hdeg : degCompat (degE L) (degE R) = true)
    (hm : degLR (degE L) (degE R) < m) (hb : l1E L + l1E R < 2 ^ W) (hpsi : psiEq m W L R = true) :
    denC L = denC R
```

Its proof is `eq_zero_of_eval_kpt`: a homogeneous polynomial of degree below `m` has pairwise distinct positions on its support (`posF_injective`, mixed-radix uniqueness), so the Kronecker value of `L - R` is a base-`2 ^ W` expansion with signed digits of absolute value below `2 ^ W`, which vanishes only if every digit vanishes. The identity is proved over `ℤ` (`sound`) and transported to `ℂ` by `MvPolynomial.map` (`soundC`, `soundN`). Values are computed as pairs of naturals (`NN`, positive and negative parts) so that every big-number step is a single kernel-accelerated `Nat` shift, multiplication, or addition; `psiRelN` evaluates relative to the minimal position of each subtree, so intermediate numbers stay proportional to the span of their subtree.

`sound4` is the sliced variant: `toC_eq_sum_sliceK` writes a tree as `Σ_k X 4 ^ k * toC (sliceK k t)`, each slice identity is between `X 4`-free expressions (`x4freeE`) and is checked at a Kronecker point on four variables (`m ^ 4` positions instead of `m ^ 5`), and the slices are reassembled with `linear_combination` using `pderiv_four_toC_sliceK` and the `pderiv_X4_pow_*` lemmas.

Auxiliary facts for certificate modules: `toC_isHomogeneous` (from `homB`), `toC_ne_zero_of_evalPT` (nonvanishing from one integer evaluation), and the witness binding below.

The `TransitiveCubic` library contains an earlier, special-purpose encoding (`TransitiveCubicKronecker`, `TransitivePolynomialDigits`) for the weighted covariant identities of the cubic family. The present module is the general six-variable, signed-coefficient, expression-level form with derivatives, scaling, and slicing.

## Use in a solve PR

Only the identity checks change; the protected `Verify` statement, the target, and the verifier are untouched. A generated `Witness.lean` helper (prepared by a script such as the capacity study's `kron_submit.py`, in the style of the existing `leangen` witness generators) contains:

1. The trees `tF`, `tA`, `tD`, `tN0`, …, `tN5`, `tQ` as `noncomputable def`s, leaves sorted by position, balanced. The witness `(D, N, Q)` is scaled to integers (`D`, `N` by `s`, `Q` by `s ^ 2`; the identities are homogeneous of degree two in the witness). If `F` or `A` has rational coefficients, write `F = C (1 / s') * toC tF'` with an integral `tF'` and check the identities multiplied by `s' ^ q`.
2. The parameters `m = (identity degree) + 1` and `W = bits(ℓ¹ bound) + 2`, where the bound is the one `l1E` computes (ℓ¹ norms multiply over products, add over sums).
3. `e1_core : denC L1 = denC R1 := soundN (m := m) (W := W) L1 R1 (by decide +kernel) × 6` and the same for `e2`, then `simp only [denC_*]` and unfolding to reach the `directional` and `divergenceNumerator` forms expected by `HodgeComputer.normalized_identity`. `KroneckerTest.lean` shows the exact tactic script on the accepted `m3-coupled` witness.
4. Homogeneity by `toC_isHomogeneous`, `D ≠ 0` by `toC_ne_zero_of_evalPT`, and legality by the accepted template's univariate specialization applied to the unfolded trees (`simp only [D, tD, toC, …]; push_cast` before `polynomial`).

To bind the verifier's protected `HodgeWitness` definitions, write the witness strings in `witness.json` as balanced parenthesized sums in the same position order as the leaves, with the coefficient conventions `X0**6`, `-X0**6`, `5*X0**6`. The verifier's translation (`scripts/proposals/check.py`) then is syntactically `denW w` for the mirrored term `w : W`, so `HodgeWitness.D = denW wD` holds by `rfl`, and `matchW_sound wD tD (by decide +kernel)` identifies it with `toC tD`. This replaces a `funext; norm_num; ring` binding whose cost is superlinear in the number of terms.

Build the module from `formal/lean` with

```sh
lake build HodgeComputer.Reflection.Kronecker HodgeComputer.Reflection.KroneckerTest
```

A proposal imports it as `import HodgeComputer.Reflection.Kronecker`; the library is part of the protected dependency cache.

## Budget

Measured on the verifier's configuration (`lean -j 2 -M 14336`, replay by `leanchecker`, 60 minutes for every module of a proposal): for degrees `d ≤ 8` the kernel cost is about 25–30 ms per witness term per identity, independent of the number of monomial products. A synthetic `e2` with 435k products checked in 58 s and 5.5 GB; a `d = 10` witness of 18.5k terms with 3.1M products in 263 s and 13.5 GB unsliced. `KroneckerTest.lean`, which checks the accepted `m3-coupled` equations (`d = 10`, 78k products in `e2`, 620 witness terms) together with homogeneity and nonvanishing, compiles in about 95 s and replays with `leanchecker` in about 95 s. Per 1000 witness terms, budget about 12 s for the verifier's own elaboration of `HodgeWitness.lean`, 17 s for the `rfl`/`matchW` binding, and 30 s per identity check, each kernel part replayed once more by `leanchecker`. A 60-minute budget therefore accommodates roughly 20–30k witness terms in total for `D`, the six numerators, and `Q`, regardless of the raw product count. Keep the witness itself small; the product count is irrelevant on this route.

Memory is governed by the degree, not by the term count. The kernel caches every intermediate value of a declaration, and each tree occurrence in an identity contributes about `depth × 2 × m⁵ × W` bits (full-span numbers at every level of its balanced tree, with their shifted copies) as soon as its leaves span the position range. The `m3-coupled` equation `e2` (`m = 21`, `W = 61`, thirteen tree occurrences) needs 8.6 GB at `lean -j 2`; degrees `d ≤ 8` (`m ≤ 17`) need about a quarter of that at any term count. For `d ≥ 10` with large witnesses, or when the peak approaches the 14 GiB limit, use `sound4`, which divides the numbers by `m` at the cost of one declaration and eight side goals per slice. Set `Elab.async false` in witness modules, as the accepted witness files do, so that concurrently elaborated identities do not add their peaks (10.9 GB instead of 8.6 GB for the test module).

## Pitfalls

- Leaves must be sorted by `pos m` and the tree balanced; unsorted leaves make every node span the whole range (one measured case: 275 s and 32 GB instead of 10 s and 0.85 GB).
- Large data definitions must be `noncomputable`; the code generator fails on 10k-leaf trees (`maximum recursion depth`) and wastes minutes before that.
- `m` must exceed the degree of the identity, not of the witness: `e1` has degree `d + 6q`, `e2` has degree `2d + 6q - 6`. `W` must exceed the ℓ¹ bound of `L - R` as computed by `l1E`; `decide +kernel` fails otherwise, it never accepts a false identity.
- `Nat.pow` is kernel-accelerated only for exponents below `2 ^ 24`; the module uses shifts (`<<<`) for Kronecker weights and `2 ^ W` only for the bound.
- Generated witness modules with thousands of terms should disable linters (`set_option linter.all false`), whose cost on huge terms is comparable to the proof itself; keep the library module itself lint-clean.
- Integer coefficients appear as `Int.cast` after unfolding `toC`; `push_cast` before `polynomial` or `norm_num` in the legality step.
- A `W` constructor name coincides with the usual bound variable `W : ℕ`; the binder shadows the type inside a theorem, which is harmless, but write `(W := 61)` for the named argument.
- The zero polynomial is a single zero leaf whose monomial carries the intended degree; `firstDeg` ignores zero leaves and `degE` returns `none` for an expression that denotes zero, which is compatible with every degree.
- A successful numerical equality without the homogeneity and ℓ¹ hypotheses proves nothing; `soundN` requires all six side goals.
