# Add the Kronecker-substitution identity checker to the shared Lean library

Maintenance PR; no solve. It adds one reusable module, its test module, and documentation. The protected target records, `HodgeComputer/Verify.lean`, the verifier scripts, workflows, and all accepted instance files are unchanged.

## What

- `formal/lean/HodgeComputer/Reflection/Kronecker.lean` (namespace `HodgeComputer.Kronecker`): explicit integer polynomials as balanced monomial trees (`PT`), expressions (`E`) with products, sums, differences, natural scalars, powers, and partial derivatives, and the soundness theorems `soundN`/`soundC` (six-variable Kronecker point, values computed as pairs of naturals with kernel-accelerated `Nat` shifts and multiplications) and `sound4` (sliced by the degree in `X 4`). `toC_isHomogeneous`, `toC_ne_zero_of_evalPT`, and the witness binding `W`/`denW`/`matchW_sound` for the verifier's generated `HodgeWitness` definitions.
- `formal/lean/HodgeComputer/Reflection/KroneckerTest.lean`: a small identity using every expression constructor, and the cleared certificate equations `e1`, `e2` of the accepted `m3-coupled` witness (scaled by 6 to integers; 78k monomial products in `e2`) proved through `soundN` in the `directional`/`divergenceNumerator` form of `normalized_identity`, with homogeneity and nonvanishing of `D`, and `#print axioms` for the soundness theorems and the applications.
- `formal/lean/HodgeComputer/Reflection/README.md`: mechanism, soundness statement, how a solve PR uses generated witness files, the budget rule (about 20–30k witness terms per 60-minute check; memory governed by the degree), and pitfalls.
- A paragraph in `formal/README.md` under reusable constructions.

The import root `formal/lean/HodgeComputer.lean` is unchanged: like the `Constructions` modules, the new modules are built by name, and `scripts/check-library.py` picks them up as changed shared modules.

## Why

The accepted template `MvPolynomial.funext; norm_num <;> ring` is superlinear in the number of monomial products and in the number of distinct monomials of the expansion. Measured on the verifier's configuration (`lean -j 2 -M 14336`), random-support identities of 20–30k products cost 200–330 s and 7.5–9.3 GB, and 185k products fail at the heartbeat limit with 12.5 GB; structured accepted identities reach about 100k. Candidate witnesses for open targets with `d ≥ 8` exceed this. The Kronecker route has kernel cost linear in the number of witness terms and independent of the product count: the same `e2` with 435k products checks in 58 s and 5.5 GB, and a `d = 10` witness of 18.5k terms with 3.1M products in 263 s and 13.5 GB. A complete proposal for `m3-coupled` generated on this route passed the exact verifier compile sequence locally in 184 s (template `Witness.lean` alone: 275 s). The soundness argument is the uniqueness of bounded signed digits in a mixed-radix expansion, with the homogeneity and ℓ¹ hypotheses decided by `decide +kernel`; the module introduces no axioms and is entirely in terms of Mathlib's `MvPolynomial`. The earlier special-purpose encodings in `TransitiveCubic` (`TransitiveCubicKronecker`, `TransitivePolynomialDigits`) remain unchanged.

## How verified

Pinned toolchain `leanprover/lean4:v4.33.0` with the committed manifest, from `formal/lean` of a fresh worktree on `origin/main` (83fbf5c9):

- `lake build --wfail HodgeComputer.Reflection.Kronecker HodgeComputer.Reflection.KroneckerTest`: completes with no warnings (the five linter findings of the study version were fixed; no linter option is disabled in the library module).
- `lake env lean -j 2 -M 14336 <file>` for both files, the verifier's compile configuration: success.
- `lake env leanchecker --verbose <module>` for both modules, as `scripts/check-library.py` does: success.
- `#print axioms` in the test module for `soundN`, `sound4`, `matchW_sound`, `KroneckerTest.small`, `KroneckerTest.e1`, `KroneckerTest.e2`: `[propext, Classical.choice, Quot.sound]` each.
- `python3 scripts/check-library.py --base origin/main` (see the final section below for the result of that run).
- `scripts/build-content.py`/`validate-records.py`/`check-content.py` were not run locally (`mistune` is not installed here); the change touches no website content, and Repository build runs them.

No accepted application imports the new module, so `check-library.py` rechecks no protected application; `HodgeComputer.Verify` is unchanged.

## Timings

Local workstation, `taskset` to two or four cores, Lean 4.33.0, mathlib cache warm.

| Step | Wall | User | Peak RSS |
| --- | --- | --- | --- |
| `lake build HodgeComputer.Reflection.Kronecker` | 21 s (28 s cold) | – | 3.6 GB |
| `lake build HodgeComputer.Reflection.KroneckerTest` | 97 s | 92 s | 8.6 GB |
| `lean -j 2 -M 14336 Kronecker.lean` | 32 s | 42 s | 3.6 GB |
| `lean -j 2 -M 14336 KroneckerTest.lean` | 97 s | 88 s | 8.5 GB |
| `leanchecker --verbose HodgeComputer.Reflection.Kronecker` | 19 s | 11 s | 2.9 GB |
| `leanchecker --verbose HodgeComputer.Reflection.KroneckerTest` | 97 s | 83 s | 8.2 GB |

The test module's memory is the full-span Kronecker numbers of the degree-10 equations (`m = 21`, `W = 61`, 31 MB per number, cached per declaration), not its 620 terms; `set_option Elab.async false`, as in the accepted witness files, keeps the peak at the largest declaration (10.9 GB with concurrent elaboration). Reference points from the capacity study on the verifier's configuration: template `e2` with 18.6k products 205 s / 7.5 GB, 28.5k products 306 s / 9.2 GB, 185k products fails; Kronecker `e1 + e2` with 435k products 58 s / 5.5 GB.

## Checklist

- Shared library change only; separate from any solve PR.
- No change to `instance.rec` files, `why.md`, `Verify.lean`, `scripts/proposals/`, workflows, website, leaderboard, or generated data.
- Attribution: developed in the Lean capacity study prepared with Claude Code assistance for `kejace`; mathematics is the classical Kronecker substitution, formalized against Mathlib's `MvPolynomial`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
