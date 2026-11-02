STATUS: in progress (started review)

# Review of maint/kronecker-checker (575e666b)

## (1) Soundness of the Lean module — read in full (lines 1-1287)

- No `native_decide`, `implemented_by`, `unsafe`, `partial`, `opaque`, `extern`, `axiom`, `sorry`, `macro`, `elab`; only `set_option autoImplicit false` (Kronecker.lean:36). Test module: `autoImplicit false`, `Elab.async false`.
- Imports: Kronecker.lean imports only `HodgeComputer.Verify`; the import root `HodgeComputer.lean` is untouched (no `Reflection` reference).
- `soundN` (945): hypotheses are exactly `homE L/R`, `degCompat`, `degLR < m`, `l1E L + l1E R < 2^W`, `psiEq m W L R`; conclusion `denC L = denC R`. Proof chain `psiE_eq_of_psiEq -> soundC -> sound -> eq_zero_of_eval_kpt` checked by hand: minimal-position digit argument with `posF_injective` (needs every exponent < m, supplied by homogeneity of degree n < m) is correct; the l1 bound via `absE`/`eval_one_absE`/`coeff_le_eval_one` is correct (products: `abs_coeff_mul_le`; `pd` via `pderivPT` on leaves).
- `sound4` (1245): adds `x4freeE L/R`; `pos4F_injective` uses `d 4 = 0` from `coeff_denZ_x4free`; correct.
- `matchW_sound` (841): only `W.add` is matched against `PT.node`; `W.sub` at a node is rejected (returns false). Leaves accept `var`, `nat`, `neg`, `mul`, `pow (var i) e` only.
- `toC_isHomogeneous` (716) and `toC_ne_zero_of_evalPT` (734): straightforward, no hidden hypothesis.
- Zero leaves: `homB` accepts `c = 0` at any degree; `firstDeg` skips them; `degE = none` handled by `denZ_eq_zero_of_degE`. Correct.
- Still to verify: `#print axioms` (test module not buildable here; will rely on PR statement + compile Kronecker.lean once), denotation shape vs Verify/accepted modules, check.py translation vs `W`.

## (2) Denotation shape vs Verify.lean / accepted modules
- `Poly := MvPolynomial (Fin 6) ℂ` (Verify.lean:12) is what `toC`/`denC` target; `directional`/`divergenceNumerator` are `HodgeTestbed` defs and the test's `e1`/`e2` match `normalized_identity`'s `h1`/`h2` shapes (with `(1 : Poly)` instead of `((1:ℕ) : Poly)`, harmless in a test).
- Accepted modules (data/m3-coupled/Witness.lean) write monomials as `C (σ := Fin 6) (c : ℂ) * X 0 ^ a * X 1 ^ b` (left-assoc, nonzero exponents only, numeral coefficient). `toC` writes `C (↑c : ℂ) * (X 0 ^ e0 * … * X 5 ^ e5)` with `Int.cast` and `^ 0` factors. Same type and vocabulary, NOT the same syntactic shape. README "Mechanism" says "the same `C c * (…)` shape as certificate modules": nit, soften (the Pitfalls section already admits `Int.cast`/`push_cast`).

## (3) Witness binding `W`/`denW`/`matchW` vs scripts/proposals/check.py
- check.py `polynomial()` (lines 166-185) emits `(MvPolynomial.X (R := ℂ) (k : Fin 6))`, `(n : HodgeComputer.Poly)`, `(-t)`, `(a ^ (k : ℕ))`, `(a op b)` for `+ - *`, and `(a * MvPolynomial.C (σ := Fin 6) (1 / n : ℂ))` for `/`. `W` mirrors the first five; there is no constructor for the `C (1/n)` case (README says scale to integers: fine).
- `matchW` only descends through `W.add`; `W.sub` is rejected at nodes and in leaves. The instance strings (instance.rec) and the accepted witness.json use binary `-`, so `matchW` can bind only witness strings deliberately emitted as sums of signed products with explicit balanced parentheses (Python's `ast.parse` left-associates unparenthesized sums: a left-deep `W.add` tree would force a left-deep `PT`, defeating `psiRelN`). README states the convention but not the two failure modes: should-fix, document "binary `-` and `/` are not matched; unparenthesized sums give a left-deep tree".
- The protected `HodgeSubmission.F`/`A` come from the fixed instance strings (with `-`), so `matchW` cannot bind them; the README's step 1 (`F = C (1/s') * toC tF'`) leaves that binding implicit: should-fix, say that `F`, `A` (small) are bound with the accepted `funext; norm_num; ring` template.
- The `rfl` claim `HodgeWitness.D = denW w` is NOT exercised anywhere in the branch (KroneckerTest only `#print axioms matchW_sound`). `(n : Poly)` for n ≥ 2 unfolds to `Nat.cast n` (instOfNat), so `rfl` is plausible; for coefficients `0`/`1` written as numerals it would not be (Zero/One instances vs `Nat.cast`), which the convention `X0**6` / `-X0**6` avoids. Testing this now with a scratch file against the exact check.py output.

## (4)-(5) Process / docs observations so far
- Verifier: sandbox.py builds `lake --no-cache build HodgeComputer.Verify *imports` by name (consistent with PR claim); policy: maxBytes 16 MiB, maxFiles 64, 3600 s, 14336 MB. README budget could mention the 16 MiB proposal cap (nit).
- check-library.py: any non-.lean change marks ALL modules affected (line 47-48) and then every accepted application is rechecked (line 70 includes Verify). PR.md says both "rechecks all accepted applications in Docker" and later "check-library.py rechecks no protected application": contradictory; should-fix the second sentence.
- MAINTAINING.md:13: "Reusable constructions belong under Constructions/ when they support concrete applications … Extract accepted machinery … check an application against the protected target before merging library changes." New `Reflection/` directory is undefined there; PR should state the placement rationale or add one sentence to MAINTAINING.md. No accepted application exists; local proposal run (184 s) is prose-only, no receipt.
- Attribution: CONTRIBUTING.md has no trailer rule; README phrasing "prepared with Claude Code assistance for `kejace`" matches the PolynomialCurves precedent ("prepared with Codex assistance for"); Co-Authored-By trailer has one precedent on main.
