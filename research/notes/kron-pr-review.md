STATUS: complete — verdict: needs fixes (docs/PR text only; no blocker in the Lean module)

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

## Verification receipts (this review)
- `lake build --no-build HodgeComputer.Reflection.Kronecker HodgeComputer.Reflection.KroneckerTest`: "All targets up-to-date" (oleans match committed sources).
- `lake env lean -j 2 -M 6000 HodgeComputer/Reflection/Kronecker.lean`: 30.4 s, maxrss 3.47 GB, zero warnings/errors. (`ulimit -v 12000000` makes Lean abort with "failed to create thread"; Lean's own `-M` cap used instead.)
- Scratch file /tmp/hodge/lean/kron-review/Scratch.lean importing both oleans: `#print axioms` for soundN, sound4, matchW_sound, toC_isHomogeneous, toC_ne_zero_of_evalPT, KroneckerTest.{small,e1,e2,denominator_ne_zero,denominator_homogeneous} = [propext, Classical.choice, Quot.sound] each.
- Binding test with the exact check.py `polynomial()` output for "-5*X0**6 + (3*X0**4*X1**2 + -X1**5*Y2)": `Dgen = denW wD := rfl` OK; `matchW_sound wD tD (by decide +kernel)` OK. `(1 : Poly) * X0^6 = denW (.mul (.nat 1) …) := rfl` FAILS (numerals 0/1 are Zero/One instances, not `Nat.cast`).
- Cited names `TransitiveCubicKronecker`, `TransitivePolynomialDigits` exist (Constructions/TransitiveCubic/Kronecker.lean, PolynomialDigitCheck.lean).

## Findings (sorted)
BLOCKER: none.
SHOULD-FIX
1. PR.md "How verified", last paragraph: "No accepted application imports the new module, so check-library.py rechecks no protected application" contradicts the earlier (correct) sentence; with README.md in the diff, scripts/check-library.py:47-48 marks all modules affected and :70 rechecks every accepted application. Reword.
2. Reflection/README.md "Use in a solve PR", binding paragraph: state the verified limits: `rfl` binds `(n : Poly)` only for n >= 2 (never emit `1*`/`0*`; `Q = "0"` is bound separately); `matchW` rejects binary `-` and `/`; unparenthesized sums are left-associated by `ast.parse`; the protected `HodgeSubmission.F`/`A` (fixed instance strings with `-`) are bound with the accepted template. Add a 3-term W-binding example mirroring check.py output to KroneckerTest.lean so CI exercises the claim (scratch example above works).
3. formal/README.md:60 + Reflection/README.md Budget + PR.md Why/Timings: capacity-study numbers (200-330 s, 185k fails, 435k/58 s/5.5 GB, d=10 3.1M/263 s/13.5 GB, proposal 184 s vs 275 s) have no receipt in the PR; link logs or a discussion post, or label as local unpublished measurements. "millions of products become checkable within the verifier's budget" should note the cited run used 13.5 GB of the 14 GiB Docker cap (use `sound4` for d >= 10).
4. Placement: MAINTAINING.md:13 defines only `Constructions/` and asks for an application against the protected target before merging. PR should justify the new `Reflection/` directory (proof infrastructure, not a construction), ideally with one sentence added to MAINTAINING.md, and acknowledge no accepted application yet (local proposal run has no artifact).
5. Kronecker.lean:747 `inductive W` shadows the bound variable `W : ℕ` used by `soundN`/`sound4` in the same namespace (README lists it as a pitfall). Rename (e.g. `WitnessTerm`) before any importer exists.
NIT
6. README "Mechanism": `toC` is not "the same shape" as certificate modules (`C (c:ℂ) * X 0 ^ a * X 1 ^ b`, no `^0`, numeral coefficients vs `Int.cast` and six factors); say "same vocabulary".
7. Docstrings: 20/165 declarations; `soundN`, `soundC`, `sound4`, `matchW_sound`, `E`, `W`, `denC` have none; some start lowercase.
8. `toZ_pdLeaf_{zero..five}` (172-236) are six copies of one `fin_cases`/`simp` proof; 66 lines > 100 columns.
9. README Budget: mention `.github/solve-policy.json` caps (16 MiB, 64 files).
10. formal/README.md: "linear in the number of witness terms" -> add "for fixed degree".
11. `#print axioms` at the end of KroneckerTest emits info on every library build; optional `#guard_msgs`.
Attribution style OK (matches PolynomialCurves README precedent; commit trailer has precedent on main).
VERDICT: needs fixes (documentation/PR text; Lean module itself is sound and lint-clean).
