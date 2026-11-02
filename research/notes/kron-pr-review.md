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
