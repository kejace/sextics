# Referee report — `diagonal-permutation`, wave-2 exclusion draft (2026-11-02 UTC)

Referee: kejace's team, Claude (Anthropic) assistance; no human review; nothing Lean-verified; no certificate claimed.
Inputs: draft post `wave2-diagonal-permutation-draftpost.md`, report `wave2-diagonal-permutation-excl.md`, packet
`results/diagonal-permutation/wave2-excl/` (+ remote receipts, copied to `results/diagonal-permutation/wave2-referee/remote-receipts/`),
verify-lens notes, aggregator status v6.3 and verification log, dossier, discussion #177.
Referee replays: fresh prime **999999937** (never used in this lane), independent code (`wave2-referee/ref_*.py`, no `hodge`/`w2lib` import);
exact rational arithmetic where stated.  Verdicts: VERIFIED = receipt located + replayed (or exact) + transfer argument checked;
PLAUSIBLE = modular evidence only / not replayed; UNSUPPORTED = no receipt or computation unfinished.

| # | claim in the draft | receipt | referee replay | verdict | wording |
|---|---|---|---|---|---|
| 1 | Lemma "E2 modulo J": `Ψ(D;u) := D div u − u(D) ∈ J` for every certificate; `Ψ mod J` independent of `u` and of `D mod F·S`; `Ψ ≡ −A·wᵀHess(F)w` for `D = w(F)` | README §1.1–1.2, `w2lib.py` docstring | rederived by hand from the contract identity (formal/README.md), incl. the cancellation of the `H`-terms and the Euler-gauge invariance | VERIFIED | keep |
| 2 | `d ≤ 4`: `ker[A] = 0` both parities (3 primes), `J_d = 0` ⇒ no denominator of degree ≤ 4 | `e1k_e1kF.json` (P1, P2), wave-1 (P3) | `ref_e1_p999999937.jsonl`: kernels 0, good reduction in degrees d, d+6; `ref_quick.json`: `A ∉ J_6` (degree 0) | VERIFIED | keep; add d = 0 |
| 3 | `d = 5,6,7`: `ker[A] = 0` both parities ⇒ every admissible D is `w(F)`, `deg w_i = d−5` (mixed parity included) | `e1k_e1kF.json` (P1, P2) | `ref_e1_p999999937.jsonl`: 0/0 at d = 5, 6, 7 with good reduction | VERIFIED | keep |
| 4 | `d = 5`: the 21 relaxation equations span all quadratic monomials ⇒ `w = 0` | `relax_d5_mixed.ms` records only the 21 monomials (the result, not the input) | `ref_hess5.json`: exact rational rank of the 21 second partials `F_ij` = 21 (fmpq_mat); same as aggregator `hessrank5.py` | VERIFIED (exact) | keep; cite the exact rank |
| 5 | `d = 6`: msolve reduced bases zero-dimensional in the Y-even (17/134), Y-odd (18/153) and mixed (35/404) sectors ⇒ no degree-6 denominator of any parity | `relax_d6{e,o,m}.{ms,gb,json}`, `singular_d6.out`, `relax_d6*_db0.sing` (P1); aggregator `hessrelax6` (P1, own construction) | `ref_relax6_p999999937.json`: rebuilt from scratch, identical counts 17/134, 18/153, 35/404; msolve `-g 2`: pure powers of all variables, colength 39/39/924; Singular `dim 0`, `vdim 39/39/924` (`ref_singular_check.jsonl`) | VERIFIED | keep ("proved"); over C by base change |
| 6 | `d = 7`: relaxation pending (relax7F) | `relax_d7e.json`: 60 unknowns / 480 independent quadrics built; msolve killed (`status −9`) after 7079 s; odd/mixed not started | none (too heavy locally) | UNSUPPORTED (open) | say "open; system built, Gröbner basis not obtained" |
| 7a | `d = 8`: `ker[A]` on `R_8^+` is one line (D8), on `R_8^-` zero | `e1k_e1kF.json` (P1, P2), wave-1 (P3), aggregator (P1) | `ref_e1_p999999937.jsonl`: 1 / 0 with good reduction | VERIFIED | keep |
| 7b | `A·D8 ∈ J` over Q with an exact `u8` | `d8_d8G.json`: `u8` lift **failed** after 4 primes | `ref_e1`: `A·D8 ∈ J mod 999999937` for the 414-term wave-1 rational D8 (6th prime) | PLAUSIBLE (screened) | "verified modulo six primes; exact lift not achieved" |
| 7c | D8 relaxation `Ψ(D8;u8) ∉ J` ⇒ no certificate with `D = D8 + gF` | `d8_d8G.json`: `in_J: false` at P1 and P2, good reduction in degree 16 | `ref_d8_p999999937.json`: Ψ(D8;u) ∉ J_16 mod 999999937 (good reduction in degrees 14 and 16; E1 nullity 954 = Koszul; rank J_16^- = 9684; u recomputed from scratch); transfer argument (README §1.3) checked by hand, see README-referee §4 for the exact scope | VERIFIED (scope: canonical representative + F-multiples only) | state the scope: excludes the canonical representative and its F-multiples, **not** the family `D8 + w(F)` |
| 7d | reduced E2' system of D8 inconsistent at three primes "with rank certification" | `d8_d8G.json`: inconsistent at P1, P2, P3 (rank 2286, nullity 159) but `certified_exact_nullvectors: 0` | none | PLAUSIBLE (screened) | "inconsistent at three primes; rank not certified" |
| 7e | full family `λ D8 + w(F)` (141 unknowns, 468 equations) "pending/screened" | no receipt (no `relax_d8*_D8` files) | none | UNSUPPORTED (open) | "not computed" |
| 8 | `d = 9`: `ker[A]` 21-dimensional in each parity; E1 stops excluding at degree 9 | `e1k_e1kF.json`: 21 / 21 at P1; sibling at P3 | `ref_e1_p999999937.jsonl`: even 21 (good reduction in degrees 9, 15); odd 21 (good reduction True); exact lower bound 21 from the sector count (`ref_sectors.json`) ⇒ equality over Q | VERIFIED | keep |
| 9 | sum principle: `Ann_6(D) = 0` (mod P1) for 49 structured factors and for the Hessian forms of the listed vector fields ⇒ products of these factors serve only the zero class | `ann_ann-local.json` + log (49 + 13 entries, rank 208 each); good reduction in degrees 8..12 from `e1k_e1kF.json` | `ref_ann_p999999937.jsonl`: 11 of them (p2X, pXY, VdmX, Hesse, sX2Y2, A, sumT, B2, p6X, H(∂0,∂0), H(E_X,E_X)) all zero with good reduction; sum principle argument checked against docs/01 (prescribed-support theorem, arbitrary parity) + weight argument for additivity | VERIFIED | keep; cite docs/01 |
| 10 | equivariant sector dimensions, `[A]: (par,σ) → (1−par,−σ)`, forced kernels (none for d ≤ 8; 21 in `(±,+)` at d = 9; 170/124/152/109 at d = 10; 275/235 at d = 11; 361/316/342/306 at d = 12) | README §4, `w2_reptheory.py`, `e1k_e1kS.json` | `ref_sectors.json`: table and forced kernels reproduced exactly from the equivariant Koszul series; `e1kS`: d = 9 `(±,−)` kernels 0 and `(±,+)` 21 at P1, d = 10 sectors 170/124/152/110 | VERIFIED | keep; the odd d = 10 total is 262 = 261 forced + 1 |
| 11 | degree-9 kernels are exactly the forced ones and `(01)`-invariant modulo J | `e1k_e1kS.json` (P1) | upper bounds from P1 and from the referee's prime, exact lower bounds from the sector count | VERIFIED | keep |
| 12 | wave-1 re-verifications (symmetry group, characters `39+31+69·std` etc., Hilbert function) | `reptheory.json`, verify lens rows 1–6, 18–19 | parity-refined Koszul series reproduced (`ref_common.koszul_series`) | VERIFIED | keep |
| 13 | "What is screened only" list in the draft (`[pending]`) | — | — | — | replace by: 7b, 7d above; families d = 8..12 and degree-9 members **not run** |

Counts: VERIFIED 12 (rows 1, 2, 3, 4, 5, 7a, 7c, 8, 9, 10, 11, 12), PLAUSIBLE 2 (7b, 7d), UNSUPPORTED/open 2 (6, 7e).  Strongest verified statement: no certificate with deg D <= 6 for any homogeneous D (all Y-parities, over C).

## Reconciliation with the verify lens (`wave2-diagonal-permutation-verify.md`, 00:40 UTC)

* Agree: the attack-1 headline ("no certificate with deg D ≤ 7") was an overclaim at that time; the excl lens's own report says so, and its new degree-5/6 closure now proves `deg D ≥ 7` for arbitrary D — the verify lens had not yet seen it.
* Agree: modular inconsistencies (18 Jacobian denominators, D8 E2') are screens, not proofs (its rows 9, 10).
* Addition not in either lens: the mod-p *relaxation* failure of D8 **does** transfer to Q without an exact lift (regular-sequence/Cramer argument, README §1.3), but only for the canonical representative and its F-multiples; the README's case (a) is misworded (it should say "the rational kernel is zero").
* Minor: the packet's degree-5 receipt records the outcome, not the inputs; the exact rank computation (aggregator, referee) is the proof.

## Gaps fixed by `README-referee.md`

Remote receipts not in the packet (copied); hard-coded paths; `fam_*`/`d8_u8_exact`/`d9_*` receipts that do not exist; `good_reduction: null` in the sector receipts (how to confirm from the data); `Singular --version` hangs; exact scope of the D8 exclusion.
