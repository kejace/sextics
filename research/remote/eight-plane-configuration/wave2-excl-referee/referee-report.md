# Referee report — `eight-plane-configuration`, wave-2 exclusion draft (2026-11-01, local time EST)

Referee: Claude (Anthropic) for `kejace`, independent code in `research/results/eight-plane-configuration/wave2-excl-referee/`
(`ref_lib.py` parses `instance.rec` with sympy; ranks with python-flint 0.9.0; msolve 0.6.5; Singular 4.3.2). Receipts: `referee.log.jsonl`.
Verdicts: **VERIFIED** = logic checked and independently recomputed (exact over ℚ, or a new prime 10⁹+9 with the no-drop check, or a
second engine); **PLAUSIBLE** = receipt present, modular only / not independently rebuilt; **UNSUPPORTED** = no receipt or wrong as stated.

## Claim table (draft post items 1–6 and the packet statements S1–S6)

| # | claim (draft) | receipt in packet | referee replay | verdict |
|---|---|---|---|---|
| 1a | E1 kernel `{D ∈ S_d : A·D ∈ J} = J_d` for d ≤ 9, both parities (S2) | `e1ladder.2147483587.jsonl`, wave-1 `e1kernel*.jsonl` (2147483647, 2147483629), `w2x_hilbert.py` no-drop check | `ref_e1.py` at p = 1000000009, d = 4..9 both parities: `dimV_p = dim J_d` (0, 3, 18, 63, 168, 378) and `rank_p J_{d+6} = exact` in every row | **VERIFIED** |
| 1b | rigor argument `dim_ℚ ker ≤ n − rank_p[A|J] + dim_ℚ J` needs and has `rank_p J_{d+6} = dim_ℚ J_{d+6}` | §5 of excl report, verify row 13 | re-derived; saturated-lattice argument also covers the d = 8 use of the mod-p kernel basis | **VERIFIED** (logic) |
| 1c | d = 10: kernel mod J = 294 (even) / 261 (odd) = forced `R_10^± − R_16^∓` | ladder d = 10 rows (ranks 9684/9714 = exact `J_16^∓`), wave-1 remote row | own Hilbert numbers: `J_16^- = 9684`, `J_16^+ = 9714`, `R_10 = 762/744`, `R_8 = 483/468`; forced bound = computed ⇒ equality over ℚ. Own d = 10 rank run at 10⁹+9: even sector kernel mod J = 294 with `rank_p J_16^- = 9684` exact (`ref_e1_d10.out`; odd sector still running at hand-off) | **VERIFIED** (logic + two primes) |
| 2a | every certificate denominator of degree ≤ 9 is `D = w(F)`, coprime iff `w ∉ {h·Euler}` (S3/R3) | excl §5 | re-derived (syzygy coefficient degree d − 5 < 5) | **VERIFIED** (logic) |
| 2b | Koszul-reduced E2 is *equivalent* to the certificate problem for `D = w(F)` (R4), μ-terms cancel | excl §5, verify row 18 | re-derived from `formal/README.md` (E = F²Q, q = 1): `N = Aw + μE/6 + K_c`, `D·div(μE/6) − (μE/6)(D) = Dμ` | **VERIFIED** (logic) |
| 2c | Hessian condition `A·Hess_F(w,w) ∈ J_{2d}`; for d ≤ 7 `Hess_F(w,w) ∈ J_{2d−6}`; gauge-invariant | excl §5 | re-derived (`w(w(F)) = (∇_w w)(F) + Hess`, `Hess(E,v) = 5v(F)`); needs `Ann_R(A)_{2d−6} = 0` both parities (= 1a) | **VERIFIED** (logic) |
| 3a | d ≤ 4 excluded (`J_d = 0`) | 1a | 1a | **VERIFIED** |
| 3b | d = 5 excluded, any parity (`Hess = 0` ⇒ `w = 0`, rank 6 per sector, 21 mixed) | `hess.jsonl` d=5, `hess_mixed.jsonl` d=5 | exact fmpq rank of the 21 quartics `F_ij` = 21 (`ref_exact.py`) | **VERIFIED** (exact) |
| 3c | d = 6 pure parity excluded: Hilbert function [1,17,18,2,0] / [1,18,18,2,0] | `hess_d6_par*_p{2147483587,2147483647}.lm`, `hess.jsonl` | own build at 10⁹+9: msolve m-primary, HF [1,17,18,2] / [1,18,18,2], colength 38/39; Singular dim 0, vdim 38/39 | **VERIFIED** |
| 3d | d = 6 mixed parity excluded (35 vars, 405 quadrics, quotient 923) | `hessmixed_d6_p2147483587.{ms,lm}`, `hess_mixed.jsonl`; aggregator at 10⁹+7 (`wave2-aggregator-hessrelax6.md`) | own build at 10⁹+9: 426 quadrics (405 independent), m-primary, HF [1,35,225,400,225,36,1], colength 923; Singular dim 0, vdim 923 | **VERIFIED** |
| 3e | transfer "(S/I)_k = 0 mod p ⇒ over ℚ" | excl §5 | Macaulay-matrix rank argument checked; requires `rank_p J_6 = 18` (checked at every prime used, incl. 10⁹+9) and p ∤ 7·11·13·17 | **VERIFIED** (logic) |
| 3f | six partials inconsistent over exact ℚ (102 vs 103) | remote `excl-epc-batch2a-partials-exact.out` (parity-restricted) ; excl-vm 210/211 | own integer-model system without parity restriction, `dF/dX0`, `dF/dY0`: 3003 × 216, rank 210, augmented 211 | **VERIFIED** (exact; superseded by 3b anyway) |
| 3g | 11 named `J_6` elements inconsistent (two primes) | lane log `a3-koszul-reduced-E2` rows | not replayed; single/two-prime modular inconsistency is evidence only (verify row 16); superseded by 3c/3d | **PLAUSIBLE** (superseded) |
| 4a | d = 7 pure: 483 independent quadrics in 60 variables per parity; full Gröbner pending | `hess.jsonl` d=7, `hess_d7_par*_p2147483587.ms` (`.lm` empty: OOM rc 137) | counts consistent with `dim R_8^+ = 483`; not rebuilt; no result exists | **PLAUSIBLE** (pending) |
| 4b | d = 7 mixed: 951 quadrics / 120 variables, not attempted | `hess_mixed.jsonl` d=7 | — | **PLAUSIBLE** (no exclusion claimed) |
| 4c | d = 7 gradient fields `∇G` excluded (linearization rank 406) | `hess_sub.jsonl` grad rows; aggregator at 998244353 | `ref_grad7.py` at 10⁹+9: rank 406 = dim Sym², both parities, `rank_p J_8 = 168` exact | **VERIFIED** |
| 4d | d = 7 subfamilies arr (29), blockX (708), blockY (439), coefX (203), coefB0 (13) excluded | `hsub_d7_par*_*.lm`, `hess_sub.jsonl` | msolve replays + independent leading-ideal reading: coefB0 (13), arr (29), coefX (203) m-primary; blockX/blockY not replayed; family encoding per `w2x_hess_sub.py` not rebuilt | **VERIFIED** (coefB0, arr, coefX) / **PLAUSIBLE** (blockX, blockY) |
| 5a | d = 8 condition `Hess ∈ Ann_R(A)_10^+ + J_10^+` (468 quadrics, 156/159 vars), Gröbner infeasible | `hess_sub8_local.jsonl` | formulation checked (1c + saturated-lattice argument) | **VERIFIED** (logic) |
| 5b | d = 8 Hamiltonian fields `C·∇G`, G ∈ {K_i,H_i,diag_i,Q_iQ_j,Q_i²} (30 cases) and random quartics: linearization full rank | `hess_sub8_local.jsonl` (p = 2147483629) | not replayed (needs the wave-1 d = 10 kernel basis); transfer sound given the wave-1 no-drop at J_16 (checked by `w2x_hilbert.py`) | **PLAUSIBLE** |
| 5c | d = 8 coefB0 (28), coefX (4006), 12 single-component fields (740) excluded | `hsub_d8_par*_{coefB0,coefX,comp-i}_p2147483629.lm` | msolve replay of `hsub_d8_par0_coefB0`: m-primary, colength 28 | **VERIFIED** (coefB0) / **PLAUSIBLE** (coefX, comp) |
| 5d | 140 individual degree-7 fields inconsistent (rank 996/997, 1056 unknowns) | `receipt-batch2a-partial.jsonl` (`rigorous_rank_full: false`) | not replayed; single prime 2147483587, no rational dual ⇒ screening only, must not be called "excluded" | **PLAUSIBLE** (screened) |
| 6a | `Ann_R(D)_6^{odd} = 0` for 36 low-degree factors, "rigorous by the rank argument" | lane log `c-factor-annihilator` (p = 2147483587; **no-drop rank not recorded**) | `ref_factors.py` at 10⁹+9: all 36 zero and `rank_p J_{6+e} = exact` recorded; `F = ΣQ_iK_i`, `K_i = diag_i + (5/p_i)H_i` re-verified exactly | **VERIFIED** |
| 6b | sum principle `served(D₁D₂) = served(D₁) + served(D₂) mod J` ⇒ no product of these factors serves A | excl R6 (no proof given) | book Ch. 1 prescribed-support theorem (q = 1: poles in V(D) ⇔ denominator `D_red`) + Ch. 3 Gysin/weight description ⇒ served depends only on V(D) and is additive over components; with `A ∉ J` (S1) the conclusion follows. Dependence on these two theorems should be stated in the post | **VERIFIED** (logic, book-dependent) |
| S1 | `A ∉ J(F)` by mod-5 Fermat + saturation; "dim J_6 = 126" | — | exact rank 18 vs 19; the number is 36 (18 per parity), argument unaffected | **VERIFIED** (typo noted) |
| H | status line "deg D ≥ 7 for D of any Y-parity" | all of the above | — | **VERIFIED** |

Totals: VERIFIED 20 (incl. logic-only items), PLAUSIBLE 7, UNSUPPORTED 0 (one numerical typo in S1).

## Reconciliation with the verify lens (`wave2-eight-plane-configuration-verify.md`, rev. 19:39 local)
* Rows 12(a), 13, 16, 17, 18, 19, 28: agree. Row 16 is applied here to item 5d (sampling = screening).
* Rows 20–21 ("d = 6 WEAKENED to pure parity; mixed sector pending"): **resolved**, not a disagreement — the mixed sector is closed by three
  independent builds (lens 2147483587; aggregator 10⁹+7; referee 10⁹+9 with msolve and Singular). "deg D ≥ 7 for arbitrary D" is rigorous.
* "No sum-principle application in this lane" (logic point 3): outdated; item 6 uses it, and the referee finds it sound at q = 1 via the book's
  prescribed-support theorem and Gysin description, once the no-drop ranks are recorded (now done).
* Row 2 (denominator primes include 17): confirmed by the referee's own parse (lcm of denominators has factor 17²).
* The lens's queued recomputations (e1rank d = 6..9 at 4294967291, hessmac) did not appear in the mirrored log; the referee's runs cover them.

## Wording changes applied to the post
"140 fields … all inconsistent" → screened (single prime); item 1 "three primes / two code paths" → four primes / three code paths;
item 6 states its dependence on the two book theorems; S1 number corrected; degree-7 full run described as pending (OOM-killed, rerun in progress).
