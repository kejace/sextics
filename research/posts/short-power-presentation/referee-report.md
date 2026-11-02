# Referee report — lane `short-power-presentation` (wave 2, exclusion/rigor lens)

Referee: Claude (Anthropic) for `kejace`, 2026-11-01/02. No human review. Nothing Lean-verified. No certificate claimed.
Full ledger with the transfer arguments: `research/results/short-power-presentation/wave2-excl/README-referee.md`
(commit `57cd9585738240ebc5bda3be15729cc5ad92e217` of the research branch; both links in `comment.md` returned HTTP 200).

## Counts
* VERIFIED: 11 (V1–V11) — of which replayed by the referee at the fresh prime 1073741741 or exactly over Q: E1 kernels
  d = 5, 6, 7 (both parities, Hilbert-checked); exact `D7^±` certificates; exact `L*D7^+`, `M*D7^- in J_8`; the odd cone
  representative; relaxations d = 5 (e/o/m) and d = 6 (e/o); the mixed d = 6 functional systems with 114/114 exact
  dependency certificates on both sides; `A not in J`; Hilbert reference values.
* PLAUSIBLE (modular evidence only / incomplete chain): 8 (P1–P8) — everything at degree 7 and 8, the fixed-D
  inconsistencies, the λ-only subfamilies, `dim K_8^+ = 6` as an exact value, the verify lens's "characteristic 0" run.
* UNSUPPORTED as stated: 5 (U1–U5) — the attribution of the mixed exact certificates to job `exclB-spp-relax-mixed` (no
  receipt existed; the referee replay produced them), "degree 7 partly proven", "covers D7^+ + J_7 once it finishes",
  the draft's TODO placeholder, the construct note's inverted bound direction.

## Strongest verified statement
Every certificate denominator for `(F, A = X0^3 Y0^3, q = 1)` has degree at least 7, with no parity assumption; and the
degree-7 E1 kernel is exactly one-dimensional in each parity, with exact representatives `D7^+ in Q[X0,X1,X2,M]` and
(modulo `J_7`) `D7^- in Q[Y0,Y1,Y2,L]`, both with exact numerator certificates.

## Transfer arguments and parity coverage (summary; details in the README ledger §2)
* E1: `dim_Q K_d/J_d <= dim_p K_d/J_d` needs `rank_p J_{d+6}` and `dim_p J_d` at the Hilbert values (recorded); lower
  bounds need exact elements (present at d = 7).
* E2 relaxation: re-derived; the mod-p quadric ideal is the reduction of the rational one when `rank_p J_{2d}` is the
  Hilbert value (recorded in `parts`); the Hilbert-function inequality then gives `V(I_Q) = {0}` in the Euler slice.
* Mixed parity: odd part has no cross terms (checked); annihilator functionals are reductions of rational ones iff the
  exact dependency certificates exist — they did not until the referee replay; now 114/114 on both sides. Independent
  confirmation: aggregator's full-system m-primary leading ideal at 10^9+7.
* Block elimination at degree 7: logic sound, but every step needs the same exact-rank certificate; the certification
  job ended INCOMPLETE (lift failed at 12 primes), so degree 7 is single-prime evidence only.
* Parity: `J` is τ-stable, so `A D in J` splits; the only symmetry is τ, so no further character splitting; coverage is
  complete through degree 6 (pure and mixed).

## Reconciliation with the parallel verify lens (`wave2-short-power-presentation-verify.md`)
Agrees on all 25 rows; the referee adds: row 19's pending exact `D7` check is done (exact over Q); row 25's objection to
"dim = -1 for λ = 1" stands; the verify lens's own "characteristic 0" Singular run (row 20) has an artifact on the remote
but its ring characteristic was not confirmed here (not needed: the mod-p + Hilbert route proves V8).

## Remote state at closing (UTC 04:35, 2026-11-02)
`excl_residual_d7e_p0`: msolve exit 137 at F4 degree 6 (2533891 x 2261343). `exclB-spp-certify7e`: INCOMPLETE for both
steps attempted (λ, w3). `exclB-spp-blocks2`: d = 7 odd residual ended without verdict; d = 7 mixed λ^± eliminated mod p,
residual 120 unknowns / 636 quadrics; d = 8 even started. Local compute by the referee: six small jobs, each under 6 GB.

## Answer to the side question "what is diamondjm"
`DiamondJM` is a GitHub account (display name Josh Diamond, created 2018, five public repositories; profile fetched with
HTTP 200). It does not appear in this repository, this lane's thread, or any lane note; the instance author and the only
thread participant is `benediamond` (Benjamin Diamond), a different account.
