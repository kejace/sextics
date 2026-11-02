# Referee report — lane `polarized-k3-remainder`, wave-2 exclusion draft (2026-11-02)

Referee: Claude (Anthropic) for `kejace`; no human review; nothing Lean-verified; no certificate claimed.  Inputs: draft
post, excl report, packet `wave2-excl/`, construct/search reports, aggregator status (lane row) and verification log
(19 entries), dossier, discussion #115 (two posts by `1337isnot1337`).  No adversarial verify-lens file exists for this lane,
so "reconciliation" below is against the aggregator's log.  Verdicts: **VERIFIED** = replayed by the referee (own code, or
the packet script re-run, exact over Q where stated) with the transfer to characteristic 0 checked; **PLAUSIBLE** =
reproducible computation that the referee did not replay, or modular-only, or a sound sketch; **UNSUPPORTED** = the
wording is not justified by the evidence.  Referee artefacts: `wave2-excl/referee/` and `wave2-excl/README-referee.md`.

## Claim table

| # | Claim (draft / report wording) | Evidence | Referee replay | Verdict |
|---|---|---|---|---|
| 1 | `[A]` injective on `(R/J)_d` for `d ≤ 9`, both parities (E1-admissible `D` of degree ≤ 9 lie in `J(F)`) | e1table (2^31−1, 2^31−19); aggregator 10^9+7 | own tensor-frame code: exact over Q (ranks 1…1246 full) + 4 primes | **VERIFIED** (exact) |
| 2 | Degree 6 for `A`: L-relaxation, block L-ranks 603/558/396/260/396/558, two-step square elimination ⇒ `D ∈ {0, cG}` | `relax-full-d6-*` | `relax.py 6 65521 --exact-constraints --squares`: 35/35 killed over Q, same L-ranks; necessity of the relaxation and the Euler gauge re-derived by hand | **VERIFIED** |
| 3 | Degree 7 for `A` replicated at p = 65521 (L-ranks 1635/1602/1635/1464/1287/1464, constraint rank 1779, 9 empty charts, 42 variables, GB DIM 0/1194/309) ⇒ `deg D ≥ 8` | lens log R13; thread DEGREE_SEVEN.md (p = 101) | not replayed (≈ 30 min under the load cap); the Z_(p)-flatness/properness transfer and the gauge-fixing checked line by line | **PLAUSIBLE** (two independent implementations, two primes; transfer sound) |
| 4 | Degree-8 relaxation for `A` regenerated: 315 parameters, rank 1050, no squares; contains `f(U)h` | `relax-full-d8-p65521-local.json` | not replayed (mod p; matches thread) | **PLAUSIBLE** |
| 5 | `D = f(U)U0^2` fails exactly for `A` in all three blocks (756/757, 828/829, 516/517) | `fixedD-koszul.jsonl` | `fixedD_koszul.py` re-run: identical ranks over Q | **VERIFIED** |
| 6 | `[A_K3]` injective on `(R/J)_d` for `d ≤ 7` | e1table-AK3 (2^31−1, 2^31−19) | own code: exact over Q, dims up to 666 | **VERIFIED** (exact) |
| 7 | `A_K3`, degree 5: relaxation kills all 6 parameters | `relax-k3-d5-*` | re-run: 6/6 over Q | **VERIFIED** |
| 8 | `A_K3`, degree 6: survivors `a_i = c_iU_i`, residual `(c0−c2)^2, (c0−c2)(c1−c2), (c1−c2)^2` ⇒ `D ~ f(U)`; `f(U)` fails exactly (136/137); control `D = G` passes | `relax-k3-d6-*`, residual .sing (mod p), `fixedD-koszul.jsonl` | re-run: 32/35 over Q; exact residual span over Q equals the stated three quadrics (`residual-k3-d6-exact.out`); 136/137 and the control reproduced | **VERIFIED** |
| 9 | Theorem: no denominator of degree ≤ 6 serves `A_K3` (over C, any parity) | 6–8 | follows from 6–8 with the hand-checked normalization/gauge | **VERIFIED** |
| 10 | `π_K3(served(f(U))) = π_K3(served(U0)) = π_K3(served(V0)) = 0` exactly over Q | `servedK3-exact.jsonl` | `servedK3.py --exact` re-run: class dim 0 for all three; the character-grading argument (μ6-homogeneous `D` ⇒ served space splits, K3 block = character 0) checked | **VERIFIED** |
| 11 | Hyperplanes `λU + μV`, `λμ ≠ 0`, fail E1 for every nonzero K3-block class | tensor identity | checked: bidegree pieces `B'⊗μV`, `λU⊗B''` vanish only if `B' = B'' = 0` | **VERIFIED** (formal) |
| 12 | Degree 7 for `A_K3` open by relaxation (rank 935/7260; 19-target 1092; no squares; 84-dim tangent spaces along `f(U)ℓ`) | `relax-k3-d7-*`, `multi-k3-d7-*`, zstruct | not replayed (mod p statements about the method) | **PLAUSIBLE** |
| 13 | `f(U)U0`, `f(U)V0` fail exactly for `A_K3` (377/378) | `fixedD-koszul.jsonl` | `f(U)U0` re-run: 377/378 over Q; `f(U)V0` not replayed | **VERIFIED** (`f(U)U0`), lens-exact for `f(U)V0` |
| 14 | Sum principle + dichotomy: no product of supports of degree ≤ 6 serves `A`; a legal octic has `X ∩ V(D)` irreducible and reduced or `(irreducible degree-7 K3-carrying support) ∪ (hyperplane section)`; the thread's `f(U)h` points are dead | logic note (R2, R4), dossier §4.3 | hypotheses listed (End_Hdg(T) = Q [thread], η ∉ NS_prim⊗C, book's prescribed-support theorem and (3,1)-additivity over components); not re-derived | **PLAUSIBLE** (conditional, theoretical) |
| 15 | Block cones `g(U)`, `g(V)`, products `g(U)h(V)` excluded in every degree (Hodge-theoretic) | R9 / draft 2b | sketch only; exceptional divisor over `C×C` gives only `[C]⊗[C]` in the Künneth `H^2⊗H^2` part (blow-up formula), but the transport through the μ6-quotient/book projectors is unwritten; the lens itself flagged it | **UNSUPPORTED** as a theorem; plausible sketch |
| 16 | Non-Jacobian E1 kernel of `[A_K3]` at `d = 8` is exactly `K4⊗K4` (81, char 4) ⊕ `{ℓ(U)ℓ'(V)(η(U)+η(V))}` (9, char 1); `dim K4 = 9`; grows to 306 at `d = 9` | R5–R7, R10–R11; aggregator K4 at 998244353 | own code: kernel 90 (d = 8) and 306 (d = 9) exactly over Q; `dim K4 = 9` exact; both blocks in the kernel by formal identities (`deg ℓη^2 = 13`) | **VERIFIED** (exact) |
| 17 | Ten tested kernel elements (8 in `K4⊗K4`, 2 in the char-1 block) fail E2 for `A_K3` (13572×2430, rank 2206, mod 2^31−1) | lens log | not replayed; modular screens | **PLAUSIBLE** (screen) |
| 18 | Random μ6-homogeneous members of `J_7`, `J_8`, `J_9` fail exactly (rank gap 1) for both targets | `fixedD-koszul-struct-d*.jsonl` | not replayed (lens exact) | **PLAUSIBLE** |
| 19 | Degree-8 μ6 sectors for `A`: char 0 → cones `g(U)`, `g ∈ J(f)_8` (6 residual quadrics); char 1 relaxation essentially empty; char 2 → `Σ U_i f_i(U) q_i(V) + g(V)` | `relax-full-d8-char*` | not replayed; sector-restricted, mod p | **PLAUSIBLE** |
| 20 | `A = A_K3 + A_mix` with both parts nonzero modulo `J(F)` (construct lens) | construct log (2^31−1) | own code, exact (bigrading; η ≠ 0 in `R_{f,6}`) | **VERIFIED** |
| 21 | `Aut(X) ⊇ D_6 = ⟨h, τ⟩`, `h(U,V) = (U, ηV)`; six planes `V = λU`, `λ^6 = −1` | construct A1–A2 | `f(ηV) = f(V)` and `f(U) + f(λU) = (1+λ^6)f(U)`: immediate | **VERIFIED** (trivial) |
| 22 | `Aut(X) = D_6` exactly; averaging-ideal exhaustion: no quotient `X/K`, `K ≤ D_6`, serves the K3 sector (ages; junior elements only `h, h^5`) | construct A1, A4–A5; probe `Aut(C_f) = 1` | not re-derived (depends on the probe's msolve length-6 scheme at 65521 + Deligne–Mumford, and on the age table) | **PLAUSIBLE** |
| 23 | Reduction: under `End_Hdg(T) = Q` and `η ∉ NS_prim⊗C`, `HasCertificate(F,1,A)` ⇔ GHC (coniveau 1) for `Λ^2T ⊂ H^4(S_2×S_2)` | dossier §4.3 (sketch), construct A7 | hypotheses explicit; ⇒ uses irreducibility of `Λ^2T` (Zarhin ⇒ MT = SO), ⇐ uses the book's existential supports for the mixed sector | **PLAUSIBLE** (reduction, not a result) |
| 24 | Construct Part B: mixed-sector chart reduction, scalar-conversion identity, F_5 seed with invertible 13×13 Jacobian, Hensel lift mod 5^6 | construct log | not refereed | **PLAUSIBLE** |
| 25 | "zero modular served-space projection is only evidence" (lens's own rigor correction) | R3 | correct (rank mod p ≤ rank over Q goes the wrong way for vanishing); exact mode used downstream | **VERIFIED** (bookkeeping) |
| 26 | E1 table for `A` re-verified at a second prime (2^31−19), both parities | lens log | own exact computation supersedes | **VERIFIED** |
| 27 | Aggregator: plain mixed-parity Hessian relaxation at `d = 6` is not m-primary (`c5_2` lacks a pure power) — no contradiction with the degree-6 exclusion | verification log | the degree-6 exclusion uses the finer L-relaxation `(G)+J^2+(R_ijk)`; consistent | **VERIFIED** (consistency) |
| 28 | `D = f(U)`'s support (join of the two copies of `C`) "carries no K3 class at all, as Chapter 3A predicts" | R6 | the exact computation is row 10; the geometric gloss is the book's | **PLAUSIBLE** (gloss) |

## Reconciliation
Aggregator verification log (19 entries at 10^9+7 / 998244353): Hilbert dims, `A ∉ J`, E1 kernels 0 for `d ≤ 9` both parities,
dense fixed-D `F_X0` inconsistent at `d = 5`, `K4` of dimension 9, the odd kernel block identity, and the non-m-primary plain
Hessian relaxation — all consistent with the packet and with this pass.  The excl lens's own corrections (modular served-space
projections demoted to evidence; exact mode added) are adopted.  Disagreements with the draft: row 15 (theorem → sketch) and the
conditional wording of row 14 (hypotheses must be displayed; degree 7 for `A_K3` is open).

## Counts
28 rows.  **VERIFIED 16** (rows 1, 2, 5, 6, 7, 8, 9, 10, 11, 13, 16, 20, 21, 25, 26, 27; of these, 25–27 are bookkeeping/consistency rows,
so 13 are mathematical statements replayed by the referee), **PLAUSIBLE 11** (rows 3, 4, 12, 14, 17, 18, 19, 22, 23, 24, 28; the `f(U)V0`
half of row 13 is lens-exact only), **UNSUPPORTED 1** (row 15 as worded).  Strongest verified statement: *the K3-sector class `A_K3 = (η(U)−η(V))/64` admits no
denominator of degree ≤ 6 over C (E1 injective through degree 7 exactly over Q; degree 5 and 6 relaxations exact; `f(U)` and the
residual quadrics exact), and `[A]` is injective on the Jacobian ring through degree 9 exactly over Q*; the thread's `deg D ≥ 8` for
the instance target rests on two independent implementations at two primes with a checked transfer argument.
