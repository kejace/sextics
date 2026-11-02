# Referee report — `rank-one-biquartic` (discussion #118), wave 2 exclusion lens

Prepared for `kejace` with Claude (Anthropic) assistance; no human review yet; nothing Lean-verified; no certificate
claimed.  Date 2026-11-01/02.  Inputs: draft post `wave2-rank-one-biquartic-draftpost.md`, exclusion report
`wave2-rank-one-biquartic-excl.md`, packet `/tmp/hodge/results/rank-one-biquartic/wave2-excl/`, adversarial notes
`wave2-rank-one-biquartic-verify.md` (written before the exclusion results existed), aggregator section 2.13 of
`wave2-status.md` and `wave2-aggregator-hessrelax6.md`, dossier, saved thread (one post by `fqlx`).

## Replays performed by the referee (own code unless stated; receipts in the packet, prefix `referee_`)

| replay | result |
|---|---|
| 6 exact degree-5 certificates vs. my own E2 columns (21+540+126) | 6/6 verified (`referee_replay.py d5`, 12 s) |
| 58 exact degree-6 certificates vs. the TRUE dense E2 columns and all 2376 Y-odd generators of `J_12` | 58/58 verified (`referee_replay.py d6`, 36 s) |
| char-0 Groebner basis regenerated from `hpipe_d6.C.json` (Singular 4.3.2) | affine dim 1, 0 of 35 coordinates outside the radical (`referee_d6_charQ.out`) |
| E1-kernel dims at the unused prime 1000000009, all parities, d = 4..7 | 0, 6, 36, 129; J block at regular-sequence rank in every degree (`referee_replay_e1.json`) |
| 3 exact `K_7` generators: `N(F) = A D` over Q; independence mod `J_7` at 1000000009 | all exact; rank `[J_7 | D1 D2 D3] = 129 = 126 + 3` (`referee_replay_k7.json`) |
| membership table at 1000000009 with regular-sequence rank check | all listed non-memberships rigorous over Q (`referee_membership.json`) |
| lane code `fam_run.py 5` at 1000000009 | rank 507, 21/21 monomials `a_k a_l` in the row space |
| aggregator `hessrelax6` input re-run with msolve 0.6.5 (`-g 1`, p = 10^9+7) | 707 leading monomials, pure power of all 35 variables, max exponent 3 |

## Claim-by-claim verdicts

Legend: VERIFIED = exact over Q, or modular with a valid transfer, and replayed/derived by the referee;
PLAUSIBLE = modular evidence only (or not replayed); UNSUPPORTED = not established by any receipt.

| # | claim | source | verdict | evidence / note |
|---|---|---|---|---|
| 1 | `F = f + g + (5/7) L^2 M`, `M = L(Y)`; `B_{eps_X F}(A X) = 8 t A L^2 M`; `B_F(A eps/6) = A F`; `f - 2g = (E_X F - 2 E_Y F)/6 in J` | excl R1, dossier 1.1/1.6/4.5 | VERIFIED | exact identities; re-asserted in every referee run (`w2lib` import assertion, `referee_membership.py`); verify S1–S9 agree |
| 2 | Euler gauge: WLOG `N(F) = A D` exactly, identity unchanged | excl Method 1 | VERIFIED | hand derivation; uses `F` irreducible (smooth) and `D` coprime; the two `F`-terms cancel since `div(phi E) = (deg phi + 6) phi` |
| 3 | Koszul parametrisation of `N - N_0`; weight-2 structure | excl Method 2 | VERIFIED | first syzygies of a regular sequence are Koszul; `F` smooth (`F == sum Z_i^6 mod 5`) |
| 4 | `V_d = 0` for `d <= 4`; `V_5 = J_5`; `V_6 = J_6` over Q, all parities (mixed included) | excl Method 5, verify_e1kernel, dossier 4.2 | VERIFIED | modular upper bound with J block at regular-sequence rank: replayed at 10^9+9 (0, 6, 36); receipts at 10^9+7 (two code paths), 2^31-1, 10007 |
| 5 | No certificate with `deg D = 5` (any parity) | excl R2 | VERIFIED | six certificates `lambda^T E2 = a_k^2` replayed against independently built columns; mod-p rank 507 at 10^9+9 |
| 6 | 58 exact certificates in the J-relaxed weight-2 ideal at `d = 6` (25 squares, 33 forms) | hpipe_d6 step C | VERIFIED | each `lambda` kills `J_12^-` and reproduces the form on the true E2 columns (sparse representatives not needed) |
| 7 | `V(I) = {D = cF}` at `d = 6`: no certificate with `deg D = 6` (any parity) | hpipe_d6 step D | VERIFIED | Singular char 0 regenerated from the receipt: all 35 coordinates of `V_6/<F>` in the radical |
| 8 | sparse mod-J column identity (*) `P_kl + P_lk == -A[Hess(G_k,G_l)+Hess(G_l,G_k)]` | w2lib | VERIFIED | hand derivation (expand `N(D)` for `N = A G`, `D = G(F)`); immaterial to the referee replay (true columns used) |
| 9 | mixed-parity degree-6 Hessian relaxation has m-primary leading ideal mod 10^9+7 (419 quadrics, 35 vars) | aggregator hessrelax6 | PLAUSIBLE (modular) | replayed with msolve; the transfer to Q is valid in principle (full-rank Macaulay matrix in degree 71; standard monomials of `(S/J)_6` mod p form a p-integral Q-basis), so this is an independent second route, but the proof of record is 5–7 |
| 10 | `V_7 = J_7 (+) K_7` over Q with `dim K_7 = 3` (Y-even, base-type, exact generators with exact E1 lifts); `V_7^- = J_7^-` | k7exact, verify_e1kernel | VERIFIED | exact lifts and independence mod `J_7` replayed at 10^9+9; upper bound 129 at 10^9+9 and 10^9+7 |
| 11 | degree 7 excluded or "screened" | (none; draft template) | UNSUPPORTED | `hpipe_d7`: 0 of 123 coordinate squares in the row space mod p; Singular mod p timed out (3600 s). Degree 7 is open |
| 12 | `dim V_d/J_d = 9+9 (d=8), 27+27 (d=9)`; `K_9 = K_7 S_2 + K_8 S_1 mod J` | wave 1, verify_e1kernel | PLAUSIBLE (modular) | 10^9+7 (d = 8, 9, both parities, rankJ_ok; generation check true) and 10007; no exact generators computed beyond degree 7 |
| 13 | `(J:A)` minimal generators 6/3/3/456 in degrees 5/7/8/10; HF of `S/(J:A)` | quotJA_p65521, dossier 4.6 | PLAUSIBLE (modular) | two primes agree (65521, 10007); the 32003 run never executed — do not cite it |
| 14 | `L, M, f, g, L^2 M, f-g, A, A L, A M, A L^2, A M^2, A L M notin J`; `fg, f L^2M, g L^2M, (L^2M)^2, f^2, g^2 notin J` | dossier 4.2/4.6, excl ledger | VERIFIED | non-membership mod p with J at regular-sequence rank is rigorous; replayed at 10^9+9 |
| 15 | `Delta' = 4 t^3 L^6 M^3 + 27 f^2 g in J` | dossier 4.6, quotJA | PLAUSIBLE (modular) | membership mod 65521 and 10007 transfers nothing to Q |
| 16 | `[A]: (S/J)_k -> (S/J)_{k+6}` surjective for `k >= 12` | dossier 4.6 | VERIFIED (logic) | modular surjectivity + `rank_p <= rank_Q` with J at regular-sequence rank; receipts at two primes |
| 17 | wave-1 fixed-denominator inconsistencies (`f, g, eps_X F, X0 F_0, ...`, all degree 6; single prime, engine v1/v2) | wave-1 log, aggregator | PLAUSIBLE as stated | now implied by 5–7; cite only the degree bound |
| 18 | attack1 §1.2 "rigorous: no degree-6 family `D_t` with `D_0 = f`" (argument) | attack1 | UNSUPPORTED (argument refuted, verify L5) | the *conclusion* for `t = 5/7` now follows from 7 by a different argument |
| 19 | "base-type ⇒ the compact-Y obstruction predicts E2 failure on the degree-7 family" | dossier 4.2 | UNSUPPORTED | verify L6: no prediction for sums base-type + polar; the degree-7 E2 failures are modular samples only |

Counts: **VERIFIED 11** (1–8, 10, 14, 16), **PLAUSIBLE 5** (9, 12, 13, 15, 17), **UNSUPPORTED 3** (11, 18, 19).

Strongest verified statement: *for the exact target (`F`, `A = X0^3 Y0^3`, `q = 1`) there is no certificate with
`deg D <= 6` — for every nonzero `D` coprime to `F`, of any Y-parity (mixed included), with `Q` arbitrary — proved over
Q by exact rational certificates and a characteristic-zero Groebner basis; the first open degree is 7, where the
E1-admissible space is exactly `V_7 = J_7 (+) K_7` (129-dimensional, `dim K_7 = 3`).*

## Transfer arguments checked

* mod p -> char 0 for E1 kernels: `nullity_p >= nullity_Q`, syzygy part identical when `rank_p J = rank_Q J`
  (regular-sequence Hilbert function checked in every run); direction correct (modular = upper bound).
* exact certificates: Farkas-type duals `lambda` verified over Q column by column (including the `J_12` columns at
  `d = 6`, which is what makes the J-relaxed certificate a certificate for the full E2 system).
* m-primary leading ideal -> no projective zero: `LT(I_p) ∋ x_i^{e_i}` for all `i` gives `(I_p)_N = (S_p)_N`,
  `N = 1 + sum(e_i - 1) <= 71`; full rank of the degree-`N` Macaulay matrix persists over Q provided the mod-p quadrics
  are reductions of the Q-quadrics — true for coordinates in the mod-p standard-monomial basis of `(S/J)_6`
  (saturated-lattice argument), given the Hilbert check.  Valid; used only as a cross-check here.
* sum principle at `q = 1` (linearity of served targets in `(A, v)`): not used in the exclusion chain; its wave-1
  use (target_image) is unaffected (verify L9).
* prescribed-support powers / `D^2` presentations: not applicable (`q = 1` is fixed; dossier gap 8).
* non-membership transfer (row 14) and modular surjectivity (row 16): correct directions.

## Reconciliation with the adversarial verify lens

No substantive disagreement.  The verify lens (22 confirmed / 1 refuted / 1 weakened) predates the exclusion results;
its refutation L5 and weakening L6 are adopted (rows 18, 19), and its L4 caution about single-prime "inconsistent"
verdicts is adopted (row 17).  Its L2 bound ("no denominator of degree <= 4 of any parity") is superseded by rows
5–7 (`<= 6`).  The aggregator's wording "Best rigorous bound: deg D >= 7 for arbitrary D" is confirmed, now by the
exact-certificate route rather than only by the modular m-primary route.

## Wording required in the post

`DMAX_PROVED = 6`; there is **no** screened degree beyond 6 (degree 7 is not even excluded mod p); the prime 32003 was
not used; the exclusion report's "exact step still running" for degree 6 is superseded by `hpipe_d6` step D; cite the
aggregator's m-primary result only as independent agreement; never claim nonexistence of a certificate.

Packet snapshot (pinned): https://github.com/kejace/sextics/tree/a6a4f461b3245f146a6eec3ed1387195b03d4f05/research/results/rank-one-biquartic/wave2-excl
