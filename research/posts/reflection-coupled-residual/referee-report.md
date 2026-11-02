# Referee report — `reflection-coupled-residual` (#154), wave-2 exclusion packet

Referee: Claude (Anthropic) for `kejace`; no human review; nothing Lean-verified; no certificate claimed by anyone.
Inputs: draft post, exclusion report (R1–R11), packet `wave2-excl/` (receipts, scripts), adversarial verify notes (23 rows),
aggregator status, dossier, saved thread. Independent replays: `wave2-excl/referee/` (own code, python-flint 0.9.0 at
p = 1000000021; Singular 4.3.2 at p = 1000003; msolve 0.6.5; exact `fractions`/`fmpq_mat` where stated). README-referee.md in the packet
gives the replay commands and the transfer arguments.

Verdict key: **VERIFIED** = replayed by the referee with independent code (or exact over Q) and transfer argument checked;
**PLAUSIBLE** = lane receipt only (single engine / single prime, or screen), transfer argument valid where applicable;
**UNSUPPORTED** = wording not backed by the receipts.

| # | claim (lane wording) | lane evidence | referee replay | verdict | downgrade / note |
|---|---|---|---|---|---|
| 1 | E1 kernel mod J is 0 for d <= 7 in all four Klein characters, 1 at d = 8 in (0,0), 0 in the others (S1, S2, S7) | e1char.jsonl, p = 2147483587 (+2147483579), good reduction every row; verify lens at 10^9+7 | ref_e1.jsonl at 10^9+21: identical table, good reduction every row | VERIFIED | zeros are theorems over C; the d = 8 value 1 is an upper bound, made exact by row 3 |
| 2 | d = 9: non-Jacobian kernels 21/0/21/0; d = 10: 186/126/136/108 (S10) | e1char.jsonl, e1char.r2.jsonl | ref_e1.jsonl: 21/0/21/0 at d = 9; d = 10 (0,0) = 186, others see receipts | VERIFIED (d = 9), VERIFIED/PLAUSIBLE (d = 10, upper bounds) | nonzero modular kernels are only upper bounds for the rational ones |
| 3 | K8 is a rational E1-admissible class outside J: A*K8 = sum N_i F_i exactly over Q; K8 not in J (S6) | K8_E1_exact.json (CRT lift), two primes | ref_k8.py: identity exact over Q with fractions; K8 not in J_8 mod 10^9+21 (good reduction) | VERIFIED | closes the verify lens's pending row 11 |
| 4 | Degree 5: the 20 classes [A F_ij] are independent in (S/J)_10, hence no degree-5 denominator (R1, S3) | hess.jsonl phi-kernel/phi-verdict | ref_hess.py rank20: F_{X2 Y2} = 0, rank 20 mod J_10, good reduction | VERIFIED | logic checked: c_i c_j = 0 for all 20 pairs including squares forces c = 0 |
| 5 | Degree 6: relaxation on W = V/T (35 unknowns, all characters) is zero-dimensional, HF 1,35,285,477,225,36,1,0 => no degree-6 denominator, mixed characters included (R2, S4) | gbcheck/hess_d6 files, Singular minAssGTZ | ref_hess export (W 35, 345 independent quadrics) + msolve (GB 709, quotient 1060) + Singular (dim 0, vdim 1060, same HF) | VERIFIED | transfer HF_Q <= HF_p checked (needs good reduction in degree 12: recorded) |
| 6 | Degree 7, characters (1,1), (1,0), (0,1): zero-dimensional (R3, S5) | gbcheck.jsonl | referee: 85 / 239 / 1109 quotient dims, HF 1,18,54,12,0 / 1,22,121,83,12,0 / 1,38,460,556,54,0 | VERIFIED | — |
| 7 | Degree 7, character (0,0): zero-dimensional, quotient 3028, HF 1,42,629,1986,359,11,0 (new since the draft; gbcheck.jsonl line 9, 4224 s) | one msolve+Singular run at 2147483587 | referee run launched (local 4 threads, and queued on the remote box); see `referee/ref_gb.jsonl` | PENDING-AT-WRITING: see addendum below | transfer argument identical to row 5; the draft post must be updated: the block is finished |
| 8 | Degree 7, pair (1,0)+(1,1): zero-dimensional, quotient 4399 (S12) | gbcheck.jsonl | referee: W 40, 300 independent quadrics, GB 2712, quotient 4399, HF 1,40,520,1543,1695,567,33,0 (64 s msolve at 10^9+21) | VERIFIED | also implied by row 11 + row 5 |
| 9 | K8 is not a certificate denominator: restricted system inconsistent, rank_Q = rank_p certified by 482 a-priori kernel vectors of rank 476 (R4, S6) | rigorfix.jsonl (two primes), full system | ref_fixedD.py (own construction from the contract identity): 22016 x 4256, a-priori rank 476, rank_p(M_sub) = 3780, rank_p([M|b]_sub) = 3781 | VERIFIED | a-priori vectors checked symbolically and mod p; the row-subset argument is spelled out in README-referee.md |
| 10 | All six coordinate sections X cap {z_i = 0} are smooth (R8, S8) | hypsmooth.jsonl (p = 32003); aggregator Singular vdim 3125 | Singular at 1000003: dim 0, vdim 3125 for all six; X itself vdim 15625 | VERIFIED | — |
| 11 | Reduction lemma served(z_i G) = served(G); Y2-parity reduction (R8, R9, S9) | "sum principle, book ch. 3 §1.1" | logic re-derived: needs Chapter 1 prescribed-support theorem (1 <=> 3, denominator D_red^q, q = 1) + Chapter 3 §3.1 Gysin description; H^2(Z_i) = Q h and h^2 vanishes off any hypersurface section | VERIFIED conditional on the book's Ch. 1 theorem (Hodge-module inputs, not Lean-verified) | cite Ch. 1, not only "ch. 3 §1.1"; state the conditionality in the post |
| 12 | Degree-9 (0,0): 7/7 tested kernel members certified non-denominators (S13) | rigorfix.r2.jsonl (rank 6006, nullity 798, two primes) | one representative (basis0) replayed with ref_fixedD.py: see addendum | VERIFIED (1 of 7 replayed) / PLAUSIBLE (other 6, same protocol) | screens of single elements, not of the 21-dimensional family |
| 13 | "No structured shape contains a member of the K8 family" (R7, S7, draft (b): "nothing outside J for d <= 8") | structscreen.jsonl | exact over Q (fmpq_mat): the family K8 + J_8 **does** contain a member with all Y-exponents even, and the Klein-invariant subring contains K8 itself; the lane's own d = 8 receipts for those two shapes say mod J = 1 | UNSUPPORTED as worded | reword: no cone with a line/plane vertex, no polynomial in squares, no coordinate multiple, no X- or Y-degree <= 4 member; the even-Y-exponent family (294 monomials, X arbitrary) and the invariant subring do meet it. No conflict with wave 1 / verify row 15: wave 1's "only even exponents" is the 126-monomial polynomial-in-squares shape, which indeed has no member (exact, referee) |
| 14 | Structured E1 screens at d = 9..12, cone families at d = 11, 12 (R7, R11, S11, S14) | structscreen.jsonl, conescreen.jsonl | not replayed | PLAUSIBLE (screens; E1 zeros are rigorous by the row-1 argument) | keep as screened |
| 15 | K8 family: ob(K8) != 0, linearised system rank 156 with ob(K8) not in the image => no first-order correction (R6) | wave-1 e2lin receipt | not replayed | PLAUSIBLE (screen; the quadratic system on 103 unknowns timed out) | keep as screened |
| 16 | Fermat positive control: pipeline returns a 19-dimensional solution cone at degree 6 (R10) | control-fermat/ | not replayed | PLAUSIBLE | control, no claim about the target |
| 17 | Draft post quotes the record SHA256 `114c5b4f…38a2` | copied from the October thread | sha256sum of the current `instance.rec` is `e694342f…bae7b6` | UNSUPPORTED (stale hash) | quote the current hash or none; F, A, q were parsed from the current record |

## Transfer arguments checked (task 2)
* mod p -> characteristic 0 for E1 kernels and for Hilbert functions: rank_p <= rank_Q plus good reduction (Fermat count) in the degrees used; the saturated lattice J_Q cap Z_(p)[z] reduces onto J_p. Correct.
* m-primary leading ideal (dim 0, finite quotient) -> no projective zero -> V(ob) = T over C: correct (HF_p(k0) = 0 => HF_Q(k0) = 0).
* Fixed-D "Farkas"/rank certification: correct as described in row 9; verified that both families of a-priori vectors are exact kernel vectors of the contract identity.
* Sum principle at q = 1: the algebraic direction served(G) subset served(z_i G) is trivial; the converse needs the prescribed-support theorem (Ch. 1, q = 1, denominator D_red) and the Gysin description (Ch. 3 §3.1); the (2,2)-classes from the smooth section (h^2 and classes of surfaces inside Z_i cap Z_G) vanish off Z_G. Correct modulo those theorems.
* Prescribed-support powers: at q = 1 existence over D is equivalent to existence over D^m (Ch. 1), so single-D exclusions (K8, degree-9 representatives) also exclude all powers of those denominators.

## Disagreements with the verify lens
None on shared rows (rows 1, 3, 4, 10 agree). The verify lens left K8's exact E1 membership pending; it is now exact (row 3). The verify lens did not examine R8/R9/R11 or the structured-shape claim; the referee downgrades row 13 and asks that the reduction lemma cite Chapter 1.

## Addendum (filled after the pending runs)
* Row 12, degree-9 representative `e1ker-d9-chi00-basis0` (referee, p = 10^9+21): system 32214 x 6804, 816 a-priori kernel vectors (all exact kernel vectors), rank 798; row subset 8400: rank_p(M_sub) = 6006 = 6804 - 798, rank_p([M|b]_sub) = 6007 -> certified inconsistent over Q: True (441.6 s). Verdict row 12: VERIFIED for this element (agrees with the lane's rank 6006 / nullity 798).
* Row 7, degree-7 character (0,0): the referee replay at 10^9+21 (local, 4 threads) and the queued remote job had not finished when this report was written; the lane's single run (quotient 3028, HF 1,42,629,1986,359,11,0, 4224 s at 2147483587) has the same valid transfer argument as rows 5, 6, 8. Verdict row 7: PLAUSIBLE (single run, not independently replayed yet).
* Degree-10 E1 kernels replayed: 186/126/136/108 with good reduction (row 2 upper bounds confirmed at a third prime).

**Final counts (17 claims):** VERIFIED 11 (row 11 conditional on the book's Chapter 1 theorem), PLAUSIBLE 4, UNSUPPORTED 2.

**Strongest verified statement.** For the protected (F, q = 1, A = X0^3 Y0^3): every certificate denominator has degree >= 7 with no parity or Klein-character restriction (E1 by character at three primes with good reduction; Hessian relaxation zero-dimensional at degrees 5, 6 with the transfer HF_Q <= HF_p); no degree-7 denominator of a single Klein character exists for (1,1), (1,0), (0,1) and none with components only in (1,0)+(1,1) [the (0,0) block: lane's run, see row 7]; the unique non-Jacobian degree-8 E1 class K8 is certified not a denominator over C.
