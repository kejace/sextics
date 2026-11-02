# Referee report — lane `one-reflection-skew`, wave-2 exclusion packet

Referee: Claude (Anthropic) for `kejace`, 2026-11-01/02. No human review; nothing Lean-verified; no certificate claimed.
Inputs: draft `wave2-one-reflection-skew-draftpost.md`, report `wave2-one-reflection-skew-excl.md`, packet
`results/one-reflection-skew/wave2-excl/` + `wave2-excl.log.jsonl`, parallel verify ledger `wave2-one-reflection-skew-verify.md`
(+ `wave2-verify.log.jsonl` on the remote mirror), aggregator §2.9, dossier, saved thread #127, search report, and the remote
receipts pulled with `rjob pull` (`remote-extradirs.log.jsonl`, `remote-families-f3c.log.jsonl`, `d7_sub_*`, `d7_5v_*`).
Referee code/logs: `results/one-reflection-skew/wave2-excl-referee/` (python-flint 0.9 only, no lane/engine imports;
Singular 4.3.2 and msolve for Gröbner steps). Pinned snapshot: commit `cdc20275d7641a9929d8da145a66bd94c5f501b2` of
`kejace/sextics` (branch `research`); packet README-referee.md is the public claim table.

Verdict scale: VERIFIED = logic checked and computation replayed (exactly, or modular with a one-sided transfer I checked);
PLAUSIBLE = consistent modular/screen evidence or receipt only; UNSUPPORTED = not established by the recorded receipts.

## Counts

15 claims examined: **9 VERIFIED, 4 PLAUSIBLE, 2 UNSUPPORTED** (one of them a wording error now withdrawn, one pending work).

## Claim table (condensed; full table in README-referee.md)

| # | Claim | Referee replay | Verdict |
|---|---|---|---|
| 1 | E1 kernels: 0 (d≤4), J_5 (6), J_6 (36) | own exact fmpz ranks over Z, 24 sectors | VERIFIED |
| 2 | C1: Hessian classes mod J_10 rank 15/15, six zero classes F_02, F_i5 ⇒ no degree-5 D | own code p=1000000007 with good-prime check (rank J_10 = 1497 = char-0 Hilbert value) | VERIFIED |
| 3 | C2: degree 6 excluded (28+2 squares, 5 survivors, 9 quadrics, Macaulay 35/35) | lane pipeline `job_degree.py 6` replayed end-to-end (45 s), counts identical; final one-sided step checked | VERIFIED |
| 4 | Exact E1 table d=5..14; 5-var (J_f:A)_d | own exact to d=8 (135/126/9, 480/336/144), mod p at 9 (1144/756/388); 5-var 0,0,5,25,84 exact, 310, 594 mod p | VERIFIED (d≤9) |
| 5 | Degree-7 form 84+25+5+15 = 129 | bookkeeping from exact dims | VERIFIED |
| 6 | Degree-7 Y2-odd mixtures excluded | exact Macaulay rank 35/35 over Z on the 21 relaxed quadrics of the 5 odd survivors (remote `d7_sub_odd.pkl`) | VERIFIED (final step) |
| 7 | Y2-free lemma (strict 5-var form) | hand proof both directions; agrees with verify lens §A.3 | VERIFIED |
| 8 | 9 non-Jacobian d=7 directions EXCLUDED-exact | own system built from `d7_fam.pkl`, exact ranks 819/820 ×9 | VERIFIED |
| 9 | Axes Y1Y2^6, Y0Y2^6, X1Y2^6 and Y1f, Y0f, X1f excluded exactly | own exact ranks 765/766, 794/795, 765/766; 561/562, 579/580, 561/562; lane items #7/#40/#106 are exactly these monomials | VERIFIED |
| 10 | Pure even classes 0010/0100/1000: "relaxed ideal supported on one coordinate axis; axis excluded ⇒ handled" | msolve at p=1000000007 (inputs rebuilt from exact quadrics) and Singular: same leading ideals; but dim 1 ≠ "axis": in 5-var systems the generators restrict to −68036492·l36² ≠ 0 on the axis; char-0 Singular for class 0100: dim 1, degree 67 (both relaxations) | UNSUPPORTED as stated; replaced by "finitely many projective points (67 for 0100 over Q; ≤495/419 mod p for 0010/1000), unidentified" |
| 11 | 129-dim family stage B/C ranks 748/757, no coordinate killed | not replayed | PLAUSIBLE (informational) |
| 12 | Structured products d=7..12: 1392 candidates, 1268 E1-fail, 93 excluded exactly, 31 pending | receipts read; exact duals not replayed | PLAUSIBLE (receipt-level) |
| 13 | d=8 non-Jacobian 5-var directions: 137 excluded exactly; "extra = 228" | receipts read; exact count of non-Jacobian directions is 310 − 175 = 135 (label error) | PLAUSIBLE; fix label |
| 14 | Degree-7 "even"/"full" leaves, class 1110 | job still running, `.lead` empty | UNSUPPORTED (pending) |
| 15 | Wave-1 single-prime 5-var screens | logic only | PLAUSIBLE (evidence, as the report says) |

## Transfer arguments

Checked and sound: Euler normalization for q = 1 (needs gcd(D,F) = 1 and F irreducible); F-gauge D ~ D + FL at fixed degree
(same N works; K_A changes by F² div w); sector restriction by group averaging; necessity of the relaxed quadrics (Koszul
syzygies exhaust ker N ↦ N(F) for a regular sequence); "full column rank mod p of exact integer quadrics ⇒ V_C = {0}";
"HF_Q ≤ HF_p ⇒ finitely many projective points over Q"; exact dual certificates (asserted over Q in `w2x.exact_dual_certificate`).
Not sound and now withdrawn: "pure powers of all coordinates but one in the leading ideal ⇒ support is the remaining axis"
(claim 10). The 6→5 reduction of the search lens (6/5-equation; lift (D², N_iD, Y2B/5, Q=0) of degree 2d) is correct in both
directions; it preserves the denominator only for Y2-free D (then it coincides with the lemma of claim 7).

## Reconciliation with the parallel verify ledger

Its 20 table rows: I agree with all CONFIRMED rows I touched (1–7, 9–17, 19, 20); its REFUTED-as-stated row 8 (dossier's list
of two-element sign-change symmetries; conclusion unaffected) and its downgrade of row 18 (single-prime screens) are correct.
Its JSON log also records a "PROBLEM" flag on the 5→6 lifting direction for a *served* target D = f_X0 at d = 5; that concerns the
search lens's numerical pipeline, not this packet, and the algebraic lift itself is correct (checked above). Its hand proof A.3
and my proof of claim 7 agree.

## Referee's own corrections

My first "second-prime" msolve attempt edited the characteristic line of the lane's `.ms` files; that is invalid (coefficients
are stored as residues of the original prime) and was discarded. All reported second-prime runs regenerate inputs from the exact
integer quadrics in the lane's pickles. A Singular mod-p multiplicity loop for the 19-coordinate classes was stopped (not needed);
the characteristic-0 Singular runs for classes 0010 and 1000 were still running at sync time.

## Strongest verified statement

For the immutable target, every certificate denominator has degree ≥ 7 (exact; `rbajaj5`'s C1/C2, re-verified by two further
independent implementations), and at degree 7 every Y2-odd denominator, mixtures included, is excluded (final step exact over Z);
moreover, for Y2-free D the certificate problem is equivalent to the strict five-variable system N'(f) = AD, D div N' − N'(D) = 0
with the same D (proved).

## Recommendations for the public comment (applied in comment.md)

Lead with the re-verified bound and the Y2-odd result; state the lemma; list single-denominator exclusions as such; withdraw
the "axis" reading and give the finite-point statement; keep classes 0010/0100/1000 attributed to `fqlx`'s addendum as the only
claimed closure; list open items and the pending jobs; fix the d = 8 label; no nonexistence claims beyond the stated degrees.
