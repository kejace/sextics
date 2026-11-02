# Referee pass over the `wave2-excl` packet (lane `one-reflection-skew`, q = 1, A = X0^3 Y0^3)

Referee: Claude (Anthropic) for `kejace`, 2026-11-01/02. No human review. Nothing here is Lean-verified; no certificate is claimed.
Inputs: `wave2-one-reflection-skew-excl.md` (report), `wave2-one-reflection-skew-draftpost.md` (draft), this packet and
`../wave2-excl.log.jsonl`, the parallel verify ledger `wave2-one-reflection-skew-verify.md` (+ its JSON log on the remote mirror),
the remote receipts pulled to `/tmp/hodge/remote/one-reflection-skew/wave2-excl/` (`remote-extradirs.log.jsonl`,
`remote-families-f3c.log.jsonl`, `d7_sub_*.pkl/.lead`, `d7_5v_*.pkl/.lead`), the dossier, the saved thread (#127) and the search report.
Referee code and logs: `../wave2-excl-referee/` (`reflib.py`, `fixedD.py`, `e1dims_more.py`, `msolve_p2.py`, `referee.log.jsonl`,
`lane-replay-d6.*`, `*_p2.ms/.lead`, `*_char0.sing`). The referee code imports neither the lane library `w2x.py` nor the engine;
it only reads the lane's pickled exact data where stated. Fresh primes: 1000000007, 998244353 (lane: 2147483647 and 1073741827).

Verdict scale. VERIFIED = logic checked and the computation replayed (exactly over Z/Q, or modular with a one-sided transfer
argument that I checked). PLAUSIBLE = consistent modular/screen evidence, no rigorous replay. UNSUPPORTED = not established by
the recorded receipts (wording must be withdrawn or weakened).

## 1. Claim table

| # | Claim (report/draft) | Lane receipt | Referee replay | Verdict |
|---|---|---|---|---|
| 1 | deg D <= 4: E1 kernel {D : A D in J} = 0; d = 5: kernel = J_5 (dim 6); d = 6: = J_6 (dim 36) | `hessian-modJ-d5`, stage A d6 (exact) | own code, exact fmpz ranks over Z, all 24 sectors: 0, 6 = J_5, 36 = J_6 | VERIFIED |
| 2 | C1 (`rbajaj5`): 21 Hessian classes [A F_ij] mod J_10, six identically zero (F_02, F_i5), the 15 others independent mod J ⇒ a = 0, no degree-5 denominator | `hessian-modJ-d5` (exact) | own code at p = 1000000007: rank J_10 = 1497 = char-0 value (good prime), rank [J_10 | 15 classes] = 1512 ⇒ 15 independent over Q (one-sided); zero classes exactly {(0,2),(0,5),(1,5),(2,5),(3,5),(4,5)} | VERIFIED |
| 3 | C2 (`rbajaj5`): degree 6 excluded — 35-param gauge-fixed family, 28+2 certified squares, 5 diagonal survivors, 9 relaxed quadrics, Macaulay deg 3 rank 35/35 | `relaxed-modF-d6` | lane pipeline replayed end-to-end locally (`job_degree.py 6`, 45 s, log `lane-replay-d6.log.jsonl`): identical counts; final step is a full-rank Macaulay matrix of exact integer quadrics mod p ⇒ V_C = {0} (one-sided, sound) | VERIFIED (pipeline replay; final step's logic checked) |
| 4 | E1-kernel table over Q, d = 5..14 (135 = 126+9, 480 = 336+144, 1144 = 756+388, ...) | `e1-kernel-exact`, `e1dims_d*.json` | own exact ranks: d = 7: 135/126/9 (non-Jacobian in fine sectors (0,0,1):1, (3,0,4):7, (3,3,1):1), d = 8: 480/336/144; d = 9 mod p: 1144/756/388; 5-var (J_f:A)_d = 0,0,5,25,84 (d = 3..7) exact, 310, 594 (d = 8, 9) mod p, Jacobian 0,0,5,25,75,175,350 | VERIFIED for d <= 8 (exact), d = 9 (mod p, consistent); d >= 10 not replayed (verify lens agrees via Singular at two primes) |
| 5 | Degree-7 E1 form 84 + 25 + 5 + 15 = 129 by t-degree (`fqlx`'s addendum form) | `d7-E1-form-check` | follows from (J_f:A)_d = 84, 25, 5, 0, 0 for d = 7..3 plus the t^5 S'_2, t^6 S'_1, t^7 blocks (15+5+1) and the six F-gauges: 135 - 6 = 129 | VERIFIED (bookkeeping checked on my exact dims) |
| 6 | Degree 7, all Y2-odd mixtures excluded (confirms `fqlx` addendum) | remote `d7-subfamily-odd`, `d7_sub_odd.pkl` | the 21 relaxed quadrics on the 5 odd survivors (18 distinct; quadric rank 9) have a degree-3 Macaulay matrix of exact rank 35/35 **over Z** (own code) ⇒ V_C = {0}; the square-killing steps that reduce 40 -> 5 survivors are lane exact-dual computations, not replayed | VERIFIED (final step exact; two independent implementations agree) |
| 7 | Y2-free lemma: for Y2-free D a certificate exists iff a Y2-free N' of degree d+1 has N'(f) = A D and D div N' - N'(D) = 0 | report (proof sketch) | hand-checked both directions (see §2); independently also derived by the verify lens (§A.3 there) | VERIFIED (proof) |
| 8 | The 9 non-Jacobian degree-7 kernel directions (lane basis) each fail the full fixed-D test ("EXCLUDED-exact") | 9 `fixed-D` rows, rank 819/820 | own sector-restricted system (714 N-columns, 123 Q-columns, 1548 rows) built independently from `d7_fam.pkl` vectors; exact ranks over Z 819/820 for all nine ⇒ inconsistent over Q | VERIFIED |
| 9 | The "axis" denominators Y1 Y2^6, Y0 Y2^6, X1 Y2^6 (= lane items #7, #40, #106 up to scalar, checked) and Y1 f, Y0 f, X1 f (strict 5-var) are excluded exactly | `fixed-D` 765/766, 794/795, 765/766; `fixed-D-5var` 561/562 | own code, exact over Z: 765/766, 794/795, 765/766; 5-var 561/562, 579/580, 561/562 | VERIFIED |
| 10 | Degree-7 pure even classes 0010/0100/1000: "the relaxed ideal is projectively 0-dimensional with support = one coordinate axis (Y1 Y2^6, Y0 Y2^6, X1 Y2^6), multiplicity 495/67/419; the axis is excluded exactly" (report; also "the Y2-free reduction does the job on the t-free part") | `d7-pure-even-class-leaf-structure`, `5v-d7-tor:*` ("OPEN") | msolve re-run at p = 1000000007 on inputs regenerated from the exact integer quadrics: same leading-ideal sizes (964/195/771) and the same single coordinate without a pure power, in both the 6-var and the strict 5-var systems. **But a pure power in the leading ideal does not put that power in the ideal.** Correct reading: Krull dimension 1 (Singular: `dim 1`, `mult 67` for class 0100), i.e. projectively finitely many points — not "support = the axis". In the strict 5-var systems the generators restricted to the axis are a nonzero multiple of l_k^2 (Singular `subst`: -68036492 l36^2 for class 0100), so the axis is not even a solution there; in the 6-var systems the axis is a genuine common zero (all l_k^2 coefficients are exactly 0) but nothing identifies the remaining points. | UNSUPPORTED as stated. Rigorous replacement (one-sided, HF_Q <= HF_p): in each of these classes the set of D satisfying the necessary relaxed quadrics is a finite set of projective points, of total length <= 67 (class 0100, replayed) and <= 495 / 419 (classes 0010 / 1000, lane values, not replayed). The classes are **not closed** by the lane's computation. |
| 11 | Degree-7 full 129-dim family: stage B 1019 quadrics rank 748, stage C 1094 quadrics rank 757, no coordinate killed | `stageB/C-modJ-d7`, `d7_stageB/C.json` | not replayed (structure only: consistent with 56/57 off-diagonal monomials in span) | PLAUSIBLE (informational; no exclusion claimed) |
| 12 | Structured products of invariants, degrees 7–12, <= 3 factors: 1392 candidates, 1268 fail E1 rigorously, 124 pass E1, 93 excluded exactly (fixed-D), rest pending | remote `families_d7-12_f3.jsonl`, `remote-families-f3c.log.jsonl` (job still running) | counts read from the receipts; E1 failures use the exact Jacobian rank (sound); exact duals not replayed | PLAUSIBLE→VERIFIED-by-receipt for the 93+1268 rows; 31 E1-passing candidates still pending |
| 13 | Degree-8 non-Jacobian 5-var directions: 137 excluded exactly (strict 5-var fixed-D) | remote `extra_dirs_d8.jsonl` (job running) | receipts read; the family label "n = 310, extra = 228" disagrees with the exact count 310 - 175 = 135 non-Jacobian directions (my exact dims; verify lens agrees) — a labeling issue in the lane's family kinds, not a mathematical claim | PLAUSIBLE (single-D statements; label to be fixed) |
| 14 | Degree-7 cases "even" and "full" of the sub-family job; class 1110 (5-var and 6-var); stage D leaf for the 129-dim family | `.lead` files empty / job running | nothing to replay | UNSUPPORTED (pending) — report as open |
| 15 | Wave-1 5-variable fixed-D screens (318 experiments, single prime, no dual) | wave-1 `log.jsonl` | logic checked (verify lens §A.4 and §2 below) | PLAUSIBLE only ("evidence, not proofs", as the report already says) |

## 2. Transfer arguments checked

* **Cleared identity and normalization.** With E = A D^2 - q D N(F) + F (D div N - N(D)) = F^{q+1} Q and q = 1: mod F,
  D (A D - N(F)) ≡ 0; gcd(D,F) = 1 and F irreducible (smooth) give A D - N(F) = g F; N' = N + (g/6) Z (Z Euler) satisfies
  N'(F) = A D and E(N') = E(N). So E1: N(F) = A D, E2: B_D(N) = F Q is lossless for q = 1. Checked.
* **F-gauge at fixed degree.** If (D_0 + F L, N, Q) is a certificate then (D_0, N, Q') is one with the same N: the fields
  N/(D_0 + F L) and N/D_0 differ by F·w with w regular away from F, and K_A changes by F^2 div w. Checked (this justifies
  quotienting the E1 kernel by F S_{d-6}, as the thread and the lane do).
* **Sector restriction** of the fixed-D system to the isotypic components of the 24-element diagonal group (D a diagonal
  eigenvector; F invariant, A of character (3,3,0) in the lane's weights) is lossless by averaging. Checked.
* **Necessary quadrics.** For D = Σ λ_a D_a with exact lifts N_a and N = Σ λ_a N_a + K (K Koszul — all syzygies of a regular
  sequence are Koszul), B_D(N) = Σ λ_a λ_b b_ab + Σ_a λ_a B_{D_a}(K); replacing λ_a·(Koszul coefficient) by free unknowns and
  pairing with exact duals y killing all correction columns and F S gives necessary quadrics y^T b(λ) = 0. Squares λ_a^2 in the
  exact Q-span of these quadrics kill coordinates. On a single direction the relaxation is strictly weaker than the fixed-D test
  (the free unknowns of the other directions are not forced to zero), which is why surviving directions must be tested exactly.
* **Modular steps.** (i) (S_p/I_p)_k = 0 for an ideal generated by exact integer forms ⇒ (S_Q/I_Q)_k = 0 ⇒ V_C(I_Q) = {0}
  (rank over Q >= rank mod p). Sound. (ii) HF_Q(k) <= HF_p(k) for all k, so "finitely many projective points mod p" transfers
  to "finitely many projective points over Q, total length <= the mod-p value". Sound. (iii) **Not sound:** inferring
  "V ⊆ coordinate axis" from "pure powers of all other variables in the leading ideal" — see claim 10. (iv) Exact dual
  certificates y^T M = 0, y^T b = 1 verified over Q (`w2x.exact_dual_certificate` asserts both) are rigorous proofs of
  inconsistency; my independent exact ranks over Z agree on all 15 denominators tested.
* **Y2-free lemma (claim 7).** (⇐) N = (N', 0), Q = 0. (⇒) average the certificate over μ6 acting on Y2 (equations linear in
  (N,Q) for fixed D; F, A, D invariant): N_i = N_i(x, s), s = Y2^6, N_5 = Y2 a(x,s), Q = Q(x,s); Euler-normalize (shape preserved);
  E1 and E2 are polynomial identities in s. At s = -f: Ñ(f) = A D + 6 f ã and, using ∂_s of E1 at s = -f (Σ f_i ∂_s N_i = -6a - 6s ∂_s a),
  B_D(Ñ) = 5 ã D. With N' = Ñ - ã E' (E' the 5-variable Euler field, E'(ã) = d·ã): N'(f) = A D, B_D(N') = 0. Checked line by line.
* **6 → 5 reduction of the search lens** (A D^2 - D N(f) + (6/5) f (D div N - N(D)) = 0 ⇔ certificate; lift (D^2, (N_i D), Y2 B/5, 0)):
  checked both directions, including the F^1 coefficient (the ∂_F ũ(f) terms cancel). The reduced denominator is the norm
  ∏ g(D) specialized at Y2^6 = -f (degree 6·deg D in general), so no degree is preserved except for Y2-free D.

## 3. Wording to use in public

* Keep: deg D >= 7 (exact; C1/C2, now re-verified by two further independent implementations); all Y2-odd degree-7 denominators
  excluded (two independent implementations, final Macaulay step exact over Z); the Y2-free lemma (proved); the exact E1 table;
  the single-denominator exclusions (exact duals / exact ranks over Z).
* Withdraw: "support = one coordinate axis", "the Y2-free reduction does the job on the t-free part", any suggestion that the
  pure even classes 0010/0100/1000 are closed by this packet. Replace with the finite-point statement of claim 10 and point to
  `fqlx`'s addendum (m-primary Gröbner basis over F_32003 on the base complement; not replayed here) as the only claimed closure.
* Open at degree 7 (unchanged from the thread): class 1110 and mixtures of >= 2 even characters; from this packet also the
  finitely many unidentified relaxed solutions in classes 0010/0100/1000 and the pending "even"/"full" leaves.
* Fix the "extra = 228" label at d = 8 (exact non-Jacobian count 135).

## 4. Characteristic-0 Gröbner attempt (referee, bonus)

See the end of this file (appended after the run): Singular `std` over Q on the exact integer quadrics of class 0100
(5-var strict and 6-var systems, 15 coordinates, 56 distinct quadrics each) to see whether the dimension-1 / finite-point
structure persists over Q and what the Q-length is. This does not close the class either; closing it needs the actual points
(e.g. a certified rational parametrization) tested with the exact fixed-D system, or a strictly stronger relaxation.

## 5. Reproduction of the referee checks

```
cd /tmp/hodge/results/one-reflection-skew/wave2-excl-referee
python3 fixedD.py 'Y1Y2^6,Y0Y2^6,X1Y2^6,Y1f,Y0f,X1f'      # exact ranks over Z (python-flint only)
python3 e1dims_more.py                                    # E1-kernel dims d = 8 (exact), 9 (mod p)
python3 msolve_p2.py                                      # second-prime msolve on inputs rebuilt from the exact quadrics
Singular -q diag_0100.sing                                # dim 1, l36 not nilpotent, generators restricted to the axis = -68036492 l36^2
```

## Appendix: result of the characteristic-0 Gröbner attempt (class 0100)

Singular 4.3.2, `std` over Q (exact), inputs = the exact integer relaxed quadrics of the lane's pickles restricted to the
alive coordinates (15 coordinates, 56 distinct quadrics in each system), files `../wave2-excl-referee/ref5v_tor_0100_char0.sing`,
`ref6v_tor_0100_char0.sing`, outputs `*_char0.full.out`:

| system (class 0100) | reduced GB size | Krull dim | multiplicity (degree) | coordinates with a pure power in the leading ideal |
|---|---:|---:|---:|---|
| strict 5-variable relaxation | 195 | 1 | 67 | 14 of 15 (all but `l36`) |
| 6-variable relaxation | 195 | 1 | 67 | 14 of 15 (all but `l40` = Y0 Y2^6) |

So **over Q (rigorous, characteristic 0): in class 0100 the denominators satisfying the necessary relaxed quadrics form a
finite set of projective points of total length 67**, in both relaxations; the structure is identical to the mod-p one.
In the 6-variable system the point Y0 Y2^6 is one of them (excluded exactly, claim 9); the others are not identified, so the
class is still open from this packet's point of view. The analogous characteristic-0 runs for classes 0010 and 1000
(19 coordinates) were still running (1200 s cap each) when this packet was synced; see `singular_char0.out` for whatever finished.
