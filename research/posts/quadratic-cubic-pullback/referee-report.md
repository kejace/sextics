# Referee report — lane `quadratic-cubic-pullback` (#175), wave-2 EXCLUSION / RIGOR packet

Referee: Claude (Anthropic) acting for `kejace`; no human review; nothing Lean-verified; no certificate claimed.
Inputs: draft post `wave2-quadratic-cubic-pullback-draftpost.md`, report `wave2-quadratic-cubic-pullback-excl.md`,
packet `results/quadratic-cubic-pullback/wave2-excl/` (README, 12 scripts, 90 master-log receipts, KH bases, Singular
receipts), the parallel verify lens `wave2-quadratic-cubic-pullback-verify.md`, the aggregator's `wave2-status.md` §2.11 and
`wave2-aggregator-hessrelax6.md`, `wave2-verification.log.jsonl`, the dossier, and the saved thread (one comment, `benediamond`).

Own replays: `results/quadratic-cubic-pullback/wave2-excl/referee/` (scripts `referee_replay.py`, `referee_hess7.py`,
`referee_identity.py`, `referee_lifts.py`, `referee_hess8.py`; receipts `referee-replay.jsonl`; Singular `referee_deg6_exact.sing`
→ `referee_deg6_singular.out`).  Fresh prime 1000000009 for every modular replay; exact (fmpq / Fraction) where stated.
Engine: python-flint 0.9.0, Singular 4.3.2, msolve 0.6.5 (msolve not used by the referee).

Legend.  **VERIFIED** = logic checked and (replayed by the referee, or two independent code bases agree).
**PLAUSIBLE** = receipt present and logic sound (including Hilbert-certified modular runs) but a single run, not replayed.
**UNSUPPORTED** = stated somewhere in the lane material but not established; must not appear in the post.

## Claim table

| # | claim (as in the draft post / excl report) | receipt | referee replay / independent check | verdict |
|---|---|---|---|---|
| 1 | `F = Σ Q_i^3` exactly (six Y-even quadrics) | reverify log 1; verify C1 | exact, 0 residual terms | VERIFIED |
| 2 | `h^{3,1}(X)^+ = 226 − 18 = 208`; `A ∉ J_6`, `δ ∉ J_6`, `A, δ` independent mod `J_6` | reverify-h31 (p = 2147483629); verify C7/C8 exact | exact fmpq ranks 18 / 19 / 19 / 20 over the 226 odd sextic monomials | VERIFIED |
| 3 | projective diagonal symmetry group `= Z/2 = {1, τ}` | reverify-symmetry; verify C4/C5 (own SNF) | own Smith normal form of the exponent-difference lattice: invariants (1,1,1,1,2), rank 5 ⇒ projective group `Z/2`; verify lens C4/C5 agree | VERIFIED |
| 4 | `F` is not a join: 21 second partials independent | reverify-notjoin (mod p); verify C10 exact; aggregator exact | exact fmpq rank 21 | VERIFIED |
| 5 | the 56 cubics in `Q_i` meet `J_6^{even}` exactly in `C·F` | reverify-cubicsQ (p₂); verify C12 (10⁹+7) | own exact fmpq ranks: cubics 56, `J_6^{even}` 18, sum 73 ⇒ intersection dimension 1, and `F` lies in both ⇒ `C·F` | VERIFIED |
| 6 | two-plane (4,4)-cancellation system has nullity 0 (exact) | reverify-twoplane44 (fmpq) | not replayed; single exact run | PLAUSIBLE |
| 7 | `K_d = J_d` for `d ≤ 6` both parities, Y-even `d ≤ 8` (primes 2³¹−1, 2147483629; J blocks at Hilbert rank) | kernel_sweep.jsonl (wave 1), smoke runs; aggregator d = 7, 8 | own code at 1000000009: `K_d = J_d` for `d = 0..7` both parities, every `J_{d+6}` block at its Hilbert value (own Koszul series); `d = 8` both parities replayed (see §Replay log) | VERIFIED |
| 8 | `K_d = J_d` for `d ≤ 9` both parities (search lens, 2³¹−1) | wave2-search.log.jsonl kernel items with `rank_J_target` = Hilbert (3441, 5070, 7129) | transfer valid because the J ranks are recorded at the Hilbert values; `d ≤ 8` covered by #7; `d = 9` single certified run | VERIFIED (d ≤ 8), PLAUSIBLE (d = 9) |
| 9 | `dim K_10/J_10 = 294 (even) / 261 (odd)` exactly | search kernel items d = 10 (`rank_J_target` 9684 / 9714 = Hilbert) | lower bound `dim R_10 − dim R_16 = 762−468 / 744−483` recomputed exactly (own Hilbert series); upper bound = one certified modular run | PLAUSIBLE (equality), VERIFIED (lower bound) |
| 10 | E2 reduction `N = A w + H z/6 + c·∇F`, `E = F(Φ(w) + L(c,w))`, all `H` terms cancel; `Φ ≡ −A (wᵀ Hess F w)` mod `J` | relax selftest; verify C35 | own derivation + exact random check at `d = 5, 6, 7` (`referee_identity.py`); Euler identities `zᵀHw = 5 w(F)`, `zᵀHz = 30F` exact | VERIFIED |
| 11 | **no certificate has `deg D = 5`** (exact relaxation `π_W(ker) = 0`, all parity mixtures; independently `K_H = 0`) | relaxL5a/relaxL5o/smokerelax (exact integer nullspace); hessL5 | own: `K_H = 0` at 1000000009 (joint 3003-column system, `J_10` blocks 750/747 = Hilbert); hand proof: `K_4 = 0` (both parities, replayed) + rank 21 of the second partials ⇒ `wᵀ Hess F w = 0 ⇒ w = 0`; verify C27 agrees | VERIFIED |
| 12 | degree 6, Y-odd polars: `dim K_H = 18` (certified), exact lifts, cone `= {0}` | hessL6 par1, hessX6 (Singular over Q: DIM 0) | own: all 18 lifted vectors lie in the condition exactly (`J_6` rref over Q); none has a diagonal entry ⇒ `w_a² = 0 ⇒ w = 0` by hand; own Singular (see §Singular) | VERIFIED |
| 13 | degree 6, Y-even polars: `dim K_H = 36` (certified at two primes), exact lifts, cone = Euler line `w = λz` (`D = 6λF`, not coprime) | hessL6/hessL6b par0, hessX6 (Singular over Q: DIM 1, one component = Euler line) | own: 36 lifted vectors exact in the condition; own exact 36-variable system restricted to the even sector; aggregator's independent m-primary verdict on the slice | VERIFIED |
| 14 | **degree 6, mixed parity: `dim K_H = 261` (certified), exact lifts, 405 quadrics, Singular over Q: cone of dimension 1, two random charts of length 923 with radical a single rational point (Euler), `dim(I + (L)) = 0` ⇒ cone = Euler line ⇒ every certificate has `deg D ≥ 7`** | hessL6 parall, hessX6_all_charts.out; aggregator hessrelax6 (m-primary at 10⁹+7, colength 923 cross-check) | own: 261 lifted vectors exact in the condition and independent over Q; own exact 426-quadric / 36-unknown system (`w ↦ wᵀ Hess F w mod J_6`, both parity parts, built from an exact `J_6` rref) with the Euler field checked to be a solution; transfer of the aggregator's modular m-primary verdict checked (principle 6 in README-referee); own Singular (see §Singular) | VERIFIED |
| 15 | T1 certification lemma (explicit fixed-D kernel `(Gz,0) ⊕ (gF, D div g − g·∇D)`, `n_str = 129 / 244`); three wave-1 degree-5 polar verdicts upgraded | excl T1; verify C34 | kernel elements re-derived (both families and the overlap); counts match own monomial counts; now superseded by #11 | VERIFIED (logic) |
| 16 | 17 structured degree-7 Jacobian-type candidates fail `A·d_w D ∈ J_14` (J block 5052 = Hilbert) | hesstest.jsonl (p = 2147483629) | own code at 1000000009: `J_14^{odd}` rank 5052 = Hilbert; `gradX0Q0`, `wQX`, `wQ3e3` (Y-odd control) all fail (stacked rank 5053) | VERIFIED |
| 17 | degree-7 Y-even polar family: `dim K_H = 1533` of 2016 (certified upper bound); cone of 483 quadrics in 63 variables undecided | hessL7 (p = 268435399, J_14 rank 5052 = Hilbert) | not replayed (877 s); no exclusion is claimed — degree 7 is **open** | PLAUSIBLE (upper bound) |
| 18 | 10 structured degree-8 candidates fail the Hessian test (J_16 block 9684 = Hilbert) | hesstest.jsonl d = 8 (hesstest8_detached.out) | remote replay (`referee_hess8.py`, 2 of 10) submitted; the J_16 block (13104 × 10164) is too large locally — see §Replay log for the outcome | PLAUSIBLE unless the remote replay lands |
| 19 | degree-6 structured candidates `E_X F`, `∇Q_k·∇F` fail the Hessian test | hesstest d = 6 | subsumed by #13/#14 (every degree-6 polar fails) | VERIFIED |
| 20 | `E_I F = Σ_{i∈I} z_i F_i ∈ J(F)` for every index set `I` (dossier §3.3 corrected) | excl report; verify C30 | trivial identity | VERIFIED |
| 21 | "every `D` of degree ≤ 6 not in `J(F)` … is excluded at every pole order (prescribed-support theorem at q = 1)" | book Ch. 1 theorem (`docs/01-certificates.md` §Prescribed support) | applicability checked: the theorem covers arbitrary homogeneous `D` (parity optional) via the squarefree equation `D_red`; its proof (mixed Hodge modules) is an external input, not re-verified | PLAUSIBLE (external theorem; wording to say "by the book's theorem") |
| 22 | dossier §3.3: "`E_I F ∈ J` only when `I` is everything" | dossier | false (see #20); the lens already corrected it | UNSUPPORTED (refuted) |
| 23 | attack-1 header: "no certificate of degree ≤ 7 exists (E1 kernel = J)" | attack1 note | `K_d = J_d` gives only `D ∈ J`; degree 7 polars are open (#17) — the correct statement is `deg D ≥ 7`, i.e. degrees 0–6 excluded | UNSUPPORTED (overclaim; not in the draft) |

Counts: VERIFIED 16 (rows 1–5, 7, 8(d≤8), 10–16, 19, 20); PLAUSIBLE 6 (rows 6, 8(d=9), 9(equality), 17, 18, 21); UNSUPPORTED 2 (rows 22, 23).

## Transfer arguments checked

1. Modular rank ≤ rational rank; a Jacobian block attaining its Hilbert value mod p has that rank over Q (regular sequence
   of six quintics, parity-split Koszul series recomputed independently: `R_6 = 218+208`, `R_10 = 762+744`, `R_12 = 885+866`,
   `R_14 = 762+744`, `R_16 = 483+468`).  Then `T ∉ J` mod p ⇒ `T ∉ J` over Q, and any kernel computed as
   `#cols − (rank[J;T] − rank J)` mod p is an upper bound for the rational kernel.  Direction of every inequality in the
   packet is right.  Note (also stated by the verify lens): without the certified J rank a modular `K_d = J_d` would be
   a **lower** bound and useless; every kernel receipt used (wave 1, excl, search) records the J rank at the Hilbert value.
2. E1 ⇒ `A D ∈ J(F)` for any certificate (reduce `E = F²Q` mod `F`, `D` coprime to `F`); parity-split because `J` is
   τ-stable and `A` is Y-odd, so mixed-parity `D` is covered.  `K_d = J_d` ⇒ `D = w·∇F` (degree-`d` part of `J` is spanned
   by `m F_i`; at `d = 6` the 36 generators are independent).
3. E2 ⇒ `A·(wᵀ Hess F w) ∈ J_{2d}` (row 10).  At `d = 6`, injectivity of `[A]` on `R_6` (= `K_6 = J_6`) turns this into
   `wᵀ Hess F w ∈ J_6` both parity parts; the exact rational quadrics are then cheap (used by the aggregator and here).
4. Exact lift pins `K_H(Q)`: lifted vectors verified exactly and independent over Q, their number equals the certified
   upper bound.  Cone decided by Singular over Q (`dim std = 0` ⇒ `{0}`; `dim 1` + chart with radical a single rational
   point + `dim std(I + (L)) = 0` ⇒ exactly one line).  Chart logic checked: a 1-dimensional cone is a finite union of
   lines through 0; every line not inside `{L = 0}` meets `{L = 1}` once; `dim(I + (L)) = 0` means no line lies in `{L = 0}`.
5. m-primary leading ideal mod p ⇒ no projective zero over C: Macaulay matrices of the exact integer quadrics reduce
   mod p to those of the mod-p quadrics (the mod-p quadrics span the reduction of the exact span when the `J_6` block is at
   Hilbert rank), so `HF(Q[w]/I_Q) ≤ HF(F_p[w]/I_p) = 0` for large degree.  Translation invariance under `w → w + hE`
   (`zᵀHw = 5w(F)`, `zᵀHz = 30F`, both in `J`) makes the slice `c_55 = 0` legitimate: the cone is a cylinder over the
   Euler direction, so `C ∩ {c_55 = 0} = {0}` gives `C = Q·E`.  The Euler line is excluded by the contract (`D = 6λF`
   is not coprime to `F`).
6. Farkas duals / sum principle: not used by this lane's packet (no served-space or positivity statement is posted).
   Prescribed-support powers: the book theorem (row 21) is the only source for "every pole order"; the kernel facts alone
   exclude `D` itself, not `D^k`, `k ≥ 2` (degree ≥ 12).

## Reconciliation with the verify lens (`wave2-quadratic-cubic-pullback-verify.md`)

* Agree on all 23 confirmations, on the refutation of dossier §3.3 (row 22), and on the weakening of the attack-1 header
  (row 23) and of "the target forces `D` Y-even" (Y-odd `D` is legal; the post states results for every parity mixture).
* The verify table (00:3x UTC) lists mixed-parity degree 6 as *not covered*; that was true of C28 at the time.  It is now
  closed three ways: the excl lens's exact Singular charts (R2'), the aggregator's m-primary leading ideal at 10⁹+7, and the
  referee's exact 426-quadric system.  No disagreement remains: `deg D ≥ 7` for arbitrary `D`.
* C31 (sector label of the 1491×201 relaxation): both sectors are in the master log (lines 12–13 parity 0, 15–16 parity 1,
  17–18 all); the post cites the all-parity run.  C32/C33 handled as above.  The verify lens's reading that the
  "Gorenstein" lower bound is plain rank–nullity is adopted.
* Two things the verify lens asked not to post were not posted: the one-prime compressed-consistency statement
  "join012 serves exactly the δ line", and the modular "served dimension 18" numbers.

## Replay log (own, fresh prime 1000000009 unless exact)

See `referee/referee-replay.jsonl`.  Summary: exact identities (rows 1, 2, 4, 10, 12–14 lifts); J block ranks at Hilbert
for degrees 6, 10, 11, 12, 14 (both parities where applicable); `K_d = J_d` for `d = 0..8` both parities; degree-5 `K_H = 0`;
three degree-7 candidates excluded; Singular over Q on the own degree-6 system — results in the lines below.

* Own degree-6 system (426 exact quadrics, 36 unknowns), Singular mod 1000000009: `FULL36 dim 1`, `SLICE_c35 dim 0`,
  `EVEN_SECTOR dim 1`, `EVEN_SECTOR_SLICE dim 0`, `ODD_SECTOR dim 0` ⇒ (transfer principles 5–6) cone = Euler line over C for
  every parity mixture; Y-odd sector `{0}`.  Singular over Q: still running at snapshot time (not needed: the mod-p run transfers).
* Own exact replays added after the table was drafted: cubics in `Q` ∩ `J_6^{even}` = `C·F` (fmpq ranks 56 / 18 / 73); symmetry
  group `Z/2` (Smith invariants 1,1,1,1,2); `K_8 = J_8` both parities (J_14 blocks 5052 / 5070 = Hilbert).
* Not replayed: structured degree-8 tests (remote stayed above the capacity gate all session; row 18 stays PLAUSIBLE);
  degree-7 full-family kernel dimension (row 17); two-plane nullity (row 6).

## Verdict summary

VERIFIED 16 · PLAUSIBLE 6 · UNSUPPORTED 2 (both are earlier internal statements already corrected; neither appears in the post).
Strongest verified statement: **every certificate for `(F, 1, X0^3 Y0^3)` has `deg D ≥ 7`** — E1 kernels `K_d = J_d` for `d ≤ 8`
(three primes, three codes), the exact degree-5 argument, and the degree-6 Hessian cone equal to the Euler line in every parity
mixture (lens Singular over Q, aggregator leading ideal, referee's exact quadrics + Singular with checked transfer).
Wording downgrades applied in the post: "every pole order" attributed to the book theorem (external); `d = 9` and the `d = 10`
equality and the degree-8 structured exclusions marked as single certified runs; degree 7 stated as open; no served-space or
compressed-consistency statements posted.
