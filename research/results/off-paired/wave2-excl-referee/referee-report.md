# Referee report — lane `off-paired`, wave-2 exclusion draft

Referee: Claude (Anthropic) for `kejace`, 2026-11-01/02. Inputs: draft post `wave2-off-paired-draftpost.md`, report
`wave2-off-paired-excl.md`, packet `results/off-paired/wave2-excl/` + lane log `wave2-excl.log.jsonl`, verify notes
`wave2-off-paired-verify.md`, status `wave2-status.md` §2.8, dossier, saved thread. Independent code and logs:
`results/off-paired/wave2-excl-referee/` (python-flint 0.9.0 only; no engine imports). Fresh prime used: 1000000009.

Verdict scale: VERIFIED = logic checked and the computation replayed (exactly, or modular with a transfer argument I checked);
PLAUSIBLE = modular/screen only; UNSUPPORTED = no valid receipt as recorded.

## Claim table

| # | Draft claim | Receipt | Referee replay | Verdict | Wording to use |
|---|---|---|---|---|---|
| 1a | $D_0=F_{X_0}/(6X_0)$: ranks 67/68, nullity 8 = explicit kernel; theorem | log `D0` fixedD at 2 primes | exact over Q: 344×75, 67/68, nullity 8 = explicit 8; identity $6X_0D_0=F_{X_0}$ exact; E1 passes 69/69 | VERIFIED | as drafted (re-derivation of `fqlx`) |
| 1b | degree-6 sectors chart-complete (msolve) | `receipts_test-local-d6` (8 charts) | inputs write squares as `t*t`; msolve 0.6.5 reads `t*t` as `t` (checked) | UNSUPPORTED | drop; degree ≤ 6 is `fqlx`'s exact theorem anyway |
| 1c | $f_6,g_6$ fail E1 (rank 176 = Koszul dim, aug 177) | log `f6` e1 | exact over Q: 216×192, 176/177 for $f_6$, $g_6$, $f_6+g_6$ | VERIFIED | theorem (also exact in the verify lens) |
| 2a | degree 7: all 36 characters chart-complete `no-solution` mod 2147483647 (171 charts) | `receipts_excl-off-paired-bilin-d7-local` | all 171 inputs affected by `t*t`; corrected rerun (`^`): 171/171 `[-1]` (69 s) | UNSUPPORTED as recorded; PLAUSIBLE after rerun | "modular screen (corrected inputs), consistent with the exact result below" |
| 2b | degree 7 exact over Q for 34 of 36 characters (relaxation pipeline) | `receipts_excl-off-paired-elim2Q-d7` (PROVEN records) | sectors (1,2),(2,2),(1,4),(2,4): colon dim 1, gauge 0, unique $D$ inconsistent exactly (188/189, nullity 18 = explicit) by independent code; Singular over Q replayed: `DIMJ2 0` for (0,3),(0,5),(1,1),(3,1); (0,0) and (2,0): `DIMC 0` for both localized relaxations, (2,0) J2 has dim 2 / one component; empty sectors (k,3): nullity_p ≥ nullity_Q is sound | VERIFIED | "exact for 34/36 characters; (1,0),(3,0) open (Singular timeouts)" |
| 3a | degree 8: all 36 characters chart-complete, **including the six (a,0)** | `receipts_excl-off-paired-bilin-d8-local` (414) | all 414 inputs affected by `t*t`; corrected rerun running (`charts_fixed_d8/results.jsonl`) | UNSUPPORTED (as recorded) | withdraw the (a,0) statement; report the corrected-rerun count as a screen only |
| 3b | degree 8 exact for 21 characters | `receipts_excl-off-paired-elim2Q-d8` | transfer sectors (1,5),(2,5),(4,1),(5,1): exact quadrics reduced mod 1000000009, msolve GB has a pure power of every variable (13/13, 14/14, 13/13, 14/14); Singular `DIMJ2 0` for (0,2),(3,3),(1,3),(4,4) | VERIFIED (sampled; same exact pipeline for the rest) | "exact for 21/36; all 21 lie inside `deanodotdotdot`'s 30, so no new degree-8 class is closed" |
| 3c | "9 attempted sectors have a positive-dimensional relaxation residual" | same file | receipts say `singular-timeout(J2)` (8) and `rational-data-failed` ((4,2)): dimension unknown | UNSUPPORTED (wording) | "8 Singular timeouts, 1 reconstruction failure; six (a,0) never attempted" |
| 4 | structured fixed denominators PROVEN (block products of degrees 7–14, wave-1 families, Hessian minors; two primes) | lane log, job `excl-off-paired-batch1-blocks` | recount: 350 distinct, 346 PROVEN (343 block products of degrees 7–13 = 347/680 of the list, + $D_0$, $f_6$, $F_{X_0}X_1Y_1$), 4 modulo p only; replays: `X1*X2*FX1` exact 245/246 (nullity 30 = explicit), `FX0*X1*Y1` exact 287/288, `X1*D0*FX2` (deg 10) at 1000000009: 920/921, nullity 126 = explicit 126 (sandwich holds), `X0*H*C` E1 exact 835/836 | VERIFIED | "343 block products of degrees 7–13 (347 tested) + 3; wave-1 (621) and Hessian (213) families NOT processed" |
| 5 | colon / non-Jacobian kernel dimensions per sector, d ≤ 12 | `kernel_dims_p2147483629.json` (d ≤ 9) + log (d = 10, 11) | modular table; totals 1, 11, 54, 177, 457, 1014, 2173, 3753 equal the exact Gröbner colon dims (dossier 1.7, verify lens #7); since nullity_p ≥ nullity_Q sector-wise and the totals agree, the per-sector values are exact | VERIFIED (d ≤ 11) | "d ≤ 11" |
| 6 | degree 9 (0,0): 27/50 charts closed | `receipts_excl-off-paired-bilin-d9-local-probe` | `t*t` inputs | UNSUPPORTED | drop |
| 7 | $\Delta_t$ full solve | none | — | no claim | "not run" |

Counts: VERIFIED 6 (1a, 1c, 2b, 3b, 4, 5), PLAUSIBLE 1 (2a after the referee's corrected rerun), UNSUPPORTED 4 (1b, 3a, 3c, 6).

## Soundness of the transfer arguments (checked by hand)
* **E1 failure mod p ⇒ over Q.** $F$ smooth (Milnor number $5^6$, verify lens), so the partials are a regular sequence and
  $\dim_{\mathbb Q} J(F)_k[\chi]$ is the equivariant Koszul value; $\operatorname{rank}_p \le \operatorname{rank}_{\mathbb Q}$ gives
  $\operatorname{rank}_{\mathbb Q}[J\,|\,AD] \ge \operatorname{rank}_p + 1 > \dim J_k[\chi]$. Sound.
* **Rank sandwich.** $(LZ,0)$ and $(FW, B_D(W))$ are kernel vectors ($B_D(LZ)=(\deg L+6-d)LD$, $B_{D}(FW)=DW(F)+FB_D(W)$);
  their intersection is $\{(FL'Z,0)\}$ because $B_D(L'Z)=0$ for $\deg L'=d-6$, so $e=s_d+\sum_i s_{d-5}-s_{d-6}$ is exact. With
  $\operatorname{rank}_p(M)=C-e$ the rational kernel is the explicit one and $\operatorname{rank}_p[M|b]=C-e+1$ forces inconsistency over
  $\mathbb Q$, hence over $\mathbb C$. Kernel-certified records verify the full reconstructed kernel exactly, same conclusion. Sound.
* **Prescribed-support theorem ($q=1$).** Each excluded $D$ is squarefree and coprime to $F$, so every denominator with the same
  reduced support is excluded; products of excluded supports are not covered (the draft says so). Sound.
* **Homogeneous $J_2$: m-primary mod $p$ ⇒ over $\mathbb Q$.** Hilbert function over $\mathbb Q$ ≤ mod $p$ (integer generator matrix,
  rank can only drop), so $(J_2)_k=S_k$ for large $k$ over $\mathbb Q$; m-primary ⇒ no nonzero complex zero ⇒ $t=0$ ⇒ $D\in F\cdot S_{d-6}$,
  illegal. Pure power of every variable in a GB ⇒ m-primary. Sound. Chart `GB = {1}` does not transfer (px−1), as the packet says.
* **Component route.** $V(J_2)=\bigcup V(P_i)$ (minAssGTZ over $\mathbb Q$), $V(P_i)\subset L_i$ (linear generators), exact
  pure-$x$ consequences of degrees 2–4 on $L_i$ m-primary ⇒ $x=0$. Sound.
* **Sector restriction / Euler gauge.** Lossless by averaging; $D\mapsto D+FL$, $N\mapsto N+AL\epsilon/6$ preserves the system. Sound.
* Farkas/constant-cofactor chart certificates (`w2_cert.py`): none succeeded (22 `no-certificate` records); nothing rests on them.
  Sum principle: not used in a rigorous claim.

## Reconciliation with the other lenses
* Verify lens (15 confirmed / 0 refuted / 1 weakened): no overlap conflict; its exact $D_0$, $f_6$, $g_6$ results agree with mine.
* Aggregator (`wave2-status.md` §2.8) suspected the chart sweeps were "mixed" by the `c*c` pitfall: in this packet they are
  uniformly affected (the 22 clean files are the other agent's `bilinear_d6` on the remote). The status's "343 tested / 340 PROVEN" is
  an earlier count; the log now gives 350 / 346.
* No disagreement on any exact statement. The draft's headline (degree-8 (a,0) sectors closed modulo p) must be withdrawn.

## README gaps fixed
`README-referee.md` in the packet: pinned engine and paths, where the batch receipts live (lane log), recomputed counts, the
`t*t` defect and corrected-rerun location, which receipts files do not exist, engine-free replay instructions.

## Late additions (end of referee pass)
* `(0,0)` degree-7 J2 over Q replayed to completion: `DIMJ2 2`, one nonlinear component (as in the receipt); together with
  `DIMC 0` for `d7_c00_comp0_local.sing` the (0,0) proof replays. Same for (2,0).
* Degree-8 corrected-syntax chart rerun (`ref_charts_fixed.py`, smallest charts first, 120 s cap, p = 2147483647): 220 of 414 charts
  finished, 220 report no solution, 0 timeouts, at the time of writing; it keeps running in the background and appends to
  `results/off-paired/wave2-excl-referee/charts_fixed_d8/results.jsonl` (copied from `/tmp/hodge/referee/off-paired/`). **Refresh the
  count in `comment.md` before posting.**
* **Sync blocker (must be fixed by a human or the orchestrator before posting).** `~/devel/sextics-research` (worktree of
  `~/devel/sextics/.git`) has been in a stopped rebase since 2026-11-01 19:10 local (`rebase-merge/`: onto `82a69660`, stopped at
  `470e9f9a` "sync brunos 2026-11-02T00:10Z`); branch `research` still points at `470e9f9a`, so every `sync.sh` since then committed
  onto a detached HEAD and `git push fork research` was rejected (non-fast-forward). The referee's local commits (`d1b251a0`,
  `589fac98`, and the one made at the end of this pass) are therefore NOT on the fork. The referee did not run any git command other
  than `sync.sh` and `rev-parse`, as instructed. The files on disk (`/tmp/hodge/results/off-paired/wave2-excl/README-referee.md`,
  `/tmp/hodge/results/off-paired/wave2-excl-referee/`) are the source of truth; after `git rebase --abort` (or resolving the conflict) and
  one `sync.sh`, replace the branch link in `comment.md` by the resulting SHA link.
* Packet links in `comment.md` use the fork's current tip `18dd4ae4` (remote box sync 02:50Z), which was checked read-only to contain
  the lane log, both `elim2Q` receipt files (byte-identical to local), the block-product receipts and the Singular scripts; its
  `README_w2.md` is the 21:02 version (without the "Exact (theorem-level)" section).
