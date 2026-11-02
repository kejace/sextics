# README-referee — replaying the `eight-plane-configuration` wave-2 exclusion packet (2026-11-01)

Written by the referee pass (Claude, Anthropic, for `kejace`); it does not modify the agent's files. It records what a stranger
needs beyond `README.md` to replay every statement S1–S6 with the pinned engine, and where the referee's independent replays live.

## Environment (pinned, as actually used)
* `python-flint 0.9.0` is installed in the **user site of `/usr/bin/python3`** (3.12); it is *not* in `~/devel/sextics-work/venv`.
  Run every script with `/usr/bin/python3` (not `-I`, which hides the user site).
* `msolve 0.6.5` (`/usr/bin/msolve`, Debian `msolve 0.6.5-1build2`); `Singular 4.3.2` (`1:4.3.2-p10`). Leading ideals: `msolve -g 1`.
  Write squares as `w^2` in `.ms` inputs (msolve silently misparses `w*w`); all inputs in this packet use `^2`.
* Engine: `~/devel/sextics-work/hodge` (`core.py`, `fixed_d.py`, `jacobian.py`); a snapshot is in `research/engine/` of the
  research branch (fork `kejace/sextics`, branch `research`). The instance is read from `~/devel/sextics/data/eight-plane-configuration/instance.rec`.
* Local scripts (`w2x_*.py`) hard-code `sys.path.insert(0, '/home/kejace/devel/sextics-work')` and the lane directory
  `/tmp/hodge/results/eight-plane-configuration` (`WD`) with the packet at `WD/wave2-excl` (`PK`). Recreate that layout or edit `w2x_lib.py`.
* Inputs that live **outside** the packet: `WD/K_decomp.json` (the quartics `K_i`, `H_i` with `F = Σ Q_i K_i`), the wave-1 kernel receipts
  `WD/e1kernel.jsonl` and `/tmp/hodge/remote/eight-plane-configuration/e1kernel_remote.jsonl` (used by `w2x_hilbert.py`), and the degree-10
  kernel bases `WD/e1kernel_basis_d10_par*.pkl` (p = 2147483629, used by `w2x_hess8_local.py`).
* Remote-only files (gcloud box, mirrored under `/tmp/hodge/remote/eight-plane-configuration/wave2-excl/`): `wave2-excl-remote-*.jsonl`,
  `excl-epc-batch2a-*.out` (exact Y-partials, the 140 degree-7 samples), `annA_d10_par*.p2147483587.pkl`.

## Statement-by-statement replay (packet `README.md` statements S1–S6)
| statement | receipt(s) in the packet | replay command | referee replay (independent code, `../wave2-excl-referee/`) |
|---|---|---|---|
| S1 `A ∉ J(F)` | (argument only; wave-1 rank 18/19) | — | `ref_exact.py`: exact fmpq rank 18 vs 19. Note: `README.md` says "dim J_6 = 126"; the correct number is 36 (18 per parity), the mod-5 saturation argument is unaffected. |
| S2 E1 kernel = `J_d`, d ≤ 9, both parities | `e1ladder.2147483587.jsonl` (= `receipt-ladder5to10.jsonl`), wave-1 `e1kernel*.jsonl`, `w2x_hilbert.py` output | `python3 w2x_hilbert.py` (exact Hilbert function + no-drop check of every wave-1 row); `python3 w2x_e1ladder.py 5 6 7 8 9 10` (normal forms, p = 2147483587; d = 9, 10 need ~10 min each) | `ref_e1.py 1000000009 4 5 6 7 8 9` (ranks, 4th prime; each row records `rank_J_p = exact_J_dim`) |
| S3 Jacobian reduction + Koszul-reduced E2 | proofs in `wave2-eight-plane-configuration-excl.md` §5 | — | re-derived in `referee-report.md` (R3–R5 checked) |
| S4 d = 5 | `hess.jsonl` rows `d=5` (linearization rank 6), `hess_mixed.jsonl` row `d=5` (rank 21) | `python3 w2x_hess.py 5 0`, `… 5 1`; `python3 w2x_hess_mixed.py 5` | `ref_exact.py`: exact rank of the 21 second partials = 21 (any parity mix) |
| S4 d = 6 pure | `hess_d6_par{0,1}_p{2147483587,2147483647}.{ms,lm,msolve.log}`, `hess.jsonl` | `msolve -g 1 -f hess_d6_par0_p2147483587.ms -o out.lm` (0.02 s); Hilbert function from the leading ideal | `ref_hess6.py 1000000009 even odd`: own build, msolve + Singular `dim 0`, colength 38 / 39 |
| S4 d = 6 mixed | `hessmixed_d6_p2147483587.{ms,lm,msolve.log}`, `hess_mixed.jsonl` (`msolve_quotient_dim` 923) | `python3 w2x_hess_mixed.py 6` (builds + runs msolve, ~10 s) | `ref_hess6.py 1000000009 mixed`: 35 vars, 426 quadrics (405 independent), m-primary, Hilbert function [1,35,225,400,225,36,1], colength 923; Singular `dim 0`, `vdim 923` |
| S5 d = 7 full pure-parity | `hess_d7_par{0,1}_p2147483587.ms` (inputs only); `.lm` files are **empty** (msolve OOM-killed, rc 137) | `bash run_msolve.sh hess_d7_par0_p2147483587.ms <threads> <tag>` (needs ≥ 60 GB) | not replayed (capacity); result still **pending** |
| S5 d = 7 subfamilies | `hess_sub.jsonl`, `hsub_d7_par*_{arr,blockX,blockY,coefX,coefB0}_p2147483587.{ms,lm}` (`grad`: linearization only) | `python3 w2x_hess_sub.py 7 <par> <family>`; `msolve -g 1 -f hsub_….ms -o out.lm` | `ref_grad7.py` (rank 406 = dim Sym², both parities, new prime); msolve replays of `coefB0`, `arr`, `coefX` with independent leading-ideal reading (`referee.log.jsonl`) |
| S5 d = 8 subfamilies | `hess_sub8_local.jsonl`, `hsub_d8_par*_{coefB0,coefX,comp-i}_p2147483629.{ms,lm}` | `python3 w2x_hess8_local.py …` (needs the d = 10 kernel bases) | msolve replay of `hsub_d8_par0_coefB0` (colength 28) |
| S6 sampling (screened) | `receipt-batch2a-partial.jsonl` (140 degree-7 fields, `rigorous_rank_full: false`), remote `excl-epc-batch2a-jsample7.out` | `python3 w2x_jsample.py 7 <n> <tag>` | not replayed; single-prime modular inconsistency = screening only (no Farkas dual recorded) |
| factor table (sum principle) | `WD/wave2-excl.log.jsonl` rows `exp = c-factor-annihilator` (36 factors; **no no-drop rank recorded**) | `python3 w2x_factors.py` | `ref_factors.py 1000000009`: all 36 annihilators 0 with `rank_p J_{6+e} = exact` recorded ⇒ rigorous over ℚ |

## Transfer arguments (why modular results are proofs here)
1. Projected kernels (E1 kernel, factor annihilators): `dim_ℚ ker ≤ n − rank_p[A·m | J] + dim_ℚ J = dim_p ker` whenever `rank_p J = dim_ℚ J`;
   the kernel contains `J_d` of known dimension, so equality of dimensions proves `ker = J_d` over ℚ. Every row used records the no-drop check.
2. Hessian relaxations (msolve/Singular): with `rank_p J_{2d−6} = dim_ℚ J_{2d−6}`, the mod-p quadric space is the reduction of the rational one;
   an m-primary leading ideal means the degree-k Macaulay matrix has full column rank mod p for some k, hence over ℚ, so `V'(ℂ) = {0}`.
   (For d = 8 the same holds with `Ann_R(A)_10 + J_10` in place of `J_10`, because the forced kernel dimension 294/261 is attained.)
3. Exact ranks (`fmpq`/`fmpz`): rank over ℚ of the Koszul-reduced system (equivalent to the full fixed-D system) decides consistency over ℂ.
4. Sum principle at q = 1: book Ch. 1 "prescribed-support theorem" (poles confined to V(D) ⇔ denominator `D_red^q = D_red`) and Ch. 3 Gysin/weight
   description (the (3,1) part of `N_{Y₁∪Y₂}` is the sum of the components' parts). Hence `served(D₁D₂) = served(D₁) + served(D₂) mod J`.
   Single-prime "inconsistent" verdicts (S6) are **not** transferable without `rank_p = rank_ℚ` or a rational Farkas dual.

## Referee files
`../wave2-excl-referee/`: `ref_lib.py` (own sympy parser, Koszul Hilbert function, sparse flint matrices), `ref_exact.py`, `ref_e1.py`, `ref_hess6.py`
(+ `ref_hess6_{even,odd,mixed}_p1000000009.{ms,lm,sing,msolve.log}`), `ref_grad7.py`, `ref_factors.py`, msolve replays `ref_replay_*.lm`,
log `referee.log.jsonl`, verdicts `referee-report.md`.

Note (referee, 22:10 local): `ref_e1.py 1000000009 10` (forced d = 10 kernel, ~1.2 GB) was still running at the time of the commit linked from the
discussion; its row lands in `referee.log.jsonl` / `ref_e1_d10.out` in a later cron sync of the research branch.
