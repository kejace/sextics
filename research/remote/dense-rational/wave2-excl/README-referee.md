# Referee notes for the `dense-rational` wave-2 exclusion packet (2026-11-01, local 21:50)

Written by the referee pass (Claude, Anthropic, for `kejace`); the agent's files in this directory were not modified.
Referee replays live in `referee/` (scripts `ref_*.py`, msolve inputs/outputs, one JSON line per check in
`referee/referee.log.jsonl`). The claim-by-claim verdict table is `referee/referee-report.md`.

## 1. What a stranger needs that the agent's README does not say

* **Paths.** The agent's scripts hard-code `/tmp/hodge/results/dense-rational` (lane dir, containing `dr.py` and
  `wave2_excl_*.py`) and `/home/kejace/devel/sextics-work` (engine package `hodge`). In the research branch these are
  `research/results/dense-rational/` and `research/engine/`. `hodge.core.load_instance` reads
  `ROOT/data/dense-rational/instance.rec` with `ROOT = /home/kejace/devel/sextics` (`core.py` line 272): a checkout of
  `hodge-computer/sextics` (any commit; `instance.rec` is immutable) must sit there, or `ROOT` must be edited.
  The scripts import `dr.py` (lane helpers) through `wave2_excl_common.py`; the engine notice about the pre-17:25 sign
  bug does not affect anything here (no `FixedDSystem` result is used; the sign of the relaxation identity was checked
  separately, see `relaxation-sign-check` in `../wave2-excl.log.jsonl` and `referee-relaxation-identity-modp`).
* **Pinned tools (as actually installed on the machine that produced the receipts).** Python 3.12.3, python-flint 0.9.0,
  msolve 0.6.5 (Ubuntu package `0.6.5-1build2`, `/usr/bin/msolve`), Singular 4.3.2 (not used by this packet).
  msolve flags used: `-g 1` (DRL leading ideal only), `-t <threads>`, input files are homogeneous with squares written
  `c^2` (msolve misparses `c*c`); the certificate is the leading-ideal block of the `.g1.out` file.
* **Prime bookkeeping.** `e1kernel_*` receipts exist here only for `d <= 6` at `p = 2147483587` (both parities). The
  `d = 7..10` sweep at that prime was queued in `wave2_excl_chain_rig.sh` but **never started**: `submit_rig.log` shows
  the remote over capacity (load 570-1970) for the whole evening. Other primes covering `d = 7, 8, 9`: wave 1 (even
  sector, 2147483647, `../e1_direct_*.log`), the aggregator (both parities, 1000000007,
  `/tmp/hodge/notes/wave2-verification.log.jsonl`), the dossier (d = 9, 101 and 10007), and the referee (both
  parities, 1000000009, `referee/referee.log.jsonl`).
* **Empty outputs are not certificates.** `fam_d7_mixed_xonly_p1073741789.g1.out` (0 bytes) and, on the remote,
  `quadsys_d7_{even,odd}_p1073741789.g1.out` (0 bytes) are unfinished msolve runs: the even job was killed inside the
  F4 degree-4 step (matrix ~870767 x 500646, ~39 GB), the odd job was still in that step after 2.5 h. The lane's
  `wave2_excl_ltcheck.py` reports `n_leading_monomials: 0` and `projectively_empty_mod_p: false` on such files; the
  referee parser reports `m_primary_mod_p: false` with reason "msolve did not finish".
* **Implicit rank cap.** `wave2_excl_e1kernel.py` checks `rank_p(J_d block) = dim J_d` but not `rank_p(N block) =
  dim J_{d+6}`, which the saturated-lattice transfer also needs (the rational solution lattice reduces *into* the mod-p
  solution space only when the generator matrix of `J_{d+6}` keeps its rank mod p). The logged numbers satisfy it:
  `rank_p(N) = dim_Kfull + rank - nD` equals `dim J_{d+6}^{1-parD}` for every `d <= 6` (63, 63, 168, 168, 378, 378,
  747, 750, 1341, 1341, 2214, 2223). The referee script `ref_e1.py` records both caps explicitly.
* **Gauge inside a family.** `wave2_excl_quadsys_fam.py` computes `W_fam ∩ gauge` mod p only. The transfer to `C`
  needs the intersection to have the same dimension over `Q`; `ref_fam7.py` computes it exactly (it is 6 for `weuler`,
  i.e. the whole gauge, and 0 for `blockdiag`, `cross`, `gradient`) and agrees with the mod-p value.

## 2. Replay (referee scripts; independent of the agent's generators except for parsing `instance.rec`)

```
cd research/results/dense-rational/wave2-excl/referee          # or /tmp/hodge/results/dense-rational/wave2-excl/referee
export PYTHONPATH=/home/kejace/devel/sextics-work               # engine (research/engine); ROOT in core.py must point at a repo checkout
python3 ref_ltparse.py ../quadsys_d6_*_p1073741789.g1.out ../fam_d7_mixed_*_p1073741789.g1.out   # re-parse the agent's certificates (seconds)
python3 ref_quad6.py 1000000009                                 # own degree-6 relaxation, 3 sectors, msolve -g 1 (seconds)
python3 ref_e1.py 1000000009 4,5,6,7,8                          # exact fmpq checks, relaxation identity, E1 sweep with both caps (minutes)
python3 ref_fam7.py 1000000009 weuler,blockdiag,cross,gradient  # own degree-7 mixed family systems + msolve -g 1 (minutes)
```
Agent's own replay commands: `README.md` section 3 (same engine; run from the lane dir).

## 3. Referee results (details in `referee/referee.log.jsonl`)

| item | agent's receipt | referee replay | verdict |
|---|---|---|---|
| `A ∉ J_6^-`, `dim J_6^- = 18` | `A-not-in-J6` (p = 2147483587: 18 -> 19) | exact `fmpq_mat`: rank 18, augmented 19 | VERIFIED |
| `Sym^2 C^6 -> S_4` injective (degree 5 excluded for every `v`) | `degree5-hessian-kernel` (exact, rank 21) | exact `fmpq_mat`: rank 21 of 126 x 21 | VERIFIED |
| `(J:A)_d = J_d`, `d <= 6`, both parities | `e1kernel_d{1..6}_par{0,1}_p2147483587.json` | implied N-cap attained (above); own sweep d = 4, 5, 6 at 1000000009 with both caps: kernel mod J = 0, N-cap and J_d-cap attained for d = 4, 5, 6 and both parities | VERIFIED |
| `(J:A)_d = J_d`, `d = 7, 8`, both parities | none at 2147483587 (chain never ran); wave 1 even at 2147483647; aggregator both parities at 1000000007 | own sweep at 1000000009 with both caps: kernel mod J = 0 for d = 7, 8 and both parities; rank_p(N) = 3441, 3441, 5052, 5070 = exact dim J_13^∓, J_14^∓; J_d caps 63, 63, 168, 168 attained (69-221 s each) | VERIFIED |
| `(J:A)_9 = J_9`, both parities | none in packet; dossier (101, 10007), aggregator (1000000007), wave 1 even (2147483647) | not replayed by the referee (three independent primes on record) | VERIFIED (by others) |
| degree-10 kernels exactly 294 / 261 | none in packet (chain never ran); wave 1 dual sweep (even, 2147483647), dossier (both, 101) | not replayed | PLAUSIBLE (consistent across two primes; counting bound + cap argument recorded in the dossier) |
| relaxation identity `B_D(Aw) ≡ -A w^T Hess(F) w (mod J)` | `relaxation-sign-check` (quadratic w, p = 2147483587) | hand derivation (see report) + own mod-p check on a linear field at 1000000009: rank J_12 = 4437, with B + A wHw still 4437 (in J), with B - A wHw 4438 (not in J) | VERIFIED |
| degree 6 excluded, even/odd/mixed (T4) => `deg D >= 7` for every certificate | `quadsys_d6_*_p1073741789` (135/153/495 leading monomials) | own parser on the agent's outputs: m-primary; own generator + msolve at 1000000009: 135/153/530 leading monomials, pure powers of all 17/18/35 variables; also verify lens (1000003) and aggregator (1000000007, mixed) | VERIFIED |
| degree-7 mixed families `weuler`, `blockdiag`, `cross` | `fam_d7_mixed_*_p1073741789` (315/562/562) | own parser: m-primary; own generator + msolve at 1000000009: 315/562/562, all variables | VERIFIED |
| degree-7 mixed family `gradient` | `fam_d7_mixed_gradient_p1073741789` (5360 leading monomials) | own parser: m-primary; own generator + msolve at 1000000009: 5267 leading monomials, pure powers of all 56 variables (276 s, 4 threads) | VERIFIED |
| degree-7 τ-pure sectors (T5) | `quadsys_d7_{even,odd}_p1073741789.ms` only; outputs empty | not replayable (msolve did not finish; even job killed, odd job running) | UNSUPPORTED |
| degree-7 mixed families `xonly`/`yonly`/pairs/triples; degree-8 families (T6 rest); T7; screens | no outputs (xonly unfinished; rest never ran) | nothing to replay | UNSUPPORTED |
