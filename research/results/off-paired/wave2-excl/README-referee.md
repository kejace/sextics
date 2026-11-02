# Referee notes for the `off-paired` wave-2 exclusion packet (`wave2-excl/`)

Written by the referee pass (Claude, Anthropic, for `kejace`), 2026-11-01/02. The agent's files in this directory were not
modified; everything the referee produced lives in `../wave2-excl-referee/` (independent replay code, logs, corrected chart reruns).
Read this file together with `README_w2.md`.

## What a stranger needs that `README_w2.md` does not say

1. **Pinned engine.** python-flint 0.9.0 (`pip install --user python-flint==0.9.0`), msolve 0.6.5 (Debian/Ubuntu package `msolve`),
   Singular 4.3.2 (`Singular -q script.sing`). The `w2*.py` scripts import the scratch engine `hodge/` from
   `/home/kejace/devel/sextics-work` (a copy is in `research/engine/` on the research branch; `sys.path.insert` lines at the top of
   `w2lib.py`/`w2_bilinear.py` must point at the directory that contains `hodge/`). Scripts also hard-code
   `/tmp/hodge/results/off-paired/wave2-excl` as the packet path (`w2lib.PK`) and `data/off-paired/instance.rec` via
   `hodge.core.load_instance` (the sextics checkout at `~/devel/sextics`). Set `W2LOG` to a scratch file before running anything,
   otherwise the scripts append to the lane log.
2. **Where the fixed-denominator receipts are.** The 1803-line lane log `../wave2-excl.log.jsonl` (synced as
   `research/results/off-paired/wave2-excl.log.jsonl`) is the only local copy of the block-product receipts
   (`job = excl-off-paired-batch1-blocks`; the per-job file `receipts_excl-off-paired-batch1-blocks.jsonl` is only in the remote
   mirror `/tmp/hodge/remote/off-paired/wave2-excl/`, which is not synced). Each denominator has an `e1` record (prime 2147483629,
   second prime only on failure) and, when E1 passes, `fixedD` records at 2147483629 and 2147483587. Denominators are given by
   `D_expr` over the block names in `w2lib.blocks()`: `M = X2^6 - Y2^6`, `S = X2^6 + Y2^6`, `H = X0^6+X1^6+Y0^6+Y1^6+X0^3X1Y0Y1/2`,
   `C = X1^3Y1/2 + X0^2X1Y0/3`, `D0 = X0^4 + X0X1Y0Y1/4 + X1X2Y0Y2/9`, `FX0..FY2` = partials of F, `f6`, `g6`, `Qx = F - f6 - g6`,
   `P2 = X2Y2`, `M2 = X2^2 - Y2^2`, `S2 = X2^2 + Y2^2`, `R` (see `w2lib.py`).
3. **Counts as of the end of the session (recomputed by the referee from the log):** 350 distinct denominators; 346 PROVEN
   (343 block products of degrees 7-13 = 347 of the 680 in `cands_blocks.jsonl`, plus `D0`, `f6`, `FX0*X1*Y1`); 4 block products
   excluded modulo both primes only (`X1*M*H`, `Y1*M*H` of degree 13, `X1*S*FY0`, `X1*M*FY0` of degree 12: kernel lift 578/580).
   Proof routes: 67 by `PROVEN-E1-fail`, 39 by `explicit-kernel`, 240 by `kernel-certified`. The 621 wave-1 denominators and the
   213 Hessian minors were **not** processed (queued, never ran). Degree 14 block products were not reached.
4. **The msolve chart sweeps are void as recorded.** `w2_bilinear.write_ms` emits squares as `t5*t5`; msolve 0.6.5 parses
   `v*v` as `v` (referee check: `x*x-4, y-x` gives the reduced GB `[y-4, x-4]`; with `x^2-4` it gives `[x-y, y^2-4]`). Every file in
   `bilinear_*/` is affected (7,969 squared terms at d = 7, 88,538 at d = 8, 14,615 in the d = 9 probe, 119 in the d = 6 test), so the
   `no-solution(modp)` receipts in `receipts_excl-off-paired-bilin-*.jsonl` and the chart tables in the report describe the wrong
   systems. The Singular scripts in `elim2_*/` write `t0*t0` too, but Singular parses that correctly (checked), and the
   `*_J2modp.ms` files use `^`, so the exact pipeline is unaffected. A corrected rerun (`v*v -> v^2`, same prime) is in
   `../wave2-excl-referee/charts_fixed_d7/` (all 171 degree-7 charts again `[-1]`) and `charts_fixed_d8/` (see `results.jsonl`).
5. **Exact sector receipts.** `receipts_excl-off-paired-elim2Q-d7.jsonl` has several records per sector for (0,0) and (1,0); the
   record that counts is the one whose `result` starts with `PROVEN`. Degree 7: 34 of 36 characters proven ((1,0), (3,0) are
   Singular timeouts); degree 8: 21 of 36 (the six (a,0) sectors were never attempted locally; (0,1), (0,5), (1,1), (2,1), (3,1), (3,5),
   (4,5), (5,5) timed out; (4,2) failed rational reconstruction). `receipts_excl-off-paired-chain-rigor-d{7,8}.jsonl` named in
   `README_w2.md` do not exist (the remote chain was never submitted).
6. **Replaying a sector proof without the engine.** `elim2_*/d{d}_c{ab}_J2.sing` contains the exact rational quadrics; `Singular -q`
   prints `DIMJ2 0` for the m-primary sectors (eight replayed by the referee, all `0`, each < 1 s) and the components otherwise;
   `d7_c00_comp0_local.sing`, `d7_c20_comp0_local.sing` print `DIMC 0` (replayed). `../wave2-excl-referee/ref_transfer.py PRIME file.sing`
   reduces the exact quadrics modulo a fresh prime and checks the pure-power criterion with msolve (`-g 2`).
7. **Replaying fixed-denominator exclusions without the engine.** `../wave2-excl-referee/ref_replay.py` (python-flint only) rebuilds the
   character-restricted Euler-normalized system `-D N(F) + F B_D(N) - F^2 Q = -A D^2` from `F`, `A` and a block expression and
   computes exact ranks over Q (fmpz_mat after column scaling) or modular ranks; its output reproduces the packet's matrix sizes,
   ranks and nullities exactly for `D0`, `X1*X2*FX1`, `FX0*X1*Y1`, `X1*D0*FX2` (degree 10, prime 1000000009), `X0*H*C` (E1, degree 11),
   and the one-dimensional degree-7 sectors (1,2), (2,2), (1,4), (2,4).
8. **Kernel-dimension table.** `kernel_dims_p2147483629.json` only holds d = 4..9; d = 10, 11 are `kernel-dims` records in the lane log.
   The report's "d <= 12" is d <= 11.
