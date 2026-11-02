# Wave 3 — degree-7/8 closure coordinator (running table)

Started 2026-11-02 ~04:00 UTC on brunos; Gröbner jobs run on hodge-sextics (128 cores / 503 GB) through `hodge.rjob`
(lane `_closure`, two jobs at a time, `--cpus 16`, `HODGE_JOB_MEM_GB=140`, `msolve -g 1 -t 16`, 6 h timeout per system).
Log: `/tmp/hodge/results/_closure/log.jsonl`. Receipts: `/tmp/hodge/results/_closure/wave3-receipts.jsonl` (remote,
pulled) and `/tmp/hodge/results/<lane>/wave2-excl/wave3/` (per lane). Runner: `_closure/wave3_run.sh` + `wave3_verdict.py`.

## 0. Rigor convention (stated once, used by every "m-primary" row)

Each input is the Hessian-relaxation system of a lane: unknowns = coordinates of the polar field `w` of a Jacobian
denominator `D = w(F)` in the stated parity/character sector, equations = the quadrics expressing
`A · Hess_F(w, w) ∈ J(F)` (mod-J reduction), reduced mod the stated prime `p` with good reduction recorded by the lens
(`good_reduction: true` in the `.json` sidecar: the Jacobian ranks used in the reduction agree mod `p` and over ℚ).
Let `I_ℚ` be the rational ideal and `I_p` its reduction. For every degree `k`,
`rank_{F_p} Mac_k(I_p) ≤ rank_ℚ Mac_k(I_ℚ)`, hence `HF(ℚ[w]/I_ℚ, k) ≤ HF(F_p[w]/I_p, k)`. If the grevlex leading ideal of
`I_p` contains a pure power of every variable, `I_p` is m-primary, its Hilbert function vanishes for `k ≥ k₀`, so the
same holds for `I_ℚ`: `V(I_ℚ)(ℂ) = {0}`, i.e. no nonzero polar field, i.e. **no denominator of that degree in that
sector over ℂ**. (Chart/dehomogenised GB = {1} mod p does NOT transfer; only homogeneous leading-ideal certificates are
used here.) Combined with the lens's E1 statement "every E1-admissible `D` of this degree is Jacobian, `D = w(F)`"
(verified exactly in the wave-2 reports), an m-primary verdict excludes the degree in that sector outright.

## 1. Inventory of OPEN degree-7/8 Gröbner systems (from wave2-status §6–7, the 16 excl reports, 10 referee reports)

| # | lane | degree / sector | unknowns / quadrics | prime | input (remote, `/tmp/hodge/results/<lane>/`) | wave-2 fate |
|---|---|---|---|---|---|---|
| 1 | reflection-coupled-residual | d7, Klein character (0,0) | 42 / 274 (426 raw) | 2147483587 | `wave2-excl/hess_d7_W00_p2147483587.ms` | `.gb` empty (killed 02:02 UTC) |
| 2 | reflection-coupled-residual | d7 mixed pairs (01,11), (00,11), (01,10), (00,10), (00,01) | 56/470, 60/468, 60/479, 64/463, 80/? | 2147483587 | `wave2-excl/hess_d7_W{0111,0011,0110,0010,0001}_*.ms` | empty `.gb` / never run; (10,11) done (`W1011.gb`) |
| 3 | klein-coupled-character | d7 non-eigen full, Y-even / Y-odd | 60 / 447, 60 / 420 (sparse) | 1000003 | `wave2-excl/local_d7full/d7_{even,odd}_all.ms` | `.gb` empty (never finished locally) |
| 4 | symmetric-sextic-residual | d7 full pure parity, par0 / par1 | 61 / 741 each | 2147483647 | `wave2-excl/relax_full_all_d7_par{0,1}_p2147483647.ms` | `.lead` empty, log header only (killed) |
| 5 | eight-plane-configuration | d7 pure parity par0 / par1 | 60 / 483 each (dreg 4) | 2147483587 | `wave2-excl/hess_d7_par{0,1}_p2147483587.ms` | OOM rc 137 at ~40 GB; `hess7b` resubmission never launched |
| 6 | dense-rational | d7 τ-even / τ-odd | 60 / 483 each | 1073741789 | `wave2-excl/quadsys_d7_{even,odd}_p1073741789.ms` | even killed in F4 degree 4; odd guard-killed 02:54 UTC at 38 GB |
| 7 | diagonal-permutation | d7 even (odd running elsewhere, `-g 2 -t 4`) | 60 / 480 | 1000000007 | `wave2-excl/relax_d7e.ms` | guard-killed 02:53 UTC at 39 GB |
| 8 | quadratic-cubic-pullback | d7 Y-even polar cone (63 vars incl. gauge) / 60-var versions both parities | 63 / 483; 60 / 483 | 268435399; 2147483629 | `wave2-excl-hess-hessL7-d7-par0-p268435399-KH.chart1.ms`; `w2_points_d7_p{0,1}_c2147483629.ms` | never finished (`.out` empty) |
| 9 | quadric-contact-residual | d7 diagonal planes `y0`, triple `000_011_101` (others running under `excl-planes7`) | 60 / 477; 48 / 451 | 2147483629 | `wave2-excl/sub_d7_y0.ms`, `sub_d7_tri_000_011_101.ms` | timed out (5406 s / 7200 s caps) |
| 10 | cyclic-triple-residual | d7 semi-invariant sectors s41–s45 (13/18 done); full d7 (126/1441) and d8 s00 running elsewhere | tiny (3–7 vars) | 999999937 | `wave2-referee/hess_d7_p999999937_X0_s4{1..5}.ms` | not finished locally |
| 11 | reflection-coupled-residual | d8: K8 family (1+102 unknowns), characters (1,1), (1,0) [(1,0) running elsewhere `-g 2`] | 103 / ?, ? | 2147483587 | `wave2-excl/hessK_d8_chi00_*.ms`, `hess_d8_W11_*.ms` | `.gb` empty |
| 12 | klein-coupled-character | d8 even k=0 eigen sector | 60-ish / 72 | 1000003 | `wave2-excl/local_radical8/d8_even_k0.ms` | `.lead` empty |

Running by other agents at 04:00 UTC (not duplicated here; results collected when they land): quadric-contact
`sub_d7_pl_even` (twice), `sub_d7_y1`, `d8_chi1100`; short-power `excl_residual_d7e_p0` (msolve + Singular),
`excl_blocks_d7o_comp`, `excl_mixed_d7_x`; cyclic-triple `relax_d7_p1000003_gauge` (126 unknowns), `relax_d8_s00`;
diagonal-permutation `relax_d7o`; reflection-coupled `hess_d8_W10`; one-reflection-skew `d7_sub_tor_1110`; dense-rational
`relax_d7_p0_q`. Out of Gröbner reach (not queued): mixed-parity d7 on dense lanes (120 vars / ~950 quadrics, dreg 5),
d8 pure cones (~160 vars / 468 quadrics), quadric-contact full mixing; polarized-k3 d7 (relaxation itself too weak).

## 2. Queue order (smallest first) and status (04:35 UTC)

Policy note: `rjob submit --force` is denied by the permission classifier ("Interfere With Workloads"), so every job
goes through the rjob capacity gate (`--wait 3600`, retried); the gate needs ≥ 110 GB available and ≤ 200 allocated
cpus (the box had 60 foreign jobs / 361 cpus / 66–97 GB available all night). One tier-1 chain (`wave3-B1`) got through;
the rest of the queue is a single sequential runner job (`wave3-Q`, `wave3_queue.sh` over `queue.txt`, one Gröbner
computation at a time, per-entry memory floor so that a 40–60 GB F4 step never drives the box under the 45 GB guard).

| job | systems | status |
|---|---|---|
| wave3-B1 (running since 04:07:52 UTC, pid 3210118, 16 cores) | klein d7 even (60/447) → klein d7 odd (60/420) → ctr s42–s45 | running |
| wave3-Q (gated submit pending) | rcr W00 → W0111 → W0011 → W0110 → ssr par0 → par1 → epc par0 → par1 → dense even → odd → dp even → qcr y0 → qcr tri → rcr W0010 → W0001 → rcr d8 W11 → K8 → klein d8 k0 → qcp Y-even/Y-odd (homogenised) | waiting for the gate |

## 3. Results (rewritten after every job)

### 3a. New wave-3 computations
(none finished yet; B1 running)

### 3b. Harvested wave-2 outputs re-parsed independently (`wave3_harvest.py`, tags harvest1/harvest2; receipts in
`/tmp/hodge/results/_closure/wave3-receipts.jsonl` and `<lane>/.../wave3/receipts.jsonl`)
All of the following finished wave-2 msolve outputs (reduced GB or leading ideal) have an m-primary leading ideal, so
by §0 the corresponding sector has no denominator of that degree over ℂ (given the lens's E1 ⇒ Jacobian statement):

| lane | degree / sector | unknowns / quadrics | prime | files |
|---|---|---|---|---|
| reflection-coupled-residual | d7 Klein characters (0,1), (1,0), (1,1) and pair (1,0)+(1,1) | 38/426, 22/426, 18/426, 40/879 (raw) | 2147483587 | `wave2-excl/hess_d7_W{01,10,11,1011}_*.gb` (complete reduced GBs, max pure power 4) |
| symmetric-sextic-residual | d7 S₃×S₃-sectors triv, sgnX, sgnY, both parities (sgn⊗sgn sectors empty) | 7/46, 3/46, 1/46, 1/20, 3/41, 7/46 | 2147483647 | `wave2-excl/relax_sector_*_d7_par*.lead` |
| cyclic-triple-residual | d7 G-semi-invariant sectors s00–s05, s20–s25, s40 (gauge-fixed) + s41–s45 (referee files, trivially all quadratic monomials) | 3–16 unknowns | 1000003 / 999999937 | `wave2-excl/relax_d7_p1000003_s*_gauge_lead.txt`, `wave2-referee/.../wave3/*.lead` |
| eight-plane-configuration | d7 subfamilies arr / blockX / blockY / coefB0 / coefX, both parities; d8 subfamilies coefB0 / coefX / comp-0..5, both parities | 9–36 unknowns | 2147483587 / 2147483629 | `wave2-excl/hsub_d{7,8}_par*_*.lm` |
| quadric-contact-residual | d7: all 28 (σ0,σ2,τ_Y)-parity-class pairs at two primes (56 files) | 24–32 / 227–258 | 2147483629 and p3 | `wave2-excl/sub_d7_*.ms.lead`, `sub_d7p3_*.ms.lead` |
| klein-coupled-character | d7 even k=0 sector | 3/15 | 1000003 | `wave2-excl/local_d7/d7_even_k0.lead` |

Not decided by this check (leading ideal positive-dimensional; the lenses closed these with second-stage exact
arguments, not re-derived here): one-reflection-skew `d7_sub_tor_{1000,0100,0010}`, `d7_5v_tor_*`; klein `local_d7`
k1, k2 sectors.

## 4. Rigorous degree bounds per lane (running)

(unchanged from wave2-status.md until a row lands in §3)
