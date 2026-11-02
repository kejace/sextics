# README-referee — replaying the wave-2 exclusion packet (klein-coupled-character, q = 1, A = X0^3 Y0^3)

Written by the referee (Claude for `kejace`, 2026-11-02). The lens's own files are untouched; everything of the referee lives in
`referee/` and in this file.  Status of each claim: `referee/referee-report.md` (claim-by-claim verdict table, soundness notes, reconciliation with the verify lens).

## Pinned engine
python-flint 0.9.0, msolve 0.6.5 (`msolve -g 1` = leading ideal, `-g 2` = reduced GB; never write `c*c`, write `c^2`), Singular 4.3.x
(4.3.2 here; `LIB "primdec.lib"` is needed for `minAssGTZ`; the lens's `LIB "poly.lib"` lines fail harmlessly on this install).
Python 3.11+.  No network needed.  The instance is `data/klein-coupled-character/instance.rec` (7F has integer coefficients).

## Two ways to replay

**A. Referee scripts (self-contained; recommended for a stranger).**  `referee/ref_lib.py` hard-codes 7F and A (checked against the
engine parser), implements the sector bookkeeping (sigma weights (1,2,4,0,0,0), Y-parity), `(J:A)` ranks, the gauge-fixed Hessian
relaxation (exact over Q or mod p), and msolve/Singular writers.  Only python-flint, Singular and msolve are required.  Run from
`referee/`:

| script | claim(s) | output |
|---|---|---|
| `Singular -q regseq_1000003.sing` (and `_1000000009`) | flatness: Jacobian ideal m-primary, vdim 15625 | `regseq_*.out` |
| `python3 ref_colon.py` | T1 kernels d = 4..10 in all 14 sectors; `(J:l)_6 = J_6` for 8 linear forms | `referee_colon_p1000000009.json` |
| `python3 ref_hess21.py` | degree-5 closure input: exact rank 21 of the second partials | stdout |
| `python3 ref_relax.py d5|d6odd|d7empty` | relaxation m-primary (mod p1) | `relax/*_p1000000009.sing`, `relax_chain.log` |
| `python3 ref_relax.py d6even|d6mixed|d7lift` | exact (Q) radical memberships, curve C, lifted curves | `relax/*_exactQ.sing`, `relax_chain.log` |
| `python3 ref_quadrics.py` | second stage at d = 6 rebuilt mod p1; exact Singular replay of the 8 lens quadrics | `referee_quadrics_p1000000009.json`, `quadrics_*.sing` |
| `python3 ref_nzd.py` | NZD lemma, four charts at p1 | `nzd_*_p1000000009.sing`, `nzd_chain.log` |
| `python3 ref_fixedD_exact.py` | six polars and 36 `Z_j F_i` inconsistent exactly over Q | `referee_fixedD_exact.json` |
| `python3 ref_stage2_exact.py d par k` | exact constant-dual ("stage-2") quadrics for a sigma-eigen sector + m-primary test mod p1 | `stage2exact_*.json/.ms/.lead`, `stage2_odd0.log` |
| `python3 ref_hs_odd0.py <lead> <tag>` | Hilbert series of a msolve leading ideal vs the Segre cone over C x P^2 | `hs_*.sing` (note: Singular `hilb(.,1)` overflows for the odd-0 ideal; use `hilb(.,2)` or dim/mult) |
| `python3 ref_stage2.py d par k` | same as `ref_stage2_exact.py` but mod p1 only (exploratory; not a proof) | |

Expected headline numbers: kernels 0 for d <= 8; d = 9 even `[2,3,3,2,10,1,3]`, odd `[0,5,0,0,9,6,4]`; d = 10 even `[24,46,39,39,54,46,46]`,
odd `[36,39,39,36,39,36,36]`; d6 odd `dim 0 vdim 57`; d6 even `gb 201 dim 2`, 12 squares in I, C = 3 quadrics (dim 2 mult 8);
d6 mixed `gb 1354`, 30 squares; 8 necessary quadrics, `(Q, l0)` dim 0 vdim 13; d7 empty sectors vdim 4/16/8; d7 even k = 1,2,4 `gb 94 dim 2 mult 8`;
d7 odd k = 0 leading ideal `dim 4 mult 72` (plain relaxation) and `dim 2 mult 64` (with exact stage-2 duals): NOT closed.

**B. The lens's scripts.**  They import `w2lib.py`, which does `sys.path.insert(0, '/home/kejace/devel/sextics-work')` (the engine; in the
research branch it is `research/engine/`) and uses absolute `/tmp/hodge/...` paths (`LANE`, `PACKET`, `LOG` at the top of `w2lib.py`;
`load()` calls `hodge.core.load_instance`, which reads `data/<id>/instance.rec` relative to the repo root).  Edit those four constants and
put the repo root / engine on `PYTHONPATH`.  `w2_job{A,C,D}.py <prime> ...` are the remote jobs; the local pieces are `w2_radical.py`,
`w2_radical_exact.py`, `w2_hfsandwich.py`, `w2_jobC_local.py` and the Singular files under `local_*`.

## Claim -> receipt map (lens files)
* T1 / kernel tables: `jobC_local.out`, log `exp = reverify-colon-kernels` (p = 1000003).
* d = 5, 6 relaxations: `local_d6/d6_{odd,even,mixed}.{ms,lead,sing}`, `local_nzd/nzd.out` (mixed d = 5), log `exp = hessian-relaxation`.
* curve C and exact radical membership: `local_d6/rad_even.sing`, `rad_mixed.sing`, `exact_proof.out` (the final `AttributeError: fmpq_mat ... nullspace`
  is harmless: the exact duals were recomputed in `exact_duals.out`), `C_quadrics_exact.json`, `exact.out`.
* second stage: `local_d6/euler_family.out` (mod p), `euler_duals_exact.json` (150 integer duals), `euler_quadrics_exact.json` (8 quadrics), `exact_duals.out`.
* NZD lemmas: `local_nzd/*.sing`, `nzd.out`, `nzd_yquadrics.out`; coordinate forms in `local_d7/eigen.out`.
* d = 7 sectors: `local_d7/d7_{even,odd}_k*.{ms,lead,sing}`, `eigen.out`; exact containments `local_radical_exact/d7_even_k{1,2,4}.*` and
  `local_radical_exact_d7_even_k*.out`; mod-p containment for odd k = 0: `local_radical/d7_odd_k0.*`, `local_radical_d7_odd_k0.out`.
* Hilbert sandwich: `w2_hfsandwich.py`, `local_hfsandwich/d7_even_k{1,2,4}.out`.
* log: `../wave2-excl.log.jsonl` (138 lines).  The four `radical-membership-lifted-EXACT` lines with "exact GB size 1" (20:15:33-20:15:53) are VOID
  (msolve `+-` formatting bug), as recorded by the lens's own `log-correction` entry; the valid entries say "exact GB size 94".

## Gaps found by the referee
1. README T4 refers to a `duals/` directory that does not exist; the exact duals for polars / `Z_j F_i` are in the verify lens's files and
   are re-proved independently in `referee/referee_fixedD_exact.json`.
2. `relaxation-verdicts.md` predates the mod-p containment for (7, odd, 0) (log 20:32:17) and the referee's finding that the sandwich cannot
   close that sector (multiplicity 72 along a degree-24 family).  `local_hfsandwich/d7_odd_k0.sing` never produced output.
3. `local_d7full/d7_even_all.gb` and `local_radical8/d8_even_k0.lead` are empty: the non-eigen degree-7 and the degree-8 claims have no receipt.
4. `jobC_local.out`: only characters 0 and 1 of the degree-9 fixed-D tests finished (12 of 132 systems).
5. `local_hfsandwich/*.out` start with a `poly.lib` error; the numbers after it are valid (replayed: `HS_L = (1+t)^3 (1-t)^9`).
6. The scripts hard-code `p = 1000003`; the README's "prime as first argument" applies only to the remote jobs.

Synced to the research branch of kejace/sextics on 2026-11-02 (referee pass complete; (7, odd, 0) left open).
