# README-referee — replaying the wave-2 exclusion packet for `cyclic-triple-residual`

Written by the referee pass (2026-11-02). The agent's files in this directory are unchanged; this file adds what a
stranger needs to replay them, records which receipts support which claim, and points to the independent replays in
`../wave2-referee/` (own code, no engine imports, fresh prime 999999937, exact over Q where feasible).

## 0. Target and conventions

`F = X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6 + (5/7)(X0^3X1Y0Y1 + X0^3X2^3 + X0^2Y0Y1Y2^2 + X0X1Y1^3Y2 + X0X2^3Y1Y2 + X1X2^3Y0Y1)`,
`q = 1`, `A = X0^3 Y0^3`, variables `X0,X1,X2,Y0,Y1,Y2`. All scripts work with `Fz = 7F` (integer); `J(Fz) = J(F)`.
Diagonal group `G` of order 18, generators `mu3 = diag(1,1,z3,1,1,1)` and `s = (w,w^5,w,w^2,w^2,1)`, `w = exp(2 pi i/6)`;
a monomial `m` has character `chi(m) = (m.(0,0,2,0,0,0) mod 6, m.(1,5,1,2,2,0) mod 6)`; `chi(A) = (0,3)`.
Certificate identity (formal/README.md): `E = A D^2 - D sum N_i F_i + F (D div N - N(D)) = F^2 Q`; its split form at
`q = 1` is E1 `sum N_i F_i = A D`, E2 `D div N - N(D) = F Q` (any certificate can be shifted by `N -> N + R Z/6` to
satisfy E1 exactly because `D` is coprime to the irreducible `F`).

## 1. Environment

* python 3 with python-flint 0.9.0 (`pip install python-flint==0.9.0`); msolve 0.6.x (binary `msolve`, options used:
  `-g 1` leading ideal, `-g 2` reduced GB, `-t T` threads); Singular 4.3.2 (`Singular -q file.sing`).
* Hard-coded paths in the agent's scripts: `/home/kejace/devel/sextics-work` (the `hodge` engine, published as
  `research/engine/`), `/tmp/hodge/results/cyclic-triple-residual` (the lane directory; `sect.py` there is the wave-1
  sector solver imported by `w2lib.py`), `/tmp/hodge/results/cyclic-triple-residual/wave2-excl` (this packet).
  To replay from a checkout of the `research` branch:
  ```sh
  R=$PWD/research                                   # checkout of kejace/sextics, branch research
  mkdir -p /tmp/hodge/results /home/kejace/devel
  ln -s $R/results/cyclic-triple-residual /tmp/hodge/results/cyclic-triple-residual
  ln -s $R/engine /home/kejace/devel/sextics-work
  ```
  The scripts append receipts to `/tmp/hodge/results/cyclic-triple-residual/wave2-excl.log.jsonl`.
* The independent replays (`../wave2-referee/ref_*.py`) need only python-flint, msolve and Singular; they hard-code
  only their own directory and `../wave2-referee/mirror/` (copies of the remote-only receipts, see section 3).

## 2. Which receipt supports which claim

| claim | receipt(s) in this packet / lane | independent replay (`../wave2-referee/`) |
|---|---|---|
| E1 kernel `Ann([A])_d mod J_d = 0` for `d <= 7`, `= 1` at `d = 8`, `120, 556, 1020` at `d = 9,10,11`; `dim J_n mod p` = c.i. value | `kerverify.py` -> `../wave2-excl.log.jsonl` (p = 1000003, 32003); `mirror/kerverify.out`, `mirror/a1_exact_8_12.out` (exact, d = 8, 9) | `ref_e1.py` -> `e1_table.json`: exact over Q for `d <= 8`, mod 999999937 for `d <= 11` |
| `D8` (25 terms, sector (0,4)): `D8 not in J_8`, `A D8 in J_14`, explicit `P` with `sum P_i d_i(7F) = A D8` | `mirror/e1_d8_exact.json` (remote-only file; the README's `relax_run.py 8` needs it) | `ref_d8.py`: all three facts exact; own exact nullspace gives `Ann_8 = 1` with `D8` in its span |
| `D8` excluded over C | `d8full2.py` -> log line `D8-sector-sealed-reverify` (sector system 2357 x 1188, nullity 33 sealed, 1155/1156) and `D8-sector-exact-fmpq-rref`; `mirror/receipts/farkas_D8.json` (exact Farkas dual, 968 nonzeros, single identity, 1041/1042 on 1188 cols) | `ref_d8.py`: own single-identity G-sector system, exact fmpz ranks; Farkas dual replayed against all 15015 unrestricted columns (`bad = 0`, `y.rhs = 1`) |
| Hessian relaxation: `d = 5` m-primary (no degree-5 denominator) | `relax_d5_p1000003*.json`, `relax0_d5_p1000003_full_gauge.json`; verify lens: Singular char 0 | `ref_hess.py d5`: exact over Q (Singular char 0, `dim 0, vdim 7`; the 21 quadrics span all quadratic forms) |
| `d = 6`, all parities/characters, gauge-fixed m-primary (no degree-6 denominator) | `relax_d6_p1000003_g_gauge.json` (`pure_powers 35/35`), `relax_d6_p1000003_g_gauge_lead.txt`, `relax6_local.err` (`dimR_total = 1751`, i.e. `rank_p J_12` = c.i. value) | `ref_hess.py d6`: p = 999999937, both the X0- and the Y2-slice, Singular `dim 0` and msolve 35/35 pure powers |
| `d = 7`, G-semi-invariant `D`, per character sector m-primary | `hsector_d7_*.local.out`, `relax_d7_p1000003_s*_gauge.json` (13 sectors (0,0)-(0,5), (2,0)-(2,5), (4,0) finished by the agent) | `ref_hess.py d7`: all 18 sectors, p = 999999937, X0-slice, Singular `dim` + msolve |
| Fermat control (`t = 0`): the same d = 6 pipeline is NOT m-primary | `relax_d6_p1000003_fermat_gauge.json` (`pure_powers 12/35`) | not replayed |
| Degree-9 kernel elements, each tested with the sector system | `mirror/kernel_exact_d9.json` (120 exact elements), `mirror/kexact9.remote.out`, `mirror/wave2-excl.log.remote.jsonl` (`e1-kernel-element-sealed-test`: 116 modular-only at p = 1000003, 2147483587; 4 exact by the fmpz fallback) | not replayed (modular evidence only for 116/120) |

## 3. Gaps found (and how they are handled)

* Remote-only files referenced by the agent's README/notes (`e1_d8_exact.json`, `kernel_exact_d9.json`, `kexact9.out`,
  `relax7.out`, the parallel session's `lib.py`, `b2_hessian.py`, `a2_farkas.py`, `receipts/`): copied to
  `../wave2-referee/mirror/` (the `/tmp/hodge/remote` mirror is not published by `sync.sh`).
* `../wave2-excl.log.jsonl` in this lane directory lacks the 120 `e1-kernel-element-sealed-test` lines written on the
  remote; the complete remote log is `../wave2-referee/mirror/wave2-excl.log.remote.jsonl`.
* Not finished at referee time (no result files exist): full `d = 7` relaxation (126 unknowns, 1441 quadrics; msolve GB
  running > 2 h on the loaded remote), full `d = 8` (337 unknowns), the degree-8 sector chain (`hsector_d8_*`),
  `families.py` (no `families.out`), the parallel session's Y-even `d = 7` msolve run. The draft post's sections 3-4
  (degree-9/10 "sector families", "structured candidates ... all sealed") have no receipts in this packet.
* The c.i. check `dim J_n mod p` was logged for `n <= 11` only; for the relaxation at degree `d` the check needed is
  `n = 2d`. It is implied by `dimR_total = 1751` in `relax6_local.err` (d = 6) and by the sector counts (d = 7), and is
  re-done explicitly in every `ref_hess.py` receipt (`rankJ_2d == ci_dimJ_2d`).
* Numbers of "quadrics" differ between scripts (raw coordinate functionals vs. a basis of their span): agent 50/675,
  referee 21/301 independent quadrics at d = 5/6. The ideals coincide (same construction); the leading-ideal sizes agree
  with the aggregator's independent run (1137 with the Y2-slice).

* Self-audit of the referee's quadric construction (`referee.log.jsonl`, `self_audit_relaxation_d6`): the Euler line
  lies in `V(I)`, the quadrics are invariant under translation along it, and their values at a random point agree
  with the direct normal form of `A G^T Hess(F) G` (1472 nonzero coordinates = the raw coordinate count).
* Full `d = 7` system built by the referee for the record only (`hessian_relaxation_d7_full_build_only`): 126 unknowns,
  921 independent quadrics (120 unknowns on the X0-slice), `rank_p J_14 = 10122` = c.i.; msolve input saved as
  `../wave2-referee/hess_d7_p999999937_X0_full.ms`; no Groebner basis was run.
* `../wave2-referee/notes-snapshot/` holds copies of the exclusion report, the draft post and the verify notes as they
  were at referee time (the canonical copies live in `research/notes/`).

## 4. Replay commands (referee scripts; each appends to `../wave2-referee/referee.log.jsonl`)

```sh
cd research/results/cyclic-triple-residual/wave2-referee
python3 ref_e1.py            # E1 table, exact d <= 8, mod 999999937 d <= 11      (~5 min)
python3 ref_d8.py            # D8 facts, own fixed-D exact ranks, Farkas replays   (~3 min)
python3 ref_hess.py d5       # exact over Q + mod p                                (seconds)
python3 ref_hess.py d6       # Singular dim + msolve, both slices                  (~1.5 min)
python3 ref_hess.py d7       # 18 character sectors                                (~10 min)
python3 ref_hess.py d7yeven  # Y-even degree-7 family (60 unknowns), msolve        (see referee-report.md)
python3 ref_hess.py d7yodd   # Y-odd degree-7 family (60 unknowns), msolve
```
