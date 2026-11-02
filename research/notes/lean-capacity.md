STATUS: in progress (restarted 2026-11-01 ~19:00 EST) — baseline + synthetic timing jobs being prepared; no timings yet.

# Lean capacity study (accepted template vs alternatives)

Work dir: /tmp/hodge/lean/capacity (log.jsonl, gen/, src/Cap/, results/). Remote jobs: lane `_cross`, names `cap-*`.

## Verifier budget (scripts/proposals/sandbox.py + .github/solve-policy.json)
- Each module: `lean -j 2 -M 14336` in Docker `--cpus 2`, then `leanchecker --verbose` per module; total wall budget 3600 s for everything (lake prep, all compiles, all replays, audit); 16 MiB source total, 64 files.
- Banned tokens: axiom sorry admit unsafe native_decide run_cmd run_tac initialize elab macro syntax #eval #exit [extern. `decide` (incl. `decide +kernel`) is allowed; axioms limited to propext/Classical.choice/Quot.sound.

## Known floor (accepted proofs that passed CI)
| instance | q | d | |F| | |D| | |N_i| | e1 raw products | e2 raw products |
|---|---|---|---|---|---|---|---|
| m3-coupled | 1 | 10 | 22 | 190 | 88,87,83,35,83,35 | 2877 | 78310 |
| contact46006 | 1 | 6 | 56 | 44 | 48,84,84,48,84,84 | 11280 | 18384 |
(raw products = number of monomial×monomial products before collection, the quantity `ring` must normalize.)

## Timings
(to be filled)
