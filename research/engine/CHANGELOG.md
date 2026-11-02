# hodge engine changelog

Engine modules live in `hodge/`; `fixed_d.ENGINE_VERSION` is bumped only when a numerical result
of an existing function changes (then every log line records `"engine": N`).  Pure additions and
speedups keep the version.

## 2026-11-02 (maintainer session on the VM)

- rjob: runs locally when invoked on the compute box itself (hostname `hodge-sextics`, or
  `HODGE_REMOTE_LOCAL=1`); `up`/`pull` become no-ops, `submit` still detaches with nohup/setsid and
  records the pid in `/tmp/hodge/jobs/<name>.pid`.  Before this the ssh-to-self failed on the VM
  (no `~/.ssh/google_compute_engine` there).
- tests: baseline `tests/run_all.py` ALL OK (jacobian 10 s, symmetry 1 s, varpro 104 s, sweep 3 s).
- engine v3 (fixed_d.ENGINE_VERSION = 3): `FixedDSystem.solve()` and `solve_compressed()` now report the rank and
  nullity of the COEFFICIENT matrix for inconsistent systems (v2 reported rank(M)+1 / nullity-1, see
  engine-bugs.md); both also return `rank_aug`.  Consistent-case numbers unchanged.
- new `hodge/sge.py` + `hodge/csrc/sge.c`: sparse modular Gaussian elimination (Markowitz-style pivoting,
  lazy column lists, exact counts in buckets) with a dense flint finish and back-substitution; verified
  against flint on random systems (rank, consistency, pivots, nullspace) and on accepted witnesses.
  `fixed_d.solve_sparse(S)` / `sge.solve_fixed_system(S)`; `solve_fixed_D_fast(method='auto')` uses it by
  default (compressed rows only as fallback).  Measured (load > 300): hesse-triangle d=6 11.6k x 2.6k
  2.4 s vs 130 s compressed; cyclic-triple d=8 40k x 7.5k 36 s vs ~900 s.
- `fixed_d.set_threads(n)` / `sge.set_threads(n)`: FLINT threads for dense rref (FLINT_NUM_THREADS env,
  which `rjob submit --cpus N` sets); 4x at 8 threads for n=4000.
- `fixed_d._compress` is chunked (memory ~ k*4096*8 bytes instead of k*nrows*8: 54 GB at d=12 before).
- new `hodge/parallel.py`: `solve_fixed_D_parallel(F, A, q, D, workers=4, method='sparse')` solves several
  primes concurrently (fork) and CRT-lifts; same return as `fixed_d.solve_fixed_D`.
- sweep records: `kind` (fixedD | e1 | skipped), `verify` (True only after core.verify_witness), `control`
  (--control flag or candidate field), `solver` = 'fast:sparse' / 'fast:compressed' / ...
- VM-side edits live as idempotent `tools/patch_*.py` (run `tools/apply_patches.sh`) because brunos-side
  `rjob submit` still rsyncs its engine copy over this one.

## 2026-11-01

- engine v2: sign fix in the `F*N(D)` term (FixedDSystem columns and exact_solution_from_pivots);
  see `/tmp/hodge/results/ENGINE-NOTICE.md`.
- modules: core, fixed_d (dense + compressed solvers, target_image, e1_solvable), jacobian,
  symmetry (diagonal group sectors), sweep (multiprocessing candidate sweep CLI), varpro (numerical
  variable-projection family search), leangen / leancheck / sandboxcheck (witness -> Lean pipeline), rjob.
