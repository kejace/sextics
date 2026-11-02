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

## 2026-11-01

- engine v2: sign fix in the `F*N(D)` term (FixedDSystem columns and exact_solution_from_pivots);
  see `/tmp/hodge/results/ENGINE-NOTICE.md`.
- modules: core, fixed_d (dense + compressed solvers, target_image, e1_solvable), jacobian,
  symmetry (diagonal group sectors), sweep (multiprocessing candidate sweep CLI), varpro (numerical
  variable-projection family search), leangen / leancheck / sandboxcheck (witness -> Lean pipeline), rjob.
