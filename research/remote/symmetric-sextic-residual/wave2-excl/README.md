# Packet: wave-2 exclusion / rigor — `symmetric-sextic-residual` (q = 1, A = X0X1X2Y0Y1Y2)

Everything here is reproducible with the engine at `~/devel/sextics-work` (python-flint, numpy), `msolve` (≥ 0.6) and
`Singular`. The library is `../w2x_lib.py`; the job drivers `../w2x_relax_job.py`, `../w2x_e1kernel_job.py`,
`../w2x_reverify_job.py`, `../w2x_farkas_job.py`, `../w2x_sector_screen_job.py` (lane dir
`/tmp/hodge/results/symmetric-sextic-residual/`). One JSON receipt per experiment is appended to
`../wave2-excl.log.jsonl`; per-run receipts are the `*.json` files here.

## Exact statements certified (see the report `/tmp/hodge/notes/wave2-symmetric-sextic-residual-excl.md` §1 for the proofs)

* **E1 kernel.** For d ≤ 6 every D of degree d with `A D ∈ J(F)` lies in `J(F)_d` (so `J_5 ∪ J_6` are the only
  E1-admissible denominators of degree ≤ 6). Modular ranks (two primes) equal the complete-intersection values
  `dim J_n^{par} = #S_n^{par} − dim R_n^{par}`, `R` the Jacobian ring (`w2x_lib.dimJ`), which makes the modular
  dimension an upper bound for the rational one.
* **Hessian relaxation.** If (D, N, Q) is a certificate and `N⁰` is any polynomial vector with `Σ N⁰_i F_i = A D`, then
  `R⁰(D) = D·div N⁰ − N⁰·∇D ∈ J(F)`. Files `relax_<kind>_<chi>_d<d>_par<par>_p<prime>.ms` are the resulting quadratic
  systems (msolve format: first line variables, second line characteristic, then the quadrics) on a basis of the
  E1-admissible family modulo the radial subspace `F·S_{d−6}`; `*.lead` is msolve's reduced grevlex leading ideal
  (`msolve -f X.ms -o X.lead -g 1`). The exclusion statement is: *every variable has a pure power in `*.lead`*
  (then V = {0} over the algebraic closure of F_p, and by the cone/torsion-freeness argument no complex certificate
  denominator exists in that family). `*.sing` is the same system for Singular (`dim(std(I))`, `vdim`).
* **Farkas certificates** (`farkas_*.pkl`): for named denominators D, an integer vector `y` indexed by the
  `S3×S3`-orbits of Y-odd monomials of degree `2d+6` such that, unfolded to the full monomial space, `y·M = 0` for
  every column `M` of the full fixed-D linear system (unknowns: all coefficients of N_0..N_5 and Q, no symmetry
  assumption) and `y·b ≠ 0` for the right-hand side `−A D²`. This is an exact rational proof that the fixed-D system
  (with F scaled by 7 and D by its denominator lcm) has no solution over C. `verify_farkas.py` replays the check.

## Replay
```
cd /tmp/hodge/results/symmetric-sextic-residual
python3 w2x_relax_job.py --kind full --d 6 --par mixed --prime 2147483647 --threads 8 --singular   # ~2 min + msolve
python3 w2x_relax_job.py --kind full --d 7 --par 0 --prime 2147483647 --threads 16
python3 w2x_relax_job.py --kind sector --chi sgnsgn --d 8 --par 1 --prime 2147483647 --threads 8 --singular
python3 w2x_e1kernel_job.py --prime 2147483629 --degrees 1,2,3,4,5,6,7,8
python3 w2x_farkas_job.py D7            # exact certificate for D7 = Δ_X·e2(X0²,X1²,X2²)
python3 wave2-excl/verify_farkas.py wave2-excl/farkas_D7.pkl
```
