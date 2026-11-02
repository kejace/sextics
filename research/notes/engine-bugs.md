# Engine bugs reported across lanes (wave 1 + wave 2) — aggregator's list with confirmation status

Updated: 2026-11-02 00:20 UTC. Engine: `/home/kejace/devel/sextics-work/hodge/{core,fixed_d,jacobian,symmetry,sweep,varpro,leangen,leancheck}.py`. Report new bugs by appending a dated bullet; the aggregator re-tests each one.

| # | module / function | symptom | reported by | status | aggregator confirmation |
|---|---|---|---|---|---|
| 1 | `fixed_d.FixedDSystem` columns and `fixed_d.exact_solution_from_pivots` | sign of the `F·N(D)` term wrong (identity `A D² − q D N(F) + F(D div N − N(D))`); every consistent/inconsistent/rank/nullity result before 2026-11-01 17:25 local unreliable | engine author (ENGINE-NOTICE.md); all lanes recomputed with `ENGINE_VERSION = 2` | **fixed** (v2 prints a banner on import) | Positive controls in the aggregator job (`cyclic-triple-residual:fermat_join_control`, `off-paired:fermat_control_f6_exact`, `rank-one-biquartic:split_control_f_exact`, `quadratic-cubic-pullback:fixedD_join012_target_delta_dense`) exercise the v2 identity end-to-end (exact lift + `core.verify_witness`); result pending. Smoke test (off-paired, local): E1/kernel/Hilbert paths OK. |
| 2 | lane script `symmetric-sextic-residual/served_scan.py` (not engine) | served-target vectors in `Sym(m)` coordinates compared with `J`/`A` vectors in `orbit_sum` coordinates (factor 36/|orbit|): spurious "2 classes beyond J served by every D" | symmetric-sextic-residual attack-1 agent (found and fixed: `check_served4.py`) | fixed in lane | not re-run by the aggregator (lane-local code); the corrected statement "the 2 classes are in J(F)" is consistent with the engine-level fact that all those D fail E1 / serve only J_6. |
| 3 | `fixed_d.solve_compressed` (design note, not a bug) | random row compression: a `consistent` verdict is probabilistic (false positive probability ≈ 1/p per random seed), an `inconsistent` verdict is rigorous mod p (R·M·x = R·b has a solution whenever M·x = b has one) | aggregator | by design | Documented here so that lanes do not quote compressed "consistent" results as proofs; `exact_solution_from_pivots` + `core.verify_witness` makes them rigorous. |
| 4 | `fixed_d.FixedDSystem.__init__` with `with_target_space=True` (used by `target_image`) | the Q-column parity uses `core.parity_of(A)` of the dummy target `A0` (first odd monomial) — correct for Y-odd targets only; `target_image` is only valid for Y-odd A (always the case for the challenge at q = 1 since deg A = 6q and A must be Y-odd) | aggregator (code reading) | not a bug for the challenge's targets | n/a |

Conventions reminder (from the engine docstrings): a mod-p E1-kernel dimension is an upper bound for the characteristic-0 dimension; `solve_fixed_D` lifts by CRT over `core.PRIMES` and verifies exactly; `sweep.py` treats E1 failure as `consistent=false` without building the full system.

## 2026-11-01 (wave2 excl, lane quadratic-cubic-pullback) — fixed_d reporting quirk (cosmetic, not a wrong answer)
`fixed_d.solve_compressed` / `solve_fixed_D_fast` report `rank` = rank of the AUGMENTED compressed matrix [M | b] and
`nullity = ncols - rank`.  For an INCONSISTENT system this is rank(M)+1, so the reported `nullity` is the true
nullity of M minus one (checked on D = dF/dX0, d = 5: reported 128, true 129 = 126 Euler elements (G z, 0) + 3 elements
(g F, D div g - g.grad D)).  Consistency verdicts are unaffected.  Workaround: when `consistent == False`, use
nullity+1 as the nullity of M.

**Aggregator confirmation of the 2026-11-01 quirk (code reading, 2026-11-02 00:28 UTC):** confirmed. Both `FixedDSystem.solve()` (dense, `M = self.dense_augmented()`; `Mr, rank = M.rref()`) and `solve_compressed()` (`Mn = [compressed M | compressed b]`; `rank` from its rref) return the rank of the *augmented* matrix and `nullity = ncols − rank`. For an inconsistent system `rank = rank(M) + 1`, so the reported nullity is one below the true nullity of the coefficient matrix; for a consistent system the numbers are correct. `solve_fixed_D` (CRT path) and `e1_solvable` (which computes `rank` and `rank_aug` separately) are unaffected. Every "rank / nullity" pair quoted from the engine for an *inconsistent* system in the wave-1 reports (e.g. `1318 / 128` for a degree-5 polar; `2359 / 243` for a degree-6 Euler field) should be read as `rank(M) = 1317, nullity(M) = 129` resp. `2358 / 244`. The aggregator's own verification log stores the raw engine numbers; its verdicts (`consistent`) are unaffected. — aggregator

## [2026-11-01 ~19:30, wave2 excl quadric-contact-residual] msolve INPUT PITFALL (not an engine bug, but affects every lane's msolve verdicts)
msolve (the /usr/bin/msolve build on brunos and on the remote) MISPARSES a repeated-variable product written as `c1*c1`:
`30*c1*c1, 3*c0*c0` gives the reduced GB `[c1, c0]` (linear!), whereas `30*c1^2, 3*c0^2` gives `[c1^2, c0^2]` (correct).
Any .ms file that writes squares as `v*c(k)*c(k)` or `v*ck*ck` (wave-1 inv.py / verdict.py in quadric-contact-residual,
the dossier's d7_system.ms, and probably other lanes' exports) was solved on a WRONG (non-homogeneous) system; outputs like
`[-1]:` ("no solution") on a homogeneous system are a symptom.  Always emit squares as `ck^2`.  Singular parses `c*c` correctly.
Also: msolve's default solving mode crashes in FGLM for primes >= ~2^31 ("Prime ... is too large") -- for zero-dimensional
(only-origin) homogeneous systems use `msolve -g 1` (leading ideal) and check that every variable has a pure power among the
leading monomials (m-primary test); see /tmp/hodge/results/quadric-contact-residual/wave2-excl/msverdict.py.

## quadric-contact-residual wave-2 agent (2026-11-01, Koszul varpro work)
- `hodge.varpro.VarProFamily`: a linearly DEPENDENT basis (e.g. raw Jacobian generators m*F_i) has directions with D(a) = 0, where b(a) = 0 and f(a) = 0 -> spurious exact minima. Only F-multiples are projected out. Workaround (lane-local w2lib.parameter_complement): also project out the null space of the family map a -> D(a).
- `hodge.varpro.VarProFamily` has no character-sector restriction of the N/Q unknowns (only Y-parity), so per-sector searches build the full system; for structured sectors this is ~50x larger than needed. Workaround: lane-local KoszulVarPro (sector-restricted Koszul formulation, Q-block eliminated by projection, `gelsy` driver ~2x faster than `gelsd` at equal cutoff).
- `hodge.rjob submit`: the rsync push of the lane dir can overwrite files that remote jobs are appending to (e.g. a shared wave2-search.log.jsonl). Remote jobs should write per-job files (we use HODGE_W2LOG=wave2-search.<job>.log.jsonl) and the puller merges them.

## rjob.up() overwrites remote notes (observed 2026-11-01 19:35 local, lane cyclic-triple-residual, excl lens)
`python3 -m hodge.rjob submit` calls `up()`, which rsyncs the LOCAL /tmp/hodge/notes/ onto the remote /tmp/hodge/notes/
without `--update`/`--delete`. Agents that run directly on the remote (Claude sessions on hodge-sextics write
/tmp/hodge/notes/wave2-<lane>-<lens>.md there) get their report files overwritten by stale local copies whenever a local
agent submits a job and the same file name exists locally. Workaround used here: check `md5sum` of the remote report before
submitting, and merge. Suggested fix: rsync notes with `--update` (or not at all) in `up()`.

**Aggregator confirmation of the msolve `c*c` pitfall (2026-11-02 00:31 UTC) and blast radius.** Reproduced locally with `/usr/bin/msolve` (`msolve -f in.ms -o out -g 2`): input `30*c1*c1, 3*c0*c0` gives the reduced Gröbner basis `[c1, c0]`; input `30*c1^2, 3*c0^2` gives `[c1^2, c0^2]`. So a repeated-variable product is parsed as the variable itself (every `v*v` becomes `v`), which turns a homogeneous quadratic system into a linear one and makes "no solution / only the origin" verdicts vacuous. Scan of every `*.ms` file under `/tmp/hodge/results` (regex `\b(\w+)\*\1\b`): **affected (71 files)** — `off-paired/` 57 (the wave-1 `msolve_d6_c*.ms`, `msolve_d7_c*.ms`, `msolve_d8_c00.ms`, `chk.ms`, `ctrl.ms`; 54 of the 55 wave-1 `.out` files say `[-1]` = no solution and are therefore **void**: the attack-1 §2.2 "degrees 6 and 7 re-derive fqlx's and deanodotdotdot's exclusions by a different method" is withdrawn — the thread results themselves are unaffected), `off-paired/wave2-excl/` 6 (the wave-2 excl lens's own `bilinear_test-local-d6/d6_c0*_sector_all_chart*.ms` charts — its table "d = 6, char [0,0]: 6 charts no-solution; [0,1]: 2 charts no-solution" is **also void** and must be re-run with `^2`), `quadric-contact-residual/` 5 (`quads_d8.ms`, `closure_quads_d8*.ms` — wave-1 verdicts already discarded by the lens), `quadric-contact-residual/wave2-excl/` 2 (its pitfall test files), `_aggregator/msolve-test/` 1 (this test). **Unaffected:** `cyclic-triple-residual/` (1 + 50 in wave2-excl), `reflection-coupled-residual/` (26 + 1), `dense-rational/` (4 + 3), `short-power-presentation/` 6, `eight-plane-configuration/wave2-excl/` 4, `diagonal-permutation/wave2-excl/` 2, `cycle-q2/wave2-excl/` 1, `symmetric-sextic-residual/wave2-excl/` 1, `off-paired/wave2-excl/` 2 — these write squares with `^` or have none. Singular is not affected. Recommendation for every lane: emit powers as `v^k`, and re-check any msolve "no solution" verdict whose `.ms` file matches the regex above. — aggregator

## varpro: `exclude_F_multiples=True` removes legitimate solutions (found 2026-11-01, lane one-reflection-skew)
`VarProFamily._complement_of_F_multiples` restricts the parameter space to the orthogonal complement of
`{a : D(a) in F*R_{d-6}}`.  That removes the spurious denominators `D = F*G`, but it also removes every solution
`D* = D0 + F*G0` whose decomposition has a nonzero F-multiple component — such a D* is generally NOT divisible by F
and is a perfectly good denominator.  Calibration on the accepted witness of `family-hesse-web` (reduced to the
5-variable problem): the true D lies in the E1 kernel (projection residual 8e-16) but not in the restricted span,
so random-restart searches in the restricted family can never find it.  Workaround: search the full family and
detect F | D afterwards (`is_spurious`); for the 5-variable strict problem there are no spurious solutions at all.
(The lane-local w2_search.py does this; engine untouched.)

## rjob submit hangs after launching the job (observed 2026-11-01 19:45 local)
Two `rjob submit` processes (excl-ctr-d8full, excl-ctr-relax78) stayed alive for 30+ min after the remote job had started
(the remote job ran normally; the local registry line was never written). The `ssh` call in submit() does not return —
probably because the detached remote process keeps the session's stdout/stderr open under heavy load. Consequence:
chained submits (`submit A; submit B`) never reach B. Workaround: wrap each submit in `timeout 600` and run them
independently; verify with `rjob status`.

**Aggregator confirmation of the `rjob.up()` notes clobber (2026-11-02 00:47 UTC):** confirmed from the source (`hodge/rjob.py`, `up()`): after pushing the engine and lane dirs it runs `rsync('/tmp/hodge/notes/', f'{USER}@{HOST}:/tmp/hodge/notes/')` with flags `-az` only (no `--update`, no `--ignore-existing`), and `submit()` always calls `up(lane)` first. Consequence: every `rjob submit` from brunos overwrites any same-named file in the remote `/tmp/hodge/notes/` with the local copy, even if the remote copy is newer. Remote-side agents already write under distinct names (`*-verify.md`, `*-excl-vm.md`, `wave2-theory-*.md`); the aggregator mirrors the remote notes into `/tmp/hodge/remote/_notes/` on every scan (22 remote-only report files at 00:45 UTC). The engine must not be edited by lane agents; the owner could add `--update` to that rsync. — aggregator

## reflection-coupled-residual wave-2 agent (2026-11-01 ~20:15) — limitations, not bugs
- `hodge.rjob pull` merges only `log.jsonl`; other per-lane JSONL logs (e.g. `wave2-search.log.jsonl`) must be merged by hand from `/tmp/hodge/remote/<lane>/` (dedupe by line).
- `hodge.varpro.VarProFamily` has no diagonal-character restriction of the unknowns and solves the full dense least squares per evaluation (≈30 s per evaluation at d=9); the lane-local `kvp.py` (Koszul-reduced, sector-restricted, sparse NF-mod-F + compression) is ~100× cheaper per evaluation for the same question.
- `fixed_d.target_image` builds the dense full fixed-D system; at d ≥ 8 use the compressed rows (`_sparse_matrix` + `_compress`) before `nullspace` (lane-local `job_served.py` does this).
- `jacobian.mult_kernel` works on Y-parity only; Klein-character sectors need the lane-local `rcr2.e1_space` (same algorithm with sector filters); results agree with the dossier/wave-1 numbers (K8 at d=8).
- zsh pitfall for agents: `set -- $spec` does not word-split in zsh; use `bash -c` for such loops.
- `hodge.rjob submit` exports OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=FLINT_NUM_THREADS=<cpus> for the job. With numpy/scipy dense linear algebra inside multiprocessing workers on the oversubscribed remote this is catastrophic: a 1002x130 lstsq took 0.3-0.4 s with 8 (or 40) BLAS threads vs 12 ms single-threaded (measured 2026-11-01, local and remote). Search scripts must set os.environ['OPENBLAS_NUM_THREADS']='1' (hard, not setdefault) BEFORE importing numpy; `hodge/__init__` only uses setdefault, so the rjob export wins. (quadric-contact-residual wave-2 agent; cost us ~1 h of remote search time.)
