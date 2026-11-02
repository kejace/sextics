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

**Aggregator confirmation of the 2026-11-01 quirk (code reading, 2026-11-02 00:40 UTC):** confirmed. Both `FixedDSystem.solve()` (dense, `M = self.dense_augmented()`; `Mr, rank = M.rref()`) and `solve_compressed()` (`Mn = [compressed M | compressed b]`; `rank` from its rref) return the rank of the *augmented* matrix and `nullity = ncols − rank`. For an inconsistent system `rank = rank(M) + 1`, so the reported nullity is one below the true nullity of the coefficient matrix; for a consistent system the numbers are correct. `solve_fixed_D` (CRT path) and `e1_solvable` (which computes `rank` and `rank_aug` separately) are unaffected. Every "rank / nullity" pair quoted from the engine for an *inconsistent* system in the wave-1 reports (e.g. `1318 / 128` for a degree-5 polar; `2359 / 243` for a degree-6 Euler field) should be read as `rank(M) = 1317, nullity(M) = 129` resp. `2358 / 244`. The aggregator's own verification log stores the raw engine numbers; its verdicts (`consistent`) are unaffected. — aggregator

## [2026-11-01 ~19:30, wave2 excl quadric-contact-residual] msolve INPUT PITFALL (not an engine bug, but affects every lane's msolve verdicts)
msolve (the /usr/bin/msolve build on brunos and on the remote) MISPARSES a repeated-variable product written as `c1*c1`:
`30*c1*c1, 3*c0*c0` gives the reduced GB `[c1, c0]` (linear!), whereas `30*c1^2, 3*c0^2` gives `[c1^2, c0^2]` (correct).
Any .ms file that writes squares as `v*c(k)*c(k)` or `v*ck*ck` (wave-1 inv.py / verdict.py in quadric-contact-residual,
the dossier's d7_system.ms, and probably other lanes' exports) was solved on a WRONG (non-homogeneous) system; outputs like
`[-1]:` ("no solution") on a homogeneous system are a symptom.  Always emit squares as `ck^2`.  Singular parses `c*c` correctly.
Also: msolve's default solving mode crashes in FGLM for primes >= ~2^31 ("Prime ... is too large") -- for zero-dimensional
(only-origin) homogeneous systems use `msolve -g 1` (leading ideal) and check that every variable has a pure power among the
leading monomials (m-primary test); see /tmp/hodge/results/quadric-contact-residual/wave2-excl/msverdict.py.
