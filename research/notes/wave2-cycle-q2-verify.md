# Wave 2 — lane `cycle-q2` — adversarial VERIFICATION (skeptical referee)

STATUS: in progress — 3 claims CONFIRMED, 0 REFUTED, 0 WEAKENED, 0 UNVERIFIED (updated 2026-11-02 00:35 UTC)

Log: `/tmp/hodge/results/cycle-q2/wave2-verify.log.jsonl` (one JSON line per check). Code: `/tmp/hodge/results/cycle-q2/verify/` (own library `vcheck.py`: engine-free, sympy parsing + python-flint linear algebra; Singular scripts in `verify/sing/`).

Conventions. "CONFIRMED" = recomputed independently (different prime and/or different tool) AND the logic transfers to characteristic zero as claimed. "WEAKENED" = the computation reproduces but the stated rigor level is too strong (e.g. single-prime inconsistency of an inhomogeneous system). "REFUTED" = wrong. "UNVERIFIED" = not yet checked / not checkable here.

## General rigor remarks (apply to many rows)

R1. `rank_{F_p} <= rank_Q` gives: zero kernel of an *integer matrix* mod p => zero kernel over Q. For kernels of maps into quotients by J(F) (E1 annihilators, "mod J" statements) one needs in addition that J has good reduction at p; this follows from smoothness of F mod p (Jacobian ideal zero-dimensional with vdim 15625 mod p), which the dossier verified at 10007 and I at 32003. With that extra step the dossier's E1 zero statements are rigorous.
R2. Inconsistency of an *inhomogeneous* system `M x = b` mod p does NOT by itself transfer to Q (counterexample: `p x = 1`). It transfers iff `rank_p M = rank_Q M`. A rigorous exclusion needs an exact Farkas dual (integer `y` with `yᵀM = 0`, `yᵀb ≠ 0`) or an exact rank over Z. The dossier's ledger ("each a zero-kernel statement mod p ⇒ over ℂ") is therefore too strong for the fixed-D exclusions of A.5/A.6 and for every "rank/aug = r/r+1" result in the thread (except fqlx's block-radial k=1 with certified integer duals). Below I recompute exact ranks over Z with flint `fmpz_mat.rank` where feasible.
R3. Row compression by a random matrix R can only turn "inconsistent" into "consistent" (false positive), never the reverse; and a Farkas dual of `R[M|b]` pulls back to one of `[M|b]`. So exact inconsistency of a compressed system is a valid proof.

## Claim table

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
| 1 | F smooth (vdim 15625) | dossier §1.4 | CONFIRMED | Singular p=32003 | transfers to Q |
| 2 | Ann_d([A]) = 0 (d≤5), 26 (d=6), 420, 825 | dossier §1.9 | CONFIRMED | Singular p=32003, identical | zero rows rigorous via R1; 26/420/825 are upper bounds for Q-dims |
| 3 | Ann_d = 0 (d≤7) for cycle-three target and siblings B, B_bal, C1; d=8 values 2/3/6/4 | dossier §1.9/A.3 | CONFIRMED | Singular p=32003 | |
