# Wave 2 — lane `cycle-q2` — lens EXCLUSION / RIGOR

STATUS: no witness (exclusion lens); in progress. Best rigorous bound so far: **deg D ≥ 6** (new; was ≥ 5 with d=5 only partially excluded).

Started 2026-11-01 ~19:05 local. Log: `/tmp/hodge/results/cycle-q2/wave2-excl.log.jsonl`.
Packet: `/tmp/hodge/results/cycle-q2/wave2-excl/` (scripts `w2lib.py`, `reverify.py`, `hess5.py`, `hess6.py`, `ann.sing`; receipts `*.receipt.json`, certificates `certs/`).
Primes used in wave 2: P1 = 1000000007, P2 = 998244353 (wave 1 used 2^31-1, 4194301, 4194287; dossier 10007, 65521). Engine v2 (post sign-fix).

## Results so far

### R1. Engine / solver cross-validation (local, 19:12)
* `core.verify_witness` (= the contract identity `E = A D² − q D Σ N_i F_i + F(D Σ ∂_i N_i − Σ N_i ∂_i D) = F^{q+1} Q`, checked against `formal/README.md`) accepts the accepted witnesses of `fermat-q2` and `block-q2`.
* engine-v2 `symmetry.restricted_system` and wave-1 `eq.py` build *identical* column sets and give identical ranks (tested on fermat-q2, block-q2 witnesses: consistent; on `D = F_X0`: inconsistent, rank(M)=442, rank([M|b])=443, 540 cols, at p = 10^9+7, 2^31−1, 10007). NOTE: the engine's reported `rank` for an inconsistent system is the *augmented* rank (443), so "rank 443 / nullity 97" in engine output = "442/443" in the dossier's notation. No real discrepancy.

### R2. Exact Farkas certificates (new rigor tool)
For an inconsistent sector system `M x = b` we compute an exact rational vector `y` supported on rank+1 rows with `yᵀM = 0`, `yᵀb ≠ 0`, verified with exact integer arithmetic (`w2lib.farkas`). This is a prime-free proof that no `(N,Q)` exists over any field of characteristic 0 (the sector restriction is lossless by averaging). Smoke test: `D = F_X0`: certificate with 443 rows, entries ≤ 166 bits, verified, 1.2 s. (Modular inconsistency alone transfers to char 0 only when `rank_p(M) = rank_Q(M)`, which a single-prime computation does not certify; the Farkas vector removes that caveat.)

### R3. Degree 5 excluded (rigorous) — `hess5.py`, receipt `hess5.receipt.json`
Chain: E1 ⇒ `A·D ∈ J`; `Ann_5([A]) = 0` (dossier A.3, re-verification by `ann.sing` pending) ⇒ `D = ∂_w F = Σ w_i F_i`, `w ∈ ℂ⁶`. With `N⁰ = (A/2) w` E1 holds exactly, and E2 becomes `H(w) + (D div K − K·∇D) = 2F²Q` with `K` Koszul (the Euler-gauge term cancels identically), where `H(w) = (∂_w F)(∂_w A) − A ∂_w²F` (degree 16). Hence `H(w) ∈ J_16` is necessary. `H(w) = Σ_{i≤j} w_i w_j H_ij`; the 21 polynomials `H_ij` lie in 21 *distinct* G-character sectors of `S_16`, and each is nonzero modulo `J_16` (rank of J-generators + H_ij = rank J + 1 in every sector, at both primes). So `[H_ij]` are independent in `R_16`, forcing `w_i w_j = 0` ∀ i,j, i.e. `w = 0`. Transfer to char 0: `dim_p J_16 = 19398 = dim_Q J_16` (regular sequence; both primes), so `rank_p ≤ rank_Q` gives independence over ℚ. **No certificate with deg D = 5 exists.** (Cross-check: `hess6.py --deg 5` independently reproduces K_5 = J_5 (dim 6) and 264 nonzero quadrics covering all 21 monomials `c_i c_j`, i.e. the ideal is `(c)²`.)

### R4. Re-verification sweep (`reverify.py`, remote job `excl-cycle-q2-reverify`, running)
190 candidates (dossier A.5 + wave-1 structural families + all `Z_j F_i`, all `F_i F_j`, squares at d=10/12): engine-v2 sector solve at P1 and P2 + exact Farkas certificate for each inconsistent case. Results pending (remote at load ~290/128, 33 jobs).

### Pending
* `ann.sing` (E1 annihilator dims d=0..9 at two primes, smoothness vdim, q=1 siblings) — waiting for remote capacity.
* `hess6.py --deg 6` (62-dim E1 kernel, 426 quadrics) + msolve/Singular Gröbner — waiting for capacity.
* sector-wise non-Jacobian E1-kernel families at d = 6,7,8 (pencils via relaxation + full tests), structured families d ≤ 12.
