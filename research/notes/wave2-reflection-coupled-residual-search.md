# Wave-2 joint search — lane `reflection-coupled-residual` (q=1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (search running; report rewritten after every batch). Started 2026-11-01 19:05 local.

Log: `/tmp/hodge/results/reflection-coupled-residual/wave2-search.log.jsonl` (one JSON line per batch; remote lines merged by pull).
Library: `/tmp/hodge/results/reflection-coupled-residual/rcr2.py` (sector-aware E1 spaces, float cokernel projectors,
Hessian-relaxation tensors + LM search, msolve export, sector-restricted VarPro, exact confirmation).

## Method (new in wave 2)
* Klein group chi(e) = (Y2-parity, (Y0+Y1)-parity); chi(F) = (0,0), chi(A) = (0,1). Residue of A: r-odd, rho-odd.
* E2-mod-J obstruction ob(D) = [D div N − N(D)] ∈ R_{2d} is independent of the E1 solution N (Koszul and F-multiple
  corrections lie in J). For Jacobian-type D = Σ g_i F_i with N = A g this is the **Hessian relaxation**
  `ob(g) = −[A · gᵀ Hess(F) g]` — a pure quadratic map, cheap to evaluate (tensor K×K×dim R_{2d}).
  For the affine family D0 + Σ g_i F_i: ob(g) = [B(D0,N0)] + [D0(g(A)+A div g) − A g(D0) − N0ᵀHg] − [A gᵀHg].
  ob is invariant under the spurious subspace Z' = {g : D(g) ∈ F·S} (incl. syzygies); searches run on its complement.
* Exact decisions mod p (msolve) for small K; numerical Levenberg–Marquardt with thousands of restarts otherwise;
  every candidate below 1e-9 is rationally reconstructed, checked exactly for the relaxation, then solved exactly
  (sector-restricted compressed fixed-D solver) and verified.

## Results so far

### R1. E1-admissible spaces (exact mod p = 2^31−1), Klein sectors of D = (Y2-parity, Y0+Y1-parity)
| d | (0,0) | (0,1) | (1,0) | (1,1) | non-Jacobian classes |
|---|---|---|---|---|---|
| 4 | 0/48 | 0/38 | 0/22 | 0/18 | 0 (rigorous: no certificate of degree ≤ 4) |
| 5 | 3 = J | 2 = J | 1 = J | 0 | 0 |
| 6 | 14 = J (R_6 = 144) | 12 = J (R_6 = 126) | 6 = J | 4 = J | 0 |
(dim R_{10}: 460/426/318/302; dim R_{11}: 505/480/363/338; dim R_12^{(0,1)} = 492, consistent with the dossier.)
| 7 | 45 = J (R_7 = 213) | 40 = J (R_7 = 198) | 23 = J (R_7 = 135) | 18 = J (R_7 = 120) | 0 — **E1 forces D ∈ J_7 in every sector** |
| 8 | 113 = J(112) + 1 (K8) | (running) | (running) | (running) | 1 Y-even (K8) |
d = 9..12: remote jobs `search-rcr-e1lo` / `search-rcr-e1hi` (queued; remote at capacity).

### R2. Hessian relaxation at d = 5, 6 — EXACT (msolve, mod p, random affine chart), spurious directions killed
Family D = Σ g_i F_i (the whole E1-admissible space at these degrees). Quadrics = coordinates of [A gᵀHg] in R_{2d}^{(0,1)} (+ R_{2d}^{(1,0)} for mixed).
| d | sector | K (g-unknowns) | polys | result |
|---|---|---|---|---|
| 5 | 00 / 01 / 10 | 3 / 2 / 1 | 23 / 19 / 8 | empty / empty / empty |
| 6 | 00 / 01 / 10 / 11 | 14 / 12 / 6 / 4 | 296 / 294 / 234 / 66 | empty / empty / empty / empty |
| 6 | Y-even mixed (00⊕11) / Y-odd mixed (01⊕10) | 18 / 18 | 418 / 473 | empty / empty |
(sector 00 at d=6 first returned the Euler direction g = z/6, D = 6F — correctly flagged spurious by the exact pipeline; killed thereafter.)
d = 7: sector (1,1): K = 18, 426 polys, **empty** (0.2 s). Sector (0,0) did not finish in 100 s locally; remaining sectors queued (remote job `search-rcr-relaxlo`, 2 h timeouts).

### R2b. Whole degree-8 Y-even E1 family K8 ⊕ J_8^{(0,0)} (K = 113, spurious 10, search dim 103, 274 obstruction coordinates)
Generalised family relaxation (pairs (D_k, N_k), tensor Π B(D_k, N_l), validated against the closed form at d=6 to 15 digits).
40 LM restarts (local): **every restart converges to the same floor rel = 2.75e-3** (median 2.77e-3, none < 1e-4): the
obstruction map has no zero on this family numerically. Exact decision (msolve on the 274 quadrics in 103 unknowns)
queued in `search-rcr-relaxlo`.

### R3. K8 (unique non-Jacobian degree-8 Y-even E1 class, Klein-invariant): full fixed-D test — remote job `search-rcr-k8full`
E1 confirmed (rank 5052 = rank_aug). Sector-restricted compressed system 22016 × 4256 (of 7497 unrestricted): **inconsistent** (rank 3781, nullity 475, 36 s).
Unrestricted compressed system 40301 × 7497: **inconsistent** (rank 6703, nullity 794, 195 s). K8 alone is not a denominator; the affine family K8 + J_8 is attacked by the affine relaxation (queued).

**Consequence:** no certificate with a pure-Y-parity denominator of degree 5 or 6 exists (mod-p evidence; the relaxation is a necessary condition and E1 forces D ∈ J there). With the rigorous d ≤ 4 exclusion: deg D ≥ 7.
