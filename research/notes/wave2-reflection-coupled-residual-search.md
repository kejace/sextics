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
d = 7..12: remote jobs `search-rcr-e1lo` / `search-rcr-e1hi` (pending capacity).

### R2. Hessian relaxation at d = 5, 6 — EXACT (msolve, mod p, random affine chart), spurious directions killed
Family D = Σ g_i F_i (the whole E1-admissible space at these degrees). Quadrics = coordinates of [A gᵀHg] in R_{2d}^{(0,1)} (+ R_{2d}^{(1,0)} for mixed).
| d | sector | K (g-unknowns) | polys | result |
|---|---|---|---|---|
| 5 | 00 / 01 / 10 | 3 / 2 / 1 | 23 / 19 / 8 | empty / empty / empty |
| 6 | 00 / 01 / 10 / 11 | 14 / 12 / 6 / 4 | 296 / 294 / 234 / 66 | empty / empty / empty / empty |
| 6 | Y-even mixed (00⊕11) / Y-odd mixed (01⊕10) | 18 / 18 | 418 / 473 | empty / empty |
(sector 00 at d=6 first returned the Euler direction g = z/6, D = 6F — correctly flagged spurious by the exact pipeline; killed thereafter.)
**Consequence:** no certificate with a pure-Y-parity denominator of degree 5 or 6 exists (mod-p evidence; the relaxation is a necessary condition and E1 forces D ∈ J there). With the rigorous d ≤ 4 exclusion: deg D ≥ 7.
