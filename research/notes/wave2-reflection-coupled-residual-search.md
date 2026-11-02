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
(filled in per batch)
