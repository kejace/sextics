# Wave-2 joint search — lane `reflection-coupled-residual` (q=1, A = X0^3 Y0^3)

STATUS: NO WITNESS FOUND (search running; report rewritten after every batch). Started 2026-11-01 19:05 local.
Remote box saturated since ~19:25 (load 240–500 on 128 cores, 22–97 jobs); submissions `search-rcr-e1lo`, `search-rcr-relaxlo`,
`search-rcr-e1hi` wait for capacity (runner gate); low-degree work done locally within the light-compute budget.
Queued next (scripts ready in the lane dir): `job_kvp_lo.sh` (full-E2 Koszul varpro d=7,8, 8 seeds), `job_kvp_d9.sh`, `job_kvp_d10.sh`,
`job_screen.sh` (structured candidates, E1 screen + exact solve), `job_served.py` (served-space profiling of 121 denominators).

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

* **Koszul-reduced full-E2 variable projection** (`kvp.py`): after solving E1 exactly, every E1 solution is N0 + K_a with
  K_a,i = Σ_j a_ij F_j (a antisymmetric, a_ij ∈ S_{d−4}) and B_D(K_a) = Σ_{i<j}[D(F_j∂_i a_ij − F_i∂_j a_ij) − a_ij(F_jD_i − F_iD_j)];
  the certificate exists iff B_D(N0) ∈ (F) + span(B_D(K_a)). Over a family D(c) = Σ c_k D_k this is
  f(c) = min_a ‖Π_F(L(c)a − quad(c))‖², with ~6× fewer unknowns than the raw numerators and the (F)-part projected out
  once; sparse index arrays make each evaluation + analytic gradient cheap (0.05–0.3 s at d = 8–9 in one sector).
  Validated: assembly vs brute force (1e-13); finite-difference gradient (6e-10); **positive controls** — Fermat sextic with
  A = X0³Y0³, d=6 J-family: residual 4e-22 → reconstruction → exact witness verified (D = X0⁶+X1⁶+X2⁶-type);
  `fermat-pair-rational-d4-strict` (d=4 kernel family, 13 classes): zeros found, exact witness recovered and verified
  (also with the sparse NF-mod-F + random-compression metric used for the large degrees: rel 1.7e-16 → verified witness).
  Two optimizers: L-BFGS variable projection (min over a inside) and a joint Levenberg–Marquardt on (c, a) with the exact
  Jacobian (controls: d=4 rel 2.6e-31 in 13 evaluations, Fermat d=6 rel 2.5e-26 in 58); at d=8 (103+538 unknowns) one
  LM restart ≈ 240 s, one L-BFGS restart ≈ 150 s (400 iterations, not converged); both end at rel 0.1–1 on our family.

* Scope: denominators of pure Y-parity (τ-eigenvectors), the engine's and the dossier's convention. A mixed-parity D of
  degree d gives the Y-even denominator D·τ(D) of degree 2d (X∖V(Dτ(D)) ⊂ X∖V(D)), so pure-parity searches up to degree
  2d cover mixed-parity denominators up to degree d. Within pure parity, both Klein-eigen (4 sectors) and Klein-mixed
  ('e' = 00⊕11, 'o' = 01⊕10) families are searched.

* **Served-space profiling** (`job_served.py`): for a denominator D, the whole space of Y-odd sextic targets it serves
  (compressed fixed-D system with target columns, nullspace mod p), reduced mod J_6 and split by Klein sector. A D whose
  served space meets the target's sector (0,1) is a building block even if it does not serve A (served(D1D2) =
  served(D1)+served(D2)). Validated against wave 1 on degree-4 cones (served = J only, 2 s each vs 20–80 s).

## Results so far

### R1. E1-admissible spaces (exact mod p = 2^31−1), Klein sectors of D = (Y2-parity, Y0+Y1-parity)
| d | (0,0) | (0,1) | (1,0) | (1,1) | non-Jacobian classes |
|---|---|---|---|---|---|
| 4 | 0/48 | 0/38 | 0/22 | 0/18 | 0 (rigorous: no certificate of degree ≤ 4) |
| 5 | 3 = J | 2 = J | 1 = J | 0 | 0 |
| 6 | 14 = J (R_6 = 144) | 12 = J (R_6 = 126) | 6 = J (R_6 = 82) | 4 = J (R_6 = 74) | 0 (recomputed, `e1_6_*.json`) |
(dim R_{10}: 460/426/318/302; dim R_{11}: 505/480/363/338; dim R_12^{(0,1)} = 492, consistent with the dossier.)
| 7 | 45 = J (R_7 = 213) | 40 = J (R_7 = 198) | 23 = J (R_7 = 135) | 18 = J (R_7 = 120) | 0 — **E1 forces D ∈ J_7 in every sector** |
| 8 | 113 = J(112) + 1 (K8), R_8 = 301 | 104 = J (R_8 = 274) | 64 = J (R_8 = 194) | 56 = J (R_8 = 182) | **1 in all sectors: K8** (sector 00) |
| 9 | 263 = J(242) + **21** (R_9 = 381 → R_15 = 360: exactly the forced minimum) | 228 = J (R_9 = 360) | 171 = J(150) + **21** (R_9 = 263) | 136 = J (R_9 = 242) | **42**: 21 in sector 00 (Y-even) + 21 in sector 10 (Y-odd); both exactly the forced minimum dim R_9 − dim R_15; dense representatives (361 terms) |
d = 10..12: remote job `search-rcr-e1hi` (queued; remote at capacity).

### R2. Hessian relaxation at d = 5, 6 — EXACT (msolve, mod p, random affine chart), spurious directions killed
Family D = Σ g_i F_i (the whole E1-admissible space at these degrees). Quadrics = coordinates of [A gᵀHg] in R_{2d}^{(0,1)} (+ R_{2d}^{(1,0)} for mixed).
| d | sector | K (g-unknowns) | polys | result |
|---|---|---|---|---|
| 5 | 00 / 01 / 10 | 3 / 2 / 1 | 23 / 19 / 8 | empty / empty / empty |
| 6 | 00 / 01 / 10 / 11 | 14 / 12 / 6 / 4 | 296 / 294 / 234 / 66 | empty / empty / empty / empty |
| 6 | Y-even mixed (00⊕11) / Y-odd mixed (01⊕10) | 18 / 18 | 418 / 473 | empty / empty |
(sector 00 at d=6 first returned the Euler direction g = z/6, D = 6F — correctly flagged spurious by the exact pipeline; killed thereafter.)
d = 7: sector (1,1): K = 18, 426 polys, **empty** (0.2 s; also empty for the generalised family relaxation, 426 polys); sector (1,0): K = 23, 427 polys, **empty**. Sectors (0,0) (45 unknowns, 429 polys) and (0,1) (40, 428) and the mixed families e/o (63, 744) did not finish in the local budget (exports `relax_d7_*.ms`; remote job `search-rcr-relaxlo`, 2 h timeouts).
Cross-check (aggregator status note v8.8, `wave2-status.md` §2.14): the excl lens closed the characters (1,1), (1,0) **and (0,1)** at p = 2147483587 (colengths 85 / 239 / 1109) and the mixed pair (1,0)+(1,1) (mult 4399), re-derived by the aggregator at 998244353 with identical quadric counts and colengths; the (0,0) block (42 unknowns after killing the spurious directions, 274 quadrics) was still running there. Our (1,1)/(1,0) closures are a third independent derivation (random affine chart, 2^31−1).

### R2a. Full E2 (Koszul varpro, sparse metric) over the complete degree-7 families — local batch B
E1 forces D ∈ J_7 in every sector (R1), so these families are everything at degree 7.
| sector | K | search dim | n_a | restarts | best rel | median rel |
|---|---|---|---|---|---|---|
| 00 | 45 | 42 | 242 | 12 | 0.495 | 0.699 |
| 01 | 40 | 38 | 238 | 12 | 0.699 | 0.943 |
| 10 | 23 | 22 | 182 | 12 | 1.048 | 1.717 |
| 11 | 18 | 18 | 178 | 12 | 1.564 | 1.81 |
| e (00⊕11) | 63 | 60 | 420 | 4 | 0.970 | 1.19 |
| o (01⊕10) | 63 | 60 | 420 | 4 | 0.633 | 0.95 | No descent below 0.1: the full-E2 landscape at degree 7 has no zero in reach of 12 restarts per
sector, consistent with the exact relaxation emptiness in sectors 10/11 (R2).

### R2b. Whole degree-8 Y-even E1 family K8 ⊕ J_8^{(0,0)} (K = 113, spurious 10, search dim 103, 274 obstruction coordinates)
Generalised family relaxation (pairs (D_k, N_k), tensor Π B(D_k, N_l), validated against the closed form at d=6 to 15 digits).
40 LM restarts (local): **every restart converges to the same floor rel = 2.75e-3** (median 2.77e-3, none < 1e-4): the
obstruction map has no zero on this family numerically. Exact decision (msolve on the 274 quadrics in 103 unknowns)
queued in `search-rcr-relaxlo`.

### R2c. Full E2 (Koszul varpro) on the degree-8 Y-even family K8 ⊕ J_8^{(0,0)} — local, preliminary
K = 113, search dim 103, Koszul unknowns n_a = 538, quotient rows 4808 (+ compression to 598): first restarts give
rel ≈ 0.08–0.12 (f ≈ 0.01) after 400 L-BFGS iterations, no descent towards zero; one evaluation ≈ 0.3 s on the loaded
local box. Hundreds of restarts queued remotely (`job_kvp_lo.sh`, 8 seeds × 80 restarts).

### R2d. Full E2 (Koszul varpro, L-BFGS, sparse metric) over the degree-8 Jacobian families in the other sectors — local batch C
| sector | family | K | search dim | n_a | restarts | best rel | median rel |
|---|---|---|---|---|---|---|---|
| 01 | J | 104 | 98 | 522 | 3 | 0.158 | 0.189 |
| 10 | J | 64 | 61 | 414 | 3 | 0.228 | 0.250 |
| 11 | J | 56 | 54 | 416 | 3 | 0.338 | 0.439 |
| 00 | K8 ⊕ J | 113 | 103 | 538 | 2 (L-BFGS, 500 it.) | 0.104 | 0.139 |
No degree-8 sector shows any descent below rel 0.1 (earlier dense-metric and LM runs on 00 JK: 0.08–1.07).

### R2e. Full E2 (Koszul varpro) over the degree-9 E1 families (21 non-Jacobian classes ⊕ J_9 per sector) — local batch D
| sector | family | K | search dim | rows | n_a | nQ | restarts | best rel | median rel | note |
|---|---|---|---|---|---|---|---|---|---|---|
| 00 | ker(21) ⊕ J_9(242) | 263 | 241 | 9306 | 1050 | 1764 | 1 | 0.0180 | 0.0180 | L-BFGS 300 it., ~8 min per restart locally |
| 00 | ker(21) ⊕ J_9(242) | 263 | 241 | 9306 | 1050 | 1764 | 1 | 0.0208 | 0.0208 | L-BFGS 300 it., seed 42 |
| 00 | ker(21) ⊕ J_9(242) | 263 | 241 | 9306 | 1050 | 1764 | 1 | 0.7676 | 0.7676 | joint LM (≤60 evaluations), seed 43: poor local minimum |
| 10 | ker(21) ⊕ J_9(150) | 171 | 161 | 9306 | 850 | 1764 | 1 | 0.0214 | 0.0214 | L-BFGS 300 it., seed 44 |
| 10 | ker(21) ⊕ J_9(150) | 171 | 161 | 9306 | 850 | 1764 | 1 | 0.0256 | 0.0256 | L-BFGS 300 it., seed 45 |
| 00 | ker(21) ⊕ J_9(242) | 263 | 241 | 9306 | 1050 | 1764 | 1 | **0.0075** | 0.0075 | L-BFGS **900 it.**, seed 41 (1288.4 s) |
| 10 | ker(21) ⊕ J_9(150) | 171 | 161 | 9306 | 850 | 1764 | 1 | **0.0131** | 0.0131 | L-BFGS **900 it.**, seed 51 (721.2 s, tag localF) |
Floors ≈ 0.02 in both kernel-bearing sectors after 300 iterations (vs 0.1–0.2 at degree 8), and both landscapes keep
descending slowly with more iterations (sector 00: 0.018 → 0.0075; sector 10: 0.021 → 0.013 at 900 iterations; not converged). Whether this approaches zero or a
positive floor is the key open numerical question; a 1500-iteration run is queued (the logged seed-41 entry with rel 4.36 is a
diagnostic 3-iteration run, not a search result).

### R3. K8 (unique non-Jacobian degree-8 Y-even E1 class, Klein-invariant): full fixed-D test — remote job `search-rcr-k8full`
E1 confirmed (rank 5052 = rank_aug). Sector-restricted compressed system 22016 × 4256 (of 7497 unrestricted): **inconsistent** (rank 3781, nullity 475, 36 s).
Unrestricted compressed system 40301 × 7497: **inconsistent** (rank 6703, nullity 794, 195 s). K8 alone is not a denominator; the affine family K8 + J_8 is attacked by the affine relaxation (queued). The excl lens made the K8 exclusion rigorous (`A·K8 ∈ J` exactly over ℚ by a CRT-lifted E1 lift; rank_ℚ = rank_p certified by 482 a-priori kernel vectors; aggregator agrees at 10⁹+7) and found that no first-order correction inside the E1 family K8 + Σ g_iF_i removes K8's obstruction (linearised obstruction rank 156, augmented 157) — consistent with our 2.75e-3 relaxation floor on that family.

**Consequence:** no certificate with a pure-Y-parity denominator of degree 5 or 6 exists (mod-p evidence; the relaxation is a necessary condition and E1 forces D ∈ J there). With the rigorous d ≤ 4 exclusion: deg D ≥ 7.

## Running summary (kept current; the final 10-line summary is derived from it)
1. NO WITNESS FOUND. The exact and numerical searches below leave every tested family without a certificate denominator.
2. deg D ≤ 4: rigorous (E1 space zero). d = 5, 6: every E1-admissible D lies in J(F) in all four Klein sectors and the Hessian
   relaxation A·gᵀHess(F)g ∈ J is EMPTY in every sector and both mixed parities (msolve mod 2^31−1, spurious directions killed)
   → no pure-parity certificate of degree ≤ 6 (mod-p evidence; extends the rigorous ≤ 4 bound).
3. d = 7: E1 forces D ∈ J_7 in all sectors (45/40/23/18); relaxation empty in sectors 10, 11 (exact, ours); sectors 00/01 and the
   mixed families exceed the local msolve budget (queued remotely) — the excl lens and the aggregator independently closed (1,1), (1,0),
   (0,1) and the mixed pair (1,0)+(1,1), leaving only the (0,0) block (and the other mixed pairs) open at degree 7; full-E2 Koszul
   varpro over all six degree-7 families (12 restarts each) never descends below rel 0.5, consistent with those closures.
4. d = 8: the only non-Jacobian E1 class in any sector is K8 (sector 00); K8 alone is inconsistent (restricted and full
   systems; made rigorous over ℚ by the excl lens); the family K8 ⊕ J_8 has a relaxation floor 2.75e-3 (40 LM restarts; the excl
   lens's linearisation shows no first-order correction exists) and full-E2 floors 0.08–0.16 (L-BFGS and joint LM); the J-only
   families in sectors 01/10/11 stay at 0.16–0.34.
5. d = 9 (new): E1 spaces computed in all sectors — 21 non-Jacobian classes in sector 00 and 21 in sector 10 (exactly the
   forced minimum dim R_9 − dim R_15), none in 01/11. Full-E2 Koszul varpro over the kernel ⊕ J_9 families: floors
   ≈ 0.02 after 300 L-BFGS iterations, descending slowly with longer runs (00: 0.018 → 0.0075 at 900 it.; 10: 0.021 → 0.013);
   a 1500-iteration run is in progress; no candidate below 1e-9, nothing to reconstruct.
6. Tools (all validated on Fermat controls that recover exact witnesses): sector-aware E1 spaces with float bases, Hessian
   relaxation (exact msolve export + LM), generalised family relaxation, Koszul-reduced full-E2 variable projection with sparse
   quotient-mod-F and random compression (L-BFGS and joint LM), served-space profiler, structured-candidate screen.
7. Remote box saturated all evening (load 240–1900, 150–230 jobs, free memory down to 0 GB): the 8-seed searches, the exact
   msolve of the degree-7/8 families, d = 10–12 kernels and the structured/served sweeps are queued (scripts ready in the lane
   dir: job_relax_lo.sh, job_kvp_lo.sh, job_kvp_d9.sh, job_kvp_d10.sh, job_screen.sh, job_served.py) behind the capacity gate.

## Auto-generated tables
Regenerated from the log by `python3 summarize_log.py` into `wave2-reflection-coupled-residual-tables.md` (same directory).
