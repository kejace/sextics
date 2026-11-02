# Wave 2 — dense-rational — SEARCH lens

STATUS: NO WITNESS (2026-11-01 20:02 local). Searches RUNNING on the remote (hodge-sextics, load ~500/128 shared with ~50 agents): search-dense-rational-nvp7 (full-E2 varpro d=7, both parities; 4 restarts each so far, best rel 3.4e-5 / 3.9e-5, no zero), search-dense-rational-phase2 (Hessian relaxations d=7,8,9 real+complex, pure+mixed, then structured candidates; hr7e: 244/3000 restarts, best f 6.9e-4, no zero), search-dense-rational-prep2 (E1 spaces: numerical dims agree with the exact table through d=8 so far), search-dense-rational-nvp9 (submitted). Pull + summarize with: bash /tmp/hodge/results/dense-rational/w2s_pull.sh && python3 /tmp/hodge/results/dense-rational/w2s_report.py

Target: q=1, A = X0^3 Y0^3, F = Fermat + (1..17)/9973 on all 230 non-pure even sextic monomials. Only symmetry: tau (Y -> -Y);
the only "character sectors" are the tau-parities of D (Y-even / Y-odd / mixed).

Log: /tmp/hodge/results/dense-rational/wave2-search.log.jsonl (one JSON line per batch). Scripts: w2s_*.py in the lane dir.

## Plan
1. E1-admissible spaces per degree 4..12, both parities (numerical SVD + exact cross-checks); kernel structure.
2. Jacobian-type degrees 7, 8, 9: Hessian relaxation (thousands of restarts, real + complex, pure + mixed sectors), then
   full normalized-E2 variable projection over the gauge-fixed w-space (inner unknowns: Koszul corrections K, rows reduced mod F).
3. Degree 10 (first non-Jacobian degree): same E2 varpro over K_full(10) = ker[A] + J_10 (F-multiples removed) and the kernel-only part.
4. Structured candidates (Hessian covariants of the cubic cover, 3x3 Hessian minors, gradient-norm forms, cones) at d=10,12: E1 then E2.
5. Rational reconstruction + exact re-solve for every candidate with residual < 1e-10.

## Results
### Infrastructure (done)
* `w2s_common.py` / `w2s_nvp.py`: normalized-E2 variable projection.  For a linear family D(w) with E1 lifts N0(w), inner unknowns are the
  Koszul corrections K (skew, degree d-4) and the quotient C is eliminated by the normal form modulo F (division w.r.t. X0^6; kernel exactly
  F*S_{2d-6}^-).  Inner system sizes: d=7: 5160 x 420; d=8: 8661 x 954; d=9: 13717 x 1890; d=10: 20736 x 3465.  Gradient = Golub-Pereyra
  formula via adjoint polynomial multiplications; checked against finite differences (ratio 1.000000).  Validation: f = 5e-29 at the accepted
  degree-6 witness of `family-block-cube` (same target), f = 5e-4 after a 1e-3 perturbation of that D.
* `w2s_hessrelax.py`: Hessian relaxation w^T Hess(F) w in (J:A)_{2d-6} with gauge removed, real (L-BFGS on the sphere) and complex (LM) modes,
  pure and mixed tau-sectors; promotion of relaxation zeros to the full E2 functional.
* `w2s_e1spaces.py`: numerical SVD bookkeeping of J_k complements (k=10..18), K_full(d) bases (d=4..12), E1 lifts at d=10, kernel structure
  (intersections with coordinate subspaces / ideal powers), exact mod-p cross-checks.
* `w2s_structured.py`: 40 structured Y-even candidates at d=10 (gradient norms, second polars w^T H w for structured cubic fields) and
  d=12 (Hessian covariants of the cubic cover a+bt+ct^2+dt^3, the 20 principal 3x3 Hessian minors, det Hess of the pure blocks, Euler cross sums).

### 1. E1-admissible denominator spaces per degree (both tau-parities) — exact
K_full(d) = {D in S_d^par : A D in J(F)} = J_d + (kernel of [A] on the Jacobian ring).  The kernel dimensions are exact: zero for d <= 9
(fqlx for d <= 8; dossier d = 9 mod 101 and 10007; wave-1 direct mod 2^31-1 for the even sector), and for d >= 10 equal to the minimal
value dim R_d - rank[A] with rank[A] = dim R_{18-d}^{opposite parity} by Macaulay duality (the dual map is injective because 18-d <= 8);
d = 10 also directly (dossier).  "search dim" = dim K_full minus the F-multiples F*S_{d-6} (spurious denominators).

| d | par | dim S_d | dim J_d | dim R_d | rank[A] (target dim) | ker mod J | dim K_full | F-mult | search dim |
|---|---|---|---|---|---|---|---|---|---|
| 4 | even/odd | 66/60 | 0 | 66/60 | 66/60 (744/762) | 0 | 0 | 0 | 0 |
| 5 | even/odd | 126/126 | 3 | 123 | 123 (843) | 0 | 3 | 0 | 3 |
| 6 | even/odd | 236/226 | 18 | 218/208 | 218/208 (866/885) | 0 | 18 | 1/0 | 17/18 |
| 7 | even/odd | 396 | 63 | 333 | 333 (843) | 0 | 63 | 3 | 60 |
| 8 | even/odd | 651/636 | 168 | 483/468 | 483/468 (744/762) | 0 | 168 | 12/9 | 156/159 |
| 9 | even/odd | 1001 | 378 | 623 | 623 (623) | 0 | 378 | 28 | 350 |
| 10 | even/odd | 1512/1491 | 750/747 | 762/744 | 468/483 (468/483) | **294/261** | 1044/1008 | 66/60 | 978/948 |
| 11 | even/odd | 2184 | 1341 | 843 | 333 (333) | **510/510** | 1851 | 126 | 1725 |
| 12 | even/odd | 3108/3080 | 2223/2214 | 885/866 | 208/218 (208/218) | **677/648** | 2900/2862 | 236/226 | 2664/2636 |

Consequences for the search: through d = 9 every admissible denominator is Jacobian-type D = w(F) (w of degree d-5, modulo the Euler gauge
w -> w + l*Euler which produces the F-multiples); the outer search spaces are 60 / 156-159 / 350 dimensional at d = 7 / 8 / 9, with inner
Koszul-correction spaces of 420 / 954-936 / 1890 unknowns and 5160 / 8661 / 13717 equations after reduction modulo F.  The Hessian
relaxation w^T Hess(F) w in (J:A)_{2d-6} has 483 / 468 / 208 equations: overdetermined at d = 7, 8 and underdetermined at d = 9 (a
>= 142-dimensional relaxation variety), vacuous at d >= 10 (dim R_20^- is tiny) — so at d = 9 and 10 only the full E2 functional can decide.
(Numerical SVD cross-checks of this table by the prep job are recorded in w2s_e1spaces.json / log events `e1-space`.)

### Batches (live at 20:02; final numbers land in wave2-search.log.jsonl as `nvp-batch` / `hessrelax-batch` events and in *.restarts.jsonl)
| search | degree | sector | outer dim | inner / equations | restarts so far | best residual | median | zeros (<1e-10) |
|---|---|---|---|---|---|---|---|---|
| full E2 varpro (jac family) | 7 | even | 60 | K 420 / 5160 rows mod F | 4 of 4000 (~500 s each under load) | rel 3.45e-5 | 3.5e-5 (all four in one basin) | 0 |
| full E2 varpro (jac family) | 7 | odd | 60 | K 420 / 5160 | 4 of 4000 | rel 3.92e-5 | 1.1e-2 | 0 |
| Hessian relaxation, real L-BFGS | 7 | even | 60 | 483 quadrics | 244 of 3000 | f 6.87e-4 | 7.8e-4 (q10 7.0e-4, q90 8.3e-4; all hit maxiter 500) | 0 |
| Hessian relaxation (odd, complex, mixed 120 unknowns / 951 eqs), d=8 (156/159 unknowns, 468 eqs), d=9 (350 unknowns, 208 eqs) | 7-9 | all | | | queued in phase2 | | | |
| full E2 varpro (jac family) | 9 | even/odd | 350 | K 1890 / 13717 | job nvp9 submitted (2 procs x 2 threads per sector) | | | |
| full E2 varpro over K_full(10) (kernel-only 294 and kernel+J 978) | 10 | even | 294 / 978 | K 3465 / 20736 | w2s_nvp10.sh ready; needs the prep2 cache (families at d=10 with E1 lifts) and a free job slot | | | |
| structured candidates (40 at d=10,12) | 10,12 | even | | numerical E1 via J_{d+6}^- complement, then exact E1, then E2 | end of phase2 | | | |

Reading of the landscape so far: at d = 7 every full-E2 restart converges to residual floors ~3.5e-5 (even) / 3.9e-5 (odd) and the
relaxation to ~7e-4 with a very narrow spread (q10..q90 = 7.0e-4..8.3e-4), i.e. a single least-squares compromise basin and no sign of a zero,
consistent with the 8:1 overdetermination (483 quadrics / 60 unknowns) and with the rigorous degree-7 msolve certificates being run by the
exclusion lens (jobs excl-dense-rational-d7even/odd).  No candidate reached the reconstruction threshold (rel < 1e-10); no rational reconstruction
was attempted; no witness.

### Numerical E1 bookkeeping (prep2, so far d = 4..8 both parities)
dim J_d, dim K_full(d), rank[A] all agree with the exact table above; SVD gaps: last nonzero singular value / s0 between 2e-7 and 0.7, first
zero ~1e-15 (clean).  Events `e1-space`, `J-complement`, `K_full` in the log.

### Rigorous side (not this lens)
The exclusion lens re-proved deg D >= 7 in all parity sectors (msolve leading ideals at d = 6) and is certifying d = 7 pure sectors; my msolve
job (w2s_ms7.sh) was therefore not submitted.

## Continuation (what the remote jobs will leave behind, and what to run next)
1. `bash w2s_pull.sh` mirrors /tmp/hodge/results/dense-rational from the remote into /tmp/hodge/remote/dense-rational and merges the remote-only
   log (wave2-search.remote.log.jsonl) into wave2-search.log.jsonl; `python3 w2s_report.py` prints the batch tables.
2. When prep2 has written w2s_cache/family_ker_d10_par0.npz and family_kfull_d10_par0.npz: submit
   `python3 -m hodge.rjob submit --lane dense-rational --name search-dense-rational-nvp10 --cpus 8 -- "bash w2s_nvp10.sh"`.
3. When phase2 has produced hr9e/hr9o promoted points (hr9e.seed*.w.npy, relaxation zeros on the >= 142-dim variety): seed the d=9 full E2 with
   `W2S_REMOTE=1 W2S_THREADS=2 python3 w2s_nvp.py --d 9 --pD 0 --kind jac --init-glob "hr9e.seed*.w.npy" --restarts 0 --procs 4 --name nvp9e-seeded`.
4. Any restart with rel < 1e-10 is reconstructed automatically (denominators 1..1e6, exact fixed-D re-solve via hodge.fixed_d, F-divisibility check)
   and written to witness-candidate.json with a WITNESS log event; then run `python3 -m hodge.leangen` + `hodge.leancheck` locally.

