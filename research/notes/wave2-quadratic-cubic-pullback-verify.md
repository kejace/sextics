# Wave 2 / VERIFY lens: `quadratic-cubic-pullback` (q = 1, A = X0^3 Y0^3, F = sum Q_i^3)

STATUS: IN PROGRESS (started 2026-11-02 00:25 UTC). confirmed 12 / refuted 0 / weakened 0 / unverified 0 (table below is filled as checks finish)

Role: adversarial referee of every claim labelled rigorous/proven/exact in
`dossier-quadratic-cubic-pullback.md`, `attack1-quadratic-cubic-pullback.md`, `wave2-quadratic-cubic-pullback-*.md`,
the lane logs (`log.jsonl`, `kernel_sweep.jsonl`, `w2_kernels.jsonl`) and (when it exists) the draft post.
Independent code lives in `/tmp/hodge/results/quadratic-cubic-pullback/verify/` (own polynomial arithmetic, NOT the
engine's `hodge.core`/`hodge.fixed_d`; different primes; Singular where feasible). Per-check JSON lines are appended
to `/tmp/hodge/results/quadratic-cubic-pullback/wave2-verify.log.jsonl`. Job logs: `/tmp/hodge/jobs/verify-quadratic-cubic-pullback-*.log`.

Verdict legend: CONFIRMED = logic checked and recomputed independently; WEAKENED = true in a narrower form than stated
(the row says exactly what survives); REFUTED = false as stated; UNVERIFIED = not checked / not checkable here.

## Claim table

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
| C1 | `F = sum_i Q_i^3` exactly with the six listed quadrics; F Y-even, 50 terms, six pure sixth powers | dossier 1.1, attack 1.1, benediamond | CONFIRMED | own sympy exact: `F - sum Q^3 = 0` | — |
| C3 | good reduction mod 5 (`F ≡ Fermat`, partials `≡ Z_i^5`) hence X smooth | dossier 1.4, checks.json | CONFIRMED | all 50 denominators coprime to 5; `F mod 5 = sum Z_i^6` (own code) | openness-of-smoothness step is standard |
| C4 | projective diagonal symmetry group exactly `Z/2 = {1,tau}`; SNF `(1,1,1,1,2)` | dossier 1.3, struct1.py | CONFIRMED | own sympy SNF of exponent-difference lattice: diag `[1,1,1,1,2,0]`, torsion `Z/2` | so no finer character than Y-parity: correct |
| C5 | monomial symmetries (perm x sign, 46080 maps): only `1, tau` | dossier 1.3 | CONFIRMED | own enumeration: 4 affine maps `{±1, ±tau}` = projectively `{1,tau}` | no reflection, no paired signs, no X<->Y exchange, no cyclic shift |
| C6 | no continuous symmetry (`hessian_symmetrizer_dimension = 1`) | checks.json, dossier 1.3 | CONFIRMED | Lie algebra `{M in gl_6 : (MZ.d)F = lam F}` has dim 1 (Euler) mod 1e9+7 (own code) | — |
| C7 | `h^{3,1}(X)^+ = 226-18 = 208`, `h^{3,1}(X)^- = 236-18 = 218`, `h^{3,1} = 426`; bigraded Hilbert function of R | dossier 1.5/1.6, benediamond, attack 1.5 table | CONFIRMED | exact rank over Q of the 18+18 degree-6 relations; own bigraded Hilbert series: `R_6=(218,208)`, `R_10=(762,744)`, `R_16=(483,468)`, total `1,6,21,...,1751,...,1` | matches every number quoted in both notes |
| C8 | `delta = det(dQ/dZ)`: 64 terms, Y-odd, `≡ -X0X1X2Y0Y1Y2 (mod 5)`; `A ∉ J_6`, `delta ∉ J_6`, `A, delta` independent mod `J_6` | dossier 1.5/1.6(4), why.md, benediamond item 2 | CONFIRMED | own sympy determinant; exact ranks over Q: `J_6^odd+A = 19`, `+delta = 19`, `+A+delta = 20` | exact over Q (stronger than the mod-p statements in the notes) |
| C9 | Fermat-limit strict witness `v = A E_X/(6(X0^6+X1^6+X2^6))` for `(Fermat,1,X0^3Y0^3)` | dossier 1.8, attack §3, witness-candidate-fermat-92047.json | CONFIRMED | own exact expansion of the cleared identity: `E = 0` | control only; says nothing about the actual F |
| C10 | F is not a join `G(l)+H(l')` in any coordinates: 21 second partials independent | attack 1.1/1.4, jointest.py | CONFIRMED | exact rank over Q = 21 | argument also excludes every split `k+(6-k)` (a join has rank `<= 21-k(6-k)`) |
| C11 | `Q_0..Q_5` regular sequence, `B=S/(Q)` Hilbert series `(1+t)^6`, tau-trace `(1-t^2)^3`, `(phi_*O)^+ = O+6O(-1)+9O(-2)`; phi finite | benediamond (1),(4); dossier 1.5 | CONFIRMED | own bigraded Hilbert function of `S/(Q)` mod 1e9+7: `1 | 3+3 | 6+9 | 10+10 | 9+6 | 3+3 | 0+1 | 0` (even+odd) | mod-p vanishing of `(S/(Q))_7` implies vanishing over Q, so no base point over C |
| C12 | the 56-dim space of cubics in `Q` meets `J_6^even` exactly in `C·F` | dossier 1.6(3), qcp_polar.py | CONFIRMED | own ranks mod 1e9+7: 56, 18, intersection 1 | — |
| C26 | book §3.4 (weighted Fano / Piola) serves only `A = delta·b(Phi)` with `wt(b) = 12 - M = 0` at q=1, i.e. the line `C·delta` | dossier 0, 2; why.md | CONFIRMED | read `docs/03-geometric-support.md` §3.4 theorem + the q=1 remark (`wt(b)=12-M`) | the dossier's reading of the theorem is correct |
