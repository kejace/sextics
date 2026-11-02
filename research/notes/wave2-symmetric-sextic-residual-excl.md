# Wave 2 — `symmetric-sextic-residual` — lens EXCLUSION / RIGOR

STATUS: no witness (in progress; last rewrite 2026-11-01 19:40). Log: `/tmp/hodge/results/symmetric-sextic-residual/wave2-excl.log.jsonl`. Packet: `/tmp/hodge/results/symmetric-sextic-residual/wave2-excl/` (library `w2x_lib.py`, jobs `w2x_*_job.py` in the lane dir).

## Best rigorous bound so far

**Every legal certificate denominator (any Y-parity, J-elements included) has degree ≥ 7.** [P]
Proof structure: (i) E1 kernel `{D : A D ∈ J}/J_d = 0` for d ≤ 6 (modular rank bookkeeping, all ranks at their
complete-intersection values, see §2.1); so an E1-admissible D of degree ≤ 6 lies in J_5 or J_6 (J_d = 0 for d ≤ 4).
(ii) **Hessian relaxation** (new, §1): for D = w·∇F a certificate forces `A·Hess_F(w,w) ∈ J`. For d = 5 (w = c constant)
this is `A·(cᵀH_F c) ∈ J_10`; `[A]: R_4 → R_10` is injective (E1 kernel at d = 4 is zero) so `cᵀH_F c = 0` as a polynomial,
which would make `c` a singular point of X — contradiction (pen-and-paper; also confirmed by Singular: dim 0, vdim 4, both parities).
For d = 6 (w linear, 36 params, 1 radial direction `w ∝ z` giving D = 6F) the relaxation system (35 unknowns, 1605 quadrics,
mixed parity) has leading ideal containing `a_k²` for every variable (msolve, grevlex) ⇒ only the trivial solution ⇒ no legal D of degree 6.
**No `S3×S3×⟨τ⟩`-semi-invariant denominator of degree 7** [P]: all eight (character, parity) sector families at d = 7 are closed
(triv: dimK 8 = J-part, V={0}, vdim 11; sgnX even: ⟨D7⟩ ⊕ J (dimK 3), V={0}; sgnY odd: ⟨D7Y⟩ ⊕ J, V={0}; the other sectors
have dimK ≤ 1 and close too; sgn⊗sgn: no E1-admissible element at all). Folded-reducer ranks at degree 14 equal the
Molien–Koszul values (276, ring 46).
Running locally now (remote saturated): d = 8 sectors, d = 7 full pure parity (61 unknowns each), E1 kernels at p = 2147483629.

## 1. The Hessian relaxation (what is proved)

Let q = 1, F smooth (hence irreducible), D coprime to F with `A D ∈ J(F)`. Fix any polynomial vector N⁰ with
`Σ N⁰_i F_i = A D`. If (D, N, Q) is a certificate, then `A D − Σ N_i F_i = F H` (irreducibility + coprimality), hence
`N = N⁰ − (H/6) z + K_B` with `K_B` a Koszul syzygy of ∇F (F_i form a regular sequence), and the cleared identity becomes
`D div N⁰ − N⁰·∇D + (D div K_B − K_B·∇D) = F Q`. Every K_B term lies in J(F). Therefore

    R⁰(D) := D·div N⁰ − N⁰·∇D  ∈ J(F)_{2d}        (necessary; independent of the choice of N⁰).

On a linear family `D(a) = Σ a_k D_k`, `N⁰(a) = Σ a_k u_k` this is a system of quadrics in a; the solution set is a cone,
invariant under `D ↦ D + F g` (radial directions, illegal denominators), so we restrict to a complement of `F·S_{d−6}`
and ask whether the only solution is a = 0. For D = w·∇F ∈ J it reads `A·Hess_F(w,w) ∈ J`.

Rigor conventions (modular computation ⇒ characteristic-0 statement): a Q̄-certificate denominator D can be rescaled and
shifted by an F-multiple so that it is p-integral with reduction D̄ ∉ (F̄) (Krull intersection in the p-adic topology).
Because `rank_p J_{d+6} = rank_Q J_{d+6}` and `rank_p J_{2d} = rank_Q J_{2d}` (both checked against the CI Hilbert function
in every run), the cokernels are p-torsion-free, so `A D̄ ∈ J̄` and `R̄⁰(D̄) ∈ J̄`: D̄ is a nonzero, non-radial solution of the
mod-p relaxation. Hence "mod-p relaxation has only radial solutions" ⇒ "no certificate denominator of that degree over C".
The check that V(I) = {0} is done by msolve's leading ideal (`-g 1`): a pure power of every variable.

## 2. Re-verification of wave-1 / dossier claims (goal a)

### 2.1 Confirmed (independent method or second prime)
| claim | wave-1 method | wave-2 check | status |
|---|---|---|---|
| F = p6(X)+p6(Y)+(5/7)QX QY(QX+2QY) | exact | exact (`core.sub == {}`) | ✓ |
| smooth (F ≡ Fermat mod 5) | checks.json | exact | ✓ |
| A ∉ J(F) | rank 18/19 mod p | **exact over Q** (fmpq rank 18 vs 19) | ✓ |
| monomial automorphism group 36 × τ | monomial_aut2.py | 72 support-preserving perms, 36 exact | ✓ |
| residue character sgn⊗sgn | argument | re-derived | ✓ |
| isotypic table (mult.py, Fermat model) | equivariant Hilbert fn | **Molien–Koszul formula** Σ_g χ(g) det(1−t⁵g)/det(1−tg): R6=(17,1,5,26…), R12 (52, 10), R18 (17, 1) | ✓ |
| age screen (ages.py) | Reid–Tai | independent numpy eigen-decomposition: min age 1 on X for all 35 non-identity residue-invariant elements | ✓ |
| E1 kernel 0 for d ≤ 6, J-ranks = CI | p = 2^31−1 | p = 2147483629 job `excl-ssr-e1kernel` (queued) | pending |
| D7, D7Y inconsistent (engine, p1) | engine v2 | engine at p2 + twisted solver at p2 + **exact rational Farkas certificate** (job `farkas`) | pending |
| invariant kernel d=8,9 reps inconsistent | inv-v1 | twisted solver at p2 with recomputed kernels (`reverify`) | pending |
| served-span: every legal D of degree ≤ 6 serves only J | charsolve p1 | p2 rerun flagging F-multiples (`reverify`) | pending |

### 2.2 Findings against wave-1 outputs (not claims)
* `scan_served.out` line 78, `1EX(F)+1EY(F)`: this is `E_X(F)+E_Y(F) = 6F` — divisible by F, illegal. It "serves" all 17
  classes, and from there on every `cum_beyondJ=17 / A_in_span=True / "!!!!! A ENTERS THE SPAN"` flag (277 lines) is
  contamination; per-D results `beyondJ=0` for the legal D stand. The attack-1 note did not claim a span hit.
* `log.jsonl` rows ≥ 366 with `"consistent": true` (families sgn3…sgnsgn8, J6-radial-pencil, quot-curve6X): these come from
  `CharSystem(..., targets=True)` whose unknowns include the target itself (A' = −A always solves it); they are not certificates.
* The dossier's "[P] every certificate denominator not in J(F) has degree ≥ 7" was correct but silent about J_5, J_6; §0 closes that gap.

## 3. Experiments (wave 2)
| exp | family / D | d | prime | result | note |
|---|---|---|---|---|---|
| relaxation | all E1-admissible D (= J_5), per parity | 5 | 2^31−1 | V={0} (Singular dim 0) | + pen-and-paper proof |
| relaxation | all E1-admissible D (= J_6), mixed parity, 35 unknowns | 6 | 2^31−1 | V={0} (msolve: a_k² ∈ LT for all k) | rigorous |
| relaxation | all E1-admissible D (D7 ⊕ J_7 / D7Y ⊕ J_7), per parity, 61 unknowns | 7 | 2^31−1, 2147483629 | queued (`excl-ssr-relax1`) | |
| relaxation | mixed parity, 124 unknowns | 7 | 2^31−1 | queued (`excl-ssr-relax7mixed`) | may not finish |
| relaxation | sectors (triv, sgnX, sgnY, sgnsgn) × parity | 7 | 2^31−1 | **all V={0} / empty** (msolve + Singular agree) | rigorous: no semi-invariant D of degree 7 |
| relaxation | sectors × parity | 8, 9 | 2^31−1 (+p2 at d=8) | running locally (d=8) / remote queue (d=9) | |
| E1 kernels, all characters, isotypic decomposition | — | 1–10 | p2 (1–10), p1 (9,10) | queued | |
| screening (twisted E1+E2, random + basis points) | sectors | 8–12 | p1, p2 | to submit | |
| Farkas exact certificates | D7, D7Y, E_X(F), E_Y(F), Δ_XΔ_Y·…, D7·e1X, Δ_X·(e2e3−7e1³e2+4e1⁵) | 6–8 | exact | to submit | |

## 4. Remote jobs
Remote saturated by other agents (load ≈ 350–390 on 128 cores) at 19:30; submissions are queued with `--wait`.
