# Referee README — wave-2 exclusion packet, `symmetric-sextic-residual` (q = 1, A = X0X1X2Y0Y1Y2)

Referee: a Claude (Anthropic) session acting for `kejace`, 2026-11-01 23:00–23:59 local (brunos). No human review yet;
nothing here is Lean-verified; no certificate is claimed; no claim of unrestricted nonexistence is made.
Referee scripts and raw outputs: `referee/` (this directory). Scratch copy: `/tmp/hodge/referee-ssr/`.

Verdict legend. **VERIFIED** = replayed by the referee with independent code or by re-running the stored input, or
confirmed by two independent implementations whose receipts the referee inspected (marked "not replayed" where so).
**PLAUSIBLE** = consistent receipts from one implementation / one prime, or a numerical or modular screen; not a theorem.
**UNSUPPORTED** = stated or drafted but not established by the packet. **REFUTED** = wrong as stated.

## 0. Counts

VERIFIED 14 · PLAUSIBLE 6 · UNSUPPORTED 2 · REFUTED 1 (claim 5, corrected in place). Strongest verified statement:
**every legal certificate denominator has degree ≥ 7 (any Y-parity, Jacobian-ideal elements included), and no
S3×S3×⟨τ⟩-semi-invariant denominator of degree 7 or 8 exists.**

## 1. Claim table

| # | claim | source | verdict | referee evidence |
|---|---|---|---|---|
| 1 | `F = p6(X)+p6(Y)+(5/7)e2(X)e2(Y)(e2(X)+2e2(Y))`, `A = e3(X)e3(Y)`, q = 1 | dossier §1.1, excl log `reverify-structure` | VERIFIED (not replayed) | two exact comparisons (dossier, excl lens); `instance.rec` read |
| 2 | X smooth (F ≡ Fermat mod 5) | checks.json | VERIFIED | 5/7 ≡ 0 mod 5 ⇒ F̄ = Σ z_i⁶, partials z_i⁵: m-primary |
| 3 | `A ∉ J(F)` | dossier, excl (exact over ℚ), aggregator (10⁹+7) | VERIFIED | referee: rank J_6^odd = 18, rank [J_6^odd | A] = 19 at p = 10⁹+9 (`referee/e1_replay.out`) |
| 4 | residue character of [A] is sgn⊗sgn, τ-character +1 | dossier §1.3 | VERIFIED | by hand: A is S3×S3-invariant and Y-odd; det(π_X,π_Y) = sgn·sgn, det τ = −1 |
| 5 | Hodge numbers of the target's sector | dossier §1.4 **(1,10,1)**; verify lens + aggregator **(17,52,17)** | **REFUTED → corrected** | referee Burnside recount in the Fermat model (`referee/burnside_sectors.out`): polynomial (triv,triv) Y-odd = residue (sgn,sgn,+) = **(17, 52, 17)**; the (1,10,1) row is the residue-(triv,triv) sector spanned at (3,1) by Δ_XΔ_Y |
| 6 | E1 kernel `{D : AD ∈ J}/J_d = 0` for d ≤ 6, both parities, with `rank_p J_{d+6} = CI` | excl (2³¹−1, 2147483629), aggregator (10⁹+7) | VERIFIED | referee at p = 10⁹+9, d = 4..7: all ranks = CI, dimK = dim J_d (`referee/e1_replay.out`) |
| 7 | d = 7 kernel beyond J: 1 per parity (`D7 = Δ_X e2(X0²,X1²,X2²)`, `D7Y`) | excl, aggregator | VERIFIED | referee: ker_R = 1 per parity at 10⁹+9; the pickled `D7` equals the referee's own `Δ_X·e2(X²)` (`referee/farkas_D_check.txt`) |
| 8 | d = 8 kernel beyond J: 6 per parity | excl (2 primes), aggregator | VERIFIED (not replayed) | two implementations, three primes, identical dims and isotypic split |
| 9 | **deg D ≥ 7 for every legal denominator** (d = 5: Hessian argument; d = 6: relaxation on J_6, 35 unknowns) | excl §1 | **VERIFIED** | §2 below: derivation re-done by hand; d = 5 argument checked by hand and by Singular on the stored files (dim 0, vdim 4, both parities); stored d = 6 msolve replayed single-threaded (identical leading ideal, 1307 tokens, pure power of every one of the 35 variables); **own-code reconstruction** of the d = 6 relaxation at p = 10⁹+9 (`referee/relax6_own.py`, 1684 quadrics / 35 unknowns) is m-primary (pure powers, exponents 2–4) |
| 10 | no `S3×S3×⟨τ⟩`-semi-invariant denominator of degree 7 (8 sectors) | excl | VERIFIED | all six nonempty sector systems replayed with msolve (identical leading ideals, pure powers for every variable); Singular on triv/even: dim 0, vdim 11; sgn⊗sgn sectors have no E1-admissible element (dimK 0, both parities; search lens's kernel table agrees); folded ring dim 46 = referee Burnside count |
| 11 | no `S3×S3×⟨τ⟩`-semi-invariant denominator of degree 8 (8 sectors) | excl | VERIFIED | same replay (triv/even 18 unknowns, 91 s single-threaded); Singular vdim 112 (triv/odd), 7 (sgnX/odd); folded ring dim 32 = referee count |
| 12 | d = 9 sectors | excl chain 3 | sgn⊗sgn both parities closed (1 unknown each, V = {0}); triv/sgnX/sgnY running at hand-off | — |
| 13 | **"no τ-semi-invariant (pure-parity) denominator of degree 7"** (full family, 61 unknowns) | draft post item (2) | **UNSUPPORTED** | both msolve runs failed: rc 124 (timeout 3601 s) and rc −9 (killed, 3348 s); `.lead` files empty; mixed parity (124 unknowns) never run |
| 14 | "screened all character sectors at degrees 9–12" | draft post | UNSUPPORTED | no such receipts in `wave2-excl.log.jsonl` (sector-screen job not finished) |
| 15 | exact Farkas certificates: `D7`, `D7Y`, `E_X(F)`, `E_Y(F)`, `E_X(F)−E_Y(F)`, `Δ_XΔ_Y` have no certificate over ℂ | excl `farkas-exact` | VERIFIED (D7, D7Y, E_X(F), Δ_XΔ_Y replayed; E_Y(F), E_X−E_Y see `referee/farkas_*.replay.txt`) | `verify_farkas.py`: unfolded `y` kills every column of the full fixed-D system (4497 resp. 2602 columns), `y·b ≠ 0`, exact integers ⇒ inconsistent over any field ⊇ ℚ |
| 16 | all 12 degree-8 E1-kernel basis elements inconsistent as single denominators | aggregator batch 2b (10⁹+7, compressed solver) | PLAUSIBLE | one prime, one implementation; modular screen of basis elements, not of the family |
| 17 | joint varpro search: nonzero stable minima in every sector searched (triv d ≤ 10, jac d ≤ 6, sign sectors, vertex cones); 358 kernel classes + 90 random combinations inconsistent mod p | search report §7–9 | PLAUSIBLE (screen) | numerical / modular; restart counts recorded; not a theorem |
| 18 | reduction `[A] ∈ served(Δ_X^aΔ_Y^b S) ⇔ [A] ∈ served(S)` for G-semi-invariant S | search report §3 | PLAUSIBLE (conditional) | referee checked the group-theoretic step: for χ = sgn⊗sgn and any transposition s, `e_χ(H⁴(X)^s) = 0` since `e_χ v = e_χ(sv) = χ(s)e_χ v = −e_χ v`; relies on the book's prescribed-support theorem and sum principle (not re-proved here). Consequence consistent with claim 15: `D7` is excluded already by the degree-4 E1 failure of `e2(X²)` |
| 19 | monomial automorphism group exactly 36 × ⟨τ⟩; age screen (min age 1) | dossier, excl reverify | PLAUSIBLE (not replayed) | |
| 20 | wave-1 log artefacts (`6F` in the served scan; `targets=True` rows) are not results | excl §2.2, search §0 | PLAUSIBLE | two lenses agree; referee did not read the wave-1 log |
| 21 | folded-reducer ranks = Molien–Koszul values (ring 46 at degree 14, 32 at degree 16) | excl receipts | VERIFIED | referee Burnside: invariant Y-odd R_14 = 46, R_16 = 32 (`referee/burnside_rings.out`) |
| 22 | mod-p ⇒ characteristic-0 transfer for the relaxation exclusions | excl §1 "rigor conventions" | VERIFIED with one addition (see §2) | needs Koszul exactness `H_1(K(F̄))_{d+6} = 0` at the prime used, not only the two saturation ranks; checked by the referee (`referee/e1_replay.out`, `referee/koszul_p2147483647.out`) |

## 2. The transfer argument, re-derived

Certificate (formal contract, q = 1): `E = A D² − D Σ N_i F_i + F(D div N − N·∇D) = F² Q`, `D` coprime to `F` over ℂ.
Reducing mod F (F irreducible): `A D − N·∇F = F H`. With any `N⁰` such that `Σ N⁰_i F_i = A D` (E1), Euler gives
`Σ (N_i − N⁰_i + H z_i/6) F_i = 0`, so `N = N⁰ − (H/6) z + K_B` with `K_B` Koszul (`F_i` regular). Substituting,
the `H` terms cancel exactly (`1 − (d+6)/6 + d/6 = 0`) and `D div K_B − K_B·∇D ∈ J` (antisymmetry of `B`, symmetry of
`F_ij`), hence **`R⁰(D) := D div N⁰ − N⁰·∇D ∈ J(F)_{2d}`**, independent of the lift. The referee also checked
`R⁰(D + Fg) ≡ R⁰(D) mod J` for every lift (cross terms `A g D − g A D`), so the solution cone is invariant under the
radial subspace `F·S_{d−6}` and restricting to a linear complement loses nothing.

Transfer (ℂ → F_p). A certificate over ℂ gives one over a number field (constructible set defined over ℚ). Normalise `D`
𝔭-integral with `D̄ ≠ 0`; if `D̄ ∈ (F̄)` replace `D` by `(D − Fg)/π` (E1 and the relaxation are shift-invariant over the
number field) — Krull's intersection theorem ends this because `D ∉ (F)`. Then: (i) `rank_p J_{d+6} = rank_ℚ` ⇒
`(S/J)_{d+6} ⊗ O` torsion-free ⇒ an integral `N⁰` exists and `A D̄ ∈ J̄`; (ii) `rank_p J_{2d} = rank_ℚ` ⇒ `R̄⁰(D̄) ∈ J̄`
for the reduced lift; (iii) **the code uses its own lift**; two lifts differ by a degree-(d+1) syzygy of `F̄`, which moves
`R̄⁰` inside `J̄` only if that syzygy is Koszul, i.e. if `H_1(K(F̄))_{d+6} = 0` at the prime used. (iii) is implicit in the
lens's "rigor conventions"; it is implied by smoothness of `F̄` but is not implied by the two rank identities alone. The
referee computed `rank ∂_2` of the Koszul complex in degree `d+6` and compared it with `6·dim S_{d+1} − rank J_{d+6}`:
equality (`H_1 = 0`) holds in degrees 11, 12, 13, 14 at p = 10⁹+9 and in degrees 12, 13, 14 at p = 2³¹−1 (the lens's
prime; see `referee/koszul_p2147483647.out`). With (i)–(iii) the mod-p statement "only radial solutions" is a
characteristic-0 exclusion for that family and degree. For the sector families the same argument applies to the
χ-isotypic summands (|G| = 36 is a unit mod p; folded ranks checked against the Molien–Koszul values, claim 21).

## 3. Replays performed by the referee (all single-threaded, niced; machine load ≈ 42 on 40 cores)

| replay | input | result | file |
|---|---|---|---|
| Farkas `D7`, `D7Y`, `E_X(F)`, `Δ_XΔ_Y`, `E_Y(F)`, `E_X−E_Y` | packet pickles + `verify_farkas.py` | 0 violating columns, `y·b ≠ 0` | `referee/farkas_*.replay.txt` |
| pickled denominators = named polynomials | own construction | `D7`, `Δ_XΔ_Y`, `E_X(F)` equal up to scalar | `referee/farkas_D_check.txt` |
| sector systems d = 7, 8 (12 nonempty) | stored `.ms`, `msolve -g 1 -t 1` | leading ideals identical to stored `.lead`; pure power of every variable | `referee/replay_sectors.out` |
| Singular dim/vdim | stored `.sing` (d5 ×2, triv d7/even, triv d8/odd, sgnX d8/odd) | dim 0; vdim 4, 4, 11, 112, 7 | `referee/singular.txt` |
| stored d = 6 mixed system | `msolve -t 1` (95 s) | identical leading ideal; pure powers 35/35 | `referee/d6mixed_replay.lead` |
| own-code d = 6 relaxation | `relax6_own.py`, p = 10⁹+9, msolve 41 s | m-primary: pure powers 35/35 (exponents 2–4) | `referee/relax6_own_p1000000009.lead` |
| E1 kernels / Jacobian ranks d = 4..7, A ∉ J | `e1_replay.py`, p = 10⁹+9 | all ranks = CI; ker_R = 0,0,0,1 per parity | `referee/e1_replay.out` |
| Koszul `H_1 = 0` in degrees 11–14 | `e1_replay.py` (10⁹+9), `koszul_only.py` (2³¹−1) | see §2 | `referee/e1_replay.out`, `referee/koszul_p2147483647.out` |
| sector Hodge numbers; invariant ring dims | `burnside_sectors.py`, `burnside_rings.py` | (17,52,17) / (1,10,1); 46, 32 | `referee/burnside_*.out` |

## 4. Corrections required in the draft post (`wave2-symmetric-sextic-residual-draftpost.md`)

* "(sgn⊗sgn, τ=+1) isotypic sector (Hodge numbers (1,10,1))" → **(17, 52, 17)**; [A] is one of 17 (3,1)-classes of its
  residue character, not "the unique" one; drop the K3-type / rank-12 language.
* item (2) "no τ-semi-invariant denominator of degree 7" → **not established** (claim 13); item (3) → degree 8 holds,
  degree 9 only for the sgn⊗sgn sectors so far; the "screened at degrees 9–12" sentence → drop (claim 14).
* "leading ideal contains `a_k²` for all k" → "a pure power (exponent ≤ 4) of every variable".
* the rigor paragraph should name the Koszul-exactness condition (§2 (iii)) next to the two rank identities.
* `D7`/`D7Y`: the Farkas certificates are exact; they are also implied by the degree-4 E1 failure via claim 18.

## 5. Open items

Full degree-7 relaxation (pure parity 61 unknowns; mixed 124 unknowns) — Gröbner computation did not finish; d = 9 sectors
triv/sgnX/sgnY; d ≥ 8 non-semi-invariant denominators (the affine cosets `D_sporadic + J_8`); exact (characteristic-0)
Gröbner versions of the sector closures; the conditional reduction of claim 18 as a stated theorem with its book citations.
