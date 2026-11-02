# Referee report — `symmetric-sextic-residual` (#117), wave-2 exclusion lane

Prepared for `kejace` with Claude (Anthropic) assistance; no human review yet; nothing Lean-verified; no certificate claimed.
Date: 2026-11-01 (local, brunos). Commit for all links: `1c42238b88abc89492f5d73698b675d1cfb9b6cd` (research branch of `kejace`'s fork).

## Counts

VERIFIED 14 · PLAUSIBLE 5 · UNSUPPORTED 2 · REFUTED 1 (the dossier's sector table; corrected value verified).

Strongest verified statement: **every legal certificate denominator has degree ≥ 7 (any Y-parity, Jacobian-ideal elements
included), and no `S3×S3×⟨τ⟩`-semi-invariant denominator of degree 7 or 8 exists** — established by E1 (Jacobian ranks at
complete-intersection values at four primes), the Hessian relaxation (necessary condition `R⁰(D) ∈ J(F)_{2d}`, re-derived
by the referee), m-primary quadric systems certified by msolve leading ideals (replayed; the degree-6 case rebuilt from
scratch by the referee at a third prime), and the transfer conditions (saturation ranks + Koszul exactness `H_1 = 0` in
degrees 11–14, checked at both primes used).

## Verdicts (full table with evidence: [README-referee.md](https://github.com/kejace/sextics/blob/1c42238b88abc89492f5d73698b675d1cfb9b6cd/research/results/symmetric-sextic-residual/wave2-excl/README-referee.md))

* VERIFIED: F/A/q structure; smoothness; `A ∉ J`; residue character sgn⊗sgn (τ = +); **sector Hodge numbers (17,52,17)**;
  E1 kernel zero for d ≤ 6; d = 7 kernel 1 + 1 (`D7`, `D7Y`); d = 8 kernel 6 + 6 (two implementations, not replayed);
  **deg D ≥ 7**; **no semi-invariant D of degree 7**; **none of degree 8**; d = 9 sgn⊗sgn sectors closed; six exact Farkas
  certificates (`D7`, `D7Y`, `E_X(F)`, `E_Y(F)`, `E_X−E_Y`, `Δ_XΔ_Y`); folded-reducer ring dimensions 46/32; the transfer
  argument (with the Koszul-exactness condition made explicit and checked).
* PLAUSIBLE (screens or conditional): the 12 degree-8 kernel basis elements inconsistent (one prime, one solver); the joint
  varpro landscape (nonzero stable minima; 358 + 90 modular single-D inconsistencies); the reduction
  `[A] ∈ served(Δ_X^aΔ_Y^bS) ⇔ [A] ∈ served(S)` (group-theoretic step checked; depends on the book's prescribed-support
  theorem and sum principle); monomial automorphism group / age screen (not replayed); wave-1 log artefacts.
* UNSUPPORTED: "no τ-semi-invariant (pure-parity) denominator of degree 7" from the full 61-unknown relaxation (both msolve
  runs timed out / killed, empty `.lead`); "sectors screened at degrees 9–12" (no receipts).
* REFUTED: dossier §1.4 "(sgn,sgn,+) has Hodge numbers (1,10,1); [A] is the unique (3,1)-class of its character". Correct:
  (17,52,17); (1,10,1) is the residue-trivial sector spanned by Δ_XΔ_Y. The draft post carried the wrong numbers and has
  been corrected in the final comment.

## Referee replays (single-threaded, niced; raw files under `wave2-excl/referee/`)

| what | result |
|---|---|
| Farkas certificates ×6 (`verify_farkas.py`, exact integers) | 0 violating columns (4497 / 2602), `y·b ≠ 0`; pickled D equal the named polynomials |
| 12 sector systems d = 7, 8 and 2 sgn⊗sgn systems d = 9 (msolve `-g 1 -t 1`) | leading ideals identical to the stored ones; pure power of every variable |
| Singular on d = 5 (×2), triv d7 even, triv d8 odd, sgnX d8 odd | dim 0; vdim 4, 4, 11, 112, 7 |
| stored d = 6 mixed system (msolve, 95 s) | identical leading ideal; pure powers 35/35 |
| own-code d = 6 relaxation at p = 10⁹+9 (`relax6_own.py`; 1684 quadrics / 35 unknowns; msolve 41 s) | m-primary (pure powers, exponents 2–4) |
| own-code E1 kernels and Jacobian ranks d = 4..7, `A ∉ J`, p = 10⁹+9 | all ranks = CI; ker_R = 0, 0, 0, 1 per parity |
| Koszul `H_1 = 0` in degrees 11–14 (p = 10⁹+9) and 12–14 (p = 2³¹−1) | all true |
| Burnside recount of the sector table; invariant ring dims 14/16 | (17,52,17) / (1,10,1); 46, 32 |

## Transfer argument (mod p → characteristic 0)

Derivation and the three conditions are in README-referee.md §2. The addition relative to the lens's "rigor conventions":
the code's E1 lift differs from the reduction of an integral characteristic-0 lift by a degree-(d+1) syzygy of `F̄`; this
moves `R̄⁰` inside `J̄` only if the syzygy is Koszul, i.e. if `H_1(K(F̄))_{d+6} = 0` at that prime. This holds for smooth
`F̄` but is not implied by the two saturation ranks; it was checked directly (rank of the Koszul `∂_2` versus
`6·dim S_{d+1} − rank J_{d+6}`).

## Corrections applied to the draft post

Sector numbers; items (2) and (3) downgraded (degree 7 full family not established; degree 9 only sgn⊗sgn); the
"screened 9–12" sentence removed; "a_k² for all k" → "a pure power of every variable"; Koszul condition named; credit and
AI-assistance lines added; no placeholders.

## Open items

Full degree-7 relaxation (pure and mixed parity); degree-9 triv/sgnX/sgnY sectors (running at hand-off); non-semi-invariant
denominators of degree ≥ 8 (affine cosets `D_sporadic + J_8`); exact characteristic-0 Gröbner versions; the conditional
reduction as a cited theorem; any construction for the 17-dimensional sector.

Links: [packet](https://github.com/kejace/sextics/tree/1c42238b88abc89492f5d73698b675d1cfb9b6cd/research/results/symmetric-sextic-residual/wave2-excl) ·
[exclusion report](https://github.com/kejace/sextics/blob/1c42238b88abc89492f5d73698b675d1cfb9b6cd/research/notes/wave2-symmetric-sextic-residual-excl.md) ·
[search report](https://github.com/kejace/sextics/blob/1c42238b88abc89492f5d73698b675d1cfb9b6cd/research/notes/wave2-symmetric-sextic-residual-search.md) ·
[final comment](comment.md).
