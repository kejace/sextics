# Draft discussion comment — `symmetric-sextic-residual` (#117) — exclusions and rigorous degree bounds (wave 2)

*(draft; numbers marked [pending] are being recomputed on the remote and will be filled in before posting; nothing is posted yet)*

**Summary.** No witness. Rigorous (characteristic-0) statements: (1) every certificate denominator has degree ≥ 7, with the
Jacobian-ideal denominators `J_5`, `J_6` now excluded as well; (2) [pending] no τ-semi-invariant (pure Y-parity) denominator of
degree 7; (3) [pending] no `S3×S3×⟨τ⟩`-semi-invariant denominator of degree 8 (and 9 in the sectors where the computation
closes). Screened (mod p, not proofs): [pending] all character sectors at degrees 9–12 at random and structured points.
Exact rational inconsistency certificates (Farkas vectors) are provided for the named denominators `D7 = Δ_X·e2(X0²,X1²,X2²)`,
`D7Y`, `E_X(F)`, `E_Y(F)`, … [pending]. No claim of nonexistence is made.

**Setup.** `F = p6(X)+p6(Y)+(5/7)·e2(X)e2(Y)(e2(X)+2e2(Y))`, `A = e3(X)e3(Y)`, q = 1. The residue class of A spans the
unique (3,1)-line of the `(sgn⊗sgn, τ=+1)` isotypic sector (Hodge numbers (1,10,1)); this and the rest of the isotypic
table were recomputed by a Molien–Koszul formula `dim R_n^{χ,±} = (1/72) Σ_g ±χ(g)·[t^n] det(1−t⁵g)/det(1−tg)`.

**Method: Hessian relaxation.** For an E1-admissible D (`A D ∈ J(F)`) pick any `N⁰` with `Σ N⁰_i F_i = A D`. A certificate
forces `R⁰(D) := D·div N⁰ − N⁰·∇D ∈ J(F)` (the gauge and Koszul freedom of N only moves `R⁰` inside J). On a linear family of
D this is a system of quadrics; its solutions form a cone containing the illegal F-multiples. We compute the family of all
E1-admissible D (modular nullspace, ranks checked against the complete-intersection Hilbert function so that the modular
dimension bounds the rational one), pass to a complement of `F·S_{d−6}`, and ask msolve for the grevlex leading ideal: a pure
power of every variable certifies that only a = 0 solves it. A Q̄-denominator can be normalized to be p-integral with
reduction not divisible by F̄; torsion-freeness of `S/J` at p in the two relevant degrees (the rank checks) then transports
the certificate condition to the reduction, so the mod-p result is a characteristic-0 exclusion.

**Results.**
* d ≤ 4: nothing is E1-admissible (`J_d = 0`, kernel 0 at p = 2³¹−1 and p = 2147483629 [pending second prime]).
* d = 5: admissible = `J_5 = {c·∇F}`; relaxation gives `A·(cᵀH_F c) ∈ J_10`, impossible for c ≠ 0 (E1-injectivity in
  degree 4 and smoothness: `cᵀH_F c ≡ 0` would make c a singular point). Also Singular: dim 0 / vdim 4.
* d = 6: admissible = `J_6` (36-dim, incl. the radial line D = 6F); relaxation: 1605 quadrics in 35 unknowns, leading ideal
  contains `a_k²` for all k. ⇒ **deg D ≥ 7 for every legal denominator.**
* d = 7: admissible = `⟨D7⟩ ⊕ J_7` (Y-even) and `⟨D7Y⟩ ⊕ J_7` (Y-odd); 61 unknowns each … [pending]
* sectors d = 8, 9 … [pending]
* screening d = 9–12 … [pending]

**Reproduction.** Scripts and receipts: `wave2-excl/` packet (library `w2x_lib.py`, msolve inputs `*.ms`, leading ideals
`*.lead`, JSON receipts with primes and ranks, Farkas pickles + `verify_farkas.py`). [link to be added]

**Caveats.** A modular `inconsistent` for a single D is a screen unless accompanied by the Farkas vector; the relaxation
exclusions are rigorous as stated (they rely on msolve's Groebner computation mod p and on the rank identities recorded in the
receipts). Nothing here bounds the degree of a *geometric* (non-explicit) support; the sum principle only reduces to
orbit products of irreducible factors, so "no semi-invariant D of degree ≤ k" gives "no irreducible P with deg(orbit product) ≤ k".
