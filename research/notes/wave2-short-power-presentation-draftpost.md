# Draft discussion comment — short-power-presentation (wave 2, exclusion/rigor lens) — DRAFT, not posted

### Low-degree denominators are excluded; the E1 kernel starts in degree 7 with a pullback from the curve-fibration base

Everything below is for the exact record `F = eX + eY + c((L+M)^6+(L−M)^6)`, `c = 5/7^7`, `A = X0^3 Y0^3`, `q = 1`,
with τ the global Y-sign. "Proven" means an exact rational computation or a modular computation with the explicit
rank argument (rank_p ≤ rank_Q, plus the Hilbert-series check that the prime is good in the degrees used);
"screened" means modular evidence only.

**1. E1 kernel (proven).** Writing `K_d^± = {D ∈ S_d^± : A·D ∈ J(F)}`, the quotient `K_d^±/J(F)_d^±` is zero for
`d ≤ 6` in both parities and is one-dimensional for `d = 7` in both parities (primes 2^31−1, 2147483629, 1000000007 and
1073741789 agree, each with the Hilbert-series check that the prime is good in the degrees used; the degree-7 classes
have exact rational representatives `D7^±` with exact numerator certificates `N(F) = A·D7` verified over Q). The even class `D7^+` is a polynomial in `X0, X1, X2, M` only (M-degrees 0, 2, 4), i.e. a cone with
vertex the line `{X = 0, M = 0}`; the odd class `D7^-` is the mirror polynomial in `Y0, Y1, Y2, L`. Geometrically
`V(D7^+) ∩ X` is the pullback of a degree-7 surface from the base `P^3_{(X0:X1:X2:M)}` of the projection of `X`
from that line, whose fibres are the plane sextics `D_s` of benediamond's model. [TODO: K_8, K_9 generation by K_7.]

**2. No denominator of degree ≤ 5, no Y-pure denominator of degree 6 (proven).** For such `D`, E1 forces
`D = w(F)` for a polynomial vector field `w`; the certificate identity then forces `A·w^T Hess(F) w ∈ J(F)_{2d}`
(after normalising `N = A w + (G/6)E + K`, `K` Koszul), a homogeneous quadratic system in the coefficients of `w`
invariant under `w ↦ w + G·E` (which changes `D` by `6GF`). Modulo the Euler directions the ideal of this system is
m-primary mod 1073741789 (Hilbert functions `[1,3,2,0]` for d = 5; `[1,17,114,358,501,225,0]` for d = 6 even;
vanishing at t = 6 for d = 6 odd, Singular `vdim` 540; mixed d = 5: `[1,6,13,8,0]`), and since the Hilbert function
of the mod-p ideal bounds the rational one from above, the only solutions over Q-bar are the Euler multiples, i.e.
`D ∝ F`. Both msolve (F4, leading ideal) and Singular (`std`, `hilb`) were used and agree.
For mixed parity `D = D^+ + D^-` the odd part of the relaxed condition is `Ξ^+(w^+) + Ξ^-(w^-) ∈ J` with no cross
terms; every linear functional vanishing on `J` and on the span of the coordinate polynomials of `Ξ^+` kills `Ξ^-(w^-)`,
which gives 36 quadrics in the 18 odd coefficients forming an m-primary ideal (Singular `vdim` 4096); the functionals
are reductions of rational ones (114 exact dependency certificates, verified with exact arithmetic), so `w^- = 0` and
then `w^+ = 0` by the Y-even case. **Hence every certificate denominator has degree at least 7.**

**3. Degree 7, Y-even (partly proven).** Here `D = λ·D7^+ + w(F)`. Block elimination by the same functional trick
(the Hessian of a sum of eight sixth powers has rank ≤ 8, so the six component blocks of `w` are coupled only through
two small cross spans) shows mod 1073741789 that `λ = 0` — no degree-7 Y-even denominator involves the kernel class —
and that the 18 coefficients `w_{k, X_i Y_j}` (`j ≠ k`) of the Y-components vanish; for Y-odd `D = λ·D7^- + w(F)`
the mirror statement holds (`λ = 0`, the X-component blocks vanish). What remains at degree 7 is a residual quadratic
system in 42 coefficients (231 independent quadrics; every random slice of codimension ≤ 11 is empty; msolve's F4 passes
degree 5 with a 959426 × 685669 matrix) — its Gröbner basis and the exact rational certificates for the eliminating
steps are still running when this is written. The kernel classes themselves are excluded at degree 8 too: the relaxed
defect ideal of `D ∈ span(K_8^±/J)` is m-primary (vdim 7, both parities), and `K_8^- = S_1·K_7`, `K_8^+ ⊋ S_1·K_7`
(two new even classes), `K_9^+/J` has dimension 18 with `S_2·K_7` spanning 15.

**4. What is screened only.** The full fixed-`D` system for `D7^+` is inconsistent mod p (search lens), but its
homogeneous solution space is 6 larger than the explicitly known one, so this single test is modular evidence; the
degree-7 relaxation covers `D7^+ + J_7` entirely once it finishes.

**5. Obstacle at degree 8.** The relaxation at `d = 8` has ~160 unknowns and ~470 independent quadrics; the systems
here are far from semi-regular (at d = 6 the Hilbert function vanishes only in degree 6 where a generic system would
vanish in degree 4), so Gröbner certification of emptiness at d = 8 is out of reach with this formulation.
Reproduction: scripts and JSON receipts in the lane directory (`excl_*.py`, `wave2-excl/README.md`,
`wave2-excl.log.jsonl`). No claim of nonexistence of certificates is made; all statements are degree bounds.
