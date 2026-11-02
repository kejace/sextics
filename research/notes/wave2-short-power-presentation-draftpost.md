# Draft discussion comment — short-power-presentation (wave 2, exclusion/rigor lens) — DRAFT, not posted

### Low-degree denominators are excluded; the E1 kernel starts in degree 7 with a pullback from the curve-fibration base

Everything below is for the exact record `F = eX + eY + c((L+M)^6+(L−M)^6)`, `c = 5/7^7`, `A = X0^3 Y0^3`, `q = 1`,
with τ the global Y-sign. "Proven" means an exact rational computation or a modular computation with the explicit
rank argument (rank_p ≤ rank_Q, plus the Hilbert-series check that the prime is good in the degrees used);
"screened" means modular evidence only.

**1. E1 kernel (proven).** Writing `K_d^± = {D ∈ S_d^± : A·D ∈ J(F)}`, the quotient `K_d^±/J(F)_d^±` is zero for
`d ≤ 6` in both parities and is one-dimensional for `d = 7` in both parities (primes 2^31−1, 2147483629 and
1073741789 agree; the degree-7 classes have exact rational representatives, checked by exact numerator
certificates). The even class `D7^+` is a polynomial in `X0, X1, X2, M` only (M-degrees 0, 2, 4), i.e. a cone with
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
[TODO: mixed parity d = 6 — 35 unknowns, 132 independent quadrics — pending; Y-pure d = 7 with the kernel class
as an extra unknown — pending.]

**3. What is screened only.** The full fixed-`D` system for `D7^+` is inconsistent mod p (search lens), but its
homogeneous solution space is 6 larger than the explicitly known one, so this single test is modular evidence; the
degree-7 relaxation covers `D7^+ + J_7` entirely once it finishes.

**4. Obstacle at degree 8.** The relaxation at `d = 8` has ~160 unknowns and ~470 independent quadrics; the systems
here are far from semi-regular (at d = 6 the Hilbert function vanishes only in degree 6 where a generic system would
vanish in degree 4), so Gröbner certification of emptiness at d = 8 is out of reach with this formulation.
Reproduction: scripts and JSON receipts in the lane directory (`excl_*.py`, `wave2-excl/README.md`,
`wave2-excl.log.jsonl`). No claim of nonexistence of certificates is made; all statements are degree bounds.
