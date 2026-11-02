# Draft discussion comment — `diagonal-permutation` (#177) — exclusion data, wave 2

(DRAFT; numbers marked [pending] are being filled in from the remote receipts. Do not post before the receipts are complete.)

**Exclusion data for `diagonal-permutation`: no certificate with `deg D ≤ 6` (any Y-parity), and the degree-7/8 picture**

Provenance: computations by `kejace` (AI-agent assisted, modular linear algebra with python-flint, Groebner bases with msolve, no independent human review claimed). Reproduction packet: `/tmp/hodge/results/diagonal-permutation/wave2-excl/` (scripts `w2lib.py`, `w2_relax.py`, `w2_e1kernels.py`, `w2_d8.py`; JSON receipts with primes and ranks; msolve inputs/outputs) — to be mirrored to the research branch.

Setup. `q = 1`; a certificate is `(D, N, Q)` with `A D² − D N(F) + F (D div N − N(D)) = F² Q`. Reducing mod F gives `N(F) = A D + F H` (E1: `A·D ∈ J(F)`), and since the partials form a regular sequence, `N = u + K + (H/6)·Euler` with `u(F) = A D` and `K` Koszul. The identity becomes `F·(Ψ(D;u) + D div K − K(D)) = F² Q` with `Ψ(D;u) := D div u − u(D)`. **Lemma.** `D div K − K(D) ∈ J(F)` for every Koszul field, so every certificate satisfies `Ψ(D;u) ∈ J(F)` ("E2 modulo J"); `Ψ mod J` is independent of the choice of `u` and of `D` modulo `F·S`. For Jacobian-type `D = w(F)` one gets `Ψ ≡ −A · wᵀ Hess(F) w (mod J)`.

What is proven (modular ranks are lower bounds for rational ranks; zero kernels mod p transfer to Q; a zero-dimensional reduced Groebner basis mod p transfers to Q by the Hilbert-function argument; good reduction of `J(F)` in every degree used is checked):

1. `d ≤ 4`: `ker([A]·: R_d → R_{d+6}) = 0` in both parities (primes 2147483647 [wave 1], 1000000007, 998244353) and `J_d = 0`, so no denominator of degree ≤ 4.
2. `d = 5, 6, 7`: `ker [A] = 0` in both parities, so every admissible `D` is `w(F)` with `deg w_i = d − 5`; the relaxation `[A]·[wᵀ Hess(F) w] = 0` is a homogeneous quadratic system in `w` (Euler multiples `gE`, i.e. `D = 6gF`, removed by a linear complement).
   * `d = 5`: the 21 relaxation equations in 6 unknowns span all quadratic monomials ⇒ `w = 0` (equivalently: `∂_w² F ≡ 0` would make `[w]` a point of multiplicity ≥ 5 on `X`).
   * `d = 6`: msolve reduced Groebner bases (p = 1000000007) are zero-dimensional in the Y-even (17 unknowns, 134 equations), Y-odd (18 / 153) and mixed-parity (35 / 404) sectors ⇒ no degree-6 denominator of any parity.
   * `d = 7`: [pending: relax7F job — even (60 / ~482), odd (60 / ~468), mixed (120 / ~950)].
3. `d = 8`: `ker [A]` on `R_8^+` is exactly one line `[D8]` (414-term rational representative; `A·D8 ∈ J(F)` verified over Q with an exact `u8` — [pending d8F]); `ker` on `R_8^-` is 0. For `D8` itself: [pending: relaxation `Ψ(D8;u8) ∈ J`? reduced E2' system at three primes with rank certification]. The full family `λ D8 + w(F)` (141 unknowns, 468 quadratic equations) is [pending/screened].
4. `d = 9`: `ker [A]` is 21-dimensional in each parity (sibling computation at 2147483647, re-verified at 1000000007 [pending]); so E1 alone stops excluding at degree 9, and the Y-even/Y-odd kernels split into (01)-sectors as [pending].

What is screened only (modular, no rank certification, or random members of a family): [pending list — structured products of degree 8–12 from the dossier list; random members of the degree-9 kernel; the degree-8 family].

Consequences via the sum principle (classes supported on `V(D1 D2)` are sums of classes supported on the components; `served(D)/J_6 ⊆ Ann_6(D) := {A' ∈ R_6^- : A' D ∈ J}`): for all 49 structured non-Jacobian factors tested (`p2X, p2Y, pXY, e2X, e2Y, X0²±Y0², X0X1±Y0Y1, (X0−X1)²±(Y0−Y1)², X0Y1±X1Y0, X2², X0X1, X2Y2, p3X, e3X, Vandermondes, sXY2, sX2Y, XYYcyc, Hesse, p3Y, e3Y, p2X², p2X p2Y, pXY², p4X, sX2Y2, p4Y, X0⁴−X0²Y0²+Y0⁴, X0⁴−X1²Y0², Vdm_X², Vdm_Y², Vdm_X Vdm_Y, ∏((Xi−Xj)²±(Yi−Yj)²), p6X, p6Y, ΣT_i, B2, A, pXY³, O1, O2, O3, F_c`) one has `Ann_6(D) = 0` (mod 1000000007; an upper bound, hence rigorous), so every product of these factors, in every degree, serves only the zero class. The same holds for the Jacobian-type denominators `w(F)` with `w` a partial, the X-Euler field, a pair-radial field, the `std`-radial field `g_A`, the swap/rotation/cyclic-shift fields: `Ann_6([wᵀ Hess(F) w]) = 0`.

Equivariant structure. The graded characters of `R` under `⟨τ⟩ × ⟨(01)⟩` give the sector dimensions; since `A` is Y-odd and `(01)`-anti-invariant, `[A]·` maps the sector `(par, σ)` of `R_d` to `(1−par, −σ)` of `R_{d+6}`. Dimension counting forces no kernel in any sector for `d ≤ 8`, and at `d = 9` forces exactly 21 in the `(01)`-invariant sector of each parity (`322 → 301`) and nothing in the anti-invariant one (`301 → 322`); the computed 21-dimensional degree-9 kernels are therefore exactly the forced ones and consist of `(01)`-invariant classes [sector computation pending]. At `d = 10` the forced sector kernels are `170 + 124 = 294` (Y-even, equal to the actual kernel) and `152 + 109 = 261` (Y-odd; actual 262).

None of this is a nonexistence claim; it says exactly which denominators are ruled out and how.
