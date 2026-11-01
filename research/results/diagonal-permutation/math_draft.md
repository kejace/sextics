## Mathematical notes (draft, resumed session)

### Hodge-theoretic meaning of the problem
`A = X0^3Y0^3 - X1^3Y1^3 ∈ R_6^-` is a `τ`-odd primitive `(3,1)`-class `ω_A = Res(A Ω/F^2)` (q = 1 ⇒ degree 6 ⇒ `H^{3,1}`).
The `τ`-odd part `H^4(X)^-` has no `(4,0)` part (the 4-form is `τ`-even), so it is a weight-4 Hodge structure of level 2 and the
generalized Hodge conjecture predicts a divisor `Y ⊂ X` with `ω_A ∈ Gysin(H^2(Ỹ)(-1))`.  A certificate `(D,N,Q)` realises this with
`Y = X ∩ V(D)`: it exhibits `A Ω/F^2` as `d(η)` + (form with pole along `D`) ... (book §3.1).  Because `H^{2,0}(Y) = 0` for every effective
divisor on `X` with rational singularities, the support `Y_D` must carry non-rational singularities: cone points of multiplicity `>= 4`
with surface base of `p_g > 0`, curves of singularities with non-rational transverse type (simple-elliptic `x^3+y^3+z^3`, `x^4+y^4+z^4`
with genus-3 base, ...), or non-normal loci giving nontrivial covers.

### Gorenstein duality for the E1 kernel
`R = S/J(F)` is Gorenstein with socle degree 24, and the socle (Hessian class) is Y-even, so the pairing respects Y-parity:
`R_d^ε × R_{24-d}^ε → R_24`.  Multiplication by the Y-odd `A` therefore satisfies
`rank(A : R_d^ε → R_{d+6}^{1-ε}) = rank(A : R_{18-d}^{1-ε} → R_{24-d}^{ε})`.  With `k_d^ε := dim Ann(A)_d^ε` this gives
`k_d^ε - k_{18-d}^{1-ε} = dim R_d^ε - dim R_{18-d}^{1-ε}`.  Consequences (exact, using the computed `k_8^+ = 1`, `k_d = 0` for `d <= 7`):
* `k_10^- = 744 - 482 = 262` (forced lower bound 261 + 1 from `D8`), `k_10^+ >= 294` (forced), `k_11^± >= 510`, `k_12^+ >= 677`.
* `k_9^+ = k_9^-` (the two parities at `d = 9` are dual to each other: `623 -> 623`), `k_8^- = k_10^+ - 294 = 0 ⇒ k_10^+ = 294` exactly;
  `k_7 = 0 ⇒ k_11^± = 510` exactly; `k_6^+ = 0 ⇒ k_12^- = 648`, `k_6^- = 0 ⇒ k_12^+ = 677` exactly, and so on.
So the complete E1 picture is known in every degree except `d = 9` (which is the only degree where the kernel is not determined by
dimension counting plus the `d <= 8` computation): **the E1-admissible denominators modulo `J(F)` are exactly: nothing for
`d <= 7`; one line (`D8`) at `d = 8` (Y-even), none Y-odd; `k_9^+ = k_9^-` unknown; `294 / 262` at `d = 10`; `510/510` at `d = 11`; ...**

### Fermat ancestry of `D8`
At the Fermat point `F_0 = ΣX_i^6 + ΣY_i^6` the Jacobian ring is monomial (`x_i^5 = 0`) and `Ann(A)_8^+` is spanned by the Y-even
degree-8 monomials lying in `(X0^2, Y0^2) ∩ (X1^2, Y1^2)` (they kill both `T_0` and `T_1`), e.g. `X0^2X1^2Y0^2Y1^2`, plus elements
like `X0^2Y0^2·m` with `m ∈ (X1^2,Y1^2)` — a large space.  Along the deformation `F_t = F_0 + t·5G` the annihilator drops to the
single line `D8` at `t = 1` (the family scan below shows whether this is generic for the 3-parameter family `Fermat + aO1 + bO2 + cO3`,
for the whole invariant family, or special to the instance).

### Vector-field reformulation for Jacobian-ideal denominators
For `q = 1` the identity forces `N(F) = A D + F H` (`deg H = d`), and then the certificate condition is exactly
`F | D·div(N) - N(D) - D·H`.  For `D = G(F) = Σ G_i ∂_i F` (the only E1-admissible shape for `5 <= d <= 8`, `D ≠ D8`),
one particular solution is `N^0 = A·G`, `H^0 = 0`, and the condition for the general `N = A G + K + F·M` (`K` a Koszul syzygy,
`M` arbitrary of degree `d-5`) reads `F | A·G(F)·div(G) + G(F)·G(A) - A·G^2(F) + (Koszul/M corrections)`.
Geometrically `Y_{G(F)}` is the tangency locus of the vector field `G` with `X`; for `G = ∂_i` it is the polar section, whose
singularities are isolated and rational, hence it cannot support a `(3,1)` class — this is confirmed by the served-target
computation (`served_dim = 18 = dim J_6^-`, i.e. only the zero class is served) for all six partials and their `S3`-combinations.
