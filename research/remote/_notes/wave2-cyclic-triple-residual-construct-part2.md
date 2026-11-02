# Part 2 — Exchange covers (book §2.7 shape) for F: all linear directions, rigorous negative

Scripts/receipts: /tmp/hodge/results/cyclic-triple-residual/wave2-construct/c1_exchange_dirs.{py,json}, c1b_fourth_power.ms (msolve),
c1b_elim.sing (Singular, char 0).  Engine not used.

## 2.1 What the §2.7 mechanism needs
F = y^6 + y·Λ + f with (after adjoining T, T^6 = Φ6(u_S)) Λ = T·Λ'(ratios, X_u) and f = T^6 + f_X(X_u): the exchange
y <-> T is then a *linear* involution of the model Z' ⊂ P^{k} × (curve {Φ6 = 1}) fixing the divisor {y = T}, and
Z'/⟨exchange, τ⟩ is a Fano fibration over the curve (holdout: fibres sextics in P(1,1,1,2,2), K = O(-1)).  Three
combinatorial requirements on F, in a coordinate system where y is a coordinate: (a) exponent set of y is {0,1,6};
(b) every monomial of Λ has S-degree exactly 1 (S = the T-block of coordinates); (c) every y-free monomial has S-degree
0 or 6, and the S-degree-6 part is Φ6 (then T^6 := Φ6).

## 2.2 Which linear forms y have exponent pattern {0,1,6}?  (coordinate-free, exact)
F has the shape c·ℓ^6 + ℓ·Λ + f in direction v (ℓ(v) = 1) iff D_v^2 F = const·ℓ^4 for a linear form ℓ.  The ideal of
{(v, a) : D_v^2 F = (a·X)^4} (126 equations, 12 unknowns) was eliminated to the v-variables with Singular over Q
(c1b_elim.sing): the v-variety has exactly two components, v ∝ e_{X1} and v ∝ e_{Y0}.  (msolve on the same system:
positive-dimensional, c1b_fourth_power.out; D_v^2 F ≡ 0 has only v = 0, c1_exchange_dirs.json.)
Exponent sets of F per coordinate (exact): X0:{0,1,2,3,6} X1:{0,1,6} X2:{0,3,6} Y0:{0,1,6} Y1:{0,1,3,6} Y2:{0,1,2,6}.
So the only candidates for the exchanged variable are y = X1 and y = Y0 (and only X2 for the cubic variant y^6+y^3Γ+Φ).

## 2.3 y = X1.  Λ1 = t·Y1·(X0^3 Y0 + X0 Y1^2 Y2 + X2^3 Y0).  Requirement (b) forces S ∌ X0, X2, Y1 and S ∋ exactly one of
{Y0,Y1} and one of {Y2}: S = {Y0, Y2}; then Φ6 = Y0^6 + Y2^6 (Fermat curve base, as in holdout) — but the X1-free
couplings t·X0^2 Y0 Y1 Y2^2 (S-degree 3) and t·X0 X2^3 Y1 Y2 (S-degree 1) violate (c).  Under Y0 = αT, Y2 = βT these
become t X0^2 Y1 αβ^2 T^3 + t X2^3 X0 Y1 β T, i.e. F = X1^6 + T^6 + X1 T Λ' + (X1,T-free) + (T^3 and T^1 terms) — not
symmetric under X1 <-> T, and no polynomial change of the cover equation repairs a T^3 + T^1 asymmetry.

## 2.4 y = Y0.  Λ0 = t·Y1·(X0^3 X1 + X0^2 Y2^2 + X2^3 X1).  (b) forces S = {Y1}; then Φ6 = Y1^6 but the Y0-free
couplings t·X0 X1 Y1^3 Y2 (S-degree 3) and t·X0 X2^3 Y1 Y2 (S-degree 1) violate (c).  Same failure mode.

## 2.5 y = X2 with the cubic pattern (x^6 + x^3 Γ + Φ, exchange x <-> T with T^3 ~ Γ or T^6 ~ Φ).  The x^3-coefficient
tC = t(X0^3 + X0 Y1 Y2 + X1 Y0 Y1) involves all five other variables, so S = {X0,X1,Y0,Y1,Y2}: T^6 = H, the base is the
sextic fourfold B̂ = {w^6 = H} ⊂ P^5 (κ = 0, a CY sextic), Z' = {x^6 + T^6 + t x^3 T^3 C(ũ)/w^3 = 0} ⊂ P^1 × B̂ is a
6:1 cover of B̂ and Z'/⟨x <-> T⟩ a 3:1 cover of B̂: not uniruled (κ ≥ 0).  The exchange exists but its quotient has no
rational curves to integrate along.  The same holds for Z = X2^3 <-> T with T^2 = H (base = the Fano double cover
{v^2 = H} of P^4, but then the exchange Z <-> T requires C(u) to transform by (Z/T)^{1/3}: not algebraic).

## 2.6 Verdict
No exchange cover of §2.7 type exists for F, for any linear form: the obstruction is the pair of couplings
X0^2 Y0 Y1 Y2^2 and X0 X2^3 Y1 Y2 (for y = X1) or X0 X1 Y1^3 Y2 and X0 X2^3 Y1 Y2 (for y = Y0) — exactly the
Y-degree-4 and Z-coupled terms that are absent from the chain/holdout sextics.  (Book boundary 8 of §2.7: "a perturbation
that destroys the exchange … is a different failure that a new conic cannot repair".)
