# wave2 dense-rational CONSTRUCT — part 4: what a mechanism must produce, audit of universal divisors, and the
# most promising candidate (the residual-tangency divisor R') with its explicit equations

## 4.1 The (2,0)-budget, sharpened: only irreducible D matter, and only the pure part of H_4(Y)
Let Y = X ∩ V(D) (reduced), Y_i its components, π_i: Y~_i -> Y_i resolutions.  The certificate exists iff res(ω_A)
lies in Im(⊕ H^2(Y~_i) -> H^4(X)) (Ch. 1).  Two refinements used below:
 (a) additivity: the budget Σ h^{2,0}(Y~_i) is a sum over components; products D = D_1 D_2 and powers add nothing
     beyond their factors.  So only IRREDUCIBLE denominators need to be designed.
 (b) the Gysin map factors through the lowest-weight part W_{-4} H_4(Y) = Im(H_4(Y~) -> H_4(Y)); classes of H^2(Y~)
     pulled back from a stratum of dimension <= 2 (exceptional divisors over curves or points, or over a surface
     through its projection) push forward to zero in H^4(X).  Hence (2,0)-forms living on exceptional divisors of
     the resolution NEVER contribute: the budget must come from H^{2,0} of the normalization-type part of Y~ that
     maps finitely onto Y.
 (c) Lefschetz for a singular ample divisor: H^2(X) -> H^2(Y) is an isomorphism for every effective ample Y
     (Andreotti–Frankel on X \ Y), so H^2(Y) = Q h and all (2,0)-classes of Y~ come from the non-normal / non-rational
     singularities of Y (attack1 §1.4 adjoint formula: h^{2,0}(Y~) = dim coker(S_d -> H^0(ω_Y/π_*ω_Y~))).
Under simplicity of V = H^4(X)^tau_prim (dim 1282, (208,866,208)) the pure part W_{-4}H_4(Y) must CONTAIN a copy of V;
so the singular locus of Y must "carry" V: a double surface Σ with h^1(I_Σ(d)) >= 208 whose non-d-normality is
transcendental, or a curve of transversally non-rational singularities whose resolution data is of V's type.
Tensor-product budgets (Fermat join: H^1(C_1) ⊗ H^1(C_2), the "100 forms" of attack1 §4.3) are impossible for a simple
V of rank 1282 (V would be a tensor product of curve Hodge structures).

## 4.2 Audit of universal tau-even divisors on W = X/tau -> P^2_X x P^2_Y (cubic a + b t + c t^2 + d t^3)
Notation a = F60, b = F42, c = F24, d = F06; on X the cubic has the root t = 1 and factors as
      (t - 1) (d t^2 + (c + d) t + (b + c + d))        [checked symbolically].
| divisor on X | equation (mod F) | degree | geometry | budget |
| R (ramification: line through [X],[Y] tangent at p) | E_X F = 6a+4b+2c  (= -(E_Y F) on X) | 6 | p = double root t = 1: b+2c+3d = 0; smooth for generic F; ample section of X | 0 by Lefschetz (c); E2 fails (packet R6/R11) |
| R' (residual tangency: the line is tangent at the OTHER pair ±p') | D_R' = (c+d)^2 + 4ad  (disc of the residual quadratic) | 12 | R' -> Δ ⊂ P^2xP^2 2:1 (p' <-> -p'), branched over Δ ∩ (V(a) ∪ V(d)); Δ is the (12,12) discriminant, cuspidal along the triple-root surface; R and R' are two different double covers of Δ | 0 for the invariant sector by the equivariance argument of 4.4 (R'/tau = Δ birationally, Δ unibranch + Lefschetz => H^2(Δ~) all (1,1) up to exceptional classes killed by (b)); the anti-invariant budget H^0(Δ~, Ω^2(log B) ⊗ O(-3,-3)) is irrelevant for A |
| full branch preimage X ∩ V(Δ~) = 2R + R' | Δ~ (degree 24, dossier §1.3) | 24 | reducible | sum of the two rows: 0 |
| Hessian/parabolic divisor | det Hess F | 24 | transversal A_1 along the corank-2 curve, ample | 0 by (c) (rational singularities) |
| smooth bihomogeneous sections f^*B, B ⊂ P^2xP^2 of bidegree (m,n), n even | any transverse b(X,Y) | m+n | ample on the smooth 4-fold W (W is smooth: checked at t = 0, ∞ via smoothness of a = 0, d = 0; elsewhere X/tau is smooth) | 0: Lefschetz on W, H^2(W) = Q^3 all (1,1) |
| ternary cones g(X), g(Y) (fibred over C_g by the S_x through C_Y) | g | k | transversal cone over C_g along C_Y; smooth fibres S_x (p_g = 10) with full Lefschetz monodromy | 0 (part 1 §1.2) and E1 fails for k <= 12 (part 3) |
| products of pieces (10 quadratic monomials in a,b,c,d) | e.g. a d, (c+d)^2 | 12 | reducible except the generic combination | E1 tested in part 3 |
Conclusion: every divisor definable from the cubic-cover geometry alone has budget 0.  This is the structural reason
mechanisms (i)-(iii) fail for a symmetry-free even sextic: all their supports are built from universal data and are
either smooth (Lefschetz) or singular along loci whose contributions are exceptional (b).

## 4.3 The most promising remaining candidate and its explicit equations
R' is the only universal, irreducible, NON-smooth, non-cone divisor of degree < 24 on an even sextic; it is also the
simplest divisor not of "polar" type (its degree 12 lies above the Jacobian-type range d <= 9 of the E1 classification).
Explicit data (variables X0,X1,X2,Y0,Y1,Y2; all pieces are the Y-degree parts of the record's F):
      D_R'  =  (F24 + F06)^2 + 4 F60 F06          (equivalently (F60 + F42)^2 + 4 F60 F06 modulo F),
      E1 required:  A D_R' ∈ J(F)_18  (Y-odd sector),   tested in part 3 (job construct-dense-rational-univ12),
      E2 (if E1 passes): the full fixed-D system at d = 12 (unknowns ~26k, equations ~162k): not run.
Hypotheses that the §2.7 recipe would need on top: a curve family on R'~ with polynomial normal and a selected source —
none can exist unless the invariant budget is nonzero, and 4.4 shows it is zero.  So R' fails for a structural reason
(independently of E1): it is a tau-stable divisor whose tau-quotient is birational to the discriminant Δ.

## 4.4 Equivariance reduction: the whole problem lives on the smooth fourfold W
The Gysin map H^2(Y~) -> H^4(X) is tau-equivariant for a tau-stable Y, and ω_A is tau-invariant, so only
H^2(Y~)^tau = H^2((Y/tau)~) can serve A: the needed budget is h^{2,0}((Y/tau)~) >= 208 for the divisor Z = Y/tau of W.
For non-tau-stable Y the prescribed-support theorem replaces Y by Y ∪ tau(Y) and the invariant classes come from
(1 + tau)·Im(H^2(Y~)), which still needs h^{2,0}(Y~) >= 208.  Consequences:
 * R (= smooth double cover of Δ), R' (= the other double cover of Δ), X ∩ V(Δ~): all have (Y/tau)~ birational to Δ~,
   hence invariant budget 0.  R' is therefore excluded for A WITHOUT the E1 test; the E1 run of part 3 is a consistency
   check of this prediction (if E1 passed, E2 would still have to fail by this argument).
 * Any Z ⊂ W that is ample and smooth has H^2(Z) = H^2(W) = Q^3 (Lefschetz): budget 0.  Pic(W) ⊗ Q = <f^*H_X, f^*H_Y, E_a>
   (E_a = exceptional P^2-bundle over C_a = the strict-transform side of the ternary cone a = 0), so "universal" divisors
   of W are the bihomogeneous pullbacks and the cone-type divisors, all audited above.
 * Hence a mechanism for a generic even sextic needs a NON-universal, singular divisor Z ⊂ W whose singular locus carries a
   1282-dimensional transcendental Hodge structure of type (208,866,208) in the pure part of H_4(Z); no construction of
   such a divisor from the equation of a symmetry-free F is known to this agent, and the search space (D of degree >= 10,
   non-Jacobian E1 kernel of dimension 294/261 at d = 10 growing to ~700) is the only place left to look (search lens).
