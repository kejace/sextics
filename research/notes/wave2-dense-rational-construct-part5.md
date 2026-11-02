# wave2 dense-rational CONSTRUCT — part 5: draft paragraphs for the discussion post (merge into the lane draft)

(Status flags in [brackets] must be resolved from part 3 before posting.  Nothing has been posted.)

## Coordinate-free exclusion of covers with an exchange, and of (quasi-)reflections
Every construction in the book that reaches a tau-invariant class on an even sextic without a diagonal symmetry adjoins a
root T^6 = g(u) and uses an exchange y <-> T (holdout family, §2.7-2.8).  Writing the pulled-back equation in the base
ring C[alpha]/(g - 1) and comparing the coefficients of y^b T^m with those of y^m T^b, the exchange forces, in the chosen
linear coordinates, F = Phi(rest, y u_1, ..., y u_r) + c y^6 + g(u) with c != 0 (part 2 §2.1).  At the point e_y this
gives F(e_y) = c and a rank-one Hessian.  For the dense sextic the locus {rank Hess F <= 1} ⊂ P^5 (120 octic minors) and
even {rank Hess F <= 2} (210 minors of degree 12) are empty: msolve (p = 1073741789, DRL) returns leading ideals
containing pure powers of all six variables, which certifies emptiness over C by the Macaulay-rank argument.  The same
generator applied to the Fermat sextic correctly returns the six coordinate points (binomials e_i^4 e_j^4).  Consequently no
exchange cover of §2.7 type exists for this F in ANY linear coordinates, for any block split and any base sextic, and no
quasi-reflection of order 3, 5 or 6 preserves F.  [V3: F is even in no linear form at all — pending/decided.]

## All ternary cones fail E1
A g ∈ J(F)_{k+6} has no nonzero solution g ∈ C[X]_k or C[Y]_k for k = 7, 8 [9..12] (ranks attain the parity-split Hilbert
caps, hence rigorous over Q).  Geometrically E_g = X ∩ {g(X) = 0} is swept by the sextic surfaces S_x = X ∩ span(x, P^2_Y)
through C_Y; by Lefschetz theory of this net the monodromy invariants of H^2(S_x) are (1,1), so the resolution of E_g has
no holomorphic 2-forms: cones serve no primitive class at any degree.

## Why no universal divisor can work (equivariance + Lefschetz)
For a tau-stable support Y the Gysin map is tau-equivariant, so the invariant class ω_A is served only by H^2((Y/tau)~):
the problem lives on the smooth fourfold W = X/tau (the cubic cover a + bt + ct^2 + dt^3 = 0 of P^2 x P^2, smooth
because {F60 = 0} and {F06 = 0} are smooth plane sextics).  H^2(W) = Q^3 is algebraic, so an ample smooth divisor of W has
no (2,0)-forms (Lefschetz), and the two universal singular candidates, the ramification divisor R = {E_X F = 0} and the
residual-tangency divisor R' = {(F24 + F06)^2 + 4 F60 F06 = 0} (discriminant of the residual quadratic after removing
the root t = 1), are double covers of the (12,12) discriminant Δ, whose resolution has H^2 algebraic (Δ is unibranch along
its cusp surface; classes on exceptional divisors over surfaces or curves push forward to zero in H^4).  The full branch
preimage X ∩ V(Δ~) = 2R + R', det Hess F, bihomogeneous sections and products of pieces are covered by the same two facts.
A mechanism for a symmetry-free even sextic therefore needs a non-universal singular divisor of degree >= 10 whose
singular locus carries the 1282-dimensional sector; none is known.  [E1 for D_R' and the ten products: pending/decided.]
