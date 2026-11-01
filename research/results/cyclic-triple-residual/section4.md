## 4. Mathematical ideas developed, partial constructions, what would finish them

**4.1 The Fermat mechanism for this class and why its deformation dies.** On F_0 = Σ Z_i^6 the class of X0^3 Y0^3
(character (4,1,1 | 4,1,1), both triples ≡ 0 mod 6) is a *join-type* class: it is served by D0 = X0^6+X1^6+X2^6 with
the strict certificate N = A·(X0,X1,X2,0,0,0)/6, Q = 0 (verified exactly, section 2). The identity behind it is
block-Euler: D0 = E_X(F_0)/6, so Σ N_i ∂_i F = A D0 (E1), and div(AG) = A, AG(D0) = A D0 (E2 = 0) because D0 is
X-homogeneous of degree 6 and A of X-degree 3. Geometrically X_0 ∩ V(D0) is the join of the plane sextics C_X, C_Y,
singular along both with transversal type the cone over a sextic curve; its resolution is a P^1-bundle over C_X × C_Y
with H^{2,0} = H^{1,0}(C_X) ⊗ H^{1,0}(C_Y) (dimension 100), and (4,1,1) is a character of H^{1,0} of each curve.
For the coupled F the same formula gives D = E_X(F) = Σ X_i ∂F/∂X_i: T = X ∩ V(D) is the **ramification divisor of
the degree-6 map X -> P^2_X × P^2_Y**, (X,Y) -> ([X],[Y]) (fibres = roots of the binary sextic F(sX, tY)). It still
contains the fixed curves Γ_X = X ∩ P^2_X = {X0^6 + X1^6 + X2^6 + tX0^3X2^3 = 0} and Γ_Y = {Y0^6+Y1^6+Y2^6 = 0}, but
only with multiplicity 2 (the couplings X0^2Y0Y1Y2^2 and X0X1Y1^3Y2 have X-degree 2; the Y-degree-2 couplings do
the same on the other side): along Γ_Y the quadratic part is t·X0·Y1Y2·(X0Y0Y2 + X1Y1^2), a rank-2 form, so the
transversal singularities are A_k (rational) and the exceptional surfaces carry no 2-forms. This is the concrete
reason the deformation of the Fermat certificate fails (exact rank 370/388, augmented 371), together with all
weighted Euler fields Σ m_i Z_i ∂_i F (= ramification divisors of the maps to weighted products): the E2 residual is
(A/6)·Σ_k k(6-k) F^{(k)} with F^{(k)} the X-degree-k part of F, i.e. (4tA/3)·(the five mixed couplings), and it is not
in F·S + B(Koszul corrections).

**4.2 Multiplicity along a plane as a linear condition inside the E1 coset.** For a plane Π ⊂ P^5 with
Γ = X ∩ Π a plane sextic (genus 10 when smooth) and a form D vanishing to order m along Π, T = X ∩ V(D) is locally
Γ × (cone over the plane curve B_P = {D_m(P;·) = 0} ⊂ P(N_{Π,P})) at smooth points of Γ (X is transverse to Π, so one
can solve one Π-coordinate). For m = 3 the transversal singularity is simple elliptic and the exceptional surface
E -> Γ is an elliptic fibration whose Leray piece H^1(Γ, R^1) ∩ F^2 has dimension >= g(Γ) - 1 = 9 (exactly
g(Γ)·g(B) = 10 if the fibration is a product, i.e. if the cubic leading form does not depend on P ∈ Γ). Since
"vanishing to order m along Π" is linear in D, it can be intersected exactly with the E1-admissible space
⟨D8⟩ + J(F)_8 (sector (0,4), 40 generators). Results (section 2.6): among the 20 coordinate planes only
Π = {X0 = X2 = Y1 = 0} admits order 3 (an 18-parameter affine family), with Γ = {X1^6 + Y0^6 + Y2^6 = 0} a Fermat
sextic; no coordinate plane admits order >= 4, and the tested members of the family are inconsistent (exact).
Non-coordinate planes (9 parameters) and the choice inside the 18-dimensional family remain open; the family is
small enough for a *parametric* treatment (the system is linear in (N, Q) with coefficients linear in the 18
parameters and right-hand side quadratic), which is the natural next computation.

**4.3 Why E1 is so restrictive here and what it says about T.** [A]·R_{18-d} = R_{24-d} for d <= 7 means every
(1,3)-class is A times a (2,2)-class (surjective Yukawa-type coupling) — the deformed Jacobian ring has no
"accidental" relations in low degree except the single sporadic one in degree 8 (and three more in degree 9).
Denominators in J(F) are tangency divisors of vector fields; denominators in the degree-8 coset D8 + J(F)_8 are the
only other option below degree 9. D8 is not congruent mod J to anything divisible by a coordinate power, to a cone
over any coordinate plane or line, or to a multiple of A, C, Z', H, H', C^2, CZ', Z'^2 (exact membership tests), so
it is not an obvious product of the structural pieces of F. Its t -> 0 limit in the family F_t (section 2.7) shows
which of the ~700 Fermat annihilators (X0^2, Y0^2)_8 survives the coupling to first order.

**4.4 What a complete argument would need.** Either (a) a D of degree >= 8 (necessarily in the E1-admissible
cosets computed here) whose divisor has non-rational singularities along curves or at many points, with
Σ h^{2,0}(T̃_i) large enough to contain the 41-dimensional (or 67-dimensional) sector — the plane/line families of
4.2 are the systematic way to generate such D compatibly with E1; or (b) a non-constructive argument as in the
companion solve (a family of Fano-fibred or rationally connected quotients whose transport hits the zeta3 sector),
for which the only structure specific to V_1 found so far is the isotrivial j = 0 elliptic family π^{-1}(L) over the
lines L ⊂ W (dossier §5.E5); or (c) a proof of non-existence of low-degree certificates: the E1 table plus a
parametric exclusion of the degree-8 coset would make "no certificate of degree <= 8" rigorous.
