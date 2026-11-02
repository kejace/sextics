# four-plus-two-residual — CONSTRUCTION lens — part 3: the core threefold T_H, contact families, and the exact failure point

Script: /tmp/hodge/work/construct42/fibration_check.py (sympy, exact), log entries step = fibration_check.

## 3.1 T_H is fibred by C_B-type genus-4 curves over P^1 x P^1 [C, exact]
For ((a:b),(c:d)) ∈ P^1 x P^1 let l = line joining (a:b:0:0) ∈ L_X = {Y0=Y1=0} and (0:0:c:d) ∈ L_Y = {X0=X1=0}. Then
   H(la a, la b, mu c, mu d) = f(a,b) la^6 + (5/7) a^3 b c d la^4 mu^2 + (10/7) a b c^3 d la^2 mu^4 + g(c,d) mu^6   [verified]
— an EVEN binary sextic in (la,mu), i.e. a binary cubic b_l(La,M) in (La,M) = (la^2,mu^2). Hence C_l := pi^{-1}(l) ⊂ T_H is
   w^3 = b_l(la^2, mu^2): a cyclic trigonal genus-4 curve with the involution s: mu -> -mu, exactly the type of C_B.
The lines meeting L_X and L_Y form a congruence of order 1, so T_H (blown up along the two genus-4 curves C_f, C_g over
L_X, L_Y) is a fibration over B = P^1 x P^1 with fibres C_l; s acts fibrewise. Consequently
   H^3(T_H) = H^2(B, R^1),  R^1 = H^1(C_l) = H^1(E_l) (s-invariant, E_l = C_l/s a j=0 elliptic curve: isotrivial family)
                                            (+) H^1(P_l)   (s-anti-invariant: the Prym of C_l, rank 6, Q(zeta3)-signature (2,1)),
   V_zeta = H^3(T_H)_zeta^{s=-1} = H^2(B, R^1_{Prym, zeta}): "H^2 of a rational surface with coefficients in a Picard-modular VHS".
The target's core factor is a class of this type; the binary factor is the SAME kind of object (H^1 of one such curve).

## 3.2 The Prym family is a genuine 2-parameter family [C]
Torus-invariants of the binary cubic, I1 = b1^3/(b0^2 b3), I2 = b2^3/(b0 b3^2), as functions of (a/b, c/d): Jacobian rank 2
at five random rational points (exact), so the period map B -> (Picard modular surface) is dominant: no isotriviality,
no fixed part over B, and generic expectation that V is Hodge-simple with End = Q(zeta3) [H]. (In the split limit K = 0
the cubic b_l has b1 = b2 = 0: the Pryms are constant up to the torus, V_0 splits into CM pieces — attack1 §1.3.)
Consequence for sweeps: a surface swept by fibres C_l over a curve Gamma ⊂ B has H^1(Sigma~)_zeta = fixed part of the Prym
family over Gamma, which is 0 for general Gamma and at most a Q(zeta3)-rank-3 piece for special Gamma: never V. [P]

## 3.3 Contact surfaces and the Abel–Jacobi triviality theorem [P]
For a quadric c ∈ H^0(P^3,O(2)) put V_c := {H = c^3} ⊂ P^3 and Sigma_c := {u = c(x)} ∩ T_H ≅ V_c (the lift; pi^{-1}(V_c) =
Sigma_c ∪ delta Sigma_c ∪ delta^2 Sigma_c). This is the 4-variable form of the book's perfect-power contact (§2.3: Q = cH^m + DK with
m = 3, D = H - c^3, K = 1) and of the Nambu/contact denominators of §2.3–2.5. Sigma_c ≅ V_c is a sextic surface in P^3:
regular (q = 0) whenever normal with rational singularities; so (part 1, consequence 3) it cannot support [a] by itself.
The natural way to manufacture irregularity is a 1-parameter family of curves on T_H and its cylinder map. Take
   E_gamma := Sigma_{c_gamma} ∩ pi^{-1}(V'_gamma), gamma ∈ Gamma  (V'_gamma ⊂ P^3 any surfaces, e.g. planes of a pencil),
or more generally any family of curves E_gamma ⊂ Sigma_{c_gamma} whose class on Sigma_{c_gamma} is a multiple of the hyperplane
class. THEOREM: the cylinder map H^1(Gamma) -> H^3(T_H) of every such family is zero. Proof: AJ(E_gamma) depends only on
c_gamma (two curves in |O_{V_c}(m)| are linearly equivalent Cartier divisors on V_c, hence rationally equivalent 1-cycles on
T_H after the lift), so AJ: Gamma -> J^2(T_H) factors through c: Gamma -> P^9 and a rational map from P^9 to an abelian
variety is constant; a family with constant Abel–Jacobi map has zero cylinder map. For a FIXED plane Pi the curves
{u = c} ∩ pi^{-1}(Pi) are the linear system |O(2)| on the triple plane pi^{-1}(Pi) (q = 0), giving the same conclusion directly.
Remark: every plane contact curve has an ambient cube root (h^0(O_F(2)) = 6 for every plane curve F), so ALL families of
plane contact curves are covered; so are all complete-intersection contact curves on quadrics (2-normal).

## 3.4 The explicit sweep family and its precise failure [P + C-ready]
Pencil l_1 - gamma l_2 of planes; binary forms lam_1, lam_2 in (l_1,l_2) of degrees r, r+3e; a conic varying with gamma,
c_hat = sum_i l_1^i l_2^{e-i} c_i(x). Denominator
   D = lam_2(l_1,l_2) H(x) - lam_1(l_1,l_2) c_hat(x)^3,   deg D = 6 + r + 3e   (core-only),
V(D) ∩ Pi_gamma = the plane sextic {H = lam(gamma) c_hat_gamma^3}, lam = lam_1/lam_2, so H is a cube on every member; the lift
needs lam^{1/3}: parameter curve Gamma: mu^3 = lam(gamma), cyclic over P^1, g(Gamma) = 2r + 3e - 2, H^1(Gamma)_zeta of dim g.
Smallest members: r = 1, e = 0: D = L_2 H - L_1 c^3 (degree 7, Gamma rational); r = 2, e = 0: D = M M' H - L L' c^3 with
L,L',M,M' in one pencil (degree 8, g = 2); r = 3: degree 9, g = 4. E1 for this family reads a·lam_1(l) c_hat^3 ∈ J(H) (the lam_2 H
term is Jacobian): at degree 7 this is 20 equations in 14 unknowns (generically violated — consistent with g = 0); at
degree 8 it is 10 cubic equations on (L,L',c) ∈ P^3 x P^3 x P^9 (expected 5-dim solution set). By the theorem of 3.3 every
member of this family, at every degree, has zero cylinder map, so no member can support [a] regardless of E1/E2 — and the
same holds for the Nambu-type realisations W = ∇D x ∇C x ∇E (4-variable Nambu bracket with an ambient cube root C).
Exact failure point of the book's §2.3/§3A contact recipe in four core variables: the source surface is a regular sextic
surface and its curve families are hyperplane-section families on it; the residue it would transport is the Gysin image of
H^1 of a regular surface, i.e. zero in the (2,1)-sector. (The 3+3 case works because the analogous object is a plane
sextic CURVE, whose H^1 is nonzero.)

## 3.5 Where a construction could still come from (not achieved)
 (a) Curves on contact surfaces Sigma_c that are NOT multiples of the hyperplane class: Noether–Lefschetz loci
     NL ⊂ P^9 of the family {V_c}. For curves with an ambient cube root and h^0(N) = 4d (contact lines, conics, (2,3)-CI
     curves: finitely many of each kind) the NL locus is a linear space (c|_E = fixed root), hence rational: AJ-constant.
     Positive-dimensional families: (2,b)-complete intersections with b >= 4 (expected dim 2 for b = 4: unknowns
     (Q, K_4, c) with H - c^3 ∈ (Q, K_4)_6: 40 equations in 42 projective parameters). Whether the AJ map is non-constant
     on an irrational curve inside such a locus is open; under simplicity of V the parameter curve must have
     zeta-genus >= 48, which forces denominators of degree far beyond 30.
 (b) Contact space curves whose cube root is NOT the restriction of a quadric (non-2-normal curves): not on any Sigma_c;
     the AJ argument does not apply. No construction attempted.
 (c) Deciding whether V is Hodge-simple would settle how large any support must be; point counts of T_H cannot reach the
     degree-48 Frobenius polynomial of V, and no cheaper invariant was found.
 (d) The secondary (normal-function) route of part 1(b): needs >= 32 monodromy-fixed eta^2-holomorphic 2-forms on a
     pencil of surfaces in Y_H; no candidate.

## 3.6 Lean / seed status
No uniruled auxiliary variety with the target's character was found, so there is no polynomial curve family, no incidence
Jacobian and no finite-field seed to report. The reduction of part 1 is itself formalisable only in the direction already
used by attack1 (reduced 5-variable certificate => 6-variable certificate, a polynomial identity); the forward direction
is Hodge-theoretic and outside the shared library.
