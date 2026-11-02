# four-plus-two-residual — CONSTRUCTION lens — part 2: the five mechanism classes, screened against the target sector

All statements use the sector bookkeeping of part 1 (h -> zeta3, r2,r5,s -> -1). "Carrier" = auxiliary variety Y with a
G-action used in a twisted quotient (X' x Y)/K (theory-permutation note §1.1); "contraction" = translated-diagonal
correspondence of Ch. 3A.3 that removes the carrier factor.

## (i) Covers creating a reflection or exchange (book §2.7 / §2.8) — [X] no usable exchange exists
§2.7 needs a variable y that occurs in F only as y^6 and linearly (y·m), so that after T^6 = (other pure powers) the
equation is a symmetric function of (y,T). In F every variable occurs with degree 5 in a mixed term
(X0^5X1, X0X1^5, Y0^5Y1, Y0Y1^5) or with degrees 4,2 (X2^4Y2^2, X2^2Y2^4): no (6,1)-pattern, for coordinates and —
since the pattern is a property of the monomial support, not of the basis of the binary block — for any linear form in
(X2,Y2) (B = b(X2^2,Y2^2) is not palindromic, so no linear change makes B a sum of two pure powers plus a symmetric
product). The only covers with an exchange-type symmetry are the trivial ones T^6 = B(X2,Y2) (giving H + T^6: a sextic
threefold, no second variable to exchange with T) and T^6 = f(X0,X1) (K couples X0,X1 individually). §2.8 (order-6
diagonal with weights (5,1,5) on (Y0,Y1,Y2)) needs F invariant under such g: here the diagonal group is G (order 24) and
no element of G has age < 1 while fixing chi_A (dossier §1.5). After the core reduction (part 1) an exchange would have to
live on T_H: H has diagonal group {1,s} only and no (6,1)-variable. [X]

## (ii) Fano / surface-fibration companions of quotients by subgroups of G — [X] carrier dichotomy
Honest quotients X/K: the target survives iff K ⊂ K0 = {1,tau,r2r5,r2s}; all K0-elements have age >= 1, the holomorphic
4-form of X is tau-anti-invariant but Omega^2 is K0-invariant, so X/K0 (and X/<tau>, X/<r2s>, X/<h^3>) have kappa = 0:
not uniruled, no surface fibration with p_g = q = 0 (the fibres of the natural maps are CY threefolds T_H-type or curves).
Twisted quotients with a carrier Y (the companion mechanism; W_H = (T_H x E)/mu3 and W_B = (C_B x E)/mu3 are the cases
Y = E). Two requirements fight each other:
 (A) support: (T_H x Y)/mu3 and (C_B x Y)/mu3 must be uniruled (middle cohomology supported); for a curve carrier the
     (dim,0)-form test is Omega_T (x) H^{1,0}(Y)_{zeta^2} = 0 and dt/w (x) H^{1,0}(Y)_zeta = 0, H^{1,0}(C_B)_zeta (x) H^{1,0}(Y)_{zeta^2} = 0:
     i.e. H^{1,0}(Y)_zeta = H^{1,0}(Y)_{zeta^2} = 0 on the binary side — then the binary companion (C_B x Y)/mu3 carries NO
     H^1(C_B)-classes at all. (For Y = E the core side passes because H^{1,0}(E) is pure zeta, but the binary side
     (C_B x E)/mu3 is the K3 W_B with 2-form dt/w ^ omega_E: the transcendental piece H^1(P) (x)_K H^1(E) is not N^1.)
     [P: any curve carrier fails on the binary side; this is the book's "direct boundary to the elliptic trick" in general form.]
 (B) dimension count of the contraction: with Gamma ⊂ W_1 x W_2 x X of dimension 3 + 1 + dim Y and classes supported on
     D_1 x D_2 (codim 2), the image in X has dimension 2 + dim Y: a divisor only for dim Y = 1. Surface carriers (e.g. the
     cyclic cubic surface S_c: v^3 = c_3(s,t,r), for which BOTH companions {H + c_3 = 0} ⊂ P(1^4,2^3) (-K = O(4)) and
     {B + c_3 = 0} ⊂ P(1^2,2^3) (-K = O(2)) are Fano and carry H^3(T_H)_zeta (x) H^2(S_c)_{zeta^2}, H^1(C_B)_{zeta^2} (x) H^2(S_c)_zeta)
     pass (A) but fail (B): the contracted class is supported on a 4-fold, i.e. nothing; and N^2-support on the core
     companion would already require N^1 H^3(T_H) ⊇ V (the core problem itself). [X]
 Conclusion: no twisted-quotient/contraction companion reaches the mu3-sector of a 4+2 split; the obstruction is on the
 binary side for curve carriers and on the dimension count for higher carriers. After part 1 this is also irrelevant:
 the binary side is automatic and the companion idea cannot manufacture coniveau on T_H.

## (iii) Galois closures of the natural finite maps + averaging-ideal test — [X]
Natural finite maps: X -> P^3 x P^1 (degree 6, Galois = <h>, intermediate quotients X/<h^2> (uniruled, serves chi(h^2)=1
sector: the complement of ours), X/<h^3> (kappa = 0, birational to (T_H x C_B)/mu3)); T_H -> P^3 (mu3); C_B -> P^1 (mu3),
C_B -> E' = C_B/sigma (mu2); S_H -> P^1 x P^1 of degree 6 (roots of the cubic b_l in (la^2:mu^2), Galois closure of type
S_3 x Z/2 — a cover of the branch surface, not of X). The target has character chi_A != 1 on every subgroup of G outside
K0, and K0-quotients are not uniruled; the averaging ideal of any family of such quotients misses chi_A (dossier §1.5,
re-verified by the verify lens, stabiliser = K0). Non-diagonal symmetries: none (coefficient patterns 5/7,10/7 vs 5/7,15/7
break every exchange; Hessian symmetriser space 2-dim). [X]

## (iv) q-raising / Gauss–Manin / deformation sub-families — [X] for this lane
§4.5 raises q from a generic q=1 certificate (direction q -> q+1, useless for q=1) and §4.4 varies the target, not F.
Deformation to the split limit H_s = f + g + sK: at s = 0 the 2a theorem serves [a] inside the 16-dim CM piece
H^1(S_f)_{eta} (x) H^1(S_g)_{eta} (attack1 §1.3); for s != 0 the three 16-dim pieces merge into the 48-dim V (period map of
the Prym family is dominant, part 3 §3.2), so the s=0 sub-Hodge structure does not deform and neither does its support
(whose denominators have degree 6 x sheet count anyway). Deforming B is irrelevant after part 1. [X]

## (v) Solved companion on the same variety — not applicable (provenance: sameVarietyInstances = []).

## Net result of the screen
Every published mechanism and its natural twists is excluded for a structural reason that is now identified precisely:
the problem is the generalized Hodge conjecture for the level-1 piece V ⊂ H^3(T_H) of a Calabi–Yau threefold (part 1),
and the book's recipes produce either uniruled quotients of X (wrong character), elliptic companions (K3 obstruction on
the binary side, irrelevant anyway) or perfect-power contacts (Abel–Jacobi trivial in four core variables, part 3).
