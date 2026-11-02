# four-plus-two-residual — CONSTRUCTION lens — part 1: framework and the core-reduction theorem

Status labels: [P] proved here (sketch given), [C] computed (script in /tmp/hodge/work/construct42/), [H] heuristic, [X] mechanism shown to fail.

## 1.1 Where the target lives (recap, all from dossier/attack1/probe, re-checked)
F = H(x) + B(y), x = (X0,X1,Y0,Y1), y = (X2,Y2); A = X0^3 Y0^3 =: a (core-only).
G = <h, r2, s> order 24; chi_A: h -> zeta3, r2, r5, s -> -1; stabiliser K0 = {1, tau, r2r5, r2s}, no junior element.
Sector of A: Hodge numbers (32,80,32), rank 144, = [H^3(T_H)_zeta^{s=-1} (x) H^1(C_B)_{zeta^2}^{sigma=-1}]
  core factor V_zeta := H^3(T_H)_zeta^{s=-1}: types (0,32,16,0); rational HS V (zeta+zeta^2 parts): dim 96, level 1,
  Q(zeta3)-type (32,16): "H^1 of a 48-dim abelian variety of Weil type".
  binary factor: the line <dt/w> in H^{1,0}(C_B)_{zeta^2}, inside the absolutely simple 3-dim Prym P (probe).
M := the sector = V (x)_{Q(zeta3)} H^1(P)' (a twist); if V and H^1(P) are "generic" M is Hodge-simple [H].

## 1.2 Pencil structure and the core-reduction theorem [P, one caveat]
Project from the plane {X2 = Y2 = 0}: Xt = Bl_{S_H} X -> P^1 = P(X2:Y2), fibre over t: Y_t = {H(x) + B(t,1) Y2^6 = 0} ⊂ P^4,
isomorphic to Y_H := {U^6 + H = 0} via U = B(t,1)^{1/6} Y2 (attack1 §1.5). Hence on U = P^1 \ {6 roots of B}:
   R^3 pi_* Q = (+)_{k=1..5} V_k (x) L^k,   V_k := H^3(Y_H)_{eta^k} (dim 104 each), L := rank-1 local system of B^{1/6}
(monodromy eta around each root; L^k pairwise non-isomorphic, nontrivial). The element h acts on Xt fibrewise as U -> eta U.
Target: [A] = [aU] (x) omega with [aU] in V_2^{(2,1)} and omega = dt/w in H^1(P^1, j_* L^2) ⊂ H^1(U, L^2): a PURE tensor in the
Leray piece Gr^1_L H^4 = H^1(U, R^3); its Gr^0 (restriction to a fibre) is 0. The joining classes (cones over S_H with vertex
at the 6 roots of B) are vertical and h-invariant, so they are not in the target sector.

THEOREM (core reduction). Let Z ⊂ Xt be a divisor supporting [A]. WLOG Z is h-invariant (replace Z by its h-orbit), i.e.
D is bihomogeneous of bidegree (d-e, e) in (x ; X2,Y2), and the Gysin preimage beta can be taken in H^2(Zt)_{eta^2}.
Then with Z_t := Z ∩ Y_t (a mu6-invariant surface of Y_t ≅ Y_H, = pi_Y^{-1}(V_t), V_t = {D(x,t,1) = 0} ⊂ P^3):
 (a) [primary] if beta ∈ L^1 H^2(Z°) (Leray over U), then for general t the fibrewise Gysin image
     W_Q := im(H^1(Zt_t)(-1) -> H^3(Y_t)) (a rational, mu6-stable sub-HS) contains [aU]; hence the single surface V_{t0} ⊂ P^3
     (general t0) supports [a] on T_H, and by attack1 §4.2 (converse direction) + prescribed-support theorem there is a
     certificate with the CORE-ONLY denominator D'(x) = D(x, t0, 1) of degree d - e <= d.
 (b) [secondary] otherwise beta has a nonzero component in H^0(U, R^2_Z)_{eta^2} (a monodromy-fixed eta^2-class on the fibre
     surfaces, necessarily killed by the fibrewise Gysin since H^4(Y_t) = Q(-2) has trivial character) and [A] arises from
     the secondary map H^0(U, R^2_Z) -> H^1(U, R^3) given by the extension class in Ext^1_U(R^2_Z, R^3_X)
     (a normal-function-type mechanism). For core-only D this cannot happen (R^2_Z|_{eta^2} = H^2(Zt_H)_{eta^2} (x) L^2 has
     no invariants). If M is simple, (b) needs >= 32 monodromy-fixed holomorphic 2-forms of character eta^2 on the fibres.
Proof of (a): the Gysin map is a morphism of Leray spectral sequences over U (both maps are over P^1); on an affine curve
base E_2 = E_oo. The image local system W ⊂ (+) V_k (x) L^k is (+) W_k (x) L^k by Schur (L^k irreducible, pairwise
non-isomorphic), so im(H^1(U,R^1_Z)) ⊂ (+) W_k (x) H^1(U, L^k) and a pure tensor [aU] (x) omega lies there iff [aU] ∈ W_2.
W_Q is the fibre of a sub-VHS, hence a Q-sub-HS of H^3(Y_H); [aU] ∈ W_Q (x) C. Descend to T_H = Y_H/mu2 (the class has
character eta^2, invariant under eta^3). The converse (core-only support on T_H => 6-variable certificate) is attack1 §4.2.

CONSEQUENCES.
 1. The binary block is cohomologically AUTOMATIC: the Prym factor dt/w is supplied by the base of the pencil
    (H^1(U, L^2)). Nothing about C_B, its Prym P, its simplicity, W_B, Kuga–Satake on W_B, or carriers for C_B can decide
    existence; benediamond's question ("does the Prym admit a correspondence supplying support?") is the wrong side.
 2. Existence (primary route) <=> the (2,1)-class [a] ∈ H^{2,1}(T_H)_zeta^{s=-1} is supported on a surface of the
    Calabi–Yau threefold T_H: u^3 = H ⊂ P(1,1,1,1,2); equivalently a core-only D works. This is the generalized Hodge
    conjecture (coniveau 1) for the level-1 piece V of H^3(T_H).
 3. If V is Hodge-simple [H], every support surface Sigma ⊂ T_H has dim_C H^1(Sigma~)_zeta >= 48, i.e. q(Sigma~) >= 48
    (Gysin images are sub-HS; a sub-HS of a simple HS containing [a] != 0 is everything). A surface in T_H with
    irregularity >= 48 is the object any construction must produce. Low-degree core-only denominators are therefore
    hopeless unless V splits; this is consistent with every exclusion so far (deg D <= 5 excluded exactly; all core
    kernel representatives up to degree 12 inconsistent in the reduced model).
 4. Any certificate of degree d yields a core-only certificate of degree <= d (primary case). Searches may be restricted
    to core-only D in the 5-variable reduced model (attack1's threefold.py) without loss, modulo the secondary case.

## 1.3 What a working mechanism must produce (translation of the book's recipes)
The accepted nonconstructive proofs = (uniruled auxiliary variety W) + (polynomial curve family with invertible incidence
Jacobian) + (integration/trace) + (correspondence to X). After 1.2 the auxiliary variety must dominate a surface
Sigma ⊂ T_H with H^1(Sigma~)_zeta ->> V_zeta, or a threefold Z ⊂ X with H^2(Z~) ⊇ M(1). The only general supply of
irregular surfaces on a threefold is a 1-parameter family of curves {E_gamma}_{gamma ∈ Gamma} (Sigma = ∪ E_gamma,
H^1(Sigma~) ⊇ H^1(Gamma)) with a NON-CONSTANT Abel–Jacobi map Gamma -> J^2(T_H) (equivalently nonzero cylinder map
H^1(Gamma) -> H^3(T_H)); if V is simple the cylinder image must be all of V, so g(Gamma)_zeta-part >= 48.
Part 3 shows that the natural families (lifted plane sections of the contact surfaces {H = c^3}) are Abel–Jacobi
trivial for ALL degrees, which is the precise 4-variable reason the 3A/2.3 "perfect-power contact" recipe fails.
