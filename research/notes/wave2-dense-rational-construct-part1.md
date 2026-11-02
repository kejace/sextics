# wave2 dense-rational CONSTRUCT — part 1: landscape and mechanism triage

(written 2026-11-02 ~02:00 local; restart after the killed attempt, which left no parts)

## 1.1 Target and the geometric statement to realize
* X = {F = 0} ⊂ P^5, F the 236-term even sextic (Fermat + (1..17)/9973 * all 230 non-pure tau-even monomials).
  K_X = O (CY fourfold), tau = (X, -Y), tau*Omega = -Omega, so H^{4,0}(X)^tau = 0.
* Target: q = 1, A = X0^3 Y0^3, residue class omega_A = A Omega/F^2 ∈ H^{3,1}(X)^tau_prim (dim 208; sector (0,208,866,208,0)).
* A certificate exists iff omega_A is supported on a divisor E = V(D) ∩ X (Ch. 1 prescribed-support theorem),
  i.e. iff omega_A ∈ Gysin image of H^2(E~) for a resolution E~ -> E. Its (3,1) part comes from H^{2,0}(E~).
* Generic expectation (dossier 5.9, Ch. 4.4): H^4(X)^tau_prim is an irreducible Hodge structure for Mumford–Tate
  generic even F, so ANY divisor supporting omega_A supports the whole sector: need h^{2,0}(E~) >= 208 with
  surjective Gysin map onto H^{3,1}^+. So a "mechanism" here = a divisor construction valid on an open set of
  even sextics (what Diamond's paper says is out of reach: ~17 of 235 dims).

## 1.2 Rigid structure of X/tau that any mechanism must use (all even sextics)
* X/tau is birational to W = {a t0^3 + b t0^2 t1 + c t0 t1^2 + d t1^3 = 0} in the P^1-bundle P(O ⊕ O(2,-2))
  over P^2_X × P^2_Y, with (a,b,c,d) = (F60, F42, F24, F06): the "natural cubic cover" of Ch. 4.7.
  W contracts to X/tau; the two sections t = 0, t = ∞ become the fixed curves C_X, C_Y (age 3/2 points).
* Fibres of W -> P^2_X: for x ∈ P^2_X, the surface S_x/iota where S_x = X ∩ P^3_x, P^3_x = span(x, P^2_Y),
  S_x: lambda^6 a(x) + lambda^4 b(x,Y) + lambda^2 c(x,Y) + d(Y) = 0 (sextic surface in P^3, even in lambda,
  p_g = 10), S_x/iota = sextic in P(1,1,1,2) with K = O(1), p_g = 3, q = 0 (Ch. 4.7 gives p_g = 20 for the
  (6,3)-surface over the OTHER passive parameters). Every S_x contains C_Y = P^2_Y ∩ X.
* Consequence for cone denominators D = g(X): E_g = X ∩ {g = 0} is the union of S_x over x ∈ C_g = {g = 0} ⊂ P^2_X,
  all glued along C_Y; blow-up of C_Y = incidence {(x,p): x ∈ C_g, p ∈ S_x} -> C_g. By the decomposition theorem
  H^{2,0}(E~_g) = (2,0)-part of the monodromy invariants of H^2(S_x) over C_g (fibres regular, R^1 = 0).
  The S_x are the codim-2 linear sections of X through the fixed plane P^2_Y; Lefschetz theory on the blow-up
  of X along C_Y gives irreducible monodromy on the vanishing part, invariants = image of H^2(X) = Q h, so
  h^{2,0}(E~_g) = 0 for every ternary g with C_g ⊂ P^2_X ample (any degree >= 1): ternary cones (X-only or
  Y-only) serve NO primitive class. This closes the gap the referee flagged in attack1 §1.4 (the
  "rational singularities elsewhere" hypothesis is replaced by the decomposition-theorem/fixed-part argument;
  singular fibres only add classes of level <= their own, which are (1,1) on surfaces with p_g from
  the smooth fibre... to be stated carefully in the final note). [theory, no computation needed]

## 1.3 Triage of the five requested mechanism classes for THIS lane
(i) Exchange on a cover (Ch. 2.7). Holdout structure: in coordinates (y; u, v) on Y-space, T^k = g(u,v),
    F must be symmetric under y <-> T in C[X][alpha,beta]/(g-1) with alpha = u/T, beta = v/T. Since
    (g - 1) has no elements of degree <= 5, the symmetry forces linear conditions on F (per X-monomial):
    all Y-monomials with (deg_{u,v}, deg_y) = (m, b), m != b, must vanish unless m + b = 6 and (m,b) ∈ {(6,0),(0,6)}.
    This is a closed condition on the pair (linear form y, complement) that is testable by linear algebra
    for ALL choices of y: see part 2 (rigorous negative expected: the needed vanishing of e.g. the
    y^1 * quintic(u,v) block is 6 conditions per X-monomial... but y is a free linear form; we test the
    exact invariant condition instead of sampling).
(ii) Fano/surface-fibration companions of quotients by subgroups of Aut: Aut_proj(F) = {1, tau} (packet R10,
    referee claim 5-6 for diagonal/permutation/infinitesimal; the full finite linear group is not classified
    but a non-diagonal automorphism of a 236-term generic form would be a miracle). X/tau: not uniruled,
    kappa >= 0 (2-canonical form Omega^{⊗2} descends; age 3/2 ⇒ canonical sing.). Twisted companions
    (X × Y)/<tau × iota>: ages >= 3/2 + 0 ≥ 1 at all fixed points (permutation-theory note O1) and
    H^{4,0}(X)^- ⊗ H^{top,0}(Y)^- gives a holomorphic top form whenever iota acts by -1 on H^{top,0}(Y)
    (forced for an involution with fixed points on a curve); so no twisted companion is uniruled.
    Surface fibrations of X/tau: over P^2_X fibres have p_g = 3 (not 0); over the Ch. 4.7 passive
    parameters p_g = 20. Fails (quantified in Ch. 4.7).
(iii) Galois closure of W -> P^2 × P^2 (S_3): V_std^K = 2, 1, 0, 0 for K = 1, C_2, A_3, S_3; the quotients
    with V_std^K != 0 (W itself and the A_3-closure's C_2-quotient = W) are not uniruled. Other natural finite
    maps: X -> P^5/tau (Veronese cone), S_x -> P^2_Y (triple covers): same cubic. No new intermediate
    quotient. Fails.
(iv) q-raising / Gauss–Manin: along F_t = F0 + t P (Fermat pencil) the Fermat certificate for A has
    first-order obstruction (deform.log: 4/12 seeds inconsistent, modular). A sector theorem on a
    sub-family (e.g. coordinatewise-even sextics, holdout chain, f(X) - f(Y)) does not propagate off the
    sub-family: certificates are not deformation-invariant (the Hodge class moves with F; the divisor
    must change). No mechanism.
(v) Solved companion on the same variety: none (sameVarietyInstances = []).

## 1.4 Decision
The only mechanism class with an honest chance for a *generic* even sextic is a universal divisor
construction on W (1.2). Part 2 tests (i) rigorously for all linear forms (cheap linear algebra), part 3
runs the E1 test of the natural universal denominators that were never tested (Delta~ of degree 24,
det Hess F, F60·F06, block products, the Y-only/X-only cones at d = 7..9 via J ∩ C[X]), and part 4
records the universal-divisor idea with its Hodge-number bookkeeping and where it fails.
