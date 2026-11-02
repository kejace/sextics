# Attack 1: quadric-contact-residual  (q = 1, A = X0^3 Y0^3)

Status: IN PROGRESS (resumed after a crash of the first run; sections are filled as experiments finish).
Work dir: /tmp/hodge/results/quadric-contact-residual/  (log.jsonl = every tested denominator)

## 0. Attack plan (written before executing)

F = S0^3 + 2 P0^2 S0 + 5 S0 Q1 + (5/7) P0 Q2 + X1^6+Y1^6+X2^6+Y2^6,
  S0 = X0^2+Y0^2, P0 = X0Y0, Q1 = P1^2 + P1 P2 + P2^2 (P_i = X_iY_i), Q2 = P1 (X2^2 + 2 Y2^2).
A = P0^3.  Visible symmetry group Gamma (order 96) = G_diag(24) x| <s0 = (X0<->Y0), r1 = (X1<->Y1)>;
the residue A*Omega/F^2 is invariant under G_aff (order 48, incl. tau_Y) and ODD under s0 and r1.

Known (fqlx, #125): deg D <= 6 impossible; (J:A)_d = J_d for d <= 7; at d = 8 a single new E1 direction
D8* = 4 G_S^2 - G_P^2 (invariant sector), which fails E2 alone; no certificate keeps (Y1/X1, Y2/X2) passive;
all quotients X/H (H in Gamma^+) have kappa = 0 and are not uniruled.

Reductions used here:
 (R1) Averaging: if (D,N,Q) certifies A then so does its image under any g in Gamma (up to the sign of A under
      tau_Y); HasCertificate.add then gives a Gamma-invariant denominator prod_g g(D).  So WLOG D is
      Gamma-invariant (at the cost of multiplying deg D by the orbit size).  All systematic closure runs below
      are over Gamma-invariant D; non-invariant D are only probed by the structural families.
 (R2) E1 kernel per V-sector: with a = X0-Y0, b = X0+Y0, u = a^2, v = b^2, F = G(u,v,X1,X2,Y1,Y2), the four
      V = <s0,sigma0>-sectors a^i b^j C[u,v,...] have Jacobian parts J_00 = (uG_u, vG_v, G_xi),
      J_10 = (G_u, vG_v, G_xi), J_01 = (uG_u, G_v, G_xi), J_11 = (G_u, G_v, G_xi); E1 for D = a^i b^j D' is
      D' in (J_ij : (v-u)^3).  Computed exactly with Singular (weighted degrees).
 (R3) E2 closure modulo J_{2d} (fqlx's lemma): for D = sum c_k D_k in the E1 kernel with lifts N_k,
      B_D(N) = sum c_k c_l B_{D_k}(N_l) must vanish in R_{2d}; Koszul and F*Q terms drop.  Quadratic system in c.

Experiments planned: (0) sanity/symmetry + D8* and its factors;  (1) E1 kernel profile, all degrees, 4 sectors;
(2) Gamma-invariant E1 kernel bases at d = 8..12 and full fixed-D solves of every basis element;
(3) quadratic E2 closure of the Gamma-invariant family at d = 8, 9, 10 (11, 12 if cheap);
(4) structural families (block pieces, derivatives, Hessian minors, cones over the (X0,Y0) P^3, products of
    invariant quadrics) filtered by E1 then solved;  (5) deformation from two-reflection-cubic-surface-d2.

## 1. Results log (chronological; final tables in later sections)

- Engine note: a sign bug in the F*N(D) term of the shared engine was fixed 2026-11-01 17:25; every number below
  was computed after that (sector solver in common.py has its own implementation, validated on collision-h3
  (d=2, consistent) and by reproducing fqlx's degree-8 E1 count exactly).
- E1 per-character profile (mod p = 2^31-1): d=7: total non-Jacobian E1 directions = 0 (dim J_7 = 126)  [matches fqlx].
- Gamma-invariant E1 kernel, d=8: dim S_8^Gamma = 57; kernel = 20 = 19 (J_8^Gamma) + 1 (D8*)  [matches fqlx].
  Quadratic E2 closure mod J_16 in sector chi_A: dim R_16[chi_A] = 52, 52 quadrics in 20 unknowns, rank 35 among
  the 210 pair-monomials, no unknown killed by the square-column test alone -> Groebner (Singular) pending.
- Odd degrees: -1 in G_aff, so Gamma-invariant polynomials exist only in even degree; invariant closure is for d even.
- V-sector 00 (D in C[u,v,X1,X2,Y1,Y2], weighted degrees, Singular p=32003): dim (J:A)_d/J_d = 0 (d<=7), 1 (d=8),
  4 (d=9: these are D8* * {X1,X2,Y1,Y2}), 176 (d=10), 324 (11), 427 (12), 472 (13), 464 (14), 412, 335, 248, 169, 104, 58, 28, 12.
  Singular reports 166 new minimal generators of the colon at weighted degree 10 (E1 stops being a filter at d >= 10).
- E1 per-character profile d=8: total non-Jacobian = 1 (sector (0,0,0,0), dim S_8[0]=87, J=32, kernel 33) [matches fqlx].
- E1 per-character profile d=9: total non-Jacobian = 4, one in each of the sectors chi(X1), chi(Y1), chi(X2), chi(Y2):
  they are D8*·X1, D8*·Y1, D8*·X2, D8*·Y2 (D8*·X0 and D8*·Y0 lie in J_9 since (X0-Y0)G_u, (X0+Y0)G_v in J).
  Consequence: at degree 9 the E1 kernel is J_9 + D8*·S_1; the mod-J_18 quadratic reduction has dim R_18 = 456
  (~30 per sector) against ~41-58 Jacobian parameters per sector, so it cannot close degree 9 by itself.
- Gauge: the family of E1-admissible D always contains s*F (s in S_{d-6}), which satisfies the identity exactly
  with N = A s Z/6, Q = A s^2; closure runs now fix the gauge (coefficients of s*F set to 0), otherwise the
  quadric variety is positive-dimensional by construction (this is why the first Singular run timed out).

## 2. Understanding of the obstruction and of the mechanism a certificate needs

(a) Hodge-theoretic position of the target.  X is a smooth sextic CY 4-fold (Omega/F is the (4,0) form).
    Under the visible group Gamma (order 96) the residue A Omega/F^2 is invariant under G_aff (incl. tau_Y, which
    negates (Y0,Y1,Y2)) and odd under both exchanges s0, r1.  Omega/F itself is tau_Y-ANTI-invariant.  Hence in the
    tau_Y-invariant sector H^4(X)^{tau_Y} = H^4(X/tau_Y) there is no (4,0)-part, and the generalized Hodge
    conjecture predicts that every tau_Y-invariant class -- in particular every Y-odd target of the challenge,
    and ours -- is supported on a divisor.  Existence is therefore expected; the problem is constructive.
(b) What kills the known constructions.  fqlx's Reid-Tai audit: every nontrivial element of Gamma^+ (stabilizer of the
    residue) has age >= 1, so X/H has canonical singularities and kappa = 0 for every H in Gamma^+: no uniruled
    quotient, hence no Jacobian/moving-curve integration on a quotient of X itself.  The companion quadric-contact-
    broken used the reflection s0 (det -1, Omega/F -> -Omega/F, A' -> -A'): its target is s0-EVEN as a form and lives on the
    Fano quotient X/s0 in P(1,1,1,1,1,2).  Our target is s0-odd as a form (A even, Omega odd), so it is in the
    anti-invariant part of the double cover X -> X/s0, equivalently forms on X/s0 with a half pole along the branch
    divisor u = (X0-Y0)^2 = 0.  The same holds for r1.  Any descent through (S0,P0) carries the Jacobian factor
    delta = X0^2-Y0^2, so A = P0^3 is never a quotient density: this is the content of why.md.
(c) The K3 fibration.  Passing to passive (alpha,beta) = (Y1/X1, Y2/X2) gives the K3 family S_{alpha,beta} in
    P(1,1,2,2) (double cover of the companion's cubic surface branched along the two plane cubics S0^2 = 4 P0^2).
    fqlx: the restriction of the target to the generic fibre is a (1,1)-class whose Gauss-Manin derivative has a
    nonzero H^{0,2}-component, so no certificate can be (alpha,beta)-passive.  A divisor W fibred in curves over the
    (alpha,beta)-plane produces (3,1)-classes whose fibrewise part is a varying combination of ALGEBRAIC classes
    [C_t] (its GM derivative stays of type (1,1)), consistent with fqlx's criterion only when the base of the
    curve family is a cover B' -> P^1 x P^1 with p_g(B') > 0 and the coefficients vary -- this is exactly how the
    companion's 27-line norm works (the pull-back of a line to the K3 is a (-2)-curve: a line meets the two branch
    cubics in 1+1 points, so the double cover of the line is rational), so its (2,0)-forms come from the degree-27
    "27 lines" cover of the parameter plane, not from the fibres.  For OUR target the analogous fibrewise class
    would have to acquire an H^{0,2} GM-component, impossible for algebraic fibre classes: so no divisor whose
    intersection with the generic K3 fibre is a curve class pulled back from the cubic surface can work, and the
    vector field must move alpha or beta (fqlx).  A serving divisor must be "transverse" to this K3 structure.
(d) The cover/deck mechanism (holdout family, book 2.7) is the only accepted technique that reaches an exchange-even
    target: a cover T^6 = Y0^6 + Y2^6 created a new exchange y <-> T whose quotient is a uniruled Fano fibration,
    and the physical target equals a quotient density up to a correction of nontrivial deck character, killed by
    the trace.  Here the pair (X1,Y1) is ALREADY a Fermat pair with exchange r1 (F depends on it only through
    X1^6+Y1^6 and P1 = X1Y1), and X/r1 in P(1,1,1,1,2,2) is a Fano 4-fold of index 2 -- but, with trivial deck group,
    it serves only r1-odd forms delta_1*B.  The Fermat variables available for a merging cover are X2, Y2 (X1, Y1):
    T^6 = X2^6+Y2^6 leaves T coupled only through T^2, T^4 (no new exchange); T^6 = X1^6+Y1^6 similarly.  The (X0,Y0)
    block S0^3 + 2 S0 P0^2 = (X0^2+Y0^2)(X0^4+4X0^2Y0^2+Y0^4) is not a sum of two sixth powers (its six roots are not a
    regular hexagon: angular gaps 35.3, 35.3, 109.4 degrees), so no linear change makes it Fermat.  Conclusion: no
    "merge two Fermat variables" cover creates a reflection entangled with the deck group here; a serving
    correspondence would have to be of a different type (see section 4).
(e) Algebraic side.  Certificates may be averaged over Gamma (the identity is linear in the field v = N/D for fixed
    (F,q,A), and HasCertificate.add multiplies denominators), so WLOG D is Gamma-invariant, of even degree.  E1
    (A D in J) is a strong filter only up to degree 9: ker(x A: R_d -> R_{d+6}) is 0 for d <= 7, is <D8*> at d = 8
    (D8* = const * G_u G_v, a universal feature of F = G((X0-Y0)^2,(X0+Y0)^2,...) with A = ((v-u)/4)^3), is D8*.S_1 at
    d = 9, and from d = 10 it is "generic" (dim R_10 - dim R_8 + 1 = 556 overall).  E2 modulo J_{2d} is decisive only
    while R_{2d} is large: dim R_16[chi_A] = 52 at d = 8, dim R_18 = 456 at d = 9, dim R_20[chi_A] = 12 at d = 10,
    dim R_24[chi_A] = 0 at d = 12.  So exhaustive algebraic exclusion ends at degree 8 (invariant) / pure sectors;
    beyond that only structured candidates can be tested.
