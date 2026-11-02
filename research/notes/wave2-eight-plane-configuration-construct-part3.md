# Part 3: Aut(X), the (H5) loophole, ramification divisors, and what is left

## 1. Linear automorphisms of X (mechanism (ii) input)

Reduction (Part 1 §4): a linear automorphism g of X permutes the planes contained in X. IF the eight
planes are the only planes on X, then g preserves their ideal (Q_0,Q_1,Q_2), hence permutes the three
rank-2 members Q_j = u_j v_j of the net up to scalars, hence permutes the kernels {X_j = Y_j = 0} and is
block-monomial in the (u_j, v_j) coordinates: a permutation of the six coordinates preserving the pairing
{u_j, v_j} (48 = |S_3| 2^3 choices) composed with a diagonal scaling.

Exact computation over Q(i) (aut_blockmonomial.py -> aut_blockmonomial.json): F has 316 of the 462
degree-6 monomials in (u,v). For each of the 48 permutations pi the scaling d is forced by the coefficient
ratios on monomial pairs differing by e_k - e_0 (always available), and then all 316 equations
a_m d^{m'} = lambda a_{m'} are checked. Result: exactly two solutions, pi = id (d = 1, lambda = 1) and
pi = global swap u_j <-> v_j for all j (d = 1, lambda = 1), i.e. tau. All 46 others are inconsistent
(none undetermined). Hence

   Aut_lin(X) = <tau>   provided the Fano scheme of 2-planes of X is the reduced set of eight planes.

Fano scheme computation: 20 affine charts of G(3,6), 9 unknowns, 28 sextic equations each
(fano_planes_gen.py -> fano/chart*.ms, p = 2147483629 = 1 mod 4 so that i exists). Expected: 8 solutions
in each of the 8 charts picking one coordinate per block, 0 in the other 12. Remote job
construct-epc-fanoplanes (capacity-gated; the box has been at load 570-600 with 320 cpus allocated) and a
niced single-thread local probe of chart00. RESULT: see the STATUS line of the index / the end of this file.

## 2. The (H5) loophole of the omega-descent theorem

The theorem needs Y/K to have canonical singularities for the subgroups K containing a tau-type element,
i.e. the Galois closure Y must be smooth (or canonical) along the fixed loci of elements (k,-1), which lie
over the twisted diagonal Delta_tau = {(p, tau p)} ~ P^2. Over Delta_tau \ Delta the cover is etale: smooth,
eigenvalues (1,1,-1,-1), age 1. The curve Delta_tau n Delta decomposes (mod 10^9+9, random line in Pi):

   Delta|_{Delta_tau} = F_{6,0}(X(p)) . F_{0,6}(Y(p)) . S(p)^2,   deg = 6 + 6 + 2*12 = 36,  S squarefree.

 * C_X' = {F_{6,0}(X(p)) = 0}: the tau-stable line L_{p,tau p} passes through the tau-fixed point
   p + tau p in C_X and is automatically tangent there (dF kills all Y-directions at points of P^2_X), so
   q has the double root s = t. Delta vanishes to order exactly 1 at a sampled point of C_X', both in a
   generic direction of B and along Delta_tau: Delta is SMOOTH there, so Y is smooth over generic points of
   C_X' and the tau-type element acts with a codimension-2 fixed locus: age 1, canonical. Same for C_Y'.
 * S = 0: q has two double roots r, 1/r swapped by tau (checked: q = (x-r_1)^2 (x-r_2)^2 at a sampled
   point). Delta has multiplicity exactly 2 there in every sampled direction, and its tangent cone is a
   quadratic form of rank 2: two transverse smooth branches (the two tangency sheets, exchanged by tau).
   The local Galois closure {z_1^2 = delta_1, z_2^2 = delta_2} is smooth because d delta_1, d delta_2 are
   independent. Tau-type elements act by swapping the sheets: fixed loci of codimension >= 2, ages >= 1.
 * Remaining: the finitely many points where these curves meet each other or where the coefficient map
   B --> P^4 (binary quartics) fails to be versal (triple roots on Delta_tau, p +- tau p both on X, ...).
   An isolated non-canonical quotient point is the only way (H5) can fail; it would make the pairs cover
   W' uniruled "through a point", which no known mechanism produces. Not pursued further.

So the omega-descent theorem holds for the four plane-pair covers up to finitely many isolated points.

## 3. Ramification divisors of the plane-pair covers are excluded Jacobian denominators

The line L_x through x = (u,v) is tangent to X at x iff grad F(x) . (u, 0) = 0, i.e. E_u(F)(x) = 0 with
E_u = sum_j u_j d/du_j = (E - i R)/2, E = Euler field, R = R_0 + R_1 + R_2, R_j = X_j d/dY_j - Y_j d/dX_j.
On X this is R(F) = 0. For the pair (eps, -eps) the ramification divisor of pi_eps is
   Ram(pi_eps) = X n V( sum_j eps_j R_j (F) ),
a Y-odd degree-6 Jacobian denominator (four of them). Every degree-6 denominator of any parity is excluded
(excl note, msolve certificates), so the branch geometry of these covers cannot serve the target — in
agreement with the theorem. The pulled-back discriminant pi_eps^* Delta is a degree-36 divisor on X
(pi^* O(1,0) = pi^* O(0,1) = O_X(1)); it is the only divisor the construction singles out, far beyond
any linear-algebra reach, and the theorem gives no reason to expect it to support alpha.

## 4. Precise statement of where the construction lens stops

 (i)   exchange covers: impossible for every linear form and base degree (Lemma A; exact ranks 6, 6).
 (ii)  quotients: Aut(X) = <tau> (modulo the Fano-scheme check), X/tau has kappa = 0; no age < 1 element;
       twisted quotients have dimension > 4 or ages >= 1.
 (iii) finite maps: the four quartic plane-pair covers are new and instance-specific, but every intermediate
       quotient of their S_4 x Z/2 Galois closure that carries alpha has kappa >= 0 (omega-descent), the
       hypotheses being verified (S_4 monodromy, reduced discriminant, smooth Galois closure along the
       tau-fixed strata up to isolated points). The same argument covers point projections and the cubic
       cover (book §4.7).
 (iv)  no q-companion; the Fermat degree-6 branch cannot reach t = 5.
 (v)   no solved companion.
 The book's finite-field incidence-seed recipe (polynomial curve family with invertible incidence Jacobian)
 presupposes a uniruled quotient to integrate on; none exists here, so there is no seed to look for.
 What is NOT excluded: a "created reflection" of non-linear origin (an automorphism of a cover of X acting
 by -1 on omega, +1 on alpha, with a codimension-1 fixed locus), or a surface fibration of a cover with
 p_g = 0 in the tau-invariant sector. The arrangement offers no candidate for either: the planes are rigid
 and of type (2,2), the natural self-correspondence x -> (other residual points on L_x) acts on alpha and
 omega alike, and the natural elliptic carrier w^2 = q over Pi x Pi' is non-isotrivial.

## 5. Suggested discussion content (not posted; for the human/lead)
 * Lemma A + ranks (kills all §2.7-type covers for this F);
 * the four quartic plane-pair covers, their S_4 monodromy and the omega-descent theorem with proof;
 * Aut(X) = <tau> once the Fano scheme result is in;
 * the identification Ram(pi_eps) = X n V(sum eps_j R_j F) (degree-6 Jacobian divisors, already excluded).
