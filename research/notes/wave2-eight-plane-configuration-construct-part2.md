# Part 2: the quartic plane-pair covers X -> Pi_eps x Pi_{-eps} and the omega-descent obstruction

This is the one mechanism that exists ONLY because of the eight planes (the book's §4.7 treats the cubic
cover X/tau -> P^2_X x P^2_Y available for every even sextic). It is mechanism (iii) of the task.

## 1. Explicit auxiliary equations

Coordinates u_j = X_j + i Y_j, v_j = X_j - i Y_j (j = 0,1,2), so Q_j = u_j v_j, F = sum_j u_j v_j K_j(u,v),
tau = simultaneous swap u_j <-> v_j. The eight planes are the coordinate 2-planes
  Pi_eps = { w_0 = w_1 = w_2 = 0 },  w_j = u_j (eps_j = +) or v_j (eps_j = -);  tau Pi_eps = Pi_{-eps}.
Fix an antipodal pair Pi = Pi_{(+++)} = {u = 0} = {(0, v)} and Pi' = Pi_{(---)} = {v = 0} = {(u', 0)}.
They are complementary (C^6 = Pi + Pi'), so every x in P^5 \ (Pi u Pi') lies on a unique line
  L_{p,p'} = { s p + t p' },  p = (0, v) in Pi,  p' = (u', 0) in Pi',   (s:t) in P^1,
giving the finite rational map  pi: X --> B := Pi x Pi' = P^2 x P^2.   On the line,
  Q_j(s p + t p') = s t u'_j v_j,   so   F(s p + t p') = s t * q_{p,p'}(s,t),
  q_{p,p'}(s,t) = sum_j u'_j v_j K_j(t u', s v) = sum_{a=0}^{4} c_a(u',v) t^a s^{4-a},
  c_a = sum_j u'_j v_j K_j^{(a,4-a)}(u',v),  K_j^{(a,b)} = part of K_j of bidegree (a,b) in (u,v),
  bideg c_a = (a+1, 5-a) in (u', v).      [In X,Y coordinates with eps = (+++):
  p = (a, i a), p' = (b, -i b), Q_j(sp+tp') = 4 s t a_j b_j, q = 4 sum_j a_j b_j K_j(s p + t p').]
So X n L = {p, p'} u {4 residual points}: pi has degree 4. The discriminant Delta = disc_{s,t}(q) is a
divisor of bidegree (18,18) on B. tau acts on B by (p, p') -> (tau p', tau p), fixed locus the twisted
diagonal Delta_tau = {(p, tau p)} ~ P^2 (dimension 2); on a tau-stable line L_{p, tau p} it acts by
s <-> t, and q is symmetric (verified). The two tau-fixed points of such a line, p + tau p = (2 X(p), 0)
and p - tau p = (0, 2 Y(p)), lie on X iff X(p) in C_X = {F_{6,0} = 0} resp. Y(p) in C_Y = {F_{0,6} = 0}
(the two tau-fixed plane sextic curves of X); for general p in Pi neither holds, so tau acts on the
four residual points of L_{p,tau p} as a double transposition without fixed points.

Galois closure. Let Y -> B be the Galois closure of pi (group G_0 = geometric monodromy), B' = B/tau,
Y' the Galois closure of X/tau -> B'. Since C(X) = C(X/tau) . C(B) (a point of X is recovered from its
tau-orbit and its image in B), C(Y) = C(Y') . C(B) and G := Gal(Y/B') is an extension of Z/2 by G_0.
With G_0 = S_4 (no outer automorphisms) G = S_4 x <tau~>, tau~ central, X = Y/S_3, X/tau = Y/(S_3 x Z/2),
B = Y/S_4, B' = Y/G.  Intermediate quotients Y/K, K <= G.

The "pairs cover" W' = Y/(S_2 x S_2 x <tau~>) (degree 6 over B': unordered pairs of residual points,
i.e. the double cover of the resolvent-cubic cover Y/D_4 -> B) is the only new 4-fold; all K with image
in a point stabiliser give covers of X/tau.

## 2. Computational verification of the hypotheses (quartic_cover.py -> quartic_cover.json)

 (H0) divisibility F(sp+tp') = s t q: identity above (exact); checked numerically at every sample.
 (H1) monodromy. Specialisation a = (3,-5,7), b = (2,11,-4), eps = (+++): q in Q(i)[s,t]; Frobenius cycle
      types at the 15 split primes 10^9+9 ... 10^9+579 (p = 1 mod 4): 4-cycle x4, (3,1) x8, (2,1,1) x1,
      (2,2) x2. A transitive group containing a 4-cycle and a 3-cycle is S_4, so Gal(q_b / Q(i)) = S_4,
      hence the arithmetic monodromy of pi over Q(i)(B) is S_4 and G_0 is normal in S_4 with transitive
      action: G_0 in {V_4, A_4, S_4}.
 (H2) Delta restricted to a random line of B (mod 10^9+9): degree 36 = 18 + 18 attained, gcd(Delta, Delta')
      = 1. So Delta is reduced and NOT a square. If G_0 were in A_4 then sqrt(Delta) would be a polynomial
      on B. Hence G_0 = S_4, and the generic point of Delta is a simple double root (local monodromy a
      transposition).
 (H3) non-isotriviality: I^3/J^2 of q takes three different values at three random points, so the
      elliptic curves {w^2 = q} over B (the natural "elliptic carrier") form a non-isotrivial family.
 (H4) twisted diagonal: q(s,t) = q(t,s) when p' = tau p (sampled), as predicted.
 Also used: X irreducible (smooth), B and B' rational.

## 3. Theorem (omega-descent obstruction for the plane-pair covers)

Let alpha in H^{3,1}(X)^{tau} be the target class and omega the holomorphic 4-form of X. For every
subgroup K <= G = S_4 x Z/2 such that the K-average of alpha is nonzero in H^4(Y)^K = H^4(Y/K), the
variety Y/K has Kodaira dimension >= 0 (assuming the canonical-singularity hypothesis (H5) below for the
one family of K where it is needed). In particular NO intermediate quotient of the Galois closure of
any of the four plane-pair covers is uniruled while carrying alpha, so the covering-curves theorem never
applies to them. The failure point is identical to the book's §4.7 result for the cubic cover.

Proof.
 (a) Representations. alpha in H^4(Y)^{S_3 x Z/2}, so its G-span is a quotient of Ind (1) = (1 + V_std) x 1;
     the trivial component is a pullback from B' = Sym^2 P^2 (rational), of type (3,1) there, hence 0.
     So span(alpha) = V_std x 1 (3-dim), and alpha corresponds to the S_3-fixed vector e_4 - (e_1+e_2+e_3)/3.
     omega in H^4(Y)^{S_3 x 1} with tau~ omega = -omega: its span is a quotient of (1 + V_std) x (1 + sgn)
     on which tau~ = -1, i.e. of (1 + V_std) x sgn; the 1 x sgn part would be an S_4-invariant 4-form,
     i.e. a 4-form on the rational B: zero. So span(omega) = V_std x sgn, with omega the same S_3-fixed vector.
 (b) Which K carry alpha. (V_std x 1)^K = V_std^{Kbar} (Kbar = image of K in S_4) is nonzero iff Kbar is
     intransitive, and the K-average of alpha is nonzero iff the projection of e_4 - (e_1+e_2+e_3)/3 to
     V_std^{Kbar} is nonzero (true for all intransitive Kbar: e.g. for S_2 x S_2 it is (e_3+e_4-e_1-e_2)/3).
     Intransitive Kbar: either Kbar fixes a point (then Y/K dominates X/tau, kappa(X/tau) = 0, and a
     dominant image of a uniruled variety is uniruled: so Y/K is not uniruled), or Kbar has orbits {1,2},{3,4}:
     Kbar in {<(12)(34)>, S_2 x S_2}. Only the latter need an argument.
 (c) Kbar in {<(12)(34)>, S_2 x S_2}, K <= Kbar x Z/2 an elementary abelian 2-group. Put K_0 = K n (S_4 x 1).
     Pick eta in V_std^{Kbar} \ 0 (exists by (b)); eta is a holomorphic 4-form on Y, K_0-invariant, and every
     element (k,-1) in K maps eta to -eta. Holomorphic forms invariant under a finite group descend to a
     resolution of the quotient, so h^{4,0}(Y/K_0) >= 1 and Y/K_0 is not uniruled. If K = K_0 we are done.
 (d) Otherwise Y/K = (Y/K_0)/iota with iota induced by some (k,-1). Fixed points of (k',-1) on Y lie over
     Fix(tau_B) = Delta_tau (dimension 2), so every such fixed locus has codimension >= 2, and iota is NOT a
     pseudo-reflection anywhere on Y/K_0: the quotient Y/K_0 -> Y/K is quasi-etale. Let R <= K_0 be the
     subgroup generated by the transpositions in K_0 (pseudo-reflections of Y along the ramification over
     Delta); their stabilisers are generated by reflections (Fix((12)(34)) = Fix((12)) n Fix((34)) for
     ordered 4-tuples), so Y/R is smooth (Chevalley), eta descends to a holomorphic 4-form on Y/R (an invariant
     form is divisible by the reflection coordinate), and K/R acts on Y/R with every nontrivial element an
     involution whose fixed locus has codimension >= 2 (age >= 1). By Reid–Tai, Y/K = (Y/R)/(K/R) has
     canonical singularities (H5), and eta^2 — invariant under K/R — descends to a nonzero section of
     omega^{[2]} on Y/K; canonical singularities let it extend to a resolution. Hence kappa(Y/K) >= 0.
 (e) Concretely, for the pairs cover W' = Y/(S_2 x S_2 x Z/2) (the only new candidate): eta =
     omega(x_1) + omega(x_2) - omega(x_3) - omega(x_4), anti-invariant under tau~, and eta^2 is a 2-canonical
     form on W'. QED.

 (H5) is the only place where the general-position input enters beyond (H1)-(H2): Y here is a smooth
 G-equivariant model of the normalisation of the fibre product; the normalisation is smooth over the
 locus where the quartic family is versal (at simple double roots this is exactly reducedness of Delta,
 (H2)); singular points can occur only over the codim >= 2 strata of Delta where the coefficient map
 B --> P^4 (space of binary quartics) is not etale. A pointwise-fixed exceptional divisor of a resolution
 over a 2-dimensional centre inside Fix((k,-1)) would turn (k,-1) into a pseudo-reflection; but kappa is
 birational and the descent of eta^2 can be done on the normalisation, where no such divisor exists, as
 long as its singularities along Delta_tau are canonical. That is the single loophole (see Part 3, §2).

## 4. Where exactly the mechanism fails, and what this means

* The failure is structural, not computational: because X is Calabi–Yau and both alpha and omega come
  from X, every element of G acts on span(alpha) and span(omega) by the same permutation matrix up to
  the sign of tau~. A quotient sees alpha iff it sees eta (= the omega-analogue), unless it contains a
  tau-type element, and tau-type elements have codimension-2 fixed loci, so they only produce eta^2.
* The same argument applies verbatim to the degree-5 projections from tau-fixed points of X
  (G = S_5 x Z/2, span(alpha) = V_std x 1, span(omega) = V_std x sgn), and to the book's cubic cover
  (S_3 x Z/2), recovering §4.7. It applies to every tau-equivariant generically finite dominant map
  X --> B with B, B/tau rational, modulo the singularity hypothesis (H5) for the tau-type elements.
* The elliptic carrier w^2 = q over B (H3) is non-isotrivial, so the twisted-product/elliptic-companion
  route (permutation-theory note, O2) would need the transcendental Hodge structure of a non-isotrivial
  family of elliptic curves over a rational 4-fold — no CM/abelian shortcut.
* What is NOT excluded by this theorem: a cover Z -> X with an automorphism acting DIFFERENTLY on alpha
  and omega with a codimension-1 fixed locus (the §2.7 "created reflection"). Part 1 Lemma A shows no
  such structure exists in the book's coordinate/linear-form form; a non-linear (e.g. Cremona-induced)
  created reflection is the only remaining shape, and nothing in the eight-plane geometry suggests one:
  the planes are rigid (h^0(N) = 0), of Hodge type (2,2), and the self-correspondence they induce
  (x -> the other three residual points on its line) acts on both alpha and omega by -1.
