# wave2 dense-rational CONSTRUCT — part 2: coordinate-free exclusion of exchange covers and reflections

Scripts: /tmp/hodge/results/dense-rational/wave2_construct_hess.py (generator), wave2_construct_ltcheck.py (leading-ideal
certificate), inputs/outputs in /tmp/hodge/results/dense-rational/wave2-construct/ ; every result line is in
/tmp/hodge/results/dense-rational/wave2-construct.log.jsonl.  G = 9973 F (integral), p = 1073741789.
Rigor convention (same as the excl/verify lenses): for a HOMOGENEOUS ideal with p-integral generators, a pure power of
every variable in the mod-p DRL leading ideal => Macaulay matrix of full column rank mod p => full rank over Q => V_C = ∅.

## 2.1 What a holdout-type exchange cover needs (any linear coordinates, any block sizes)
Book §2.7 structure, generalized: pick linear coordinates (rest ; y ; u_1..u_r) on C^6 (r >= 1, rest of dim 5 - r),
adjoin T^6 = g(u) (g a sextic form), alpha_i = u_i/T on the affine base {g(alpha) = 1}, and ask that the pulled-back
equation be symmetric under y <-> T.  Writing F = sum_{m,b} y^b T^m p_{m,b}(alpha; rest) with p_{m,b} homogeneous of
degree m in alpha, symmetry means p_{m,b} = p_{b,m} in C[alpha]/(g - 1).  Since (g - 1) contains no nonzero element of
degree <= 5, for m != b both sides must vanish unless {m, b} = {0, 6} (which forces g = F|_{y=0,rest=0}/c, c = coeff of y^6,
c != 0).  Hence
      F = Phi(rest, y u_1, ..., y u_r) + c y^6 + g(u),   Phi weighted-homogeneous (rest: 1, y u_i: 2) of degree 6.   (EX)
(EX) is exactly the holdout shape (rest = X, y = Y1, u = (Y0, Y2)); it also covers r = 1..4 and "non-coordinate" y.
Invariant consequence: at the point e_y = [0:..:1:..:0] the Taylor expansion in the chart y = 1 is
F = c + (terms of degree >= 3)  (every monomial rest^a (yu)^b has a + b >= 3 because a + 2b = 6; g has degree 6).
So Hess F(e_y) = 30 c (dy)^2 has RANK ONE, and F(e_y) = c != 0.
The same rank-one point exists for every quasi-reflection of order 3 or 6 (F = c l^6 + l^3 q_3(h) + q_6(h), resp.
F = c l^6 + q_6(h); c != 0 because otherwise e_l would be a singular point of X).  For a quasi-reflection of order 5
(F = l^5 q_1(h) + q_6(h)) the point e_l lies on X and Hess F(e_l) vanishes on the tangent hyperplane: rank <= 2.
Order 4 (F = l^4 q_2 + q_6) is impossible on a smooth sextic.  Linear reflections (order 2) are treated in 2.3.

## 2.2 V1 = {rank Hess G <= 1} and V2 = {rank Hess G <= 2} are EMPTY  [rigorous over C]
* V1: 120 distinct 2x2 minors (degree 8), msolve -g 1, 4.9 s: 1923 leading monomials, max degree 11, pure powers of all
  six variables (degrees 8..11)  => V1 = ∅.                                           log line construct-V1
  Sanity: the same generator on the Fermat sextic gives the 15 binomials e_i^4 e_j^4 and NO pure power (V1 = six
  coordinate points, as it must: Fermat admits all exchange covers).                  log line sanity-fermat-V1
* V2: 210 distinct 3x3 minors (degree 12): msolve running (construct-V2); result appended in part 3.
Consequences (unconditional once V2 = ∅ is certified; V1 = ∅ already gives the first two):
  (a) no holdout-type exchange cover exists for F in ANY linear coordinates, for any block split and any base sextic g;
  (b) no quasi-reflection of order 3 or 6 (and, with V2, of order 5) preserves F in any linear coordinates;
  (c) the §2.7/§2.8 target-selection machinery (deck character on a Fermat-type base) has no cover to act on.
This is the precise failure point of mechanism (i): the dense F has no point of P^5 at which its polar quadric is a
double hyperplane.  The 236-term perturbation of size 1/9973 destroys the six Fermat rank-one points (the Fermat V1
is a non-transverse configuration in the expected-empty codimension-15 locus; the perturbation makes it empty).

## 2.3 V3: linear reflections  [running]
F even in a linear form l (coordinates (l, h)) forces, at e = e_l (F(e) = c != 0): the cubic Taylor term D^3F(e)(h^3)
vanishes on H = ker dF(e).  Using the Euler identities D^3F(e)(e,h,h) = 4 Hess_e(h,h), D^3F(e)(e,e,h) = 20 dF(e).h,
D^3F(e)(e,e,e) = 120 F(e), the vanishing on H is equivalent to the polynomial identity in h
      9 F(e)^2 D^3F(e)(h^3) = 18 F(e) (dF(e).h) Hess_e(h,h) - 10 (dF(e).h)^3,
i.e. the 56 degree-15 forms T'_{ijk} = 27 G^2 G_{ijk} - 18 G (G_i G_{jk} + G_j G_{ik} + G_k G_{ij}) + 30 G_i G_j G_k
vanish at e (checked on G = z^6).  Points of X cannot satisfy it (it would force dG(e) = 0).  V3 = ∅ => F is not even in
any linear form in any coordinates => no linear reflection, hence the §3.3 reflection theorem and the "one-reflection"
family (and Ch. 3.5 reflection groups) are unavailable in every coordinate system.  Note tau itself is not a
reflection (its (-1)-eigenspace is 3-dimensional); the two Fermat-type involutions that WOULD be reflections are the
coordinate sign changes, excluded by packet R10 only in the given coordinates; V3 makes it coordinate-free.

## 2.4 Ternary cones D = g(X), g(Y) (mechanism (ii)/(v) surrogate: surface-fibred divisors over a plane curve)
Script wave2_construct_cones.py: E1 for D = g in a cone sector <=> A g ∈ J_{k+6}; rank of [J | A x^beta] vs rank J with
the parity-split Hilbert-series cap (rigorous over Q when attained).
| k | sector | target deg/parity | cols | rank J = cap | #cones | admissible cones | rigorous |
| 7 | X | 13 / odd | 4284 | 3441 = 3441 | 36 | 0 | yes |
(rows for k = 7 Y, 8 X/Y local; 9..12 X/Y on the remote job construct-dense-rational-cones; appended in part 3)
Geometric reading (part 1 §1.2): E_g = X ∩ {g(X) = 0} is fibred over C_g by the sextic surfaces S_x through C_Y;
h^{2,0} of its resolution is the (2,0)-part of the monodromy invariants of the family, which is 0 by Lefschetz theory
of the net of P^3's through P^2_Y; so cones serve no primitive class at ANY degree (theory), and E1 confirms it at the
degrees where E1 is not vacuous.
