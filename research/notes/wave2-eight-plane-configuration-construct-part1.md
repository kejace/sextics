# Wave 2 — eight-plane-configuration — CONSTRUCTION lens — Part 1: framework and mechanism inventory

Lane: `eight-plane-configuration`, q = 1, A = X0^3 Y0^3, F = sum_i Q_i K_i, Q_i = X_i^2 + Y_i^2.
Target class alpha = Res(A Omega/F^2) in H^{3,1}(X): tau-INVARIANT (A odd, Omega odd), sector dim 208.
Inherited rigorous facts (excl/verify notes): deg D >= 7 for every certificate; d in {7,8,9} forces
D = w(F) Jacobian type; d = 10 is the first degree with non-Jacobian E1 kernel (294 even / 261 odd).
Work dir: /tmp/hodge/results/eight-plane-configuration/wave2-construct/ ; log wave2-construct.log.jsonl.

## 1. What a construction must do (book Ch. 2-3 distilled)

The book's only support engine is the covering-curves theorem (Ch.3 §2.1): if a smooth projective
4-fold W is uniruled, ALL of H^4(W) is divisor-supported, and the divisor pulls back along any
correspondence to X. Every accepted construction therefore produces a uniruled 4-fold W with a
correspondence Gamma: W -> X whose image contains alpha, in one of three ways:
 (Q) W = X/G for a finite G with alpha G-invariant and X/G uniruled (reflection quotient, Fano weighted
     hypersurface, or an element of age < 1 — Kollar–Larsen);
 (C) W = Z/<sigma, tau~> for a cover Z -> X (base change along a curve, §2.7) with a CREATED reflection
     sigma on Z; fibrewise Fano quotient; polynomial-curve integration + trace;
 (S) a surface fibration W -> B with p_g = q = 0 fibres in the relevant sector (Ch.3 §2.3).

Key invariant: X is Calabi–Yau (K_X = 0) with holomorphic 4-form omega = Res(Omega/F), and
tau^* omega = -omega while tau^* alpha = +alpha. Any W that "sees" alpha but still carries omega
(or a descendant omega^2 as a 2-canonical form) has kappa >= 0 and is NOT uniruled. So every
construction must act DIFFERENTLY on omega and on alpha. On X itself only tau does so, and tau's
fixed locus has codim 3 (age 3/2), so omega^2 descends to X/tau: kappa(X/tau) = 0 (book §4.7).

## 2. Mechanism inventory for this lane

| # | mechanism | status | where it fails |
|---|---|---|---|
| (i) | §2.7-type cover creating an exchange y <-> T, T^m = beta(w), any linear forms | EXCLUDED (Lemma A + exact rank computation, §3) | F_{4,2} and F_{2,4} coefficient quadrics span 6 dims (max); every one has Gram rank 3 |
| (ii) | quotients X/G, twisted quotients (X x Y)/G, Fano/surface companions | EXCLUDED for G = <tau> (book §4.7: kappa(X/tau) = 0); twisted quotients have dim > 4 or inherit age >= 1; Aut(X) = <tau> expected, Fano-scheme check running (Part 3) | no element acts with age < 1; tau is not a reflection |
| (iii) | Galois closures of the natural finite maps: (a) book's cubic cover of X/tau (fails, §4.7), (b) NEW: the four quartic covers X -> Pi_eps x Pi_{-eps} from antipodal plane pairs | (b) EXCLUDED by the omega-descent theorem (Part 2), hypotheses verified computationally | every intermediate quotient that carries alpha has kappa >= 0 |
| (iv) | q-raising / Gauss–Manin; deformation from the coordinatewise-even sub-family (Fermat) along F_t = Fermat + t sum Q_i H_i / p_i | NOT USABLE: the sub-family is codim 160+ in the arrangement family; the degree-6 Fermat witness cannot continue to t = 5 inside degree 6 (dossier §2.3, excl: no deg-6 D at t=5) | no q-companion exists (q = 1 only instance on this F) |
| (v) | modify a solved companion's support geometry | N/A | no solved companion on this variety |

## 3. Lemma A (necessary conditions for an exchange cover) and its verification

Setting (book §2.7 generalised). Split the 6 linear coordinates as S (3 "spectators") + ell + W
(dim W = 2), let T^m = beta_m(w) define the base curve, and suppose F = Phi(S, ell, T) with
Phi symmetric under ell <-> T and coefficients in K_b = C(base curve). Write
F = sum_a S^a sum_j ell^j P_{6-a-j}(w) (P homogeneous of degree 6-a-j in w).

Claim. Polynomiality + symmetry force, for every spectator degree a and j < (6-a)/2:
 P_{6-a-j}(w) = T^{6-a-2j} P_j(w) as functions on the base, hence (irreducible base of degree m)
 * m >= 3: P_{6-a-j} = P_j = 0 unless 6-a-2j = 0 mod m;  * m = 2: P_{6-a-j} = beta_2^{(6-a-2j)/2} P_j.
Proof: the coefficient c of S^a ell^j T^k (k = 6-a-j) must satisfy c T^k in C[w]_k, so c = P_k(alpha,beta)
is a form of degree k in the base coordinates (alpha,beta) = w/T; symmetry c_{j,k} = c_{k,j} on the
base curve means P_k - gamma^{k-j} P_j vanishes on {beta_m(alpha,beta) = gamma^m}; a form of degree k < m
in the ideal (beta_m - gamma^m) is zero, and for m = 2 the relation is multiplication by beta_2 powers.

Consequence for a = 4 (the (S^4, (ell,W)^2) part of F): every quadric q_j(ell,w) attached to a spectator
quartic monomial lies in  C(ell^2 + beta_2(w)) + ell W*  (m = 2, dim 3)  or in  ell W*  (m >= 3, dim 2).
So the 15 quadrics span at most 3 dimensions, and each is DEGENERATE (contains the factor ell when
m >= 3; rank <= 3 is not enough for m = 2 but span <= 3 is).

tau-compatibility. For the tau-invariant target the book's quotient is by <sigma, tau~> with tau~ a lift of
tau commuting with sigma; this forces S and ell + W to be tau-eigenspaces: S = X-block (exchange inside the
Y-block) or S = Y-block (exchange inside the X-block). Both are covered by the two tests below.
(Non-tau-compatible splittings cannot use tau~ to kill omega; they are left as a formal caveat.)

Computation (exact over Q, mech_i_ranks.py, result mech_i_ranks.json):
 F_{4,2}: 15 quadrics in Y (coefficients of the X-quartic monomials): span rank 6 = dim Sym^2(C^3)^*;
          Gram ranks all 3 (each quadric is irreducible, no linear factor).
 F_{2,4}: 15 quadrics in X (coefficients of the Y-quartic monomials): span rank 6; Gram ranks all 3.
Verdict: no exchange cover of §2.7 type exists for ANY choice of linear forms (ell, W) in either block and ANY
base degree m. This is rigorous (ranks over Q, Lemma A).

## 4. Mechanism (ii): quotients and twisted quotients

* Diagonal symmetries mod scalars = <tau> (SNF, dossier §1.3, verified). Non-diagonal: a linear automorphism
  of X permutes the planes on X; if the 8 planes are ALL planes on X (Fano-scheme computation, Part 3),
  then g preserves the ideal (Q_0,Q_1,Q_2), hence permutes its three rank-2 members Q_i, hence is block-
  monomial in the (u_i, v_i) = (X_i +- i Y_i) coordinates (48 candidate permutations x torus); a finite
  exact check then decides Aut(X). Expected: Aut(X) = <tau>.
* X/tau: kappa = 0 (omega^2 descends; tau has age 3/2), not uniruled (book §4.7). No other honest quotient.
* Twisted quotients (X x Y)/<tau x iota> (permutation-theory framework): dimension 4 + dim Y > 4, so the
  covering-curves theorem does not address H^4; and the age of tau x iota at any fixed point is
  >= age_X(tau) = 3/2 >= 1, so Kollar–Larsen never makes them uniruled. The book's elliptic companions
  need a 3-dimensional "variety part" (split F), absent here.

## 5. Mechanism (iv)

The only sub-family of the 162-dimensional arrangement family {sum Q_i K_i} where a sector theorem holds is
the coordinatewise-even locus (reflections present), e.g. Fermat at t = 0 of F_t = Fermat + t G. Known:
at t = 5 no degree-6 denominator of any parity exists (excl note, rigorous), so the Fermat degree-6 witness
cannot be continued at fixed degree. A Gauss–Manin/first-order test at t = 0 would only confirm where the
degree-6 branch dies, not produce a certificate; there is no q-companion to raise from. Not pursued further.

(Part 2: the quartic plane-pair covers and the omega-descent theorem. Part 3: Fano scheme / Aut(X) and the
remaining loophole analysis.)
