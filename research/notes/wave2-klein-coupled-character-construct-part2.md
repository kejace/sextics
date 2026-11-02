# wave2 construct — klein-coupled-character — part 2: linear exchanges, fibrewise automorphisms, Aut(X)

Scripts: /tmp/hodge/results/klein-coupled-character/w2c/ (w2c_sigma_gen.py, parse_msolve*.py, w2c_monomial_aut2.py,
smooth_check.sing); outputs sigma*_*.out (msolve), log wave2-construct.log.jsonl. Prime for mod-p runs: 1000003;
"exact" = msolve in characteristic 0 with rational parametrization (-P 1), elim-polynomial degree = #points over Qbar.

## 2.1 The book-2.7 linear exchange, all linear forms and complements: a degenerate-polar condition
2.7 needs coordinates (s | u' | u'') with F = c s^6 + s * sum_j u''_j b_j(u') + C'(u') + N(u'')  (s <-> T swap on the cover
T^6 = N(u''), base {N = c}). Necessary: for the point p = e_s and the hyperplane H = ker(d_p^5 F), the polars d_p^k F
vanish on H for k = 2,3,4. By the polarization identity d_p^4 F(Z) ~ Z^T Hess F(p) Z and d_p^5 F(Z) ~ grad F(p).Z, the
k = 4 condition is exactly:   Hess_F(p) restricted to ker(grad F(p)) == 0   with F(p) != 0.            (DP)
The same condition (DP) with F(p) = 0 is "the tangent hyperplane section of X has multiplicity >= 3 at p".
Computation (unknowns p in a chart, auxiliary m with H_ij(p) = F_i(p) m_j + m_i F_j(p); 22 equations, 11-12 unknowns):
 * F(p) != 0 (Rabinowitsch): charts x,y,z EMPTY; charts t,r,s: ONE point each mod p, i.e. p in {e_t, e_r, e_s}.
   (mod 1000003; char-0 replay submitted to the remote: construct-kcc-sigmaXQ.)
 * For p = e_t: B_5 = (5/7)(xy^3 r + 3 yz^3 s) is linear in u'' = (r,s) as required, but the remainder
   C = K + (10/7) x^3 z r s + r^6 + s^6 does not split as C'(x,y,z) + N(r,s): the third coupling term is quadratic in u''.
   Same for e_r, e_s by the cyclic structure. Invariant form: M = (5/7)(xy^3 tr + 2x^3z rs + 3yz^3 ts) is a quadratic form
   in Y over k(X) with Gram determinant 2abc = 2(5/14)(5/7)(15/14) x^4 y^4 z^4 != 0 (rank 3), so M is divisible by no
   linear form in Y; 2.7's "mixed term = s * (linear in u'')" fails for EVERY linear form s(Y). On the X side the pure
   block K has no sixth power of a linear form at all (K has no point p with F(p) != 0 satisfying (DP): charts x,y,z empty).
 => RIGOROUS (mod p now, exact pending): no 2.7-type cover with a linear exchange exists for F, for any coordinate pairing
    or non-coordinate linear form, with any complement and any base form N.

## 2.2 Fibrewise automorphisms of the Y-ratio fibrations (base-change / dihedral mechanism of part 1.3)
A base-change cover Z = X x_B B' with an exchange not defined over k(B) requires Aut(geometric generic fibre) to be
strictly larger than Aut_{k(B)} = <sigma,tau>. Fibres: X -> P^2_Y has sextic surfaces S_Y c P^3 (general type, K = O(2),
Bir = Aut = linear); X -> P^1_{[r:s]} (and [t:r], [t:s]) has sextic threefolds V c P^4 (Pic = Z by Lefschetz, Bir = Aut = linear).
Tool: the (DP)-locus Sigma(V) = {triple points} is Aut-invariant and finite; points are separated by the multiplicity
of the tangent section (3 at e_x,e_y,e_z where the cubic term is z r s-type; 4 at the Fermat hyperflex points where the
cubic term vanishes and the quartic coupling survives). An automorphism permuting {e_x,e_y,e_z} and their tangent
hyperplanes {y=0},{z=0},{x=0} is block-monomial: g = (monomial on x,y,z) + (M on the remaining coordinates); the pure
Y-part (Fermat sextic / binary form t^6 + cT^6 / e lambda^6) forces M monomial too (Aut of a Fermat curve of degree >= 4
is (mu_n)^2 x| S_3). Monomial automorphisms: exponent-lattice + permutation test (w2c_monomial_aut2.py):
   S_(1,2,3): diagonal projective group Z/14 (hand check: (a,b,c)=(u,u,u)+k(1,2,4)N/7, d = u or u+N/2), the two
   support-preserving 3-cycles fail the coefficient test; V_(1,2): Z/14, no permutations; X: Z/14, 3-cycles fail.
Results:
 * SURFACE FIBRE S_(1,2,3) (smooth, Singular): Sigma(S) = {e_x,e_y,e_z} EXACTLY over Qbar (charts x,y,z: elim degree 1;
   chart lambda: empty). Hence Aut(S_(1,2,3)) = Z/14 = <sigma,tau>|_S. Matsusaka-Mumford (injective specialization of
   automorphisms for non-ruled smooth fibres) => Aut(S_eta-bar) = <sigma,tau> for the geometric generic fibre of X -> P^2_Y.
   THEOREM: the P^2_Y fibration carries no exchange (linear or birational) over any base cover; mechanism 1.3 is impossible there.
 * THREEFOLD FIBRE V_(1,2) (smooth): mod p: chart x: {e_x}; y: {e_y}; z: {e_z}; chart t: 6 points (t^6 = -65 T^6 Fermat
   points, each of ideal length 4); chart T: running; char-0 replay running (sigmaQP_V_*). Expected Sigma(V) = 9 points,
   giving Aut(V) = Z/14 and (M-M) no exchange on X -> P^1_{[r:s]}. Pairings [t:r], [t:s] (V2, V3): running over Q.
 * FOURFOLD X: Sigma(X) mod p = {e_x,e_y,e_z} (length 4 each) + 18 Fermat hyperflex points {x=y=z=0, one of t,r,s = 0,
   ratio a 12th root of unity with 6th power -1} (12 points per chart t,r,s, length 22 each). With the multiplicity
   separation (3 vs 4) this gives Aut(X) = Z/14 once the char-0 replay confirms Sigma(X)_Qbar = these 21 points
   (semicontinuity already gives finiteness over Q; only extra points reducing to the same residue classes mod p remain).
   This closes dossier 1.5's open item "complete classification of finite non-diagonal linear automorphisms".

## 2.3 Hodge-theoretic facts used for the fibrations (hand, from the Jacobian rings of the fibres)
 * S_Y: H^{2,0} = R_2 = all 10 quadrics in (x,y,z,lambda); sigma-weights x:1,y:2,z:4,lambda:0 => the zeta^3-eigenspace of
   H^{2,0}(S_Y) is span(xy), ONE-dimensional (tau-even). So the zeta^3-part of R^2 pi_* has a (2,0) piece: book 3.2's
   surface-fibration theorem cannot apply to this sector (dossier 5 item 4 closed negatively).
 * V (threefold fibres over P^1): H^{3,0} = R_1 = linear forms; no linear form has weight 3 => H^{3,0}(V)[zeta^3] = 0;
   the sector of H^3(V) is of type (2,1)+(1,2) only. A relative support theorem for this sector would need the zeta^3-part
   of H^3 of a sextic CY threefold to be supported on divisors (Abel-Jacobi type); no book mechanism gives that.
 * Dossier 1.5 correction (not load-bearing): sigma^3 tau has age 1/2 at the coordinate points (weights (2,1,1,1)/14 on
   the tangent directions (z,t,r,s)/x), so G = <sigma,tau> DOES contain junior elements; that is why W = X/G is Fano
   (-K = O(3)) while X/sigma and X/tau have kappa = 0. The audit's "juniorElementsFixingTarget: []" stays correct.
