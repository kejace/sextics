# wave2 construct — klein-coupled-character — part 3: quotient companions, Galois closures, deformation, and the
# best remaining candidate (Geiser-twisted fibre product) with the exact place where it fails

## 3.1 (ii) Quotients by subgroups of the full symmetry group (Aut(X) = Z/14, part 2.2; nothing to compute remotely)
Subgroups: 1, <tau>, <sigma>, <sigma^k tau> = <sigma,tau> (k != 0, order 14). Fixed loci: sigma^k: the plane Fermat curve
{X=0, G(Y)=0} (normal weights (1,2,4)k/7, age 1) and e_x,e_y,e_z (ages >= 1); tau: the curves {Y=0,K=0} and {X=0,G=0}
(age 3/2); sigma^3 tau: age 1/2 at e_x,e_y,e_z (junior).  Kollar-Larsen:
 * X/<sigma>: no junior element -> canonical, K ~ 0, kappa = 0, NOT uniruled; h^{3,1} = dim R_6[weight 0] = 38 + 19 = 57 != 0.
   Serves only weight 0 anyway (averaging).
 * X/<tau>: kappa = 0, h^{3,1} = dim R_6^- = 208; H^4(X)^tau contains the whole target sector: supporting it IS the problem.
 * X/<sigma,tau> = W (book 4.3): junior sigma^3 tau -> uniruled; Fano model in P(2,2,2,1,1,1), -K = O(3), Hodge (19,128,19);
   averaging ideal = weight-0 projector only. Dossier 2.4 stands.
 => no quotient companion (Fano, surface-fibration or otherwise) of X itself can reach weight 3.

## 3.2 (iii) Galois closures of the natural finite maps and the averaging-ideal (constituent) test of part 1.1
 (a) X -> W: already Galois (Z/14); intermediate quotients = 3.1.
 (b) Projection from a point c on X (degree 5 onto P^4; from e_x it is <sigma,tau>-equivariant): Galois closure Z with group
     Gamma <= S_5 x Z/14, Delta = S_4 x 1, Ind_N^Gamma(chi~) = (1 + std_4) (x) chi. The std (x) chi part needs K with
     (std (x) chi)^K != 0, i.e. K meets S_5 in a subgroup with >= 2 orbits on the 5 sheets (S_4, S_3 x S_2, S_3, ...),
     and Z/K supported. Those quotients are the (components of the) fibre powers X x_{P^4} X minus diagonals and their
     quotients: K_{fibre product} = pi^* O(5) + (diagonal correction) on the residual component, i.e. general type; the A_5
     double cover of P^4 is branched along the degree-30 discriminant (K = O(10)). Rational quotients only for K >= S_5-type,
     where (std (x) chi)^K = 0. => the test fails for every intermediate quotient.
 (c) Fibrations (non-finite maps): Stein factorisation reduces to connected-fibre fibrations; for the Y-ratio fibrations
     (sigma trivial on the base) part 2.2 shows Aut(geometric generic fibre) = <sigma,tau>, so every automorphism of the
     Galois closure that is fibrewise descends to X; the constituent test then reduces to 3.1 (Mackey). For bases on which
     sigma acts (X-ratios, mixed invariants), an exchange e = s'ds' would descend to an element acting on the base, so it
     cannot be sigma^k tau^j restricted to fibres; no such construction reaches a weight != 0 class.
 (d) Polynomial Piola maps (book 3.4): sigma-eigen Phi of weight profile (2,2,2,1,1,1)-type cannot exist (each sigma-eigen
     quadric in x,y,z is a monomial; xy^3, x^3z, yz^3 need six distinct monomials); dossier item 7 (non-eigen, mixed) remains
     open as before — not a character-changing mechanism in any case (delta.b(Phi) targets are tied to the finite map).

## 3.3 (iv) Deformation, q-raising, Gauss-Manin
 * The only sub-family of the 4.3 family with a sector theorem for weight 3 is the split locus Q_0=Q_1=Q_2=0 (three
   Y-reflections, codim 9 in the 49-dim family). No explicit certificate for x^3 t^3 exists there (atlas: Y-degree-1 targets
   only), so the first-order continuation system (ORD1 of wave2-theory-deformation) has no t=0 input; divisor-supportedness
   of a fixed class is not preserved along F_s = K + G + sM (Hodge-locus type condition), and 4.5 q-raising needs a
   certificate at the same member. Nothing transfers. (v) no companion on the same variety.

## 3.4 The single most promising mechanism: a non-lifting involution of W and the twisted fibre product (part 1.2)
 Set-up. W -> P^2_Y has cubic-surface fibres C_eta: p(U,V,W) + T(aUV + bUW + cVW) + e T^3 = 0 over K = k(P^2_Y),
 T = Y0^2, a = (5/7)Y1/Y0, b = (10/7)Y1Y2/Y0^2, c = (15/7)Y2/Y0, e = G/Y0^6. K-rational points: e_U,e_V,e_W (p vanishes
 there). GEISER INVOLUTION from e_U (gradient of P at e_U is (0,0,1,0)):
     gamma: (U:V:W:T) -> (-(UW + V^2 + aTV + bTW)/W : V : W : T),  i.e.  U -> -U - (V^2 + (5/7)V tr + (10/7)W rs)/W
 in weighted coordinates (fixes V,W,Y). Check: gamma^2 = id. (Same with e_V, e_W by the cyclic structure.)
 Non-lifting (RIGOROUS): k(X) = k(W)(x), x^14 = h := U^5 W^3 / V (adjugate of the exponent matrix). gamma lifts iff
 gamma^* h = h^j . (14th power); div(gamma^*h / h^j) = 5[Q] - 5j[U] + 3(1-j)[W] - (1-j)[V], [Q] = {UW + V^2 + aTV + bTW = 0}
 irreducible and different from [U],[V],[W]; coefficient 5 of [Q] is never 0 mod 14. So gamma is a birational involution
 of W = X/<sigma tau> that does not lift: exactly the ingredient 1.2 asks for, with g = sigma tau (sigma-part != 1).
 Explicit auxiliary equations. Z_gamma = X x_{W,gamma} X = { (P,P') : pi(P') = gamma(pi(P)) }; as a cover of X it is the
 Kummer cover  w^14 = phi^5 x^3 y^9,  phi := z^5 x + (15/7) y z^3 t s + t^6 + r^6 + s^6  (so pi^*(gamma U) = phi/(x^3 z));
 branch divisor on X: {phi = 0} (mult 5), {x = 0} (3), {y = 0} (9). Group: d (deck), e (lift of sigma tau), swap s',
 r = tau x tau. Constituent test (1.1): passes with K = <s'> or <r,s'>, even enlarged by <d^7,e^7>.
 WHERE IT FAILS: the hypothesis "Z/K divisor-supported".
   - K = <s'>: omega = p_1^* Omega_X + p_2^* Omega_X is a nonzero swap-invariant holomorphic 4-form on a resolution
     (the two pullbacks have d-eigenvalues +1, -1); an involution quotient is canonical, so h^{4,0}(Z/<s'>) >= 1: NOT
     supported (coniveau). RIGOROUS.
   - K = <r,s'> (r kills that form, as in holdout): the Esnault-Viehweg eigensheaf count H^0(X, O(3j - 6[5j/14] - [3j/14]
     - [9j/14])) (j = 0..13; degrees 0,3,5,2,4,5,2,4,6,3,4,6,3,5) decomposed by the (Z/14)^2-characters (j,k) with the
     r-sign gives 948 (r = (tau, w -> w)) resp. 1053 (r = (tau, w -> -w)) dimensions of r-invariant off-diagonal pieces
     V_(j,k), j != k; each pair V_(j,k) + V_(k,j) contributes a swap-invariant form. Script w2c_geiser_h40.py. Caveat: the
     count is taken on X (branch divisor not snc), so it is an upper bound on h^{4,0}(Z~); but the sheer size (vs. the
     19-dimensional obstruction-free holdout picture where the fibrewise 3-forms are swap-ANTI-invariant) leaves no room
     for a Fano-type quotient: the quotient is not a Fano fibration over a curve (Z_gamma has no fibration over a curve
     compatible with the swap: the base P^2_Y is 2-dimensional and the fibres of Z_gamma/<r,s'> over it are degree-98
     covers of the rational surface C_eta/gamma, not rational).
 Why this is structurally different from holdout: there the exchange acted on the fibres of X -> P^1 (a linear swap of
 two fibre coordinates, invariants S,P polynomial, quotient fibre a weighted Fano 3-fold). Here part 2.2 proves no
 exchange acts on any Y-ratio fibre of X (Aut = Z/14 exactly), so the exchange has to live on a cover of X branched
 along a sextic-type divisor; such covers have p_g >> 0 and the swap cannot be made to act by -1 on all of it.
 Consequence for the Lean interface: no auxiliary uniruled/Fano quotient reaching weight 3 exists in this class, hence no
 finite-field incidence seed (polynomial curve family with invertible Jacobian) can be produced for this target.

## 3.5 Status of the construction lens (summary of rigor)
 RIGOROUS: no 2.7 linear exchange (mod p; char-0 replay pending for the e_t,e_r,e_s uniqueness); Aut of the geometric generic
 surface fibre = Z/14 (exact); Aut of the threefold fibres V_(1,2) = Z/14 (exact in charts x,y,z; t,T pending);
 Aut(X) = Z/14 (exact in charts x,t; y,z,r,s pending — mod p complete); Geiser gamma non-lifting; Z_gamma/<swap> not supported.
 EVIDENCE: Z_gamma/<r,swap> not supported (eigensheaf count on X).
 OPEN: non-base-change covers of other shapes; birational (Bertini-type) involutions of W from K-rational curves give the
 same obstruction pattern (Kummer cover of X branched along a positive-degree divisor => many (4,0)-forms).

## 3.6 Addendum: a rigorous lower bound for the <r,swap> quotient (r = tau x tau, i.e. r = (Y -> -Y, w -> w))
 The j = 1 eigenspace V_1 = { s . w^{-1} Omega_Z : s in S_3 } (deg_1 = 3, no strict-transform floors since 5,3,9 < 14) consists
 of forms that are regular on the normalised cover (vanishing orders along the ramification divisors: 13-5 = 8 over phi,
 13-3 = 10 over x, 13-9 = 4 over y) and extend over every exceptional divisor E of a log resolution of the branch divisor
 (EV condition floor(m_E/14) <= a_E: m_E = 8,14,12 for pairwise intersections with a_E = 1, m_E = 17 for the triple locus with
 a_E = 2). r acts on s w^{-1} Omega_Z by (-1)^{1+deg_Y s}: 28 r-invariant forms (Y-degree 1 or 3), of which 25 lie in (1,k)
 pieces with k != 1 ((1,3):3, (1,5):3, (1,7):10, (1,9):3, (1,11):3, (1,13):3). Each gives the nonzero <r,swap>-invariant
 form v + swap^* v; (Z/2)^2-quotient singularities are canonical (only reflections are junior among involutions), so
   h^{4,0}( Z_gamma / <r,swap> ~ ) >= 25  (and >= 948 if the full count of part 3.4 extends).
 Hence Z_gamma/<r,swap> is NOT divisor-supported, for the symmetric lift r = tau x tau; the other lift (tau,1) does not
 commute with the swap and generates <tau x 1, 1 x tau, swap>, whose invariants are a subspace of the above with the same
 j = 1 analysis (28 r'-invariant forms, 27 off-diagonal).  The Geiser candidate is dead for every admissible K.
