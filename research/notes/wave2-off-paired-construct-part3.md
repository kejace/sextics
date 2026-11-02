# Wave 2 — off-paired — CONSTRUCTION lens — part 3: the one new sufficient criterion, and the no-hub theorem

## E. Descent through a reflection (variant of book §4.6; proof included)
Setting: X smooth projective, xi ∈ H^4(X,Q), rho ∈ Aut(X) with rho^* xi = −xi (here rho = sigma2 or any mirror sigma2 h^k).
Lemma.  Let p: Z → X be proper, dominant, generically finite, Z smooth projective, e ∈ Aut(Z) with p∘e = p (deck),
sigma ∈ Aut(Z) an involution, and suppose d := sigma e sigma satisfies p∘d = rho∘p.  Let R ⊂ Aut(Z) be a group commuting
with sigma whose elements act trivially on p^*xi (e.g. lifts of elements of G).  If a smooth projective model of Z/<sigma, R>
has divisor-supported H^4, then xi ∈ N^1 H^4(X).
Proof.  a := p^*xi, V := H^4(Z)/N^1.  a + sigma^*a is <sigma,R>-invariant, hence in N^1 (quotient transfer, book §3.1.4), so
sigma^*a ≡ −a.  e^*a = a exactly.  d^*a = sigma^* e^* sigma^* a ≡ sigma^* e^*(−a) = −sigma^* a ≡ a.  But d^*a = (p∘d)^*xi =
(rho∘p)^*xi = p^*(−xi) = −a.  So 2a ∈ N^1, hence a ∈ N^1 and xi = p_*a/deg p ∈ N^1.  ∎
Remark.  Unlike §4.6 no supported quotient X/<g> is needed: the sign of rho does the averaging.  For off-paired this is the
only descent-type criterion compatible with the facts of part 2 §D (no quotient of X by G-elements is uniruled).

Realisation.  W := X/sigma2 = {H + C m + s^6 − 6 s^4 m + 9 s^2 m^2 − 2m^3 = 0} ⊂ P(1,1,1,1,1,2), s = X2+Y2, m = X2Y2 (Fano,
−K = O(1); K(X) = K(W)(√B), B = s^2 − 4m = (X2−Y2)^2).  For a birational involution phi of W put
   Z := X ×_W X^phi = {(x1,x2) ∈ X × X : pi(x1) = phi(pi(x2))},  e = (1, sigma2), d = (sigma2, 1), sigma(x1,x2) = (x2,x1),
   p = pr_1.   Then p∘e = p, sigma e sigma = d, p∘d = sigma2∘p.  ✓ (all hypotheses of the lemma except support)
   r = (tau, tau) (tau = g^3 acts trivially on the target class) is admissible for R when phi commutes with the image of tau.
Obstruction for phi an automorphism of W (the centraliser of sigma2 in <G,sigma2> is <g^2, h^3>; these lift to X so X^phi ≅ X and
Z is reducible — useless).  For a hypothetical linear involution phi of P(1,1,1,1,1,2) preserving W and not lifting to X:
   K_Z = pi^*(K_W + B/2 + phi^*B/2) = pi^*O(−1 + 1 + 1) = pi^*O(1), p_g(Z) = 2 + h^0(W,O(1)) = 7, and the sigma-symmetric
   combination of the two Ω/F-pullbacks survives in Z/<sigma>; with r it dies, but the (−,−)-type forms O(1)·(√B√phi^*B)^{-1}
   include r-invariant ones (r acts on O(1) with eigenvalue +1 on X0,X1,X2-type sections) ⇒ p_g(Z/<sigma,r>) ≥ 1: NOT uniruled.
   So phi must be genuinely birational and must make the branch loci "small" (the holdout analogue works because the
   (y,T)-Klein quotient has weight-2 invariants S,P and the fibre equation is a sextic in P(1,1,1,2,2) over a CURVE).
Status: no phi found.  Precise requirement for a search: a birational involution phi of the genus-4 fibration W → P^3_x
   (fibres s^2(s^2−3m)^2 = −(H+Cm−2m^3), genus 4 = sigma2-quotients of the plane-sextic fibres; W → Y_- is the h^2-quotient)
   not induced by <g^2,h^3>, commuting with tau, such that the degree-4 cover Z/<sigma> → W/phi (tower W/phi ← (W/phi)(√(B·phi^*B))
   ← Z/<sigma>) has a uniruled quotient by r.  Fibrewise this asks for an involution of the genus-4 curve W_b (or of a cover
   of it) exchanging the two "sheets" X2 ↔ Y2 with a deck-like transformation — the genus-4 curve W_b = X_b/sigma2 has the
   extra involution h^3 only (s ↦ −s), which lifts to X.  Verdict: FAILS for every automorphism; open for birational phi.

## F. No-hub theorem (why every §2.7-type cover fails; exact statement)
Write F = ℓ^6 + ℓ·L + M in a coordinate system containing a linear form ℓ.  The §2.7 mechanism needs (i) no ℓ^k terms for
2 ≤ k ≤ 5, (ii) M − f6(X) a binary sextic FORM in two further coordinates (so T^6 := M − f6 defines a constant base curve),
(iii) L divisible by T after the substitution.  For off-paired the Y-quadratic part Q_X(Y) = (a/2)Y0Y1 + (b/3)Y0Y2 + (c/2)Y1Y2
has rank 3 over C(X), so (i) forces ℓ to be Q_X-isotropic, i.e. ℓ ∈ {Y0,Y1,Y2} up to the X-dependent isotropic cone, and then
M − f6 = (sum of the other two Y^6) + (the mixed term not containing ℓ) is inhomogeneous of degrees (6,2): (ii) fails.
Concretely the obstruction to each hub is one mixed monomial: hub Y0 ← X1^3X2Y1Y2/2 (second link), hub Y1 ← X0^2X1X2Y0Y2/3
(extra term), hub Y2 ← X0^3X1Y0Y1/2 (first link).  X-side hubs fail because F has an X1^3 term.  The exchange Y0 ↔ T with
T^6 = Y1^6 + Y2^6 + (c/2)Y1Y2 makes F formally symmetric but does not preserve the cover relation
alpha^6 + beta^6 + (c/2) alpha beta T^{−4} = 1 (verified by direct substitution), so sigma is not an automorphism of Z.
This is sharper than "no coordinate pairing works": no linear hub at all.

## G. Galois closures / averaging ideal (precise)
X → X/D6 = P(1,1,1,1,2) is Galois with group D6 = <h, sigma2>; the target is the det character, trivial exactly on the
rotation subgroups <h^k>; the quotients X/<h>, X/<h^2>, X/<h^3> all have p_g = 1 (det h = 1 on Ω/F), so no intermediate
quotient of this closure serves the target (averaging ideal of §3.1 empty).  The τ-cubic cover X/τ → P^2×P^2 is book §4.7.2.
