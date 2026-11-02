# Wave 2 — off-paired — CONSTRUCTION lens — part 1: target, structure, mechanism inventory

(Started 2026-11-01 ~22:00 local; continues the killed attempt, which had produced only a header.)

## 1. Target and symmetry (from dossier, verified by the referee note)
F = Fermat + X0^3X1Y0Y1/2 + X1^3X2Y1Y2/2 + X0^2X1X2Y0Y2/3, A = X0^3Y0^3, q = 1.
Block form (x = (X0,X1,Y0,Y1)):  F = H(x) + C(x)·X2Y2 + X2^6 + Y2^6,
  H = X0^6+X1^6+Y0^6+Y1^6 + X0^3X1Y0Y1/2 (deg 6),  C = X1^3Y1/2 + X0^2X1Y0/3 (deg 4).
Diagonal group G = <g,h> = Z/6 x Z/6, g=(0,2,0,1,3,3), h=(0,0,1,0,0,5); sigma2: X2<->Y2.
Residue of A: G-invariant, sigma2-anti-invariant. Sector H^4(X/G)_prim = (22,50,22).
Under the dihedral group D6 = <h, sigma2> (order 12) acting on the (X2,Y2)-block the target carries the
determinant character (trivial on rotations h^k, -1 on the six reflections sigma2 h^k).

## 2. Structural facts derived in this session (to be verified exactly in part 2)
S1. X/D6 = P(1,1,1,1,2) (coordinates x, m = X2Y2): the D6-invariants of the block are m and
    I6 = X2^6+Y2^6, and F = 0 is LINEAR in I6.  So the D6-invariant sector is trivial, as expected.
S2. Klein cover.  Ŷ := X/<h^2> is the (Z/2)^2-cover of P(1,1,1,1,2) given by
      t1^2 = H + C m − 2 m^3,   t2^2 = H + C m + 2 m^3,      t1 = i(X2^3+Y2^3), t2 = i(X2^3−Y2^3).
    Its three intermediate double covers are
      Y_- = X/<h^2, sigma2>      : t1^2 = H + Cm − 2m^3  ⊂ P(1,1,1,1,2,3)   (degree 6, −K = O(3): FANO),
      Y_+ = X/<h^2, sigma2 h^3>  : t2^2 = H + Cm + 2m^3  ⊂ P(1,1,1,1,2,3)   (FANO),
      X/<h>                      : v^2 = (H+Cm)^2 − 4m^6 ⊂ P(1,1,1,1,2,6)   (K = 0, p_g = 1),
    and the target sector is exactly the (−,−) character of this Klein cover (the sigma2-anti-invariant
    part of H^4(X/<h>)).  Both branch threefolds B± = {H + Cm ± 2m^3 = 0} ⊂ P(1,1,1,1,2) are sextic
    Calabi–Yau threefolds = (reflection-fixed threefolds X ∩ {X2 = ±Y2}) / (z ↦ −z).
S3. Elliptic fibrations.  Over P^3_x the Fano fourfolds Y_± are WEIERSTRASS elliptic fibrations
      t^2 = ∓2 m^3 − C m − H   (m weight 2, t weight 3),  discriminant ∝ 2C^3 − 27H^2 (resp. with H→−H),
    j-invariant a function of kappa = C^3/H^2 only.  The fibre of Ŷ → P^3 is the genus-4 curve
    E+ ×_{P^1_m} E− (étale Z/3 quotient of the plane-sextic fibre of X), with Jacobian ~ E+ × E− × J(E_b),
    where E_b: v^2 = (H+Cm)^2 − 4m^6 is the genus-2 fibre of X/<h> → P^3 (the thread's genus-2 fibration).
    The target lives in the J(E_b)-part: H^3(P^3, V), V = R^1 of the genus-2 family (rank 4).
S4. The genus-2 family depends only on kappa.  Its total space over the kappa-line,
      Σ = {v^2 = (H+Cu)^2 − 4u^6} ⊂ P(1,2,3,3) (coordinates u,C,H,v), K_Σ = O(−3),
    is a RATIONAL surface: with w = Cu/H, lambda = 1/kappa, the branch locus (1+w)^2 = 4 lambda^2 w^6 is the pair
    of graphs lambda = ±(1+w)/(2w^3), and v/(1+w) =: v'', z := 2 lambda w^3/(1+w) satisfy v''^2 + z^2 = 1.
    Rational parametrisation: z = (1−s^2)/(1+s^2), v'' = 2s/(1+s^2); the genus-2 fibre over lambda is the
    (3,2)-curve  E_lambda : (1+w)(1−s^2) = 2 lambda w^3 (1+s^2)  in P^1_w × P^1_s, hyperelliptic involution s ↦ −s,
    and g acts by (lambda, w, s) ↦ (−lambda, w, 1/s).  Hence H^1(P^1_kappa, V) is pure (1,1): the "transcendental
    part of the genus-2 surface" that the thread's currents use is algebraic; all (3,1)-classes of the sector
    involve the pencil of surfaces S_kappa = {C^3 = kappa H^2} ⊂ P^3 (normalisation: the (4,6) complete
    intersection {C = kappa s^2, H = kappa s^3} ⊂ P(1,1,1,1,2), K = O(4), p_g ≈ 45) or monodromy-twisted
    base changes Σ_n = Σ ×_{P^1} (kappa = rho^n), which are (n,3,2)-surfaces in (P^1)^3 with K = O(n−2,1,0),
    so p_g(Σ_n) = 2(n−1) if the singularities are rational, every 2-form in a NONTRIVIAL deck character.
S5. X/<h> is birationally the fibre product Σ ×_{P^1_kappa} B̃ (B̃ = P^3 with the pencil resolved).

## 3. Mechanism inventory (task items i–v) — verdicts so far
(i)  Cover creating an exchange (book §2.7).  A hub variable ℓ needs F = ℓ^6 + ℓ·L + M with M − f6 a
     binary sextic form in the two other Y's; Q_X(Y) has rank 3 (no hub), so every coordinate hub leaves an
     inhomogeneous quadratic term (hub Y0: X1^3X2 Y1Y2/2 = second chain link; hub Y1: the extra term;
     hub Y2: the first link), and a non-coordinate linear hub would have to be Q_X-isotropic over C(X)
     AND leave a binary form — impossible for a rank-3 form.  FAILS (precisely: the exchange Y0 ↔ T does
     not preserve the cover relation alpha^6+beta^6 + c·alpha·beta·T^{-4} = 1).
     NEW variant (descent through a reflection, part 2 §D): a cover Z → X with deck e and an involution
     sigma with sigma·e·sigma = lift of a reflection acting by −1 on the target makes the target supported
     as soon as Z/<sigma, r> is supported — no supported quotient of X is needed.  Realisation
     Z = X ×_W X^phi (W = X/sigma2 Fano, phi a birational involution of W not lifting to X).  Open: no phi found.
(ii) Quotient companions.  NEW theorem (Reid–Tai): every element of G has age ≥ 1 at every fixed point
     (dossier 1.5, confirmed by the referee) and det(gamma)^2 = 1, so for every K ≤ G the quotient X/K has
     canonical singularities and (Omega/F)^{⊗2} descends: P_2(X/K) ≥ 1, hence X/K is NOT uniruled.
     Twisted companions (X × Y)/K with P_2(Y) ≥ 1 (elliptic curves, K3, abelian) likewise.  FAILS rigorously
     for every honest or elliptic-twisted quotient; only p_g = q = 0 surface fibrations remain (fqlx excluded all
     monomial ones).
(iii) Galois closures.  The natural Galois map is X → X/D6 = P(1,1,1,1,2) with group D6; the target is the
     det character, served only by rotation quotients X/<h^k>, all with p_g = 1.  The cubic cover X/tau → P^2×P^2
     is the book's §4.7.2 (V_std unserved).  FAILS (averaging-ideal test negative on every intermediate quotient).
(iv) Deformation/q-raising.  Sub-family with a sector theorem: the chain (t = 0 in F_t = chain + t·extra),
     every odd target served by §4.6; the mechanism dies at first order (unmatched T^2 term) and §4.5 only
     raises q on a fixed F.  No certificate to continue.  FAILS as a construction; gives no obstruction either.
(v)  No solved companion on this variety.  N/A.
