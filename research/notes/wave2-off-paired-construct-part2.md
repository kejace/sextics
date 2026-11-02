# Wave 2 — off-paired — CONSTRUCTION lens — part 2: the fibre-product picture and what it forbids

All identities below were verified exactly (c1_identities.py) and the surface data by Singular (c2*, c4, c6*);
scripts and outputs in /tmp/hodge/results/off-paired/wave2-construct/.

## A. The three fibrations of the (X2,Y2)-block (exact)
With m = X2Y2, on X:
  (X2^3+Y2^3)^2 = −(H + Cm − 2m^3),  (X2^3−Y2^3)^2 = −(H + Cm + 2m^3),  (X2^6−Y2^6)^2 = (H+Cm)^2 − 4m^6,
  (H+Cu)^2 − 4u^6 = (H + Cu − 2u^3)(H + Cu + 2u^3),   disc_u(2u^3 − Cu − H) ∝ 2C^3 − 27H^2.
Consequences (all over P^3_x, x = (X0,X1,Y0,Y1)):
  Y_- = X/<h^2,sigma2> : t^2 = 2m^3 − Cm − H   Weierstrass elliptic fibration, Fano (sextic in P(1,1,1,1,2,3));
  Y_+ = X/<h^2,sigma2 h^3>: t^2 = −2m^3 − Cm − H  likewise;  Ŷ = X/<h^2> = Y_+ ×_{P(1,1,1,1,2)} Y_-;
  X/<h> : genus-2 fibration v^2 = (H+Cm)^2 − 4m^6 (the thread's), Jac(fibre of Ŷ) ~ E_+ × E_- × J(E_b).
The target sector = sigma2-anti-invariant part of H^4(X/<h>) = the (−,−) character of the Klein cover Ŷ → P(1,1,1,1,2);
the other two characters are H^4 of Fano fourfolds (served); the trivial character is H^4 of P(1,1,1,1,2).
Sector Hodge numbers (c3_sectors.py, Fermat-ring count): h-invariant Y-odd primitive = (56,144,56), ENTIRELY
sigma2-anti-invariant; split by g-character: g=0: (22,50,22) [target], g=2 and g=4: (17,47,17) each.

## B. The genus-2 surface Σ is rational; its base changes
Σ = {v^2 = (H+Cu)^2 − 4u^6} ⊂ P(1,2,3,3) has exactly 2A1 + 2A2 (c2c_types), K = O(−3): rational.  Explicit
parametrisation (w = Cu/H, lambda = 1/kappa = H^2/C^3):
   lambda = (1+w)(1−s^2) / (2w^3(1+s^2)),   v/(H+Cu) = 2s/(1+s^2);
fibre E_lambda = the (3,2)-curve (1+w)(1−s^2) = 2 lambda w^3 (1+s^2), hyperelliptic involution s ↦ −s, g: (lambda,w,s) ↦ (−lambda,w,1/s).
Hence H^1(P^1_kappa, V) (V = R^1 of the genus-2 family, rank 4) is pure (1,1): NO (2,0)-class comes from the genus-2
variation alone.  Base changes kappa = rho^n (double-plane models v^2 = (rho^b + rho^a u)^2 − 4u^6, 3a−2b = n; c6b_types):
   n : 1  2  3  4  5  6        p_g(Σ̃_n) : 0  0  0  1  2  1
   (n=2: branch sextic with a J10-type quadruple point (mu=15) + A5 ⇒ 1−1 = 0;  n=3: quadruple point mu=19 ⇒ 0;
    n=4: A5 + A11 ⇒ K3;  n=5: octic branch, J10 + A13 + A5 ⇒ 3−1 = 2;  n=6: A17 ⇒ K3.)
   These match the Hodge-bundle formula deg f_*ω(Σ̃_n) = 2n − ⌊n/3⌋ − ⌊2n/3⌋ − ⌊n/2⌋, p_g = deg − 2, i.e. local exponents
   (1/3,2/3) at kappa=0 (orbifold point of P(3,2), fibre smooth) and (1/2,0)+unipotent at kappa=∞ (H=0, nodal fibre).
   So twisted (2,0)-classes H^0(P^1, Ω^1(log) ⊗ V^{1,0} ⊗ chi) exist only for characters chi of order ≥ 4 ramified at 0,∞
   (first: order 4 and 6, one each; order 5: two), never for order 2 or 3.

## C. Leray analysis of the sector (X/<h> ≅ Σ ×_{P^1_kappa} B̃ birationally; W := H^2(S̃_kappa), S̃_kappa quasi-smooth
   (4,6) c.i. in P(1,1,1,1,2), p_g = 45, q = 0 — c4_Skappa)
   H^4(X/<h>)^− ≈ H^1(P^1_kappa, V ⊗ W) (+ boundary terms).  Divisor supports, by type:
   (V1) vertical D = pullback of a surface Λ ⊂ P^3: Gysin image = H^1(P^1, V ⊗ [components of Λ∩S_kappa]) — only the
        ALGEBRAIC classes of the pencil surfaces enter, twisted by the permutation local system of the components.
        kappa constant on Λ (H, C, Δ-components, every divisor the thread used): image 0 — this is the geometric content
        of deanodotdotdot's "every component must move the modulus" and of benediamond's contact obstruction.
        kappa non-constant with Λ ∩ S_kappa irreducible: H^1(P^1,V)⊗[Γ] is pure (2,2) (Σ rational): Tate only.
        kappa acquiring an n-th root on Λ (n ≥ 4 needed by §B): the conjugate components must be NON-homologous on S̃_kappa.
        For every Λ of the shape {f(C,H) = mu(x)^n} (C = mu^2, H = nu^3, Δ = mu^2, C^3H^2 = mu^4, ...) the components are
        hypersurface sections mu = ζ·(function of s) of S̃_kappa ⊂ P(1,1,1,1,2), hence linearly equivalent: image 0.
        What survives: curves on the generic S̃_kappa not cut by hypersurfaces, i.e. a Picard local system of the pencil
        with non-trivial finite monodromy.  After the base change kappa = mu^{−2} the pencil is {H = mu s^3} ⊂ T, T = {s^2 = C}
        the (singular) quartic double solid; Lefschetz-type pencils have irreducible monodromy on the vanishing part, so
        this is expected to fail, but it is the one remaining loophole of type V1 (not decided here).
   (H1) horizontal D pulled back from a curve Λ ⊂ Σ (block divisors a X2^6 + b X2^3Y2^3 + c Y2^6, mirrors, X2^6 ± Y2^6,
        torsion multisections of Y_±): Gysin image = cup of the (1,1) normal-function class of the multisection in
        H^1(P^1, V ⊗ chi^{-1}) with H^0(P^1, W ⊗ chi) — needs FLAT chi-twisted (2,0)-forms of the pencil (a finite-monodromy
        transcendental sub-VHS of W).  None expected (big monodromy); consistent with the search lens: all such D serve 0.
   (M)  mixed D (fibrewise a surface D_kappa ⊂ E_kappa × S_kappa dominating both): needs H^{1,0}(D̃_kappa) beyond H^{1,0}(E)
        pairing into H^1(E) ⊗ H^2_tr(S), i.e. a correspondence making J(E_kappa) ⊗ (something) a sub-VHS of W_tr.
        This is the only type compatible with the GHC prediction that the (22,50,22) sector IS supported; it is exactly the
        kind of object the holdout cover produces (same block family F = H + C X2Y2 + X2^6 + Y2^6, there with a hub).
   Control: on the Fermat (C ≡ 0) the cone over C_X is of type (M)-like through the Ext^1 component (join classes), which is
        why the horizontal/vertical dichotomy is NOT a theorem: the Leray graded map vanishes but the extension class does not.
   NOTE the earlier draft claim "horizontal divisors cannot serve (−)-classes" is FALSE (Fermat join is a counterexample);
        the correct statements are the three bullets above.

## D. Reid–Tai: no quotient of X by a subgroup of G is uniruled (rigorous)
   Every gamma ∈ G∖1 has age ≥ 1 at every fixed point on X (dossier 1.5, referee-confirmed), G has no pseudo-reflections,
   det(gamma) ∈ {±1}.  Hence X/K has canonical singularities for all K ≤ G and (Ω/F)^{⊗2} is K-invariant and extends to any
   resolution: P_2(X/K) ≥ 1, kappa(X/K) ≥ 0.  The same holds for (X × Y)/K whenever Y carries a K-invariant bicanonical form
   (elliptic curves, abelian/K3/Enriques carriers).  Only p_g = q = 0 surface fibrations of X/K remain; fqlx excluded all
   monomial ones.  For K containing a reflection the target is not K-invariant.  (Replaces "not known uniruled" by "not uniruled".)

## B'. Hodge-bundle prediction for all base changes kappa = rho^n (formula of §B; verified n ≤ 6 by double planes)
| n | deg f_*omega(Σ̃_n) | p_g(Σ̃_n) | new (2,0) in characters of exact order n |
|---|---|---|---|
| 1 | 2 | 0 | 0 |
| 2 | 2 | 0 | 0 |
| 3 | 2 | 0 | 0 |
| 4 | 3 | 1 | 1 |
| 5 | 4 | 2 | 2 |
| 6 | 3 | 1 | 1 |
| 7 | 5 | 3 | 3 |
| 8 | 5 | 3 | 2 |
| 9 | 5 | 3 | 3 |
| 10 | 6 | 4 | 2 |
| 11 | 7 | 5 | 5 |
| 12 | 6 | 4 | 2 |
Reading: a vertical divisor on which kappa acquires an n-th root can only carry twisted (2,0)-classes for n ∈ {4,5,6,7,8,...}
(never n = 2, 3), and even then only if the conjugate fibre components are non-homologous on S̃_kappa (part 2 §C).
