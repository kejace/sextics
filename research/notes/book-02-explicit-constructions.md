# Technique map: Chapter 2, "Explicit constructions: transport, contact, and algebraic cycles"

Source: `/home/kejace/devel/sextics/docs/02-explicit-constructions.md` (969 lines, read completely).
Supporting context: `formal/README.md`, `docs/01-certificates.md`, `data/catalogue.json`, `data/<id>/instance.rec`.

Conventions used below (from Chapter 1 and the formal contract):

- Variables `Z=(X0,X1,X2,Y0,Y1,Y2)`, `tau(X,Y)=(X,-Y)`, `F` Y-even smooth sextic, target `A` Y-odd of degree `6q`.
- Field `v=N/D`, `deg D=d`, `deg N_i=d+6q-5`, `gcd(D,F)=1`, component parities Y-odd for `v_0..v_2`, Y-even for `v_3..v_5`.
- Cleared equation `E = A D^2 - q D N(F) + F (D div N - N(D)) = F^{q+1} Q`.
- Euler-normalized form: `E1: q N(F) = A D`, `E2: D div N - N(D) in (F^q)`. Strict means `E=0`.
- Addition of certificates: `(D,N,Q)+(E,M,R) = (DE, EN+DM, E^2Q+D^2R)` (Lean: `HodgeComputer.CertificateOperations`).
- Descent: normalized trace `Tr_{K/k}(v)/[K:k]` over any finite extension; denominator becomes the norm.

Legend for status tags used in this file:

- **[Thm]** stated and proved in the chapter (I checked the argument's shape, not every line).
- **[Verified]** polynomial identity I re-derived with sympy 1.12 in this session (script runs recorded at the end).
- **[Stated]** identity/number quoted from the chapter that I did not independently recompute.
- **[Evidence]** computational seed / specific numeric data the chapter reports (finite-field seeds, roots in intervals).
- **[Limit]** explicit scope statement or boundary the chapter itself asserts.
- **[My inference]** my own cross-reference, not a claim of the chapter.

Global framing the chapter gives: every explicit certificate is a *transport identity*: change variables so `F` becomes an equation with accessible geometry, build a rational field downstairs, pull it back by cofactors; the Jacobian determinant of the map decides which target is served. Geometry (cycle classes) is needed only to prove that the served targets span.

---

## 2.1 Cofactor (Piola) transport **[Thm]**

Hypotheses: `Phi` a polynomial or rational map between affine spaces of the same dimension (rational allowed on its domain; *every* coordinate, including parameters written as ratios, is differentiated). `J = dPhi`, `Delta = det J`, `w` a target-space rational field, `G` a target function.

Key identities:

```
V = adj(J) · w(Phi)
V(G(Phi)) = Delta · w(G)(Phi)
div V     = Delta · (div w)(Phi)
```

Proof shape: `Phi^*(iota_w Theta_t) = iota_V Theta` (minor expansion), then `d` of both sides; plus `J adj(J)=Delta I` and the chain rule. So a divergence-free `w` serving `B` downstairs serves `Delta · B(Phi)` upstairs ("determinant-density target").

**Spectral polynomial constructor** (q=1): give `Phi_i` positive weights `m_i`, `G` weighted degree 6, `w` a polynomial derivation of weighted degree 0, density `b` of weight `12 - sum m_i`. Suppose constants `mu, nu` with

```
D0 = w(G),   w(D0) = (mu+nu) D0 - mu nu G,   div(b w) = (mu+nu) b .
```

Then

```
F = G(Phi),  A = b(Phi) Delta,  D = D0(Phi),  N = adj(dPhi) · (b w)(Phi),
N(F) = A D,   D div N - N(D) = mu nu A F     (so E2 holds with quotient mu nu A, q=1).
```

Degrees: `deg A = 6`, `deg D = 6`, `deg N_i = 7`. Legality (`D != 0`, `gcd(D,F)=1`) is a separate check.

**Nilpotent branch**: `w` of weighted degree `e`, `w^2(G)=0`, `div(b w)=0` gives `N(F)=AD`, `div N = N(D) = 0`, `deg D = 6+e`, `deg N = 7+e`. For `e != 0` a constant-spectrum quadratic relation on `G` can only be of this nilpotent form (weights of the three terms differ).

**[Limit]** "These spectral hypotheses are sufficient conditions for explicit fields. They are not a characterization of certificates." Nonlinear polynomial `Phi` has a broader weighted-Fano interpretation (Chapter 3, §3.4); rational maps can introduce boundary terms; scope/necessity belong to Chapter 4.

## 2.2 Euler and Nambu fields

**Subset Euler [Thm]**: `E_I = sum_{i in I} Z_i d_i`. If `E_I^2(F) = 6 E_I(F)` and the monomial target satisfies `deg_I A + |I| = 6`, then `v = A E_I / E_I(F)` is strict at `q=1`; `A E_I / (q E_I(F))` works at any `q` under the same condition. Denominator `E_I(F)` (degree 6) must be legal. Finite sums serve spans.

**Nambu fields [Thm]**: homogeneous `P, h_1..h_4` with `sum (deg h_i - 1) = 2`,

```
iota_W (dZ0 ^ ... ^ dZ5) = dP ^ dh_1 ^ ... ^ dh_4 .
```

Closed five-form => `div W = 0`; `W` kills `P, h_1..h_4`. If `F = Gamma(h_1..h_4) + P K` then `A = W(K)`, field `W/P`, strict: `W(F) = P W(K)`, `div(W/P) = 0`. `deg W = deg P + 1`, so the field has degree 1. This is cofactor transport through `Phi = (K, P, h_1..h_4)` with target equation `t0 t1 + Gamma` and derivation `d_0`. Rational first integrals: same proof, but polynomiality of `A` and pole cancellation must be checked.

**Adjunction / residue interpretation [Thm]**: if `E = V(P,h_1..h_4)` is a zero-dimensional complete intersection and `d = deg P` (the five degrees sum to `d+6`), the Koszul resolution twisted by `d` shows that restriction of degree-`d` polynomials to `H^0(E, O_E(d))` has one-dimensional cokernel. The four `deg h_i` sum to 6 so `Gamma` has no linear term, hence `dF = K dP` on `E` and smoothness makes `K` a unit on `E`. Euler column reduction: `(A/K)|_E = (deg K) J_E` with `J_E` the affine CI Jacobian, and the residue functional has `Res(J_E) = length(E)`. So `A/K` spans the missing adjoint direction.
**[Limit]** "the contact selects a specific target through a residue pairing; it does not serve an arbitrary prescribed A."

## 2.3 Perfect-power contact and its source surface

Hypotheses: split ternary `F = P(U) + Q(V)` (`U = X`-block, `V = Y`-block), with homogeneous `D, H, K`:

```
Q = c H^m + D K,  m in {2,3,6},  c != 0,  k = deg H = 6/m,  W = grad D x grad H  (cross product in V).
```

Multiplier `B(U,T)` of weighted degree `3-k` with weights `(1,1,1,k)`. Field and target **[Thm]**:

```
v = B(U,H) W / D,    A = B(U,H) W(K),
W(D) = W(H) = div W = 0,  W(Q) = D W(K)   (strict, q=1).
```

Support `V(F,D)` maps to the cyclic surface `S = {P(U) + c T^m = 0} subset P(1,1,1,k)`, `T = H(V)`. Smoothness of `F` forces smoothness of `P=0` and of `S`; `K_S = O_S(3-k)`, whose sections are exactly the multipliers `B`: 1 constant (m=2), 3 linear (m=3), 10 quadratic (m=6). Source identification: `iota_W vol_V = dD ^ dH`; residue of `iota_W vol_V / D` along a reduced component of `D` is its multiplicity times `dH`; after contraction with the Euler field one gets `Res_S(B Omega_weighted/(P + cT^m))` up to orientation/multiplicity; Gysin map of the support. **[Limit]** "Nonvanishing still requires checking that target class."

Rational version: `H, K` rational with `Q = cH^m + DK` gives `v(F)=A`, `div v=0` in the function field; a certificate only if `A` is the required polynomial and all remaining poles are legal. **[Limit]** A root `H` on the normalization of `D` need not descend to an ambient polynomial (distinct nodal branches may have distinct root values); enforcing an ambient polynomial root excludes such contacts; allowing rational roots without checking extra poles admits false candidates.

Picard model **[Thm]**: cubic surface `S`, twisted cubic `Gamma`, `Gamma^2 = 1`, `r = [Gamma] - H_S`, `r^2 = -2`, `r·H_S = 0` (root in `E_6`); for a cyclic triple-plane quotient `S -> P^2` the invariant rational Picard space is `<H_S>`, so `Gamma` cannot have invariant class. Moral: an ambient divisor has invariant class; the correspondence can vanish in the required nontrivial character. "A rational Nambu field may need a second divergence-free field to cancel its poles away from the intended denominator" — an additional rational identity, not a consequence of the perfect-power relation.

## 2.4 Cubic-surface theorem (complete target map downstairs) **[Thm]**

Statement: `G` a smooth homogeneous cubic in four variables over characteristic-zero `K`. For every `q >= 1` and `B` homogeneous of degree `3q-1` there is a rational field `w_q` with components of degree `3q-3`, `q w_q(G) = B`, `div w_q = 0`. One common denominator of active degree 45, coprime to `G`, serves every `B` and `q`.

Inputs: cubic surface = blowup of `P^2` in 6 points, 27 lines, 45 tritangent planes (5 through each line); Villaflor Loyola Thm 1.1: for `G = l1 Q1 + l2 Q2`, `L = V(l1,l2)`, `Res(det d(l1,Q1,l2,Q2) Omega/G^2)` is a fixed nonzero scalar times `H_S - 3[L]`; Griffiths `R(G)_2 = H^2(S)_prim`.

Line cofactor field (reusable):

```
Phi_L = (l1, Q1, l2, Q2),  J_L = det(dPhi_L)  (quadratic),
w_L = adj(dPhi_L) e_2 / l1,   w_L(G) = J_L,  div w_L = 0   (degree 0 field).
```

(Downstairs `G = ab + cd`, field `(1/a) d_b`.) Six mutually skew lines: `r_i = [L_i] - H_S/3` have Gram matrix `-I_6 - (1/3) 1 1^T`, eigenvalues `-3, -1^5`, det 3, independent; together with the 4 partials `G_i` (independent, else `G` is a cone) they span the 10 quadratics. A `10x10` solve gives `B = sum a_i G_i + sum b_j J_{L_j}` and the strict field `sum a_i d_i + sum b_j w_{L_j}`.

Descent: `w-bar = Tr_{E/K}(w)/[E:K]`; the active derivatives commute with trace. Trace via multiplication matrices: component `N/D` has normalized trace `tr(adj(M_D) M_N) / ([E:K] det M_D)`. Legality: every conjugate denominator is a product of linear forms, which cannot contain the irreducible cubic.

Higher `q` recursion (reusable, **[Verified]** by hand: value and divergence both check):

```
B = sum b_i G_i (deg b_i = 3q-3; possible since K[u]/(G_0..G_3) has Hilbert series (1+t)^4),
c = div b,  choose w with (q-1) w(G) = c, div w = 0,
v_q = [ b + G w - c E_4 / (3(q-1)) ] / q .
```

Common denominator: pick `p in K^4` off the 45 tritangent planes, normalize each plane `l(p)=1`, `T = prod l` (degree 45, Galois-stable, in `K[u]`); use a tritangent plane through each chosen line as `l1`; numerator degree `45 + 3q - 3`.

### Finite-algebra realization of line contact **[Thm]**

`R` char-0 domain of passive coefficients, `E` commutative finite free `R`-algebra of rank `m` (defining polynomial may be reducible / have repeated roots). A line-contact field `n/l` over `E` with `n(G) = l J`, `div(n/l) = 0`, `l, n_i` linear in the 4 active variables, `J` quadratic. With multiplication matrices in an `R`-basis:

```
D = det M_l,   N_i = lambda( M_{n_i} adj(M_l) e ),   lambda: E -> R any R-linear functional,
N(G) = D lambda(J),   D div N - N(D) = 0,   D and N_i of active degree m.
```

Certification trick: localize the passive ring at polynomials nonzero at one specified fiber; an exactly nonzero specialized `10x10` minor is a unit there. **[Limit]** "A zero minor at a different fiber means that this chart fails, not that the target has no certificate." Coprimality of `D` with `G` and with the physical `F` after substitution are separate obligations.

Instances credited: `quadric-contact-broken` and `holdout-q2` (`HeadOfNaps`). In `holdout-q2`: monic degree-27 line algebra, six coefficient functionals + four gradients span all ten quadratic sources.

## 2.5 From sextics to cubic surfaces

### Coordinatewise-even sextics, all targets, all q **[Thm]**

Hypothesis: `F` even in each coordinate separately (`F = f(Z_0^2,...,Z_5^2)`). For a monomial `A`, pick the odd-exponent coordinates `x` (from X) and `y` (from Y); remaining `(a,b)` in X, `(c,d)` in Y. Map (reusable):

```
alpha = b/a,  beta = d/c,
Psi = (x^2, a^2, alpha, y^2, c^2, beta),   det(dPsi) = 16 x y   [Verified: each triple (x,a,b)->(x^2,a^2,b/a) has Jacobian 4x]
G(u,s,w,t) = f(u, s, alpha^2 s, w, t, beta^2 t)   cubic over K = k(alpha,beta), smooth (specialize alpha=beta=0).
```

For `A = x^r a^m b^n y^j c^h d^l` (`r,j` odd, `m+n`, `h+l` even):

```
B = alpha^n beta^l u^{(r-1)/2} s^{(m+n)/2} w^{(j-1)/2} t^{(h+l)/2} / 16   (times permutation sign), active degree 3q-1,
v = adj(dPsi) w_q(Psi)  with w_q from 2.4 (zero passive components appended),
q v(F) = A,  div v = 0,  physical component degree 6q-5.
```

Parity: `Psi(SZ) = Psi(Z)` with `S = diag(1,1,1,-1,-1,-1)` gives `adj(dPsi)(SZ) = -S adj(dPsi)(Z)`, so `v(SZ) = -S v(Z)` — exactly the required parities. Legality: the square map restricts `F=0` dominantly onto `G=0`, so a denominator coprime to `G` cannot pull back to a multiple of `F`; passive poles and chart boundaries cannot contain `F`. Constant block-linear coordinate changes preserve everything.

**Two-reflection extension [Thm]**: `F` globally Y-even and even separately in one X coordinate `x` and one Y coordinate `y` serves every allowed `A` odd in `x` and `y`. Smoothness of the generic cubic via `F_b = G_alpha / a`, `F_a = 2a G_s - alpha G_alpha / a`, `F_x = 2x G_u`, `F_y = 2y G_w`. **[Limit]** partial-target theorem; "global Y-evenness alone does not supply the selected reflections."

Denominator: active degree 45 over `K` but **[Limit]** "not a numerical degree bound in the six original coordinates." One finite physical support serves every `q` for fixed `F` and pair (fix the line basis and decompositions of the 56 degree-five active monomials into the partials); nine coordinate pairs still finite support; no bound on physical degree or pole multiplicity. "It gives no practical bound on the size of that calculation."

**Passive multipliers [Thm]**: if a strict `v` kills every passive ratio `r_j`, then for rational `h(r)`, `q(hv)(F) = hA`, `div(hv)=0`. **[Limit]** preserves the density sector; "It does not give arbitrary targets merely because they lie on the same sextic."

### Binary-cubic Hessian construction (`one-reflection`, `fqlx`; hidden reflection due to `fringewidth`)

Physical order `(x,a,c,y,b,d) = (X0,X1,X2,Y0,Y1,Y2)`,

```
F_R = x^6+a^6+c^6+y^6+b^6+d^6 + x^3 a y b + a^3 c y b        (= instance one-reflection)
u = x^2, p = y b, v = d^2, z = (y^3-b^3)/d, r = a/x, s = c/x
F_R = f(u,p) + v^3 + v z^2,   f = alpha u^3 + beta u^2 p + 2 p^3,  alpha = 1 + r^6 + s^6,  beta = r + r^3 s   [Verified]
Jacobian (that order) = 12 (y^3+b^3)/x                                                                        [Verified]
```

So density `u^2/24` serves `x^3 (y^3+b^3)/2`; the symmetry `(a,b) -> (-a,-b)` gives the other half; sum is `x^3 y^3`.

Binary-cubic Hessian data (reusable, all **[Verified]**):

```
delta = 324 alpha^2 + 24 beta^3
h = -beta^2 u^2 + 18 alpha u p + 6 beta p^2
H = (18 alpha u + 12 beta p) d_u + (2 beta^2 u - 18 alpha p) d_p,     k = H(f)/3,   div H = 0
H(h) = 0,   H(k) = 3 delta f,   delta f^2 - k^2 = 4 h^3 / 9,
delta u^2 = 108 alpha f_u + 12 beta^2 f_p - 12 beta h .
```

Over `K(sqrt delta)`: `f_± = (f ± k/sqrt delta)/2` are cubes of linear forms (`f_+ f_- = h^3/(9 delta)` **[Verified]**). With `lambda^2 = 3`: `v^3 + v z^2 = (1/2)(v + lambda z/3)^3 + (1/2)(v - lambda z/3)^3` **[Verified]**; `g = (v + lambda z/3)^3 / 2`. The Euler field `E_+ + E_g` on the cube-root coordinates has divergence 2, weight 1 on `h`, weight 3 on `f_+ + g`, so `h (E_+ + E_g) / (3 (f_+ + g))` is divergence-free and sends the full cubic to `h`. Average with the `sqrt delta`-conjugate **including the denominator**; combine with `d_u, d_p` through the `delta u^2` identity to serve `u^2/24`.

Denominator **[Stated]**: two half-denominators of degree 30, product degree 60. Both half-fields annihilate `a/x, c/x`; numerators share an `x^3` factor, so replacing it by any cubic `P_3(x,a,c)` serves the family `P_3 y^3` with the same denominator.

### Endpoint cubics and pole order two (`holdout-q2`, `HeadOfNaps`)

```
F_H = x^6+a^6+c^6+y^6+b^6+d^6 + x^3 a y b + a^3 c b d       (= instances holdout*, same variable order)
kappa^2 + kappa + 1 = 0,  alpha = y/x, beta = b/a (passive),
U = x^2, V = a^2, S = c^2 + kappa^2 d^2, T = kappa c d (active),
G = S^3 - 3 S T^2 + (1+alpha^6) U^3 + (1+beta^6) V^3 + alpha beta U^2 V + beta kappa^2 V^2 T,   G = F_H   [Verified]
```

Coefficient fiber `(1+alpha^6, 1+beta^6, alpha beta, beta kappa^2) = (2,2,1,1)` certifies the inverses needed by the §2.4 finite-algebra solver **[Stated]**.

`q=2` step (reusable):

```
U^3 V^2 = b(G) with b polynomial of component degree 3,  h = div b,  solve w_1(G) = h, div w_1 = 0,
w_2 = (b + G w_1)/2 - h E_4/6,   2 w_2(G) = U^3 V^2,  div w_2 = 0 .
```

Horizontal lift (all **[Verified]**: divergence 0, kill `alpha, beta`, act on `(U,V,S,T)` as `2 Delta_kappa` times the coordinate derivatives), `Delta_kappa = c^2 - kappa^2 d^2`:

```
C_U = (Delta/x) d_x + (Delta y/x^2) d_y,   C_V = (Delta/a) d_a + (Delta b/a^2) d_b,
C_S = c d_c - d d_d,                       C_T = -2 kappa d d_c + 2 kappa^2 c d_d .
```

Lifting `(alpha^3 beta^2 / 2) w_2` gives a strict pole-two field for `P Delta_kappa`, `P = x^3 a^2 y^3 b^2`, and

```
P c^2 = (2+kappa)/3 · P Delta_kappa + (2+kappa^2)/3 · P Delta_{kappa^2}     [Verified]
```

gives `holdout-q2` (`A = X0^3 X1^2 X2^2 Y0^3 Y1^2`). Passive multiplication serves the 34 monomials `x^{6-i} y^i a^{4-j} b^j c^2` and `... d^2`, `0<=i<=6`, `0<=j<=4`, `i+j` odd. **[Limit]** "a list of targets, not a claim that their residue classes are linearly independent."

### The density boundary **[Thm, Limit]**

If a linear reflection `rho` fixes the map (`Phi rho = Phi`) then `det(dPhi)(rho Z) = -det(dPhi)(Z)`: every density target `det(dPhi) B(Phi)` is `rho`-odd, and `rho`-fixed passive multipliers keep it so. Endpoint example: `rho_kappa(c,d) = (kappa d, kappa^2 c)` fixes `U,V,S,T` and negates `Delta_kappa` **[Verified]**. The central `q=1` target `x^3 y^3` is `rho_kappa`-invariant, so it is outside the span of all endpoint density sectors (averaging over the group generated by the endpoint reflections kills the sectors and fixes `x^3y^3`). "It does not exclude a different map or support. The distinction between a missing density sector and a missing cubic solver is central to the frontier in Chapter 4."

## 2.6 The conic specialization and physical degree bounds **[Thm]**

Hypothesis: `Phi_0 = (Q1,Q2,Q3,Q4,r,s)`, `Q_i` quadrics, `r,s` linear, `F = H(Phi_0)` smooth, `H` of weights `(2,2,2,2,1,1)`. Then for every `q>=1` and density `B` of weight `6q-4`, `A = det(dPhi_0) B(Phi_0)` has a strict certificate `N_q / D^q` over an extension of degree `<= 27`, `deg D = 10`, `deg N_q = 16q-5`; for a fixed chosen line the same `D` serves all `B` and `q`.

Proof data: `Phi_0` finite and surjective (common zero would be singular on `F`), `H` quasismooth, the cubic surface `H(Q,0,0)=0` smooth. 27 lines form a finite étale scheme (normal bundle `O(-1)`), so a line exists over `K/k` of degree `<= 27`. In adapted quadratic coordinates `(x,y,U,V)` with the line `U=V=0`:

```
G = a x^2 + h x y + b y^2 + c x + d y + e,  coefficients in (U,V,r,s) of physical degrees 2,2,2,4,4,6,
Q = [[a, h/2, c/2],[h/2, b, d/2],[c/2, d/2, e]],   det Q != 0 (else singular point of H),
D = -8 det Q   (physical degree 10).
```

`q=1` solver (reusable; coefficient rows **[Verified]**, `det M = -8 det Q` **[Verified]**):

```
B_1 = kappa x + lambda y + mu  (spectator coefficients),
M = [[2a, h, -c],[h, 2b, -d],[c, d, -2e]],   M (p1, p2, rho)^T = (kappa, lambda, mu)^T  (solve by adjugate, det M = D),
w = (rho x + p1) d_x + (rho y + p2) d_y,   w(G) - G div w = B_1,   div w = 2 rho.
```

Cofactor transport gives physical `V` with `V(F) - F div V = A`, degree 1; `s = div V` has degree 0; Euler correction `v_1 = V - s E_6 / 6` has `v_1(F) = A`, `div v_1 = 0`; numerator degree 11; `D` retained since `div V = 2 rho det(dPhi)`. Legality: `D` independent of `x,y` while `H` has positive active degree. Norm denominator degree `<= 270`.

**Coefficient raising identity** (reusable, **[Verified]** by hand for value, divergence, and the numerator form): for a family `F_t = F + tT`, `q v_q(t)(F_t) = A` with `A` independent of `t`, `div v_q(t) = 0`, dot = `d/dt` at 0,

```
v_{q+1} = ( q T v_q - F v_q-dot ) / (q+1) - v_q(T) E_6 / [6(q+1)],
(q+1) v_{q+1}(F) = A T,   div v_{q+1} = 0,
N_{q+1} = [ q (T D + F D-dot) N_q - F D N_q-dot - (D/6) N_q(T) E_6 ] / (q+1)     (for v_q = N_q/D^q).
```

Coefficient differentiation raises the power of `D` by at most one and adds no new divisor. Derivatives exist at every smooth member (27-line scheme finite étale over the smooth locus; Grassmannian chart with nonzero minors at the initial line; no successive extensions). Weight bookkeeping: a weight-`6q-4` monomial is a product of `3q-2` weight-two units; one unit is the initial density and the rest group into `q-1` sextic directions `T`.

Parity projection `v -> (v(Z) - S v(SZ))/2`; denominator divides `D(Z)^q D(SZ)^q`.

| Conditions | degree over line field | after trace |
|---|---:|---:|
| `Phi_0` equivariant for involutions + required parity | `<= 10q` | `<= 270q` |
| general `Phi_0`, then parity projection | `<= 20q` | `<= 540q` |

Coordinatewise-even `F` at `q=1`: square four coordinates, keep two linear; `det = const · C_I` (product of the four squared coordinates); density weight 2 gives `C_I (sum b_j Z_j^2 + b_rs r s)`; 6 carriers x 6 monomials + the all-six-odd monomial = 37 of the 226 allowed Y-odd sextic monomials with a degree-10 denominator. **[Limit]** the six-line construction serves the other 189 two-odd monomials "but has no stated physical degree bound." Arbitrary quadrics go beyond coordinatewise-even equations "while retaining its determinant-density target restriction."

## 2.7 Creating a reflection on a finite cover (`holdout` family; `fqlx` PR 128)

Chain sextic (instances `holdout`, `holdout-*`, `holdout-q2`, `holdout-bundle`):

```
F_H = sum_i (X_i^6 + Y_i^6) + X0^3 X1 Y0 Y1 + X1^3 X2 Y1 Y2 .
```

Smooth: at a nonzero critical point the max-modulus coordinate contributes a pure term of modulus `6M^5` versus at most `4M^5` from mixed terms. **[Limit]** "The endpoint cubic construction of §2.5 does not reach the reflection-invariant target X0^3 Y0^3."

Cover coordinates (reusable):

```
T^6 = Y0^6 + Y2^6,  alpha = Y0/T,  beta = Y2/T,  y = Y1,   K_b = Frac( C[alpha,beta]/(alpha^6+beta^6-1) ),
S = y^2 + T^2,  P = y T   (symmetric in y,T: the created reflection),
H = sum X_i^6 + P (alpha X0^3 X1 + beta X1^3 X2) + S^3 - 3 S P^2,   weights (1,1,1,2,2),   H = F_H   [Verified]
```

Residue divergence with normal coordinate `X0`, `g = H_{X0}`:

```
E_i = d_i - (H_i/g) d_{X0},   div_g(w) = sum_i ( E_i w_i - w_i E_i g / g ),
r_B = d_{X0} B / g - B H_{X0X0} / g^2 .
```

**Moving-conic primitive [Thm]**: for every polynomial `B in K_b[X0,X1,X2,S,P]` there is `w in K_H^4` with `div_g(w) = -r_B`. Mechanism: endpoint on `X1=0`; polynomial curve `phi(t)` with `X1 = kt`, other `X` quadratic in `t`, `S,P` quartic; `H(phi(t)) = 0` gives 12 positive-degree coefficient equations in 12 unknowns; exact complex seed on an actual fiber of the Fermat coefficient curve with invertible Jacobian; localize at that determinant (étale, flat) to reach the generic endpoint; evaluation embeds `K_H` into `E(t)`. Same matrix builds a polynomial normal `psi(t)` with `grad H(phi) · psi = 1`; with the curve tangent and three endpoint tangents, the 5x5 matrix has nonzero polynomial determinant `J(t)`. Pull back `r_B` by the residue Piola identity (polynomial in `t`), integrate termwise `sum c_n t^n -> sum c_n t^{n+1}/(n+1)`, divide by `J`, normalized trace back to `K_H` (trace is required, not an arbitrary functional: derivations commute with trace). **[Limit]** the primitive exists for every polynomial source on this quotient, but which physical targets it serves is a separate identification.

**Changing the ternary part** (`holdout-ternary-deformation`, `hodgepodgeconjecture` PR 151): replacing `sum X_i^6` by any ternary sextic `f(X)` keeps the exchange, the Fermat base and the target-selection identities, but the conic must be re-solved. Example `f = sum X_i^6 + (1/14)(X0^5 X1 + 2 X1^5 X2 + X2^5 X0)`; at `(alpha,beta) = (1,0)` in characteristic 5 the curve `X0 = 1+2t^2, X1 = t, X2 = 4+2t+2t^2, S = 4+t+t^2+3t^3+3t^4, P = 4+t+4t^2+4t^3+3t^4` lies on `14G = 0` **[Evidence]**; 13 unknowns (3 coefficients of `X0`, 2 of `X2`, 4 each of `S,P`) vs 13 coefficients of the degree-12 substituted equation, invertible Jacobian; a normal for `14G` must be multiplied by 14. **[Limit]** a perturbation that destroys the exchange or makes the base depend on an active coordinate is a different failure that a new conic cannot repair.

**Target selection by the deck character** (reusable; all identities below **[Verified]**):

```
J_c = 2 alpha^5 (T^2 - y^2) / T        (physical-to-quotient Jacobian, base coordinate beta)
Sufficient source pair:  J_c B = A - C t,   J_c d_{X0} B = d_{X0} A - (d_{X0} C) t,   Tr(t) = 0 over the physical field.
Central target:  B = X0^3 S / (2 alpha^2),   J_c B = X0^3 Y0^3 - X0^3 Y0^3 Y1^4 T^{-4}.
Deck action: physical coordinates fixed, (T, alpha, beta) -> (zeta T, zeta^{-1} alpha, zeta^{-1} beta), zeta^6 = 1.
Monomial criterion: A = M(X) Y0^a Y1^b Y2^c,  m = a+c+1,  e = b-m,
   B = M(X) alpha^{a-5} beta^c / 2 · (y^m T^b - T^m y^b)/(y^2 - T^2)   (polynomial in S,P over K_b),
   J_c B = A - C T^e,  C = M(X) Y0^a Y2^c Y1^m;  trace kills the correction iff e != 0 mod 6.
```

Useful exponents for degree-six targets: `e in {-4,-2,2,4}`; this covers `holdout` and the nine other single-monomial `q=1` holdout entries (per-target identities recorded by `DiamondJM`, discussion 44). **[Limit]** `e = 0`: the numerator of `B` vanishes; `e = -6`: the correction `T^{-6} = (Y0^6+Y2^6)^{-1}` is physical, deck-invariant, not killed by trace. "These are limitations of the selected source map, not nonexistence results."

**Equal coupling coefficients** (`one-reflection` sextic `F_eq = sum(X_i^6+Y_i^6) + Y0 Y1 (X0^3 X1 + X1^3 X2)`): same cover, `H_eq = sum X_i^6 + alpha P (X0^3 X1 + X1^3 X2) + S^3 - 3SP^2`; the base is still `alpha^6 + beta^6 = 1` (do **not** impose `2 alpha^6 = 1`). New conic seed needed (`DiamondJM`, PR 144): equality of coupling coefficients in `fqlx`'s conic family reduces to the octic

```
g(r) = r^8 + 7r^7 - 58r^6 + 149r^5 - 178r^4 + 115r^3 - 44r^2 + 10r - 1 = 0,  g(-13) = 35869443, g(-12) = -35000185,
```

root in `(-13,-12)` where seed denominators and both Jacobian factors are nonzero **[Evidence]**; `beta_c^6 = 1 - alpha_c^6` places it on the Fermat base. Transverse target (`one-reflection-y2-sector`): `A = X0^3 Y0 Y1 Y2`, `B = beta X0^3 P / (2 alpha^4)`, `J_c B = A - X0^3 Y0 Y1^3 Y2 T^{-2}` **[Verified]**. **[Limit]** For `F_eq` geometric existence for every odd target at every `q` follows "more directly" from §3.3 (in `U = Y0+Y1`, `V = Y0-Y1` the equation is separately even in `U,V,Y2`); the checked construction and the geometric theorem have different scopes.

**Returning to a homogeneous certificate** (reusable; `fqlx` PR 128): residue primitive in `K_F` with representatives `w_i = n_i/d`, `d != 0 mod F`, `g_i = d_i F`:

```
D = g_0 d,   N_0 = A d - sum_{i>=1} g_i n_i,   N_i = g_0 n_i (i>=1),
N(F) = A D exactly,   D div N - N(D) = F Q   =>  cleared defect F^2 Q (q=1), numerators of degree deg D + 1.
```

**[Limit]** "a finite algebraic existence proof ... It does not provide an expanded denominator or a numerical degree bound."

## 2.8 A complementary quotient (`holdout-bundle`, `fqlx` PR 161)

Target missed by 2.7: `A_* = X0^3 Y0 Y1^2` (deck exponent `e = 0`). Order-6 diagonal `g` fixing `X` and acting on `(Y0,Y1,Y2)` with weights `(5,1,5)` mod 6; equation and residue invariant (numerator character `zeta`, volume `zeta^5`). Quotient coordinates (reusable; **[Verified]**):

```
p = Y0 Y1,  v = Y1 Y2,  lambda = Y1^3 / (Y0^3 + i Y2^3),   k = C(lambda),
K = sum X_j^6 + X0^3 X1 p + X1^3 X2 v + a p^3 + b v^3,   a = lambda + 1/lambda,  b = i (lambda - 1/lambda),
theta = (X0,X1,X2,p,v,lambda):  K(theta) = F_H,  J_theta = -6 lambda Y1,  rho = -X0^3 p/(6 lambda),  J_theta rho = A_*,
Y1^6 = lambda (p^3 + i v^3).
```

Curve/normal/integration: endpoint `(c0,0,c2,p0,v0)` with `c0^6 + c2^6 + a p0^3 + b v0^3 = 0`; curve `Phi_1 = t`, `Phi_0, Phi_2` quadratic, `Phi_p, Phi_v` quartic; 12 unknowns, equations `[t^j] K(Phi) = 0`, `1<=j<=12`. Seed over `F_5` **[Evidence]**: `(a,b) = (0,3)`, `(c0,c2,p0,v0) = (1,1,4,1)`, `u = (1,0,4,1,4,2,2,4,4,1,1,2)` (order: 2 coeffs of `Phi_0`, 2 of `Phi_2`, 4 of `Phi_p`, 4 of `Phi_v`), `det Delta = 3`. Integral endpoint relation at `(a,b) = (0,-2)`, realized by `lambda = i`; §5.4 lifting gives a nonsingular complex point on the actual coefficient locus. No inversion of `a` or `lambda^2+1` at the seed. Normal vector starts as `e_0/(6 c0^5)` and is corrected with `Delta^{-1}`; evaluation determinant has constant term `1/(6 c0^5)`. Injectivity of `K_quot -> E(t)` from independent endpoint coordinates, `Phi_1 = t`, irreducibility. Physical transport through `theta` includes derivatives of `lambda` (not a constant on `F`). Bundle `X0^3 Y0^3 + X0^3 Y0 Y1^2` by the Chapter 1 addition formula; quotient remainders retained.

Source boundaries for `A = M(X) Y0^a Y1^b Y2^c`, `e = b - a - c - 1`:

| `e` | `A / J_theta` in quotient coordinates |
|---|---|
| `0` | `-M(X) p^a v^c / (6 lambda)` (polynomial; same curve/integration serves it) |
| `-6` | `-M(X) p^a v^c / (6 lambda^2 (p^3 + i v^3))` (fiber denominator; needs a different source representative or primitive; rational integration can meet nonzero residues) |

**[Limit]** "The displayed source maps need not exhaust their closure under automorphisms"; for `e=-6` geometric existence still follows from §4.6. Closing sentence of 2.8: "None establishes that a suitable map or supporting family exists for a general Y-even sextic."

## 2.9 Moving cubic maps and Jacobian residue transport

### Change of normal density (reusable) **[Thm]**

```
ResDiv_g(w) = sum_i d_i w_i - (sum_i w_i d_i g)/g = g sum_i d_i (w_i/g),
ResDiv_g( (g/h) det(M) N w ) = (g/h) det(M) ResDiv_h(w),   M = coordinate Jacobian, N = M^{-1}, h = new normal derivative.
```

Balanced-numerator normal-chart identity (Lean `ResidueTransport.NormalChange`): `f = F_p`, `g = F_q`, `A = B p q` with `B` independent of `p,q`,

```
G = p f + q g,   H = p^2 F_pp + 2 p q F_pq + q^2 F_qq,   D = d_q - (g/f) d_p,   C = -B p q^2 / G,
A_p/f - A F_pp/f^2 = A (3G - H)/G^2 + D C - C (D f)/f      (f, G != 0).
```

### First variation of a Weierstrass map (reusable) **[Thm]**

`eta^2 = xi^3 + a xi + b`, `Da = Db = 0`, `D xi/(2 eta) = m/g` (`m,g,eta != 0`), dot = compatible first variation fixing the tangent coordinate and `m`:

```
2 eta eta-dot = (3 xi^2 + a) xi-dot + a-dot xi + b-dot,
D xi-dot/(2 eta) - D xi eta-dot/(2 eta^2) = -m g-dot/g^2,
g-dot/g = (a-dot xi + b-dot)/(2 eta^2) - (g/m) D( xi-dot/(2 eta) ).     <-- exact-differential correction; omitting it changes the source
```

Established at a single coefficient value with dual numbers `K[eps]/(eps^2)` (Lean `ResidueTransport.MovingMap`; avoids assuming `d/d(coefficient)` is a derivation of the specialized function field).

### Cyclic paired-sign cubic (`cycle-three`, `hodgepodgeconjecture` PR 165)

```
F_C = sum_i (X_i^6 + Y_i^6) + (1/2) sum_{i mod 3} X_i^3 X_{i+1} Y_i Y_{i+1},   q=1, A = X0^3 Y0^3.
r = X2^2,  v = Y0/X0,  c_i = (Y_i/X_i)/v,  x = X0^2/r,  y = X1^2/r   (chart z = 1),
U = a x^3 + b y^3 + c z^3 + d x^2 y + e y^2 z + f z^2 x,
(a,b,c) = (1+v^6, 1+c_1^6 v^6, 1+c_2^6 v^6),   (d,e,f) = (1/2) v^2 (c_1, c_1 c_2, c_2).
s_4 = a e^2 f + b d f^2 + c d^2 e,
s_6 = 27 a^2 b^2 c^2 + 4 a^2 c e^3 + 4 a b^2 f^3 + 6 a b c d e f + 4 b c^2 d^3 - d^2 e^2 f^2,
syzygy on U=0 (covariants H, Theta, J of the construction):  J^2 = 4 Theta^3 + 4 s_4 Theta H^4 - s_6 H^6,
xi = Theta/H^2,  eta = -J/(2H^3):  eta^2 = xi^3 + s_4 xi - s_6/4,   D_y xi/(2 eta) = 3/U_x.
T = r v^6,  X = r^2 xi,  Y = r^3 eta:   V = Y^2 - X^3 - (1/8) a_8 X - (1/256) b_256 = 0,  a_8 = 8 r^4 s_4,  b_256 = -64 r^6 s_6
   (degrees 4 and 6 in (r,T); consistent with eta^2 relation — checked by hand).
```

Polynomial curve on `V=0`: degrees `(2,4,6)` in `(r,X,Y)`, parameter `T`; 13 coefficient equations; seed over `F_5` at `(c_1,c_2) = (2,3)` **[Evidence]**: `r = 1+2T+4T^2`, `X = 3+T+3T^2+4T^3+4T^4`, `Y = 2+3T+3T^2+4T^3+4T^4+T^5+T^6`; endpoint parameters `r(0), X(0)` retained, other 13 coefficients unknown. Physical recovery: vary `a = alpha + v^6` at `alpha = 1`, source `P_J = -V_alpha/(144 r)` (polynomial: `V_alpha` has the factor `r`), residue-volume multiplier `F_{X0} det M / (2Y) = -144 r v^3`; pull back the primitive, add the moving-map correction, subtract the normal-chart correction. **[Limit]** "no bound on deg D, and no expanded witness, is asserted."

### Transitive paired-sign cubics (`mixed-biweight`, `fqlx` PR 170; `holdout-moving-branch`)

```
U = a x^3 + b y^3 + c z^3 + d x^2 y + e y^2 z + f x^2 z     (transitive triangle of mixed monomials)
s_4 = -3 b c d f + c d^2 e - (1/3) e^2 f^2,
s_6 = 27 a^2 b^2 c^2 + 4 a^2 c e^3 + 4 b^2 c f^3 + 4 b c^2 d^3 - 4 b c d e f^2 + (8/3) c d^2 e^2 f + (8/27) e^3 f^3,
same syzygy / invariant-differential form; Lean `TransitiveCubic`.
mixed-biweight:   (a,b,c) as above,  d = (1/2) c_1 v^2,  e = (1/2) c_1 c_2 v^2,  f = (1/3) c_2 v^4;
   V = Y^2 - X^3 - A_4 X - B_6,  A_4 = r^4 s_4,  B_6 = -r^6 s_6/4;  seed at c_1=c_2=1 over F_5: r = 1+T^2, X = 3+4T^2+3T^4, Y = 2+3T^6; 13x13 Jacobian invertible  [Evidence]
holdout-moving-branch:  d = c_1 v^2,  e = c_1 c_2 v^2,  f = (1/7) c_2^2 v^4;  integral Jacobian equation scale 111132;
   seed at c_1=c_2=1 over F_5: r = 1+2T^2, X = 4+4T+T^2+3T^3+T^4, Y = 1+3T^2+4T^4+2T^6; fix r(0)=1, X(0)=4  [Evidence]
```

Both use `a = alpha + v^6`, source `-V_alpha/(144 r)`, multiplier `-144 r v^3`, target `X0^3 Y0^3`.

### Reducing equations and source variations (shear; `fqlx`, discussion 58; Lean `TransitiveCubic.Shear`)

```
(U + eps Q)(x + eps v(x)) = U + eps (Q + sum_i v_i U_{x_i}),  v(x) = M x, tr M = 0, det(1 + eps M) = 1   (eps^2 = 0):
   variations may be reduced modulo Jacobian directions of trace-zero linear fields.
Scalar correction:  delta((h - c g)/g) = delta(h/g) - delta(c);  if delta c = 0 the complete normal expression is preserved.
Example U = a x^3 + b y^3 + c z^3 + d x^2 y + e y^2 z + f x y^2,  a != 0,  d^2 - 3 a f != 0:
   choose k with 3 a k^2 + 2 d k + f = 0 (finite extension), shear x = X + k y (det 1):
   U = a X^3 + b' y^3 + c z^3 + d' X^2 y + e y^2 z,  b' = b + a k^3 + d k^2 + f k,  d' = d + 3 a k,  (d')^2 = d^2 - 3 a f.
   Source x y z -> (X + k y) y z;  infinitesimal shear X -> X - eps z/(2 d') gives
   (X + k y) y z - z U_X/(2 d') = k y^2 z - (3a/(2 d')) X^2 z   (two coefficient directions of the transitive family).
   Normal derivative d_X on z=1 (before restriction to U=0) fixes the correction coefficient 1/(2 d').
   For the xyz variation:  I_4-dot = (9 a b c - c d f)/2,   I_6-dot = a c d e^2,   with I_4 = s_4, I_6 = -s_6/4.
```

Derivations fixing `a,d,f` fix `k`; parameter derivations may move `k` and must enter the physical coordinate matrix. Covariant denominators proved nonzero by the §5.4 specialization criterion. **[Limit]** "These steps change the presentation of a source; they do not enlarge the supported representations of the quotient."

### Coefficient and source scope **[Thm]**

Symmetry hypotheses: global Y-evenness + the three paired signs `(X_i,Y_i) -> (-X_i,-Y_i)` + diagonal action with weights `w = (2,1,0,3,2,1)` mod 6. Equation space: 26 monomials (6 pure sixth powers, 2 balanced `u0u1u2` monomials with `1, v^6`, 9 per edge orientation). Paired chart (`(u_0,u_1,u_2)` weights `(4,2,0)`, sum 0 mod 6):

```
U = sum_i (a_i + b_i v^6) u_i^3 + (g_0 + g_1 v^6) u_0 u_1 u_2
    + v^2 (d_0 u_0^2 u_1 + d_1 u_1^2 u_2 + d_2 u_2^2 u_0) + v^4 (e_0 u_0 u_1^2 + e_1 u_1 u_2^2 + e_2 u_2 u_0^2).
```

Classical `SL_3` invariants `I_4 = s_4`, `I_6 = -s_6/4` (coefficient degrees 4, 6; Fisher): `I_d(U(v))` has `v`-degree `<= 6d` in powers of `v^6`, so `r^d I_d(U(v))|_{v^6 = T/r}` is homogeneous of degree `d` in `(r,T)` — the Jacobian model `Y^2 = X^3 + A_4 X + B_6` exists on the whole 26-dimensional coefficient space wherever the generic cubic is smooth; no differential-character twist (determinant-one action). The genus-one surface can have multiple fibers; its Jacobian need not be birational to it.

Polynomial sources: 22 odd sextic monomials are paired-sign-invariant with diagonal-invariant residue (6 of Y-degree 1, 10 of degree 3, 6 of degree 5). Writing such a target divided by `r^3` as `c v^j M(u)`, `j in {1,3,5}`, `k_M` the coefficient of `M` in `U`:

```
P_{d,M,j} = r^{d-1} v^{j-3} (dI_d/dk_M)(U(v)),  d = 4,6,  polynomial in (r,T);  Jacobian source = (c/144)(P_{4,M,j} X + P_{6,M,j}).
```

Character formula: a monomial residue with exponent vector `e` has character `zeta^{<w,e> + 9}` under the diagonal action (9 = determinant character of the projective volume), *independently of the pole order*; averaging to the invariant quotient retains exactly character one. Hodge numbers of the invariant primitive sector `(22,74,22)` versus `(208,866,208)` for the global involution **[Stated]**. **[Limit]** Chapter 4 proves the stronger geometric statement (smooth `(3,3)` surface => whole `Gamma`-invariant `H^4` has divisor support, on a nonempty coefficient open); extraction by the Jacobian method still needs the covariant differential, its first variation, a serving curve family with polynomial normal, and physical transport; "a particular seed need not survive every specialization"; "an elliptic correspondence alone does not integrate an arbitrary rational differential."

---

## Reusable explicit formulas (index)

1. Cofactor transport `V = adj(dPhi) w(Phi)`, `V(G(Phi)) = det(dPhi) w(G)(Phi)`, `div V = det(dPhi) (div w)(Phi)` (2.1). **[Thm]**
2. Spectral constructor `(mu,nu)`: `D0 = w(G)`, `w(D0) = (mu+nu)D0 - mu nu G`, `div(bw) = (mu+nu) b` => `E2` quotient `mu nu A` (2.1); nilpotent branch `w^2(G) = 0`, `div(bw) = 0` (2.1).
3. Subset Euler field `A E_I/(q E_I(F))` under `E_I^2(F) = 6 E_I(F)`, `deg_I A + |I| = 6` (2.2).
4. Nambu field `iota_W vol = dP ^ dh_1 ^ ... ^ dh_4`, `F = Gamma(h) + PK`, `A = W(K)`, field `W/P` (2.2).
5. Perfect-power contact `Q = cH^m + DK`, `W = grad D x grad H`, `v = B(U,H) W/D`, `A = B W(K)` (2.3).
6. Line cofactor field `w_L = adj(d(l1,Q1,l2,Q2)) e_2 / l1`, `J_L = det`, and the `q`-recursion `v_q = [b + Gw - cE_4/(3(q-1))]/q` (2.4). **[Verified by hand]**
7. Finite-algebra line realization `D = det M_l`, `N_i = lambda(M_{n_i} adj(M_l) e)`; trace `tr(adj(M_D) M_N)/([E:K] det M_D)` (2.4).
8. Coordinatewise-even map `Psi = (x^2,a^2,b/a,y^2,c^2,d/c)`, `det = 16xy`, density `B = alpha^n beta^l u^{(r-1)/2} s^{(m+n)/2} w^{(j-1)/2} t^{(h+l)/2}/16`, `v = adj(dPsi) w_q(Psi)` (2.5). **[Verified: det]**
9. Binary-cubic Hessian package: `delta, h, H, k`; `H(h)=0`, `H(k)=3 delta f`, `delta f^2 - k^2 = 4h^3/9`, `delta u^2 = 108 alpha f_u + 12 beta^2 f_p - 12 beta h`; `f_±`; `v^3 + vz^2` as two cubes with `lambda^2 = 3` (2.5). **[Verified]**
10. `F_R` rational coordinates `(u,p,v,z,r,s)`, Jacobian `12(y^3+b^3)/x` (2.5). **[Verified]**
11. Endpoint cubic `(U,V,S,T)` with `kappa`, `G` formula, `w_2 = (b + G w_1)/2 - h E_4/6`, horizontal fields `C_U, C_V, C_S, C_T`, `Pc^2` decomposition into `Delta_kappa, Delta_{kappa^2}` (2.5). **[Verified; w_2 is the q=2 case of item 6: (b+Gw_1)/2 - hE_4/6, checked by hand]**
12. Density-boundary reflection law `det(dPhi)(rho Z) = -det(dPhi)(Z)` (2.5). **[Thm]**
13. Conic solver: `Q`, `D = -8 det Q`, matrix `M`, field `w = (rho x + p1) d_x + (rho y + p2) d_y`, `w(G) - G div w = B_1`, Euler correction `v_1 = V - sE_6/6` (2.6). **[Verified: rows and det]**
14. Coefficient-raising identity `v_{q+1} = (qTv_q - F v_q-dot)/(q+1) - v_q(T)E_6/[6(q+1)]` and its numerator form (2.6). **[Verified by hand]**
15. Parity projection `(v(Z) - S v(SZ))/2` with denominator `D(Z)^q D(SZ)^q` (2.6).
16. Fermat cover `T^6 = Y0^6 + Y2^6`, `alpha, beta, S = y^2+T^2, P = yT`, `H`, `J_c = 2 alpha^5 (T^2-y^2)/T`, central source `B = X0^3 S/(2 alpha^2)`, monomial criterion `B = M alpha^{a-5} beta^c (y^m T^b - T^m y^b)/(2(y^2-T^2))`, `J_c B = A - C T^e` (2.7). **[Verified]**
17. Residue divergence `div_g`, source `r_B = d_{X0}B/g - B H_{X0X0}/g^2`, moving-conic ansatz (12x12), normal `psi`, determinant `J(t)`, termwise integration, trace (2.7).
18. Deformed ternary `f` with `(1/14)(X0^5X1 + 2X1^5X2 + X2^5X0)` and its char-5 seed curve (2.7). **[Evidence]**
19. Equal-coupling octic `g(r)` and transverse source `B = beta X0^3 P/(2 alpha^4)` (2.7). **[Verified: source identity]**
20. Homogeneous extraction `D = g_0 d`, `N_0 = Ad - sum_{i>=1} g_i n_i`, `N_i = g_0 n_i` (2.7).
21. Complementary quotient `p = Y0Y1`, `v = Y1Y2`, `lambda = Y1^3/(Y0^3 + iY2^3)`, `K`, `a = lambda + 1/lambda`, `b = i(lambda - 1/lambda)`, `J_theta = -6 lambda Y1`, `rho = -X0^3 p/(6 lambda)`, `Y1^6 = lambda(p^3 + iv^3)`, the `e = 0 / -6` table, and the `F_5` seed (2.8). **[Verified except seed]**
22. `ResDiv_g` Piola change of normal density; balanced-numerator correction `C = -Bpq^2/G` (2.9).
23. Weierstrass first-variation formula (boxed) and the dual-number method (2.9).
24. Cyclic cubic data: `(a,b,c,d,e,f)`, `s_4, s_6`, syzygy, `xi = Theta/H^2`, `eta = -J/(2H^3)`, `D_y xi/(2 eta) = 3/U_x`, `T = rv^6`, `X = r^2 xi`, `Y = r^3 eta`, `a_8 = 8r^4 s_4`, `b_256 = -64 r^6 s_6`, source `-V_alpha/(144 r)`, multiplier `-144 r v^3`, `F_5` seed (2.9). **[Stated; a_8/b_256 scaling consistent with the eta^2 relation]**
25. Transitive cubic data: `s_4, s_6`, the two coefficient specializations (`mixed-biweight`: `d,e,f = c_1v^2/2, c_1c_2v^2/2, c_2v^4/3`; `holdout-moving-branch`: `c_1v^2, c_1c_2v^2, c_2^2v^4/7`, scale 111132) and seeds (2.9). **[Stated/Evidence]**
26. Shear reduction `3ak^2 + 2dk + f = 0`, `x = X + ky`, `b', d'`, `(d')^2 = d^2 - 3af`, infinitesimal shear `X -> X - eps z/(2d')`, `I_4-dot = (9abc - cdf)/2`, `I_6-dot = acde^2` (2.9). **[Stated]**
27. General paired-sign/diagonal family `U` (26 monomials), `P_{d,M,j} = r^{d-1} v^{j-3} dI_d/dk_M`, source `(c/144)(P_4 X + P_6)`, character `zeta^{<w,e>+9}` with `w = (2,1,0,3,2,1)` (2.9). **[Stated]**

Note on the task's named items: the "Hessian/adjugate E1 lift" in this chapter is item 1 (adjugate transport) plus item 9 (binary-cubic Hessian); the "Koszul" material in this chapter is the adjunction/residue-pairing argument of 2.2 (Koszul resolution of a zero-dimensional complete intersection) and the Hilbert series `(1+t)^4` of 2.4. A "Koszul correction map" as a computational device is in §5.2 of Chapter 5 ("The full correction map"), not in this chapter.

## Lean modules the chapter points to

`HodgeComputer.Constructions.HoldoutCover` (+ `GeneralizedCover`, `certificate_of_source_pair`), `HodgeComputer.CertificateOperations` (addition), `ResidueTransport.NormalChange`, `ResidueTransport.MovingMap` (jets and commuting-derivation forms), `TransitiveCubic`, `TransitiveCubic.Shear`, `HigherResidues` (third-order poles, split-ternary), `PolynomialCurves` (integral lifting, polynomial-family integration), `AlgebraicSpecialization` (§5.4 nonvanishing at algebraic roots), `AveragingProjectors` (§3.1).

## Instances the chapter attributes to each technique (all currently **solved**, i.e. `proof.json` present)

| Technique | Instance(s) named in the chapter | Credit |
|---|---|---|
| 2.4 finite-algebra line realization | `quadric-contact-broken`, `holdout-q2` | `HeadOfNaps` |
| 2.5 binary-cubic Hessian + hidden reflection | `one-reflection` (half-denominators 30+30, product 60; family `P_3 y^3`) | `fqlx`, after `fringewidth` |
| 2.5 endpoint cubic, q=2 | `holdout-q2` (+34 listed monomials) | `HeadOfNaps` |
| 2.7 Fermat cover + moving conic | `holdout` and the nine other single-monomial holdout entries (`holdout-p-y14`, `holdout-p3`, `holdout-x1-y04-y1`, `holdout-x12-p-y02`, `holdout-x13-y13`, `holdout-x0-y03-y12`, `holdout-all-six`, `holdout-x0-x14-y0`, `holdout-x04-p`; `e` values `-4,-2,2,4` — I checked each `e` against its record) | `fqlx` PR 128; per-target identities `DiamondJM` disc. 44 |
| 2.7 changed ternary part | `holdout-ternary-deformation` | `hodgepodgeconjecture` PR 151 |
| 2.7 equal-coupling seed | `one-reflection-y2-sector` | `DiamondJM` PR 144 (seed from `fqlx`'s family) |
| 2.8 complementary quotient + bundle | `holdout-bundle` (`A = X0^3Y0^3 + X0^3Y0Y1^2`) | `fqlx` PR 161 |
| 2.9 cyclic Jacobian transport | `cycle-three` | `hodgepodgeconjecture` PR 165 |
| 2.9 transitive covariants | `mixed-biweight`, `holdout-moving-branch` | `fqlx` PR 170; second specialization |
| 2.9 shear reduction | (discussion 58 only; no instance anchor in this chapter) | `fqlx` with Codex |
| 2.6 conic theorem | no instance anchor in this chapter | — |

---

## Where this chapter says the open frontier instances fail

**Fact:** the chapter names no *open* instance by ID. Every `#case/...` anchor in it (`quadric-contact-broken`, `holdout-q2`, `one-reflection`, `holdout`, `one-reflection-y2-sector`, `cycle-three`, `holdout-moving-branch`) is a solved instance. The chapter instead states boundaries of each technique in general terms and defers the frontier to Chapter 4 ("The distinction between a missing density sector and a missing cubic solver is central to the frontier in Chapter 4"). The verbatim boundary statements are:

1. (2.1) "These spectral hypotheses are sufficient conditions for explicit fields. They are not a characterization of certificates."
2. (2.2) "the contact selects a specific target through a residue pairing; it does not serve an arbitrary prescribed A."
3. (2.3) "A root H on the normalization of D need not descend to an ambient polynomial ... Enforcing an ambient polynomial root would exclude such contacts. Conversely, permitting rational roots without checking their extra poles would admit many false candidates."
4. (2.4) "A zero minor at a different fiber means that this chart fails, not that the target has no certificate."
5. (2.5) "global Y-evenness alone does not supply the selected reflections." / "This argument bounds neither its physical degree nor its pole multiplicities." / passive multipliers: "It does not give arbitrary targets merely because they lie on the same sextic."
6. (2.5, density boundary) "every density target det(dPhi)B(Phi) is odd under rho ... The central q=1 target x^3y^3 is invariant. Even combining these endpoint constructions does not reach it ... Thus this target is outside the linear span of those density sectors. It does not exclude a different map or support."
7. (2.6) "The six-line construction serves the other 189 two-odd monomials too, but has no stated physical degree bound. Arbitrary quadratics enlarge the bounded construction beyond coordinatewise-even equations while retaining its determinant-density target restriction."
8. (2.7) "The endpoint cubic construction of §2.5 does not reach the reflection-invariant target X0^3 Y0^3." (resolved there by the cover) / "When e=0, the numerator defining B is zero. When e=-6, the correction is deck-invariant ... Neither case is removed by this trace argument. These are limitations of the selected source map, not nonexistence results for certificates." / "If a perturbation destroys the exchange or makes the coefficient curve depend on an active coordinate, replacing the conic alone does not repair that different failure." / "It does not provide an expanded denominator or a numerical degree bound."
9. (2.8) "The e=-6 densities have a fiber denominator and require a different source representative or primitive construction: rational integration can encounter nonzero residues." / "None establishes that a suitable map or supporting family exists for a general Y-even sextic."
10. (2.9) "no bound on deg D, and no expanded witness, is asserted." / "These steps change the presentation of a source; they do not enlarge the supported representations of the quotient." / "A nonsingular incidence persists on an open, but a particular seed need not survive every specialization." / "Removing the extra diagonal symmetry, or the paired signs themselves, can also change the surface's geometric genus ... an elliptic correspondence alone does not integrate an arbitrary rational differential."

**[My inference]** Mapping these boundaries onto the currently open instances (from `data/catalogue.json`; open = no `proof.json`), using only the chapter's own criteria. These are my cross-references, not statements of the chapter:

| Open instance | `(q, A)` | Which chapter boundary applies |
|---|---|---|
| `cycle-q2` | `q=2`, `A = X0^2X1^2X2^3Y0^2Y1^2Y2` on `F_C` | Character under `w=(2,1,0,3,2,1)`: `<w,e>+9 = 17+9 = 26 = 2 mod 6`, so the residue has character `zeta^2`, not the invariant character the 2.9 Jacobian quotient retains ("independently of the pole order"). The §2.9 sector serves only character one; a different serving quotient (§3.1 averaging ideal) would be needed. Also the chapter's q-raising tool (2.6) needs a `t`-family with `A` fixed and is not applied to `F_C`. |
| `off-paired` | `q=1`, `A = X0^3Y0^3`; `F_C` with `X0X2^3Y0Y2/2` replaced by `X0^2X1X2Y0Y2/3` | The new term is odd under the paired sign `(X1,Y1)->(-X1,-Y1)` and has weight `9 != 0 mod 6` under the diagonal action, so it lies outside the 26-monomial paired-sign/diagonal space of 2.9; boundary 10 ("Removing ... the paired signs themselves ..."). |
| `one-reflection-skew` | `q=1`, `A = X0^3Y0^3`; `F_eq + (X0X1Y0Y1)(Y0^2-Y1^2)/7` | The extra term is antisymmetric under the exchange `Y0 <-> Y1` (`y <-> b`) that the 2.5 Hessian coordinates `(u,p,v,z,r,s)` and the 2.7 cover (`y,T` exchange) rely on; `(y^2-b^2)` is not a function of `p=yb`, `z=(y^3-b^3)/d`, `v=d^2`. Boundary 8 ("a perturbation [that] destroys the exchange"). |
| `quadric-contact-residual` | `q=1`, `A = X0^3Y0^3` on the `quadric-contact-broken` sextic | The 2.4 finite-algebra realization served the broken-reflection target `X0^3X1^2Y0 - X0X1^2Y0^3`; by boundary 6 every determinant-density target is odd under any reflection fixing the map, and `X0^3Y0^3` is the invariant ("residual") sector. |
| `quadratic-cubic-pullback` | `q=1`, `A = X0^3Y0^3`; `F = H(Q_1..Q_4, r, s)`-type | Boundary 7: the conic theorem reaches only determinant-density targets `det(dPhi_0)B(Phi_0)`; the instance title "residual density" indicates `X0^3Y0^3` is outside that sector. |
| `dense-rational`, `short-power-presentation`, `eight-plane-configuration` | `q=1`, `A = X0^3Y0^3`, full or near-full coefficient support | Boundary 9: no technique in this chapter applies to a general Y-even sextic; Chapter 1 adds that for full support all Y-odd targets of one degree lie in a single character coset (no fine character decomposition). |
| `cyclic-triple-residual`, `reflection-coupled-residual`, `klein-coupled-character`, `symmetric-sextic-residual`, `four-plus-two-residual`, `polarized-k3-remainder`, `rank-one-biquartic`, `diagonal-permutation` | various | Not addressed by any construction of this chapter (their solved companions `cyclic-triple-fano`, `reflection-coupled-fano`, `polarized-contact-*`, `split-ternary-q2` are Chapter 3/3a material). The only applicable general statement is boundary 9. |

What the chapter does **not** claim: it never asserts nonexistence of a certificate for any target; every boundary is phrased as a limit of a specific map, source, seed, or density sector ("not nonexistence results for certificates").

---

## Verification log (this session)

sympy 1.12, `python3 -I`, all exact:

- `det` of `(x^2, a^2, b/a)` w.r.t. `(x,a,b)` = `4x`.
- `F_R = f + v^3 + vz^2` under `(u,p,v,z,r,s)`; Jacobian `= 12(y^3+b^3)/x`; `H(h)=0`; `H(k)=3 delta f`; `delta f^2 - k^2 = 4h^3/9`; `delta u^2 = 108 alpha f_u + 12 beta^2 f_p - 12 beta h`; `div H = 0`; `v^3+vz^2` two-cube split; `f_+ f_- = h^3/(9 delta)`.
- `G = F_H` for the endpoint cubic (reduction modulo `kappa^2+kappa+1`); `C_U, C_V, C_S, C_T` have zero divergence, kill `alpha, beta`, and act as `2 Delta_kappa` times coordinate derivatives; `rho_kappa` negates `Delta_kappa` and fixes `S, T`; `P c^2` identity.
- `J_c = 2 alpha^5 (T^2 - y^2)/T`; `H = F_H` on the cover; `J_c B` for the central target; the monomial criterion identity for every odd Y-degree `<= 5` exponent triple with `e != 0` (observed `e in {-6,-4,-2,2,4}`); transverse-target identity on `F_eq`.
- `K(theta) = F_H`; `J_theta = -6 lambda Y1`; `J_theta rho = X0^3 Y0 Y1^2`; `Y1^6 = lambda(p^3 + i v^3)`.
- Conic solver: `w(G) - G div w` has vanishing quadratic part and linear/constant coefficients equal to the rows of `M`; `det M = -8 det Q`.
- Not machine-checked: the cubic-surface line-class spanning argument, the 45-plane denominator, the degree-30/60 denominators of `one-reflection`, all finite-field seeds, the ternary-cubic covariant syzygies and `s_4, s_6` formulas, the shear `I_4-dot, I_6-dot`, the Hodge numbers `(22,74,22)` / `(208,866,208)`.
