# Technique map: Chapter 3A "Split sextics, elliptic contraction and the K3 remainder"

Source: `/home/kejace/devel/sextics/docs/03a-split-and-k3.md` (349 lines, read completely).
Section anchors (line numbers in that file): 3A.1 l.7 · 3A.2 l.43 · 3A.3 l.72 · "Algebraic contraction with moving poles" l.106 · "Higher-pole residue calculus" l.135 · "Involution-stable ternary decompositions" l.185 · 3A.4 l.191 · "The algebraic split decomposition and its boundaries" l.215 · 3A.5 l.241 · "A polynomial test for the constructed part" l.279 · 3A.6 l.293 · "An unconditional order-three family" l.331.

Conventions used below (from Chapter 1 / `formal/README.md`): variables `Z=(X0,X1,X2,Y0,Y1,Y2)`, `tau(X,Y)=(X,-Y)`, target `A` is Y-odd of degree `6q`, certificate `v=N/D`, `deg D=d`, `deg N_i=d+6q-5`, `gcd(D,F)=1`, cleared identity
`A D^2 - q D sum_i N_i dF/dZ_i + F (D sum_i dN_i/dZ_i - sum_i N_i dD/dZ_i) = F^(q+1) Q`.
Normalized form after Euler shift: `E1: q N(F) = A D`, `E2: D div N - N(D) in (F^q)`.

Status labels: **[proved in book]** = stated with proof in the chapter; **[delegated]** = relies on a theorem of another chapter or an external citation; **[computational evidence]** = checked by me in this session with sympy; **[my cross-reference]** = linkage to repo instances that the chapter itself does not make.

---

## 0. External inputs the chapter relies on (so hypotheses can be tracked)

| Input | Where | Exact statement used |
|---|---|---|
| Quotient-transfer lemma | Ch.3 §3.1 (`docs/03-geometric-support.md` l.27) | `M` smooth projective with faithful finite `G`-action, `N` smooth projective model of `M/G`: `H^k(N)` divisor-supported iff `H^k(M)^G` divisor-supported. Proof includes all graph and exceptional boundaries. |
| Covering-curve support theorem | Ch.3 §3.2 | `W` smooth projective `n`-fold, actual morphism `e: P^1 x S -> W` dominant generically finite, `S` smooth integral affine of dim `n-1` ⇒ all `H^n(W)` supported on one divisor computable by finite algebraic operations (graph closure boundary images ∪ positive-dim-fibre locus ∪ non-étale locus, enlarged to a divisor `E`). |
| Bounded Fano input | Ch.3 §3.2 (Debarre Thm 7.5, 7.7; Miyaoka–Mori) | smooth Fano `n`-fold: covered by rational curves with `-K.C <= n+1`; normal projective `n`-fold with ample Q-Cartier `-K`: `-K.C <= 2n`. Search complete Hom schemes `P^1 -> W` up to that `L`-degree; if `W` singular, resolve and lift. |
| Reflection theorem | Ch.3 §3.3 (l.129) | `F` even in one coordinate `z` ⇒ every target odd in `z` has divisor support in its residue class (keep global Y-odd condition). |
| From support to every `q` | Ch.3 §3.6 | enlarge computed support to homogeneous `D` with `F ∤ D`, take union of Galois conjugates and of `tau`-image (product of equations, rational, Y-even); for `1<=q<=4` solve the complete linear system `q N(F)=A D^q`, `D^q div N - N(D^q) = F^q R`, `deg N_i = q(deg D+6)-5`; for `q>=5` write `A = sum B_i F_i` (socle degree 24) and induct `v=(B+F w)/q`, `K_q(A;v)=(F/q) K_{q-1}(div B; w)`. No strictness, no degree bound, no runtime bound. |
| Prescribed-support theorem | Ch.1 | support on `V(D)` ⇒ certificate with denominator `D_red^q` (hence `D^q`); at `q=1` any power of `D` is equivalent to `D`; at `q=2` the definitive test on a support is over `D^2`. |
| Trace descent / character averaging | Ch.1 | normalized trace over finite extensions; `(T_gamma v)_i = mu_gamma (v_i∘gamma)/(nu_gamma gamma_i)` averaging. |
| Shioda–Katsura split blowup | external (Tohoku Math J. 31, Thm 1.7, Rem 1.10) | blowup of `S_6 x S_6` along `C x C`, diagonal `mu_6` quotient = blowup of split fourfold along two disjoint pure-block curves. |
| Kuga–Satake model | Floccari §5.2, Rem 5.5 | `W=H^1(A)=C^+(T)`, `j: T(1) -> End(W)`, `j(v)(w)=v w e`, Hodge embedding; `End(W)(-1) ≅ W⊗W ⊂ H^2(A x A)`. |
| Algebraic KS for order-3 non-symplectic K3 with fixed curve | Bolognesi–Laterveer (arXiv 2310.11981, Thm 4.4, §6) | supplies the algebraic cycle for the "unconditional order-three family". |
| Weighted Griffiths residue grading for `S_6` | "stated in Chapter 3" per the text (no section anchor given; a grep for "Steenbrink"/"residue grading" over `docs/*.md` hits only 3A itself) | Hodge pieces of `H^2(S_6)_prim` sit in total Jacobian degrees 2, 8, 14 of `R_f ⊗ k[u]/(u^5)`. |

---

## 1. §3A.1–3A.3 The split (block-sign) theorem — elliptic companion + contraction

### Statement [proved in book, modulo Ch.3 inputs]
**Split theorem.** `F = f(X) + g(Y)` with `f, g` ternary homogeneous sextics, `X = V(F) ⊂ P^5` smooth, `tau[X:Y] = [X:-Y]`. Then all of `H^4(X)^tau` is divisor-supported, and over an effective characteristic-zero field the supporting divisor and legal certificates **for every positive `q`** can be constructed by finite algebraic operations. Includes `f = g`.

### Hypotheses on F, precisely
- Block structure: `F = f(X0,X1,X2) + g(Y0,Y1,Y2)`; no cross terms; `f, g` arbitrary ternary sextics (no Fermat/diagonal condition, no coefficient condition beyond smoothness).
- Smoothness of `X` ⟺ smoothness of both plane curves `C_f = {f=0}` and `C_g = {g=0}` (partials of `F` are exactly the two groups of partials of `f` and `g`).
- Symmetry group used: block sign `tau`; auxiliary `mu_3` deck actions on triple planes and the Fermat elliptic curve. No further symmetry of `f, g` is required.

### Targets served
Every Y-odd `A` of degree `6q`, every `q >= 1` (geometric all-pole statement). Characters: the invariant residue sector decomposes as `V_{f,zeta}⊗V_{g,zeta^2} ⊕ V_{f,zeta^2}⊗V_{g,zeta}` plus three Tate classes.

### The denominator `D`
- **Not explicit.** `D` is the computed support of the covering-curve algorithm on a resolution `W~_f` of the Fano companion (see §2 below), transferred through the quotient maps (graph + exceptional boundaries included), intersected/projected via the translated diagonal, transferred to `X/tau` and `X`, then enlarged to a homogeneous hypersurface, multiplied by Galois conjugates and by its `tau`-image. Result: rational, globally Y-even, `F ∤ D`.
- Degree: **no bound stated** ("supplies no bound on the degree `e` of the support"). Numerator then solved by the Ch.3 §3.6 fixed-denominator linear system over `D^q`.
- For the rational-function realizations actually accepted (crossed-ternary q=1, split-ternary q=2) the denominators are products of powers of the affine ternaries `p = f(x,a,1)`, `q0 = g(r,z,1)`, the chart variables `t, lambda`, integer constants, and norms of companion-field denominators after trace; **neither accepted proof exports `deg D`** (existential certificates). Numerator degree `deg D + 1` (q=1) and `deg D + 7` (q=2).

### Key geometric objects and numerical facts (all [proved in book])
- Triple plane `T_f = {u^3 = f(z0,z1,z2)} ⊂ P(1,1,1,2)`, `rho(z,u) = (z, zeta u)`, `zeta^3=1≠zeta`; avoids the weighted vertex; smooth iff `C_f` smooth.
- `pi_* O_{T_f} = O ⊕ O(-2) ⊕ O(-4)` on `P^2` ⇒ `H^1(T_f) = 0`, `H^3(T_f)=0` by duality; `T_f/rho = P^2` ⇒ invariant part of `H^2(T_f)` is the hyperplane line.
- Branch curve `C_f`: genus 10, `e(C_f) = -18`; `e(T_f) = 3·3 - 2·(-18) = 45`, `b_2 = 43`, `dim V_f = dim H^2(T_f)_prim = 42`, each nontrivial eigenspace has dim 21.
- `K_{T_f} = O(1)` ⇒ `h^{2,0} = 3`; affine residue denominator `u^2` ⇒ every holomorphic 2-form has `rho`-eigenvalue `zeta`. Hodge triples `(h^{2,0},h^{1,1},h^{0,2})`: `V_{f,zeta} = (3,18,0)`, `V_{f,zeta^2} = (0,18,3)`.
- Fermat elliptic curve `E = {a^3+b^3+c^3=0}`, `epsilon[a:b:c] = [a:b:zeta c]`; on `c=1`, `epsilon(a,b) = (zeta^{-1}a, zeta^{-1}b)`; holomorphic form `da/b^2` has eigenvalue `zeta`.
- `U_f = [V_f ⊗ H^1(E)]^{diag mu_3}` has only types `(2,1),(1,2)`, each of dim 21: no invariant `(3,0)` because both holomorphic factors carry the same character `zeta`. (This is the motivation; support is proved geometrically.)

### §3A.2 Fano elliptic companion
- `W_f = {f(z) + s^3 + t^3 = 0} ⊂ P(1,1,1,2,2)`.
- Rational map `T_f x E ⇢ W_f`, `([z:u],[a:b:c]) ↦ [c z : u a c : u b c]`; substitution gives `c^6 f(z) + u^3 c^3 (a^3+b^3) = 0`; on `c=1`: `s = u a`, `t = u b`. Deck transformation `(u,a,b) ↦ (zeta u, zeta^{-1} a, zeta^{-1} b) = (rho, epsilon)`. Inverse: adjoin `u^3=f`, `a = s/u`, `b = t/u`. Generically degree 3; `W_f` birational to `(T_f x E)/diag mu_3`.
- Singularities: cone gradient `(∇f, 3s^2, 3t^2)` smooth off vertex; ambient singular locus is the weighted line `z=0`; it meets `W_f` in the 3 points `s^3+t^3=0`; local type `1/2(1,1,1)` (age 3/2 > 1, terminal). These are the only singularities.
- `-K_{W_f} = O(1)`; `O(2)` is very ample: 6 quadratic `z`-monomials + `s,t` give a closed embedding in `P^7`; every sextic `z`-monomial is a product of three quadratic ones so the equation lives in that presentation.
- Normal-Fano bound: covering curves of anticanonical degree `<= 6` (=`2n`, `n=3`) ⇒ embedding degree `<= 12`. Algorithm: search complete Hom schemes through degree 12 in the `P^7` embedding, **retain curves through the three quotient points**, lift generic curve of a dominant component to a resolution `W~_f`, spread out, choose smooth affine parameter surface with generically finite `P^1 x S -> W~_f`; Ch.3 gives a divisor supporting all `H^3(W~_f)`.
- Quotient transfer ⇒ divisor `D_f ⊂ T_f x E` supporting `U_f = H^3(T_f x E)^{diag mu_3}` (uses `H^1(T_f)=H^3(T_f)=0` and `H^1(E)^epsilon = 0`).

### §3A.3 Contraction of the elliptic factors
- `D_f x D_g` has dim 4 in `(T_f x E) x (T_g x E)`; intersect with pullback of a **general translated diagonal** `Delta_a ⊂ E x E`, project to `T_f x T_g`: dim `<= 3`. Properness: apply the difference map `E x E -> E` to each component of `D_f x D_g`; a general value avoids the image or drops fibre dimension by one. The translation is selected by a finite dimension test over a finite extension.
- Translation acts trivially on cohomology ⇒ cup with `[Delta_a]` then pushforward: `(U_f ⊗ U_g)(1) -> [V_f ⊗ V_g]^{diag mu_3}` (H^6 of sixfold → H^4 of fourfold: +2 from cup, −4 from projecting along `E x E`). Image supported on the proper projection of the 3-dimensional intersection.
- Over C: `U_f = V_{f,zeta}⊗H^1(E)_{zeta^2} ⊕ V_{f,zeta^2}⊗H^1(E)_zeta`; the two elliptic eigenlines pair nondegenerately (Poincaré duality) ⇒ contraction surjects onto `V_{f,zeta}⊗V_{g,zeta^2} ⊕ V_{f,zeta^2}⊗V_{g,zeta}` = the whole invariant tensor. All maps rational ⇒ support descends.
- Quotient comparison (cover coordinates): on `X0 Y0 ≠ 0`, `x_i = X_i/X0`, `y_i = Y_i/Y0`, `t = Y0/X0`; `f(1,x) + t^6 g(1,y) = 0`, `tau(t) = -t`; invariant `r = t^2` with `r^3 = -f/g`; on `T_f x T_g` (`u^3=f`, `v^3=g`) the diagonal invariant `r = -u/v` satisfies the same equation ⇒ `X/tau` birational to `(T_f x T_g)/diag mu_3`.
- Künneth: `H^4(T_f x T_g)^{diag mu_3} = Q(-2)^3 ⊕ [V_f ⊗ V_g]^{diag mu_3}`. Tate terms algebraic; tensor term supported by contraction. Transfer to a smooth model of the quotient, then to `X/tau`, then `X`. Then Ch.3 §3.6 for the certificate at every `q`. **No Hodge conjecture for abelian varieties is used.**

### Stated limits
- Finite but not practical: no degree bound on `D`, no runtime bound, no strictness promised.
- The geometric all-pole theorem is distinct from a formal proof of a particular higher-pole target (see §3 below).

---

## 2. "Algebraic contraction with moving poles" — the rational-function realization (q=1 template)

[proved in book at the level of identities; full formalization in PR #149 / `HodgeComputer.Constructions.EllipticCompanions`]

Setting: Fermat elliptic function field with invariant derivation `delta`; two companion constructions give scalar identities
```
delta f + div_x C = alpha            (first block; coefficients constant for the second block's derivatives)
delta g + div_y G = beta * theta     (second block; conversely)
```
- A generic elliptic **translation** separates the finite poles of the first primitive from the second block's pole set. Must be an actual field map commuting with the derivations (numerical avoidance of finitely many points is not enough).
- Finite residue functional `R` with `R(delta h) = 0`, `R(d_i h) = d_i R(h)`, `R(g) = 0` (the last by proved pole separation).
- Pairing: `R(f theta) = 3 alpha - div_x xi`, where `xi` = finite theta-pole evaluations + an **infinity correction** built from a finite rational approximation beyond the pole order (omitting it gives an incorrect finite-pole formula).
- Product rule:
```
U_i = ( R(C_i g) + beta xi_i ) / 3,   V_j = R(f G_j) / 3,   div_x U + div_y V = alpha beta.
```
- Descent: finite algebraic extensions by normalized trace; the translation parameter removed by trace to its rational coordinate line followed by the constant Laurent coefficient (both commute with coefficient derivatives). Then coordinate transport and legal clearing.

Explicit data from the accepted crossed-ternary proof (`data/crossed-ternary-asymmetric/solution.md`) that realizes this template — **[from the accepted solution, not the chapter]**:
- `a_e^3 + b_e^3 + 1 = 0`, `delta a_e = b_e^2`, `delta b_e = -a_e^2`; Weierstrass presentation `y_e^2 = x_e^3 - 432`, `a_e = -(y_e+36)/(6 x_e)`, `b_e = (y_e-36)/(6 x_e)`.
- On the first triple plane `S = u a_e`, `T = u b_e`; from companion vector `(V_x,V_a,V_S,V_T)`: `f = 3(b_e V_S - a_e V_T)/u^2`, `C_x = -3V_x/u`, `C_a = -3V_a/u`; scalar equation `delta f + d_x C_x + d_a C_a = B/u^4 - (a_e B_S + b_e B_T)/u^3`.
- `theta = b_e/a_e^2`; sign-fix identity `delta(b_e^2/a_e) = b_e/a_e^2 - a_e b_e`.
- Theta kernels for `h = h0(x_e) + h1(x_e) y_e`: `K0 = -36 x_e (x_e^3+3456)/(x_e^3-1728)^2`, `K1 = 3888 x_e^4/(x_e^3-1728)^2`; theta fibres `r^3 = 1728`; local residue `Res_r(h0 K0 + h1 K1) = -12 h0'(r) + 432 h1'(r) + (r^2/2) h1(r) = -delta h(r,-36)`.
- Infinity functional `B` with `B(1)=0`, `B(d_i h) = d_i B(h)`, `InfCoeff Tr(f theta omega) = B(delta f)`; triangular recurrence with denominators `5+2n`; `xi_i = sum_{r^3=1728} C_i(r,-36) + B(C_i)`.
- Weighted divergence after multiplying by `u^4 v^2` and trace: `sum_i ( d_i W_i - (4 d_i p/(3p)) W_i - (2 d_i q/(3q)) W_i ) = s`, `s = x^4 a r`, `alpha = x^4 a/u^4`, `beta = r/v^2`.
- Physical chart: `X0 = t x, X1 = t a, X2 = t lambda, Y0 = t, Y1 = t lambda r, Y2 = t lambda z`; `F = t^6 (p + lambda^6 q)`, `A = t^6 lambda s`; `P_d = sum W_i d_i p`, `Q_d = sum W_i d_i q`; `w_i = -lambda W_i/(3p)`, `w_lambda = (3 p s + P_d + lambda^6 Q_d)/(18 p lambda^4 q)`; physical field `V = (t w_x, t w_a, t w_lambda, 0, t(lambda w_r + r w_lambda), t(lambda w_z + z w_lambda))`; chart Jacobian adds `2 w_lambda/lambda`; identities `V(F) = A`, `div V = F · ( -(3ps + P_d)/(9 p^2 q lambda^5 t^6) )`.
- Legality: track fractions simultaneously in the ambient chart field and on `F=0` via `lambda_0^6 = -p/q`; irreducibility of `F` by specialization `X0=X1=Y1=0, X2=T-1, Y2=i` → `Y0^6 + (T-1)^6 - 1`, Eisenstein at `T=0`.

Stated limits (chapter): the realization "does not turn a fixed pole-one primitive into an arbitrary higher-pole certificate by multiplication"; must preserve companion primitives, pole separation, all residue boundaries, derivative-compatible descent, exact target.

---

## 3. "Higher-pole residue calculus" — third-order pole (q=2) via second normal residue

[proved in book as identities; formalized in PR #159 / `HodgeComputer.Constructions.HigherResidues`]

- `H` hypersurface equation in five affine coordinates, normal coordinate `x0`, `g = d_0 H ≠ 0` in the function field. For polynomial numerator `B`:
```
R_1(B) = d_0 (B/g),
R_2(B) = (1/2) d_0 ( R_1(B)/g )
       = (1/2) ( B_00/g^2 - (3 B_0 H_00 + B H_000)/g^3 + 3 B H_00^2/g^4 ).
```
  Second residue density is `R_2(B)/g` in tangent coordinates. Not by itself a polynomial-curve primitive (denominators along the family, varying Jacobian).
- Polynomial-curve family with polynomial normal `psi`, `dH(psi) = 1`; thicken `c(t) + eps psi(t)`:
```
H(c + eps psi) = eps + a eps^2 + b eps^3 + ...,
B(c + eps psi) J(eps) = K_0 + K_1 eps + K_2 eps^2 + ...,
P = K_2 - 3 a K_1 + (6 a^2 - 3 b) K_0,
P = J(0) ( R_2(B) + ResDiv_g Z ),   ResDiv_g Z = g sum_i d_i (Z_i/g)   (Z = tangent correction).
```
  Integrate `P` in the curve parameter, Piola-transport, subtract `Z` ⇒ primitive for `R_2(B)`; normalized trace descends through the curve family's finite extension. Polynomiality of numerator, normal and coordinate derivatives is essential.
- Split sextic, target of block bidegree `(7,5)`: on `H = 7(f + S^3 + T^3)`, `h` of degree 7 in the ordinary block variables, `B = 49 h S T`, normal coordinate `T`, `g = 21 T^2`:
```
(1/2) R_1(B) = -7 h S/(6 T^2),    R_1( (1/2) R_1(B) ) = R_2(B).
```
  (Accepted proof adds: tangent equation right side `-2hS/(9T^5)`; `T ≠ 0` proved in the companion field.)
- Weighted radial descent: `(T,x,a,y,S) = (rho^2 T0, rho x0, rho a0, rho, rho^2 S0)`, cone weights `(1,1,1,2)`, normal `21T^2` weight 4, source weight −1 ⇒ affine coefficient divided by `rho`; radial residue gives the affine primitive. (Specific weights justify the descent; not a general assertion.)
- Elliptic identity converting to the theta character: `delta(a^2/b + b^2/a) = b/a^2 - a/b^2` (Fermat coordinates `a_e^3+b_e^3+1=0`); scalar source of first companion `-h a_e/(3 u^5 b_e^2)` → theta source `(h/u^5)(b_e/a_e^2)`; second companion gives ordinary source `r^3 z^2/v^4`.
- After contraction and normalized trace the four-vector `W` satisfies
```
div W - 5 W(p)/(3p) - 4 W(q0)/(3 q0) = s     (p, q0 affine ternaries, s affine target)
```
  Substitution `X = t(x,a,1)`, `Y = t lambda (r,z,1)` gives `F = t^6 (p + lambda^6 q0)`; the vertical component and the Jacobian contribution `2 V_lambda/lambda` must both be retained; legal homogeneous clearing gives a **strict** certificate with numerator degree `deg D + 7`, `Q = 0`.
- Instance: `split-ternary-q2` (q=2, `F = f(X)+g(Y)` with `c=5/7`, `f = x^6+a^6+y^6 + c(x^5 a + x^2 a^3 y + 3 x y^5 + 2 a^5 y)`, `g = r^6+z^6+w^6 + c(2 r^5 z + 2 r^3 z w^2 + r w^5 + 3 z^5 w)`, `A = X0^4 X1^3 Y0^3 Y1^2`, block degrees (7,5)). Accepted (fqlx, PR #159). Moving-conic seed modulo 5, 13 coefficient equations with Jacobian determinant 1 in F_5, 12-dim free/normal block determinant 4.

Stated limits: "reusable conditional steps, not a theorem producing the required family on every hypersurface"; the geometric all-pole statement does not remove scalar-weight, descent, or denominator obligations of another formal application.

---

## 4. "Involution-stable ternary decompositions" — the crossed case

[proved in book, nonconstructive; constructive realization = PR #149]

- Hypothesis: `F = p(U) + q(V)` for complementary 3-dimensional subspaces `U, V` preserved by `tau`. If all three negative eigendirections lie in one summand → block-sign split theorem. Otherwise multiplicities (1,2) after exchange; normal form `p(X0,X1,Y0) + q(X2,Y1,Y2)`; global Y-evenness forces `p` even in `Y0` and `q` even in `X2`.
- Let `r` negate `Y0`, `s` negate `X2`, `rho` negate the entire second block; `tau = rho r s`. Reflection theorem supports `r`-invariant and `s`-invariant cohomology; split theorem supports `rho`-invariant cohomology. On a simultaneous eigenspace with `tau = +1`: at least one of `r, s` is `+1`, or `r = s = -1` and so `rho = +1`. Three supported spaces exhaust the invariants.
- Serves: every Y-odd target, every positive `q` (existence). Character bookkeeping matters: a target even under both `r` and `s` (residue anti-invariant under each) lies only in the `rho`-sector and needs the split/elliptic mechanism, not reflections.
- Instance **[my cross-reference]**: `crossed-ternary-asymmetric` (q=1, `F = P(X0,X1,Y0) + P(Y1,Y2,X2)`, `P(u,v,w) = u^6+v^6+w^6 + (5/7)(u^5 v + u^3 v w^2 + u v w^4 + 2 u v^5)`, `A = X0^4 X1 Y1`, block degrees (5,1)); solved by `1337isnot1337` (PR #149) using exactly §2 above.

---

## 5. §3A.4 Block swap: `F_f = f(X+Y) + f(X-Y)`

### Hypotheses
`f` arbitrary smooth ternary sextic; `U = X+Y`, `V = X-Y`; the global involution becomes the **block swap** `(U,V) ↦ (V,U)`, not the block sign. The split theorem does not apply.

### Sextic cover and its eigenspaces [proved in book]
- `S_6 = {u^6 = f(z)} ⊂ P^3`, `delta(z,u) = (z, eta u)`, `eta` primitive 6th root; `H_k = (H^2(S_6,C)_prim)_{delta^* = eta^k}`, `1<=k<=5`.
- Jacobian ring `R_f ⊗ k[u]/(u^5)`, `Hilb(R_f,t) = (1+t+t^2+t^3+t^4)^3`; residue volume contributes `eta`, so `u^{k-1}` has character `eta^k`; Hodge pieces in total Jacobian degrees 2, 8, 14; read `R_f` at degrees `2-(k-1), 8-(k-1), 14-(k-1)`:

| k | (h^{2,0},h^{1,1},h^{0,2}) | rank |
|---|---|---|
| 1 | (6,15,0) | 21 |
| 2 | (3,18,0) | 21 |
| 3 | (1,19,1) | 21 |
| 4 | (0,18,3) | 21 |
| 5 | (0,15,6) | 21 |

- `H_3` is rational and is the pullback of `V_2 = H^2(S_2)_prim` from the K3 double plane `S_2 = {w^2 = f} ⊂ P(1,1,1,3)` via `w = u^3`.

### Algebraic split decomposition [proved in book + Shioda–Katsura]
- `alpha^6 = -1`, `alpha^2 = eta`. `Phi: S_6 x S_6 ⇢ {f(U)+f(V)=0}`, `((z,u),(z',v)) ↦ [v z : alpha u z']`; generic fibres = diagonal `mu_6`-orbits; indeterminacy `C x C`, `C = {f=0}`. Blow up `C x C`; chart `v = u lambda`: action `u ↦ eta u`, `lambda` fixed; eliminating `u^6 = f(z)` gives quotient chart `f(z') = lambda^6 f(z)`, induced map `[lambda z : alpha z']` = blowup chart of the target along its pure-block copy of `C`; chart `u = v lambda'` gives the other copy. So the smooth quotient `B` = blowup of the split fourfold along its two disjoint pure-block curves.
- Blowup bookkeeping: `S_6 x S_6` blown up along `C x C` adds `H^2(C x C)(-1)` = 2 Tate + `(H^1(C)⊗H^1(C))(-1)`; original invariant `H^4` = 3 Tate + primitive invariant tensor; each pure-block `C` (codim 3) adds `H^2(C)(-1) ⊕ H^0(C)(-2)` = 2 Tate per curve. Cancelling:
```
H^4({f(U)+f(V)=0})_prim ≅ [H^2(S_6)_prim ⊗ H^2(S_6)_prim]^{diag mu_6} ⊕ [H^1(C)⊗H^1(C)](-1).
```
  Maps/projectors algebraic (blowup maps, finite averaging, hyperplane projectors, projective-bundle formula). The joining term is supported on the **joining divisor** swept by lines connecting points of the two copies of `C` (image of the exceptional `P^1`-bundle over `C x C`).
- Swap lift `j = (delta^{-1}, id) ∘ flip` (`Phi ∘ j` = swap of `Phi` after multiplication by `alpha`); `j^2 = (delta^{-1},delta^{-1})` trivial on the diagonal quotient. On `H_k ⊗ H_{6-k}`, `j^* = eta^{-k} ·` (tensor flip); at `k=3` it is `-flip` ⇒ swap-fixed K3 contribution is `Λ^2 V_2`. On the joining term the graded sign gives fixed part `Λ^2 H^1(C)(-1)`. Exterior squares are forced by these signs.

---

## 6. §3A.5 Second elliptic companion and supported ranks

### Pairs `H_2⊗H_4`, `H_4⊗H_2` [delegated to §1]
Pull back the triple-plane support along `S_6 -> T_f`, `u ↦ u^2`.

### Pairs `H_1⊗H_5`, `H_5⊗H_1` [proved in book modulo Ch.3]
- `E_23 = {a^2 + b^3 + c^6 = 0} ⊂ P(3,2,1)`, `epsilon(a,b,c) = (a,b,eta c)`; holomorphic differential `db/a` on `c=1` transforms by `eta`.
- `W_{23,f} = {f(z) + s^2 + t^3 = 0} ⊂ P(1,1,1,3,2)`.
- Map `S_6 x E_23 ⇢ W_{23,f}`, `((z,u),(a,b,c)) ↦ [c z : u^3 a : u^2 b]`; substitution gives `u^6 (c^6 + a^2 + b^3)`; on `c=1`, `a = s/u^3`, `b = t/u^2`; diagonal `mu_6` quotient of degree 6.
- Smoothness: cone gradient `(∇f, 2s, 3t^2)` has no nonzero critical zero; ambient singularities = isolated weight-3 and weight-2 vertices, excluded by the equation; where `s,t ≠ 0` weights coprime ⇒ trivial stabilizer. `-K = O(2)`, `O(1)` ample.
- Curve bound (smooth threefold, `-K.C <= 4`): weighted map degree `e <= 2`. Incidence: three binary `e`-forms, one binary `3e`-form, one binary `2e`-form: `8e+5` coefficients, `6e+1` equations, `6e x 8e` freeness matrix; basepoint exclusion applies to all five forms collectively (the three weight-one forms may vanish together where the other two are nonzero).
- Support for `H^3(W_{23,f})` from Ch.3; quotient transfer supports `[(H_1 ⊕ H_5) ⊗ H^1(E_23)]^{diag mu_6}` (other characters pair with no elliptic eigenspace); translated-diagonal contraction as in §3A.3 supports `H_1⊗H_5` and `H_5⊗H_1` — the whole tensor, not just leading Hodge piece.

### Complete invariant decomposition for `F_f` [proved in book]

| Invariant summand | Rank | (h^{3,1},h^{2,2},h^{1,3}) | Support |
|---|---|---|---|
| Primitive sixth-root pair (`H_1⊗H_5`) | 441 | (90,261,90) | second elliptic companion |
| Primitive third-root pair (`H_2⊗H_4`) | 441 | (54,333,54) | triple-plane companion |
| Joining-curve exterior square `Λ^2 H^1(C)(-1)` | 190 | (45,100,45) | joining divisor |
| K3 exterior square `Λ^2 V_2` | 210 | (19,172,19) | **Unresolved in general** |

Checks given: `Λ^2 V_2` middle dim `1·1 + C(19,2) = 172`; `Λ^2 H^1(C)` extremes `C(10,2) = 45`. Supported primitive rank 1072 with (189,694,189); +ambient = 1073; full primitive (208,866,208). No generic Picard-rank assumption; extra algebraic K3 classes only add already-supported exterior-square terms.

### Polynomial test for the constructed part [proved in book]
- For globally Y-odd `A` of degree `6q`: `Â(U,V) = A((U+V)/2, (U-V)/2)` (antisymmetric under swap). `P_0` = projector retaining monomials whose total `V`-degree is `≡ 0 mod 6`.
```
P_0 Â = 0  ⇒  the residue target has the constructed divisor support.
```
- Proof: `h(U,V) = (U, eta V)`; a monomial of `V`-degree `b` has residue character `eta^{b+3}` (3-coordinate volume contributes `eta^3`); `h` lifts under `Phi` to `(delta, id)`; the K3 summand has character `eta^3`, so receives exactly `b ≡ 0 mod 6`. Finite-group actions commute with residue and pole reduction ⇒ test on the full class at arbitrary `q`.
- At `q=1`: the excluded terms are exactly `B(U) - B(V)` with `B(U) = Â(U,0) = 2^{-6} A(U,U)`, so the test is `A(X,X) = 0`. Counts: 198 mixed antisymmetric sextic monomials; 9 independent mixed Jacobian directions `f_i(U) V_j - f_i(V) U_j` (partials form a regular sequence ⇒ no degree-one syzygies); `198 - 9 = 189` = supported (3,1) count.
- For `q > 1`, intermediate bidegrees such as (6,6) can occur; `A(X,X)=0` alone is **not** the stated test — use `P_0 Â` directly.
- **[computational evidence]** sympy: for a random Y-odd sextic `A`, `P_0 Â == B(U)-B(V)` and `B(U) == A(U,U)/64` (True); mixed sextic monomials 406, swap-fixed 10, antisymmetric dimension 198 (matches).

---

## 7. §3A.6 Algebraic Kuga–Satake criterion

### Exterior-power theorem [proved in book, given the KS model]
If a projective K3 `S` has an **algebraic** Kuga–Satake correspondence (algebraicity of the standard Clifford embedding `j`, not mere existence of a KS Hodge structure), then `Λ^r H^2(S) ⊂ H^{2r}(S^r)` is supported in codimension `>= r-1` for every `r >= 1`. For `r=2` it supports `Λ^2 V_2` and completes all targets at every positive `q` on `F_f`.

- Clifford calculation: `j(v_1)⋯j(v_r) = L_{v_1⋯v_r} R_{e^r}`; antisymmetrize; symbol map `Alt: Λ^r T -> C(T)` injective; `J_r = L_Alt R_{e^r}` injective (`r` even: `e^r` nonzero scalar, left multiplication faithful on `C^+`; `r` odd: evaluate at 1, `Alt·e` injective since `e` invertible). `J_r: Λ^r T(r-1) ↪ H^2(A x A)`. For `r=2`: `[j(v),j(w)] = q(e) L_{vw - wv}`.
- Algebraicity ⇒ `J_r` induced by `Γ ∈ CH^{r+1}(S^r x B)_Q`, `B = A x A` (external products, contractions, antisymmetrizer algebraic; projectors: `H^1(A)` via multiplication-by-`n` graph polynomial; polarization via isogeny; degree-2 projector on `S` from diagonal and point classes; transcendental projector subtracts the NS projector via inverse intersection matrix; permutation graphs for exterior power).
- Smooth complete-intersection surface `Σ ⊂ B`; weak Lefschetz ⇒ `Λ^r T(r-1) ↪ H^2(Σ)`; polarization adjoint + exterior-power projector = algebraic correspondence `Σ -> S^r` surjecting onto `Λ^r T` (Gram operator nondegenerate; inverse by Cayley–Hamilton; polarizations algebraic by Lefschetz decomposition on a surface). Adjoint cycle codim `r+1` in `Σ x S^r`, support dim `<= 2 + 2r - (r+1) = r+1`, projection to `S^r` has codim `>= r-1`. Terms with `p` NS factors gain `p` extra codimensions.

### Unconditional order-three family [delegated to Bolognesi–Laterveer]
- `f(z0,z1,z2) = g_6(z0,z1) + z2^3 g_3(z0,z1) + a z2^6` smooth. On `S = {w^2 = f}`, `z2 ↦ zeta_3 z2`, `w ↦ w` is non-symplectic of order 3 (multiplies the 2-form by `zeta_3`); fixed locus contains `z2 = 0`, the genus-2 curve `w^2 = g_6`. External theorem ⇒ algebraic KS ⇒ exterior-square theorem completes `F_f` for every allowed target and positive `q` — **existence only**.
- Smoothness criterion [proved in book]: `a ≠ 0`, `g_6` squarefree, and `Δ = g_3^2 - 4 a g_6` squarefree. (On `z2=0`: `g_6` squarefree. Off it, `t = z2^3`: singular point ⟺ `t = -g_3/(2a)`, `Δ=0`, `dΔ=0`; repeated root of `Δ` with `g_3 ≠ 0` gives a singular point; if `g_3` vanishes there, `g_6` already has a repeated root. `[0:0:1]` excluded iff `a ≠ 0`.)

### Stated limits (verbatim gist)
- "For arbitrary smooth `f`, algebraicity of the needed Kuga–Satake cycle is not established here. The general polarized-split case therefore retains its K3 remainder, and the full globally Y-even problem remains open."
- "The argument here does not turn that cycle into an explicit certificate." Fano companion sectors have the finite support construction; the K3 sector does not.
- **Direct boundary to the elliptic trick**: the double-plane involution acts by `-1` on K3 transcendental cohomology; for any curve with nonzero anti-invariant `H^1`, the anti-invariant `H^{1,0}` is nonzero; the diagonal invariant tensor contains `H^{2,0}(S) ⊗ H^{1,0}(C)_-`; a smooth quotient model has a holomorphic 3-form and cannot be uniruled. The third- and sixth-root companion selections avoid this piece; the K3 exterior square needs its antisymmetry inside a larger correspondence.

---

## 8. Catalogue of explicit formulas reusable computationally

Cover / cone coordinates
1. Triple plane: `T_f = {u^3 = f(z)} ⊂ P(1,1,1,2)`; companion `W_f = {f(z)+s^3+t^3=0} ⊂ P(1,1,1,2,2)`; map `([z:u],[a:b:c]) ↦ [cz : uac : ubc]`, inverse `a = s/u, b = t/u`; deck `(u,a,b) ↦ (zeta u, zeta^{-1}a, zeta^{-1}b)`.
2. Quotient chart of split `X/tau`: `x_i = X_i/X0, y_i = Y_i/Y0, t = Y0/X0`, `f(1,x) + t^6 g(1,y) = 0`, `r = t^2`, `r^3 = -f/g`, and `r = -u/v` on `T_f x T_g`.
3. Sextic cover: `S_6 = {u^6 = f} ⊂ P^3`, `delta = (z, eta u)`; K3 `S_2 = {w^2=f} ⊂ P(1,1,1,3)`, `w = u^3`.
4. Swap-split map: `Phi((z,u),(z',v)) = [v z : alpha u z']`, `alpha^6=-1, alpha^2=eta`; blowup chart `v = u lambda`, quotient chart `f(z') = lambda^6 f(z)`, induced `[lambda z : alpha z']`; swap lift `j = (delta^{-1}, id)∘flip`.
5. Second companion: `E_23 = {a^2+b^3+c^6=0} ⊂ P(3,2,1)`, `W_{23,f} = {f+s^2+t^3=0} ⊂ P(1,1,1,3,2)`, map `[cz : u^3 a : u^2 b]`, inverse `a = s/u^3, b = t/u^2`.
6. Physical chart for split/crossed transport: `X = t(x,a,1)`, `Y = t lambda (r,z,1)` (crossed variant `X0=tx, X1=ta, X2=t lambda, Y0=t, Y1=t lambda r, Y2=t lambda z`), `F = t^6(p + lambda^6 q0)`; Jacobian contribution `2 V_lambda/lambda` (resp. `2 w_lambda/lambda`).
7. Weighted radial descent on the companion cone: `(T,x,a,y,S) = (rho^2 T0, rho x0, rho a0, rho, rho^2 S0)`.

Elliptic identities (Fermat cubic `a^3+b^3+1=0`, `delta a = b^2`, `delta b = -a^2`)
8. `delta(a^2/b + b^2/a) = b/a^2 - a/b^2` (used for the (7,5) theta conversion).
9. `delta(b^2/a) = b/a^2 - a b` (sign fix, crossed case); `theta = b/a^2`.
10. Weierstrass model `y^2 = x^3 - 432`, `a = -(y+36)/(6x)`, `b = (y-36)/(6x)`; theta kernels `K0 = -36x(x^3+3456)/(x^3-1728)^2`, `K1 = 3888x^4/(x^3-1728)^2`; theta fibres `r^3 = 1728`; `Res_r(h0 K0 + h1 K1) = -12 h0'(r) + 432 h1'(r) + (r^2/2) h1(r) = -delta h(r,-36)`; infinity recurrence denominators `5+2n`.

Residue / contraction identities
11. `delta f + div_x C = alpha`, `delta g + div_y G = beta theta`, `R(f theta) = 3 alpha - div_x xi`, `U_i = (R(C_i g) + beta xi_i)/3`, `V_j = R(f G_j)/3`, `div_x U + div_y V = alpha beta`.
12. Weighted divergence equations after clearing cube roots: q=1 crossed `div W - (4/3) W(p)/p - (2/3) W(q)/q = s`; q=2 split (7,5): `div W - (5/3) W(p)/p - (4/3) W(q0)/q0 = s`. (Pattern: weights `(deg_p-side)/3`, i.e. `u^{-k}` sources give `-k W(p)/(3p)`.)
13. Second normal residue: `R_1(B) = d_0(B/g)`, `R_2(B) = (1/2)(B_00/g^2 - (3B_0 H_00 + B H_000)/g^3 + 3 B H_00^2/g^4)`; thickened jets `H = eps + a eps^2 + b eps^3`, `B·J = K0 + K1 eps + K2 eps^2`, density `P = K2 - 3a K1 + (6a^2 - 3b) K0 = J(0)(R_2(B) + ResDiv_g Z)`, `ResDiv_g Z = g sum_i d_i(Z_i/g)`.
14. (7,5) source: `H = 7(f+S^3+T^3)`, `B = 49 h S T`, `g = 21 T^2`, `(1/2)R_1(B) = -7hS/(6T^2)`, `R_1((1/2)R_1(B)) = R_2(B)`; scalar source `-h a_e/(3u^5 b_e^2)` → `(h/u^5)(b_e/a_e^2)`; second block `r^3 z^2/v^4`.

Character / projector formulas
15. `rho`-character of holomorphic 2-forms on `T_f`: `zeta`; on `E`, `da/b^2 ↦ zeta`; on `S_6`, `u^{k-1}` has character `eta^k` (volume contributes `eta`); on `E_23`, `db/a ↦ eta`; on `H_k ⊗ H_{6-k}`, `j^* = eta^{-k}·flip`.
16. K3-sector projector: `Â(U,V) = A((U+V)/2,(U-V)/2)`; `P_0` keeps monomials with `V`-degree `≡ 0 (mod 6)`; monomial of `V`-degree `b` has residue character `eta^{b+3}` under `h = (U, eta V)`; at `q=1`, `P_0 Â = B(U) - B(V)`, `B(U) = A(U,U)/64`.
17. Hilbert series `Hilb(R_f) = (1+t+t^2+t^3+t^4)^3`; Hodge degrees 2, 8, 14 of `R_f ⊗ k[u]/(u^5)`; eigenspace table in §5.
18. Order-three K3 family and smoothness test: `f = g_6 + z2^3 g_3 + a z2^6`, smooth ⟺ `a ≠ 0`, `g_6` squarefree, `Δ = g_3^2 - 4 a g_6` squarefree.

Search bounds (companion curve searches)
19. `W_f`: `-K = O(1)`, `O(2)` very ample to `P^7` (6 quadratic z-monomials + s,t), curves with `-K.C <= 6`, embedding degree `<= 12`, must retain curves through the three `1/2(1,1,1)` points.
20. `W_{23,f}`: `-K = O(2)`, `-K.C <= 4` ⇒ weighted degree `e <= 2`; three binary `e`-forms + one `3e`-form + one `2e`-form; `8e+5` coefficients, `6e+1` equations, `6e x 8e` freeness matrix.

Not present in this chapter (contrary to what one might expect from the task's examples): no Hessian/adjugate `E1` lift, no Koszul correction map, no Jacobian-of-cubic normalization. Those belong to Chapters 2/2a/4.

---

## 9. Boundaries summary table

| Construction | F hypotheses | Targets / q | Characters served | D | Limit |
|---|---|---|---|---|---|
| Split theorem (3A.1–3A.3) | `f(X)+g(Y)`, `C_f, C_g` smooth | all Y-odd A, all `q>=1` | both `V_{f,ζ}⊗V_{g,ζ²}`, `V_{f,ζ²}⊗V_{g,ζ}` + Tate | computed support, degree unbounded, Y-even rational | finite but impractical; nonconstructive degree |
| Moving-pole contraction (rational realization) | same, one target at a time | q=1 template (PR #149) | one scalar product `alpha·beta` | from clearing; `deg D` not exported | not a pole-multiplication device |
| Higher-pole residue calculus | same; needs polynomial curve family with polynomial normal | q=2 (third-order pole), block bidegree (7,5) done | theta × ordinary | strict; `deg N = deg D + 7`; `deg D` not exported | no family theorem on every hypersurface |
| Crossed decomposition | `p(X0,X1,Y0)+q(X2,Y1,Y2)`, `p` even in `Y0`, `q` even in `X2`, smooth | all Y-odd A, all `q` | r-, s-, rho-sectors | nonconstructive | none beyond split theorem |
| Block swap `f(X+Y)+f(X−Y)` | `f` smooth ternary | all `q` for the supported 1072-dim part; test `P_0 Â = 0` | `H_1⊗H_5`, `H_2⊗H_4`, `Λ²H^1(C)` | joining divisor + companion supports (not explicit) | `Λ²V_2` (210; (19,172,19)) unresolved |
| KS exterior-power criterion | `f` with algebraic KS for `w²=f`; unconditional for `g_6 + z2³g_3 + a z2^6` | completes all A, all `q` | `Λ²V_2` | none (existence of cycle only) | no explicit certificate; arbitrary `f` open |

---

## 10. Where this chapter says the open frontier instances fail

Instance names literally mentioned in the chapter:
- The only `#case/` anchor is **"Split ternaries · third-order pole" (`split-ternary-q2`)** — cited as the accepted implementation (fqlx, PR #159) of the higher-pole residue calculus. It is *solved*, not open.
- PR #149 (`1337isnot1337`, crossed-ternary construction) is cited for the kernels of the moving-pole contraction; this corresponds to the *solved* instance `crossed-ternary-asymmetric` (the chapter does not write the instance id).
- **No open instance id is named anywhere in the chapter.**

What the chapter does say fails, and which open instances that covers **[my cross-reference via `data/*/why.md`]**:

1. **Polarized split / block swap, K3 remainder** — the chapter's table marks the K3 exterior square `Λ²V_2` (rank 210, Hodge numbers (19,172,19)) "Unresolved in general"; §3A.6 ends: "For arbitrary smooth f, algebraicity of the needed Kuga–Satake cycle is not established here. The general polarized-split case therefore retains its K3 remainder, and the full globally Y-even problem remains open." It adds that even where the KS cycle exists (order-three family) "the argument here does not turn that cycle into an explicit certificate", and gives the structural obstruction: the elliptic trick cannot reach the K3 sector because the invariant tensor contains `H^{2,0}(S)⊗H^{1,0}(C)_-` (holomorphic 3-form on any smooth quotient model ⇒ not uniruled).
   - Open instance `polarized-k3-remainder` (q=1, `A = X0 X1 X2 Y0 Y1 Y2`) is exactly this family: **[computational evidence]** `F = (f(X+Y)+f(X-Y))/2` with `f(z) = z0^6 + (5/7) z0^5 z1 + (5/7) z0^2 z1^3 z2 + (15/7) z0 z2^5 + z1^6 + (10/7) z1^5 z2 + z2^6` — the same `f` as the first block of `split-ternary-q2`; its Jacobian ring Hilbert function is `1,3,6,10,15,18,19,18,15,10,6,3,1` = `(1+t+…+t^4)^3`, so `C_f` is smooth. The chapter's test fails here: `A(X,X) = X0^2 X1^2 X2^2 ≠ 0`, and `P_0 Â = (U0^2U1^2U2^2 - V0^2V1^2V2^2)/64 ≠ 0`, so the target genuinely reaches the K3 summand (agrees with its `why.md`). In the displayed coordinates `f` is **not** in the order-three normal form (its `z2`-degrees are {0,1,5,6}, not all `≡ 0 mod 3`), so the unconditional KS family does not apply as written; a hidden coordinate change is not excluded by the chapter or by `why.md`.
2. Other open instances in the catalogue (`cycle-q2`, `cyclic-triple-residual`, `dense-rational`, `diagonal-permutation`, `eight-plane-configuration`, `four-plus-two-residual`, `klein-coupled-character`, `off-paired`, `one-reflection-skew`, `quadratic-cubic-pullback`, `quadric-contact-residual`, `rank-one-biquartic`, `reflection-coupled-residual`, `short-power-presentation`, `symmetric-sextic-residual`) are **not addressed by Chapter 3A at all**: none has the ternary-split, crossed-ternary, or polarized-swap block structure the chapter requires, so the chapter neither supports nor obstructs them. (E.g. `symmetric-sextic-residual` has `S_3 x S_3` symmetry without a two-ternary splitting; its `why.md` points to Chapter 3 §§3.3–3.4, not 3A.)
3. Generic caution repeated in the chapter: the geometric all-pole theorems are existence statements; "the geometric all-pole statement does not remove the scalar-weight, descent, or denominator obligations of another formal application", and the rational realization "does not turn a fixed pole-one primitive into an arbitrary higher-pole certificate by multiplication". So for a *new* split or crossed target (any `q`), the proven route is: companion polynomial-curve family (finite-field seed → lift), scalar sources of the right weights (`u^{-k}` ↔ `-k W(p)/(3p)`), elliptic conversion to theta/ordinary characters, translated contraction with infinity correction, trace descent, chart transport with the Jacobian term, and legality by simultaneous ambient/hypersurface tracking.
