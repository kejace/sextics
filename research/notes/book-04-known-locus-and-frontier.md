# Technique map: Chapter 4, "The known locus and the open frontier"

Source: `/home/kejace/devel/sextics/docs/04-known-locus-and-frontier.md` (416 lines, read completely).
Cross-referenced definitions pulled from `docs/02-explicit-constructions.md` (§2.1, §2.4, §2.7, §2.9), `docs/03-geometric-support.md` (§3.1, §3.2, §3.6), `docs/05-computation.md` (§5.2), `docs/01-certificates.md`, `formal/README.md`, and `data/holdout/solution.md` ("Geometric scope of the Fermat cover"). Anything taken from outside Chapter 4 is marked `[§x.y]` or `[holdout/solution.md]`.

Status tags used below: **PROVEN** (theorem with proof in the book), **EVIDENCE** (finite computation, reproduced here), **CONJECTURE/OPEN** (stated as hypothesis or missing premise), **NOT DONE** (the chapter says explicitly that the step has not been carried out), **CHECKED** (I re-verified independently; script `/tmp/hodge/notes/check-04-counts.py`, output reproduced in §12 below).

Conventions. Variables `Z = (X0,X1,X2,Y0,Y1,Y2)`; `τ(X,Y) = (X,−Y)`; `F` smooth sextic with `F(X,−Y)=F(X,Y)`; target `A` is `Y`-odd of degree `6q`; certificate `v = N/D`, `deg D = d`, `deg N_i = d + 6q − 5`, `gcd(D,F)=1`, with
`K_A(v) = A − q·v(F) + F·div(v)` divisible by `F^{q+1}`. Cleared form (the Lean contract): `E = A D² − q D Σ N_i ∂_iF + F (D Σ ∂_iN_i − Σ N_i ∂_iD) = F^{q+1} Q`. Normalized system for fixed `D`: `E1: q·N(F) = A·D`, `E2: D·div N − N(D) ∈ (F^q)`. [Ch.1]

Chapter 4 is a *measurement* chapter: it (a) measures which coefficient families and which cohomological representations are covered by the constructions of Chapters 2–3a, (b) proves two geometric support theorems on coupled families (paired-sign surface quotient §4.2, cyclic-seven weighted Fano §4.3), (c) proves the generic-target reduction (§4.4), the `q`-raising identity (§4.5) and a deck/exchange descent criterion (§4.6), and (d) identifies the transcendental local systems that block the natural quotients of a *general* even sextic (§4.7). Almost all certificates it produces are existential (denominator = equation of a divisor support, not written down); the only explicit certificate recipe inside the chapter is the strict `q=1` cubic-surface construction of §4.3.

---

## 1. §4.1 Coefficient spaces, target spaces, orbit dimensions

### 1.1 Counting formulas  **PROVEN / CHECKED**

- `S = C[X0,X1,X2,Y0,Y1,Y2]`, `±` = total `Y`-parity. Summand of `Y`-degree `j` in `S_6` has dimension `C(8−j,2)·C(j+2,2)`.
  `dim S_6^+ = 28+90+90+28 = 236`, `dim S_6^- = 63+100+63 = 226`. Parameter space of even sextics: `P^235`.
- Smooth `F` ⇒ the six quintic partials form a regular sequence (homogeneous s.o.p.). Koszul resolution with 3 even + 3 odd generators gives, for `R_F = S/J(F)`:
  `H(t) = (1+t+t²+t³+t⁴)^6`, `Θ(t) = (1+t²+t⁴+t⁶+t⁸)^3` (τ-trace; `Θ = (1−t^5)^3(1+t^5)^3/((1−t)^3(1+t)^3)`),
  `dim (R_F)_d^± = ([t^d]H ± [t^d]Θ)/2`.
  Values: `d=6: (426 ± 10)/2 → R^+ = 218, R^- = 208`; `d=12: (1751 ± 19)/2 → 885 / 866`; `d=18: 218 / 208`.
- Griffiths ⇒ invariant primitive Hodge numbers `(h^{3,1}, h^{2,2}, h^{1,3}) = (208, 866, 208)`; the ambient `h²` adds one algebraic middle class.
- Degree-6 check: `226 − 18 = 208`; the 18 odd Jacobian directions are `Y_j F_{X_i}` and `X_i F_{Y_j}` (no linear syzygies among the quintic partials). At Fermat these have unique leading monomials `6 Z_i^5 Z_j`, so the 208 `Y`-odd sextic monomials with all exponents ≤ 4 are a basis of the primitive target space [holdout/solution.md].

### 1.2 Orbit-dimension recipe  **PROVEN**

- `G = GL3(X) × GL3(Y)`, `dim G = 18`. For a linear family `V` preserved by a connected `H ⊂ G`: `dim(G·V) ≤ dim V + 18 − dim H`. Equality is proved by a differential-rank computation at one point.
- At Fermat `F0 = Σ Z_i^6` the 18 orbit derivatives are the distinct monomials `6 Z_i^5 Z_j`, `i,j` in the same parity block (9 + 9). For a *monomial* family containing `F0`, count the `Z_i^5 Z_j` outside `V`: that number added to `dim V` is the saturated affine dimension (no generic-rank assumption needed).
- Generic affine orbit dimension 18, projective 17. Reported triple: (projective normal-form dim, projective saturation dim, quotient dim = saturation − 17). Ambient quotient: `235 − 17 = 218`.

### 1.3 The all-target table (every target, stated pole orders)  **PROVEN (dimensions) / cites other chapters (constructions)**

| Family (exact hypothesis on F, smooth) | Construction (where) | q | normal form | saturation | quotient |
|---|---|---|---|---|---|
| Coordinatewise even: `F = f(X0²,X1²,X2²,Y0²,Y1²,Y2²)`, `f` cubic (56 coeffs) | constructive **strict** cubic-surface fields [Ch.2 §§2.4–2.5] | every `q>0` | 55 | 67 | 50 |
| Separately even in each of `Y0,Y1,Y2` (or `X/Y` mirror) (119 coeffs = 28+45+36+10) | finite covering-curve construction [Ch.3] | every `q>0` | 118 | 124 | 107 |
| Split `f(X)+g(Y)`, arbitrary ternary sextics (56 coeffs) | finite elliptic-companion construction [Ch.3a] | every `q>0` | 55 | 55 | 38 |
| Crossed `p(X0,X1,Y0)+r(X2,Y1,Y2)`, globally `Y`-even (16+16 coeffs) | finite elliptic-companion construction [Ch.3a] | every `q>0` | 31 | 39 | 22 |
| Polarized split `f(X+Y)+f(X−Y)` with `f = g6(z0,z1) + z2³ g3(z0,z1) + a z2^6` (7+4+1 = 12 coeffs) | K3-correspondence existence [Ch.3a §3A.6] | every `q>0` | 11 | 24 | 7 |
| Binary blocks `F = P0(X0,X1) + P1(Y0,Y1) + P2(X2²,Y2²)` (7+7+4 = 18 coeffs) | explicit binary-block construction [Ch.2a] | 1, 2 | 17 | 25 | 8 |

Proof data for saturations (all CHECKED): coordinatewise-even: torus dim 6, all 12 off-diagonal Fermat monomials outside ⇒ `56+12 = 68` affine. Three `Y`-reflections: preserving group `GL3(X) × (Y torus)` dim 12, 6 Fermat monomials outside ⇒ `119 + 6 = 125`. Split: already `G`-stable. Crossed: `GL2(X0,X1)×Gm(X2)×Gm(Y0)×GL2(Y1,Y2)` dim 10, 8 outside ⇒ 40. Binary: dim-10 subgroup, 8 outside ⇒ 26. Polarized: with `U = X+Y`, `V = X−Y`, block-linear change `δU = CU + EV`, `δV = CV + EU`, `δF = 6 Σ_i [U_i^5 (CU+EV)_i + V_i^5 (CV+EU)_i]`; the 9 `E` directions miss the family tangent, exactly 5 `C` directions (`i,j ∈ {0,1}` and `i=j=2`) meet it ⇒ `12 + 18 − 5 = 25` affine; matching upper bound from `GL2(z0,z1) × Gm(z2)` (dim 5). Witnesses smooth (polarized reduces to Fermat by the invertible `U,V` change).

**Records (chapter's own summary).** Largest all-target, all-`q` families: `(118,124,107)` for finite constructive *legal* certificates; `(55,67,50)` for constructive *strict* certificates. "Finite construction" = terminating algebraic operations, not necessarily practical. Explicitly: "These are dimensions of established all-target families, not upper bounds on the full solvable locus." The §2.7–2.9 cover/Jacobian constructions do not enlarge these records; §4.6's chain application proves all targets for one fixed sextic; the §4.2 surface quotient supports one sector `(22,74,22)` on a 26-coefficient family.

---

## 2. §4.2 Scope of the partial-target constructions

### 2.1 Two-reflection theorem [Ch.2]  **sector theorem**
- Hypothesis: `F` even in one selected `X_i` and one selected `Y_j` (plus globally `Y`-even): 104 coefficients (CHECKED for `X0,Y0`).
- Serves: all targets odd in both selected coordinates: 34 raw `q=1` monomials; exactly 2 odd Jacobian directions in that character (`X_i F_{Y_j}`, `Y_j F_{X_i}`) ⇒ primitive sector dim 32 (CHECKED). Count also = coefficients of `Π_i (1 − t ε_i)^{-1}` in the group algebra of the three sign characters.
- Boundary: "This is a sector theorem, not all-target coverage on that coefficient family."

### 2.2 Conic construction [Ch.2 §2.6]  **bounded scope**
- Weighted sextic `H` in `P(2,2,2,2,1,1)`: 77 coefficients (CHECKED). Four arbitrary quadrics (4·21) + two linears (2·6) ⇒ coefficient image has affine dim ≤ `77 + 84 + 12 = 173` ⇒ proper intersection with `S_6^+` (dim 236).
- Serves the specified Jacobian-density targets, with the denominator bounds of Ch.2; "it does not assert that those densities fill the target space."

### 2.3 Weighted-Fano theorem [Ch.3 §3.4]
- Serves `A = det(dΦ)·b(Φ)` for nonlinear polynomial Piola maps `Φ`; covers every primitive `q=1` nonlinear polynomial Piola atom under its hypotheses, removes the spectral condition of §2.1; the *density condition* remains. "A theorem that transfers all cohomology from a smaller model need not transfer all cohomology of the covering fourfold."

### 2.4 Polarized split, arbitrary smooth ternary `f` [Ch.3a]
- Supported: rank-1073 invariant substructure (incl. ambient class), primitive Hodge numbers `(189,694,189)`; its `q=1` sector contains every odd `A` with `A(X,X) = 0`.
- Remainder: rank 210, Hodge `(19,172,19)`, = exterior square of primitive K3 cohomology (K3 = double plane `w² = f`). Ch.3a proves the decomposition and a *sufficient* algebraic Kuga–Satake criterion, applied to the order-three family. **OPEN**: "No unconditional support theorem for that remaining piece for arbitrary `f` is asserted."

### 2.5 Paired-sign sector theorem (surface quotient)  **PROVEN (geometric existence); certificates existential**

Group data.
- `P = ⟨s_0, s_1, s_2⟩`, `s_i: (X_i,Y_i) ↦ (−X_i,−Y_i)`; projectively `P ≅ (Z/2)²` (`s_0s_1s_2` scalar).
- `h = diag(ζ², ζ, 1, ζ³, ζ², ζ)`, `ζ^6 = 1`, i.e. weights `w = (2,1,0,3,2,1) mod 6`. `P ∩ ⟨h⟩ = 1`, `τ = s_1 h³` (CHECKED: `h³ = diag(1,−1,1,−1,1,−1)`). `Γ = P × ⟨h⟩`, order 24, contains `τ`.
- `Γ`-invariant sextics: the 26-coefficient space of §2.9 (CHECKED): six pure sixth powers, two balanced monomials (`X0X1X2Y0Y1Y2`-type, i.e. `u0u1u2` and `v^6 u0u1u2` in the chart), nine monomials for each edge orientation. In the paired chart the general member is the ternary cubic [§2.9]
  `U = Σ_i (a_i + b_i v^6) u_i³ + (g0 + g1 v^6) u0u1u2 + v²(d0 u0²u1 + d1 u1²u2 + d2 u2²u0) + v⁴(e0 u0u1² + e1 u1u2² + e2 u2u0²)`.

Statement. If `F` is in this space, `X = V(F)` smooth, and the surface `S` (below) is geometrically smooth over the passive field, then `H^4(X,Q)^Γ ⊆ N^1 H^4(X,Q)`. Surface-smoothness holds on a nonempty open of the 26-coefficient space. Supported primitive sector: Hodge `(22,74,22)`, rank 118 (+1 ambient = 119). Every residue in this sector has a certificate at every positive `q` (via Ch.1 prescribed-support theorem / §3.6: denominator `D^q` for `q ≤ 4` where `D` is an equation of the support; `D` is **not explicit** in the chapter).

Quotient field (cover coordinates).
- `[u0:u1:u2] = [X0²:X1²:X2²]`, `v = Y0/X0`, `c_i = (Y_i/X_i)/v` (`i = 1,2`); these generate the `P`-invariant function field modulo `U = 0`.
- `h` fixes `c_1, c_2`; `v ↦ ζ v`; `(u0,u1,u2) ↦ (ζ⁴u0, ζ²u1, u2)`.
- Quotient by `τ`: `b = v²`. Over `K = C(c_1,c_2)`, homogenize `(b:r)`:
  `S = {H = 0} ⊂ P²_u × P¹_{(b:r)}`, `H = Σ_{j=0}^3 b^j r^{3−j} C_j(u)`, bidegree `(3,3)`;
  `C_0, C_3 ∈ span{u0³,u1³,u2³,u0u1u2}`, `C_1 ∈ span{u0²u1, u1²u2, u2²u0}`, `C_2 ∈ span{u0u1², u1u2², u2u0²}`
  (`C_0 = Σ a_i u_i³ + g0 u0u1u2`, `C_3 = Σ b_i u_i³ + g1 u0u1u2`, `C_1` = `d`-terms, `C_2` = `e`-terms).
- Residual order-three action, `ω = ζ²`: `φ: (b:r; u0:u1:u2) ↦ (ωb : r ; ω²u0 : ωu1 : u2)`, preserves `H`. `X/Γ` is birational over the passive surface to the family `T = S/⟨φ⟩`.

Surface cohomology (CHECKED numerics).
- `K_S = O_S(0,1)`, `q(S) = 0` (Künneth on `O(−3,−3)`), `p_g(S) = 2` (canonical sections = residues with numerators `1, b` on `r = 1`). `det` of the action on `u` is 1, `db` has character `ω` ⇒ the two holomorphic 2-forms have characters `ω, ω²`.
- Action is free: a fixed point lies over `b = 0` or `r = 0` at a coordinate vertex of `P²`; none is on `S` because the six pure sixth-power coefficients are nonzero (else the coordinate vertex is a singular point of `X`). Hence `T` smooth, `p_g(T) = q(T) = 0`, `χ(O_T) = 1`, `b_2(T) = 10` (from `χ(O_S) = 3`, `K_S² = 0`, degree-3 free cover; `e(S) = 36`, `b_2(S) = 34`).
- Surface-fibration theorem [§3.2: fourfold with dominant rational map to a surface whose geometric generic fibre has a smooth projective model with `p_g = q = 0`, geometrically connected ⇒ `H^4 = N^1H^4`] + quotient transfer [§3.1] ⇒ the theorem.
- Not rationality: `S → T` étale ⇒ `κ(T) = 1`; genus-one fibration with multiple fibres of order 3 over the images of `b = 0, ∞`; the rational elliptic surface from §2.9's (degree-4, degree-6) Weierstrass coefficients is its Jacobian. "The algebraicity of `H²(T)`, rather than uniruledness of `T`, is the relevant support input."

Openness (explicit witness). Fermat member: `H = Σ_{i=0}^2 (r³ + c_i^6 b³) u_i³`, `c_0 = 1`; choose nonzero `c_i` with distinct sixth powers ⇒ the three coefficient zero-sets on `P¹` are disjoint and reduced ⇒ `S` smooth (at a singular point the `u`-partials force all but one coordinate to vanish; the derivative of the remaining coefficient excludes singularity; fibres at `b = 0, ∞` are smooth diagonal cubics). Properness + openness of smoothness ⇒ nonempty coefficient open. The theorem is stated *with* this condition; smoothness of `X` alone does not imply it.

Character count. Projective volume has `h`-weight `Σw = 9`, so a numerator with exponent vector `e` gives a `Γ`-invariant residue iff all `e_i + e_{i+3}` even and `⟨w,e⟩ ≡ 3 (mod 6)`. Counted in the Fermat ring `C[Z]/(Z_i^5)` (Koszul series constant on the smooth locus):

```python
from itertools import product
w = (2, 1, 0, 3, 2, 1)
def invariant(e):
    return (all((e[i]+e[i+3]) % 2 == 0 for i in range(3))
            and sum(a*b for a, b in zip(e, w)) % 6 == 3)
assert [sum(sum(e) == d and invariant(e) for e in product(range(5), repeat=6))
        for d in (6, 12, 18)] == [22, 74, 22]
```
(CHECKED.) Dimensions of this family: normal form 25, saturation 37 (torus dim 6; 12 off-diagonal Fermat derivatives outside), quotient 20.

Application named: `mixed-polar` ([#case/mixed-polar], solved): its paired surface is smooth over the passive field, its balanced numerator has weight `9 ≡ 3 (mod 6)` ⇒ `Γ`-invariant residue; the instance's covariant construction realizes the existence statement as an algebraic certificate. Attribution: cyclic instance of the surface quotient introduced by `hodgepodgeconjecture` (discussion 55, comment 19077807); the chapter isolates family-wide hypotheses and adds the relative-divisor and finite-descent step via §3.2.

---

## 3. §4.3 Coupled cyclic-seven family  **PROVEN; strict explicit q=1 fields**

Family. `σ` acts on `(x,y,z) = (X0,X1,X2)` with weights `(1,2,4) mod 7`, trivially on `Y`. Complete `σ`-invariant `Y`-even family (50 coefficients, CHECKED):
`F = a x^5 y + b y^5 z + c z^5 x + d x²y²z² + x y³ Q0(Y) + x³ z Q1(Y) + y z³ Q2(Y) + g(Y)`, `Q_i` arbitrary ternary quadratics, `g` arbitrary ternary sextic.

Theorem. Every smooth member has its entire `(τ,σ)`-invariant `H^4` divisor-supported (so certificates for all `σ`-invariant `Y`-odd targets at every positive `q`, existential); all 19 primitive `q=1` directions in this sector admit **strict** cubic-surface constructions.

Cover coordinates.
- `U = x²y/z`, `V = y²z/x`, `W = z²x/y`; `p = aU²W + bV²U + cW²V + dUVW`; `P = p + UV Q0 + UW Q1 + VW Q2 + g`; then `F = P` identically.
- Weighted fourfold `𝒲 = {P = 0} ⊂ P(2,2,2,1,1,1)`.
- Exponent matrix of `(x,y,z) ↦ (U,V,W)`: rows `(2,1,−1), (−1,2,1), (1,−1,2)`, `det = 14` (CHECKED); kernel = `⟨σ, −1⟩` (`−1` = `τ` projectively on a chart `Y_j = 1`) ⇒ degree-14 map identifying `k(𝒲)` with the invariant function field of `X`. `F = 0` is irreducible (reducible hypersurface of this dimension is singular along the intersection of components).

Quasismoothness (stratum-by-stratum, all PROVEN in text): smooth `F` ⇒ `F(X,0)` and `g` smooth plane sextics and `abc ≠ 0`. On `UVW ≠ 0` the coordinate map is étale (Jacobian `14xyz`). `U = 0, VW ≠ 0`: `P_V = W(cW + Q2)`, `P_W = V(2cW + Q2)` cannot both vanish; cyclic for other strata. `U ≠ 0, V = W = 0`: reduces to `g = 0`; `Y ≠ 0` use `∇g ≠ 0`, `Y = 0` use `P_W = aU² ≠ 0`. `U = V = W = 0`: smoothness of `g`.
Ambient singular locus meets `𝒲` along the smooth plane cubic `p = 0` at `Y = 0`; transverse type `½(1,1,1)`, age `3/2`, terminal. Weighted adjunction `−K_𝒲 = O(3)`. Normal-Fano support theorem [Ch.3] ⇒ divisor-supported `H^4` on a resolution; quotient transfer ⇒ `H^4(X)^{τ,σ}` supported.

Strict `q=1` construction (explicit recipe).
- Project to `P²_Y`; generic fibre over `k(P²_Y)`: `p(U,V,W) + T[UV Q0 + UW Q1 + VW Q2] + T³ g = 0`, a smooth cubic surface (affine part smooth by generic smoothness; `T = 0` smooth by smoothness of `p`).
- Normalizations: `T = Y0²`, `α = Y1/Y0`, `β = Y2/Y0`. Full six-coordinate Jacobian of `(U,V,W,T,α,β)` w.r.t. `(x,y,z,Y0,Y1,Y2)`: `28 xyz / Y0`.
- Required numerator character: `σ`-trivial, `τ`-odd. Numerator space: `xyz · k[Y0,Y1,Y2]_3 ⊕ span{x³y², y³z², z³x²} · k[Y]_1`, `10 + 9 = 19` directions (CHECKED by enumeration).
- Cubic-surface densities (divide by the Jacobian): `T² α_i α_j α_k / 28` and `U T α_i/28, V T α_i/28, W T α_i/28` with `α_0 = 1, α_1 = α, α_2 = β`.
- Then Ch.2's cubic-surface ("six-line") theorem [§2.4], normalized trace and full Piola (cofactor) transport [§2.1] give strict fields. Coordinates `Y`-even, Jacobian `Y`-odd ⇒ correct field parity; chart divisors (`x,y,z,Y0`) and line divisors avoid `F`.
- Denominator shape: pullback of the degree-45 tritangent-plane product `T_45 = Π ℓ` of §2.4 (in cubic coordinates, after trace over the field of the 27 lines) times chart monomials; explicit in principle, degree not stated in Ch.4.
- Jacobian-triviality excluded: `Y_j F_{X_i}`, `X_i F_{Y_j}` have `σ`-characters `−w_i, +w_i ≠ 0`.

Equivariant Koszul series (reusable): in `Z[u,v]/(u^7 − 1, v² − 1)`, the coefficient of `t^{6q} u^0 v^1` in
`Π_i (1 − t^5 u^{−w_i} v^{e_i}) / (1 − t u^{w_i} v^{e_i})`, `w = (1,2,4,0,0,0)`, `e = (0,0,0,1,1,1)`,
is 19 at `q=1` and `395 − 276 + 9 = 128` at `q=2` (polynomial, first Koszul, second Koszul terms); Hodge `(19,128,19)` (CHECKED both via series and via Fermat-ring enumeration).
Dimensions: projective family 49; saturation ≤ `49 + 6 = 55` (diagonal `X` changes and all `Y` changes preserve it), equality at the smooth split Klein member (six off-diagonal `X` orbit derivatives have the six distinct nontrivial `σ`-characters; Klein smoothness: on the torus multiplying the three derivative equations gives `1 = (−5)³`). "This remains partial-sector coverage on a family of projective dimension 49, saturated dimension 55."

---

## 4. §4.4 One generic nonzero target suffices  **PROVEN (conditional reduction)**

"Generic" = geometric generic point of the full even-sextic parameter space `P^235`, equivalently a finite-type certificate/support family dominating it. "An isolated certificate, however difficult to find, does not satisfy this hypothesis."

**Theorem.** If a dominating finite-type family of reduced divisor supports on smooth even sextics carries one nonzero invariant primitive `(3,1)` class at the geometric generic point, then its Gysin image carries the entire invariant primitive `H^4` there. In particular a *generic* `q=1` certificate for any one nonzero primitive target (e.g. `X0³Y0³`) suffices for every required target (`q = 1, 2` via Ch.1; all `q` via §4.5/§3.6). If `E_i` are the resolved support components then `Σ_i h^{2,0}(E_i) ≥ 208`.

Inputs: (i) infinitesimal residue variation = multiplication by the deformation polynomial in the Jacobian ring (differentiate `A Ω/F²`: next numerator `−2gA`); (ii) images of morphisms of polarizable VHS are polarizable subvariations, orthogonal complements are subvariations, morphisms are strict; applied fibrewise over a shrunk base of a resolved support family.

Linear-algebra step at Fermat (reusable computation).
- `R = C[Z]/(Z_0^5,…,Z_5^5)`, `V = R_6^-` (dim 208), `G = R_6^+` (dim 218), `W = R_12^-` (dim 866).
- Higgs map `M_g: V → W, u ↦ gu`. Diagonal Fermat automorphisms act by distinct characters on monomials ⇒ Hodge metric diagonal ⇒ for monomial `g = Z^a`, `M_g† M_g` is diagonal with nonzero entry on `Z^u` iff `u + a ≤ (4,…,4)` coordinatewise.
- Three combinatorial facts (CHECKED by running the chapter's script): support signatures distinguish all 208 `u`; every `u` has a common nonzero product with the anchor `b = (1,1,1,1,1,1)`; `G·V = W`.

```python
from itertools import product
def check(extra, anchor, expected):
    def parity(u): return sum(u[3:]) % 2
    box = list(product(range(5), repeat=6))
    G = [u for u in box if sum(u)==6 and parity(u)==0 and extra(u)]
    V = [u for u in box if sum(u)==6 and parity(u)==1 and extra(u)]
    W = {u for u in box if sum(u)==12 and parity(u)==1 and extra(u)}
    signatures, images = [], []
    for u in V:
        image = {tuple(x+y for x,y in zip(u,a)) for a in G
                 if all(x+y<=4 for x,y in zip(u,a))}
        images.append(image)
        signatures.append(tuple(all(x+y<=4 for x,y in zip(u,a)) for a in G))
    assert (len(G),len(V),len(W)) == expected
    assert len(set(signatures)) == len(V)
    anchor_image = images[V.index(anchor)]
    assert all(image & anchor_image for image in images)
    assert set().union(*images) == W
    return expected
assert check(lambda u: True, (1,1,1,1,1,1), (218,208,866)) == (218,208,866)
assert check(lambda u: u[5]%2==0, (1,1,1,1,2,0), (144,126,492)) == (144,126,492)
```
- Distinct signatures ⇒ distinct joint eigenvalue tuples (independent of the actual metric values) ⇒ polynomial interpolation gives all monomial projectors `P_u`. For a common product exponent `w ≥ u, b`, put `a = w − u`, `c = w − b` (allowed even sextic exponents): `P_b M_{Z^c}† M_{Z^a} P_u` is a nonzero matrix unit `u`-line → `b`-line; adjoints and composition give all matrix units ⇒ the `M_h† M_g` generate `End(V)`.

Geometric step: a polarizable subvariation `U` of the invariant primitive variation has `U^{3,1}` stable under all `M_h†M_g` ⇒ `U^{3,1} = 0` or `V`. Using Fermat even if the support family lives on another open: finitely many matrix-unit words have nonzero determinant in a real-analytic neighbourhood of Fermat; the smooth locus of a dominant finite-type map meets it with surjective tangent map. If `U^{3,1} = V`, Higgs stability and `G·V = W` force `U^{2,2} = W`; conjugation gives `(1,3)`. Strictness gives `⊕ H^{2,0}(E_i) ↠ H^{3,1}(X)^τ` ⇒ bound 208; components contracted to codim ≥ 2 contribute no `(3,1)`.

Remark: a certificate over the algebraic closure of the universal coefficient field has finitely many algebraic coefficients ⇒ lies in a finite extension ⇒ denominator and support spread over a nonempty open of a finite-type base change ⇒ theorem applies. `X0³Y0³` is nonzero at Fermat (exponents ≤ 4) hence generically.

### Intermediate one-reflection locus  **CONJECTURE/NOT DONE**
- Retain `τ` and the `Y2`-reflection `r`. Residual sector: `τ = +1`, `r = −1`; numerators `Y`-odd and `Y2`-even. Hodge `(126,492,126)`; allowed deformation space 144-dimensional (`extra = u[5] % 2 == 0`, anchor `(1,1,1,1,2,0)`; CHECKED).
- Same proof verbatim: one dominating support with a nonzero residual `(3,1)` class carries the whole residual sector; the `r`-invariant complement already has Fano support [Ch.3]. Would force `Σ h^{2,0}(E_i) ≥ 126`.
- Coefficient dim 158 (CHECKED), preserving group `GL3(X) × GL2(Y0,Y1) × Gm(Y2)` dim 14, four outside Fermat monomials `Y2^5 Y0, Y2^5 Y1, Y0^5 Y2, Y1^5 Y2` ⇒ `158 + 18 − 14 = 162` affine; conditional projective dims `157 / 161 / 144`.
- Verbatim: "A support with the required nonzero residual image has not been constructed."

---

## 5. §4.5 Raising `q` by differentiating a generic certificate  **PROVEN (identity; hand-CHECKED)**

**Proposition.** `F_t = F + tC`, `C` an even sextic, `A` independent of `t`, `K_q(F_t, A; v_t) = F_t^{q+1} R_t`, where `K_q(F,A;v) = A − q v(F) + F div v`. Then
`w = (q C v − F v̇)/(q+1)` satisfies `K_{q+1}(F, A·C; w) = −F^{q+2} Ṙ/(q+1)`.
Proof identity: `K_{q+1}(F, AC; w) = C·K_q(F,A;v) − F·∂_t K_q(F_t,A;v_t)|_0 /(q+1)`; differentiating `F_t^{q+1}R_t` cancels the `C F^{q+1} R` term.
(My check: `∂_t K_q|_0 = −q v̇(F) − q v(C) + C div v + F div v̇`; expanding both sides, the difference is `F v̇(F) − CF div v + F[(q+1)C div v − (q+1)v̇(F)]/(q+1) = 0`.)
- Degrees/parity: `deg w = 6(q+1) − 5`; `C, F` even ⇒ same parity as `v`. For `v = N/D`, `w` has denominator dividing `D²`; reduced support and coprimality preserved.
- Over the universal coefficient field, `Ḟ = C` is a derivation; extends uniquely to finite separable extensions; apply, then normalized trace [Ch.1] ⇒ certificate over the original field. "An analytic choice of solutions is unnecessary, but a generic family is necessary."
- Monomial factorization lemma: every odd monomial of degree `6q` is a product of one odd sextic and `q−1` even sextics (`Y`-degree `s`, pick `b ∈ {1,3,5}`, `b ≡ s mod 6`, take `b` `Y`-factors and `6−b` `X`-factors; group the rest into all-`Y` and all-`X` sextics). Iterating + linear combinations ⇒ all positive `q` from generic `q=1`.
- Boundary: "it cannot be applied to an isolated successful tuple without such an extension."

---

## 6. §4.6 Descent through a deck action and an exchange  **PROVEN**

**Proposition.** `p: Z → X` proper dominant generically finite between smooth complex projective varieties; `d, g` finite-order automorphisms of `Z, X`; `r, σ` commuting involutions of `Z`; `τ` involution of `X`; `e = σ d σ`; suppose `p d = p`, `p e = g p`, `p r = τ p`. If smooth projective models of `Z/⟨r,σ⟩` and `X/⟨g⟩` have divisor-supported `H^k`, then `H^k(X,Q)^τ ⊆ N^1 H^k(X,Q)`.
Proof sketch: `a = p^*ξ` in `V = H^k(Z)/N^1`; `d^*a = r^*a = a`; `p^*ξ + σ^*p^*ξ` is `⟨r,σ⟩`-invariant hence 0 in `V` ⇒ `σ^*a = −a` ⇒ `e^*a = a`; order of `g` divides order of `e` (dominance) ⇒ `a = Av_e a = p^*(Av_g ξ) = 0` in `V`; `p_* p^* = deg p`.
Requirement: "apply the criterion after exhibiting a finite symmetry group and resolving its action together with the cover projection equivariantly. Finiteness is a hypothesis to check; individual birational symmetries do not by themselves supply the required models."

Chain application (the `holdout` sextic, [#case/holdout]; data from §2.7 and holdout/solution.md):
- `F_H = Σ_i (X_i^6 + Y_i^6) + X0³ X1 Y0 Y1 + X1³ X2 Y1 Y2` (smooth by the max-modulus estimate `6M^5 > 4M^5`).
- Cover: `T^6 = Y0^6 + Y2^6`, `α = Y0/T`, `β = Y2/T`, `y = Y1`; `K_b = Frac(C[α,β]/(α^6+β^6−1))`; `S = y² + T²`, `P = yT`; quotient `H = Σ X_i^6 + P(α X0³X1 + β X1³X2) + S³ − 3SP²`, weights `(X0,X1,X2,S,P) = (1,1,1,2,2)`.
- `r: (y,T) ↦ (−y,−T)`, `σ: (y,T) ↦ (T,y)`; `d: (y,T,α,β) ↦ (y, ζT, ζ^{-1}α, ζ^{-1}β)`, `e = σdσ: (ζy, T, ζ^{-1}α, ζ^{-1}β)`, `ζ^6 = 1`; `r = d³e³`; group generated is finite (quotient of `(μ6×μ6) ⋊ Z/2`). Physical `g` on `X`: diagonal with weights `(0,0,0,5,1,5) mod 6`; `g` fixes `F_H`; fixed point with only `Y0,Y2 ≠ 0`, `Y0^6 + Y2^6 = 0`, tangent weights `(1,1,1,2)/6`, age `5/6` ⇒ finite-quotient theorem [§3.3] supports `H^4(X)^g`. Cover quotient fibres: degree-6 hypersurfaces in `P(1,1,1,2,2)`, quasismooth, three `½(1,1,1)` points, `−K = O(1)` ⇒ normal Q-Fano threefolds ⇒ uniruled total space ⇒ support.
- Conclusion (Ch.4 §4.6 + §3.6): **every odd target at every positive pole order for `F_H`** (existential). Character bookkeeping [holdout/solution.md]: for `A = M(X) Y0^a Y1^b Y2^c`, `e = b − a − c − 1`; cover criterion serves `e ≢ 0 (mod 6)` (158 of 226 raw monomials, 149 of 208 basis), `g`-quotient the complement (68 / 59). Attribution: consequence of `fqlx`'s cover construction.

---

## 7. §4.7 Scope of the support mechanisms  **PROVEN obstructions (for the named mechanisms)**

Three conditions for any construction: (1) auxiliary geometry supplies divisor support; (2) the correspondence reaches the required representation (averaging ideal of §3.1: a simple `Q[G]`-module `ρ` is served by quotients `M/K_j` iff `ρ^{K_j} ≠ 0` for some `j`); (3) a computed certificate binds the primitive to the exact source, normal derivatives and legal denominator (polynomial-density theorem or chosen curve family, only for sources satisfying its hypotheses).
Key negative statements: "Neither an additional incidence seed on the same quotient nor a higher-order integration identity enlarges that quotient's averaging ideal." "Complete coverage on a special variety or symmetry locus does not establish [the §4.4] hypothesis."

### 7.1 Transcendental part of the surface fibres
- Paired signs + global `Y`-evenness: 74 coefficients, `(1/2) Σ_{k0+k1+k2=3} (Π_i (2k_i+1) + 1) = 74` (CHECKED); saturation `74 + 12 = 86` vs 236 ⇒ "No construction confined to these symmetries can dominate the full coefficient space."
- General even sextic `F = Σ_{j=0}^3 F_{6−2j,2j}(X,Y)`. Fix `c_i = (Y_i/X_i)/(Y0/X0)`, `c_0 = 1`, `b = (Y0/X0)²`: surface
  `Σ_{j=0}^3 b^j r^{3−j} F_{6−2j,2j}(X, cX) = 0` in `P²_X × P¹_{(b:r)}`, bidegree `(6,3)`, general member smooth, `K = O(3,1)`, `q = 0`, `p_g = 20` (CHECKED `h^0(O(3,1)) = 20`); coefficient space surjects onto the full `(6,3)` system.
- Paired signs (`u_i = X_i²`) ⇒ bidegree `(3,3)`, `p_g = 2`; extra order-three quotient ⇒ `p_g = 0`. Sequence `20 → 2 → 0`.
- Even on the smooth `(3,3)` surface: `e = 36`, `b_2 = 34`, free quotient `b_2 = 10` ⇒ each nontrivial `Z/3`-eigenspace has dim 12 with Hodge `(1,10,1)` (CHECKED numerics). Hence "the surface-fibration theorem does not extend to a nontrivial character merely by replacing the numerator or increasing its pole order."
- Without the extra diagonal symmetry: paired cubic has coefficients of degree ≤ 3 in `b`; Weierstrass coefficients degrees ≤ 12, 18; minimal model with fundamental bundle `O_{P¹}(n)`, `n = 3` ⇒ `χ(O) = 3`, `K` from `O(1)`, Kodaira dim 1 ⇒ not dominated by rational curves with passive parameters fixed. Minimality criterion: coprime Weierstrass coefficients + exact degree 18 (excludes nonminimality at finite places and at infinity). Ref. Schütt–Shioda §§6.6, 8.2.

### 7.2 Natural cubic cover of a general even sextic and its standard representation
- `B = P²_X × P²_Y`, `(a,b,c,d) = (F_{6,0}, F_{4,2}, F_{2,4}, F_{0,6})`. Along the line `[λX : μY]` the quotient `Q = X/τ` has equation `a + bt + ct² + dt³ = 0`, `t = (μ/λ)²`.
- Discriminant, section of `O_B(12,12)`: `Δ = b²c² − 4ac³ − 4b³d − 27a²d² + 18abcd`.
- General `F`: monodromy `S_3` (connected, simple ramification). Over the étale locus: `Q ⊕ V_std` (rank 2 standard rep). Constant summand contributes weight-4 classes from `H^4(B)` ⇒ algebraic; the primitive problem needs `V_std`. Discriminant double cover carries sign; `V_std^{A_3} = 0`; pulling to the double cover splits `V_std` into two nontrivial cubic characters, not constant.
- `Q` not uniruled: `τ`-fixed locus = two smooth plane sextics `X ∩ {X=0}`, `X ∩ {Y=0}`, eigenvalues `(1,−1,−1,−1)`, age `3/2`, terminal `½(1,1,1)`, torsion `K` ⇒ Kodaira dim 0 (Kollár–Larsen Thm 2).
- Subgroup table for a smooth equivariant model `Z` of the `S_3`-Galois closure: `K = 1`: `dim V_std^K = 2`, `Z` (dominates `Q`, not uniruled); `C_2`: 1, model of `Q` (not uniruled); `A_3`: 0, discriminant double cover; `S_3`: 0, rational base `B`.
- Conclusion (verbatim): "no collection of uniruled intermediate quotients of this natural Galois closure can serve `V_std` through the averaging ideal of §3.1. This excludes a specific extension of the Fano strategy, not divisor support itself." A general construction must supply support for the transcendental `V_std` directly or via a larger correspondence.

---

## 8. Reusable formulas (computational inventory)

Inside Chapter 4:
- **R1** Dimension counts: `dim S_6^{(j)} = C(8−j,2)C(j+2,2)`; `H(t) = (1+…+t⁴)^6`; `Θ(t) = (1+t²+…+t⁸)^3`; `dim(R_F)_d^± = ([t^d]H ± [t^d]Θ)/2`.
- **R2** Saturation bound `dim G·V ≤ dim V + 18 − dim H`; Fermat test monomials `Z_i^5 Z_j` (same block); quotient = saturation − 17.
- **R3** Paired-sign quotient coordinates `u_i = X_i²`, `v = Y0/X0`, `c_i = (Y_i/X_i)/v`, `b = v²`; `h`-action `v ↦ ζv`, `u ↦ (ζ⁴u0, ζ²u1, u2)`; `(3,3)` surface `H = Σ b^j r^{3−j} C_j(u)`; `φ = (ωb : r; ω²u0 : ωu1 : u2)`.
- **R4** Fermat `(3,3)` witness `H = Σ (r³ + c_i^6 b³) u_i³` with the distinct-sixth-power smoothness recipe (a template for verifying surface smoothness of any member: check `u`-partials at coordinate vertices and the fibres over `b = 0, ∞`).
- **R5** `Γ`-sector membership test: paired parities + `⟨(2,1,0,3,2,1), e⟩ ≡ 3 (mod 6)` (script above); residue character of a monomial `ζ^{⟨w,e⟩ + 9}` [§2.9].
- **R6** Cyclic-seven coordinates `U = x²y/z, V = y²z/x, W = z²x/y`; `p, P`; exponent matrix `det 14`; `T = Y0², α = Y1/Y0, β = Y2/Y0`; Jacobians `14xyz` (3 × 3) and `28xyz/Y0` (6 × 6); the 19 numerators and their cubic-surface densities `T²α_iα_jα_k/28`, `UTα_i/28`, `VTα_i/28`, `WTα_i/28`; equivariant Koszul series `Π_i (1 − t^5 u^{−w_i} v^{e_i})/(1 − t u^{w_i} v^{e_i})`.
- **R7** Fermat Higgs-map generation test (script above) — template for any symmetry-restricted sector: filter `G, V, W` by a character, check signature separation, anchor connectivity, `G·V = W`.
- **R8** `q`-raising: `w = (qCv − F v̇)/(q+1)` serves `A·C` at pole order `q+1`, with `K_{q+1}(F,AC;w) = C K_q(F,A;v) − F ∂_t K_q(F_t,A;v_t)|_0/(q+1)`.
- **R9** Descent criterion relations `pd = p`, `pe = gp`, `pr = τp`, `e = σdσ`; chain data: `d, e, g = diag weights (0,0,0,5,1,5) mod 6`, age-`5/6` fixed point.
- **R10** General `(6,3)` surface `Σ b^j r^{3−j} F_{6−2j,2j}(X, cX)`; natural cubic cover `a + bt + ct² + dt³`, `Δ = b²c² − 4ac³ − 4b³d − 27a²d² + 18abcd ∈ O_B(12,12)`.
- **R11** Paired-sign coefficient count `(1/2)Σ_{k0+k1+k2=3}(Π(2k_i+1) + 1) = 74`.

Referenced by Chapter 4 but defined elsewhere (pointers; formulas copied for convenience):
- **R12 Cofactor/adjugate transport [§2.1]** — the "adjugate E1 lift": for `Φ` with `J = dΦ`, `Δ = det J`, `V = adj(J)·w(Φ)` satisfies `V(G(Φ)) = Δ·w(G)(Φ)` and `div V = Δ·(div w)(Φ)`. Spectral constructor: `Φ_i` of weights `m_i`, `G` weighted degree 6, `w` weight-0 derivation, density `b` of weight `12 − Σm_i`, with `D0 = w(G)`, `w(D0) = (μ+ν)D0 − μν G`, `div(bw) = (μ+ν)b` ⇒ `F = G(Φ)`, `A = b(Φ)Δ`, `D = D0(Φ)`, `N = adj(dΦ)(bw)(Φ)`, `N(F) = AD`, `D div N − N(D) = μν A F` (a `q=1` certificate with `deg D = 6`, `deg N = 7`). Nilpotent branch: `w²(G) = 0`, `div(bw) = 0` ⇒ `div N = N(D) = 0`.
- **R13 Cubic-surface ("six-line") theorem [§2.4]**: smooth cubic `G` in 4 variables over char-0 `K`; for every `q ≥ 1` and `deg B = 3q − 1` there is `w_q` with `q w_q(G) = B`, `div w_q = 0`, components of degree `3q − 3`; one common active denominator `T = Π ℓ` over the 45 tritangent planes (normalized `ℓ(p) = 1`), degree 45, coprime to `G`. Line field: `G = ℓ1 Q1 + ℓ2 Q2`, `Φ_L = (ℓ1, Q1, ℓ2, Q2)`, `J_L = det dΦ_L`, `w_L = adj(dΦ_L) e_2/ℓ1`, `w_L(G) = J_L`, `div w_L = 0`; six skew-line `J_L` + four `G_i` span the 10 quadrics. Trace descent `w̄ = Tr_{E/K}(w)/[E:K]`, computed as `tr(adj(M_D)M_N)/([E:K] det M_D)`. Higher `q`: `B = Σ b_i G_i`, `c = div b`, `v_q = [b + G w − c E_4/(3(q−1))]/q`.
- **R14 Jacobian-of-cubic normalizations [§2.9]** (used by the paired-sign family): cyclic sextic `F_C = Σ(X_i^6+Y_i^6) + ½Σ_{i mod 3} X_i³X_{i+1}Y_iY_{i+1}`; chart `r = X2²`, `x = X0²/r`, `y = X1²/r`, `z = 1`; `U = ax³+by³+cz³+dx²y+ey²z+fz²x` with `(a,b,c) = (1+v^6, 1+c_1^6v^6, 1+c_2^6v^6)`, `(d,e,f) = ½v²(c_1, c_1c_2, c_2)`; `s_4 = ae²f + bdf² + cd²e`, `s_6 = 27a²b²c² + 4a²ce³ + 4ab²f³ + 6abcdef + 4bc²d³ − d²e²f²`; covariant syzygy on `U = 0`: `J² = 4Θ³ + 4s_4 Θ H⁴ − s_6 H⁶`; `ξ = Θ/H²`, `η = −J/(2H³)`, `η² = ξ³ + s_4ξ − s_6/4`, `D_yξ/(2η) = 3/U_x`; `T = rv^6`, `X = r²ξ`, `Y = r³η`; `V = Y² − X³ − (1/8)a_8 X − (1/256)b_256 = 0`, `a_8 = 8r⁴s_4`, `b_256 = −64r^6 s_6` (degrees 4, 6 in `(r,T)`); polynomial curve degrees `(2,4,6)` in `(r,X,Y)`; seed over `F_5`, `(c_1,c_2) = (2,3)`: `r = 1+2T+4T²`, `X = 3+T+3T²+4T³+4T⁴`, `Y = 2+3T+3T²+4T³+4T⁴+T⁵+T⁶`; source `P_J = −V_α/(144r)`; residue-volume multiplier `F_{X0} det M/(2Y) = −144 r v³`. General family: `I_4 = s_4`, `I_6 = −s_6/4`, `r^d I_d(U(v))|_{v^6 = T/r}` homogeneous of degree `d` in `(r,T)`; sources `P_{d,M,j} = r^{d−1} v^{j−3} ∂I_d/∂k_M (U(v))`, Jacobian source `(c/144)(P_{4,M,j} X + P_{6,M,j})`.
- **R15 Computed support → all `q` [§3.6]**: `δ = D^q` (`1 ≤ q ≤ 4`), solve the complete linear system `q N(F) = Aδ`, `δ div N − N(δ) = F^q R`, `deg N_i = q(deg D + 6) − 5`; for `q ≥ 5` use Jacobian socle degree 24: `A = Σ B_i F_i`, `c = div B`, `v = (B + F w)/q` with `w` serving `(q−1, c)`, `K_q(A; v) = (F/q) K_{q−1}(c; w)`.
- **R16 Koszul correction map [§5.2]** (the "E1 lift freedom"): all `E1` lifts are `N0 + K`, `K = Σ_{i<j} c_ij (F_j ∂_i − F_i ∂_j)`, `deg c_ij = d + 6q − 10`; `T_D(c) = B_D(K) = Σ_{i<j} ( D[(∂_i c_ij)F_j − (∂_j c_ij)F_i] − c_ij (D_i F_j − D_j F_i) )`; `E2 ⇔ −B_D(N0) ∈ F^q S_{2d−6} + im T_D` inside `S_{2d+6q−6}`.
- **R17 Ch.1 operations**: addition `D_{A+B} = DE`, `N_{A+B} = EN + DM`, `Q_{A+B} = E²Q + D²R`; character averaging `(T_γ v)_i = μ_γ (v_i∘γ)/(ν_γ γ_i)`; trace descent `v_Q = Tr_{K/Q}(v)/[K:Q]`; gauge `(D,N) ↦ (D + tF, N + tAε/6)` at `q=1`, `deg D = 6`.
- The "binary-cubic Hessian construction" lives in §2.5 (lines 258–315 of `02-explicit-constructions.md`); not read in full here — pointer only.

---

## 9. Stated limits and boundaries (chapter's own words, condensed)

1. §4.1: all-target records `(118,124,107)` legal / `(55,67,50)` strict are "dimensions of established all-target families, not upper bounds on the full solvable locus."
2. §4.2 two-reflection: sector theorem (32 primitive dims), not all-target. Conic: coefficient image ≤ 173, densities don't fill the target space. Weighted-Fano: density condition remains. Polarized split: rank-210 K3-exterior-square remainder has no unconditional support theorem for arbitrary `f`. Paired-sign: "leaves the other representations of the involution-invariant cohomology untreated"; smoothness of `S` is an extra open condition, not implied by smoothness of `X`; Hodge dims `(22,74,22)` vs `(208,866,208)`.
3. §4.3: partial-sector coverage (19 of 208 `(3,1)` directions), family of projective dim 49, saturation 55.
4. §4.4: needs a *dominating* family; "An isolated certificate, however difficult to find, does not satisfy this hypothesis." Any support for a nonzero `(3,1)` class needs `Σ h^{2,0}(E_i) ≥ 208` (126 for the one-reflection residual sector) — so supports must have nonnormal or non-rational-singularity components [§3.1]. One-reflection residual support "has not been constructed."
5. §4.5: differentiates a coefficient family; inapplicable to an isolated tuple.
6. §4.6: finiteness of the symmetry group and an equivariant resolution are hypotheses to verify.
7. §4.7: paired-sign symmetric families saturate to dim 86 ≪ 236; `p_g` sequence `20 → 2 → 0`; nontrivial `Z/3` characters on the `(3,3)` surface carry `(1,10,1)` — not reachable by changing numerator or pole order; without the diagonal symmetry the genus-one surface has `κ = 1`; for a general even sextic the needed local system is `V_std` of the `S_3` cubic cover, and no uniruled intermediate quotient of the Galois closure serves it; `X/τ` itself has `κ = 0` (not uniruled). "This excludes a specific extension of the Fano strategy, not divisor support itself."

---

## 10. Where this chapter says the open frontier instances fail

**Instances named by the chapter itself** (only two, both solved): `mixed-polar` (`#case/mixed-polar`, §4.2, as an application of the paired-sign sector theorem) and `holdout` (`#case/holdout`, §4.6, chain-sextic application of the deck/exchange descent). No open instance is named in Chapter 4. The chapter's frontier statements are structural; the verbatim ones are:

- §4.2: "This is a sector theorem, not all-target coverage on that coefficient family." / "No unconditional support theorem for that remaining piece for arbitrary `f` is asserted." / "It leaves the other representations of the involution-invariant cohomology untreated."
- §4.3: "This remains partial-sector coverage on a family of projective dimension 49, saturated dimension 55."
- §4.4: "An isolated certificate, however difficult to find, does not satisfy this hypothesis." / "A support with the required nonzero residual image has not been constructed."
- §4.5: "This differentiates a coefficient family; it cannot be applied to an isolated successful tuple without such an extension."
- §4.7: "Complete coverage on a special variety or symmetry locus does not establish that hypothesis." / "No construction confined to these symmetries can dominate the full coefficient space." / "Thus the surface-fibration theorem does not extend to a nontrivial character merely by replacing the numerator or increasing its pole order." / "Consequently no collection of uniruled intermediate quotients of this natural Galois closure can serve `V_std` through the averaging ideal of §3.1."

**My cross-reference to the catalogue's open instances** (`data/catalogue.json`; status = no `proof.json`/`Solution.lean` in `data/<id>/`, checked 2026-11-01). The mapping below is mine, based on each instance's `why.md`; the instances marked * cite Chapter 4 explicitly.

| Open instance | `q`, `A` | Which Chapter-4 boundary applies |
|---|---|---|
| `cycle-q2` * | 2, `X0²X1²X2³Y0²Y1²Y2` on the cyclic three-link `F_C` (a member of the 26-coefficient `Γ`-family) | Residue character `ζ6^{17+9} = ζ6² ≠ 1` ⇒ outside the §4.2 paired-sign sector; exactly the nontrivial order-three eigenspaces of §4.7 (`(1,10,1)` each, holomorphic forms survive) — "does not extend to a nontrivial character merely by replacing the numerator or increasing its pole order." |
| `klein-coupled-character` * | 1, `X0³Y0³` on a member of the §4.3 cyclic-seven family | `A` has `σ`-character `ζ7³ ≠ 1`; §4.3 supports only the `σ`-trivial sector (19 directions). |
| `one-reflection-skew` | 1, `X0³Y0³`, `F` even in `Y2` | Target is `Y2`-even ⇒ lies in the §4.4 one-reflection *residual* sector `(126,492,126)`, whose required dominating support "has not been constructed"; the §2.7 cover does not survive the skew term. |
| `polarized-k3-remainder` | 1, `X0X1X2Y0Y1Y2`, `F = ½(f(U)+f(V))` | Target reaches the rank-210 K3 exterior-square remainder of §4.2 — "No unconditional support theorem … asserted"; `f` is not in the order-three family of the §4.1 table. |
| `off-paired` | 1, `X0³Y0³` | Breaks two of the three paired signs ⇒ leaves the 26-coefficient `Γ`-family (§4.2) and the §2.7 cover; §4.7's `p_g` sequence (20 → 2 → 0) shows which quotient geometry is lost. |
| `quadric-contact-residual`, `reflection-coupled-residual`, `cyclic-triple-residual` | 1, `X0³Y0³` on sextics with a *solved* density companion | Complementary-character companions: §4.7 "changing the numerator can change its representation even when the equation and pole order are fixed. Neither an additional incidence seed on the same quotient nor a higher-order integration identity enlarges that quotient's averaging ideal." |
| `symmetric-sextic-residual`, `diagonal-permutation` | 1, `X0X1X2Y0Y1Y2` (sign ⊗ sign of `S3 × S3`) and `X0³Y0³ − X1³Y1³` (standard rep of `S3`) | Representation outside the averaging ideal of the invariant quotients (§3.1 criterion, as used in §4.7); the standard-representation obstruction of §4.7 is the same shape. |
| `rank-one-biquartic`, `four-plus-two-residual` | 1, `X0³Y0³` | Controlled departures from the split row (`f(X)+g(Y)`) of the §4.1 table / reflection sectors; no all-target family of §4.1 contains them. |
| `quadratic-cubic-pullback` | 1, `X0³Y0³`, `F = Σ Q_i³` | §4.2 weighted-Fano serves the determinant density `det(∂Q_i/∂Z_j)` only; the target is outside that density sector. |
| `dense-rational`, `short-power-presentation`, `eight-plane-configuration` | 1, `X0³Y0³`, dense (all 236 monomials) | No symmetry family applies; the relevant chapter-4 analysis is §4.7 (natural `S3` cubic cover, `V_std` not served by uniruled quotients; `X/τ` has `κ = 0`) and §4.4 (an isolated certificate would not by itself give a dominating family). |

---

## 11. Separation of status

- **PROVEN in the chapter**: all dimension counts and saturations (§4.1); two-reflection/conic/weighted-Fano/polarized scope statements (§4.2, citing Ch.2/3/3a); paired-sign sector theorem incl. surface geometry, freeness, openness (§4.2); cyclic-seven support + strict `q=1` sector (§4.3); generic-target reduction and one-reflection variant (§4.4); `q`-raising proposition (§4.5); deck/exchange descent (§4.6); `(6,3)`/`(3,3)` surface genera, `κ(T) = 1`, `κ(Q) = 0`, `S3`-cover representation analysis (§4.7). External inputs used without proof: Griffiths residues, polarizable VHS strictness, Kollár–Larsen, Schütt–Shioda minimal Weierstrass models, Bertini/adjunction.
- **EVIDENCE (finite computations, reproduced)**: the two embedded Python checks; `(19,128,19)`; all coefficient counts. See §12.
- **CONJECTURE / missing premise**: existence of a dominating support family with a nonzero generic `(3,1)` class (§4.4); its one-reflection residual analogue; unconditional support of the K3 remainder for arbitrary `f` (§4.2).
- **NOT DONE / explicitly open**: one-reflection residual support (§4.4); support for `V_std` of the natural cubic cover of a general even sextic (§4.7); any construction dominating `P^235`.
- **Not in this chapter**: explicit denominators for the existential theorems (§4.2, §4.3 all-`q`, §4.6); degree bounds on supports (Ch.1 remark: no bound on `e = deg D_red` is supplied); the Hessian/binary-cubic construction (§2.5).

## 12. Verification log

Script: `/tmp/hodge/notes/check-04-counts.py` (plain Python, `python3 -I`, ~1 s; no CAS needed). Output:
```
dim S6^+ = 236  dim S6^- = 226
d=6: [H]=426 [Theta]=10  R^+=218  R^-=208
d=12: [H]=1751 [Theta]=19  R^+=885  R^-=866
d=18: [H]=426 [Theta]=10  R^+=218  R^-=208
Fermat Y-odd counts d=6,12,18: [208, 866, 208]
coordinatewise-even coeffs: 56  three-Y-reflection coeffs: 119 (both re-derived by enumeration)
one-reflection (even in Y2) coeffs: 158
two-reflection (even in X0,Y0) coeffs: 104 ; raw q=1 targets 34 ; primitive 32
paired-sign coeffs: 74 (enumeration and closed form) ; Gamma-invariant coeffs: 26 ; cyclic-7 family coeffs: 50
weighted sextic coeffs P(2,2,2,2,1,1): 77
paired-sign sector Hodge dims: [22, 74, 22]
cyclic-7 dims (Fermat ring and Koszul series): [19, 128, 19] ; exponent matrix det 14 ; 19 numerators = 10 (Y-deg 3) + 9 (Y-deg 1)
4.4 checks: (218, 208, 866) and (144, 126, 492) pass
(6,3) surface p_g = 20 ; (3,3) surface p_g = 2, chi(O) = 3, e = 36, b2 = 34 ; free Z/3 quotient e = 12, b2 = 10
```
Hand-verified: `τ = s_1 h³`; the `h`-weights of `u_i, v, c_i`; `C_j` spans from the §2.9 cubic; the §4.5 identity (expansion in §5 above); degree/parity/denominator bookkeeping in §4.5. Not re-derived: the geometric theorems themselves (support, quasismoothness strata, VHS arguments), taken as stated.
