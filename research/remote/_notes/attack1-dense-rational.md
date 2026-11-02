# Attack 1 — `dense-rational` (q = 1, A = X0^3 Y0^3, F = Fermat + (small integers)/9973 on all 236 even monomials)

**Status: NO WITNESS FOUND.** No certificate, no Lean claim. This note records (1) the obstruction
picture and the geometric mechanism any certificate must realise, (2) every experiment run (exact,
mod p = 2147483647 unless stated; engine v2 after the 2026-11-01 sign fix), (3) —, (4) the
mathematics developed, (5) ranked next steps.  Raw log: `/tmp/hodge/results/dense-rational/log.jsonl`;
scripts in `/tmp/hodge/results/dense-rational/`.

(Sections are filled in below as results arrive.)

## 1. The obstruction and the geometric mechanism a certificate must realise

### 1.1 Set-up and what is known
* `F = F0 + (1/9973) F1`, `F0 = Σ Z_i^6` (Fermat), `F1` has integer coefficients 1..17 on the 230 non-pure even
  monomials. `A = X0^3 Y0^3`, `q = 1`. Variables `(X0,X1,X2,Y0,Y1,Y2)`.
* `A ∉ J(F)` (checked: rank 18 vs 19), so the class `[A] ∈ R_6^- ` (dim 208) is a nonzero invariant (3,1)-class.
  At Fermat `X0^3Y0^3` has all exponents ≤ 4, hence is nonzero there too; the dense perturbation keeps it nonzero.
* No hidden linear structure: catalecticant ranks of `F` are the generic `(6, 21, 56, 21, 6)` (a linear change of a
  split form `P(U)+Q(V)` would give second-partial rank ≤ 12); the stabiliser Lie algebra `{M ∈ gl_6 : M·F = λF}`
  is the Euler line only; full monomial support ⇒ the diagonal symmetry group is `{scalars} × ⟨τ⟩` and the exponent
  lattice `L(F)` has rank 5 (Ch. 1: all Y-odd sextic targets lie in one character coset — no sector decomposition).
* Prior work (fqlx, discussion #59, 2026-10-20): `deg D ≥ 7` for every legal denominator (complex coefficients,
  arbitrary `Q`); every Euler-normalised degree-7 `E1` solution has `D = ∇G·w` (`w` quadratic) — i.e. `D ∈ J(F)` —
  and a necessary condition is the quadratic system `wᵀ Hess(G) w ∈ J_8` (951 equations in 120 unknowns after
  removing the Euler gauge `w ↦ w + ℓ·ε`); 31 numerical starts, two unfinished Gröbner runs; nonlinear homogeneous
  polynomial Piola maps are excluded; the `H_s`-pole module is excluded.

### 1.2 Where E1 can be satisfied: the kernel of `[A]` on the Jacobian ring (new, exact)
`E1: N(F) = A·D` is solvable iff `[A]·[D] = 0` in `R = S/J(F)`. `D ∈ J_d` always works trivially
(`N = A·w` for `D = w(F)`). The new information is the kernel of `[A] : R_d^+ → R_{d+6}^-` (Table in §2.1):
it is **zero for d ≤ 8** (direct Macaulay-matrix nullspace), [d = 9: see table], and for `d ≥ 10` it is as large as
the dimension count allows (`[A]` is *injective* on `R_e^-` for every `e ≤ 8`, so `ker_d = dim R_d^+ − dim R_{18−d}^-`):
`294, 510, 677, 720, 702, 595, 474, 330, 218` for `d = 10..18`.
Consequence: **every certificate denominator of degree ≤ 8 (≤ 9 if the d = 9 entry is 0) lies in `J(F)`, i.e. is a
"polar divisor" `D = w(F) = Σ w_i ∂_iF`** (`w` a polynomial vector field of degree `d − 5`). Its support
`Y_w = {z ∈ X : w(z) ∈ T_zX}` is the locus where the field `w` is tangent to `X`. Degree 10 is the first degree where
non-polar denominators can satisfy `E1` (a 294-dimensional space modulo `J_10^+`, which itself has dim 750).

### 1.3 The geometric mechanism: a (2,0)-budget on the support
Write `Y = X ∩ V(D)` (reduced support; only the support matters by the prescribed-support theorem) with components
`Y_i = X ∩ V(G_i)` and resolutions `π_i : Ỹ_i → Y_i`. A certificate exists iff `res(ω_A)` lies in the Gysin image
`W = Im(⊕ H^2(Ỹ_i, Q) → H^4(X, Q))`, a sub-Hodge structure. Since the Gysin map has type (1,1), the (3,1)-part of `W_C`
is the image of `⊕ H^{2,0}(Ỹ_i)`. Hence, **unconditionally**, some `Ỹ_i` has `h^{2,0} ≥ 1`: the support must have
non-rational singularities (smooth sections, ADE/ordinary-double loci, ordinary triple points, etc. are all excluded).

Let `V = H^4(X,Q)^τ_prim` (rank 1282, Hodge numbers `(208, 866, 208)`). If `V` is a *simple* Hodge structure
(expected for this `F`: it has no symmetry beyond `τ`, and over the full even-sextic family the Higgs maps generate
`End(V)` — Ch. 4 §4.4), then `W ∩ V ≠ 0 ⇒ V ⊂ W`, so
```
      Σ_i h^{2,0}(Ỹ_i) ≥ 208   and   ⊕ H^{2,0}(Ỹ_i) ↠ H^{3,1}(X)^τ          (the (2,0)-budget).
```
Equivalently (Ch. 4 §4.4's bound for a dominating family, here forced for the single instance by simplicity):
any legal `D` for *any* target serves *all* 208 primitive targets. A testable corollary: for every `D` passing `E1`,
`target_image(D)` must have dimension 18 (= `J_6^-`, always served) or 226 (everything); an intermediate value would
disprove simplicity.

### 1.4 Exact formula for the budget and what it excludes
For reduced `Y = X ∩ V(D)` of degree `d` (Cohen–Macaulay, `ω_Y = O_Y(d)` because `K_X = 0`), Grauert–Riemenschneider
and `H^1(Y, ω_Y) = H^2(Y, O_Y)^∨ = 0` give the exact, unconditional formula
```
   h^{2,0}(Ỹ) = h^1(Ỹ, ω_Ỹ) = dim coker( H^0(Y, O_Y(d)) → H^0(Y, ω_Y / π_*ω_Ỹ) ),
   H^0(Y, O_Y(d)) = S_d / (F·S_{d−6} + C·D),   s(d) := dim = C(d+5,5) − C(d−1,5) − 1 .
```
The sheaf `ω_Y/π_*ω_Ỹ` ("adjoint quotient") is supported on the non-rational locus and is computable singularity by
singularity: a (transversally) ordinary `k`-fold point of a threefold contributes jets of order ≤ `k−4`
(`C(k,4)` conditions); a curve `C ⊂ Sing Y` with transversal type "cone over a smooth plane curve of degree `k`"
contributes `⊕_{j=0}^{k−3} Sym^j(N*_{C/X}) ⊗ O_C(d)` (rank `C(k,3)`); an ordinary double surface `Σ` contributes
`O_Σ(d)`; a transversally ordinary `k`-fold surface contributes `⊕_{j≤k−2} Sym^j(N*_{Σ/X})(d)`.
**A certificate requires the adjoint conditions to be dependent (cokernel ≠ 0), and ≥ 208-fold dependent under
simplicity.** Consequences for this `F` (details in §4):

1. *Isolated singular points.* Degree-`d` forms surject onto jets of order `≤ d−4` at a point, so a point of
   multiplicity `m ≤ d` imposes independent conditions: contributes nothing. (Several points can only help if they
   are in special position; the cost/payoff count `C(m+4,5)` vs `C(m,4)` makes this hopeless.)
2. *The τ-fixed plane sextics `C_X = X ∩ {Y=0}`, `C_Y = X ∩ {X=0}` (genus 10).* `N_{C_Y/X} = O(1)^3`, so a
   transversal degree-`k` cone along `C_Y` contributes `Σ_{j=0}^{k−3} C(j+2,2)·h^0(O_{C_Y}(d−j))` conditions
   (`= 540` for `d = k = 7`, `945` for 8, `2376` for 10). The restriction map is *surjective* because the
   `X`-degree-`j` part of a form is `Σ_{|α|=j} X^α f_α(Y)` and `C_Y ⊂ P^2` is projectively normal. Hence
   **ternary cones `D = g(X0,X1,X2)` (vertex `P_Y`) and `D = g(Y)` give `h^{2,0}(Ỹ) = 0`** whenever the rest of
   `Y` has rational singularities — unconditional exclusion at every degree (this is exactly what fails relative to
   Fermat, where the second singular curve `C_X` made the 460 → 560 restriction map non-surjective and gave
   `h^{2,0} = 10·10 = 100`, the (3,3)-block sector). The counts also show that *both* curves together can never
   reach the budget: `D` vanishing to order `k` along `C_Y` and `k'` along `C_X` forces `k + k' ≤ d + 2`
   (bidegree bookkeeping using `F ≡ d(Y) mod (X)^2`, `F ≡ a(X) mod (Y)^2`), and
   `a(k) + a(k') < s(d) + 208` for all `d` (e.g. `d = 6: 90 ≪ 668`; `d = 8: 325 ≪ 1473`; asymptotically `d^4/32 ≪ d^5/120`).
   This kills the natural "structural" denominators `a(X)`, `d(Y)`, `a+b`, `c+d`, `E_X(F)`, `E_Y(F)`, products
   `a·d`, etc. at every degree, and the whole family `D ∈ (X)^k + (Y)^{k'} + (F)`.
3. *Other curves `C ⊂ X`.* Vanishing to order `k` along a curve of degree `δ` costs ≈ `(δd+1−g)·C(k+3,4)`
   conditions on `S_d` and pays ≈ `(δd+1−g)·C(k,3) − (2g−2)C(k,4)` adjoint conditions; the ratio is `≈ 4/k`, so a
   payoff exceeding `s(d) ≈ dim S_d` is impossible for `k ≥ 4`, and `k ≤ 3` gives transversal `p_g ≤ 1`.
   Curves (lines, conics, the 720 τ-invariant lines, plane sections) cannot carry the budget.
4. *Surfaces.* After normalisation `ν : Ȳ → Y`, `h^2(O_Ȳ) = dim ker(H^2(Σ, ν_*O_Ȳ/O_Y) → H^3(Y, O_Y))`; for an
   ordinary double surface `Σ` this is `dim coker(S_d → H^0(Σ, O_Σ(d))) = h^1(P^5, I_Σ(d))`, the failure of
   `d`-normality of `Σ`. Complete-intersection surfaces (`X ∩ P^3`, `X ∩ H ∩ Q`, …) are projectively normal ⇒ 0.
   **The only mechanism left open for a general even sextic is a support that is non-normal along a surface
   `Σ ⊂ X` that fails `d`-normality by ≥ 208 (e.g. an isomorphic projection, or a ruled surface swept by a
   1-parameter family of rational curves on `X`), with `D ∈ I_Σ^{(2)}` of degree `d`.** This is consistent with
   Ch. 4 §4.7 (the needed local system is `V_std` of the `S_3`-cover; no uniruled intermediate quotient serves it).

## 4. Mathematical ideas developed (partial constructions; what would be needed)

### 4.1 Degrees ≤ 8 (≤ 9) are "polar divisors": the E2 problem on `J(F)`
By §1.2 every candidate of degree `d ≤ 8` is `D = w(F)`, `w ∈ (S_{d−5})^6`, modulo the gauge `w ↦ w + ℓ·ε`
(`ε` = Euler field, `deg ℓ = d − 6`), which changes `D` by `6ℓF` and `N` by `Aℓε/6` and preserves E1, E2 mod F.
With `N = A·w + K∇F` (`K` skew, entries of degree `d − 4`) one computes
```
  div N = w(A) + A div w + Σ_j (∂·K)_j F_j ,      N(D) = A·w(D) + Σ_{ij} K_ij F_j D_i ,
  w(D) = w(w(F)) = (w·∇w)·∇F + wᵀ Hess(F) w ,
```
so `B_D(N) = D div N − N(D) ≡ −A·wᵀHess(F)w  (mod J)`; the necessary condition `A·wᵀHess(F)w ∈ J_{d+8}` and the
injectivity of `[A]` on `R_{d+1}^+` (true for `d+1 ≤ 8`, §2.1) give fqlx's `wᵀHess(F)w ∈ J_{d+1}` — i.e.
**`w²(F) ∈ J`**: the field `w` must be "nilpotent on `F` modulo the Jacobian ideal". The same computation for
`d = 8` (`w` cubic, 6·56 − 21 = 315 unknowns) gives `dim R_9^+ = 623` quadratic equations. In both cases the system is
vastly overdetermined; a solution needs a structural reason.

Geometrically `Y_w = X ∩ V(w(F)) = {z ∈ X : w(z) ∈ T_zX}` is the locus where the field `w` is tangent to `X`.
It is singular at `z` iff `∇(w(F)) = (Dw)ᵀ∇F + Hess(F)·w` is proportional to `∇F`. The (2,0)-budget (§1.3–1.4)
demands a non-rational, highly dependent singular locus — for `d = 7, 8` necessarily a surface `Σ ⊂ X` along which
`w(F)` vanishes doubly modulo `F`, with `h^1(I_Σ(d)) ≥ 208`. The natural polar fields fail:
`w = E_X` (or `E_Y`) gives `Y = X ∩ V(E_Y F) = ` the ramification divisor of `X → P²_X × P²_Y`, singular along
`C_X ∪ C_Y` with transversal `A_1` (because `E_Y F ≡ −6a − 4b − 2c` has `(X)`-order exactly 2) — rational;
`w = m(Z)·v` gives the reducible `V(m) ∪ V(∂_vF)` (hyperplane-type section ∪ smooth polar); `w = ℓ·E_X`
gives `V(ℓ) ∪ V(E_XF)`. I found no field `w` of degree 2 or 3 whose tangency locus is non-normal along a surface of
`X`; this is the concrete open problem for degrees 7–8, and the budget says a random search (fqlx's 31 starts, mine)
cannot succeed.

### 4.2 Deformation from Fermat (why.md's suggestion) — framing and the formal/exact gap
At `F0` the degree-6 denominators for `A = X0³Y0³` include the 12 subset sums `D_I = Σ_{i∈I} Z_i^6`, `|I| = 3`,
`I` containing exactly one of `X0, Y0` (strict fields `A E_I/(6 D_I)`), their pencils `D_I + tF0`, and products.
`F_ε = F0 + εF1`, `F1 = 9973(F − F0)` (integer coefficients). Order-`ε` equation for `(D1, N1, Q1)` with the
canonical `N0 = A E_I/6`, `Q0 = 0`:
```
  A D_I D1 + F0[A D1 − N0(D1)] − D_I N1(F0) + F0[D_I div N1 − N1(D_I)] − F0² Q1 = D_I · A·E_I(F1)/6 .
```
This is linear; with the full numerator freedom `N0 ↦ N0 + K0` (Fermat kernel) the order-`ε` problem becomes
bilinear in `(K0, D1)`. Artin approximation makes "formal solution to all orders at degree `d`" equivalent to
"degree-`d` certificates exist for all but finitely many `ε`"; since `deg 6` fails at `ε = 1/9973` (fqlx), the
degree-6 formal deformation must be obstructed at a finite order unless `1/9973` is exceptional. Results in §2.4.
The limit `ε → 0` of a hypothetical family of degree-`d` certificates may degenerate to a *spurious* Fermat point
(`D0 ∈ (F0)`), so a complete deformation analysis must start from the whole Fermat solution variety, not only
from the legal subset sums.

### 4.3 Fermat join versus dense: where the 100 (2,0)-forms came from and why they die
At Fermat `Y = X ∩ V(X0^6+X1^6+X2^6) = join(C_X, C_Y)` (singular along both plane sextics, transversal type
"cone over a plane sextic", transversal `p_g = C(6,3) = 20`). The adjoint quotient has `h^0 = 280 + 280 = 560` while
`h^0(O_Y(6)) = 460`, so `h^{2,0}(Ỹ) ≥ 100`, and `= 100 = g(C_X)·g(C_Y)` (Ỹ is a `P¹`-bundle over `C_X × C_Y`):
exactly the `(3,3)`-bidegree sector `{P_3(X)·Q_3(Y)}` (10·10 targets) that the split field serves.
For dense `F`, `a(X)` restricted to `X` has `(Y)`-order 2 along `C_X` (`a ≡ −b − c − d`, `b ∈ (Y)²`), so the second
singular curve degenerates to an ordinary double curve and the budget collapses to `280 − 280 = 0`. Restoring two
non-rational curves needs `k + k' ≤ d + 2` and never reaches the simple-HS budget (§1.4 item 2). This is the precise
sense in which "the geometry of the Fermat certificate does not survive the perturbation".

### 4.4 Consequences of the budget for search strategy
* All isolated-singularity, curve-singularity, complete-intersection-double-surface and `C_X ∪ C_Y`-type supports
  are excluded at every degree (unconditionally where surjectivity of the adjoint restriction is proved: ternary
  cones with vertex `P_X`/`P_Y`, any plane as vertex, smooth sections; conditionally on simplicity for the rest).
* The polar divisors of degree 7–8 are excluded unless their tangency locus is non-normal along a surface.
* What is needed: a surface `Σ ⊂ X` with `h^1(I_Σ(d)) ≥ 208` for a `d` such that `I_Σ^{(2)}` has a form of degree
  `d` (modulo `F`); or a correspondence-type existence proof (Kuga–Satake / `V_std` of the natural `S_3`-cover,
  Ch. 4 §4.7) which no current construction supplies.
