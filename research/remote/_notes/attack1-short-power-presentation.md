# Attack 1 — `short-power-presentation` (q = 1, A = X0^3 Y0^3, F = eight sixth powers)

**Status: NO WITNESS FOUND (work in progress; resumed after a crash of the first run).**
No certificate, no Lean claim. Engine v2 (post sign-fix). Raw log:
`/tmp/hodge/results/short-power-presentation/log.jsonl`; E1-kernel table `e1kernel.jsonl`; scripts in the same dir.

## 0. Attack plan (written 2026-11-01, resume)
F = eX + eY + c((L+M)^6 + (L-M)^6), c = 5/7^7, L = X0+2X1+3X2, M = Y0+3Y1+5Y2 (verified exactly by exp0.py).
So F = sum of 8 sixth powers of linear forms l_1..l_8 = X0,X1,X2,Y0,Y1,Y2, g(L+M), g(L-M) (g^6 = c): X is the
codimension-2 linear section of the Fermat sextic sixfold V ⊂ P^7 by P^5 = {z7+z8 = 2gL, z7-z8 = 2gM}.
Steps:
 1. E1 kernel of [A] on the Jacobian ring, both Y-parities of D, d = 1..8 (first run did even d ≤ 6: all zero mod J).
 2. If kernel = 0 up to d0: every denominator of degree ≤ d0 is a polar divisor D = w(F); the E2 condition mod J is
    the quadratic system  w^T Hess(F) w ∈ J  (Koszul-syzygy argument, see §4). Solve it exactly for d = 5, 6
    (18 unknowns) with Singular/msolve → nonexistence theorem for deg D ≤ 6, then attempt d = 7.
 3. Structured denominators from the eight-power geometry: subset sums of the 8 powers (mixed parity handled by
    D·τ(D) or a parity-free builder), partial sums of F, E_X F, E_Y F, hyperplane-pencil denominators
    L^{6a} P(eX/L^6, M/L) (benediamond's (e,m) base), cones over the special curves/surfaces of X.
 4. Geometric (2,0)-budget analysis of the supports (as in attack1-dense-rational §1.3–1.4) for this F.

## 1. The obstruction and the geometric mechanism a certificate must realise

### 1.1 Set-up
* `F = Σ_{k=1}^{8} l_k^6` with `l_1..l_6 = X0,X1,X2,Y0,Y1,Y2`, `l_7 = g(L+M)`, `l_8 = g(L−M)`, `g^6 = c = 5/7^7`,
  `L = X0+2X1+3X2`, `M = Y0+3Y1+5Y2`.  `A = X0^3 Y0^3 = l_1^3 l_4^3`, `q = 1`.  `A ∉ J(F)` (D = 1 fails E1: rank 18 vs 19).
* The eight points `[l_k] ∈ P^5*` are NOT in general position: `l_1,l_2,l_3` span the plane `Π_X`, `l_4,l_5,l_6` span
  `Π_Y`, and `l_7, l_8` lie on the line joining `L ∈ Π_X` and `M ∈ Π_Y`.  The two linear relations are
  `R_1 = l_7 + l_8 − 2g(l_1+2l_2+3l_3)`, `R_2 = l_7 − l_8 − 2g(l_4+3l_5+5l_6)`.
* Jacobian ideal: `F_j = 6 Σ_k l_k^5 ∂_j l_k`, so `J(F)_5 = {Σ_k λ_k l_k^5 : λ ∈ U}` where `U ⊂ C^8` is the 6-dim
  space of evaluation vectors `(l_k(c))_k`, `c ∈ C^6` (= annihilator of `R_1, R_2`).  In particular
  `J(F)_d = { w(F) = 6 Σ_k l_k^5 w(l_k) : w a polynomial vector field of degree d−5 }`, and the second partials span the
  8-dim space `⟨l_k^4⟩`: `w^T Hess(F) w = 30 Σ_k l_k^4 w(l_k)^2`.  This is the whole content of "short power
  presentation" at the level of the Jacobian ring, and it is what makes the E2 analysis below tractable.
* Geometry: `X = V ∩ P^5`, `V ⊂ P^7` the Fermat sextic sixfold.  Divisors of `X` are `W ∩ P^5` for divisors `W ⊂ V`
  (all complete intersections).  The Fermat "subset-sum" divisors `W_I = V ∩ {Σ_{k∈I} z_k^6 = 0} = join(Fermat_I, Fermat_{I^c})`
  restrict to supports whose only non-rational singularities are: (|I| = 3) a plane sextic curve `S^3_{I^c} ∩ P^5` with
  transversal type cone over the Fermat sextic curve (transversal `p_g = C(6,3) = 20`) plus isolated points;
  (|I| = 4) finitely many isolated points (6+6) with transversal type cone over a sextic surface in `P^3` (`p_g = C(6,4) = 15`).
  A single such curve gives `≤ 280` adjoint conditions against `h^0(O_Y(6)) = 460` forms and the restriction is
  surjective (Taylor coefficients along the vertex plane + projective normality of plane curves), so `h^{2,0}(Ỹ) = 0`:
  **none of the eight-power subset sums can serve a (3,1)-class** (the Fermat mechanism needs TWO vertex planes spanning
  `P^5`, which forces `F` to be a Thom–Sebastiani sum in those coordinates — excluded by the Hessian-symmetrizer test).
  The plane `P^2_I ⊂ P^7` meets `P^5` in a point unless `I = {1,2,3}` or `{4,5,6}` (then a line: `{Y = 0, L = 0}` resp. `{X = 0, M = 0}`).

### 1.2 The (3,1)-support constraint (book ch. 3 §3.1) and the (2,0)-budget
Same as in `attack1-dense-rational.md` §1.3–1.4: `Y = X ∩ V(D)` must have a component whose resolution carries
`h^{2,0} ≥ 1` (non-normal or non-rational singularities), and
`h^{2,0}(Ỹ) = dim coker( S_d/(F S_{d−6} + C D) → H^0(Y, ω_Y/π_*ω_Ỹ) )`.  For this `F`:
* ternary cones `D = g(X)` / `g(Y)` (vertex `P_Y = {X=0}` resp. `P_X`): along `C_Y' = X ∩ {X = 0} = {eY + 2cM^6 = 0}` (smooth plane
  sextic, genus 10; note `F ≡ eY + 2cM^6 mod (X)^2`) the transversal type is the cone over `{g = 0}`; the adjoint
  restriction is surjective ⇒ `h^{2,0} = 0` (e.g. 540 conditions vs 785 forms for `d = 7`).  Same for `g(Y)`.
* `D ≡ F mod (binary form in L, M)`: support = union of hyperplane sections through the `P^3 = {L = M = 0}` ⇒ excluded.
* `D = eX + 2cL^6` (= `F|_{M=0}` extended): support singular along `C'' = {Y=0, eX + 2cL^6 = 0}` with transversal type
  `z^2 = (sextic)` (`p_g = 3`): 81 conditions vs 460 ⇒ excluded.
* Two-vertex-plane mechanisms are impossible for `d ≤ 8` (`l F ∈ (u)^6 + (v)^6` forces `F` to be a direct sum in `(u,v)`).
Remaining mechanisms (as for dense-rational): supports non-normal along a surface `Σ ⊂ X` that fails `d`-normality,
or a curve with a special (non-projectively-normal) position.  `X` contains the 36 lines `{L = M = 0, X0 = r_i X1, Y0 = s_j Y1}`
(roots of the binary sextics `b_X = eX|_{L=0}`, `b_Y = eY|_{M=0}`), the 6+6 cones `{L = 0, X0 = r_i X1, eY + 2cM^6 = 0}`,
`{M = 0, Y0 = s_j Y1, eX + 2cL^6 = 0}` over genus-10 curves, and the 136080 lines `P^3_{pairing,ζ} ∩ P^5` coming from
the standard `P^3`'s of the Fermat sixfold; none of these is obviously non-`d`-normal.
