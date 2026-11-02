# Research dossier: `cyclic-triple-residual` (OPEN)

Prepared 2026-11-01 from the repository at `/home/kejace/devel/sextics` (branch `main`, HEAD `83fbf5c`),
the saved discussion thread `/tmp/hodge/discussions/cyclic-triple-residual.md` (GitHub discussion #155,
two posts), the companion discussion #112 (fetched read-only), the companion's accepted `solution.md`,
`fqlx`'s two linked research packets (fetched read-only into `/tmp/hodge/ext/cyclic-residual-ext/`),
and the book chapters 01, 02, 03, 04. Nothing in the repository was modified; nothing was posted.

Status labels used below:

* **[P]** proven / exact computation (over Q, hence over C when the system is rational).
* **[C]** computational evidence only (modular rank, screened ansatz, bounded search).
* **[H]** heuristic, conjecture, or structural reading not proved.
* **[X]** attempted route that failed or stalled.
* **[N]** new in this dossier (computed here; scripts in `/tmp/hodge/scratch-residual/`).

---

## 1. Exact target

### 1.1 Statement

* `id: cyclic-triple-residual`, title "Cyclic triple cover · complementary character", `q = 1`,
  variables ordered `X0, X1, X2, Y0, Y1, Y2`.
* Target fingerprint (as quoted by `fqlx`, "the repository's target fingerprint, not the raw SHA256 of
  `instance.rec`"): `c5713bc06b37db0a64123364f8687a8ebc82e9722cc39ba88357659768eed0bc`.
* `A = X0^3 * Y0^3` (one monomial, degree 6 = 6q, Y-degree 3: Y-odd as required).
* `F` (12 terms, degree 6, Y-even):

```
F = X0^6 + X1^6 + X2^6 + Y0^6 + Y1^6 + Y2^6
  + (5/7) * ( X0^3*X1*Y0*Y1 + X0^3*X2^3 + X0^2*Y0*Y1*Y2^2
            + X0*X1*Y1^3*Y2 + X0*X2^3*Y1*Y2 + X1*X2^3*Y0*Y1 )
```

Monomial support (exponent vectors): (6,0,0,0,0,0), (0,6,0,0,0,0), (0,0,6,0,0,0), (0,0,0,6,0,0),
(0,0,0,0,6,0), (0,0,0,0,0,6) with coefficient 1, and (3,1,0,1,1,0), (3,0,3,0,0,0), (2,0,0,1,1,2),
(1,1,0,0,3,1), (1,0,3,0,1,1), (0,1,3,1,1,0) with coefficient 5/7.

### 1.2 Sparse presentation (verified by hand [P])

With `t = 5/7`, `Z = X2^3`,

```
C = X0^3 + X0*Y1*Y2 + X1*Y0*Y1                                   (cubic, X2-free)
H = X0^6 + X1^6 + Y0^6 + Y1^6 + Y2^6
    + t*(X0^3*X1*Y0*Y1 + X0^2*Y0*Y1*Y2^2 + X0*X1*Y1^3*Y2)          (sextic, X2-free)
F = Z^2 + t*Z*C + H.
```

Completing the square: `Z' = Z + t*C/2 = X2^3 + (5/14)*C`, `H' = H - t^2*C^2/4`, and `F = Z'^2 + H'`.

Per-variable exponent sets in `F`: X0 in {6,3,2,1,0}, X1 in {6,1,0}, X2 in {6,3,0}, Y0 in {6,1,0},
Y1 in {6,3,1,0}, Y2 in {6,2,1,0}. **X2 is the only variable in which F is a polynomial in a proper power**
(`X2^3`). `F` is linear in `Y0` outside `Y0^6`, but not in `Y1` (term `X0*X1*Y1^3*Y2`) or `Y2`
(`X0^2*Y0*Y1*Y2^2`); this matters for the cover trick of §2.7 (see §2.2 below).

Geometric shape:

* `pi : X -> W = {Z^2 + t Z C + H = 0} ⊂ P(1,1,3,1,1,1)` is the cyclic triple cover by `mu3: X2 -> zeta3 X2`,
  totally ramified along `R = {X2 = 0} ∩ X ≅ {H = 0} ⊂ P^4` (a sextic threefold; `H` has the mod-5 Fermat
  reduction, so it is smooth).
  `W` is the weighted Fano fourfold of index 2 (`-K_W = O_W(2)`), i.e. a double cover of `P^4`
  branched along the sextic threefold `{H' = 0}`; the deck involution of `W -> P^4` is `Z -> -Z - tC`.
* The composite `X -> P^4` (forget `X2`; the centre `[0:0:1:0:0:0]` is not on `X`) is finite of degree 6
  with fibres the roots of `x^6 + tC x^3 + H`; its Galois closure has group `Z/3 ≀ Z/2` of order 18
  (§1.7).
* `det(dPhi) = 3*X2^2` for `Phi = (X0, X1, X2^3, Y0, Y1, Y2)`, weights (1,1,3,1,1,1), sum 8.

### 1.3 Smoothness and nontriviality of the class [P]

`checks.json`: reduction mod 5 kills the coefficient 5/7 and gives the Fermat sextic, whose partials are
`Z_i^5`; so `X` is smooth mod 5 and hence over Q. The target is nonzero in `F_5[Z]/(Z_i^5)` (all exponents
<= 4); since `J(F)_6` has the same dimension 36 over Q and over F_5, the lattice `J(F)_6 ∩ Z_(5)[Z]` is
saturated and `A ∉ J(F)` over Q: the residue of `A` is a **nonzero primitive (3,1)-class** on `X`
(no trivial certificate). `K_X = O_X` (sextic fourfold in `P^5`), so `X` is a Calabi–Yau fourfold.

### 1.4 Projective diagonal symmetry group [P, recomputed N]

`G = {g ∈ (Z/6)^6 : g·m ≡ 0 (mod 6) for every exponent m of F} / <(1,1,1,1,1,1)>`, order **18**,
`G ≅ Z/6 × Z/3` (element orders: 1, one element of order 2, eight of order 3, eight of order 6).
Generators (weights mod 6 on `(X0,X1,X2,Y0,Y1,Y2)`):

```
mu3 : (0,0,2,0,0,0)    X2 -> zeta3*X2
tau : (3,3,3,0,0,0) ≡ (0,0,0,3,3,3)    Y -> -Y
g3  : (2,4,0,4,4,0)    (X0,X1,X2,Y0,Y1,Y2) -> (z*X0, z^2*X1, X2, z^2*Y0, z^2*Y1, Y2),  z = zeta3
```

Residue character of a monomial numerator `Z^e` under `g`: `chi_res(e,g) = Σ_i (e_i + 1) g_i mod 6`
(the `+1` is the determinant of the volume form `Omega`). For `A = X0^3 Y0^3`:
`chi_res(A,g) = 4 g0 + g1 + g2 + 4 g3 + g4 + g5 (mod 6)`.

Full table (weights; order; `chi_res(A)`; `chi_res` of the companion target `X0 X2^2 Y0^3`; ages of the
fixed components on `X`, Reid–Tai convention, computed as in `scripts/audit-frontier.py`):

| g | ord | chi(A) | chi(fano) | fixed blocks and ages |
|---|---|---|---|---|
| (0,0,0,0,0,0) | 1 | 0 | 0 | – |
| (0,0,2,0,0,0) | 3 | 2 | 0 | {X2=0}: age 1/3 (**junior**) |
| (0,0,4,0,0,0) | 3 | 4 | 0 | {X2=0}: age 2/3 (**junior**) |
| (1,5,1,2,2,0) | 6 | 2 | 2 | {X0,X2}-block 11/6; {Y0,Y1}-block 17/6 |
| (1,5,3,2,2,0) | 6 | 4 | 2 | {Y0,Y1}-block 13/6 |
| (1,5,5,2,2,0) | 6 | 0 | 2 | {Y0,Y1}-block 5/2; {X1,X2}-block 3/2 |
| (2,4,0,4,4,0) | 3 | 2 | 4 | {X2,Y2}-block 7/3; {X1,Y0,Y1}-block 4/3 |
| (2,4,2,4,4,0) | 3 | 4 | 4 | {X0,X2} 5/3; {X1,Y0,Y1} 5/3 |
| (2,4,4,4,4,0) | 3 | 0 | 4 | {X1,X2,Y0,Y1}-block: age exactly 1 |
| (3,3,1,0,0,0) | 6 | 4 | 0 | {Y0,Y1,Y2} 7/6; {X0,X1} 13/6 |
| (3,3,3,0,0,0)=tau | 2 | 0 | 0 | {Y}-block 3/2; {X}-block 3/2 |
| (3,3,5,0,0,0) | 6 | 2 | 0 | {Y0,Y1,Y2} 11/6; {X0,X1} 11/6 |
| (4,2,0,2,2,0) | 3 | 4 | 2 | {X2,Y2} 5/3; {X1,Y0,Y1} 5/3 |
| (4,2,2,2,2,0) | 3 | 0 | 2 | {X1,X2,Y0,Y1}: age exactly 1 |
| (4,2,4,2,2,0) | 3 | 2 | 2 | {X1,Y0,Y1} 4/3; {X0,X2} 7/3 |
| (5,1,1,4,4,0) | 6 | 0 | 4 | {X1,X2} 5/2; {Y0,Y1} 3/2 |
| (5,1,3,4,4,0) | 6 | 2 | 4 | {Y0,Y1} 11/6 |
| (5,1,5,4,4,0) | 6 | 4 | 4 | {Y0,Y1} 7/6; {X0,X2} 13/6 |

Consequences:

* **mu3 character of the residue of A is zeta3^1** (`chi_res = 2 mod 6`): the polynomial `A` is
  mu3-invariant, the volume form has character zeta3, so the residue is not invariant. The companion target
  `X0 X2^2 Y0^3` has numerator character zeta3^2 and invariant residue. This is the obstruction stated in
  `why.md`.
* **Stabiliser of the residue character of A** = the cyclic group `<s>` of order 6,
  `s = (1,5,5,2,2,0)`, i.e. `(X0,X1,X2,Y0,Y1,Y2) -> (w X0, w^5 X1, w^5 X2, w^2 Y0, w^2 Y1, Y2)`, `w = zeta6`,
  with `s^3 = tau`. **Every nontrivial element of `<s>` is senior** (ages 3/2, 5/2 for the order-6 elements,
  exactly 1 for `s^2, s^4`). Hence no quotient `X/K` by a subgroup `K ⊂ G` fixing the class is uniruled
  (Kollár–Larsen: the quotient of a Calabi–Yau by a group without junior elements has Kodaira dimension
  zero); the "finite quotients beyond reflections" mechanism of Chapter 3 is unavailable.
  This is exactly `checks.json`'s `juniorElementsFixingTarget: []` (the companion has
  `[(0,0,2,0,0,0),(0,0,4,0,0,0)]`).
* The residue character is also nontrivial under `g3` (`chi = 2`), so `A` is not even in the
  "mu3-invariant plus g3-invariant" sector; the full G-character of the residue is
  `(mu3, tau, g3) -> (2, 0, 2)` in the additive mod-6 notation of the scripts.

### 1.5 Non-diagonal symmetries, blocks [P]

* No coordinate reflection preserves `F` (`coordinateReflections: []`).
* Signed permutations of the variables preserving `F` (brute force over all 720·64 cases [N]): only the
  identity, `tau`, `-1`, and `-tau`. No variable permutation symmetry.
* `hessianSymmetrizerDimension = 1`: no linear direct-sum decomposition, even over C.
* `subsetEulerAtoms: []`: no coordinate-block Euler construction serves `A`.
* The full group `Aut(X) ⊂ PGL_6` has **not** been computed by anyone (gap); only the diagonal part (18)
  and signed permutations (trivial) are known.

### 1.6 Hodge-sector dimensions [N, exact]

Because the Koszul complex of the partials is `G`-equivariantly exact and `F` is `G`-invariant, the graded
`G`-character of `S/J(F)` equals that of the Fermat sextic (same subgroup `G ⊂ mu6^6`). Griffiths' residue
map is equivariant with the determinant twist. Counting Fermat-ring monomials (all exponents <= 4) in
degrees 6, 12, 18 therefore gives the exact `(h^{3,1}, h^{2,2}_prim, h^{1,3})` of every sector:

| sector | (3,1) | (2,2) | (1,3) |
|---|---|---|---|
| all of `H^4_prim(X)` | 426 | 1751 | 426 |
| `tau`-invariant (Y-odd numerators) | 208 | 866 | 208 |
| `V_0`: mu3 residue character 1 (served by the Fano quotient `W`) | 32 | 184 | 32 |
| `V_1`: mu3 residue character zeta3 (**contains A**) | 107 | 341 | 69 |
| `V_2`: mu3 residue character zeta3^2 (complex conjugate of `V_1`) | 69 | 341 | 107 |
| exact `G`-sector of `A` | 41 | 119 | 26 |
| its complex conjugate `G`-sector | 26 | 119 | 41 |
| exact `G`-sector of the companion target | 10 | 59 | 9 |

The tau-invariant part splits into 9 `G`-sectors with (3,1)-dimensions 41, 33, 33, 26, 22, 21, 13, 10, 9.
`V_1 = (41,119,26) ⊕ (33,112,21) ⊕ (33,110,22)`, `V_0 = (10,59,9) ⊕ (9,59,10) ⊕ (13,66,13)`.

The 41 degree-6 monomials in the exact `G`-sector of `A` are:
`X2^3*Y2^3, X2^3*Y1^3, X2^3*Y0*Y1^2, X2^3*Y0^2*Y1, X2^3*Y0^3, X1*Y1^2*Y2^3, X1*Y0*Y1*Y2^3, X1*Y0*Y1^4,
X1*Y0^2*Y2^3, X1*Y0^2*Y1^3, X1*Y0^3*Y1^2, X1*Y0^4*Y1, X1^2*X2^3*Y1, X1^2*X2^3*Y0, X1^3*Y2^3, X1^3*Y1^3,
X1^3*Y0*Y1^2, X1^3*Y0^2*Y1, X1^3*Y0^3, X0*Y1*Y2^4, X0*Y1^4*Y2, X0*Y0*Y2^4, X0*Y0*Y1^3*Y2, X0*Y0^2*Y1^2*Y2,
X0*Y0^3*Y1*Y2, X0*Y0^4*Y2, X0*X1*X2^3*Y2, X0*X1^2*Y1^2*Y2, X0*X1^2*Y0*Y1*Y2, X0*X1^2*Y0^2*Y2, X0*X1^4*Y2,
X0^2*X1*Y1*Y2^2, X0^2*X1*Y0*Y2^2, X0^3*Y2^3, X0^3*Y1^3, X0^3*Y0*Y1^2, X0^3*Y0^2*Y1, X0^3*Y0^3,
X0^3*X1^2*Y1, X0^3*X1^2*Y0, X0^4*X1*Y2`.

**Reading [H].** The smallest Q-Hodge substructure of `H^4(X,Q)` that can contain the class of `A` is
contained in the Galois-orbit sector `(41,119,26) ⊕ (26,119,41)`, Hodge numbers `(67, 238, 67)`. If that
Q-Hodge structure is simple (plausible for the "generic" coupling 5/7 within the `G`-symmetric family, not
proved), then by the argument of Chapter 4 §4.4 and Chapter 3 §3.1, **any divisor supporting the class of
`A` has resolved components `E_i` with `Σ h^{2,0}(E_i) >= 67`**, versus `>= 19` for the companion sector
`(10,59,9) ⊕ (9,59,10)`. Smooth or rationally-singular normal divisors contribute no (3,1)-class at all
(Chapter 3 §3.1). This is the quantitative reason why no low-degree denominator can exist (§4 item 2, §5.E1).

### 1.7 The Galois closure and the §4.6 exchange criterion [N, computed; conclusions H]

Let `Xbar2` be a cube root of the other root `Zbar = -tC - X2^3` of `y^2 + tCy + H`. The Galois closure of
`X -> P^4` is

```
Ztilde = { X2^3 + Xbar2^3 + t*C = 0 ,  X2^3 * Xbar2^3 = H } ⊂ P^6_(X0,X1,X2,Xbar2,Y0,Y1,Y2),
```

a complete intersection of type (3,6) in `P^6`, with group `Gamma = <d, e, sigma> ≅ Z/3 ≀ Z/2` (order 18;
`d: Xbar2 -> z Xbar2`, `e: X2 -> z X2`, `sigma: X2 <-> Xbar2`, `sigma d sigma = e`). `p: Ztilde -> X` is
the quotient by `<d>`, `p∘e = mu3∘p`, and the lift `r` of `tau` commutes with `sigma`. This is precisely
the data of the Chapter 4 §4.6 descent criterion, with `X/<mu3> = W` Fano (supported). The permutation
representation on the six roots is `1 ⊕ sign ⊕ rho ⊕ rho'` (two distinct 2-dimensional irreducibles);
`V_0` is the `sign`-isotypic part, `V_1`, `V_2` are the `rho`, `rho'` parts. Subgroups `K` with
`rho^K ≠ 0`: `{1}`, `<d>` (quotient `X`), `<e>` (quotient `X' ≅ X`), and the conjugates of `<sigma>`.

* `K_{Ztilde} = O(3+6-7) = O(2)` on the smooth locus.
* `Ztilde/sigma = { s1^3 - 3 s1 s2 + tC = 0, s2^3 = H } ⊂ P(1,1,1,1,1,1,2)` (`s1 = X2 + Xbar2`,
  `s2 = X2 Xbar2`), a weighted complete intersection (3,6) in weights summing to 8: `K = O(1)`.
  Modulo a singularity check this quotient is of general type, in any case not rationally connected.

So, exactly as in the book's §4.7 table for the natural `S3`-closure, **no uniruled intermediate quotient
of the natural Galois closure serves `V_1`**, and the holdout-style exchange trick (§2.7/§4.6) fails with
this cover. Other covers with an exchange are not excluded, but none is visible: the §2.7 shape
`y^6 + T^6 + yT·Λ(X) + f` requires the exchanged coordinate to enter the couplings linearly and paired with
a single partner, which fails for every `Y_i` of `F` (see §2.2).

### 1.8 The passive-ratio curve fibration used by `fqlx` [P, verified by hand N]

Put `X0 = h, X1 = v h, X2 = W, Y0 = a y, Y1 = b y, Y2 = y, t = y^2` (so `(a,b,v)` are passive ratios).
Then, with `c = 5/7`,

```
F = C(W,h,t) = W^6 + c*(h^3 + B*h*t)*W^3 + P*h^6 + Q*h^4*t + R*h^2*t^2 + S*t^3,
B = b(1+av),  P = 1+v^6,  Q = c*a*b*v,  R = c*(a*b + v*b^3),  S = 1 + a^6 + b^6 ,
```

a weighted plane curve of degree 6 in `P(1,1,2)[W,h,t]`: genus 4, `K = O(2)` with sections
`W^2, Wh, h^2, t`. The `mu3`-characters of `H^{1,0}` are `W^2 -> 1`, `Wh -> z^2`, `h^2 -> z`, `t -> z`:
the quotient `u = W^3`, `u^2 + c(h^3+Bht)u + (...) = 0 ⊂ P(1,2,3)`, is an **elliptic curve**, and the
Prym (dimension 3) carries the two conjugate rank-3 local systems of Hodge signature (2,1) and (1,2).
`fqlx`'s normalised form `W^6 + c(h^3+hT)W^3 + L h^6 + r h^4 T + T^3/x` (invariants `L, r, x`) is obtained by
`t = alpha T + beta h^2`, `W -> rho W`, with `L = (27PS^2 - 9QRS + 2R^3)/(3A^2)`, `r = (3SQ - R^2)/(BA)`,
`x = 3B^3/A`, `A = 3S - BR` (their notation). The physical source descends exactly to
`(a^3/2) * h^4 t^2 * Omega_(h,W,t) / C^2 ∧ da db dv` (affine Jacobian `h y^2`).

---

## 2. Companion and related solved instances

### 2.1 `cyclic-triple-fano` (same `F`; `A = X0*X2^2*Y0^3`) — SOLVED, PR #153, `benediamond`, accepted 2026-10-20

Method (nonconstructive, no witness, no degree bound; all in Lean):
`rbajaj5`'s free weighted line `[X0:X1:Z:Y0:Y1:Y2] = [0:s:0:is:t:it]` on `W` (discussion #112, 2026-10-17;
twisted normal map `H^0(O^5 ⊕ O(2)) -> H^0(O(5))` of rank 6, minor `-22500/2401`); anchored
six-coefficient deformation through an endpoint on `Y1 = 0` with Jacobian determinant `-22500/2401`;
étale lifting over the integral endpoint ring (specialisation `xi^6 + (T+1)^6 - 1`, Eisenstein);
polynomial normal; dominance via five commuting derivations; polynomial-family integration of the
quotient source `B = X0*Y0^3`; normalised finite trace; legal clearing; homogeneous extraction; cubic
transport `D = phi(d)`, `N_j = X2^2 phi(N_j)` (`j ≠ X2`), `N_{X2} = phi(N_Z)/3`. The denominator is
`phi(d)`, a polynomial in `X2^3` of unknown degree.

Why it fails here (`why.md`, §3.4 of the book, confirmed by §1.4): the transport multiplies every source
by `det(dPhi) = 3X2^2`; reachable sources are `X2^2 * B(X0,X1,X2^3,Y)` with X2-exponent ≡ 2 (mod 3), i.e.
residues in `V_0`. `A` has X2-exponent 0, residue in `V_1`. "No change of polynomial B repairs this
mismatch"; dividing by `X2^2` leaves the polynomial-source theorem and introduces the covering-field
character and poles. Also: the line `[0:s:0:is:t:it]` lies in the branch locus `Z = 0`.

### 2.2 Other solved instances with the same numerator `X0^3*Y0^3`

* `holdout` (PR #128, `fqlx`): chain sextic `F_H = Σ X_i^6 + Σ Y_i^6 + X0^3X1Y0Y1 + X1^3X2Y1Y2`.
  Cover `T^6 = Y0^6 + Y2^6`, exchange `y = Y1 <-> T` (a reflection on the cover), quotient
  `H = Σ X_i^6 + P(αX0^3X1 + βX1^3X2) + S^3 - 3SP^2` with Fano-threefold fibres over the Fermat curve
  `α^6 + β^6 = 1`; moving-conic primitive; target selection by the deck character `e = b - a - c - 1`
  (`e ≢ 0 mod 6` serves). Precedent for "reach a complementary character on a finite cover".
* `holdout-bundle` (PR #161, `fqlx`): second summand `X0^3Y0Y1^2` served by a **complementary diagonal
  action** `(Y0,Y1,Y2) -> (z^5Y0, zY1, z^5Y2)` (zeta6) with invariant residue and a junior fixed point
  (age 5/6), quotient coordinates `p = Y0Y1, r = Y2/Y0, R = Y0^6`, `det = 6Y0^5`, relative hypersurface
  over `C(λ)`, moving conics, then `HasCertificate.add`.
* `cycle-three` (PR #165, `hodgepodgeconjecture`): paired-sign sextic (even in every variable),
  reduction to a plane cubic over passive ratios, covariants `H, Θ, J`, Jacobian residue primitive
  (`ResidueTransport` library).

Why none transfers: (i) no coordinate reflections and no paired signs on our `F` (cycle-three needs
them); (ii) the stabiliser of our residue character has no junior element (holdout-bundle mechanism);
(iii) the §2.7 exchange needs `F = y^6 + T^6 + yT·Λ + (rest)` after adjoining `T`: for `Y0` (the only
`Y`-coordinate entering the couplings linearly) one would need the remaining `Y1,Y2`-terms to be
`T·(…)` with `T ∝ Y1` or `Y2`, but `X0^2Y0Y1Y2^2` and `X0X1Y1^3Y2` are not of that form; for `Y1`, `Y2` the
cubic/quadratic exponents already break linearity; (iv) the natural Galois closure is of general type
(§1.7).

### 2.3 Other related instances

* `reflection-coupled-fano` (solved, PR #150) / `reflection-coupled-residual` (open) and
  `quadric-contact-residual` (open): the same pattern "density sector solved, complementary character
  open" on other sextics. Nothing from their discussions applies directly; the notes in
  `/tmp/hodge/notes/accepted-group5.md` §4.2 recommend for all three: "look for a finite cover on which the
  complementary character becomes a nontrivial deck character, rather than a new curve on the already
  known quotient".
* `crossed-ternary-asymmetric` (PR #149, `1337isnot1337`) and `split-ternary-q2`: elliptic companions
  `u^3 = p` and `EllipticCompanions`/`HigherResidues` libraries; require a ternary split (absent here).

### 2.4 What the book says specifically

* Ch. 3 §3.4 (last paragraph): "A target lacking this density character is not obtained by substituting a
  different polynomial B. The free line solves the quotient integration problem; it does not remove the
  character restriction of the finite map."
* Ch. 3 §3.3: residue character formula `Σ (e_j+1) w_j`; junior-element criterion (Kollár–Larsen).
* Ch. 4 §4.7: "Neither an additional incidence seed on the same quotient nor a higher-order integration
  identity enlarges that quotient's averaging ideal."
* Ch. 4 §4.4: one generic nonzero certificate on a dominating family serves the whole invariant sector;
  `Σ h^{2,0}(E_i) >= 208` for a generic support. The `G`-symmetric family has 49 invariant sextic
  monomials (projective dimension 48 [N]); the analogous statement for the `(67,238,67)` sector is
  not in the book (see §5.E4).
* `provenance.json`: "No existence conclusion for this target follows from the screened
  quotient-density constructions."

---

## 3. Prior work, comment by comment

Discussion #155 (`https://github.com/hodge-computer/sextics/discussions/155`) contains exactly **two
posts, both by `fqlx` (prepared with OpenAI Codex agents)**. There are no other contributors.

### 3.1 Comment `DC_kwDOVLOOJc4BI59p`, `fqlx`, 2026-10-22T00:26Z (edited 20:27Z) — "No certificate yet"

Packet: gist `https://gist.github.com/fqlx/3e90a700d3631f9218680c27857f7741/b28663e84c720043e7200f8aa0a9e8ee58b6782b`
(README, `model_and_next_steps.md`, `rootfree_rr_bound.md`, `trace_dual_bound.md`, 15 Python/Sage audit
scripts, `source.json` (135/379/47-term polynomials `Tpoly, Rnum0, Delta0` in `eta, V, s`),
`original_scalar_input.json.gz.b64` (SHA256 `78219df9…`), `replay_receipt.json`, `sha256.json`).
Tools: SymPy, python-flint, SageMath 10.9.

Claims (all on the "worked invariant level" `I = 1/44, J = -121/12`, one normalised curve of the
elliptic-quotient family; **the packet does not independently derive the scalar reduction from the
six-variable target**):

1. **Degree-84 field model** [C→P in scope]: a compact polynomial `f` of degree 84 in `H` (leading
   coefficient `9 c^96 x^48`) with a universal sparse-resultant identity `Res_V(g1, g2) = 9 N^12 f`
   (70 terms, replayed over Q), rational inverse formulas for `eta, V, s`, audited in the Q quotient field
   at `b = 2`. The `b = 2` fibre is irreducible over Q; monic regular specialisation gives arithmetic
   irreducibility of the `Q(b)` model (not geometric irreducibility; ACTIVE.md later claims geometric
   irreducibility via finite-field factorisation, a smooth rational point and exact reduction).
   Integral coordinate: `W_int = t N A V`, `W_int^2 - t H D^2 W_int + t^2 N A^3 = 0`, clearing polynomial
   `t(b)` of degree 22, coefficients in `Q[b,Y]` of `(Y,b)`-degrees (3,18), (8,48). (Here `N, D, A` are
   packet-internal polynomials, not certificate data.)
2. **Charged scalar ODE** [X, unsolved]: `[(δ - a_κ)^2 - A0 (δ - a_κ) + C0] u = B_src`,
   `κ^3 = s/M^2`, `a_κ = δ log κ`, `δ = v(b) d/db` with explicit `v, A0, C0`; `rootfree_rr_bound.md` gives a
   complete characteristic-zero pole bound (`u ∈ L(E)`, `E = Pole(B_src) + 2 Pole(ω) + Σ h_P P`);
   `trace_dual_bound.md` gives a finite polynomial superset `u = Σ p_i(b) Y^i / (T f'(Y))` via the monogenic
   trace dual, avoiding a maximal order. **The resulting linear system was never assembled or solved.**
   Earlier runs: Q maximal order stopped after 35m26s, GF(101) saturation after 6m59s, GF(101)
   rational-function inverse after 5m40s, a Magma GF(101) run timed out at 60 s.
3. **Conjugate-character differential machinery** [C]: generic rank-3 "conjugate connection /
   polarisation", a complete horizontal frame on the worked level, explicit rational `C_t` correction;
   five scalar replay stages, the characteristic-zero integral-model audit and seven conjugate-frame
   stages passed. Essential gap: `G·W_t = G·W_W + d(G C) - dG ∧ C`; the `dG ∧ C` term and the actual `κ`
   charge are unmatched; "the horizontal frame alone supplies no source primitive".

Stated next steps: bounded specialisation/interpolation for the inverse in the `H` basis, polynomial HNF
of `[Y^i, W_int Y^i]`, normalisation of the enlarged order or direct assembly of the bounded scalar
system; compute the source coefficients and all `dG ∧ C` terms. "Any successful single-level construction
must still be lifted in the invariant parameters and transported to legal homogeneous certificate data."

### 3.2 Reply `DC_kwDOVLOOJc4BI_xu`, `fqlx`, 2026-10-23T01:06Z — "Final research handoff: research stopped"

Packet: `https://github.com/fqlx/sextics/tree/79bd381f50fbcbdf385835982afbac440c99b145/research-handoff`
(738 files; `README.md`, `manifest.json` with SHA-256 of every file, `research/ACTIVE.md`, rounds 5–11).
Six agents; "Do not resume the searches without a new instruction." Pinned repository commit
`cb551854693127a2542fb405f7ceb47ae0b35a09`. Python 3.12.14, SymPy 1.14.0, python-flint 0.9.0; some audits
need Singular.

Claims and status:

* **Exact fixed-support exclusion [P]** (`research/round10/fermat_cone_fixed/`): for
  `D = X2^6 + Y0^6 + Y1^6 + Y2^6` (reduced, legal), the complete normalised `q=1` system
  `Σ N_i F_i = A D`, `D div N - N(D) = F Q` (deg `N_i` = 7, deg `Q` = 6), projected onto the
  `mu3 × (global Y-sign)` character, has 2,760 rows, 1,248 columns (1,129 numerator + 119 quotient);
  modulo 1,000,003 the ranks are 1,236 / 1,237; twelve exact homogeneous kernels
  `N = F V - V(F) ε/6`, `Q = D div V - V(D)` for the mixed linear fields `V = Y_j ∂_{X_i}` (`i = 0,1`),
  `X_i ∂_{Y_j}` (`i = 0,1`, `j = 0,1,2`) are independent mod p, so the rational rank is exactly 1,236 and the
  system is inconsistent over C. By the prescribed-support theorem (Ch. 1), every multiplicity of this
  support is excluded. Independently rebuilt by a second agent (`round10/portfolio/audit_fermat_fixed.py`).
  **Replicated here [N]**: with the full 18-element reduction the system has 968 rows, 440 columns
  (395 + 45), exact Q-rank 436, augmented 437, kernel 4 (the four of the twelve fields with the right
  full character); `fqlx`'s 1129/119/2760 counts are reproduced exactly.
* **Second cone [C in packet → P here]** (`round10/source_center_fermat/result.json`):
  `D = X1^6 + X2^6 + Y1^6 + Y2^6`, ranks 1,232 / 1,233 mod 1,000,003 with the same twelve kernels, so the
  rank upper bound 1,236 is not sharp: "finite-prime inconsistency only". **Upgraded here [N]**: exact
  Q-rank 435, augmented 436 over the 440-column full-character system: excluded over C.
* **Generic one-direction scalar route excluded [P in its scope, X as a route]**
  (`round8/off_paired_B_resonance/resonance_obstruction.md`, ACTIVE.md): a two-place valuation
  obstruction over `C(eta, V, s)` — the necessary univariate equation forces a polynomial of degree <= 10
  to be divisible by one of degree 12. "This is not an exclusion of unrestricted certificates, and is not
  an exclusion on every special invariant level."
* **Rank-three transverse Prym connection** [C] (`round9/cyclic_full/`): curve
  `C = W^6 + (5/7)(h^3 + h t)W^3 + L h^6 + r h^4 t + t^3/x`; actual-character basis
  `(h^2 Ω/C, t Ω/C, t^4 Ω/C^2)`, conjugate basis `(Wh Ω/C, Wh t^3 Ω/C^2, W h^3 t^2 Ω/C^2)`; exact flatness
  over `Q(L,r,x)`; clean global frame with only the geometric poles `x`, `D1 = 27L^2 + 4r^3x`, `D2`
  (elliptic discriminant); complete free module `Der(-log(x D1 D2))` with determinant `-120050 x D1 D2`
  (Saito criterion); source identity `Res(h^4 t^2 Ω/C^2) = -∇_r ω_1`; coupled divergence equation
  `Σ_j [V_j(Z_j) + N_j^T Z_j + (V_j log H) Z_j] = -N_r^T e1 + 2z N_r^T e0 - z^2 N_L^T e0` with
  `dlog H = 3da/a + (13/3)dA/A - 5db/b - 4d(1+av)/(1+av) - (4/3)dS/S - dΔ47/Δ47`,
  `Jac = (1215/343) b^2 (1+av) S^2 Δ47 / A^5`.
* **Rank-five marked connection** [C] (`round10/horizontal5/`, `round9/open_curve/`): adds the
  graph-pole forms `σ0 = Res(h^4 Ω/((t - s h^2) C))`, `σ1 = Res(W^3 h Ω/((t - s h^2) C))`; all six
  curvature pairs vanish over `Q(L,r,x,s)`; poles `x, D1, D2, P = Lx + rxs + s^3`,
  `Q = 196Lx + 196rxs - 25xs^2 - 50xs - 25x + 196s^3`; third-kind projection
  `θ_χ = (1 - z^2) w^2 ((-q-u) σ0 - σ1)`; physical graph `s = T_eff = N1/N0`. Independent audits passed.
* **Bounded current searches, all negative [C, screened]**: constant five-component currents 30/31;
  degree-two physical currents 300/301; both-frame degree-two module 2400/2401; marked-log constant
  module 60/61; closed-seed module 156/157; old-mark trace-current modules 144/145, 144/145, 288/289
  (all mod 1,000,003); two nullflow trials with nonzero third-jet defects at four points. "These do not
  exclude the unrestricted horizontal route."
* **Effective-graph divisor `D32`** [C, geometry unfinished] (`round9/abel_prym/effective_graph.pickle`):
  the physical support candidate `Y2^2 N0 - X0^2 N1`-type divisor of degree 32, 558 terms, independent
  of `X2`, X-degrees 17 and 19 (hence **Y-odd**), multiplicity 17 along `X0 = X1 = 0`; six conductor
  surfaces `Y1 = 0, Y0^6 + Y2^6 = 0`; exact decomposition `D32 = q20 S6^2 + Y1 r25 S6 + Y1^2 s30`,
  normalisation generator `θ = q20 S6/Y1` (degree 25) with `θ^2 + rθ + qs = 0`, quotient `O_T(-25)`;
  generic collision with `Γ = {Y = 0} ∩ X` (`w^6 + (5/7)w^3 + 1 + v^6 = 0`, genus 10): plane degree-13
  cone with six collinear nodes, normalisation genus 60, graded adjoint lattices of ranks
  60,49,39,30,22,15,10,6,3,1 (total 235; 31,995 target slots) with explicit degree-32 lifts;
  Fermat/conductor corner: quasismooth weighted CI of degrees (12,12) in weights (3,2,6,6),
  discrepancy -8, transverse geometric genus 10, adjoint slots `1, eps, eps^2, b, r, W^3`,
  R^1 O characters (2,2,6), Y-even slots `1, eps^2, b, r`; coarse exceptional elliptic surface with
  `chi = 13, p_g = 12` whose two-forms have the **wrong** (mu3, Y) character ("actual source lives in
  conormal layers", first allowed normal-layer index `m ≡ 1 mod 6`); 36 infinity conductor lines with
  Newton vertices (0,6),(2,3),(6,1),(14,0), adjoint ideal `(u^2, b, r)`, transverse genus 2,
  homogeneous condition `P ∈ (S6, Y1^2, Y1 X0^2)`; actual infinity vertex of weights
  `(X0,X1,X2,Y1) = (7,5,5,16)`, `ord F = 30`, `ord D = 155`, discrepancy 32, canonical numerator must have
  `wt(P|X) >= 123`; main coarse residual component rational; special locus `196 v^6 + 171 = 0` weighted
  surface with three A1 and one A3 points. **Unresolved**: complete normalisation, endpoint/thickening
  corrections, global adjoint/Leray map, Gysin image of the source. "Neither nonzero global serving H20
  nor a primitive has been proved." No fixed-`D` linear test of `D32` is reported anywhere.
* Elliptic Fano companion / Jacobian-character counts (`round8/cyclic_companion/`): "the whole
  companion's even sector contains H41, so no unsupported coniveau-two inference is permitted"; special
  fixed line divisors lie over a source-empty hyperplane [C].
* Explicitly **withdrawn/invalid** items (handoff README): the Chebyshev one-dimensional interpolation
  shortage (wrong trace-dual filtration), `global_adjoint_with_RES.*`, the special-P0 adjoint sufficiency
  and the T0-surface-to-threefold rationality inference (unaudited).

### 3.3 Companion discussion #112 (`cyclic-triple-fano`)

* `rbajaj5`, 2026-10-17 (edited 10-19): the free line, the 6×8 normal map with minor `-22500/2401`
  (hence `H^1(φ^*T_W(-1)) = 0`), a seven-coefficient family with Jacobian `-135000 i/2401`, a 4×4 evaluation
  minor equal to 1 (dominant generically finite evaluation). "The source density `b = X0Y0^3/3` vanishes
  on the seed, while its `a`-derivative is `-i s^4/3`." Packet:
  `https://github.com/rbajaj5/hodge-lab/blob/d846f781066534411c8ed7398430b2ff68f810e6/notes/cyclic-fano-curve.md`.
* `benediamond`, 2026-10-19: independent reconstruction, `M M^* = diag(72, 25/49, 25/49, 25/49, 25/49, 72)`
  (rank 6); suggests computing boundary/branch/exceptional images; "a larger blind denominator search"
  is explicitly discouraged.

Both posts concern `V_0` only. Nobody in either thread discusses a divisor candidate for `V_1` other than
`fqlx`'s `D32` and the two Fermat cones.

### 3.4 Summary of rigorous exclusions versus screened failures

| item | status |
|---|---|
| `D = X2^6+Y0^6+Y1^6+Y2^6`, any multiplicity | excluded over C [P] (fqlx; replicated N) |
| `D = X1^6+X2^6+Y1^6+Y2^6` | modular only in packet; exact exclusion here [P, N] |
| all 63 Fermat-type cones `Σ_{i∈S} Z_i^6`, 2 random-coefficient full cones | excluded over C [P, N] |
| `C`, `Z' = 2X2^3 + tC`, `H`, `H'`, `X0Y0`, `X0^3Y0^3`, five polar quintics `∂F/∂Z_i`, `C Z'`, `H C`, `H Z'` | excluded over C [P, N] |
| random `G`-semi-invariant `D` of degree 1..5, every character sector (68 sectors) | excluded over C [P, N] |
| generic one-direction scalar route on the `(eta,V,s)` field | excluded in that scope [P]; not a certificate exclusion |
| bounded current modules (ranks k/k+1 mod p) | screened failures [C] |
| `D32` geometric audits | inconclusive |
| diagonal quotients `X/K` (junior-element method) | unavailable [P] |
| natural Galois-closure exchange (§4.6) | unavailable: `Ztilde/sigma` has `K = O(1)` [N, H] |

---

## 4. What remains unexplored

1. **No fixed-denominator test of `D32`** (degree 32, 558 terms, available in
   `research/round9/abel_prym/effective_graph.pickle`). All of `fqlx`'s effort on `D32` was geometric
   (normalisation, adjoints) and inconclusive. The definitive test is linear algebra (§5.E2).
2. **Degrees > 9 never swept**; the only degree-6 supports tested were two Fermat cones (now all 63,
   plus the natural divisors above). No `mu3`-invariant denominators of the form
   `d(X0,X1,X2^3,Y0,Y1,Y2)` beyond those listed, no products of polar divisors, no Hessian
   (degree 24), no discriminant-type divisors of the line family on `W`.
3. **The companion's own support divisor was never computed.** The §3.2 support theorem gives an
   algorithm (graph closure of the line family on `W`, boundary images, branch and positive-dimensional
   fibre loci). Its pullback `π^{-1}(E)` is a natural `mu3`-invariant candidate for `V_1` too; nobody has
   written it down or tested it.
4. **The Fano threefold of lines `𝓕` on `W` is unknown** (only one free line and its local structure).
   Its Hodge numbers, Albanese, and fibrations control `V_1 ⊕ V_2` through the elliptic family of §5.E5.
5. **Symmetry sectors never used beyond `mu3 × tau`.** The full order-18 group (`g3` in particular) cuts
   the systems by a further factor ~3 and splits `V_1` into (41,119,26) ⊕ (33,112,21) ⊕ (33,110,22). No
   search used it; no one computed the sector dimensions or the `Σ h^{2,0} >= 67` bound.
6. **No §4.4-type theorem for the `G`-symmetric family** (projective dimension 48): "one generic
   certificate for `X0^3Y0^3` on the cyclic-triple family serves the whole `(67,238,67)` sector", with the
   `Σ h^{2,0}(E_i) >= 67` consequence. The Higgs-field/matrix-unit computation of §4.4 can be repeated
   verbatim on the `G`-invariant deformation space (49 monomials) [not done].
7. **No nonlinear joint `(D, N)` search** at any degree, structured or not (expected empty at low degree
   by the `h^{2,0}` bound, but never run).
8. **No deformation/continuation**: at `t = 0` (Fermat) the target is served by the reflection
   `X0 -> -X0` (which the couplings destroy); nobody studied the one-parameter family `F_t` or the
   `G`-symmetric family for a certificate that spreads.
9. **No use of the isotrivial `j = 0` elliptic curves `E_L = π^{-1}(L)` over lines `L ⊂ W`**, the one
   structure of `X` that carries the twisted characters `V_1, V_2` with rank-one (finite-monodromy)
   coefficients (§5.E5).
10. **Full `Aut(X)`** not computed; the `H^{3,1}`-sector Hodge structure's simplicity not examined.
11. `fqlx`'s remaining concrete items: the charged scalar linear system on the worked level (never
    assembled), the `dG ∧ C` source matching, HNF normalisation of the degree-84 order — all parked, and all
    confined to one invariant level (a positive result would still need lifting in `(L, r, x)` and transport).

---

## 5. Candidate attack plan

Sizes below are for `q = 1`: unknowns are six numerators of degree `d + 1` and `Q` of degree `2d - 6`
(`deg E = 2d + 6 = deg F^2 Q`); equations are the coefficients of `E1 = Σ N_i F_i - A D` (degree `d + 6`)
and `E2 = D div N - N(D) - F Q` (degree `2d`). With a `G`-semi-invariant `D` the system splits into
characters: `chi(N_i) = chi_A + chi_D + e_i`, `chi(Q) = chi_A + 2 chi_D`, rows `chi_A + chi_D` (E1) and
`chi_A + 2 chi_D` (E2). Exact counts for `chi_D = 0` [N] (full reduction; unreduced sizes in the last
two columns):

| d | numerator cols | Q cols | E1 rows | E2 rows | total cols | total rows | unreduced cols | unreduced rows |
|---|---|---|---|---|---|---|---|---|
| 3 | 76 | 0 | 173 | 45 | 76 | 218 | 757 | 2,464 |
| 6 | 395 | 45 | 484 | 484 | 440 | 968 | 5,214 | 12,376 |
| 9 | 1,376 | 484 | 1,150 | 2,386 | 1,860 | 3,536 | 24,206 | 49,153 |
| 12 | 3,693 | 2,386 | 2,386 | 7,972 | 6,079 | 10,358 | 85,057 | 152,404 |
| 15 | 8,436 | 7,972 | 4,530 | 21,047 | 16,408 | 25,577 | 240,849 | 390,412 |
| 18 | 17,098 | 21,047 | 7,972 | 47,432 | 38,145 | 55,404 | 579,656 | 868,153 |
| 24 | 55,082 | 95,396 | 21,047 | 176,088 | 150,478 | 197,135 | 2,388,975 | 3,194,317 |
| 32 | 186,423 | 422,059 | 59,857 | 670,384 | 608,482 | 730,241 | 10,040,499 | 12,201,111 |
| 40 | 499,716 | 1,332,552 | 143,580 | 1,931,048 | 1,832,268 | 2,074,628 | 30,762,039 | 35,150,577 |

Exact Q-ranks with dense python-flint are instantaneous to d = 6, ~50 s at d = 9, and stall at d = 12;
beyond that a sparse modular solver (structured elimination / Wiedemann) is needed. Note that
`E1` alone is automatically solvable once `d + 6 > 24` (Jacobian socle degree), so all obstruction sits in
`E2` for `d >= 19`.

### E1. Complete the low-degree screen [done here, N]

All 63 Fermat-type cones, the natural divisors of §3.4, and random semi-invariant `D` of degree <= 5 in
all 68 character sectors are excluded exactly. Chance that *any* `D` of degree <= 9 works: essentially
zero (the `h^{2,0}` bound; every tested system has the right-hand side independent of the column space).
Remaining cheap items: random semi-invariant `D` of degree 6–9 in each sector (minutes each), products
of polar quintics, `F_{X2}/X2^2 · (…)`. Scripts: `/tmp/hodge/scratch-residual/fixedD_cones.py`
(`test(D)` builds the full-character system and returns exact ranks).

### E2. The definitive test of `fqlx`'s `D32` [not done; the single most informative computation]

`D32` is Y-odd (X-degrees 17, 19) and `mu3`-invariant; use `D32` itself as a (reduced, odd) denominator
with numerators of the transformed parity (Ch. 1 allows an odd reduced `D`), or `D32^2`. Full-character
system: ~608k columns × 730k rows (`D32` itself; for `D32^2`, degree 64, numerators of degree 65 and `Q` of
degree 122 give roughly 2×10^7 columns after reduction, so test the reduced odd `D32` first). Column density
~1,100 nonzeros (E2 part). Modular rank mod two large primes with a sparse solver; a modular solution
→ rational reconstruction → exact verification. A modular inconsistency is only evidence unless the
kernel is sealed. Expected: the geometric audits found several `p_g`-carrying pieces (12, 10, 60, 2) but
with wrong characters or unverified Gysin images; chance of success **~10–15 %**. Cost: a day of sparse
linear algebra engineering; the data are in the handoff (`effective_graph.pickle`).

### E3. Compute and test the companion's support divisor [not done]

Run the §3.2 algorithm on `rbajaj5`'s line family: closure of the graph of the evaluation map
`P^1 × S -> W`, images of the boundary, branch locus and positive-dimensional-fibre locus; take a
homogeneous equation `E_W(X0,X1,Z,Y)` of the union, pull back by `Z = X2^3`. This `D` is `mu3`-invariant and
supports **all** of `H^4(W̃)`; the question is whether `H^3(𝓕, L_χ)`-type classes also die on its
complement. Degree unknown (possibly large). Chance **~10–20 %**; a negative is rigorous only with the
rank sandwich. Cost: elimination in Singular (local `Singular 4.x` and the `ddg-cas` image), then E2-type
linear algebra.

### E4. A §4.4 theorem for the `G`-symmetric family and the `h^{2,0} >= 67` bound [theory, high chance]

Repeat Ch. 4 §4.4 verbatim with the deformation space `Def = R_6^+ ∩ (G-invariant sextics)` (49
monomials) in place of its `G`, `V` = the `(41)`-sector, `W` = its `(2,2)`-part (119), checking (i) distinct
support signatures, (ii) anchor connectivity, (iii) `Def·V = W`, at the Fermat point of the family
(`t = 0`). If the finite checks pass, then: a certificate
for `X0^3 Y0^3` at the generic point of the 48-dimensional `G`-symmetric family serves the whole sector,
and **every** support of the class on a generic member needs `Σ h^{2,0}(E_i) >= 67` (41 + 26 after
adding the conjugate sector). Chance of establishing the statement: **~80 %** (same finite computation).
Value: it rules out blind low-degree searches once and for all and tells what geometry a support must
have; it also says that one generic solve settles 41 + 33 + 33 + (conjugates) targets at once.

### E5. The isotrivial elliptic family over the lines of `W` [new structural route, H]

For a line `L ⊂ W` not in the branch locus, `E_L = π^{-1}(L) = {X2^3 = c_L(s,t)}` is a plane cubic with the
order-3 automorphism `X2 -> z X2` fixing three points: **`j = 0` for every `L`** (isotrivial family), with
`H^1(E_L) = chi ⊕ chibar` (no invariant part). Over the 3-dimensional family `𝓕` of lines
(`-K_W·L = 2`, `h^0(N_L) = 3` at a free line, dominant by `rbajaj5`), the universal curve `𝓔 -> X` is
generically finite and dominant, so `H^4(X) ↪ H^4(𝓔)`. The Leray pieces `H^4(𝓕)`, `H^2(𝓕) ⊗ R^2` and the
classes supported over the discriminant are divisor-supported on `𝓔` (hard Lefschetz on the threefold `𝓕`,
the 3-section `R ∩ 𝓔`), while `R^1` has no `mu3`-invariant part. Hence, modulo `N^1`:

```
V_0 ⊂ N^1 H^4(X)   (recovering the Fano theorem),
V_1 ⊕ V_2 ↪ IH^3(𝓕, L_chi ⊕ L_chibar),   L_chi = (R^1)_chi  a rank-one local system of finite monodromy,
```

trivialised on a finite cover `𝓕' -> 𝓕` (root ordering of `c_L` and the `mu6` ambiguity). Thus
`V_1 ⊕ V_2` is (a sub-Hodge structure of) `H^3(𝓕')_{chi'} ⊗ H^1(E_0)`, which forces
`h^{2,1}(𝓕')_{chi'} >= 107`, `h^{0,3}(𝓕')_{chi'} >= 69` (the target class is a `(2,1)`-class of `𝓕'` times
the holomorphic form of `E_0`). Divisor support on `𝓕' × E_0` of such classes cannot come from
`D × E_0` (that would need `N^1 H^3(𝓕')`); it must come from "graph-type" divisors, i.e. from maps of
(subvarieties of) `𝓕'` to `E_0` — the CM structure. **Experiments**: (a) compute `𝓕` (incidence of
`(5 linear forms, 1 binary cubic)` with 7 coefficient equations modulo `GL_2`), its singularities,
`h^{3,0}`, `h^{2,1}`, irregularity and the monodromy of `L_chi`; (b) look for fibrations of `𝓕'` with fibres
mapping to `E_0`. Chance of a *proof* within the challenge: **~10 %**; chance of clarifying why the sector
is hard: high. This is the only route in this dossier that uses structure specific to `X` rather than to
the quotient.

### E6. Degree-6/9 `mu3`-invariant denominators built from the line family [cheap, low chance]

Natural `mu3`-invariant divisors: the union of lines of `𝓕` tangent to the branch locus `{Z = 0}`
(discriminant of the elliptic family), the locus swept by lines through the six special points
`Y0^6 + Y2^6 = 0` etc., the surface-cones of §1.4's fixed loci. Each needs an elimination to get its
equation (Singular), then E1-type tests. Chance **< 5 %** individually (degrees are small, `h^{2,0}`
bound), but they document the obvious candidates.

### E7. Continue `fqlx`'s coupled-current / Prym-connection route [X-prone]

Assemble the charged scalar system on the worked level with the trace-dual superset; or extend the
rank-five current search beyond the bounded modules (all negative so far). Even a positive on one invariant
level needs generic `(L, r, x)` dependence and physical transport. Chance **~5 %**; cost high; the packet
is reproducible (`scalar_replay.py`, `conjugate_replay.py`, Sage 10.9 for some audits — Sage is not
available locally or in `ddg-cas`).

### E8. Nonlinear joint `(D, N)` search in a structured family [low chance]

Fix the shape `D = polynomial in (X0, X1, X2^3, Y0, Y1, Y2)` of weighted degree 6 (246 coefficients) or
degree 9, and treat consistency of the augmented system as a determinantal condition on `D`; solve by
Newton/homotopy on the bilinear system `E(D, N, Q) = 0` numerically, then reconstruct. Expected empty by
E4's bound (a sextic divisor on `X` cannot carry 67 holomorphic two-forms unless pathologically singular).
Chance **~2 %**; only worth running after E4 is settled negatively.

### E9. Full `Aut(X)` and hidden correspondences [cheap theory, low chance]

Compute the automorphism group of `F` beyond diagonal/permutation elements (e.g. via the Hessian or the
polar map); any non-diagonal automorphism not commuting with `mu3` would be a correspondence mixing
`V_0` into `V_1`. Expectation: `Aut = G` (order 18). Chance of a surprise **< 5 %**.

### E10. Other exchange covers [theory]

Look for covers `Z -> X` with an involution conjugating a deck transformation into the `mu3` generator
other than the natural Galois closure (`Ztilde/sigma` has `K = O(1)`). The structure
`F = Z'^2 + H'` suggests only the natural one. Chance **< 5 %**.

---

## 6. Feasibility assessment

* **Explicit witness.** Every tested denominator of degree <= 9 is excluded exactly, and the
  `h^{2,0} >= 67` heuristic (rigorous if E4 succeeds) says a support must contain components whose
  resolutions carry at least 67 independent holomorphic two-forms in the right characters. Cones over
  degree-`m` surfaces in `P^3` have `p_g = C(m-1,3)` (10 at `m = 6`, 165 at `m = 12`), so denominators of
  degree < 12 are essentially hopeless and a working `D` is likely to have degree in the 20s–30s
  (`fqlx`'s `D32` is in this range). Such a witness would have ~10^5–10^6 numerator coefficients after
  symmetry reduction; **kernel verification of the cleared identity for a degree-32 witness is far beyond
  what the accepted witness proofs (`d <= 4`) have done**, so even a found `D32` certificate would need a
  structured Lean proof (a decomposition of the identity along the construction), not `decide`.
* **Nonconstructive argument.** All accepted nonconstructive proofs made a geometric theorem effective
  through an explicit curve family with a polynomial normal, étale lifting, finite trace, legal clearing
  and transport (`HoldoutCover`, `PolynomialCurves`, `EllipticCompanions`, `ResidueTransport`,
  `HigherResidues`, `TransitiveCubic`). None of these handles a nontrivial deck character on the *same*
  quotient; the only character-changing precedents are (i) a different quotient with a junior element
  (unavailable here, §1.4) and (ii) an exchange on a finite cover (unavailable with the natural cover,
  §1.7). A proof for `V_1` therefore needs a **new** mechanism: most plausibly a family of curves on `X`
  whose parameter space carries the `chi`-character (E5), integrated with an elliptic/CM-type residue
  calculus. The Lean cost would be comparable to the crossed-ternary or cyclic-cubic proofs (thousands of
  lines) plus new library material for twisted traces.
* **What a Lean proof must supply** in either case: the protected statement
  `HodgeComputer.HasCertificate F A 1` via `hasCertificate_of_normalized` (exact `N(F) = A D`,
  `D div N - N(D) = F Q`, homogeneity, `IsRelPrime F D` over C — irreducibility of `F` is available by a
  monic Eisenstein specialisation as in the companion's `CyclicDomain.lean`, which can be reused since
  `F` is unchanged), and the coordinate-order binding. The companion's `CubicLift.lean` transport cannot be
  reused because its source transport is exactly the `X2^2` density.
* **Bottom line.** This target is a genuine Hodge-theoretic question (coniveau of a `(107,341,69)`
  twisted sector attached to a cyclic triple cover of a Fano fourfold). The screened constructions are
  exhausted; the discussion contains one exact exclusion (now extended to all low-degree natural
  candidates) and a large, inconclusive differential-geometric apparatus. The most valuable next steps are
  E4 (settles what any support must look like), E2 (the one untested concrete candidate), and E5
  (the only `X`-specific structure).

---

## Appendix A. Reproduction of the new computations [N]

Scripts (pure Python 3 + python-flint 0.9.0; local Singular 4.x available; `ddg-cas` docker image has
Singular 4.4.1, giac, python-flint, no SymPy/Sage/M2):

* `/tmp/hodge/scratch-residual/chars.py` — group `G`, element orders, residue characters, ages, sector
  dimensions, the 41 sector monomials (output reproduced in §1.4, §1.6).
* `/tmp/hodge/scratch-residual/fixedD_cones.py` — `test(D)`: full 18-character fixed-`D` system, exact
  Q-ranks via `fmpz_mat.rank()`; sweep of the 63 Fermat-type cones (29 s total).
* `/tmp/hodge/scratch-residual/batch2.py` — natural divisors (§3.4 table; degree-9 cases ~50 s each).
* `/tmp/hodge/scratch-residual/perm_and_random.py` — signed-permutation symmetries; random
  semi-invariant `D` of degree 1–5 in every character sector.
* `/tmp/hodge/scratch-residual/sizes.py` — DP table of §5.
* `/tmp/hodge/scratch-residual/deg12_modular.py` — modular (p = 1,000,003) screen of `H·C·Z'` (degree 12);
  the system (10,358 × 6,079) was built in 8 s but the dense `nmod_mat` rank had not finished after ~25 min
  when this dossier was closed, so **no result is recorded** (modular evidence only even when it finishes;
  rerun with a sparse solver).

Key numeric cross-checks: `fqlx`'s `mu3 × tau` column counts 1129 (numerator) / 119 (quotient) / 2760 rows
for `d = 6` are reproduced exactly; their exclusion of `X2^6+Y0^6+Y1^6+Y2^6` is reproduced over Q with
rank 436/437 on 440 columns.

## Appendix B. Links

* Instance folder: `data/cyclic-triple-residual/` (`instance.rec`, `why.md`, `provenance.json`, `checks.json`).
* Discussion #155: https://github.com/hodge-computer/sextics/discussions/155
* Companion discussion #112: https://github.com/hodge-computer/sextics/discussions/112
  (`rbajaj5` comment: https://github.com/hodge-computer/sextics/discussions/112#discussioncomment-19002605)
* Companion solve PR #153: https://github.com/hodge-computer/sextics/pull/153 (`data/cyclic-triple-fano/solution.md`)
* `fqlx` gist (2026-10-22): https://gist.github.com/fqlx/3e90a700d3631f9218680c27857f7741/b28663e84c720043e7200f8aa0a9e8ee58b6782b
* `fqlx` research handoff (2026-10-23): https://github.com/fqlx/sextics/tree/79bd381f50fbcbdf385835982afbac440c99b145/research-handoff
  — key files: `research/ACTIVE.md`; `research/round10/fermat_cone_fixed/{README.md,search.py,result.json,packet.pickle}`;
  `research/round10/source_center_fermat/result.json`; `research/round10/portfolio/audit_fermat_fixed.py`;
  `research/round8/coupled_source/README.md`; `research/round8/off_paired_B_resonance/resonance_obstruction.md`;
  `research/round9/cyclic_full/{README.md,CLEAN-FRAME.md}`; `research/round9/open_curve/README.md`;
  `research/round10/horizontal5/README.md`; `research/round10/horizontal5/fresh/boundary/README.md`;
  `research/round10/conductor_collision/README.md`; `research/round11/d32_corner_elliptic/README.md`;
  `research/round11/klein_coupled_handoff/VERIFIED_CYCLIC_PACKET_INDEX.md`;
  `research/round9/abel_prym/effective_graph.pickle` (the `D32` polynomial).
* `rbajaj5` curve packet: https://github.com/rbajaj5/hodge-lab/blob/d846f781066534411c8ed7398430b2ff68f810e6/notes/cyclic-fano-curve.md
* Book: `docs/01-certificates.md` (prescribed-support theorem, character averaging),
  `docs/02-explicit-constructions.md` §2.7–2.8 (reflection on a finite cover, complementary quotient),
  `docs/03-geometric-support.md` §3.1–3.4, 3.6, `docs/04-known-locus-and-frontier.md` §4.4, 4.6, 4.7.
* Fetched copies of the external material (untrusted data, read only):
  `/tmp/hodge/ext/cyclic-residual-ext/` (`disc112.json`, `gist.json`, `handoff_README.md`, `handoff/`).
