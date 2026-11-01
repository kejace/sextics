# Research dossier: `quadric-contact-residual` (OPEN)

Prepared 2026-11-01 from the repository snapshot at `/home/kejace/devel/sextics` (commit `83fbf5c`),
the saved discussion thread `/tmp/hodge/discussions/quadric-contact-residual.md` (GitHub discussion
#125: exactly **one** top-level comment, by `fqlx`, 2026-10-24, no replies; thread #119 of the solved
companion is empty), the companion and ancestor solutions, the book chapters cited by `why.md`
(§§1, 2.4–2.9, 3.1, 3.3, 4.2, 4.6, 4.7), the two curation scripts `scripts/audit-cubic-frontier.py`
and `scripts/audit-jacobian-scope.py`, and independent computations of my own (python-flint mod
`p = 1000000009`, sympy, msolve/Singular; scripts and logs in `/tmp/hodge/scratch/qcr-dossier/`).

Conventions. Variables are ordered `X0, X1, X2, Y0, Y1, Y2`; `τ(X,Y) = (X,−Y)`; `S = C[Z]`;
`J = (∂_0F,…,∂_5F)`; `R = S/J` (Jacobian ring, Hilbert function `1,6,21,56,126,246,426,666,951,1246,
1506,1686,1751,1686,…`); for a six-tuple `N`, `N(P) = Σ N_i ∂_iP`, `div N = Σ ∂_iN_i`,
`B_D(N) = D·div N − N(D)`. A certificate for `(F,q,A)` is `(D,N,Q)` with `D ≠ 0` homogeneous,
`gcd(D,F) = 1`, `deg N_i = deg D + 6q − 5`, and the cleared identity
`E := A·D² − q·D·N(F) + F·B_D(N) = F^{q+1}·Q`. Degree bookkeeping at `q = 1`, `d = deg D`:
`deg N_i = d + 1`, `deg E = 2d + 6`, `deg Q = 2d − 6`. Euler normalization (book §1) makes this
equivalent, for fixed `D`, to the two linear conditions `E1: q·N(F) = A·D` and `E2: B_D(N) ∈ (F^q)`.

Labels: **[proved]** (argument given/cited and checked), **[computed]** (exact or modular computation
by me, this session; mod-`p` ranks are lower bounds for ranks over `Q`, so "independent"/"not in the
ideal" conclusions transfer to `Q` and `C`, while "in the ideal" conclusions are modular evidence unless
stated), **[replicated]** (fqlx's claim recomputed by me and found identical), **[reported]** (claimed in
the thread or in accepted solutions; not replayed), **[inference]** (my reasoning), **[conjecture]**,
**[failed]**, **[pending]**.

---

## 0. Status at a glance

* Target: `q = 1`, `A = X0³Y0³ = P0³`, on the "quadric-contact" sextic `F_b` (the undeformed
  two-reflection sextic `F*` plus the perturbation `(5/7)·X0X1Y0Y1(X2² + 2Y2²)`). No certificate,
  no existence proof, no Lean theorem.
* **Rigorous exclusions over `C`** (`fqlx`, replicated where marked): every denominator of degree
  `≤ 6` [reported; E1 part replicated]; at degree 7 every `D` lying in one of the 8 parity sectors
  (hence every pure character sector) [reported]; at degree 8 the one-dimensional non-Jacobian E1
  direction `D8* = 4G_S² − G_P²` taken alone (plus the `F·S_2` gauge) [reported]; at every degree,
  certificates whose field keeps the ratios `t1 = Y1/X1`, `t2 = Y2/X2` passive (weighted-K3
  Gauss–Manin argument) [reported; Hilbert-series input replicated]; uniruledness of `X/H` for every
  `H ⊆ Γ⁺` (Reid–Tai) [reported].
* **Why the solved companion's method fails [proved]:** the companion `quadric-contact-broken`
  (same `F_b`, target `A' = X1²·X0Y0·(X0² − Y0²)`) is solved by a cubic-surface reduction through
  the quotient map `(X0,Y0) ↦ (S,T) = (2(X0²+Y0²), 4X0Y0)`, whose Jacobian `16(X0² − Y0²)` is odd
  under the exchange `s0: X0 ↔ Y0`. Every density target of that map is `s0`-odd; `X0³Y0³` is
  `s0`-even. More strongly (fqlx §4), *no* certificate at any degree can keep `t1, t2` passive, which
  the companion's does.
* **New data in this dossier [computed]:** (i) the full symmetry bookkeeping (projective monomial
  group of order 96, residue character `ε` of the target); (ii) the Hodge numbers of the target's
  character sector, `(h^{4,0},h^{3,1},h^{2,2},h^{1,3},h^{0,4})_ε = (0, 22, 51, 22, 0)` (two
  independent methods); (iii) E1 colon dimensions `dim (J:A)_d` for `d ≤ 11`, reproducing fqlx for
  `d ≤ 8` and showing that E1 stops being a useful filter from `d = 10` on (`556` non-Jacobian
  directions at `d = 10`, `1020` at `d = 11`, and `(J:A)_d ⊇` a subspace of codimension
  `dim R_{d+6}`); (iv) the degree-7 problem *with sectors mixed* is a genuine quadratic system
  (120 unknowns, 1496 quadrics, linearized rank 950) on which fqlx's square-column elimination removes
  nothing; a Gröbner computation is **[pending]** (see §3.5); (v) natural candidate denominators
  `S0^k, P0^k, δ^k` (`k ≤ 6`, `S0 = X0²+Y0²`, `δ = X0²−Y0²`) and `G_S², G_P², (2G_S ± G_P)²` all
  fail E1; `D8*·(any form)` passes E1 trivially; full fixed-denominator linear solves for
  `D8*` and `D8*·{S0, P0, δ, X0², Y0², S1, P1}` are **[pending]**.
* **Conjectural status [inference]:** the `ε`-isotypic part `V_ε ⊂ H⁴(X,Q)` is a rational
  sub-Hodge structure with no `(4,0)` part, so the generalized Hodge conjecture predicts
  `V_ε ⊆ N¹H⁴(X)`, i.e. a certificate for every one of its 22 `(3,1)`-directions, including the
  target. No geometric mechanism for that support is known: in both natural fibrations of `X`
  (K3 fibers over `(t1,t2)`, `(3,3)`-surfaces over `(c1,c2)`) the restricted class is
  transcendental.

---

## 1. The exact target

### 1.1 Data (`data/quadric-contact-residual/instance.rec`)

```
F = X0**6 + 5*X0**4*Y0**2 + 5*X0**2*X1**2*Y1**2 + 5*X0**2*X1*X2*Y1*Y2 + 5*X0**2*X2**2*Y2**2
  + 5*X0**2*Y0**4 + 5*X0*X1*X2**2*Y0*Y1/7 + 10*X0*X1*Y0*Y1*Y2**2/7 + X1**6 + 5*X1**2*Y0**2*Y1**2
  + 5*X1*X2*Y0**2*Y1*Y2 + X2**6 + 5*X2**2*Y0**2*Y2**2 + Y0**6 + Y1**6 + Y2**6
q = 1
A = X0**3*Y0**3
```

Monomial support of `F` (16 monomials; exponent vector in the order `X0,X1,X2,Y0,Y1,Y2`, coefficient):

| # | exponent | coeff | | # | exponent | coeff |
|---|---|---|---|---|---|---|
| 1 | (6,0,0,0,0,0) | 1 | | 9 | (2,2,0,0,2,0) | 5 |
| 2 | (4,0,0,2,0,0) | 5 | | 10 | (2,1,1,0,1,1) | 5 |
| 3 | (2,0,0,4,0,0) | 5 | | 11 | (2,0,2,0,0,2) | 5 |
| 4 | (0,0,0,6,0,0) | 1 | | 12 | (0,2,0,2,2,0) | 5 |
| 5 | (0,6,0,0,0,0) | 1 | | 13 | (0,1,1,2,1,1) | 5 |
| 6 | (0,0,0,0,6,0) | 1 | | 14 | (0,0,2,2,0,2) | 5 |
| 7 | (0,0,6,0,0,0) | 1 | | 15 | (1,1,2,1,1,0) | 5/7 |
| 8 | (0,0,0,0,0,6) | 1 | | 16 | (1,1,0,1,1,2) | 10/7 |

Coefficient set `{1, 5, 5/7, 10/7}`; all non-Fermat coefficients are `≡ 0 (mod 5)`.

### 1.2 Structured presentation [proved]

Put `S_i = X_i² + Y_i²`, `P_i = X_iY_i` (`i = 0,1,2`), `δ = X0² − Y0²`,
`B = P1² + P1P2 + P2² = X1²Y1² + X1X2Y1Y2 + X2²Y2²`. Then

```
F_b = S0³ + 2·S0·P0²                                  (= X0⁶ + 5X0⁴Y0² + 5X0²Y0⁴ + Y0⁶)
    + 5·S0·B
    + (5/7)·P0·P1·(X2² + 2·Y2²)                       (the deformation term, t = 5/7)
    + (X1⁶ + Y1⁶) + (X2⁶ + Y2⁶)                       (= S1³ − 3S1P1² + S2³ − 3S2P2²)
```

`F* := F_b − (5/7)·P0P1(X2²+2Y2²)` is the sextic of the solved ancestor
`two-reflection-cubic-surface-d2` (even in `X0` and in `Y0` separately). The deformation term is odd
in `X0`, odd in `Y0`, and breaks `X2 ↔ Y2` (coefficients `5/7` vs `10/7`). The `(X0,Y0)`-block enters
`F_b` only through `S0` and `P0`; this is the basis of the companion's quotient map.

In fqlx's notation `G_S = ∂F/∂S0 = 3S0² + 2P0² + 5B`, `G_P = ∂F/∂P0 = 4P0S0 + (5/7)P1(X2²+2Y2²)`,
so that `∂_{X0}F = 2X0G_S + Y0G_P`, `∂_{Y0}F = 2Y0G_S + X0G_P` [replicated].

Paired structure: every monomial has even degree in each pair `(X_i,Y_i)`, so `F_b` is invariant under
the three paired sign changes `σ_i: (X_i,Y_i) ↦ −(X_i,Y_i)`. Fixing `t_i = Y_i/X_i` and writing
`u_i = X_i²` gives the plane cubic over `C(t0,t1,t2)` [computed from the presentation above]:

```
U_t(u) = (1 + 5t0² + 5t0⁴ + t0⁶)·u0³ + (1 + t1⁶)·u1³ + (1 + t2⁶)·u2³
       + 5(1 + t0²)·u0·(t1²u1² + t1t2·u1u2 + t2²u2²) + (5/7)·t0t1(1 + 2t2²)·u0u1u2 .
```

Its support is `{u0³, u1³, u2³, u0u1², u0u1u2, u0u2²}` — no `u0²u_i`, no `u1²u2`, no `u1u2²`.
The map `(X_i,Y_i) ↦ (u_i, t_i)` has constant Jacobian `2` per pair (`8` total), so *every* target is a
polynomial density on this quotient; `A = t0³u0³`. The sextic is **not** in the 26-coefficient
`Γ`-family of book §2.9/§4.2: the weight-`(2,1,0,3,2,1)` action does not preserve it (e.g.
`X0⁴Y0²` has weight `14 ≡ 2 (mod 6)`), as `why.md` states.

### 1.3 Parity and degrees

`F_b` is `Y`-even, `A` is `Y`-odd (`checks.json`: `involutionParities [1,−1]`). `deg F = deg A = 6`.
At `q = 1`: `deg N_i = d + 1`, `deg Q = 2d − 6`. Component parities: `N_0,N_1,N_2` `Y`-odd,
`N_3,N_4,N_5` `Y`-even when `D` is `Y`-even (a common denominator can always be made `Y`-even).

### 1.4 Symmetry group and characters

**Diagonal symmetries [computed, agrees with fqlx and `checks.json`].** The difference lattice of the
support has Smith form `Z^6/M ≅ Z ⊕ (Z/2)² ⊕ Z/6`, so the projective diagonal group is
`G_diag ≅ (Z/2)² × Z/6`, order 24. The exact-invariance group `{λ : λ^e = 1 ∀ e ∈ supp F}` has Smith
form `(Z/2)² × (Z/6)²`, order 144 `= 48 · 3`: it is fqlx's 48-element group
`G_aff = ⟨σ0, σ2, h1, τ_Y⟩` times the scalar `μ3`. Same projective image. Generators:

```
σ0: (X0,Y0) ↦ (−X0,−Y0)        σ2: (X2,Y2) ↦ (−X2,−Y2)
h1: (X1,Y1) ↦ (ζ6 X1, ζ6⁻¹ Y1)  τ_Y = τ: (Y0,Y1,Y2) ↦ −(Y0,Y1,Y2)
```

(Note `σ1 = σ0σ2·(−1)` projectively, and `−I = σ0σ2h1³`.) A monomial `Z^e` has fqlx's character
`χ(e) = (e0+e3 mod 2, e2+e5 mod 2, e1−e4 mod 6, e3+e4+e5 mod 2) ∈ (Z/2)²×Z/6×Z/2`, with
`χ(F) = 0`, `χ(A) = (0,0,0,1)`.

**Non-diagonal symmetries [computed].** Among all `720` coordinate permutations combined with diagonal
scalings, exactly four preserve `F_b` up to scalar: `id`, `s0: X0 ↔ Y0`, `r1: X1 ↔ Y1`, `s0r1`.
(`X2 ↔ Y2` is broken by the deformation; no other exchange or cyclic permutation survives.) The
visible projective symmetry group is therefore `Γ = G_diag ⋊ ⟨s0, r1⟩`, order 96 (fqlx; replicated).
As a linear group with scalars `{±1, ±ζ6…}` it has order 192 in my enumeration. Hidden non-monomial
linear symmetries are not excluded by this check; the companion's `checks.json` records a scalar
Hessian symmetrizer space (no hidden linear block split).

**Residue character of the target [computed, agrees with fqlx].** For `g ∈ Γ` with `F(gz) = c_g F(z)`
(`c_g = 1` for all 192 elements), the residue `Res(PΩ/F^{q+1})` transforms by
`ψ_P(g)·det(g)/c_g^{q+1}`, where `P(gz) = ψ_P(g)P(z)`. For `A = X0³Y0³`:

| generator | `ψ_A` | `det` | residue character `ε` |
|---|---|---|---|
| `σ0` | +1 | +1 | +1 |
| `σ2` | +1 | +1 | +1 |
| `h1` | +1 | +1 | +1 |
| `τ_Y` | −1 | −1 | +1 |
| `s0` | +1 | −1 | **−1** |
| `r1` | +1 | −1 | **−1** |

So the residue is trivial on all of `G_diag` (fqlx: "all fixing `AΩ/F²`") and odd under both
exchanges. `Γ⁺ := ker ε = G_diag ⋊ ⟨s0r1⟩` (order 48), `Γ⁻ = Γ \ Γ⁺` contains the 8 reflections
`s0, s0σ0, r1h1^k (k = 0..5)` (all of Reid–Tai age `1/2`), each fixing `A` and negating `Ω`.

**Hodge numbers of the `ε`-sector [computed, two methods agree].** Taking the `ε`-isotypic part of
the residue representation on `H⁴_prim(X)` (numerator character `ψ_A·c^{q−1} = ψ_A` in degree `6q`):

```
(h^{4,0}, h^{3,1}, h^{2,2}, h^{1,3}, h^{0,4})_ε = (0, 22, 51, 22, 0)       (rank 95)
```

(degree 6: `dim S_6^ψ = 25`, `dim J_6^ψ = 3`, so 22 independent targets in the sector; degree 12:
`191 − 140 = 51`; confirmed by the equivariant Koszul trace formula
`tr(g|R_d) = Σ_k (−1)^k e_k(μ⁻¹) h_{d−5k}(μ)`.) For comparison the full primitive cohomology is
`(1, 426, 1751, 426, 1)` and the `τ`-invariant residues form `(0, 208, 866, 208, 0)` (book §2.9).
The `(4,0)` class `Res(Ω/F)` has character `det`, which differs from `ε` on `τ_Y`; hence `h^{4,0}_ε = 0`.

**The companion's target** `A' = X0³X1²Y0 − X0X1²Y0³ = X1²·P0·δ` is *not* `G_diag`-invariant as a
residue: under `h1` it picks up `ζ6² = ζ3` (`X1²`), and `r1` swaps the `ζ3` and `ζ3⁻¹` sectors, so
`A'` lies in a 2-dimensional irreducible `Γ`-representation. Under the subgroup `G_diag ⋊ ⟨s0⟩` its
character sector has `h^{3,1} = 3` [computed]. The two targets on the same variety live in
different, non-interacting isotypic components.

### 1.5 Smoothness [proved; `checks.json`]

Reduction mod 5 kills every non-Fermat coefficient (`5, 5/7, 10/7` all have numerator divisible by
5 and denominator prime to 5), leaving `Σ(X_i⁶ + Y_i⁶)`, whose partials `6Z_i⁵` have no common
projective zero over `F_5`. Smoothness being open in flat families, `X = V(F_b) ⊂ P⁵` is smooth over
`Q` and `C`; `(∂_0F,…,∂_5F)` is a regular sequence (used throughout by fqlx and below).

### 1.6 Block / split structure

No linear split: `F_b` is not of the form `f(U) + g(V)` or a function of fewer quadrics; the
`(X0,Y0)` block is coupled to `(X1,Y1,X2,Y2)` through `5·S0·B` and `(5/7)·P0·P1·(X2²+2Y2²)`.
The only "block" structure is the paired-sign structure of §1.2 and the `s0`-exchange on the first
pair. `checks.json`: `subsetEulerAtoms = []`, `coordinateReflections = []` (no single-coordinate sign
symmetry survives the deformation), `degreeSixJacobianRanks [36, 37]` (the 36 products
`Z_i·∂_jF` are independent and `A ∉ J_6`, i.e. `A` is a nonzero primitive class).

---

## 2. Companion and related solved instances

### 2.1 `quadric-contact-broken` (same `F_b`; `A' = X1²·X0Y0·(X0²−Y0²)`) — SOLVED

* Solver `HeadOfNaps`, PR #122, accepted 2026-10-19 (`acceptance.json`: 56 Lean modules, strict
  certificate, `witnessChecked: false`, no `witness.json`; axioms `propext, Quot.sound, Classical.choice`).
* **Method ("finite-algebra contact certificate", book §2.4).** Coordinates `(x,a,c,y,b,d) =
  (X0,X1,X2,Y0,Y1,Y2)`; passive ratios `α = b/a = t1`, `β = d/c = t2`; active quadratic coordinates
  `S = 2(x²+y²)`, `T = 4xy`, `B = a²`, `C = c²`. Then `G = 16·F_b` with
  `G = 2S³ + ST² + S(pB² + qBC + sC²) + kTBC + uB³ + vC³`,
  `(p,q,s,k,u,v) = (40α², 40αβ, 40β², (20/7)α(1+2β²), 16(1+α⁶), 16(1+β⁶))`,
  a cubic surface over `K = Q(α,β)`. A line `S = ξB + ηC, T = ζB + ωC` on it is constructed over the
  degree-27 finite free algebra `E = O[r]/(P)` (`O` = localization of `Q[α,β]` at the fiber `(1,2)`;
  `P` monic of degree 27 from a resultant). A divergence-free tangent field along the line sends
  `G ↦ L·J` (`J` quadratic); an invertible `8×8` minor of the source matrix gives a functional sending
  `J ↦ TB`. The regular representation of `E` and the adjugate of the `27×27` pencil `ρ(L)` give the
  cubic-level denominator `D0 = det ρ(L)` (degree 27 in `(S,T,B,C)`; a product of 27 conjugate
  "plane" linear forms over the algebraic closure) and numerators.
* **Explicit denominator?** Only semi-explicit: after clearing passive denominators and homogenizing
  `α = b/a, β = d/c` by `(ac)^n` (exponent `n` from an existential lemma), `D = a²c²·Q` of degree
  `2n + 58`, numerators of degree `2n + 59`, `Q = 0` (strict). Legality is proved on the fiber
  `(a,c,y,b,d) = (1,1,0,1,2)` where `F_b ↦ x⁶ + 35x² + 67` and `D ↦ κ·det(2x²I − M)`, with a matrix
  Bézout identity.
* **Density identity.** `det ∂(S,T)/∂(x,y) = 16(x² − y²) = 16δ` and `A' = (δ/4)·BT`. The target is
  exactly `(Jacobian density) × (polynomial source BT)`.
* **Why it fails here [proved; `why.md`, book §2.5 "density boundary"].** The quotient map
  `Φ = (S,T,B,C,α,β)` is fixed by the linear reflection `s0: X0 ↔ Y0`; differentiating `Φ∘s0 = Φ`
  gives `det(dΦ)(s0 z) = −det(dΦ)(z)`. Hence every density target `det(dΦ)·b(Φ)` — and every
  `h(α,β)`-multiple of one (passive multipliers, §2.5) — is `s0`-odd. `A = X0³Y0³` is `s0`-even.
  Averaging over `⟨s0⟩` kills all reachable targets and fixes `A`, so `A` is not in the linear span of
  that density sector. (`why.md` also notes that `r1: X1 ↔ Y1` fixes `A`; `r1` does not fix `Φ`
  (`α ↦ 1/α`), so for `r1` this is an observation about the character, not a density obstruction.)
* **Stronger statement [reported, fqlx §4].** No certificate at any degree keeps `t1 = α`, `t2 = β`
  passive — see §3.4. The companion's field is exactly of that passive type.

### 2.2 `two-reflection-cubic-surface-d2` (the undeformed `F*`; same `A'`) — SOLVED (`benediamond`)

`D = 4(X0² + Y0²)`, `N = (−Y0X1², 0, 0, X0X1², 0, 0)` (rotation field; strict). Uses
`F*_u − F*_w = −2(u−w)(u+w)` for `u = X0², w = Y0²`. On `F_b` the same `(D,N)` leaves the defect
`N(F_b) − A'D = (5/7)·X1³Y1·(X0²−Y0²)(X2²+2Y2²)` (degree 8, not divisible by `D`), and `why.md` of the
companion proves that **no** numerator works with `D = X0² + Y0²` (rank `336/337` of the 336 cubic
coefficients mod 101; full column rank ⇒ characteristic-zero inconsistency) [reported].

On `F*` the residual target `X0³Y0³` is odd in `X0` and `F*` is even in `X0`, so the reflection
theorem (book §3.3) gives geometric existence of a certificate for `(F*, 1, X0³Y0³)` through the
Fano quotient `X*/⟨X0 ↦ −X0⟩ ⊂ P(2,1,1,1,1,1)`. The deformation term is odd in `X0`; this reflection
does not survive on `F_b`. No instance `(F*, X0³Y0³)` is in the catalogue.

### 2.3 Other related constructions (why they do not apply)

* **Reflection quotients (§3.3).** `F_b` has no single-coordinate reflection. The 8 reflections in
  `Γ⁻` (`s0`, `s0σ0`, `r1h1^k`) do give Fano quotients (`fqlx` identifies `X/⟨s0⟩` with a degree-3
  hypersurface in `P(1,1,1,1,1,2)`) and support every residue *odd* under the chosen reflection; the
  target is even under all eight [proved by the character table above].
* **Paired-sign / Jacobian transport (§2.9, §4.2).** Needs the extra weight-`(2,1,0,3,2,1)` order-6
  symmetry (absent) to kill the two holomorphic 2-forms of the `(3,3)` surface; `why.md` and
  `scripts/audit-jacobian-scope.py` show the generic auxiliary elliptic surface here (over the
  `T = t0²` line at `c1 = 2, c2 = 3`) has coprime Weierstrass coefficients of degrees `11, 18`, hence
  is minimal with `χ(O) = 3`, Kodaira dimension 1, not uniruled: a new incidence seed with the
  passive parameters fixed cannot give a dominant rational-curve family.
* **Fermat cover (§2.7) / complementary quotient (§2.8).** Both are built for the chain sextic's
  specific couplings; neither the `T⁶ = Y0⁶ + Y2⁶` cover nor the weight-`(5,1,5)` quotient
  preserves `F_b`. The §2.7 mechanism ("create a reflection on a finite cover whose deck character
  selects `X0³Y0³`") is the closest template for a transfer serving the complementary sector — see
  §5, experiment 7.
* **Weighted-Fano / conic constructions (§2.6, §3.4)** serve determinant densities of polynomial
  Piola maps; no such map is visible for `F_b` beyond the companion's, which has the wrong character.

---

## 3. Prior work in the discussion (#125), comment by comment

**Only one comment exists.** `fqlx`, node id `DC_kwDOVLOOJc4BJFYx`, created and last updated
2026-10-24T03:37:33Z, title "Exact Low-Degree Exclusion (`deg D ≤ 6`), Degree-7 Sector Analysis,
Degree-8 Non-Jacobian Generator `D8*`, and Weighted K3 Gauss–Manin Obstruction". No replies.
AI disclosure: computed with the Antigravity coding agent on behalf of `fqlx`. Reproduction: a
standalone script `verify_quadric_contact_residual.py` is embedded in the comment (exact rational
arithmetic with `fractions`, no external packet URL); it covers the E1 colon dimensions for
`d ≤ 7`, the identity `N*(F) = A·D8*`, and (per its docstring) the Reid–Tai and K3 checks.
Discussion #119 (companion) has no comments. The companion's PR is
<https://github.com/hodge-computer/sextics/pull/122>.

### 3.1 §1 of the comment: Koszul reduction and exclusion of `deg D ≤ 6` over `C`

**Argument [proved; I re-derived it].** For `D ∈ J_d`, write `D = Σ C_i ∂_iF`. Then `N0 = A·C`
satisfies E1, and any E1 solution is `N = A·C + K` with `K` a degree-`(d+1)` syzygy of `∇F`, hence
Koszul (regular sequence): `K = Σ_{r<s} Q_rs(∂_sF e_r − ∂_rF e_s)`, `deg Q_rs = d − 4`. Modulo
`J_{2d}`: `B_D(K) ∈ J` (`div K ∈ J` because mixed second partials cancel, `K_ℓ ∈ J`), `F·Q ∈ J`
(Euler), and `B_D(AC) ≡ −A·Σ_{i,j} C_iC_j ∂_i∂_jF`. So a **necessary** condition for E2 is the
quadratic "Hessian congruence"

```
H(C) := A · Σ_{i,j} C_i C_j ∂_i∂_j F  ≡ 0   (mod J_{2d}).
```

The gauge `(D,N) ↦ (D + tF, N + tAZ/6)` corresponds to `C ↦ C + tZ/6` and lies in the kernel of the
quadratic map (`Σ Z_iZ_j∂_ijF = 30F`, cross terms `Σ_i Z_i∂_ijF = 5∂_jF`).

**Table (fqlx; E1 column replicated by me mod `p`):**

| `d` | `dim S_d` | `dim (J:A)_d` | `dim J_d` | status |
|---|---|---|---|---|
| 0–4 | 1, 6, 21, 56, 126 | 0 | 0 | excluded at E1 [replicated] |
| 5 | 252 | 6 | 6 | E2: the 21 columns `A·∂_i∂_jF` have rank 21 in `S_10/J_10` (dim 1506) ⇒ all `c_ic_j = 0` [reported; linearized full rank, rigorous if correct] |
| 6 | 462 | 36 | 36 | E2: gauge-fix `m_55 = 0` (35 parameters `m_ij`, `C_j = Σ m_ij Z_i`); 630 columns of rank 334 in `S_12/J_12` (dim 1751); step 1: 28 off-diagonal parameters killed by independent square columns; step 2: on `(m00, m03, m11, m22, m30, m33, m44)` the 28 columns have rank 15 with explicit relations `m22m30 = 0`, `m22² − 3m30² = 0`, `m03² − m30² = 0`, `m22m33 − (19/6)m30² − m33² = 0`, then `m00² = 0`, `(3/2)m11m44 = 0`, `m11² + m44² = 0` ⇒ all zero over `C` [reported] |

Conclusion (fqlx): **every certificate denominator has `deg D ≥ 7` over `C`.** The E1 part is
rigorous and replicated; the E2 closures at `d = 5, 6` are exact computations I did not replay but
whose logic is sound (independent square columns force vanishing over any field; the explicit
rational relations are stated as exact).

### 3.2 §2 of the comment: degree 7 (sectors) and the degree-8 generator `D8*`

* **Degree 7, E1 [replicated]:** `dim (J:A)_7 = 126 = dim J_7` (no non-Jacobian E1 direction).
  (The comment's "`dim S_13 = 27,132`" is a typo — that is the dimension in seven variables;
  `dim S_13 = 8568`, `dim J_13 = 6882` [computed], the latter as stated.)
* **Degree 7, E2 [reported]:** with the 6-dimensional gauge `F·S_1` removed, 120 parameters; split
  into the 8 parity sectors `(b0, b2, y) ∈ (Z/2)³` (all six `h1`-characters mixed within a sector).
  The two sectors with `(b0,b2) = (1,1)` (24 parameters) die by independent square columns mod
  `J_14`. In each of the other six sectors (16 parameters) 14 die the same way, leaving the family
  `D = Z_k(u·X1∂_{X1}F + v·Y1∂_{Y1}F)`; using the exact identity
  `X1²∂²_{X1}F − Y1²∂²_{Y1}F = 5(X1∂_{X1}F − Y1∂_{Y1}F) ∈ J_6` (because `(X1,Y1)` enter only through
  `X1⁶ + Y1⁶` and `P1`), the `u²` and `v²` columns coincide mod `J_14` while the `uv` column is
  independent, forcing `u² + v² = 0`, `uv = 0`, so `u = v = 0`. **Hence: no degree-7 denominator in
  any single parity sector (a fortiori any pure character sector).** Mixtures across parity sectors
  are explicitly left open ("120 total parameters").
* **Degree 8, E1 [replicated]:** `dim (J:A)_8 = 337 = dim J_8 + 1`; the extra direction lies in the
  invariant sector `χ = (0,0,0,0)` (`dim S_8[χ] = 87`, `dim J_8[χ] = 32`, colon 33 — all three
  numbers reproduced). Closed form:

  ```
  D8* = 4G_S² − G_P² = (2G_S − G_P)(2G_S + G_P) = ((∂_{X0}F)² − (∂_{Y0}F)²)/(X0² − Y0²),
  2X0Y0·D8* = (2Y0G_S − X0G_P)∂_{X0}F + (2X0G_S − Y0G_P)∂_{Y0}F,
  N* = ½X0²Y0²·(2Y0G_S − X0G_P, 0, 0, 2X0G_S − Y0G_P, 0, 0),   N*(F) = A·D8*.
  ```
  [replicated: `A·D8* ∈ J_14` mod `p`, and the factorization is immediate.]
  Geometric reading [inference]: `∂_{X0}F ± ∂_{Y0}F = (X0 ± Y0)(2G_S ± G_P)`, so `2G_S ± G_P` is the
  polar quartic of the `s0`-fixed point `p± = [1:0:0:±1:0:0]` with the hyperplane `X0 ± Y0 = 0`
  removed; `V(D8*) ∩ X` is the union of the two reduced apparent-contour divisors of the projections
  from `p+` and `p−`. These are the natural "contact" divisors attached to the exchange symmetry.
* **Degree 8, E2 on `D8*` [reported]:** in `(S_16/J_16)[(0,0,0,1)]` (dimension `1040 − 988 = 52`),
  `B_{D8*}(N*) ≢ 0` (12 nonzero remainder terms). As in §3.1 one checks that for
  `D = c·D8* + F·L` (`L ∈ S_2`) and `N = c·N* + A·L·Z/6 + K`, `B_D(N) ≡ c²·B_{D8*}(N*) (mod J_16)`
  (all cross terms land in `J` because `A·D8* ∈ J`, `F ∈ J`, `div K, K_ℓ ∈ J`) [proved], so this
  excludes `D ∈ C·D8* + F·S_2`. It does **not** exclude `D = c·D8* + Σ C_i∂_iF` with `C ∈ S_3^6`
  nonzero (316 effective parameters) — see §4.

### 3.3 §3 of the comment: symmetry audit and Reid–Tai ages on `Γ⁺`

[reported] `Γ = G_diag ⋊ ⟨s0, r1⟩` of order 96, `Γ⁺ = G_diag ⋊ ⟨s0r1⟩` of order 48 (replicated in
§1.4). The 8 reflections of age `1/2` all lie in `Γ⁻`. Over all 47 non-identity elements of `Γ⁺`,
146 projective eigenspace cases on `P⁵` and 90 tangent cases on `X`, the Reid–Tai ages are
`{1: 21, 3/2: 16, 2: 53, 5/2: 48, 3: 6, 4: 2}`, minimum 1. By Reid–Tai / Kollár–Larsen, `X/H` has
canonical singularities and Kodaira dimension 0 for every `H ⊆ Γ⁺`: **no uniruled quotient of `X` by a
residue-preserving subgroup** (so book §3.3 "finite quotients beyond reflections" cannot apply).

### 3.4 §4 of the comment: all-degree weighted-K3 Gauss–Manin obstruction for passive `(t1, t2)`

[reported; Hilbert-series and `F_sub` identity replicated] In the `⟨σ1,σ2⟩`-invariant coordinates
`(X0, Y0, U = X1², V = X2², t1 = Y1/X1, t2 = Y2/X2)` (constant Jacobian 4), `F_b` becomes

```
F_sub = X0⁶ + 5X0⁴Y0² + 5X0²Y0⁴ + Y0⁶ + (1+t1⁶)U³ + (1+t2⁶)V³
      + 5(X0²+Y0²)(t1²U² + t1t2·UV + t2²V²) + (5/7)·t1(1+2t2²)·X0Y0·UV ,
```

weighted-homogeneous of degree 6 in weights `(1,1,2,2)`; for generic `(t1,t2)`,
`S_{t1,t2} = {F_sub = 0} ⊂ P(1,1,2,2)` is a quasismooth K3 (`6 = 1+1+2+2`, `3A_1` points on
`X0 = Y0 = 0`), the double cover of the companion's cubic surface branched along `S² − 4T² = 0`.
A certificate whose field has `v(t1) = v(t2) = 0` restricts to each fiber as a fiberwise primitive,
forcing the residue of `AΩ_S/F_sub²` to vanish in `H²(S_η)`; then its Gauss–Manin derivatives
vanish, in particular their Kodaira–Spencer projections to `H^{0,2} ≅ R(F_sub)_12` (one-dimensional
socle), i.e. `A·∂_{t_i}F_sub ∈ J_12(F_sub)`. At `(t1,t2) = (2,3)`: weighted Hilbert function
`(1,2,5,8,12,14,16,14,12,8,5,2,1,0)` [replicated exactly], `dim S_12 = 140`, `dim J_12 = 139`
[replicated], and the normalized socle functional `λ` (`λ(U⁶) = 1`) has
`⟨λ, A∂_{t1}F_sub⟩ ≠ 0`, `⟨λ, A∂_{t2}F_sub⟩ ≠ 0` [reported]. **Conclusion: no rational certificate
at any denominator degree, over `C(t1,t2)` or any finite extension, keeps `(t1,t2)` passive.**
Geometric gloss [inference]: `s0` acts on each K3 fiber as a non-symplectic involution (it negates
`Ω_S`); the companion's class `A'Ω_S/F_sub²` is `s0`-invariant, hence lies in the invariant lattice
`H²(S)^{s0} = π*H²(cubic surface)`, which is algebraic — that is why `A'` was solvable with passive
`(t1,t2)`. The residual class is `s0`-anti-invariant and fqlx's computation shows it is not even
Gauss–Manin flat, so it is transcendental on the generic fiber.

### 3.5 §5 of the comment: scope summary and open directions (fqlx's own list)

Unconditionally ruled out over `C`: (1) all `deg D ≤ 6`; (2) all degree-7 `D` in a pure character
or parity sector; (3) `D8*` on its own (plus `F·S_2` gauge); (4) any certificate keeping
`(t1,t2)` passive; (5) uniruled `X/H`, `H ⊆ Γ⁺`.
Open per fqlx: degree 7 with mixed parity sectors (120 parameters); degree `≥ 8` coupling `D8*`
(and higher analogues in `(J:A)_d`) with active `(t1,t2)` components; a non-obvious correspondence or
fibration in which `t1, t2` vary along the fibers.

**My addition on the degree-7 mixed problem [computed].** With the gauge `(Z_kZ_5, e_5)` dropped,
the 120 unknowns `C ∈ S_2^6/(F·S_1)` and the 7260 product columns `H_ab = A·m_am_b·∂_{i_ai_b}F`
reduced modulo `J_14` (`dim S_14/J_14 = 1506`) give: no product vanishes mod `J_14`; the linearized
rank is **950**; and the iterated square-column elimination removes **no** unknown once all sectors are
allowed to mix (every square column `c_a²` is dependent on cross products `c_bc_d` of the same
character). The residual system is 1496 homogeneous quadrics in 120 unknowns (files
`d7_system.ms`, `d7_system.sing` in the scratch directory). Gröbner computation with msolve/Singular
mod `p`: **[pending at the time of writing — see the addendum at the end of this file]**. If the
projective zero set mod `p` is empty, then by upper semicontinuity of fiber dimension the zero set over
`Q̄` is empty too, and degree 7 is rigorously closed (the Hessian congruence being necessary).

### 3.6 Credits recorded in the thread

Problem formulation: Benjamin E. Diamond and the `hodge-computer/sextics` contributors; companion
solution: `HeadOfNaps` (with OpenAI Codex assistance, per its `solution.md`); Griffiths–Dolgachev
residue theory for weighted hypersurfaces; Reid–Tai / Kollár–Larsen.

---

## 4. What remains unexplored (concrete gaps)

1. **Degree 7, mixed sectors.** Exactly the 120-unknown quadratic system above. Never solved
   (pending here).
2. **Degree 8, the full E1-admissible family.** `(J:A)_8 = J_8 ⊕ C·D8*` (337 dims), gauge `F·S_2`
   (21): 316 unknowns `(C ∈ S_3^6, c)`. The mod-`J_16` necessary condition is quadratic:
   `H(C) + c·[D8*·div(AC) − A·C(D8*) − Σ_j C_j ∂_j(A·D8*)] + c²·B_{D8*}(N*) ≡ 0 (mod J_16)`
   (I derived the cross term using `N*(∂_jF) ≡ ∂_j(A·D8*) (mod J)`; it should be re-derived before
   use), 951 equations. Never run.
3. **Degree 9.** `dim (J:A)_9 = 760 = dim J_9 + 4` [computed]; the four non-Jacobian directions are
   `X1·D8*, X2·D8*, Y1·D8*, Y2·D8*` (sectors `(0,0,1,0), (0,1,0,0), (0,0,5,1), (0,1,0,1)`), while
   `X0·D8*, Y0·D8* ∈ J_9`. Nothing tested at degree 9.
4. **Degrees `≥ 10`.** E1 collapses: `dim (J:A)_10 = 2053 = 3003 − 950` (556 non-Jacobian directions
   spread over 24 sectors), `dim (J:A)_11 = 3702 = 4368 − 666` (maximal), and `(J:A)_d = S_d` for
   `d ≥ 13` since `R_{19} = 0` [computed]. From here on only full fixed-denominator linear solves or
   structural (geometric) choices of `D` are informative. No fixed-`D` solve of degree `≥ 7` is
   recorded anywhere for this instance.
5. **Natural denominator families never tested** (beyond `D = X0²+Y0²` on the companion): powers
   of the old quadric `S0^k` (why.md's own suggestion), `P0^k`, `δ^k`, polars/contours
   `(2G_S ± G_P)`, products `D8*·(form)`, `D8*²`, norms of linear pencils over finite algebras in the
   `ε`-sector. My E1 screen (§3.5 data, `cand.log`): `S0^4, P0^4, δ^4, S0²δ², S0²P0², G_S², G_P²,
   (2G_S±G_P)²` fail at degree 8; `S0^5, P0^5, δ^5` fail at degree 10; `S0^6, δ^6` fail at degree 12;
   every `D8*·(quadratic form)` passes E1 (trivially). Full solves for the `D8*`-multiples are
   pending (§6 addendum).
6. **Symmetry sectors never used constructively.** No construction has been attempted on the
   `Γ⁺`-quotient or on any cover of `X`; the `ε`-sector `(0,22,51,22,0)` is a level-2 sub-Hodge
   structure and the GHC predicts its support, but no divisor family with the needed `h^{2,0}`
   (non-normal or non-rational singularities, book §3.1) has been proposed.
7. **Nonlinear joint `(D,N)` searches**: none beyond fqlx's sector analysis.
8. **Deformation/continuation** from `t = 0` (`F*`, where the `X0`-reflection quotient serves
   `X0³Y0³` geometrically) to `t = 5/7`: never attempted; no explicit certificate at `t = 0` is in
   the catalogue either.
9. **Active-`(t1,t2)` transport** (fqlx's own suggestion): the Jacobian-transport machinery of §2.9
   applied to the genus-one fibration `U_t = 0` over `C(t0,t1,t2)` with moving-map corrections —
   not attempted; the cohomological obstruction (fibers with `p_g > 0` in both natural slicings) means
   a new idea about *which* surfaces carry the class is required first.

---

## 5. Candidate attack plan

Linear-system sizes at `q = 1` for a fixed `D` of degree `d`: unknowns `6·dim S_{d+1}` (numerators,
roughly halved by parity) `+ dim S_{2d−6}` (`Q`); equations `dim S_{2d+6}`. Numerically:

| `d` | numerator coeffs `6·dim S_{d+1}` | `dim S_{2d−6}` (`Q`) | equations `dim S_{2d+6}` |
|---|---|---|---|
| 7 | 12012 | 1287 | 53130 |
| 8 | 18018 | 3003 | 80730 |
| 9 | 26208 | 6188 | 118755 |
| 10 | 36288 | 11628 | 169911 |
| 12 | 51408 | 33649 | 324632 |

(Engine: `~/devel/sextics-work/hodge/fixed_d.py`, compressed modular solve + rational lift.)

**E1. Close degree 7 (running).** Gröbner basis mod `p` of the 1496 quadrics in 120 unknowns. If
only the trivial solution: rigorous exclusion of all degree-7 denominators (post as a reply to fqlx's
comment with the script). Chance that it closes: high (the linearized image already has rank 950 and
all pure sectors are dead); chance that it yields a certificate: negligible (the condition is only
necessary; a nontrivial solution would then need a full fixed-`D` solve).

**E2. Degree 8, full family.** Set up the 316-unknown quadratic system of §4.2 (after re-deriving the
cross term), reduce mod `J_16` sector-wise, Gröbner mod `p`. Size: 951 equations, 316 unknowns —
feasible for msolve. If solutions exist, each gives a candidate `D` for a fixed-`D` solve (sizes row
`d = 8`). Chance of a certificate: low (`< 5%`); value: closes degree 8 or produces the first
genuinely new candidate.

**E3. Fixed-`D` solves for structured degree-9/10 candidates.** `D = Z_i·D8*` (`i = 1,2,4,5`) at
degree 9; `D8*·{S0, P0, δ, X0², Y0², S1, P1, S2, P2, X_iY_j}` at degree 10; `D8*·G_S`, `D8*·δ²`,
`D8*·S0²` at degree 12; `D8*²` at degree 16 (`Q` of degree 26, `dim S_38 = 962598` equations —
heavy). Each is one linear solve (minutes to tens of minutes). Rationale: `V(D8*)` is the only
low-degree E1-admissible support and is geometrically meaningful (contour divisors of `p±`); if the
`ε`-class is supported on a divisor built from these contours, a low multiple would show it.
Chance: low (`~5%`); the contour divisor of a smooth hypersurface from a general point is smooth,
hence (book §3.1) cannot carry a `(3,1)` class, and `p±` are only mildly special.

**E4. Geometric `(3,1)`-support search.** Enumerate non-normal/non-rationally-singular divisors
compatible with the symmetry: (a) cones with vertex an `s0`-fixed line over surfaces of degree
`≥ 4`; (b) unions of hyperplane sections through `s0`-stable `P³`'s; (c) the preimages in `X` of
Noether–Lefschetz curves of the K3 fibration `X/⟨σ1,σ2⟩ → A²_{(t1,t2)}` (threefolds fibred in K3s on
which the restricted class becomes algebraic; their `(3,1)` Gysin image depends on their
singularities). For each candidate compute `h^{2,0}` of a resolution (or at least test E1 on the
defining equation) before any linear solve. Chance: low-to-moderate for finding *a* supporting
divisor family, but each test is cheap. This is the only route that attacks the actual obstruction
(the class is transcendental on all natural surface fibers).

**E5. Finite cover with a new reflection (template §2.7).** Seek a cyclic cover `Z → X` (adjoin a
root of a `Γ`-semi-invariant function, e.g. `T² = δ = X0² − Y0²` or `T⁶ = Y0⁶ + …`) on which the
pulled-back residue of `A` becomes *odd* under some reflection of `Z` whose quotient is Fano or
has algebraic `H²` of its surface fibers, and such that the deck character of the correction is
nontrivial (so trace recovers `A`). Concretely: with `T² = δ`, the exchange `s0` lifts to
`(X0,Y0,T) ↦ (Y0,X0,±iT)`; test whether `AΩ/F²` pulled back is odd under a lifted reflection. The
companion's double cover of the cubic surface is exactly `T² = S² − 4T'²` in the other direction, so a
further root of `δ` is the natural next step. Chance: moderate as an idea, low as a complete
construction (needs a new incidence seed and target identity, cf. §2.7's caveats).

**E6. Jacobian transport with active `(t0,t1,t2)`.** Use the genus-one fibration `U_t = 0` of §1.2
(target density `t0³u0³/8`, polynomial), compute the classical invariants `I_4, I_6` of the cubic
`U_t` (Fisher), the Weierstrass model over `C(t)`, and the first-variation source for the target;
then look for a *different* two-dimensional slicing of the base `A³_t` over which the Jacobian
elliptic surface is rational (`χ(O) = 1`) or K3 with the class algebraic. The audit script shows the
slicing "fix `(c1,c2)`, vary `t0`" gives `χ = 3`; the slicings "fix `(t1,t2)`" give the K3s of §3.4
(class transcendental). Other pencils through the base (e.g. fix `t0` and one ratio) have not been
computed. Chance: low (`~5%`); cost: moderate (sympy/Singular invariants in 3 parameters).

**E7. Continuation from `t = 0`.** On `F*`, construct an explicit certificate for `X0³Y0³` via the
`X0`-reflection quotient (book §3.3 incidence; or search low degrees directly — `F*` has many more
symmetries, so small denominators are plausible), then test the same `D_0` on `F_b` (one fixed-`D`
solve) and the first-order deformation system. Chance that `D_0` survives at `t = 5/7`: very low
(the obstruction is a character, not a coefficient size); value: a certificate on `F*` is itself a
useful data point for the catalogue's `X0`-reflection sector.

**E8. Averaging-ideal bookkeeping for covers.** Before any construction, compute for a proposed
auxiliary variety `Z → X` with finite group `G_Z` the averaging ideal of book §3.1 and check that the
`ε`-isotypic representation of `Γ` is reached; this is a finite character computation and rules out
most cover proposals in minutes (as the diagonal character test already rules out the companion's
quotient).

**E9. Numerical sanity on the Hodge side (optional).** Compute the period/Gauss–Manin behaviour of the
`ε`-sector along the one-parameter family `F_t` (Picard–Fuchs in `t`) to see whether the 22-dimensional
`(3,1)` part has any decomposable sub-structure (e.g. a CM or isotrivial piece) that could be
supported by a special surface family. This is exploratory and does not by itself give a certificate.

---

## 6. Feasibility: explicit witness vs. nonconstructive argument; Lean requirements

* **Explicit witness.** Rigorous lower bound `deg D ≥ 7` (fqlx), with degree 7 pending and degree 8
  essentially one-dimensional beyond `J`. If a polynomial certificate exists at degree `d ≤ 12`, the
  fixed-`D` linear algebra (sizes above) finds and lifts it in minutes once `D` is known, and the
  repository's witness pipeline (`hodge/leangen.py` → `Witness.lean`/`Solution.lean`, verified by
  `HodgeComputer.Verify` through kernel replay) handles witnesses of that size (a regenerated accepted
  witness passed the official sandbox in ~5 minutes). The bottleneck is purely the choice of `D`.
  Given that every "natural" low-degree support is excluded and that the class is transcendental on
  the natural surface fibers, I estimate the chance of a witness of degree `≤ 12` at well under 10%.
* **Nonconstructive argument.** The companion-style route (finite coefficient algebra, trace, legal
  clearing) requires an auxiliary geometry serving the `ε`-character; none is known. A purely
  cohomological existence proof would need divisor support of the `ε`-sector `(0,22,51,22,0)` — true
  conditionally on the generalized Hodge conjecture for `X`, but no supporting divisor family is
  identified, and the project's Lean contract requires an algebraic certificate in any case.
* **Lean.** Either (a) an explicit `(D, N, Q)` bound to `HodgeSubmission.F/A/q` through
  `HodgeComputer.HasCertificate` (as in every witness-based accepted solution), or (b) a
  finite-algebra construction in the style of `quadric-contact-broken` (56 modules, `maxHeartbeats 0`
  in the data-heavy files, matrix Bézout legality on a fiber), which presupposes a solved geometric
  transfer. (a) is routine once `D` exists; (b) is a multi-week formalization even with the shared
  `HodgeComputer.Constructions` library.

---

## 7. Files

* Instance: `/home/kejace/devel/sextics/data/quadric-contact-residual/{instance.rec, why.md, provenance.json, checks.json}`
* Companion: `/home/kejace/devel/sextics/data/quadric-contact-broken/{solution.md, acceptance.json, why.md, checks.json, *.lean}`
* Ancestor: `/home/kejace/devel/sextics/data/two-reflection-cubic-surface-d2/{solution.md, why.md}`
* Book: `/home/kejace/devel/sextics/docs/01-certificates.md` (key equation, Euler gauge),
  `02-explicit-constructions.md` (§2.4 finite-algebra line contact, §2.5 density boundary, §2.7 cover,
  §2.9 paired cubics), `03-geometric-support.md` (§3.1 averaging ideal, §3.3 reflections),
  `04-known-locus-and-frontier.md` (§4.2 paired-sign surface, §4.6 descent, §4.7 scope)
* Curation scripts: `/home/kejace/devel/sextics/scripts/audit-cubic-frontier.py`, `audit-jacobian-scope.py`
* Thread: `/tmp/hodge/discussions/quadric-contact-residual.md` (#125), `/tmp/hodge/discussions/quadric-contact-broken.md` (#119, empty)
* My scripts/logs: `/tmp/hodge/scratch/qcr-dossier/{qcr.py, sym.py, e1.py, hodge.py, trace.py, k3.py, cand.py, d7.py, fixedD.py, *.log, d7_system.ms, d7_system.sing}`

---

## Addendum: results of the pending computations

(To be filled in below as the background jobs finish.)
