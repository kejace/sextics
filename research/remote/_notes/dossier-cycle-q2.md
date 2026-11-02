# Research dossier: `cycle-q2` (OPEN) — "Cyclic coupling · pole order two"

Compiled 2026-11-01 from: `data/cycle-q2/{instance.rec, why.md, provenance.json}` (no `checks.json` exists for this instance), the complete GitHub discussion #61 (12 posts by `benediamond`, `fqlx`, `DiamondJM`, 2026-10-22 … 2026-10-28, saved at `/tmp/hodge/discussions/cycle-q2.md`), the accepted companion `data/cycle-three/{solution.md, proof.json, Solution.lean, acceptance.json}` (nonconstructive, no witness), the accepted `q=2` instances `holdout-q2` (nonconstructive-style finite algebra, no witness), `fermat-q2` and `block-q2` (explicit witnesses), the book chapters `docs/01` (certificates), `docs/02` (§2.6–§2.9), `docs/03` (§3.1, §3.2, §3.6), `docs/04` (§4.2, §4.6, §4.7), `formal/README.md`, the earlier notes in `/tmp/hodge/notes/accepted-group*.md`, and my own computations (scripts and logs in `/tmp/hodge/scratch/cq2/`, summarized in Appendix A).

Conventions. Variables `(X0,X1,X2,Y0,Y1,Y2)`; `ζ = ζ6 = e^{2πi/6}`, `ω = ζ²`. A diagonal element is written by its weight vector `w ∈ (Z/6)^6` (acting by `Z_i ↦ ζ^{w_i} Z_i`), modulo the global scalar, normalized by `w(X2)=0`. "Rigorous" = proved by a cited argument or an exact/modular computation whose logic gives a characteristic-zero conclusion; "evidence" = modular/numerical computation that does not by itself transfer; "claim" = asserted in the thread and not independently re-derived here. A failed search in a bounded ansatz is never a nonexistence result.

---

## 0. Status in one paragraph

Open. `q = 2`, `A = X0²·X1²·X2³·Y0²·Y1²·Y2` on the same cyclic sextic as the accepted `cycle-three`. Nobody has a certificate, an existence proof, or a Lean theorem. The thread (≈81 KB, three contributors, 12 posts, ~15 reproduction packets) has established, rigorously or with exact modular evidence: (i) the class `α = Res(AΩ/F³)` has a **nonzero `H^{3,1}` component** (fqlx; interval pairing bound), so no "`H^{2,2}`-only" construction and no support whose resolution has `H^{2,0} = 0` in the right character can serve it; (ii) the residue has the **nontrivial character `ω²` of the order-72 diagonal group `G`** (kernel `K` of order 24), so every quotient/average used by the accepted constructions on this variety (paired-sign sector theorem, cycle-three Jacobian transport) erases it; (iii) the obvious transplants are dead: cycle-three's rational-elliptic-surface model becomes, after the needed cube-root cover, an elliptic surface with `p_g = 2` (DiamondJM) on which the fibrewise source is **not** exact (benediamond, exact `H^{0,2}` test), the holdout-q2 cubic-surface/27-line method has no analogue (DiamondJM), every monomial surface fibration of `X/K` has `p_g > 0` and the whole non-trivial character sector carries holomorphic 2-forms (benediamond, five-point and lattice-width arguments), and no subgroup quotient of the visible order-216 group `G ⋊ Z/3` is uniruled (Reid–Tai); (iv) the problem has been reduced several times to a concrete **rational-current / Gauss–Manin integrability equation** over a two-dimensional base (rank-2 elliptic, rank-12 target sector on the triple cover, rank-18 proper-K3 module) — each equation is explicit, each has been searched in bounded rational ansätze over finite fields, and none has been solved or shown unsolvable. New facts established in this dossier (Appendix A): the exact Hodge numbers of the character sectors (`G`-sector of the target: `(0,6,26,6,0)`); the `q=2` E1 annihilator is zero through `deg D = 5` (so `deg D ≤ 4` is impossible and `deg D = 5` forces `D ∈ span(∂_iF)`), 26-dimensional at `deg D = 6`; a **validated sector-reduced fixed-denominator solver** (≈150× smaller than the dense system that DiamondJM could not run) that excludes all degree-5/6 polar-type denominators, the block-radial and `B6`/`K10` supports, the natural monomials, and `F_{X2}²` (results in A.5–A.6); and an elementary **`q`-raising identity**: a `q=1` certificate for `B = X0X1²Y0Y1²` on the one-parameter cyclic family `F_λ = F + (λ−½)X0X2³Y0Y2`, rational in `λ` and regular at `λ=½`, differentiates to a `q=2` certificate for `A` with denominator `D²` (§5, P2) — this "`q=1` sibling" problem has never been attacked at certificate level.

---

## 1. Exact target

### 1.1 Statement (`data/cycle-q2/instance.rec`)

```
q = 2
F = X0^6 + X0^3*X1*Y0*Y1/2 + X0*X2^3*Y0*Y2/2 + X1^6 + X1^3*X2*Y1*Y2/2 + X2^6 + Y0^6 + Y1^6 + Y2^6
A = X0^2*X1^2*X2^3*Y0^2*Y1^2*Y2
```

Certificate (`formal/README.md`, `docs/01`): nonzero homogeneous `D` with `gcd(D,F)=1` over `ℂ`, six numerators `N_i` homogeneous of degree `deg D + 6q − 5 = d + 7`, and the exact polynomial identity

```
E := A·D² − 2·D·Σ_i N_i ∂_iF + F·(D·Σ_i ∂_iN_i − Σ_i N_i ∂_iD) = F³·Q .
```

Degrees (`d = deg D`): `deg E = 2d + 12`, `deg Q = 2d − 6`; for `d ≤ 2` the certificate must be strict (`Q = 0`). Normalized form (ch. 1, Euler gauge): `E1: 2·N(F) = A·D`, `E2: D·div N − N(D) ∈ (F²)`. Equivalent rational statement: `F³ | A − 2·dF(v) + F·div(v)` for `v = N/D`. Parities (ch. 1): WLOG `D` is `Y`-even (replace `D` by `D·τ(D)`), `N_0,N_1,N_2` `Y`-odd, `N_3,N_4,N_5` `Y`-even, `Q` `Y`-odd. No denominator-degree bound and no requirement `Q = 0`.

### 1.2 Structure and sparse presentation

```
F = Σ_i (X_i^6 + Y_i^6) + ½·( X0^3·X1·Y0·Y1 + X1^3·X2·Y1·Y2 + X2^3·X0·Y2·Y0 )      (9 monomials; 2F ∈ Z[Z])
```

Support (exponent vectors): six pure sixth powers and the three oriented couplings `(3,1,0,1,1,0)`, `(0,3,1,0,1,1)`, `(1,0,3,1,0,1)`, i.e. `X_i³·X_{i+1}·Y_i·Y_{i+1}` for `i ∈ Z/3`. Each coupling has degree 4 in its "source" pair `(X_i,Y_i)` and degree 2 in its "target" pair, so `F` is even in each pair `(X_i,Y_i)` separately (paired signs) and `Y`-even. The chain sextic `F_H` of the holdout family is the same equation with the closing edge `X2³X0Y2Y0` removed (and coefficients 1 instead of ½ in the catalogue's `holdout-q2`); the family `F_λ = F + (λ − ½)·X0X2³Y0Y2` interpolates between the chain (`λ = 0`) and the cycle (`λ = ½`).

Target: `A = X0²X1²X2³Y0²Y1²Y2`, degree 12, `Y`-degree 5 (odd), every coordinate present; pair degrees `(4,4,4)`; the unequal exponents in the third pair (`X2³Y2` vs `X2²Y2²`) break the cyclic symmetry and single out one position of the cycle. Useful factorizations (all monomial, 122 of degree-6 type per fqlx):

| factorization | note |
|---|---|
| `A = M1·B`, `M1 = X0X2³Y0Y2` (closing coupling), `B = X0X1²Y0Y1²` | `M1` is a monomial of `F` and `G`-invariant ⇒ `B` has the same residue character as `A`; this is benediamond's "support in a family" identity `α = −½ ∇_λ Res(BΩ/F_λ²)` at `λ=½` and the basis of P2 in §5 |
| `A = m1·m2`, `m1 = X0X1X2Y0Y1Y2`, `m2 = X0X1X2²Y0Y1` | DiamondJM: `2AΩ/F³ = ∂_{s1}∂_{s2}(Ω/(F+s1m1+s2m2))|_0` — a *mixed second variation*; `m1` has trivial `G`-character (fqlx), `m2` has character `ω²` |
| `A = B1·C1`, `B1 = X1²Y0²Y1²` (G-invariant), `C1 = X0²X2³Y2` | fqlx: the only two factorizations with a `G`-invariant factor are `B1` and `B2 = M1` |

### 1.3 Degrees, parity, sector sizes

`A` is `Y`-odd and `F` is `Y`-even as the contract requires. The `Y`-parity involution `τ = (0,0,0,3,3,3)` lies in `G`, and `χ_A(τ) = ζ^{3·(3+3+2)} = 1`: the residue is `τ`-invariant (this is the statement that `AΩ` is `τ`-even). Sector-restricted linear-system sizes are in §1.11.

### 1.4 Smoothness (rigorous)

`F` is smooth: the Jacobian ideal of `2F` is zero-dimensional with `vdim = 5^6 = 15625` modulo `10007` (my Singular run, Appendix A.3); the Lean proof of cycle-three (`PhysicalField.lean`) proves irreducibility via a monic Eisenstein specialization, and the catalogue's acceptance of cycle-three presupposes smoothness. benediamond gives a direct **Newton-nondegeneracy** proof valid for every face: at a critical point of a face polynomial, multiplying the partials by the variables and summing absolute values gives `6S ≤ 3Σ|m|` with `S = Σ|Z_i|⁶` over active variables, while weighted AM–GM and the coupling exponent column sums `(4,4,4,2,2,2)` give `6|m| ≤ Σ e_i(m)|Z_i|⁶`, hence `6S ≤ 2S`, contradiction. Newton nondegeneracy is what makes the toric/Khovanskii arguments of §3.5 below legitimate.

### 1.5 Symmetry group (rigorous; Appendix A.1)

* **Projective diagonal group `G`.** `w ∈ (Z/6)^6/scalars` preserves `F` iff the three coupling weights vanish: `3w0+w1+w3+w4 ≡ 3w1+w2+w4+w5 ≡ 3w2+w0+w5+w3 ≡ 0 (mod 6)`. Order **72**, `G ≅ Z/2 × Z/6 × Z/6` (fqlx: `SNF(F/X2⁶) = diag(1,1,2,6,6)`; element orders: 7 of order 2, 8 of order 3, 56 of order 6 — consistent). Generators (checked to span `G`): `g = (2,0,0,2,4,2)`, `h = w0 = (2,1,0,3,2,1)`, the paired signs `s0 = (3,0,0,3,0,0)`, `s1 = (0,3,0,0,3,0)` (and `s2 = s0s1` modulo scalars), and `τ = (0,0,0,3,3,3)`.
* **Residue character of `A`.** `χ_A(w) = ζ^{Σ_i (e_i+1) w_i}` with `e+1 = (3,3,4,3,3,2)` (numerator weight plus the volume weight `Σ w_i`). Values on `G`: `{1, ζ², ζ⁴} = {1, ω, ω²}`; kernel `K` of **order 24**, `K ≅ (Z/2)² × Z/6` (element orders 7/2/14 of orders 2/3/6), containing `τ` and all three paired signs. `G/K ≅ Z/3` is generated by `g`, on which `χ_A(g) = ζ⁴ = ω²` (benediamond's "target has character `ω²`"). On `why.md`'s element `h = (2,1,0,3,2,1)`: numerator weight 17, volume 9, `χ_A(h) = ζ^{26} = ζ²` (the `why.md` statement). On `γ = (ζ3,ζ3,1,ζ3²,1,1) = (2,2,0,4,0,0)`: `χ_A(γ) = 1` (so `γ ∈ K`; the earlier note in `accepted-group1.md` that cycle-three's residue is moved by `γ` concerns the *other* target `X0³Y0³`).
* **Book subgroup.** `Γ_book = P × ⟨h⟩` (order 24, the group of the paired-sign sector theorem, docs/04 §4.2) is **not** `K` (`h ∉ K`); `K = P × ⟨k⟩` for an order-6 `k ∈ K`, e.g. `k = g+h = (4,1,0,5,0,3)`, which moves the passive coordinates (`c1 ↦ ζ^{-2}c1`, `c2 ↦ ζ²c2`), so `X/K` is *not* a fibrewise quotient of the cubic family (relevant to §3.5).
* **Non-diagonal symmetries.** Exhaustive check over the 720 variable permutations: the only permutations preserving the monomial support (up to diagonal rescaling) are the identity and the **cyclic pair permutation** `R: (X_i,Y_i) ↦ (X_{i+1},Y_{i+1})` and `R²`. `R` preserves `F` exactly (coefficients all ½) and normalizes `G`; `Γ = G ⋊ ⟨R⟩` has order **216** (benediamond; not claimed to be the full automorphism group). There are **no coordinate reflections** and no `X ↔ Y` exchanges. `A ↦ R(A) = X0³X1²X2²Y0Y1²Y2²`, `R²(A) = X0²X1³X2²Y0²Y1Y2²`: the three translates have **pairwise distinct `G`-characters** (`χ_{R(A)} − χ_A` takes the values `{0,2,4}` on `G`, so it is not identically 0). Consequence (chapter-1 character averaging): a certificate for the `R`-invariant sum `A + R(A) + R²(A)` **does** yield one for `A` alone, by projecting onto `χ_A` (denominator becomes a product of `G`-conjugates). `why.md`'s caveat "must still separate the other components" is therefore automatic once one works `G`-equivariantly.
* **`Y`-parity / general facts from chapter 1.** Since `F` contains all six pure powers, `L(F)` (lattice of exponent differences) has index 72 in the degree-0 lattice modulo scalars; `G = Hom(Z⁶/L(F), μ6)`. Numerator and denominator of a character component are both `G`-eigenvectors (support-coset rule); averaging may multiply `deg D` by up to 72.

### 1.6 Equivariant Hodge numbers of `H⁴_prim(X)` (rigorous; Appendix A.2)

Computed from the `G`-graded Koszul/Hilbert series of the Jacobian ring `R = S/J(F)` (`J` is generated by a regular sequence; `R_k` has dimensions `1, 426, 1751, 426, 1` in degrees `0, 6, 12, 18, 24`; Griffiths: numerators of degree `6k−6` at pole order `k` give `H^{5−k, k−1}_prim`). "Sector of `χ`" = classes on which the group acts by `χ`.

| group | sector | `(h^{4,0}, h^{3,1}, h^{2,2}, h^{1,3}, h^{0,4})` | rank |
|---|---|---|---|
| all of `H⁴_prim` | — | `(1, 426, 1751, 426, 1)` | 2605 |
| `⟨τ⟩` | `τ`-invariant (= `Y`-odd numerators) | `(0, 208, 866, 208, 0)` | 1282 |
| paired signs `P` (4) | `χ_A` (= invariant) | `(1, 126, 479, 126, 1)` | 733 |
| `Γ_book` (24) | invariant | `(0, 22, 74, 22, 0)` | 118 — matches docs/04 §4.2 |
| `Γ_book` (24) | `χ_A|Γ_book` | `(0, 18, 78, 18, 0)` | 114 |
| `K` (24) | `χ_A` (= invariant) | `(0, 18, 78, 18, 0)` | 114 |
| **`G` (72)** | **`χ_A` (target sector)** | **`(0, 6, 26, 6, 0)`** | **38** |
| `G` (72) | invariant | `(0, 6, 26, 6, 0)` | 38 |
| `G` (72) | the `h^{4,0}` sector (numerator weight `0`) | `(1, 6, 30, 8, 0)` | 45 |

So the target class lives in a rank-38 `G`-isotypic piece with a 6-dimensional `(3,1)` part; fqlx's checkpoint says the `(3,1)` component of `α` is nonzero there (§3.2). Among the 90 degree-12 `Y`-odd monomials of character `χ_A`, the Jacobian-ring image is the 26-dimensional `h^{2,2}` piece; `A ∉ J_{12}` (Singular). The `K`-sector is three `G`-sectors (`1, ω, ω²`); its Hodge numbers `(18,78,18)` are the `(1,10,1)`-per-surface data of DiamondJM/the book seen from the fourfold. (These numbers refine `why.md`; nothing in the thread computed them.)

### 1.7 Block / split structure (what geometry is visible)

1. **Paired-sign cubic.** On the chart `Y_i = r_i X_i`, `u_i = X_i²`: `F = (1+r0⁶)u0³ + (1+r1⁶)u1³ + (1+r2⁶)u2³ + ½(r0r1·u0²u1 + r1r2·u1²u2 + r2r0·u2²u0)` — a plane cubic in `u` with pure cubes and the three directed terms `u_i²u_{i+1}` (`why.md` of cycle-three, DiamondJM §1). With `v = r0`, `c_i = r_i/v`, `R = X2²`, `x = X0²/R`, `y = X1²/R`: `U = ax³+by³+cz³+dx²y+ey²z+fz²x`, `(a,b,c) = (1+v⁶, 1+c1⁶v⁶, 1+c2⁶v⁶)`, `(d,e,f) = ½v²(c1, c1c2, c2)` (cycle-three's chart). The target downstairs is `A = r0²r1²r2·u0²u1²u2²` (DiamondJM): a degree-6 source on a plane cubic, hence zero in the cubic's Jacobian ring — the class is carried entirely by base derivatives (second variation), not by a fibrewise source.
2. **Weighted K3 fibration** (benediamond, final comment). On `X1 = 1`, `λ = Y0/X0`, `μ = Y2/X2`, `U = X0²`, `V = X2²`, `D = Y1`, homogenizing coordinate `C`, weights `(2,2,1,1)`: `W = (1+λ⁶)U³ + (1+μ⁶)V³ + C⁶ + D⁶ + ½(λU²CD + μVC³D + λμUV²)`, a degree-6 K3 surface in `P(2,2,1,1)` over the `(λ,μ)`-plane; the source is `𝒞 = Res(λ²μ U²V²C²D² Ω_W / (4W³))`.
3. **Fermat-type cover** (holdout construction; lost here): adjoin `T⁶ = Y0⁶ + Y2⁶`; the closing term becomes `½αβX0X2³T²` with no matching `y²` term, so the `y ↔ T` exchange that produces the weighted quotient and its moving conic is absent (`why.md` of cycle-three; docs/04 §4.6 discussion).
4. **Deformation structure.** Two-parameter family `F + s1·X0X1X2Y0Y1Y2 + s2·X0X1X2²Y0Y1` (DiamondJM): both directions add an `xyz` term to the cubic with coefficients `g1 = v³c1c2`, `g2 = v²c1` — the `v²` is the character obstruction in coordinates.

### 1.8 Why the residue character kills the quotient constructions (rigorous)

Averaging a certificate over a finite group preserves certificates and projects onto the target character (docs/01). The paired-sign sector theorem (docs/04 §4.2) supports exactly `H⁴(X)^{Γ_book}`, rank 119 (118 primitive), and the cycle-three transport serves polynomial sources on the `Γ_book`-quotient Jacobian model — all of trivial `Γ_book`-character. `χ_A(h) = ζ² ≠ 1`, so the target is orthogonal to everything those mechanisms reach, "independently of the pole order" (`why.md`, docs/04 §4.7: "the surface-fibration theorem does not extend to a nontrivial character merely by replacing the numerator or increasing its pole order"). The same holds for the §2.8 complementary quotient (weights `(5,1,5)` on `Y`: the closing monomial has `Y`-weight `5+5 = 10 ≡ 4`, so `F` is not even invariant).

### 1.9 The `q=2` E1 annihilator (new; rigorous for the zero statements; Appendix A.3)

`E1` says `A·D ∈ J(F)_{12+d}`. Let `Ann_d = {[D] ∈ (S/J)_d : [A·D] = 0}`. Singular, `p = 10007`:

| `d` | `dim S_d` | `dim (S/J)_d` | `dim Ann_d([A])` (cycle-q2, q=2) | `dim Ann_d([X0³Y0³])` (cycle-three, q=1) |
|---|---|---|---|---|
| 0–4 | 1,6,21,56,126 | same | 0 | 0 |
| 5 | 252 | 246 | **0** | 0 |
| 6 | 462 | 426 | **26** | 0 |
| 7 | 792 | 666 | 420 | 0 |
| 8 | 1287 | 951 | 825 | 2 |
| 9 | 2002 | 1246 | 1190 | 52 |

Since rank over `F_p` ≤ rank over `Q`, a zero kernel mod `p` is a zero kernel over `Q` (and over `ℂ`, the system being defined over `Q`). Consequences:
* **`deg D ≤ 4` is impossible** (E1 alone; `J_d = 0` there).
* **`deg D = 5` forces `D ∈ span(F_{X0},F_{X1},F_{X2},F_{Y0},F_{Y1},F_{Y2})`** (with `D` `Y`-even: `D = c0F_{X0}+c1F_{X1}+c2F_{X2}`; with `D` additionally `G`-semi-invariant, a single `F_{X_i}`). All three `F_{X_i}` are excluded by the full solve (A.5). The non-semi-invariant combinations have not been tested (would need the dense `d=5` system: 37,254 unknowns × 80,730 equations, or a `K`-sector reduction).
* **`deg D = 6`**: at most 26 classes mod `J_6` (dimension 36), i.e. a ≤ 62-dimensional subspace of `S_6`; per `G`-character and `Y`-even the kernels are small (trivial sector: 6 among 12 monomials; most sectors 0–2; A.6).
* For comparison, the `q=1` cycle-three target has `Ann_d = 0` through `d = 7`, strengthening rbajaj5's `deg D ≥ 7` (E1 forces `D ∈ J` for `d ≤ 7`; rbajaj5's E2 analysis excluded `d ≤ 6`); the unknown cycle-three denominator thus satisfies `D ∈ J_7` or `deg D ≥ 8`.
* The `q=1` "siblings" (§5 P2): `Ann_d([B]) = Ann_d([B_bal]) = Ann_d([C1]) = 0` for `d ≤ 7` (`B = X0X1²Y0Y1²`, `B_bal = X0X1X2Y0Y1Y2`, `C1 = X0²X2³Y2`), so any `q=1` certificate for them has `deg D ≥ 7` with `D ∈ J_7` or `deg D ≥ 8` (A.3).

### 1.10 Geometric constraint on any denominator (rigorous given fqlx's checkpoint 1)

docs/03 §3.1: on a smooth sextic fourfold a divisor supporting a nonzero `(3,1)` class must have a component that is non-normal or has non-rational singularities (`Pic X = Z`, Kodaira vanishing, `H²(𝒪_D) = 0`). Since `α` has a nonzero `(3,1)` component (§3.2, fqlx), **every** legal `D` must satisfy this, and moreover the `(2,0)`-forms of the resolution must carry the character `χ_A` (restricted to the stabilizer of `D`). This excludes: monomial denominators (coordinate hyperplane sections are smooth sextic threefolds, their unions are non-normal only along smooth complete-intersection surfaces with `H¹ = 0`), products of smooth transversal divisors, and generic hypersurface sections. It favours: cones with vertex a line over surfaces with `p_g > 0` (quartics are excluded by E1; sextic cones in four variables appear in the `d = 6` E1 kernel, A.6), cones with vertex a plane over curves of genus ≥ 1, divisors singular along curves with non-rational transversal type, and polar/Hessian-type divisors whose singularities are not understood.

### 1.11 Sector-restricted fixed-`D` system sizes at `q = 2` (new; Appendix A.4)

For `D` `Y`-even and `G`-semi-invariant, the chapter-1 averaging argument shows that a certificate exists iff one exists with `N_i` in the character `χ_A + w_i + χ_D` (additive weights) with the stated parities and `Q` in character `χ_A + 2χ_D`, `Y`-odd. Sizes for `χ_D = 1`:

| `d` | unknowns `N` | unknowns `Q` | equations | dense system (N+Q unknowns vs equations) |
|---|---|---|---|---|
| 5 | 505 | 1 | 1,053 | 37,254 vs 80,730 |
| 6 | 729 | 7 | 1,678 | 51,870 vs 118,755 |
| 7 | 969 | 18 | 2,388 | 71,055 vs 169,911 |
| 8 | 1,305 | 41 | 3,323 | 96,027 vs 237,336 |
| 10 | ≈2,200 | 165 | ≈6,100 | ≈3.9·10⁵ vs ≈4.4·10⁵ |

The reduction factor is ≈144 (72 characters × 2 parities). This is why DiamondJM's dense degree-6 run exceeded 15 GB while the sector-restricted one takes seconds in pure Python (A.5).

---

## 2. Companion and related solved instances

### 2.1 `cycle-three` (same variety; `q = 1`, `A = X0³Y0³`; `hodgepodgeconjecture`, PR #165, accepted 2026-10-21; discussion #55)

**Method** (book §2.9 "cyclic paired-sign cubic"; solution.md): nonconstructive. Paired-sign descent to the ternary cubic `U` of §1.7 over the passive field `ℂ(v, c1, c2)`; classical covariants `H, Θ, J` with the syzygy `J² = 4Θ³ + 4s4ΘH⁴ − s6H⁶` (`s4 = ae²f+bdf²+cd²e`, `s6 = 27a²b²c²+4a²ce³+4ab²f³+6abcdef+4bc²d³−d²e²f²`), giving the Jacobian `η² = ξ³ + s4ξ − s6/4` with tangent identity `D_yξ/(2η) = 3/U_x`; polynomial Weierstrass model `V = Y² − X³ − a8X/8 − b256/256` in `(R,X,Y;T=Rv⁶;c1,c2)` with `a8, b256` of degrees 4, 6 in `(R,T)` — a **rational elliptic surface** over each `(c1,c2)` (`p_g = 0`), on which *every polynomial source has a rational residue primitive* (polynomial curve family of degrees `(2,4,6)`, 13 coefficient equations, `F5` seed at `(c1,c2,R(0),X(0)) = (2,3,1,3)`, arithmetic lifting, finite trace); first-variation (dual numbers) identity `ġ/g = (ȧ_Wξ+ḃ_W)/(2η²) − (g/3)D_y(ξ̇/(2η))`; Piola transport with residue-volume multiplier `F0·det M/(2Y) = −144Rv³`; moving-map correction field `W = (λhp1, 0, −λh1p3, λhp4, 0)`; normal-chart correction `C = −p0³p3⁴/h`; residue clearing and homogenization. **Denominator: not exhibited** (choice/trace); rbajaj5 excluded `deg D ≤ 6`; §1.9 shows E1 forces `D ∈ J_7` at `d = 7`. Lean: 60 files, 14.6k lines; the covariant syzygy is checked as 190 kernel-`decide` coefficient identities.

**Why it fails here** (`why.md`, DiamondJM, benediamond, docs/04 §4.7): (a) **character**: the model is the `Γ_book`-quotient; the target has character `ζ²` under `h` and `ω²` under `G/K`; the transported `q=2` source carries the factor `v² = T3/R` (DiamondJM §3; benediamond's factor `t`), a function on the **triple cover** `t = v⁶ = τ³`, not on the model; (b) **pole order**: `q=2` needs field components of degree 7 and divisibility by `F³`; the `q=1` first-variation identity is replaced by a *mixed second variation* with second-order moving-map and normal-derivative corrections (the `A4·V` term in benediamond's `N` "must be retained"); (c) **geometry of the cover**: pulling the Weierstrass family back along `(r:t) = (ρ³:τ³)` gives coefficients of degrees `(12,18)`, minimal, non-constant `j`, an elliptic surface `S3` with `χ = 3`, `p_g = 2` (DiamondJM) — the "every polynomial source has a primitive" theorem is false there, and benediamond proved the fibrewise source is **not** exact (§3.1). Also, by (a), Gauss–Manin derivatives of the cycle-three class in the coupling-coefficient directions cannot produce `α` (benediamond, 10-24).

### 2.2 `holdout-q2` (`q = 2` on the chain `F_H = ΣX_i⁶+ΣY_i⁶ + X0³X1Y0Y1 + X1³X2Y1Y2`, `A = X0³X1²X2²Y0³Y1²`; `HeadOfNaps`, accepted)

**Method**: strict (`Q = 0`) finite-algebra construction. Endpoint pair `(c,d) = (X2,Y2)` kept *active* through `S = c² + k²d²`, `T = kcd` (`k² + k + 1 = 0`), giving the cubic surface `G = S³ − 3ST² + 𝖺U³ + 𝖻V³ + 𝖼U²V + 𝖽V²T` over four coefficient parameters; a line on it via a degree-27 eliminant (`AdjoinRoot` algebra); `q=1` solver for *every quadratic source* (invertible 10×10 matrix on the good fiber `(2,2,1,1)`); explicit quintic reduction (`R = 27𝖺²𝖻 + 4𝖼³`, `H = R(R⁴ − 2916𝖺⁸𝖽⁶)`, `H = 515198970880` at the good fiber); pole-two correction with the Euler field `W2 = (N+GW1)/2 − hE/6`; horizontal lifts `C_U, C_V, C_S, C_T`; two endpoint targets `PΔ_k`, `PΔ_{k²}` combined over a product denominator. **Denominator**: cleared representation `6D0` with `D0 = det(SI − UM_e − VM_f)` of degree 27 in the cubic-surface variables, further multiplied by `x²a²` and the clearing exponent; legality via a matrix Bezout identity against `c⁶ + 5`. No witness file.

**Why it fails here** (DiamondJM, 10-23 §2; `accepted-group5.md`): the construction needs the endpoint pair to enter `F` only through `cd` and `c⁶ + d⁶`. In the cycle every pair enters through both `ZW` and `Z³W` (pair 2 through `X1³X2Y1Y2` and `X0X2³Y0Y2`); `c²` alone is not a rational function of `(S,T)` (`c²`, `k²d²` are the roots of `z² − Sz + T²`), so the closing coupling `X0X2³Y0Y2 = αx²c³d` is not polynomial in `(U,V,S,T)`. Keeping `(X2², X2Y2, Y2²)` gives `(2,3)` complete-intersection fibres in `P⁴` — K3 surfaces, no cubic surface, no 27-line algebra. The transferable part is only the generic pole-two Euler correction (fqlx's conditional Lean adapter, §3.4).

### 2.3 `fermat-q2` (`q = 2` on Fermat, same `A` as holdout-q2; explicit witness)

`D = E_I(F) = 6(X1⁶ + X2⁶)` for the subset `I = {X1,X2}`, `N = (A/2)E_I`, strict: works because `deg_I A + |I| = 4 + 2 = 6`. The subset-Euler trick needs `F` separated in `I`; the couplings kill it here (the solution itself warns the "shortcut does not survive a coupled sextic").

### 2.4 `block-q2` (`q = 2`, binary block `P(X0²,Y0²)` plus spectators; explicit witness, `deg D = 6`)

`D = P(u,w) = u³ + 2u²w − uw² + w³`, numerators of degree 13 found by a complete fixed-denominator solve; the target coefficient `−610/87` is forced by the binary cubic's discriminant. Relevance: a template for a *strict explicit* `q=2` witness with `deg D = 6` and a worked Lean file pair; the cyclic `F` has no binary block.

### 2.5 The chain's all-`q` support theorem (docs/04 §4.6, docs/03 §3.6) and the `q`-raising tools (docs/02 §2.6)

On the chain `F_H`, the Fermat cover `T⁶ = Y0⁶ + Y2⁶` with the deck action `d`, the exchange `σ: y ↔ T`, the physical diagonal `g` (fixed point of age 5/6) and `τ` satisfy `pd = p`, `pe = gp`, `pr = τp`, so `H⁴(X)^τ ⊆ N¹H⁴(X)`: **every odd target at every pole order** is supported, and §3.6 turns a support `D` into a certificate with denominator `D^q` by a complete fixed-denominator solve (numerators of degree `q(deg D + 6) − 5`). For the cycle the exchange is lost (§1.7.3), the §2.8 quotient is lost (§1.8), and no replacement has been found. §3.6 also records the pole-reduction identity `v = (B + F·w)/q`, `K_q(A;v) = (F/q)K_{q−1}(c;w)` for `q ≥ 5` (needs `A ∈ J`, false at `q = 2`: `A ∉ J_{12}`).

### 2.6 Siblings by structure (open): `off-paired`, `one-reflection-skew`, `quadric-contact-residual`, `cyclic-triple-residual`

All are "complementary-character" lanes (memory note): the accepted mechanism on the same or neighbouring variety serves one averaging-ideal character and the open target sits in another. `cycle-q2` is the only open `q = 2` lane and the only one where the target is a *second variation* of accepted-type classes.

---

## 3. Prior work in discussion #61, post by post

Packet URLs are listed in Appendix B. Post IDs are the GitHub `DC_…` node IDs from the saved thread.

### 3.1 `benediamond`, 2026-10-22 02:02 (`DC_kwDOVLOOJc4BI6Qv`): "A character-cover reduction and an explicit remaining connection problem" — partial

* **Source transport (claim with symbolic checks).** `A = B·H_def`, `B = X0X1X2Y0Y1Y2`, `H_def = X0X1X2²Y0Y1`; `−½∂_s(B/(F+sH_def)²)|_{s=0} = A/F³`. Cubic coordinates `v = Y0/X0`, `u = (Y1/X1)/v`, `w = (Y2/X2)/v`, `R = X2²`, `t = v²`, `T = Rt³`; `U = ax³+…` with `(a,b,c) = (1+t³, 1+u⁶t³, 1+w⁶t³)`, `(d,e,f) = t(u,uw,w)/2`; deformation `U + τxyz`, `τ = sut`. Fisher-normalized short Jacobian `η² = ξ³ + I4(τ)ξ + I6(τ)` with jets `I4'(0) = L = (9abc − def)/2`, `I4''(0) = 0`, `I6'(0) = M = abef² + acde² + bcd²f`, `I6''(0) = −I4(0)/6`. With `A4 = R⁴I4(0)`, `B6 = R⁶I6(0)`, `V = Y² − X³ − A4X − B6`, `P = R³LX + R⁵M`, the target transports modulo rational exact forms to `N·Ω_J/V³` with **`N = u²wtR/1728 · (A4·V − 12P²)`**, `Ω_J = dY∧dT∧du∧dw∧dR∧dX`; multiplier 3 (`P ↦ 𝒪(3P − H)`), density `−144Rv³`. The factor `t` carries the cubic character: the source lives on `t³ = T/R`.
* **Certified nonzero vertical remainder (rigorous, 5-adic with error bound).** Along `R(T) = r2T² + r1T + r0` the cover `t³ = T/R(T)` has genus 2 (`r0r2Δ ≠ 0`, `Δ = r1² − 4r0r2`); Hermite reduction `tQ(T)dT = d(tp(T)) + t(c0/T + c1/R(T))dT`; at the nonsingular 5-adic branch through the cycle-three seed `(u,w,r0,X(0)) = (2,3,1,3)`: `v5(c0) = v5(c1) = −2`, `25c0 ≡ 2`, `25c1 ≡ 4 (mod 5)`; 60- and 100-digit computations, independent jet implementation, `v5(error) ≥ N − 30`. Excludes vertical exactness on that branch and splitting only.
* **Next equation (open).** Hyperelliptic model `z² = 1 − 2bt³ + Δt⁶` (`b = r1`, `p = r0r2`); `ω`-eigenspace basis `θ0 = dt/z`, `θ1 = t³dt/z`; connection `M_D = [[bh−γ, h],[−Δh, −bh−4γ]]`, `γ = DΔ/(6Δ)`, `h = (2p·Db − b·Dp)/(6pΔ)`; normalize `r0 = 1`, `ξ = X(0)`, `c_i = C_i`; the rank-two inhomogeneous equation **`DH + M_D H = (C0,C1)ᵀ`** for a rational `H` over the incidence field (`Dξ = 1`, `Du = Dw = 0`). No solution found.
* Packet: `cycle-q2-research-2026-10-22.zip` (8 scripts, pinned deps).

### 3.2 `fqlx`, 2026-10-22 02:21 (`DC_kwDOVLOOJc4BI6Uj`): "Paused research checkpoint" — exclusions

1. **Nonzero `H^{3,1}` component (rigorous interval arithmetic, claim of correctness of the transport coefficients).** With `α = Res(AΩ/F³)`, `β = Res(X0X1²Y0Y1²Ω/F²)`, `h(v,w) = Q(v, w̄)`: `−51/200 < h(α,β)/((2π)⁴/6⁵) < −49/200`. Hence no `H^{2,2}`-only construction. (Note `β` is exactly the `q=1` sibling class of §5 P2; `A = M1·B`.)
2. **Supports `B6` and `K10` excluded, with arbitrary powers, each separately** (proper models have `H^{2,0} = 0`): `B6 = 3ΣX_i⁶ + ΣX_i³X_{i+1}Y_iY_{i+1} = ½ΣX_iF_{X_i}` (degree 6), `K10 = 2(Y2⁵∂_{Y1} − Y1⁵∂_{Y2})F` (degree 10). Mixed/additional supports not excluded. (My full solves for `B6`, `K10` agree: A.5.)
3. **Several complete Hecke traces vanish** (fixed-K3-base channel, complete scaling-even channel, differential-normalized target inverse).
4. **Normalized coefficient inverse**: degree 168, geometric monodromy `S168`, irreducible over the degree-2016 source field; the degree-nine target torsor unresolved.
5. **Full-normal current defect** `ε³𝒯⁴(f0)/6`, `f0 = 1/(λxy²)`; degree-14 pole support with 3312 rank-zero torus ODPs outside `K10 = 0`; rank-one tangencies and proper `H^{2,0}` open.
* Routes named: degree-three coset of the degree-nine torsor via a Hessian-triangle adapter; a globally regular eigenform on the adapter cover with nonzero source pairing; closing the full-normal current defect.
* Packet: gist `fqlx/74567249a07cb927125e1f785782c22d` (frozen research packet + decoder).

### 3.3 `fqlx`, 2026-10-22 03:14 (`DC_kwDOVLOOJc4BI6ec`): "New exact partials" — reply to 3.2

* Homogeneous part of 3.1's rank-two equation has **no nonzero algebraic horizontal section** (finite residues nilpotent, exponents `±1/3` at infinity) ⇒ a rational solution, if any, is unique; exact arithmetic circuit for `C`.
* For fixed `(u,w)` the retained cover `s = t³` of the Jacobian elliptic surface has **36 nodal fibres and `p_g = 2`** (no dominating-rational-curve transplant of cycle-three).
* Scoped Gauss–Manin test: a second-kind relative primitive would solve a rank-two equation with a 14-coefficient complete rational ansatz; at `(u,w) = (2,3)` mod 5: rank 14 vs augmented 15 ⇒ generically inconsistent (third-kind residues on moving divisors not excluded).
* One target-torsor coset excluded (Hessian-triangle line cubic stays irreducible over the minimal field; conditional on finite-flat/`S168` results).
* Packet: gist `fqlx/6c64c2750df88fbe9f4572aa8b734347`.

### 3.4 `fqlx`, 2026-10-23 03:03 (`DC_kwDOVLOOJc4BJAEn`): toric reduction, coupled consistency, exclusions — reply to 3.2

* **Divergence-free lift (exact).** In the chart `R = X2²`, `u = X0²/R`, `z = X1²/R`, `v`, `w = v²`, `c1`, `c2`, the five fields `∂_u,∂_z,∂_w,∂_{c1},∂_{c2}` scaled by `8X2⁶/v` become homogeneous degree-7 rational fields `D̃_c` on `P⁵` with `div = 0` and `D̃_c(X2²) = 0`; hence any rational solution of the chart equation `B − 2V(U) + U·div(V) = U³Q` lifts to the physical `(D,N,Q)` identity *with zero normal/Euler correction*. After two normal reductions (`V_u^{(0)} = B/(2g)`, `V_u^{(1)} = R1/(2g)`) the problem is a tangent divergence equation modulo `U³` on `{U2 = 0}`.
* **Toric structure (exact).** `SNF(F/X2⁶) = diag(1,1,2,6,6)` (`G ≅ Z/2×Z/6×Z/6`), `SNF(U) = diag(1,1,1,3,3)` (`(Z/3)²`); unimodular change to `(m0,m1,m2) = (c1wu²z, c1c2wz², c2wu)`; `U` depends on `(u,z)` only through `(u³,z³)`; index-9 invariant lattice `Λ` with explicit monomial relations; `g̃ = uU_u = 3(x0+y0) + m0 + ½m2`, `B = m0²m2/(8x0)`, remainders `R1 = B(m0+m2)/g̃²`, `R2 = (B/2)(3(m0+m2)²/g̃⁴ − (4m0+5m2)/g̃³)`.
* **Closed-form partial solutions.** Exact one-operator solution on the 6-variable subsystem `{y1 = m1 = 0}`: `F_{c1} = m2(3m0²(m0+m2) − 3m0²g2 + 3m0g2² − g2³)/(24x0m0g2³)`; on all 8 variables the pole-reduction primitive with `u_{c1} = m0 + m1 + 12y1` collapses the `g2^{−5..−2}` poles identically, leaving a simple pole `−m0²m2(u_{c1}² − 300u_{c1}y1 + 5400y1²)/(8x0u_{c1}⁴g2)`; in the free 8-variable ring all generators of `im(T_k) + U2R` vanish on the critical plane `𝒵_crit` where `RHS_base = −186624·y2(3y2−2y1)²(2y2−y1)² ≢ 0` (the `y2 = 0` slice is 100% solvable).
* Of the 122 monomial factorizations `A = B·C`, only `B1 = X1²Y0²Y1²` and `B2 = X0X2³Y0Y2` are `G`-invariant; both leave `C` with the nontrivial `μ3` character; `C_bal = X0X1X2Y0Y1Y2` has trivial character.
* **Coupled surface connection mod 1009 (evidence).** All 12 irreducible factors of the discriminant `D126(2,w)`: residue matrices of rank 1, no local bivariate residue obstruction. Univariate slices `u0 ∈ {2,3,5}`: uncoupled (`H_u = 0` or `H_w = 0`) inconsistent (`189 < 190`, `1053 < 1054`, `1701 < 1702`); the **coupled** slice equation is consistent (`rk(A) = 323 + 23·deg D`, `dim ker = 55 + deg D`); `(u−2)`-jets consistent to order 3 with nullity growth 50/51 per order; truncated polynomial-in-`(u−2)` solutions inconsistent ⇒ a global rational solution needs `u`-denominators.
* Conic-pencil divisors `g2(λ) = 0` (character `ζ3`, right sector, vertical Hermite collapse) vs `g0(λ) = 0` (wrong sector); origin-free internal cubic `Ψ_W` irreducible of degree 3 mod 103 (excludes descending non-coset components); universal conic primitives `B_a,B_b,B_c` with `2cB_a − 8bB_b − 88aB_c + 15(4ac−b²)I = 0`; 5-adic moving-amplitude defects `≡ 2 (mod 5)`.
* **Block-radial denominators excluded**: `D = G_rad,i^k`, `G_rad,i = (X_i∂_{X_i} + Y_i∂_{Y_i})(2F)`, over `Q` for `k = 1` (degree 6, certified integer duals) and mod 1009 for `k = 2` (degree 12). (My solver agrees for `k = 1`: A.5.) Fifteen five-border scalar-adjoint modules (≤315 unknowns/326 points) and five 12-D surface-connection trials (≤1152 unknowns/1260 equations): full column rank, augmented `+1`.
* Packet: gist `fqlx/d0cdc991a111dcd9b19ea37764f533b2`.

### 3.5 `DiamondJM`, 2026-10-23 01:26 (`DC_kwDOVLOOJc4BI_1U`): holdout-q2 does not transfer; exact character table

Items 1–4 as summarized in §1.7.1, §2.2, §1.5: paired-ratio cubic; no cubic surface; **character table by enumeration** (order 72; values `1, ζ6², ζ6⁴`; index-3 kernel of order 24 containing `τ` and the paired signs); dense degree-6 fixed-denominator test **not run** (>15 GB); the `q=1` sweep over ≈1300 natural degree-6 supports for the other open instances found nothing. Repo `DiamondJM/hodgeHoldout`, branch `claude/inspiring-turing-tx2pbu`, `cycleq2/structure.py`, `cycleq2/NOTES.md`, `sweep/RESULTS.md`.

### 3.6 `DiamondJM`, 2026-10-23 15:15 (`DC_kwDOVLOOJc4BJC7I`): the cycle-three Jacobian route at `q=2`

1. `2AΩ/F³ = ∂_{s1}∂_{s2}(Ω/F_s)|_0`, `F_s = F + s1m1 + s2m2` (exact); both directions add `xyz` with coefficients `g1 = v³c1c2`, `g2 = v²c1`.
2. **Invariants with an `xyz` term** (exact, normalization of cycle-three): `48s4 = 48(ae²f+bdf²+cd²e) + 216abc·m − 24def·m − m⁴`; `216s6 = 5832a²b²c² + 864a²ce³ + 864ab²f³ + 1296abcdef − 540abc·m³ − 864abef²m − 864acde²m + 72ae²fm² + 864bc²d³ − 864bcd²fm + 72bdf²m² + 72cd²em² − 216d²e²f² − 36def·m³ − m⁶` (point counts agree mod 101–109).
3. **Bookkeeping** (exact): `∂_{s1}∂_{s2}A_W = 0`, `∂_{s1}∂_{s2}B_W` carries `v⁵`; the Gauss–Manin formula `∂1∂2[dξ/η] = [3(A1ξ+B1)(A2ξ+B2)/(4η⁵) − (A12ξ+B12)/(2η³)]dξ` with the multiplier `−144Rv³` leaves a source carrying `v² = T3/R` — the triple cover.
4. **Triple-cover model** `S3`: coefficients of degrees `(12,18)` in `(ρ,τ)`, minimal, `χ = 3`, `p_g = 2` (vs `p_g = 0` for cycle-three's rational elliptic surfaces).
5. **Hodge-level test** (exact over `Q` at `(c1,c2) = (2,3),(3,5),(1,4)`): weighted Jacobian ring of `S3 ⊂ P(1,1,6,9)`: `dim R1 = 2 = p_g`, `dim R37 = 2 = h^{0,2}`, socle `R38`; target sector `(1,10,1)`, invariant `(0,8,0)`; transported source `Res(PΩ/G³) + Res(P'Ω/G²)` with `deg P = 33 ≤ 37`; **`P ∉ J37`** (nonzero fibrewise transcendental `(0,2)` part; no rational GM solution with denominators `t^jΔ^k`, `j ≤ 4`, `k ≤ 2`, numerator degree ≤ 40) **but `P ∈ J37 + ∂_{c1}G·S19 + ∂_{c2}G·S19`** — the base Kodaira–Spencer images span `R37`, so on the fourfold the `(0,2)` part can be absorbed by `c`-derivatives: *no obstruction at this level*; transcendental sub-variation ≥ `(1,2,1)`.
* Open: global two-variable rational integrability of `Σ` in the rank-12 target-sector GM system over `ℂ(c1,c2)`; caveats: re-derive `−144Rv³` with the `xyz` term; write the second-order moving-map correction.

### 3.7 `DiamondJM`, 2026-10-23 18:28 (`DC_kwDOVLOOJc4BJD3V`): the transcendental part is the whole target sector

Derivative spans of `s = Res(ρΩ/G)` in the target sector: dimensions `1, 3, 6, 10, 12, 12` for orders `≤ 0..5` at three points mod `1000003` ⇒ the rank-12 sector is entirely transcendental (no algebraic classes to discard). Engine checks: ring dimensions, exact flatness, agreement with item 5. Consequence: the class must vanish in de Rham cohomology of the full rank-12 system over `ℂ(c1,c2)` or be supported over a Noether–Lefschetz curve; connection entries along `c2 = 3` have `c1`-degree > 130 (frame poles) ⇒ direct rational reconstruction impractical; a good global frame or an NL search is needed. An unposted rank-4 computation (Jacobian ring missing `−x³`) is withdrawn. Scripts `gm_engine.py`, `gm_jets.py` (`python3 gm_jets.py 2 3 6`), `sigma_engine.py`, `degree_probe.py`.

### 3.8 `benediamond`, 2026-10-22 17:09 (`DC_kwDOVLOOJc4BI-IE`, reply to 3.1): "Two stronger obstructions" — rigorous (modular rank + good reduction)

1. **Fixed-parameter source not rationally exact.** For generic fixed `(u,w)`: `G = y² − x³ − a(t³)x − b(t³)`, `p = L(t³)x + M(t³)`, `ω = t³(a/y³ + 9p²/y⁵)dx∧dt`; proper model the quasi-smooth degree-18 surface `G_h = 0` in `P(1,1,6,9)`; pole-three numerator `N_h = t³z⁴(a_hG_h − 12p_h²)`, `deg N_h = 37 = 3·18 − 17`; weighted Griffiths residue: `H^{0,2}` component `[−12t³z⁴p_h²]` in degree 37 of the Jacobian ring is **nonzero** (ranks `68 → 69` mod 5 at `(u,w) = (2,3)`; character block exact rational ranks `23 → 24`; explicit 24-entry mod-5 functional). Deleting divisors cannot kill an `H^{0,2}` component ⇒ **no rational one-form primitive with arbitrary poles**, also after finite extensions (trace). Hence the endpoint equation `DH + M_DH = C` of 3.1 has no solution over its incidence field or any finite extension (normalized incidence evaluation is generically finite onto the surface retaining `t`; `(T,ξ)`-Jacobian is 1 at `T = 0`).
2. **`K10` normal-current repair obstructed at every order** (assuming fqlx's `h^{2,0}(K10) = 0`): `ω0 = dx∧dz/(λxy²)` on `1 + x⁶ + y⁶ + z⁶ = 0`, `λ⁵ = z³/x²`, is a nonzero holomorphic form on a smooth proper model (Fermat-30 cover); `W_2H²` of an open of `K10` is of type `(1,1)` (Deligne strictness), so no rational repair `Ω_rel − Ψdx∧dz` can restrict to `ω0`.
* Scope: closes the fixed-`(u,w)` endpoint and the prescribed `K10` extension; not the instance. Packet `cycle-q2-two-obstructions-2026-10-22.zip` (three proof notes, dependency-free modular checker, exact rational calculation).

### 3.9 `benediamond`, 2026-10-24 01:00 (`DC_kwDOVLOOJc4BJFBR`, reply to 3.1): "Constructive reductions"

1. **Three base variables.** On `X2 = 1`: `p = Y0/X0`, `q = Y1/X1`, `r = Y2`, `u = X0²`, `v = X1²`; `U_m = (1+p⁶)u³ + (1+q⁶)v³ + (1+r⁶) + ½pq·u²v + ½qr·v² + ½pr·u + m·uv`; `E_m: η² = ξ³ + a(m)ξ + b(m)`; `ω0 = dξ/(2η)`, `ω1 = ξω0`; `Σ = (p²q²r/24)∇_m²ω0|_{m=0} ∧ dp∧dq∧dr`; sufficient equation **`s = ∇_pR_p + ∇_qR_q + ∇_rR_r`, `R_• ∈ K²`, `K = ℂ(p,q,r)`** (rank-two proper connection; relative differentials with extra residues are another channel).
2. **Support differentiates in a family.** `M1 = X0X2³Y0Y2`, `B = X0X1²Y0Y1²`, `F_λ = F + (λ−½)M1`, `β_λ = Res(BΩ/F_λ²)`: **`α = −½∇_{∂λ}β_λ|_{λ=1/2}`**. A smooth proper algebraic family of Gysin maps containing `β_λ` near `λ = ½` would support `α` (relative Gysin commutes with GM); support at a single fibre is insufficient; cycle-three's class has a different character so its GM derivatives in coupling directions cannot produce `α` or `β_λ`.
3. Seek a family/correspondence whose variation matches the transported source (higher Gysin maps of varying horizontal curve families).

### 3.10 `benediamond`, 2026-10-23 17:15 (`DC_kwDOVLOOJc4BJDdB`): five-point obstruction — rigorous (convex geometry + Khovanskii)

Weight vector `β = (3,3,4,3,3,2)`, kernel `K` (24), `M_G ⊂ M_K` of index 3, `g = diag(ζ²,1,1,ζ²,ζ⁴,ζ²)` generates `G/K`, target character `ω²`. Dehomogenize at `X2 = 1`, `P = 6Δ5 ⊂ M_K`; a monomial map to a 2-torus is a saturated rank-2 passive lattice `N ⊂ M_K`, `L = M_K/N`, `Q = π(P)`; the residual cubic group acts fibrewise iff `N ⊂ M_G`. **Claim**: if `N ⊂ M_G`, `int(Q) ∩ L` contains a point of nontrivial `G/K`-character ⇒ the sum of the two nontrivial character sectors in the generic proper surface fibre is not of type `(1,1)`. Proof: the five `K`-invariant points `p = (1,0,1,1,2,1)` (colour 1), `q = (2,1,1,0,1,1)` (2), `r = (0,2,2,2,0,0)` (2), `s = (0,2,0,2,0,2)` (1), `t = (1,3,0,1,1,0)` (1) cannot all lie on `∂Q`: otherwise `Q` is a tetrahedron forcing two of the vertex images `a2, b1, b2` to coincide, killing one of `X2²/Y1², X2²/Y2², Y1²/Y2² ∈ M_K ∖ M_G`, contradicting `N ⊂ M_G`. Newton nondegeneracy of every face (§1.4) + Khovanskii (holomorphic 2-forms ↔ interior lattice points, with colour = character) gives actual forms. **Broader exclusion**: every rank-3 lattice image of `P` is nonhollow (lattice width of `P` in `M_K` is 3; of the 24 dual residue classes, 13 have width 3, 10 width 4, trivial class ≥ 6; Blanco–Haase–Hofmann–Santos reduce a hollow image to `M3,M5,M6,M9,M10`; exhaustive `24³` assignment enumeration: compatible assignments 164/12/76/78/74, with full convex hull **0**) ⇒ **no monomial surface fibration of `X/K` or `X/H` (`H ⊂ K`) has `p_g = 0`**. Not excluded: nonmonomial fibrations, correspondences with verified source action, absolute coupled-current arguments, cancellation of the specific class. Packet `cycle-q2-toric-2026-10-23.zip` (`cycle_kernel_hollow.py`, `cycle_character_projection_check.py`).

### 3.11 `benediamond`, 2026-10-23 17:16 (`DC_kwDOVLOOJc4BJDdU`): the order-216 group has no uniruled quotient — rigorous (Reid–Tai, Kollár–Larsen)

`Γ = G ⋊ ⟨R⟩`, order 216. For every `H ⊂ Γ`, `X/H` has canonical singularities and torsion canonical class; a resolution has Kodaira dimension 0, not uniruled. Diagonal elements: 110 fixed-eigenspace cases, all ages ≥ 1 (three of age exactly 1); for `gR`, `gR²` (eighteenth-root weights `s_X, s_X+6, s_X+12, s_Y, …`, normal weight `−6b`): ages `4/3:36, 5/3:36, 11/6:72, 2:36, 13/6:72`; no quasi-reflections. Standard-library Python check included. Says nothing about divisor-supported `H⁴` on a non-uniruled quotient, non-monomial fibrations, auxiliary covers, or correspondences.

### 3.12 `fqlx`, 2026-10-24 15:25 (`DC_kwDOVLOOJc4BJHNH`): fixed-current obstruction and scalar-curl repair

1. **Conditional Lean adapter** converting a homogeneous tangent primitive of the second normal residue into the `q=2` identity with legal denominator and degrees (Lean 4.33, three allowed axioms, kernel replay) — the *generic clearing* obligation is discharged; the source-closing primitive is the missing input.
2. **One fixed horizontal-lift/old-current system excluded** (proper K3 route): `Ba = r`, `a_μ − La_λ + Za = 0` with `B` the `17×34` proper source map; transport `z_ν = P_νz + h_ν`; curvature `𝒦 = (P_λ)_μ − (P_μ)_λ + [P_λ,P_μ]`; at `(λ,μ) = (2,3)` over `F101` determinants `79, 49, 94` ⇒ unique rational candidate; its `λ`-PDE residual `97 ≠ 0` (two routes agree in 17 coordinates). Excludes that fixed-current system in characteristic 0.
3. **Scalar-curl repair** `c = κ(C_src)/Y_{μλ}`, `f0 = c_λ`, `cθ_λ + f0ω = ∇_λ(cω)`, forcing `C_src − ∇_μ∇_λ(cω)`: three third-derivative certificates and the `λ`-compatibility test **pending**; a general scalar curl leaves 14 complementary compatibility identities.
* 2026-10-27 16:41 (`DC_kwDOVLOOJc4BJU92`): recovered packet `cycle-q2-recovered-2026-10-27.zip` (SHA-256 `5a08eb15…`), 13 arithmetic identity checks replay; the repair scripts are not standalone.

### 3.13 `benediamond`, 2026-10-27 17:21 (`DC_kwDOVLOOJc4BJVKG`, updated 10-28): proper-K3 rational currents and a scalar closure equation

Setup of §1.7.2. Over `K = Q(λ,μ)`, `v_ij = ∇_λ^i∇_μ^jω`; the 18 vectors indexed by `S = {i+j ≤ 4} ∪ {(0,5),(1,4),(2,3)}` form a basis of the proper residue module (determinant 16 mod 101 at `(2,3)` ⇒ generic invertibility); `G = (v_α)`, `a = G⁻¹𝒞`; explicit integration by parts gives **`𝒞 = ∇_λJ_λ + ∇_μJ_μ + Hω`, `H = Σ_α(−∂)^αa_α ≠ 0`** (formulas `I_λ(a;i,j)`, `I_μ(a;i,j)`); border set `ℬ = {(3,2),(4,1),(5,0),(0,6),(1,5),(2,4),(3,3)}`, `P_β = ∂^β − Σ_α b_{β,α}∂^α` with `P_βω = 0`; **closure equation `H + Σ_{β∈ℬ} P_β^*(z_β) = 0`, `z_β ∈ K`** (captures *all* rational currents in the proper module). Attempts: direct currents in all 18 classes, derivatives of the source, scalar-adjoint corrections with natural poles — no closing current (bounds recorded). Complementary ansatz `h_ν = μ^{−2}Σ_{j≤12}(μ−3)^j h_{ν,j}(λ)`: at `λ = 2` over `F1009` the frozen equation `Γ_λh_λ + (∂_μ + Γ_μ)h_μ = 𝒞` has a 26-dimensional solution space (a seed, not a family primitive; omits `∂_λh_λ`). Packet `cycle-q2-rational-current.zip`.

### 3.14 `benediamond`, 2026-10-27 17:21 (`DC_kwDOVLOOJc4BJVKO`, reply to 3.12): ten-class connection on fqlx's cubic

Depressed cubic `θ³ + pθ + q = 0`, trace-zero basis `e1 = θ`, `e2 = θ² + 2p/3`, `d` the monic normalization in `b` of `Δ = 4p³ + 27q²`, `K = Q(P,U,Z)`; representatives `e1, be1, b²e1, b³e1, e2/d, …, b⁵e2/d` (four affine elliptic classes incl. infinity + six logarithmic classes above the finite branch divisor, retaining homogeneous residue freedom); polynomial currents `J1 = (−3q/4, −p/4)ᵀ`, `J2 = (4p²/15, −3q/5)ᵀ`, `det[J1 J2] = Δ/60`; moving matrices `Γ_j`; equation **`Σ_{j=P,U,Z}(∂_j + Γ_j)u_j = C`, `u_j ∈ K^{10}`**, with reconstruction `v_j = Ru_j`, `v_b = v0 − ΣE_ju_j`. `F1009` searches with degree-3 polynomial factors and simple poles on four moving base planes: no closing current. Packets `cycle-q2-finite-branch-currents.zip`, `cycle-q2-rational-current.zip`.

### 3.15 `fqlx`, 2026-10-27 17:53 (`DC_kwDOVLOOJc4BJVUs`): handoff — strategy atlas, exact source bridge, two closure routes

* Public files: strategy atlas (40 files; interactive graph, five routes; 955 files inventoried / 676 assets recovered), current-and-remainder packet, general rational-integration packet, gist `fqlx/6e65179d28dfbd05fd0e94e831a52323` with `PUBLIC_README.md`, `PUBLICATION_REPLAY.json`.
* **Toric–K3 source bridge (exact, SymPy-verified):** on `X1 = C = 1`: `P = λd`, `u = d³`, `Z = μ`, `a = U_K/(V_Kd)`, `b = 1/(V_Kd)`; `G = (d³/V_K³)W1`, Jacobian `3d/V_K³`; with the Euler factor 2 reproduces the K3 numerator `λ²μU_K²V_K²d²/4`.
* **Mixed current leaves one scalar channel:** `Σ_{j=P,U,Z,b}(∂_j + B_j)u_j = −(ℛ, 0)ᵀ`, `ℛ = −PZbA2/(81T⁶)`; the identity still contains a nonzero `F²ℛ_phys` term that clearing/`Q` cannot absorb; six-sheet synthesis `ζξ³ + ξ + 1 = 0` (`ξ = (p/q_c)θ`, `ζ = q_c²/p³`) gives a scalar adjoint equation over the degree-6 extension — unsolved.
* **Proper-K3 presentation reduced from 7 to 6 adjoint inputs** (`P33` is generated: `P33 = ∂_μP32 + b_{32,31}P32 + …`), `ℬ6 = {(3,2),(4,1),(5,0),(0,6),(1,5),(2,4)}`.
* Corrections: `elliptic-boundary-repair.md` withdraws an earlier obstruction (indexing error); `rank13-connecting-map.md`: keep homogeneous residue freedom (compact class must lie in the image `H¹(M5) → H²(M8)`); formal jets are local evidence only; Macaulay2 attempts stopped during source normalization.
* Next deliverable named: a literal rational solution of the full cubic remainder equation or of the six-adjoint proper-K3 equation, with reconstruction, then legal `(D,N_i,Q)` and the protected Lean theorem.

### 3.16 Rigorous vs. screened — the ledger

**Rigorous (characteristic zero):** `G`, `K`, characters (§1.5); nonzero `(3,1)` component of `α` (interval bound; relies on the posted transport coefficients); fibrewise non-exactness of the transported source on the fixed-parameter surface, hence no solution of the 3.1 endpoint equation (modular rank with good reduction + Deligne/weighted residue theory); no monomial surface fibration of `X/K` with `p_g = 0`, and the whole nontrivial-character sector of fibrewise-quotient fibrations carries 2-forms; no uniruled quotient by subgroups of the 216-group; `K10`-repair obstruction at every order (conditional on `h^{2,0}(K10) = 0`); `B6`, `K10` supports excluded (`H^{2,0} = 0` of proper models, posted as checked but not re-audited here); block-radial `G_rad,i` excluded over `Q`; E1 degree bounds (§1.9) and the fixed-`D` exclusions of A.5 (each a zero-kernel statement mod `p` ⇒ over `ℂ`).
**Screened / bounded searches (no nonexistence content):** the 14-coefficient rank-two ansatz; `F1009`/`F101` current searches with stated pole/degree windows; the fixed horizontal-lift/old-current system (excludes one candidate current, not the module); mod-`p` jet consistency; DiamondJM's `t^jΔ^k` denominators with numerator degree ≤ 40.
**Claims not independently audited:** the exact transport formulas (3.1, 3.6, 3.13, 3.15 bridges) agree with each other pairwise in several places (fqlx's bridge reproduces benediamond's K3 numerator; DiamondJM's `s4,s6` reduce to cycle-three's) but no one has re-derived the full second-order moving-map correction and the `−144Rv³` multiplier with the `xyz` term (DiamondJM's own caveat).

---

## 4. What remains unexplored (concrete gaps)

1. **No sector-reduced fixed-denominator sweep at `q=2` had been run** before this dossier (DiamondJM's dense attempt was abandoned). Appendix A.5 starts it (degrees 5, 6, 10, ≈30 natural denominators); nothing of degree 7–9, nothing of degree 11–14, and no sweep over the E1-admissible families (§1.9) beyond the cone-type elements of A.6.
2. **The `q=1` sibling `(F_λ, B = X0X1²Y0Y1²)` has never been attacked at certificate level** — neither at `λ = ½` nor as a family. Its annihilator is zero through `d = 7` (A.3), so the search starts at `d = 7` with `D ∈ J_7` (18-dimensional `Y`-even semi-invariant families of the form `Σ ℓ_i F_i`) and at `d = 8`. The sizes at `q=1` are tiny (numerators of degree `d+1`).
3. **Non-semi-invariant degree-5 denominators** `c0F_{X0} + c1F_{X1} + c2F_{X2}` (forced by E1) were not tested (needs the dense `d=5` solve or a two-parameter determinantal computation).
4. **Cones** (vertex a line, four-variable sextics with `p_g = 10`; vertex a plane, three-variable forms) have not been systematically tested; the trivial-sector E1 kernel at `d = 6` contains such cones (A.6: e.g. two independent ones in `(X0,X1,Y0,Y1)`, two in `(X0,X2,Y0,Y2)`, two in `(X1,X2,Y1,Y2)`), i.e. exactly the `(3,1)`-compatible shape of §1.10; the full E2 test on them is running (A.6).
5. **Noether–Lefschetz search** on DiamondJM's triple-cover family `S3 → ℂ(c1,c2)`: curves in the `(c1,c2)`-plane where the rank-12 transcendental sector acquires algebraic classes (extra sections / isotriviality / CM). DiamondJM named it; nobody ran it. Likewise the proper-K3 family of 3.13: K3 surfaces with Picard jumps along curves in `(λ,μ)`.
6. **Nonlinear joint search** `(D, N, Q)` with `D` in the E1-admissible sector families (bilinear system; determinantal conditions on the `D`-parameters) — never run at any degree.
7. **Deformation/continuation from solved neighbours**: (a) the chain `F_H` (`λ = 0`), where every odd target at every `q` is supported through the Fermat cover — no attempt to continue the *support divisor* `D_λ` (as an algebraic family of hypersurfaces) from `λ = 0` to `λ = ½`; (b) cycle-three's own certificate as a function of the three coupling coefficients (it is defined over the 26-coefficient family) — its coupling-direction derivatives cannot give `α` (character), but the question whether `β_λ` is supported on cycle-three's unknown denominator for *some* `λ` was never asked.
8. **Non-monomial fibrations and correspondences** (explicitly left open by 3.10/3.11): e.g. quotients by `K = P × ⟨k⟩` with `k` moving the passive coordinates (§1.5), the Hessian/dual correspondences of the cubic family, the Prym/Kuga–Satake-type structures suggested by `p_g = 2` elliptic surfaces (the `K`-sector `(18,78,18)` and its `G`-sectors `(6,26,6)` are small enough for an explicit Hodge-structure study).
9. **The second-order transport identities themselves** (multiplier with `xyz`, second-order moving map) have not been written down and checked; fqlx's conditional Lean adapter is waiting for exactly that input.
10. **Primes/independence**: every modular exclusion in the thread is at a single prime (5, 101, 1009, 10007, 1000003 in different places); the zero-kernel results of A.5 were replayed at `p = 65521` for the completed cases (A.5), the others only at one prime.

---

## 5. Candidate attack plan

Degrees at `q=2`: `deg N_i = d+7`, `deg Q = 2d−6`, `deg E = 2d+12`. Sector sizes from §1.11.

**P1. Finish the sector-reduced fixed-`D` sweep (cheap, decisive for low degree).**
Idea: for each `Y`-even `G`-semi-invariant `D` of degree `d ≤ 12` in a "natural" list, solve the complete sector-restricted linear system mod two primes; any `SOLVABLE` hit is then lifted to `Q` by exact linear algebra and checked for `gcd(D,F) = 1`. Families: all `Z_jF_i` products (done, `d=6`), `F_iF_j` (`d=10`, partly done), `Z_j Z_k F_i F_l`, Hessian-type minors of `F`, the fqlx toric polynomials (`g2(λ)`, `g0(λ)`, `Prim_0..20` after homogenization), the E1-admissible families of §1.9 at `d = 6,7` sampled by random combinations and by their cone-type members (A.6), and cones over `G`-semi-invariant sextic surfaces in each of the 15 four-variable subsets and over plane curves in three-variable subsets. Size: `d=6`: 736 × 1,678 (seconds); `d=10`: ≈2,350 × 6,100 (minutes in pure Python; use a compiled mod-`p` kernel for `d ≥ 11`). Chance of a hit: low (the `(3,1)` constraint says the support must be geometrically special; the natural polar/toric candidates are all excluded so far), but the cone family is the one shape that is both E1-admissible and `(3,1)`-compatible, and it has never been looked at.

**P2. Reduce to the `q=1` sibling via the parameter derivative (new, exact).**
Identity (checked by hand, §A.7): let `F_λ = F + (λ−½)M1`, `M1 = X0X2³Y0Y2`, `B = X0X1²Y0Y1²` so `A = M1·B`. If `v_λ = N_λ/D_λ` is a `q=1` certificate for `(F_λ, B)` with coefficients in `Q(λ)` (or any field of functions of `λ`) regular at `λ = ½`, i.e. `B − v_λ(F_λ) + F_λ div v_λ = F_λ² Q_λ`, then differentiating in `λ` and putting
```
w := ½·( M1·v_λ − F_λ·v̇_λ )          (degree 7, denominator D_λ²)
```
gives `A − 2·w(F_λ) + F_λ·div w = −½·F_λ³·Q̇_λ`, i.e. **`w|_{λ=½}` is a `q=2` certificate for `(F, A)` with denominator `D_{½}²` and quotient `−½Q̇`** (provided `D_{½} ≠ 0` and `gcd(D_{½}, F) = 1`). This is the certificate-level form of benediamond's "support differentiates in a family" and of the book's `q`-raising tool; it needs no geometry. The sibling problem `(F_λ, q=1, B)` is a degree-6 `Y`-odd target of the *same* character `ω²` on the *same* one-parameter cyclic family; its class `β_λ` is exactly the one fqlx paired with `α`. Plan: (i) run P1 for `(F, q=1, B)` at `λ = ½` (sizes ≈ `d=7`: ≈100 unknowns; `d=8`: ≈150) over the E1-admissible families (`D ∈ J_7`, all of `S_8`-sector monomial supports, cones); (ii) if a certificate is found at `λ = ½`, re-solve with `λ` a formal parameter over `F_p(λ)` or by interpolation at several `λ` to get a rational family, then differentiate. Variants: any factorization `A = M'·B'` with `M'` `Y`-even gives a family `F + sM'` and a sibling `(F, q=1, B')` at `s = 0` (122 choices); with `M' = m2 = X0X1X2²Y0Y1` the sibling `B' = m1 = X0X1X2Y0Y1Y2` has **trivial `G`-character** (served, nonconstructively, by the cycle-three machinery on the 26-coefficient family), but the family `F + sm2` leaves the `Γ_book`-invariant space (`m2` has `h`-weight 2), so the family certificate must be built on the triple cover — this is precisely DiamondJM's second-variation setup; the `M1`-family keeps all symmetry and the `ω²` character. Chance: the sibling is strictly smaller than the `q=2` problem and has not been tried; success at `λ = ½` for some low-degree `D` is uncertain (same `(3,1)` obstruction, same character), but the degree thresholds (`d ≥ 7`) are now known and the systems are tiny.

**P3. Nonlinear joint search in the E1-admissible sector families (`d = 6, 7`).**
Idea: parametrize `D = Σ_k c_k D_k` over the E1 kernel of one `(χ_D, Y`-even`)` sector (dimensions 0–6 at `d=6`, A.6), form the sector system `M(c)·(N,Q) = A·D(c)²` whose matrix is linear in `c` and right-hand side quadratic; solvability is a determinantal condition on `c` (vanishing of all `(r+1)`-minors of the augmented matrix, `r` the generic rank). Compute over `F_p` by random evaluation to get the generic rank, then search the `c`-variety by Gröbner/elimination on a few minors or by Newton/alternating projection. Size: ≤ 6 parameters, 736 unknowns. Chance: low-to-moderate as a *closure* of P1 at `d ≤ 7` (if nothing is found, every `Y`-even semi-invariant `D` of degree ≤ 7 is excluded).

**P4. Solve or kill the proper-K3 closure equation `H + Σ_β P_β^*(z_β) = 0` (3.13/3.15).**
Idea: it is a linear inhomogeneous system of PDEs in 6 (after fqlx's reduction) rational unknowns `z_β ∈ Q(λ,μ)`; the natural next step is a *two-variable* rational ansatz with denominators restricted to the discriminant and branch divisors of the K3 family (which must be computed: the loci in `(λ,μ)` where `W` is singular or acquires Picard jumps), with Laurent windows in *both* variables (fqlx's `(u−2)`-jet result says `u`-denominators are necessary; benediamond's seed at `λ = 2` omits `∂_λ`). If instead one wants a nonexistence proof in the proper module, compute the de Rham cohomology of the rank-18 connection `∇` on the open base (Picard–Fuchs system of the K3 family, which is a Fuchsian system in two variables) and test whether the class of `𝒞` is zero in `H²_dR(base, ∇)`; this is a finite computation once the connection matrices and singular locus are explicit (the degree-130 frame poles DiamondJM saw must first be removed by a better frame, e.g. a basis adapted to the lattice polarization of the K3 fibres). Chance: moderate as a decision procedure; the equation is explicit and the ansätze so far were narrow.

**P5. Noether–Lefschetz / special-fibre search on `S3 → ℂ(c1,c2)` and on the K3 family.**
Idea (DiamondJM 3.7, never executed): the class dies on the generic fibre only if it is globally exact; otherwise it must be supported on a curve `C` in the base where the fibres acquire algebraic classes in the `ω`-sector (`H⁰(C, V_ω|_C) ≠ 0`). Compute, for the elliptic surface `S3(c1,c2)` with `p_g = 2`, the loci where the Mordell–Weil rank jumps or the fibration becomes isotrivial/CM (point counts mod several primes along candidate curves such as `c1c2 = const`, `c1 = ζ3^k c2`, `c1⁶ = c2⁶`, the `R`-fixed locus) and test the pairing of the transported source with the new algebraic classes. If a curve `C` is found with the right character, the certificate's denominator is the pullback of `C` (a polynomial in `c1, c2` = rational in `Y_i/X_i`, homogenized) times the fibre-level denominator — a *predicted denominator family* that P1 can then test directly. Chance: speculative but this is the one route that would explain a high-degree, structured denominator.

**P6. Continuation of the chain support (`λ = 0 → ½`).**
Idea: at `λ = 0` every odd target is supported via the Fermat cover; the §3.6 recipe produces, for the `q=1` sibling `B` on `F_0`, a denominator `D_0` (union of conjugate cover-pushforward divisors, `Y`-even). Compute `D_0` explicitly (fixed-`D` solve at `λ = 0`, which is the chain where the support theorem guarantees consistency), then test whether the *same shape* of divisor, deformed with `λ` (e.g. the pushforward of the moving conic family of the cover for `F_λ`, which still exists as a correspondence even though the exchange symmetry is lost), stays a support for small `λ` by solving the sector system over `F_p[λ]/(λ^k)` (jets) — exactly the "support in a family" criterion of 3.9 made computational. Chance: low; the thread's expectation (and `why.md`) is that the lost exchange is essential. But a jet computation would make that expectation a theorem for the specific divisor family.

**P7. The `K`-quotient with moving passive coordinates and non-monomial fibrations.**
Idea: `X/K` (order 24, `K ≠ Γ_book`) is a quotient by a group that fixes `χ_A`; `k ∈ K` of order 6 acts on `(v, c1, c2)` by `(ζ, ζ^{−2}, ζ²)`, so the quotient is not a fibrewise quotient of the cubic family — outside the scope of the five-point argument (fibrewise) but inside the lattice-width exclusion only for *monomial* projections. Study `X/K` by a non-monomial map: e.g. the Jacobian map of the cubic family itself (`(u,v) ↦ (ξ,η)` is non-monomial) composed with `K`-invariant functions of the base; or the projection to the `R`-quotient. Determine `p_g` of the resulting surface fibration or whether `X/K` has a rational-curve covering in the right sector. Chance: low (3.11 shows `X/K` is not uniruled; the only remaining hope is a surface fibration with `p_g = q = 0` that is non-monomial, which seems unlikely given 3.10's structure).

**P8. Prove nonexistence of low-degree certificates and publish the ledger.**
Even if no certificate is found, P1+P3 at `d ≤ 7` and the E1 bounds (§1.9) give a clean, reproducible statement "every `Y`-even `G`-semi-invariant denominator of degree ≤ 7 is excluded" (zero-kernel mod `p` ⇒ over `ℂ`), which nobody has for this instance (DiamondJM's `d = 6` test never ran). This is the analogue of rbajaj5's cycle-three bound and should be posted with scripts.

---

## 6. Feasibility assessment

**Explicit witness.** Everything known points to a denominator that is geometrically special (must carry `(2,0)`-forms in character `ω²` on its resolution) and of degree at least 6 (E1), realistically ≥ 8–10 by analogy with the siblings (`Ann` of the `q=1` siblings vanishes through `d = 7`; the accepted cycle-three denominator has `deg ≥ 7` and is unknown). At `d = 10–12` the sector-restricted systems have 2–4 thousand unknowns — trivial for a compiled solver — so *if* the right family of `D` is guessed (cones, pullbacks of NL curves, cover pushforwards), an explicit witness is computationally easy and would give a short Lean proof (`ring`-level identity plus coprimality by specialization, as in `block-q2`). The hard part is the guess; the geometric constraints (§1.10, §3.10–3.11) and the E1 families (§1.9) are the only guidance available, and the natural first candidates have now all been excluded (A.5–A.6).

**Nonconstructive argument.** The thread's three reductions (rank-2 elliptic connection over three base variables; rank-12 triple-cover sector over `ℂ(c1,c2)`; rank-18 proper-K3 module over `Q(λ,μ)`) are all of the form "decide whether an explicit class is exact for an explicit Gauss–Manin system over a two-dimensional base". That is a finite, if large, computation (P4) and is the most promising route to a *decision*. A positive answer would give a primitive as rational currents on the auxiliary family, which must then be transported back (second-order moving-map and normal corrections — not yet written) and cleared; fqlx's conditional Lean adapter covers the clearing. A negative answer in the proper module would not end the problem (boundary/third-kind currents, other supports) but would be strong evidence that the class is not of coniveau 1 in any of the visible geometries — which, under the generalized Hodge conjecture, cannot be the final answer, so it would point to a genuinely new correspondence.

**What a Lean proof would require.** (a) Witness route: `proof.json` with `witness.json` (`D`, six `N_i`, `Q`), `Solution.lean` proving `HodgeComputer.HasCertificate HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q` via the exported identity (kernel-checked `ring`/`decide` on polynomials of degree `2d+12 ≈ 32–36`, feasible as in `block-q2` but large: the identity has ~10⁵ monomials dense, far fewer in the sector) plus `IsRelPrime D F` over `ℂ` by a monic specialization argument (`AlgebraicSpecialization`). (b) Nonconstructive route: reuse `HoldoutCover.*`, `ResidueTransport.*` (cycle-three), `PolynomialCurves.ArithmeticLifting`, and fqlx's pole-three clearing adapter; new obligations are the second-order jet calculus (dual numbers of order 2 or two independent first-order jets, `TrivSqZeroExt` twice), the triple-cover function field (adjoining `v²`), the existence of a rational primitive on the auxiliary family (this is where P4/P5 must deliver an *algebraic* construction — a polynomial curve family with nonsingular seed, as in cycle-three, is the proven template), and the transport identities with the `xyz` term. Expect the cycle-three scale (10⁴ lines) or more. (c) Via P2: a `q=1` sibling certificate over `Q(λ)` (explicit or by a family version of the cycle-three construction) plus the ten-line differentiation identity of §A.7 — this would be by far the cheapest Lean route if the sibling is found explicitly.

---

## Appendix A. Computations performed for this dossier (`/tmp/hodge/scratch/cq2/`)

All Python run with `python3 -I`; Singular 4.x; primes `p = 10007` (and `65521` where stated). Zero-kernel / "no certificate" conclusions transfer to characteristic 0 because the systems are defined over `Z[½]` and `rank_{F_p} ≤ rank_Q`.

**A.1 `group_char.py`, `sectors2.py` — symmetry and characters (rigorous).** `|G| = 72`; `χ_A` values `{0,2,4}` (exponents of `ζ6`); `|K| = 24`; `τ, s0, s1, s0s1 ∈ K`; `χ_A(h) = 2` with `h = (2,1,0,3,2,1)` (numerator weight 17 + volume 9); `χ_A(g) = 4` for `g = (2,0,0,2,4,2)`; `γ = (2,2,0,4,0,0) ∈ K`; element orders `G`: `{1:1, 2:7, 3:8, 6:56}`, `K`: `{1:1, 2:7, 3:2, 6:14}`; generators `g, h, s0, s1, τ` span `G`; `Γ_book ≠ K`, `Γ_book ⊄ K`. Only three variable permutations preserve the support of `F` (identity, `R`, `R²`). `χ_{R(A)} ≠ χ_A`, `χ_{R²(A)} ≠ χ_A` as `G`-characters. 90 degree-12 monomials share the `G`-character of `A`.

**A.2 Equivariant Hodge numbers (rigorous).** Trace of `w ∈ G` on `R_k` is `[t^k] Π_i(1 − t⁵ζ^{−w_i}) / Π_i(1 − tζ^{w_i})` (Koszul resolution of the regular sequence `∂_iF`, `∂_iF` of weight `−w_i`); sector dimensions by character orthogonality; results in §1.6 (sum over all 72 characters reproduces `(1,426,1751,426,1)`; the `Γ_book`-invariant sector reproduces the book's `(22,74,22)`).

**A.3 `ann.sing`, `sib.sing` (Singular, `p = 10007`).** `dim std(jacob(2F)) = 0`, `vdim = 15625` (smooth). `A ∉ J_12`; `B, B_bal, C1 ∉ J_6`. Annihilator dimensions `dim (J:A)_d/J_d` as in §1.9; for the `q=1` siblings `Ann_d = 0` for `d ≤ 7`, `(3, 6, 4)` at `d = 8`, `(60, 60, 82)` at `d = 9` for `(B, B_bal, C1)`. (`(J:A)` has 27 minimal generators of degree 6 outside `J`, all `Y`-odd in the printed sample; the `Y`-even semi-invariant part is the one relevant for denominators and is computed in A.6.)

**A.4 `group_char.py` — sector sizes.** Table in §1.11 (`N_i` in character `χ_A + w_i + χ_D` with parities odd/odd/odd/even/even/even; `Q` in `χ_A + 2χ_D`, `Y`-odd).

**A.5 `fixedD_q2.py`, `batch10.py` — sector-reduced complete fixed-`D` tests at `q=2` (rigorous exclusions).** Solver validated on a synthetic certifiable target (`A_syn = 2M(F) − F div M` with random sector-correct `M`: reported `SOLVABLE` for `D = F_{X2}` and `D = X2²Y0²`, `validate.py`). Results, `p = 10007` (rank of the homogeneous system vs. augmented; "+1" means no certificate):

| `D` | `deg` | unknowns (`N`;`Q`) | rows | rank / aug | result |
|---|---|---|---|---|---|
| `F_{X0}`, `F_{X1}`, `F_{X2}` | 5 | 537–540;1–3 | ≈1,050–1,090 | 442/443, 427/428, 427/428 | excluded |
| `B6 = ΣX_iF_{X_i}` | 6 | 729;7 | 1,645 | 628/629 | excluded (agrees with fqlx) |
| `G_rad,i = (X_i∂_{X_i}+Y_i∂_{Y_i})F`, `i=0,1,2` | 6 | 729;7 | ≈1,580–1,610 | 609/610, 610/611, 611/612 | excluded (agrees with fqlx, `k=1`) |
| `X_jF_{X_i}` (9 cases) | 6 | 706–736 | ≈1,530–1,650 | all `+1` | excluded |
| `Y_jF_{Y_i}` (9 cases) | 6 | 705–736 | ≈1,570–1,640 | all `+1` | excluded |
| monomials `X0²X1²X2²`, `Y0²Y1²Y2²` | 6 | 729;7 | ≈1,250–1,280 | 519/520, 520/521 | excluded |
| `K10 = (Y2⁵∂_{Y1} − Y1⁵∂_{Y2})F` | 10 | 2,184;165 | 6,108 | 1965/1966 | excluded (agrees with fqlx) |
| `F_{X2}²` | 10 | 2,199;165 | 6,108 | 1690/1691 | excluded |
| `F_{X0}F_{X1}` | 10 | 2,196;165 | 6,108 | 1957/1958 | excluded |
| `F_{Y2}²`, `F_{Y0}F_{Y1}`, `F_{X2}F_{X0}` | 10 | — | — | *(still running when the dossier was closed; see `batch10_10007.log`)* | — |

Replay at `p = 65521` (`fixedD_65521.log`): **identical ranks** for all 27 degree-5/6 cases above (`F_{X_i}`, `B6`, `G_rad,i`, all 18 `Z_jF_i`, the two monomials). Note: `D` is never checked for coprimality by the solver; it only matters for `SOLVABLE` hits, of which there were none.

**A.6 `e1sector.py`, `cones.py` — E1-admissible `Y`-even semi-invariant denominators at `d = 6` (and `d = 7`) and the full test on them.** Per `G`-character sector at `d = 6` (`e1sector_d6.log`): trivial sector 12 monomials, kernel 6, with kernel elements supported on every four-variable subset (dimension 1 on most, 2 on `(X0,X1,Y0,Y1)`, `(X0,X2,Y0,Y2)`, `(X1,X2,Y1,Y2)`); other sectors have kernels of dimension 0–2 and none supported on ≤ 4 variables among those printed. `cones.py` extracts the kernel vectors (per sector and per 3-/4-variable subset) and runs the full sector solve on each basis element and on a random combination of the full kernel (`cones_d6.log`). Results at `d = 6`, `p = 10007`, **all "no certificate"** (rank vs. augmented rank always `+1`):
* Trivial sector (6-dim kernel): the four-variable kernel elements are *diagonal sextic cones* `X_a⁶ + c_1 Z_b⁶ + c_2 Z_c⁶ + c_3 Z_d⁶` (one per 4-subset; e.g. `X0⁶ − 3Y1⁶ − Y0⁶ + X1⁶`, `X0⁶ + 4000Y0⁶ + 2002X2⁶ + 6004X1⁶` mod `p`) and the three *coupling-plus-powers* elements `X_i³X_{i+1}Y_iY_{i+1} + 18Y_{i+1}⁶ − 6X_{i+1}⁶` (equivalently `X_i³X_{i+1}Y_iY_{i+1} + 6Y_{i+1}⁶ + 6Y_i⁶ − 6Y_{i+2}⁶` etc.); 18 cone-type elements, 6 basis elements and one random combination tested — all excluded (ranks 577–619 / `+1`; random combination 629/630).
* Sector `χ_D = (0,1,0,3,0)` (weights on `g,h,s0,s1,τ`): kernel `X0X1X2²Y0Y2 + 4X1X2⁵ + 3336X1⁴Y1Y2` — excluded. Sector `(0,2,0,0,0)`: 2-dim kernel (dense 8-monomial combinations) — both basis elements and a random combination excluded (631/632). Sector `(0,2,3,0,0)`: `X0X2³Y0² + 12Y0Y2⁵ + X1³X2Y0Y1` — excluded. Sector `(0,3,0,3,0)`: 1-dim — excluded. Sector `(0,4,0,0,0)`: 1-dim — excluded (642/643). Sector `(0,4,3,0,0)`: 2-dim (`X0³X1Y1Y2 + 12Y0⁵Y2 + X0X2³Y2²` and a dense one) — excluded. Sector `(0,5,0,3,0)`: `X0³X2Y0Y1 + 3X1²X2²Y1Y2 + 12X1⁵X2` — excluded. The remaining sectors of the `d = 6` run and the whole `d = 7` run (`e1sector_d7.log`: trivial sector kernel 10 among 14 monomials, cone-type elements on every 4-subset again; sector `(0,0,3,0,0)` kernel 7; sector `(0,1,0,3,0)` kernel 10) were still running when this dossier was closed; see `cones_d6.log`, `cones_d7.log`.

Interpretation: E1-admissible denominators of degree 6 are either (i) `J_6`-type (`Z_jF_i`, `F`), (ii) diagonal Fermat-type cones (vertex a line, `p_g = 10` surfaces — exactly the `(3,1)`-compatible shape of §1.10) or (iii) sparse 2–3-term semi-invariants built from one coupling and pure powers; **every one of them fails E2**. Any `SOLVABLE` line in a future run must be followed by an exact lift and a `gcd(D,F)` check (elements of `J_6`, in particular `F` itself and `Z_jF_i`, are in the E1 kernel and some are not legal). Summary statement supported by A.3 + A.5 + A.6 (mod `p`, hence over `ℂ` for each tested `D`): **no certificate with `deg D ≤ 4`; none with `deg D = 5` and `D` semi-invariant; none with `deg D = 6` for any of the tested E1-admissible semi-invariant `D` (which exhaust the trivial sector up to the 6-parameter family, of which only basis elements and one random point were tested — E2 is quadratic in `D`, so the family is not fully excluded; see P3).

**A.7 The `q`-raising identity (hand computation, §5 P2).** With `K¹(v) = B − v(F_λ) + F_λ div v = F_λ²Q`, differentiate in `λ` (`∂_λF_λ = M1`, `A = M1B`): `−v̇(F) − v(M1) + M1 div v + F div v̇ = 2FM1Q + F²Q̇`. For `w = c(F v̇ − M1 v)`: `K²(w) = A + c[−F v̇(F) + F² div v̇ + 2M1 v(F) − FM1 div v − F v(M1)] = A + c[2M1(v(F) − F div v) + 2F²M1Q + F³Q̇] = A + c[2A + F³Q̇]`, so `c = −½` gives `K²(w) = −½F³Q̇`. Degrees: `v` has degree 1, `F v̇` and `M1 v` degree 7 = `6·2 − 5`; denominator `D²`.

**A.8 Not done.** No exact lift over `Q` was needed (no hits). The second-order transport formulas of the thread were not re-derived. The `H^{3,1}` pairing bound of fqlx was not reproduced. The posted packets were not downloaded (links in Appendix B).

## Appendix B. Links

* Discussion #61: https://github.com/hodge-computer/sextics/discussions/61
* benediamond packets: `cycle-q2-research-2026-10-22.zip` https://github.com/user-attachments/files/33796383/cycle-q2-research-2026-10-22.zip ; `cycle-q2-two-obstructions-2026-10-22.zip` https://github.com/user-attachments/files/33834787/cycle-q2-two-obstructions-2026-10-22.zip ; `cycle-q2-toric-2026-10-23.zip` https://github.com/user-attachments/files/33883422/cycle-q2-toric-2026-10-23.zip ; `cycle-q2-rational-current.zip` https://github.com/user-attachments/files/34087110/cycle-q2-rational-current.zip ; `cycle-q2-finite-branch-currents.zip` https://github.com/user-attachments/files/34100104/cycle-q2-finite-branch-currents.zip
* fqlx gists/packets: https://gist.github.com/fqlx/74567249a07cb927125e1f785782c22d ; https://gist.github.com/fqlx/6c64c2750df88fbe9f4572aa8b734347 ; https://gist.github.com/fqlx/d0cdc991a111dcd9b19ea37764f533b2 ; recovered packet https://github.com/user-attachments/files/34031189/cycle-q2-recovered-2026-10-27.zip (SHA-256 `5a08eb15b134e73ab467cf7116631050b6a43f05b91cc98628d04394163da267`) ; strategy atlas https://github.com/user-attachments/files/34034680/q2-strategy-research-bundle.zip ; current-and-remainder https://github.com/user-attachments/files/34034681/cycle-q2-current-and-remainder.zip ; general rational integration https://github.com/user-attachments/files/34034682/cycle-q2-general-rational-integration.zip ; handoff gist https://gist.github.com/fqlx/6e65179d28dfbd05fd0e94e831a52323/5643a96366bee5f659c4f99dd1eca4605c941e21
* DiamondJM: https://github.com/DiamondJM/hodgeHoldout branch `claude/inspiring-turing-tx2pbu`, directory `cycleq2/` (`structure.py`, `model_check.py`, `triple_cover.py`, `invariants.py`, `inv_check.py`, `second_variation.py`, `gm.py`, `step2.py`, `solve_rational.py`, `jacobian_ring.py`, `ks_test.py`, `sectors.py`, `gm_engine.py`, `gm_jets.py`, `sigma_engine.py`, `degree_probe.py`, `NOTES.md`), `sweep/RESULTS.md`
* Accepted companions: `data/cycle-three/solution.md` (PR #165), `data/holdout-q2/solution.md`, `data/fermat-q2/`, `data/block-q2/`; book: `docs/02-explicit-constructions.md` §2.8–§2.9, `docs/03-geometric-support.md` §3.1, §3.2, §3.6, `docs/04-known-locus-and-frontier.md` §4.2, §4.6, §4.7
* References cited in the thread: Fisher, *genus one invariants* (https://www.dpmms.cam.ac.uk/~taf1000/papers/g1inv.pdf) and *Testing equivalence of ternary cubics*; Dolgachev, weighted residue theory; Deligne, Hodge II/III (strictness); Khovanskii, holomorphic forms on Newton-nondegenerate hypersurfaces; Blanco–Haase–Hofmann–Santos (arXiv:1607.00798) hollow lattice polytopes; Kollár–Larsen, quotients of Calabi–Yau varieties (arXiv:math/0701466).
