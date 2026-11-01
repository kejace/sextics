# Attack 1 — `polarized-k3-remainder` (q = 1, A = X0X1X2Y0Y1Y2, F = ½(f(U)+f(V)))

**STATUS: NO WITNESS FOUND.** Every tested denominator is inconsistent (see §2). The main result of
this session is structural (§1): under the Hodge-theoretic facts already established in the discussion
(`End_Hdg(T) = Q`), a certificate for this instance exists **if and only if** the generalized Hodge
conjecture (coniveau 1) holds for the exterior square Λ²T of the transcendental lattice of the K3 double
plane `w² = f`, unless the specific class η = [z0²z1²z2²] ∈ R_{f,6} ≅ H^{1,1}_prim(S) happens to be a
C-combination of algebraic classes (impossible if ρ = 1). No known theorem gives that coniveau statement for
this K3, so **no denominator of any degree can be found by search unless a new case of the generalized
Hodge conjecture is true**; and every denominator that works must cut X in a *singular* divisor whose
resolution carries a weight-2 Hodge structure of adjoint type so(T)(1).

Work directory: `/tmp/hodge/results/polarized-k3-remainder/` — `log.jsonl` (one JSON per tested
denominator; `"engine": 2` = computed with the engine after the sign fix of 2026-11-01), `uvring.py`
(Jacobian ring R_F = R_f ⊗ R_f in the (U,V) frame), `e1sweep*.py`, `pfsolve.py` (parity-free served-space
solver), `runD.py`/`batch.sh` (engine runs), `candidates_d8.json`, `run-*.out`, `pf-*.json`,
`target-image-*.json`, `k3proj.py`.

Conventions: `U = X + Y`, `V = X − Y`; then `F = ½(f(U) + f(V))` with

    f(z) = z0^6 + 5/7 z0^5 z1 + 5/7 z0^2 z1^3 z2 + 15/7 z0 z2^5 + z1^6 + 10/7 z1^5 z2 + z2^6

(verified exactly: `F((U+V)/2,(U−V)/2) − ½(f(U)+f(V)) = 0`), τ: Y ↦ −Y is the block swap U ↔ V, and

    Â(U,V) := A((U+V)/2,(U−V)/2) = 2^{-6} ∏_i (U_i² − V_i²)       (8 monomials).

Y-even/Y-odd in the (X,Y) frame = swap-symmetric/antisymmetric in the (U,V) frame. In the (U,V) frame
`J(F) = (∂_i f(U), ∂_i f(V))` and the Jacobian ring is the tensor product `R_F = R_f ⊗ R_f`,
`Hilb(R_f) = (1,3,6,10,15,18,19,18,15,10,6,3,1)` (verified mod p = 2^31−1; this also re-verifies
smoothness of f mod p, hence of X: the six partials form a regular sequence).

---

## 1. The obstruction and the geometric mechanism a certificate would need

### 1.1 Where the target lives

By the book's sector decomposition (Ch. 3A §3A.4–3A.5), via the rational map
Φ: S₆ × S₆ ⇢ X, S₆ = {u⁶ = f(z)} ⊂ P³,

    H⁴(X)_prim^τ = [primitive sixth-root pair, 441] ⊕ [primitive third-root pair, 441]
                 ⊕ [Λ²H¹(C)(−1), 190]  ⊕  [Λ²V₂, 210],     V₂ = H²(S₂)_prim, S₂ = {w² = f}.

The first three are supported on divisors by the book's companion constructions (existentially). The K3
sector receives exactly the monomials of V-degree ≡ 0 (mod 6); at q = 1 these are `B(U) − B(V)`,
`B(U) = Â(U,0) = U0²U1²U2²/64`. In Jacobian-ring terms the K3 component of the residue is

    r = η ⊗ ω − ω ⊗ η ∈ (Λ²V₂)^{3,1} = ω ∧ V₂^{1,1},   η := [z0²z1²z2²] ∈ R_{f,6} ≅ H^{1,1}_prim(S₂),

ω the K3 two-form (↔ 1 ∈ R_{f,0}). **Verified:** z0²z1²z2² is a *standard monomial* of R_{f,6}
(dim 19), so η ≠ 0 and the K3 component of the target is nonzero (independent confirmation of the
why.md claim). Since the residue of a pole-order-2 form lies in F³H⁴ and the τ-invariant part has no
(4,0), the class is exactly its (3,1) part, i.e. the Jacobian-ring class is the whole story.

### 1.2 Equivalence with the generalized Hodge conjecture for Λ²T

Write `V₂ = T ⊕ NS_prim` (rational Hodge structures), so

    Λ²V₂ = Λ²T ⊕ (T ⊗ NS_prim) ⊕ Λ²NS_prim.

Facts used: (a) by the discussion (1337isnot1337, Frobenius/Galois argument) `End_Hdg(T) = Q`, hence
by Zarhin `MT(T) = SO(T)` and Λ²T = so(T) is an **irreducible** Hodge structure (adjoint representation,
rank T ≥ 19 > 4); (b) `Hom_HS(Λ²T, T) = Hom_HS(Λ²T, Q(−2)) = 0` (dimension); (c) the sector projectors are
algebraic correspondences (book §3A.4), and algebraic correspondences preserve the coniveau filtration N¹.

**Proposition.** Let Y = X ∩ V(D) for a legal D, N_Y ⊂ H⁴(X,Q) the image of H²(Ỹ,Q)(−1) → H⁴(X,Q)
(Ỹ = disjoint union of resolutions of the components). If the K3 component r of the target has nonzero
projection to Λ²T, then

    r ∈ N_Y ⊗ C   ⟹   Λ²T ⊂ π_{K3}(N_Y) ⊂ N¹H⁴(X,Q).

*Proof.* W := smallest Q-Hodge substructure of Λ²V₂ with r ∈ W_C. Its projection to Λ²T is a nonzero
sub-HS, hence all of Λ²T; by (b) a sub-HS surjecting onto Λ²T contains it. Since π_{K3}(N_Y) ∋ r is a
Q-Hodge substructure, W ⊂ π_{K3}(N_Y). ∎

Transporting along Φ and the pullback S₆ → S₂ (algebraic correspondences, injective on this sector),
Λ²T ⊂ N¹H⁴(X) ⟺ Λ²T ⊂ N¹H⁴(S₂ × S₂). Conversely the book's prescribed-support theorem turns a
support into a certificate. Hence:

> **Theorem (conditional on `End_Hdg(T) = Q` and on η ∉ NS_prim ⊗ C).** `HasCertificate` for this
> instance ⟺ the generalized Hodge conjecture GHC(S₂ × S₂, 4, 1) holds for Λ²T(S₂), i.e. Λ²T is
> supported on a divisor of S₂ × S₂.

η ∉ NS_prim ⊗ C is automatic if ρ(S₂) = 1 (then NS_prim = 0); the discussion's bound is 1 ≤ ρ ≤ 3, and
ρ = 1 is what one expects for this sextic. If ρ ∈ {2,3}, the only loophole is that the rational class
η = [z0²z1²z2²] is a C-combination of the ≤ 2 primitive algebraic classes — a transcendental coincidence
that can only be tested with period computations (Lairez–Sertöz style), not by linear algebra.

The GHC for Λ²T of a K3 is known only when the Kuga–Satake correspondence is algebraic (book §3A.6:
exterior-power theorem), i.e. for Kummer/Shioda–Inose type (ρ ≥ 17), Paranjape's six-line double planes
(ρ ≥ 16), the Bolognesi–Laterveer order-three family and its relatives (non-symplectic automorphism ⟹
`End_Hdg(T) ⊋ Q`), RM/CM cases — all excluded here by ρ ≤ 3 and `End_Hdg(T) = Q`. So the instance is
**equivalent to an open case of the generalized Hodge conjecture**; a search can only succeed if that
case is true, and any witness would be a theorem about this K3.

### 1.3 Two consequences for searching (used to prune the experiment menu)

**(i) Dichotomy.** For every legal D, the K3-sector projection of the space of (3,1)-targets served by D
is either ⊆ ω ⊗ NS_prim (dimension ≤ ρ − 1 ≤ 2; = 0 if ρ = 1) or ⊇ (Λ²T)^{3,1} (dimension ≥ 17).
There is no partial progress inside the K3 sector: a D serves either (essentially) none of it or all of it.
Computationally: `K3_projection_dim` of `target_image(D)` (my `k3proj.py`/`pfsolve.py`) is the right
diagnostic; so far it is 0 for every D tested.

**(ii) Sum principle for reducible denominators.** If D = D₁D₂ with coprime factors then
served(D) = served(D₁) + served(D₂) (modulo J(F)). Reason: at q = 1, served(D) = N_{V(D)} ∩ F³ by the
prescribed-support theorem, and N_{Y₁∪Y₂} = N_{Y₁} + N_{Y₂} because the image of H⁴_{Y}(X) → H⁴(X) is
pure of weight 4 and equals the image of ⊕ H²(Ỹ_i)(−1) (Deligne: W_{−4} of H₄^{BM}(Y) is the image of
the resolution of its components). Combined with the discussion's exclusion of all degrees ≤ 7,
**every product of factors of degree ≤ 7 is excluded**: f(U)·h, (f(U)−f(V))·h, P_i·(cubic), all
products of linear forms, all "cone" products. **Only denominators with an irreducible factor of degree
≥ 8 can work.** (Computational check of the principle: §2.2.)

### 1.4 What the divisor must look like

A *smooth* hypersurface section Y of X has H²(Y) ≅ H²(X) = Q (Lefschetz), so its Gysin image is Q·h².
Hence the divisor V(D) ∩ X must be singular, and the needed classes come from the singularities:

* vertex-a-line cones D ∈ k[l₁,…,l₄]: Y is singular at the 6 points of L ∩ X; blowing up, the exceptional
  surface is the base surface S' = V(D) ⊂ P³ (p_g(S') = C(d−1,3) = 35 for d = 8). A (2,0)-form on S'
  lifts to Ỹ iff its image in W₂H³(Y) vanishes, and then pushes to H^{3,1}(X). The Hodge structure the
  Gysin image can carry is a sub-quotient of ⊕_{6 points} H²(S')(−1); to reach the K3 sector it would have
  to contain so(T)(1) — i.e. S' itself would have to be a "Kuga–Satake surface" for this K3 (the surface
  Σ ⊂ A × A of the book's exterior-power theorem), which no octic in P³ written down by hand is.
* vertex-a-plane cones D ∈ k[l₁,l₂,l₃]: Y is singular along the plane sextic Γ = Π ∩ X; the exceptional
  divisor is a C'-bundle over Γ (C' = base curve), h^{2,0} = g(Γ)g(C') = 10·21 = 210 for octics; the Hodge
  structure is H¹(Γ) ⊗ H¹(C')(−1), of "product of curves" type, which contains so(T)(1) only if J(Γ) or
  J(C') involves the Kuga–Satake abelian variety (dimension 2^{rank T − 2}) — impossible for plane curves.
  Note Π = {Y = 0} gives Γ = C = {f = 0}, the branch sextic itself (D = g(Y) for ternary octics g); we
  tested this family (§2) — inconsistent, as predicted.
* The divisors on S₆ × S₆ that a D of degree d defines are the μ₆-invariant (d,d)-divisors
  {Σ c_{ab} z^a z'^b u^{|b|} v^{|a|} = 0}; the joining divisor is D = f(U) (the exceptional divisor of the
  blow-up along C × C), hyperplane sections U_i are {ℓ(z) = 0} × S₆, bilinear D = UᵀMV are the
  (1,1)-correspondences {zᵀMz' = 0}, cubics c(U) − κc'(V) are {u³c'(z') = κ v³c(z)}, etc. Graph-type
  divisors (graphs of automorphisms/correspondences) have Gysin image in H²⊗H⁴ ⊕ H⁴⊗H², families of
  curves S₆ × C' give T ⊗ NS: neither reaches Λ²T. Only a divisor whose resolution has "new" (2,0)-forms of
  adjoint type does.

### 1.5 E1 (Jacobian) structure — exact facts

Multiplication by [Â] on R_F = R_f ⊗ R_f (mod p = 2^31 − 1), dimensions
`dim R_F,d = 1,6,21,56,126,246,426,666,951,1246,1506,1686,1751,…` (palindromic, socle degree 24):

| d | dim R_d | dim R_{d+6} | rank | kernel | sym / anti |
|---|---|---|---|---|---|
| 0–6 | 1,6,21,56,126,246,426 | 426,…,1751 | = dim R_d | 0 | — |
| 7 | 666 | 1686 | 666 | 0 | — |
| 8 | 951 | 1506 | 951 | 0 | — |
| 9 | 1246 | 1246 | 1246 | 0 | — |
| 10 | 1506 | 951 | 951 | 555 | 294 / 261 |
| 11 | 1686 | 666 | 666 | 1020 | 510 / 510 |
| 12 | 1751 | 426 | 426 | 1325 | 677 / 648 |

So [Â] has **maximal rank in every degree**: injective for d ≤ 9, surjective for d ≥ 10. Consequences:
every E1-admissible denominator of degree ≤ 9 lies in J(F) (extends the discussion's "injectivity in the
screened degrees" through d = 9 and justifies the Koszul parametrization D = Σ a_i ∂_iF, N = Aa + Koszul
up to degree 9); from d = 10 on E1 is no constraint beyond dimension counting (the symmetric E1 kernel
at d = 10 is 294-dimensional modulo J(F)_{10}). Everything is decided by E2.

(Sections 2–3 below are filled in from the experiment log.)

## 4. Mathematical ideas developed, partial constructions, what would be needed to finish

4.1 **The sharp reformulation (§1.2)** is the main output: the instance is equivalent to GHC(coniveau 1)
for Λ²T of the degree-2 K3 `w² = f`, with the single loophole η ∈ NS_prim ⊗ C (needs ρ ≥ 2).
Consequently the only two routes to a solve are:
 (A) prove Λ²T ⊂ N¹H⁴(S₂×S₂) for this K3 — e.g. by proving algebraicity of its Kuga–Satake
     correspondence (no known criterion applies: ρ ≤ 3, End_Hdg(T) = Q, no non-symplectic automorphism,
     no RM/CM, not a Kummer/Shioda–Inose or six-line double plane), or by a new geometric construction of
     a divisor on S₂ × S₂ whose resolution carries so(T)(1) in its H²;
 (B) show ρ ≥ 2 and η ∈ NS_prim ⊗ C (then r lies in the T ⊗ NS sector, which the book supports) —
     testable only by certified period computations (Lairez–Sertöz) of the K3 `w² = f` and of the
     class of z0²z1²z2²; the discussion's finite-field data (one splitting line mod 3, obstructed lift mod 9)
     make ρ = 1 the likely truth, which closes this route.

4.2 **Sum principle (§1.3 ii).** served(D₁D₂) = served(D₁) + served(D₂) mod J(F). Combined with the
discussion's degree ≤ 7 exclusion it shows that every reducible candidate with all factors of degree ≤ 7
is dead, which removes the entire "cone product / structural product" part of the experiment menu and
explains why the discussion's degree-8 relaxation found `D = f_U h` as spurious necessary points: these are
products of excluded pieces. Any future degree-8/9 exclusion only needs to treat irreducible D (or D with
an irreducible factor of degree ≥ 8).

4.3 **Dichotomy (§1.3 i)** gives a cheap, decisive diagnostic for any proposed D: compute the served
target space (`target_image`), reduce to R_6 = R_f ⊗ R_f and project to the K3 block (bidegrees (6,0),(0,6)).
Nonzero projection of dimension ≥ 17 ⟺ D solves every K3-sector target at once (in particular this
instance after adding the book's supports for the other sectors via the product rule). This replaces
"consistency for one target" by a sector-level statement.

4.4 **Maximal rank of [Â] in every degree (§1.5)**: the Koszul parametrization D = Σ a_i ∂_iF is forced up
to degree 9; for d ≥ 10 the E1 kernel is the whole dimension excess (sym part 294 at d = 10, 510 at 11,
677 at 12). Sparse rational representatives of the d = 10 kernel do not exist in the standard-monomial
basis (supports ≥ 121 and no small-height rational lift), so "the E1 kernel" does not single out natural
candidates — it is just a dimension count.

4.5 **Where a construction would have to come from.** By the exterior-power theorem (book §3A.6) an
algebraic Kuga–Satake class produces a divisor of S₂ × S₂ swept by a 1-parameter family over a
Kuga–Satake surface Σ ⊂ A × A (dim A = 2^{rank T − 2} ≥ 2^{17}); the corresponding D(U,V) would have
astronomically large degree. A degree-8 or 10 denominator would require a *different* mechanism producing
a weight-2 Hodge structure of adjoint type from a low-degree hypersurface section of X. For singular
hypersurface sections the Hodge structures available are H²(S') of the base surface of a cone (vertex a
line) or H¹(Γ)⊗H¹(C') (vertex a plane) or the H² of exceptional divisors over worse singularities — all
"geometric" Hodge structures of small Mumford–Tate type unrelated to SO(T). We found no candidate
mechanism; this is consistent with the experiments.

## 5. Ranked next steps (with compute estimates)

1. **Decide ρ and test the loophole (B)** — certified periods of `w² = f` (Lairez–Sertöz `PeriodSuite`
   or Sertöz's K3 period code), LLL for NS, then check whether the period vector of η = [z0²z1²z2²] is in
   the NS_prim span. CPU: hours; needs software not present here. Outcome: either closes the instance
   conditionally (ρ = 1 ⟹ equivalent to open GHC) or exhibits the only possible cheap route.
2. **Independent replication of the degree ≤ 7 exclusion** using the Koszul parametrization in the
   (U,V) frame (my `uvring.py` gives D = Σ a_i∂_if(U) + Σ b_i ∂_if(V) directly) — degrees ≤ 6 are
   minutes, degree 7 ≈ the discussion's 261 s replay; worth doing because the whole "degree ≥ 8" frontier
   and the sum principle rest on it.
3. **Degree-8 exclusion for irreducible D** via the sum principle: the discussion's 315-parameter
   relaxation only needs to be decided off the reducible locus; a Gröbner/msolve attack restricted to the
   complement of {D = f_U h, D = (f_U − f_V)h, P_i·cubic} might finish where the unrestricted one timed out.
   CPU: hours–days (msolve), memory several GB.
4. **Sector-level diagnostics for any future candidate** (minutes each): `pfsolve.py` / `k3proj.py`
   compute the K3-block projection of the served space; by the dichotomy a single nonzero value would be
   decisive, so this is the right first test for any D anybody proposes, before trying a Lean proof.
5. **Degree-10 structural candidates** (polar pairing Σ ∂_if(U)∂_if(V) etc.): ≈ 20k × 20k dense modular
   solves, 30–60 min and 3–5 GB each on an idle machine; expected inconsistent by §1; low priority.
6. **Discussion post**: the equivalence theorem (§1.2), the dichotomy and the sum principle (§1.3), the
   maximal-rank table (§1.5), and the negative table (§2) are publishable progress; they sharpen
   1337isnot1337's "degree ≥ 8" statement to "only irreducible degree ≥ 8, and in fact only a divisor
   carrying so(T)(1)".
