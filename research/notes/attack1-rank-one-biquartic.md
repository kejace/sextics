# Attack 1 — `rank-one-biquartic` (q = 1, A = X0³Y0³)

**Status: NO WITNESS FOUND (so far).** All results below are exact (python-flint mod 2^31-range primes, rational lifts where stated).
Work directory: `/tmp/hodge/results/rank-one-biquartic/` (log of every tested denominator: `log.jsonl`; engine = `/home/kejace/devel/sextics-work/hodge`, "engine v2" after the 17:25 sign fix — every fixed-D consistency result quoted here was (re)computed with v2; E1 and E1-kernel results use my own code and are unaffected).

## 1. The obstruction and the geometric mechanism a certificate needs

### 1.1 Exact structure of (F, A)
```
F = f(X) + g(Y) + t·L(X)²·M(Y),   t = 5/7,
f = X0⁶ + (5/7)X0⁵X1 + (5/7)X0²X1³X2 + (15/7)X0X2⁵ + X1⁶ + (10/7)X1⁵X2 + X2⁶      (the split-ternary-q2 X-block)
g = Y0⁶ + (10/7)Y0⁵Y1 + (10/7)Y0³Y1Y2² + (5/7)Y0Y2⁵ + Y1⁶ + (15/7)Y1⁵Y2 + Y2⁶   (the split-ternary-q2 Y-block)
L = X0² + X0X1 + 2X1X2 + 3X2²,   M = Y0² + Y0Y1 + 2Y1Y2 + 3Y2²            [verified exactly: F − f − g − tL²M = 0]
A = X0³Y0³ (block bidegree (3,3)).
```
Jacobian ring Hilbert function (mod p) = (1+t+…+t⁴)⁶ → smooth; A ∉ J(F) (class nonzero); D = 1 fails.
Only symmetry: τ itself (projective diagonal group of order 2, no reflections) → no character sectors to exploit.

### 1.2 What the t = 0 certificate is, and why it does not deform
At t = 0 the target has block degrees (3,3), so the **elementary split Euler field is strict**:
```
D₀ = f(X),  N₀ = (A·X0/6, A·X1/6, A·X2/6, 0, 0, 0),  Q = 0
```
(check: N₀(F₀) = (A/6)·ε_X f = A f = A D₀; div N₀ = (3A + ε_X A)/6 = A; N₀(D₀) = A f; so E ≡ 0).
Geometrically V(f) ∩ X₀ = {f = 0 = g} is the **join of the two plane sextics C_f, C_g**, a ℙ¹-bundle over C_f × C_g with
H^{2,0} = H^{1,0}(C_f) ⊗ H^{1,0}(C_g) ≅ (R_f)₃ ⊗ (R_g)₃ (100-dim): exactly the bidegree-(3,3) sector of the invariant (3,1)-classes, and A = x0³ ⊗ y0³ is one of them.

For t ≠ 0 the radial fields give ε_X F = R_X := 6f + 4tL²M, ε_Y F = R_Y := 6g + 2tL²M, R_X + R_Y = 6F. The naive deformation
v = (A/R_X) ε_X has v(F) = A but div v = 8tA·L²M/R_X² ≠ 0 (mod F), so it is not a certificate. First-order deformation theory of
(D₀,N₀) forces D₁ = (2/3)L²M mod the t=0 E1-kernel, i.e. the only degree-6 continuation is D = R_X/6 — which fails (below).
Rigorous consequence of the degree-6 E1-kernel computation (§2.3): **no family of degree-6 denominators D_t with D₀ = f exists**; any certificate for t ≠ 0 has a denominator that does not specialise to f, i.e. the support must jump.

### 1.3 Quotient geometry (the triple cover)
With S = |X|², U = |Y|² the fibre equation over B = ℙ²_X × ℙ²_Y is the binary cubic a S³ + b S²U + d U³ = 0,
a = f, b = tL²M, d = g. X/τ → B is a **non-Galois triple cover** (cyclic exactly at t = 0) with branch divisor
{g = 0} ∪ Δ', Δ' = {4t³L⁶M³ + 27f²g = 0} (bidegrees (0,6)+(12,6)), Tschirnhausen bundle O(−2,−4) ⊕ O(−4,−2).
The ramification divisor of X ⇢ B is exactly R = X ∩ V(R_X) (R_X = ε_X F): a degree-6 divisor, birational to Δ'.
The map X̃ → B contracts 24 ruled surfaces (joins of the 12 points C_f ∩ {L=0} with C_g, and of C_f with the 12 points C_g ∩ {M=0}).
Natural divisors and their (τ-invariant) holomorphic 2-forms (which is what a supporting divisor needs — Gysin H^{2,0}(Z̃) → H^{3,1}(X)):
- Z_f = X ∩ V(f) = {f = 0, g + tL²M = 0}: Z_f/τ ≅_bir C_f × Σ with Σ = {w² = M·g} the double plane branched along conic ∪ sextic (p_g = 3).
  Prediction "3 classes" — **refuted by target_image: Z_f serves only J₆** (the Gysin map kills them; see §2.4).
- Z_g = X ∩ V(g): Z_g/τ ≅_bir ℙ²_X × C_g, no 2-forms (confirmed: serves only J₆).
- Base-locus divisors X ∩ V(L), X ∩ V(M) (= X₀ ∩ V(L), X₀ ∩ V(M), independent of t): quasi-smooth weighted hypersurfaces, no 2-forms (confirmed).
- R (ramification divisor): singular model Δ' with transversal A₅ along the 12 planes {p_i} × ℙ²_Y (p_i ∈ C_f ∩ V(L)); the A₅ chain is
  reversed by the monodromy of √(−M³/g), so two of its exceptional divisors are ℙ¹-bundles over Σ (3 forms each). Prediction ≤ 72 classes; see §2.4 for the computed value.
Pull-backs π⁻¹(Γ) of smooth threefolds Γ ⊂ B carry no 2-forms (Lefschetz in the ℙ¹-bundle), so a support must come from
reducible/singular pull-backs: the cubic splits over Γ iff Γ ⊂ {g p³ + tL²M p q² + f q³ = 0} with p, q bihomogeneous of bidegrees (κ+2, π), (κ, π+2) — the smallest such Γ has degree 12 (p ∈ S(2,0), q ∈ S(0,2)).

## 2. Experiments  (see §2.x below; tables are updated as runs finish)
