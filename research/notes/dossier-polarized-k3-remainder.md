# Research dossier — `polarized-k3-remainder` (OPEN)

Compiled 2026-11-01 from: `data/polarized-k3-remainder/{instance.rec,why.md,provenance.json,checks.json}`,
the complete GitHub thread https://github.com/hodge-computer/sextics/discussions/115 (saved at
`/tmp/hodge/discussions/polarized-k3-remainder.{md,json}`: exactly ONE top-level comment and ONE reply, both by
`1337isnot1337`), the full reproduction packet of that comment (local copy
`/tmp/hodge/downloads/k3-remainder-gist/packet/`, 12 reports + scripts + logs), book chapters 3A, 3, 1, 4
(`docs/03a-split-and-k3.md`, `docs/03-geometric-support.md`, `docs/01-certificates.md`,
`docs/04-known-locus-and-frontier.md`), the accepted companions `split-ternary-q2`, `crossed-ternary-asymmetric`,
`polarized-contact-coupled-d2`, `polarized-contact-deformation-d2`, `cone-ladder-linear`, the external sweep
`/tmp/hodge/ext/hodgeHoldout/sweep/` (DiamondJM, not posted in #115), the crashed-session notes
`/tmp/hodge/notes/attack1-polarized-k3-remainder.md` and `/tmp/hodge/results/polarized-k3-remainder/`, and
independent exact re-computations done while writing this file (python-flint / sympy; listed in Appendix A).

Labels used below: **[PROVEN]** = rigorous argument (human-checked, usually with exact computer algebra),
**[COMP]** = exact computation whose scope is exactly as stated, **[MOD-p]** = modular evidence only,
**[CONJ]** = conjecture/heuristic, **[FAILED]** = attempt that produced no conclusion, **[UNVERIFIED]** = claim
from an internal note not yet independently checked, **[NEW]** = derived while writing this dossier.

---

## 0. Executive summary

* Target: `q = 1`, `A = X0 X1 X2 Y0 Y1 Y2`, `F = ½(f(X+Y) + f(X−Y))` with the ternary sextic
  `f = z0^6 + z1^6 + z2^6 + (5/7)(z0^5 z1 + 2 z1^5 z2 + 3 z2^5 z0 + z0^2 z1^3 z2)`.
  The involution `τ: Y ↦ −Y` is the **block swap** `U ↔ V` (`U = X+Y`, `V = X−Y`), not a block sign, so the
  all-target split theorem (Ch. 3A.1–3A.3) does not apply.
* In the (U,V) frame the target is `Â = 2^{-6} ∏_i (U_i² − V_i²)` (8 monomials). Its two pure-block monomials
  `(U0²U1²U2² − V0²V1²V2²)/64` land in the **K3 exterior-square sector** `Λ²V₂`, `V₂ = H²(S₂)_prim`,
  `S₂ = {w² = f} ⊂ P(1,1,1,3)` (rank 210, Hodge numbers (19,172,19)); its six mixed monomials (U/V bidegrees
  (4,2),(2,4)) land in the **primitive sixth-root pair** sector (rank 441), which the book supports only
  *existentially* via the weighted companion `W_{23,f} = {f + s² + t³} ⊂ P(1,1,1,3,2)` (no explicit or Lean
  construction of that companion exists in the repository).
* **[PROVEN, thread]** Every legal denominator has degree ≥ 8 (all degrees 0–7 excluded over ℂ, non-strict
  certificates and complex coefficients allowed). Degree 8 is undecided; the necessary relaxation used for
  degrees 6–7 is *nonempty* at degree 8 (spurious reducible points `D = f(U)·h`), so that method cannot finish.
* **[PROVEN, thread]** For the K3 `S₂`: smooth reduction at 3; Frobenius polynomial `(z−1)^4 Q₁₈(z)/3`;
  `1 ≤ ρ(S₂,ℂ) ≤ 3`; `End_Hdg(T) = ℚ` (Galois group of `Q₁₈` is `(C₂)^9 ⋊ S₉`); no non-symplectic order-3
  automorphism, no rational Hodge similarity to any K3 with RM/CM or with a finite-order transcendental
  isometry ≠ ±1; no splitting lines/conics, no class `E² = 0, H·E = 3`, no sum-of-three-cubic-squares
  presentation, Waring rank ≥ 10; the fixed sextic is not linearly conjugate to the order-three normal form.
* **[UNVERIFIED but carefully argued, crashed-session note]** Under `End_Hdg(T) = ℚ` and `ρ = 1` (or more
  generally `η = [z0²z1²z2²] ∉ NS_prim ⊗ ℂ`), `HasCertificate` for this instance **⇔** the generalized Hodge
  conjecture (coniveau 1) for `Λ²T(S₂) ⊂ H⁴(S₂ × S₂)`, i.e. algebraicity-type input (e.g. an algebraic
  Kuga–Satake correspondence) for this specific low-Picard K3. No known theorem supplies it.
* Consequently: an explicit low-degree witness is **not expected to exist at all** (any supporting divisor must
  carry a weight-2 Hodge structure of adjoint type `so(T)(1)`, `h^{2,0} ≥ 17`), and a nonconstructive Lean proof
  would need new mathematics. Realistic contributions are: sharpening the obstruction (degree 8/9 for
  irreducible D), deciding `ρ` and the loophole, a constructive certificate for the *mixed* (sixth-root) part of
  the target via the `(2,3)`-weighted companion, and publishing the equivalence theorem.

---

## 1. Exact target

### 1.1 Statement (verbatim data, `instance.rec`)

```
q = 1,   A = X0*X1*X2*Y0*Y1*Y2
F = X0^6 + 5/7 X0^5 X1 + 15 X0^4 Y0^2 + 25/7 X0^4 Y0 Y1 + 50/7 X0^3 X1 Y0^2 + 5/7 X0^2 X1^3 X2 + 15/7 X0^2 X1^2 Y1 Y2
  + 15/7 X0^2 X1 X2 Y1^2 + 15 X0^2 Y0^4 + 50/7 X0^2 Y0^3 Y1 + 5/7 X0^2 Y1^3 Y2 + 10/7 X0 X1^3 Y0 Y2 + 30/7 X0 X1^2 X2 Y0 Y1
  + 25/7 X0 X1 Y0^4 + 30/7 X0 X1 Y0 Y1^2 Y2 + 15/7 X0 X2^5 + 150/7 X0 X2^3 Y2^2 + 10/7 X0 X2 Y0 Y1^3 + 75/7 X0 X2 Y2^4
  + X1^6 + 10/7 X1^5 X2 + 15 X1^4 Y1^2 + 50/7 X1^4 Y1 Y2 + 5/7 X1^3 X2 Y0^2 + 100/7 X1^3 X2 Y1^2 + 15/7 X1^2 Y0^2 Y1 Y2
  + 15 X1^2 Y1^4 + 100/7 X1^2 Y1^3 Y2 + 15/7 X1 X2 Y0^2 Y1^2 + 50/7 X1 X2 Y1^4 + X2^6 + 75/7 X2^4 Y0 Y2 + 15 X2^4 Y2^2
  + 150/7 X2^2 Y0 Y2^3 + 15 X2^2 Y2^4 + Y0^6 + 5/7 Y0^5 Y1 + 5/7 Y0^2 Y1^3 Y2 + 15/7 Y0 Y2^5 + Y1^6 + 10/7 Y1^5 Y2 + Y2^6
```

42 monomials; Y-degree distribution 7 / 14 / 14 / 7 in Y-degrees 0 / 2 / 4 / 6 (so `F(X,−Y) = F`). `A` has
degree 6, is Y-odd (`A(X,−Y) = −A`), and `deg A = 6q` as required. Variables ordered `X0,X1,X2,Y0,Y1,Y2`.

### 1.2 Sparse presentation **[COMP, re-verified exactly]**

With `U = X + Y`, `V = X − Y`:

```
F = ½ ( f(U) + f(V) ),
f(z0,z1,z2) = z0^6 + z1^6 + z2^6 + (5/7)·( z0^5 z1 + 2 z1^5 z2 + 3 z2^5 z0 + z0^2 z1^3 z2 )
            = z0^6 + 5/7 z0^5 z1 + 5/7 z0^2 z1^3 z2 + 15/7 z0 z2^5 + z1^6 + 10/7 z1^5 z2 + z2^6     (7 monomials)
Â(U,V) := A((U+V)/2, (U−V)/2) = 2^{-6} ∏_{i=0}^{2} (U_i² − V_i²)
        = 2^{-6} [ U0²U1²U2² − (U0²U1²V2² + U0²U2²V1² + U1²U2²V0²) + (U0²V1²V2² + U1²V0²V2² + U2²V0²V1²) − V0²V1²V2² ]
```

(`sympy`: `expand(F − ½(f(U)+f(V))) == 0`.) The book writes `F_f = f(U) + f(V)`; the record's `½` is a harmless
scalar. The same ternary `f` is the X-block of the accepted `split-ternary-q2` (there `F = f(X) + g(Y)` with a
different `g`), and the related 5/7-coefficient family appears in `crossed-ternary-asymmetric`.

### 1.3 Degrees of a certificate

`deg D = d` (legal: `D ≠ 0`, `gcd(D,F) = 1` over ℂ), `deg N_i = d + 6q − 5 = d + 1`, `deg E = 2d + 6`,
`deg F^{q+1} = 12`, so `deg Q = 2d − 6` (strict certificates forced for `d ≤ 2`). Parity: with `D` Y-even one may
take `N_X` Y-odd, `N_Y` Y-even (field `v = N/D` Y-odd), `Q` Y-odd; averaging over `τ` makes this harmless
(`formal/README.md`). Unknown counts for a fixed `D` (parity-restricted):

| d | N unknowns (6·#deg d+1 /2) | Q unknowns | equations (deg 2d+6 monomials /2) |
|---|---|---|---|
| 6 | 2376 | 226 | ≈ 16 800 |
| 7 | 3861 | 643 | ≈ 26 565 |
| 8 | 6006 | 1491 | ≈ 40 365 |
| 9 | 9009 | 3003 | ≈ 59 000 |
| 10 | 13104 | 5814 | ≈ 85 000 |
| 12 | 25704 | 16824 | ≈ 162 000 |

(Engine run at d = 8 built exactly `rows=40326 cols=7497 (N 6006, Q 1491)`.)

### 1.4 Symmetry **[COMP, re-verified]**

* **Projective diagonal group of F:** order 2, generated by `τ = diag(1,1,1,−1,−1,−1)`. (Smith normal form of
  the exponent-difference lattice of F: invariant factors `(1,1,1,1,2,0)` → torus part = scalars only, finite
  part `ℤ/2`.) Matches `checks.json`: `projectiveDiagonalGroupOrder: 2`, `coordinateReflections: []`,
  `juniorElementsFixingTarget: []`.
* **Diagonal symmetries of f:** none (Smith form `(1,1,0)`: scalars only). No coordinate permutation preserves
  `f` (coefficients 5/7, 10/7, 15/7 on distinct monomial shapes).
* **Residue character of A:** `A` is τ-odd; the volume form `Ω = Σ ± Z_i dZ_0∧…` is τ-odd (three negated
  coordinates), so `ω_A = AΩ/F²` is **τ-invariant**: the residue lies in `H⁴(X,ℂ)^τ`, which has no (4,0) part
  (the (4,0) class is τ-anti-invariant), so the residue class of a q = 1 target is purely of type (3,1). The
  τ-invariant primitive Hodge numbers are (208, 866, 208) (Ch. 3A.5).
* **Non-diagonal structure noted in the sources:** (i) the block swap `U ↔ V` (this *is* τ in the (U,V) frame);
  (ii) the μ₆-grading `h(U,V) = (U, ηV)` on the auxiliary cover `S₆ × S₆` (not a symmetry of F) which grades the
  (U,V) monomials by V-degree `b` with residue character `η^{b+3}` (Ch. 3A.5); (iii) **no** other linear
  decomposition exists: the Hessian symmetrizer space has dimension 2 (`checks.json`,
  `hessianSymmetrizerDimension: 2`; why.md: "the U/V splitting is the only nontrivial linear decomposition").
  **[NEW, conditional]** If `ρ(S₂) = 1` then `Aut(S₂) = {1, ι}` (covering involution), every automorphism of the
  plane sextic `C = {f=0}` lifts to `S₂`, hence `Aut(C) = 1` and the full projective linear automorphism group
  of F is `{1, τ}` (any automorphism must preserve the unique splitting up to swap). Not stated in the thread.
* **Character decomposition of the target** (Ch. 3A.5 "polynomial test"): monomial of V-degree `b` ↦ character
  `η^{b+3}`; `b ≡ 0 (6)` → K3 sector `Λ²V₂`; `b ≡ 2, 4` → primitive sixth-root pair `H₁⊗H₅ ⊕ H₅⊗H₁`;
  `b ≡ 1, 5` → primitive third-root pair `H₂⊗H₄ ⊕ H₄⊗H₂`; `b ≡ 3` → joining `Λ²H¹(C)(−1)` and Tate.
  For `Â`: `b ∈ {0, 6}` (two pure monomials) → **K3**; `b ∈ {2, 4}` (six mixed monomials) → **sixth-root pair**.
  No third-root or joining component. (Book: `P₀Â = (U0²U1²U2² − V0²V1²V2²)/64 ≠ 0`, equivalently
  `A(X,X) = X0²X1²X2² ≠ 0`; `checks.json` `diagonalTarget: X0^2*X1^2*X2^2`.)
* **Sector dimensions at q = 1 [COMP, re-derived]:** swap-antisymmetric sextics in (U,V): `(462 − 10)/2 = 226`
  (= Y-odd sextics in (X,Y)); of these 28 are pure (`B(U) − B(V)`), 198 mixed. `J(F)_6` has dimension 36
  (no linear syzygies), antisymmetric part 18 (9 pure, 9 mixed). Hence `h^{3,1}(X)^τ = 226 − 18 = 208 =
  19 (K3: 28 − 9) + 189 (mixed: 198 − 9)`, matching Ch. 3A.5.

### 1.5 Smoothness **[PROVEN]**

All non-pure coefficients of `F` (15, 5/7, 25/7, 50/7, 10/7, 30/7, 150/7, 75/7, 100/7) are divisible by 5 in
`ℤ_(5)` (7 is a unit), so `F ≡ Σ_i X_i^6 + Σ_i Y_i^6 (mod 5)` — the Fermat sextic — and the reduced partials are
`X_i^5, Y_i^5` (`checks.json`: `smoothReductionPrime: 5`, `reducedPartials`). The Fermat fourfold is smooth,
so the discriminant of F is a 5-adic unit, hence nonzero: X is smooth over ℚ and ℂ. Equivalently (packet
`REPORT.md`): the six partials form a homogeneous regular sequence; `Hilb(S/J(F)) = (1+t+…+t⁴)^6`
(1, 6, 21, 56, 126, 246, 426, 666, 951, 1246, 1506, 1686, 1751, …, socle 24). Likewise `7f ≡ 2(z0^6+z1^6+z2^6)
(mod 5)`, i.e. `f ≡ Fermat (mod 5)` [COMP], so `C = {f = 0}` is smooth and `S₂` is a smooth K3.
`Hilb(R_f) = (1,3,6,10,15,18,19,18,15,10,6,3,1)` [COMP, exact over ℚ]; `η = [z0²z1²z2²] ≠ 0` in `R_{f,6}`
(rank of `J(f)_6` is 9, rank with η appended is 10) [COMP, exact over ℚ] — this is the why.md claim "the
exponents two in B show its projection survives", here checked directly in characteristic 0.

### 1.6 Where the target lives (Ch. 3A.4–3A.5, **[PROVEN in the book]**)

Via `Φ: S₆ × S₆ ⇢ {f(U)+f(V)=0}`, `S₆ = {u⁶ = f(z)} ⊂ P³`, `((z,u),(z',v)) ↦ [v·z : α·u·z']`, blown up along
`C × C` and quotiented by diagonal μ₆ (Shioda–Katsura split blow-up), one gets algebraic correspondences and
projectors with

```
H⁴(X)_prim^τ = [sixth-root pair, 441, (90,261,90)]  ⊕  [third-root pair, 441, (54,333,54)]
             ⊕  [Λ²H¹(C)(−1), 190, (45,100,45)]      ⊕  [Λ²V₂, 210, (19,172,19)]
```

The first three are supported by: second elliptic companion `W_{23,f} = {f + s² + t³} ⊂ P(1,1,1,3,2)` with
`E₂₃ = {a² + b³ + c⁶} ⊂ P(3,2,1)` (sixth-root pair); triple-plane companion `W_f = {f + s³ + t³} ⊂ P(1,1,1,2,2)`
with the Fermat cubic (third-root pair); the joining divisor swept by lines between the two copies of `C`
(joining sector). "Support" is Chapter-3 Fano covering-curve existence, i.e. finite algebraic operations in
principle, **not** an explicit certificate. The K3 summand `Λ²V₂` is "Unresolved in general"; it is supported
iff (sufficient criterion) `S₂` has an algebraic Kuga–Satake correspondence (exterior-power theorem, Ch. 3A.6),
known for the order-three family `f = g₆(z0,z1) + z2³g₃(z0,z1) + a z2⁶` (Bolognesi–Laterveer) — which this `f`
is not (§2, §3).

In Jacobian-ring terms (`R_F = R_f ⊗ R_f` in the (U,V) frame), the K3 component of the residue is
`r_K3 = η ⊗ ω − ω ⊗ η ∈ (Λ²V₂)^{3,1} = ω ∧ V₂^{1,1}` with `ω ↔ 1 ∈ R_{f,0}` the K3 two-form and
`η = [z0²z1²z2²] ∈ R_{f,6} ≅ H^{1,1}_prim(S₂)`.

---

## 2. Companion / related solved instances and why their methods fail here

| instance (status) | F shape | method | denominator | why it fails for `polarized-k3-remainder` |
|---|---|---|---|---|
| `split-ternary-q2` (accepted, `fqlx`, PR #159, q = 2) | `f(X) + g(Y)` with **the same f**, `A = X0⁴X1³Y0³Y1²` | block-sign split theorem made constructive: polynomial companion curves on `7(f+S³+T³)`, second normal residue `R₂`, radial descent, Fermat-elliptic contraction, normalized trace, physical transport, homogeneous clearing (Lean, existential; no expanded D) | unspecified finite degree, `Q = 0` | the involution here swaps the blocks instead of negating one; the all-target split theorem is for `f(X)+g(Y)` only (Ch. 3A.1). The companion machinery (μ₃, `s³+t³`) serves the *third-root* sector; our target has no third-root component, and the K3 sector is outside every companion. |
| `crossed-ternary-asymmetric` (accepted, `1337isnot1337`, PR #149, q = 1) | `P(X0,X1,Y0) + P(Y1,Y2,X2)`, `A = X0⁴X1Y1` | involution-stable ternary decomposition (Ch. 3A.3) + explicit moving-conic companion on `P + S³ + T³`, elliptic kernels with infinity correction, translated-diagonal contraction, Lean existential | finite degree, non-strict (`Q ≠ 0` allowed) | same obstruction: requires a block-sign (ρ-invariant) sector; nothing in the K3 sector. Its reusable library `EllipticCompanions` is relevant for the *mixed* part of our target only after replacing the cubic companion by the `(2,3)` weighted one. |
| `polarized-contact-coupled-d2`, `polarized-contact-deformation-d2` (accepted, d = 2 strict) | `f(U)+f(V)+1000003·r(U)r(V)` with `f = r² + z2·Q` (square contact along the line `z2 = 0`), `A = g(U) − g(V)` | Hamiltonian contact field `W = (−r_{z1}, r_{z0}, 0)` with two first integrals; `D = U2·V2 = X2² − Y2²`, `N = (V2 W(U), −U2 W(V))`, `div N = 0`, `N(D) = 0` | **degree 2**, `D = X2² − Y2²` | works because the K3 `w² = f` has a **splitting line** (`z2 = 0` lifts to two components `w = ±r`), i.e. extra Néron–Severi classes and a pure-block target built from the contact data. For our `f`, splitting lines and smooth splitting conics are **excluded over ℂ** (packet `SMALL_GEOMETRY.md`), and `1 ≤ ρ ≤ 3` with no low-degree special curves; no contact presentation `f = r² + ℓ·Q` exists. |
| `cone-ladder-linear`, `contact46006` (accepted, d = 1 / cubic contact) | `f(X)+f(Y)` with `f = c·e³ + P·k` | contact ladder fields `∇P × ∇e` times passive `Y0` | degree 1 (`X1 − X2`) | block-sign presentation plus a perfect-power contact; neither available here. |
| order-three polarized family (book Ch. 3A.6, no catalogue instance) | `f = g₆ + z2³g₃ + a z2⁶` | Bolognesi–Laterveer algebraic Kuga–Satake ⇒ exterior-power theorem ⇒ existence for every target and q | astronomically large (via Kuga–Satake surfaces), never explicit | **[PROVEN, packet]** our `f` is not linearly conjugate to that normal form (Hessian has no rank-≤1 point: six 2×2 minors generate all degree-14 forms, rank 120 mod 101 → ℚ), and more strongly `S₂` has **no** non-symplectic order-3 automorphism and no rational Hodge similarity to such a K3. |
| Fermat fourfold (`fermat-*`, Fermat joins) | `Σ Z_i^6` | Fermat join `D = X0^6+X1^6+X2^6` etc. | degree 6 | `F ≡ Fermat (mod 5)` only; over ℂ the Fermat double plane has `ρ = 20`, `T` of rank 2 (CM), so its K3 (3,1) sector is `ω ∧ NS_prim` (supported by curves × K3). Deforming to our `f` (`ρ ≤ 3`) destroys exactly that support; mod-5 solvability would not lift. |

Why.md's own statement: "A linear change cannot turn the block swap into a block-sign presentation"; the
remainder `Λ²H²(w²=f)_prim` "reaches the K3 summand" because `B(U) − B(V) ≠ 0` modulo J; "no such
[Kuga–Satake] correspondence for this particular ternary sextic is supplied".

---

## 3. Prior work in discussion #115, comment by comment

The thread contains only `1337isnot1337`'s material (prepared "with OpenAI Codex assistance"), all in one
top-level comment (node `DC_kwDOVLOOJc4BIp5k`, created 2026-10-19T16:57:05Z, last edited 2026-10-19T23:44:08Z,
URL https://github.com/hodge-computer/sextics/discussions/115#discussioncomment-19045988) and one reply
(node `DC_kwDOVLOOJc4BIr7M`, 2026-10-20T00:00:41Z; decoded URL
https://github.com/hodge-computer/sextics/discussions/115#discussioncomment-19054284 — decoded from the node
id, unverified). Reproduction packet: https://gist.github.com/1337isnot1337/d504e828c4f50d8dade0687792ea9ec0
(local copy `/tmp/hodge/downloads/k3-remainder-gist/packet/`). Upstream pinned at
`b56894e8…` → `bfb13916…` → `4faa55ef…`; target unchanged. Instance was chosen by `secrets.choice` from the
29 then-open ids (`sextics-choice.json`); no reroll.

The comment was edited through the day; the packet's dated reports let one reconstruct the chronology.

### 3.1 `REPORT.md` — degrees 0–5 excluded **[PROVEN]**

Normalization: from `E = F²Q`, reduce mod F; `gcd(D,F) = 1` and F irreducible give `AD − N(F) = FH`,
`deg H = d`; Euler shift `N' = N + Hε/6` gives
```
N'(F) = A·D,      B_D(N') := D·div N' − N'(D) ∈ (F).
```
Necessary for non-strict certificates too; in particular `A·D ∈ J(F)`.

* Degrees 0–4: original-coordinate matrices mod 101: `Σ F_i N_i` (deg N = 5): 3003 × 1512, rank 1497 (= char-0
  rank via 15 Koszul relations); with all `A·D`, deg D = 4: 3003 × 1638, rank 1623 → the 126 columns `A·D` are
  independent mod J ⇒ no quartic `D` with `AD ∈ J`; lower degrees by multiplying with monomials.
* Degree 5: ranks 2682 / 2928 (2772 Jacobian + 252 denominator columns, 4368 rows, mod 101); increment 246 =
  `dim (S/J)_5` ⇒ multiplication by A on `(S/J)_5` injective in char 0 ⇒ `D = Σ c_i F_i`.
  Then with `G = 7F`, `P_i = ∂_iG`, `H_ij = ∂_i∂_jG`, every degree-6 numerator is `M = A·c + Σ s_{jkl} z_l (P_k e_j − P_j e_k)`
  (90 syzygy parameters, Koszul exactness). Relaxing products `c_i s_jkl` to independent unknowns, let `W ⊂ S_10`
  be spanned by 540 polynomials `B_{P_i}(z_l(P_k e_j − P_j e_k))` and 126 `G z^m`. The quadratic defect is
  `B_D(M₀) = Σ c_i² T_ii + Σ_{i<j} c_i c_j T_ij`, `T_ii = P_i ∂_iA − A H_ii`, `T_ij = P_i∂_jA + P_j∂_iA − 2A H_ij`.
  Exact ℚ ranks: `dim W = 411`, `dim(W + ⟨T_ij⟩) = 432` ⇒ all 21 quadratic monomials `c_i c_j` forced to 0 ⇒ `c = 0`.
  Same ranks mod 101, 103. Separately, 12 fixed quintics (`G_{Z_i}`, `G_{X_i} ± G_{Y_i}`) had been tested
  (ranks 210/211, 195/196).
* Side facts: multiplication by `A` on the Jacobian quotient is injective in grades 0–9 mod 101 (grade 9:
  rank 1246 mod 101 but 1244 mod 103 — "a rank drop at one prime is not a char-0 kernel"); at grade 10 the map
  `R_10 (1506) → R_16 (951)` is surjective, so E1 no longer forces `D ∈ J` there.
  Files: `sextics-direct-degree-four.py`, `sextics-jacobian-screen.py`, `sextics-all-polar-screen.py`,
  `sextics-polar-screen.py`, `sextics-jacobian-101/103.jsonl`.

### 3.2 `DEGREE_SIX.md` — degree 6 excluded **[PROVEN]**, plus Hessian screen

* Grade 6: `A·D ∈ J ⇒ D ∈ J` (dimension and rank both 426) ⇒ `D = Σ P_i a_i`, `a_i = Σ_j C_ij z_j` (unique C,
  no linear syzygies). Split coordinates, `G = 14F̂`, `T = 64Â = ∏(U_i² − V_i²)`. Every normalized numerator is
  `M = T a + S`, `S = Σ_{j<k} q_jk (P_k e_j − P_j e_k)`, `q_jk` quadratic.
* Enlarged correction ideal `L = (G) + J(G)² + (R_ijk)`, `R_ijk = P_k H_ij − P_j H_ik`; `B_D(S) ∈ L` always, so a
  certificate needs `B_D(Ta) ∈ L_12`. Degree-12 matrix: 3927 columns (54 primitive `R` × 56 cubics, 21 `P_iP_j`
  × 21 quadrics, `G` × 462 sextics), split by V-degree mod 6 (rows 966/990/1071/1100/1071/990; correction ranks
  603/558/396/260/396/558; augmented 713/715/617/470/617/715). 666 quadratic monomials in the 36 entries of C.
* Euler gauge `C_55 = 0`. Exact ℚ row reduction: each of the 30 off-diagonal squares `C_ij²` lies in the
  constraint span ⇒ off-diagonal entries vanish over ℂ; then rank 15 on the 15 quadratic monomials of the 5
  remaining diagonal entries ⇒ all vanish ⇒ `C = λI`, `D = 6λG` — zero or ∝ F, illegal. (≈ 89 s.)
* Fixed-denominator script also excluded 7 structured sextics; illegal control `D = G` passes.
* **Hessian screen [PROVEN]:** `Hess(f)` has no rank-≤1 point over ℂ (the six primitive 2×2 minors' sextic
  multiples span all degree-14 forms: 120 × 168, rank 120 mod 101 ⇒ ℚ). A smooth order-three normal form has
  Hessian rank 1 at `[0:0:1]`. Hence `f` is not linearly conjugate to the order-three family.
* Notes the withdrawn OpenAI preprint "rational Hodge conjecture for products of K3 surfaces" (Oct 4, 2026;
  withdrawn Oct 6 for a sign error) — not used.

### 3.3 `DEGREE_SEVEN.md` — degree 7 excluded **[PROVEN]**

* Grade 7: injective (666 = 666) ⇒ `D = Σ P_i a_i` with quadratic `a_i` (126 coefficients), unique.
  Euler gauge `D' = D + Gh`, `M' = M + Thz/6` (`h` linear) ⇒ zero the 6 coefficients of `a_5` divisible by `z_5`
  → 120 parameters `c_0..c_119`, not all zero for legal D.
* Necessary quadrics: `B_D(Ta) ∈ L_14`; 17 541 correction columns, 7 260 quadratic monomials; by V-degree
  residue: rows 1926/1926/1926/1935/1980/1935; correction rank over ℚ = over F₁₀₁ = 1635/1602/1635/1464/1287/1464;
  augmented ranks over F₁₀₁ 1926/1920/1926/1934/1975/1934; membership equations 2 528 rows, rank 1 779.
* Projective emptiness over `F̄₁₀₁`: a global Gröbner basis in 120 variables exceeded a 2800 MiB guard
  (**[FAILED]**, proves nothing); instead charts `c_0=…=c_{k−1}=0, c_k=1` for `k = 0..8` all give the unit ideal
  (Singular); then squares in the row span force 69 more zeros, leaving 42 coordinates whose homogeneous ideal
  has affine dimension 0 (basis 1 194, quotient dimension 309) ⇒ empty projective set.
* Transfer to ℂ: equal ranks over ℚ and F₁₀₁ ⇒ unit maximal minors over `ℤ_(101)`, free cokernel, so the
  reduced equations are the reduction of the membership scheme ⊂ `P^119_{ℤ_(101)}`; empty special fiber +
  properness ⇒ empty generic fiber ⇒ no complex point. Independent replay passed in 261 s
  (`sextics-degree-seven-projective-exclusion.py`, logs `degree-seven-replay-*.log`).

### 3.4 `DEGREE_EIGHT.md` — assembled, **undecided**; relaxation provably nonempty

* Grade 8 injective (951 = 951) ⇒ `D = Σ P_i a_i`, cubic `a_i`; Euler gauge removes 21 → **315 parameters**.
  `L` has Gröbner bases of 563 elements over ℚ and F₁₀₁ with identical Hilbert function
  `1,6,21,56,126,252,461,786,1266,1894,2580,3162,3417,3192,2541,1714,1050,666,426,246,126,56,21,6,1,0`;
  `dim (S/L)_16 = 1050`. All 49 770 quadratic monomials assembled; modular row rank 1050; 6 277 983 nonzero
  coefficients; 46 sampled columns verified against direct differentiation. No square appears in the row span.
* **Explicit spurious points [PROVEN]:** `f_U := 7f(U)`, `h` any quadratic eigenvector of the block sign
  `V ↦ −V` (12-dim even / 9-dim odd), `D = f_U h`, `a_i = U_i h/6 (i<3)`, `a_i = 0 (i≥3)`. Then
  `B_D(Ta) = (f_U h²/6)(E_U T − 3T)`, `E_U = Σ_{i<3} U_i ∂_{U_i}`; `f_U(E_UT − 3T) ≠ 0` in grade 12 (126-term
  normal form) but its products with all 66 quartic monomials of even V-degree vanish in grade 16 (22 nonzero /
  38 zero among the other 60). So the relaxation contains positive-dimensional rational families and **cannot**
  exclude degree 8. (`sextics-degree-eight-block-euler.py`.)
* **One fixed octic excluded [PROVEN]:** `D = f_U U0²`: full normalized search (all quartic Koszul
  coefficients, all degree-10 multiples of G) is inconsistent — V-degree-0 block coefficient rank 756,
  augmented 757, over ℚ and F₁₀₁ (3366 × 888). Illegal control `D = G U0²` passes (`B_D(Ta) = G T U0⁴`; block
  ranks 486/486, 540/540, 486/486). No other octic excluded.
* **[FAILED]** affine chart `c_0 = 1` with Singular `slimgb` (memory guard 2800 MiB, 26 s); msolve 0.6.5 (7.9 GB
  exponent array; 2800 MiB and 10 GiB runs failed at import / hash expansion); patched msolve (storage-only
  patch, byte-identical on a degree-7 control) ran degree-3 rounds then hit a 600 s timeout; Singular `satstd`
  saturation by `c_0` exceeded 2800 MiB (29 s). No conclusion.
* **Limit:** in degrees ≥ 17 the relaxed ideal has the same graded quotients as `J(G)`, so higher-degree versions
  of this relaxation have many solutions regardless of certificates — "continuing this relaxation alone cannot
  decide unrestricted HasCertificate".

### 3.5 `PICARD_BOUND.md` — `1 ≤ ρ(S₂,ℂ) ≤ 4` **[PROVEN]** (van Luijk method)

Integer twist `w² = 7f`; mod 3 the branch is `a⁶+b⁶+c⁶+2a⁵b+b⁵c+2a²b³c`, smooth (unit ideals on all charts over F₃).
Exact point counts over `F_{3^n}`, `n = 1..11` (C++ enumeration, independent tuple check for n ≤ 3):
`#S(F_q) = 16, 88, 865, 6796, 59941, 533197, 4794148, 43056604, 387466417, 3486963493, 31381582303`;
traces on H²: `6, 6, 135, 234, 891, 1755, 11178, 9882, 45927, 179091, 522693`. Newton identities + functional
equation (middle coefficient −354294 ≠ 0 fixes the sign) give
```
P(3z)/3^22 = (z−1)^4 · Q₁₈(z)/3,
Q₁₈ = 3z^18+6z^17+11z^16+13z^15+15z^14+14z^13+14z^12+12z^11+12z^10+11z^9+12z^8+12z^7+14z^6+14z^5+15z^4+13z^3+11z^2+6z+3,
```
irreducible, non-cyclotomic (monic normalization non-integral) ⇒ exactly 4 root-of-unity eigenvalues ⇒
`ρ(S₂,ℂ) ≤ 4`. Also notes the Sept-30-2026 OpenAI "K3 quadratic locus" preprint needs an isometric embedding of
`T` into an 8-dimensional space — impossible here (`dim T ≥ 18`).

### 3.6 `PICARD_LIFT_OBSTRUCTION.md` — `ρ ≤ 3` **[PROVEN]** (Elsenhans–Jahnel)

On the line `a + b = 0`: `f₇(−b,b,c) = 9b⁶+7c⁶+15b⁵c−15bc⁵ ≡ c⁶ (mod 3)`, so the special fiber contains the
(−2)-curve `R: a+b=0, w=c³` (`H·R = 1`). First-order lifts `a = −(1+3λ)b − 3μc`, `w = c³ + 3(αb³+βb²c+γbc²+δc³)`
(parameters in `F̄₃`): coefficients of `b⁶, b⁵c, b⁴c²` give `−λ = 0`, `2 − 2λ − μ = 0`, `−2μ = 0` —
inconsistent (rank 6 vs augmented 7). Riemann–Roch (`h⁰(O(R)) = 1`) + Nakayama ⇒ the line bundle `O(R)` does
not lift mod 9. If `ρ(ℂ) = 4`, specialization would be an isomorphism with torsion-free cokernel
(Elsenhans–Jahnel Thm 1.4, `e = 1 < p − 1`), so `R` would lift — contradiction. Hence `ρ ≤ 3`.

### 3.7 `ORDER_THREE_EXCLUSION.md` — no non-symplectic order-3 automorphism **[PROVEN]**

If `ρ(ℂ) = 4`: Artin–Tate (valid under the temporary equal-rank hypothesis) gives
`|Br|·|disc NS(S_{F₃})| = Q₁₈(1) = 211` (prime), so `disc NS ≡ −211 (mod squares)`. Artebani–Sarti lattice
theory: an order-3 non-symplectic automorphism forces even `ρ` (so 2 or 4), invariant lattice `U` or `U(3)`;
`ρ = 2`: `NS ∈ {U, U(3)}` both incompatible with an ample `H`, `H² = 2`; `ρ = 4` trivial action: 3-elementary
disc square class `−1` or `−3` ≠ `−211`; nontrivial action: `disc ≡ −3`. All excluded.

### 3.8 `FROBENIUS_ORDER_THREE.md` — finite-order Hodge isometries of `T` are ±1 **[PROVEN]**

`Q₁₈(1) = 211`, `Q₁₈(−1) = 9`; on the 18-dim non-Tate summand `W` the alternating matrix `A(g − g^{-1})` has
`det(g − g^{-1}) = 211`, so `disc W ≡ 211` (a square in ℚ₅). A fixed-point-free order-3 isometry would force
`disc ≡ 3^9 ≡ 3` in ℚ₅ (non-square) — contradiction; orders 7 and 19 (`18 = (p−1)m`) via p-adic cohomology give
`disc ≡ p^m`, parity mismatch; order 4 via ℚ₂ (211 ≡ 3 mod 8 non-square). Transfer from a rational Hodge
isometry uses Buskin's theorem (every rational Hodge isometry of K3s is algebraic) + specialization. `T` is a
simple Hodge structure. Conclusion: no rational Hodge similarity to a K3 with a non-symplectic order-3
automorphism, or with any finite-order transcendental isometry ≠ ±1.

### 3.9 `FROBENIUS_GALOIS.md` + `HODGE_ENDOMORPHISMS.md` — `End_Hdg(T) = ℚ` **[PROVEN modulo cited theorems]**

`Q₁₈(z) = z⁹P(z + 1/z)`, `P(y) = 3y⁹+6y⁸−16y⁷−35y⁶+19y⁵+56y⁴+3y³−23y²−5y+1`; `Q₁₈ mod 2` irreducible of degree 18;
`P mod 23` square-free with factor degrees (4,5) ⇒ element (4)(5) ⇒ primitive ⇒ (Jordan, via Jones) `A₉ ⊂ G`,
odd ⇒ `S₉`; `P mod 61` factors with pair-cycle lengths (1,3,5) and flip parities (even, even, odd) ⇒ 15th power
flips exactly 5 pairs ⇒ full kernel ⇒ `Gal = (C₂)^9 ⋊ S₉`, order 185 794 560; 9 real roots of `P` in (−2,2).
No quotient of two distinct roots is a root of unity ⇒ centralizer of every `g^f` on `W_ℓ` is `F ⊗ ℚ_ℓ`
(`F = ℚ(λ)`, degree 18). `F` has no quadratic subfield; any homomorphism `G → S_d`, `d ≤ 7`, factors through
`C₂ × C₂`. Then (Zarhin: `E = End_Hdg(T)` is totally real or CM; Charles/DMOS: `E ⊗ ℚ_ℓ` commutes with a
finite-index Galois subgroup) a case analysis over `(ρ, n = 22 − ρ, e = [E:ℚ], m = n/e)` — totally real needs
`m ≥ 3`; `m ≥ 2` cases (`e ∈ {3,7}` for ρ=1; `{2,4,5,10}` for ρ=2; none for ρ=3) killed by Chebotarev + the
small-quotient lemma + no quadratic subfield; `m = 1` (ρ = 2, CM, e = 20) killed by a permutation-character
argument (`χ_E ≥ χ_F` with equal averages but `20 ≠ 18`) — gives `End_Hdg(T) = ℚ`. Consequences stated: every
rational Hodge isometry of `T` is ±1; no RM/CM; the *ordinary* rational Hodge conjecture for `S₂ × S₂` holds
(divisor products, Tate summands, transcendental projector) but this is **not** the coniveau statement needed.

### 3.10 `SMALL_GEOMETRY.md` (the reply of 2026-10-20) **[PROVEN]**

On the reduction at 3 all geometric NS classes are Frobenius-fixed (no cyclotomic factor besides `(z−1)^4`), so
every (−2)-curve and every line bundle descends to F₃. All 13 F₃-lines: exactly one splitting line `a + b = 0`
(`f₃|line = c⁶`), components `R_±` with `R_+ + R_- = H`, `R_±² = −2`, `R_+R_- = 3`. All 234 smooth F₃-conics:
none splits. All 29 525 cubics `g` mod sign: `f₃ − g²` never has an irreducible cubic factor ⇒ no class
`E² = 0, H·E = 3` on the reduction (uses RR, adjunction, the reducible cases `iR_+ + (3−i)R_-` with squares
`−18+30i−10i²`), hence none over ℂ (specialization). Independently at 101: the splitting-line incidence
(7 equations, 6 unknowns, three dual charts) has unit Gröbner bases ⇒ no complex splitting line (properness).
Corollary: `f ≠ g₀² + g₁² + g₂²` for cubics (would give rulings `E, E'` with `E + E' = 3H`, `H·E = 3`).
Waring rank ≥ 10: middle catalecticant has rank 9 with apolar cubic `b·c²` = `z1 z2²` [COMP, re-verified
exactly], and a 9-term decomposition would need all 9 points on `b = 0` or `c = 0`, contradicting the
`b⁵c`, `a²b³c` terms. Complete replay 43 s (`sextics-small-geometry.py`).
Open question posed by the author: "is there a correspondence or support criterion applicable to this
low-Picard, `End_Hdg(T) = ℚ` K3 that does not require the book's order-three structure, real/complex
multiplication, or the excluded low-degree geometry? A direct construction would need a nonzero Gysin image on
the selected K3 exterior-square component."

### 3.11 `VERIFICATION.md` (repository hygiene; nothing mathematical beyond the above)

Lean 4.33.0 pinned toolchain built (`lake build HodgeComputer.Verify` passed; generated target compiles; "does
not verify a certificate"); content checks pass on main (76 instances, 56 accepted, 20 open); no files changed,
no PR. The new shared `HoldoutCover` library was read and "uses a different fixed cover and does not supply the
source and derivative identities required here".

### 3.12 External, **not posted in #115**: DiamondJM fixed-denominator sweep **[MOD-p, screened]**

`https://github.com/DiamondJM/hodgeHoldout`, `sweep/RESULTS.md` + `sweep/logs/polarized-k3-remainder.log`:
100 "natural" degree-6 denominators in the (X,Y) frame — sums of 2, 3, 4 sixth powers of coordinates
(`pow01 … pow2345`) and every coordinate sub-equation of `F` on 2–5 coordinates containing a coupling
(`Fsub01 … Fsub12345`) — full cleared system mod `p = 268435399`, all **inconsistent** (ranks 2304–2358 out of
2602 columns, 37–60 s each). Superseded by the thread's rigorous degree-6 exclusion; the sweep itself is only
modular evidence and only for Y-even D.

### 3.13 Crashed local session (2026-11-01, this machine; **not posted**) — see §4.3 for the structural claims

Experiments actually completed (`/tmp/hodge/results/polarized-k3-remainder/log.jsonl`, engine v2 after the
`F·N(D)` sign fix at 17:25):
* `D = 1`: inconsistent (i.e. `A ∉ J(F)`, rank 19 — the 19-dim K3 (3,1) block) **[MOD-p, p = 2³¹−1]**.
* `D = f(U) − f(V)` (degree 6, τ-odd): inconsistent, rank 2304, nullity 298, 138 s **[MOD-p]** — consistent
  with the thread's degree-6 theorem.
* Parity-free served-space solver at `D = U0` (degree 1): served dimension 35 but all served targets lie in
  `J(F)` (class dimension 0; K3-block projection 0) **[MOD-p]** — consistent with the degree ≤ 4 exclusion.
* E1 multiplication table `[Â]: R_d → R_{d+6}` (mod 2³¹−1): full rank in every degree `d ≤ 12`; injective for
  `d ≤ 9` (ranks 1, 6, 21, 56, 126, 246, 426, 666, 951, 1246), surjective for `d ≥ 10` with kernels 555 (sym
  294 / anti 261), 1020 (510/510), 1325 (677/648) at `d = 10, 11, 12`. Since a full modular rank bounds the
  char-0 rank from below, **injectivity through d = 9 holds over ℚ/ℂ [PROVEN]** (the thread had this only mod
  101, with a mod-103 drop at d = 9). Kernel representatives at d = 10 have supports ≥ 121 and no small-height
  rational lift [COMP].
* Started but **crashed/unfinished [FAILED]**: fixed-D runs for the octics `f(U)U0² + swap`, `f(U)U0V0 + swap`
  (E1 solvable, rank 5052 of 6006; full 40326 × 7497 systems built; rref never completed), a 24-candidate
  degree-8 list (Fermat octics, cones over octic Hessian-minor curves, polar octics `Σ f_i(U)V_i³ ± swap`,
  perturbations of `f(U)U0²`), two degree-10 polar pairings `Σ f_i(U)f_i(V)`, `Σ f_i(U)f_{i+1}(V) + swap`, the
  parity-free served-space batch for `f(U)`, `f(U)−f(V)`, `f(U)U0`, `f(U)(U0−V0)`, `U0−V0`, and an exact Koszul
  degree-5 Singular exclusion (96 parameters, 1505 equations; never run). None of these produced a result.

---

## 4. What remains unexplored / open sub-problems

### 4.1 Rigorous gaps in the obstruction
1. **Degree 8** is undecided for all D except `D = f(U)U0²` (and, if the sum principle of §4.3 holds, all
   reducible octics). The 315-parameter necessary relaxation cannot decide it; the full bilinear system
   (D-coefficients × N-coefficients) has never been attacked at degree 8.
2. **Degrees 9–12** have never been swept in any form (no fixed-D candidates, no universal argument). The E1
   constraint stops being informative at `d = 10` (surjectivity), so everything is in E2.
3. The thread's degree ≤ 7 exclusion has **not been independently replicated** by anyone else (it is a single
   Codex-assisted packet; the chart-by-chart Gröbner step and the ℤ_(101) flatness argument deserve a second
   implementation).
4. **No τ-odd / non-τ-symmetric denominators** were separately considered in the fixed-D sweeps (the universal
   degree ≤ 7 arguments are parity-free, so that gap is only in the screened sweeps).

### 4.2 Geometry of the specific K3 never examined
5. Exact `ρ ∈ {1,2,3}` unknown. Only the prime 3 was used (`ρ(S̄_{F₃}) = 4` by Tate for K3s). A second good
   prime with `ρ(S̄_{F_p}) = 2` and a discriminant comparison would prove `ρ = 1`. `p = 5` is useless
   (`f ≡ Fermat`, supersingular reduction since `5 ≡ −1 mod 6`); `p = 7` is bad (denominators); `p = 11, 13, …`
   need p-adic point counting (`controlledreduction`-type software), not available locally.
6. Whether `η = [z0²z1²z2²] ∈ H^{1,1}_prim(S₂)` is a ℂ-combination of algebraic classes (only relevant if
   `ρ ≥ 2`) — a period computation (Lairez–Sertöz) never attempted.
7. Classes `E` with `H·E = 4` (totally tangent quartics, `E² = 2`), automorphisms of the plane sextic `C`,
   bitangent/tritangent configurations beyond lines and conics — never screened (would only matter via item 6).
8. Hodge-theoretic *type* of any candidate divisor: nobody has written down a divisor on `S₂ × S₂` or on `X`
   whose resolution could carry `so(T)(1)`; the only known mechanism is the exterior-power theorem from an
   algebraic Kuga–Satake class (Ch. 3A.6), none known for this K3.

### 4.3 Structural claims from the crashed session, **[UNVERIFIED]** (should be checked, then posted)
9. **Equivalence theorem.** Assume `End_Hdg(T) = ℚ` (proven) and `η ∉ NS_prim ⊗ ℂ` (automatic if `ρ = 1`).
   Then `HasCertificate(F, 1, A)` ⇔ `Λ²T ⊂ N¹H⁴(S₂ × S₂, ℚ)` (GHC, coniveau 1, for `Λ²T`).
   Sketch: Zarhin ⇒ `MT(T) = SO(T)`, so `Λ²T = so(T)` is an irreducible ℚ-Hodge structure (`rank T ≥ 19`);
   `Hom_HS(Λ²T, T) = Hom_HS(Λ²T, ℚ(−2)) = 0` by dimension. If `r = Res ω_A ∈ N_Y ⊗ ℂ` for `Y = X ∩ V(D)`, the
   algebraic K3 projector gives `r_K3 ∈ π_K3(N_Y) ⊗ ℂ`; the smallest ℚ-sub-HS `W` with `r_K3 ∈ W_ℂ` projects onto
   a nonzero sub-HS of `Λ²T`, hence onto `Λ²T`, and by semisimplicity + the Hom-vanishing `W ⊃ Λ²T`; transport
   along `Φ` and `S₆ → S₂`. Converse: Ch. 3A.4 projectors + prescribed-support theorem + `HasCertificate.add`
   with the book's supports for the mixed sector. Uses: book §3.1 (N¹ generated by Gysin images), §3A.4
   (algebraic projectors), prescribed-support theorem (Ch. 1).
10. **Dichotomy.** For every legal D, the K3-sector projection of the served (3,1)-space is either
    `⊆ ω ⊗ NS_prim` (dim ≤ ρ − 1 ≤ 2; = 0 if ρ = 1) or `⊇ (Λ²T)^{3,1}` (dim `20 − ρ ≥ 17`). There is no partial
    progress inside the K3 sector; the right diagnostic for any proposed D is the K3-block projection of its
    served space (`k3proj.py` / `pfsolve.py`), not consistency for the single target.
11. **Sum principle at q = 1.** `served(D₁D₂) = served(D₁) + served(D₂) (mod J(F))` for coprime factors, because
    `served(D) = (N_{V(D)∩X})^{3,1}` by the prescribed-support theorem and, by the §3.1 weight argument, the
    (3,1) part of `N_{Y₁∪Y₂}` is the sum of the (3,1) parts of `N_{Y₁}`, `N_{Y₂}` (exceptional components over
    `Y₁ ∩ Y₂` contribute only (2,2)). With the degree ≤ 7 exclusion this kills **every reducible denominator all
    of whose factors have degree ≤ 7** — `f(U)h`, `(f(U) − f(V))h`, `P_i·cubic`, products of linear forms,
    cone products — and explains the spurious degree-8 relaxation points. Only D with an irreducible factor of
    degree ≥ 8 can work. The planned computational check (`served(f(U)U0)` vs `served(f(U)) + served(U0)`)
    never ran.
12. **Shape of a supporting divisor.** Smooth hypersurface sections have `H²(Y) = ℚ` (Lefschetz) ⇒ no (3,1)
    support; the book's §3.1 forces a non-normal or non-rational-singular component. Cones with vertex a line
    over an octic surface `S' ⊂ P³` (`p_g = 35`) would need `H²(S')_tr ⊇ Λ²T(1)` (Hodge numbers (19,172,19)):
    `S'` would be a "Kuga–Satake surface" for this K3; cones with vertex a plane over a curve give
    `H¹(Γ) ⊗ H¹(C')(−1)` (product-of-curves type), incompatible with adjoint type. No low-degree mechanism found.

### 4.4 Constructive gap independent of the K3 sector
13. Even the *supported* part of this target — the six mixed monomials in the sixth-root sector — has **no
    explicit or Lean-existential certificate anywhere in the repository**: all accepted elliptic-companion
    proofs use the μ₃ companion `f + s³ + t³`; the sixth-root sector needs `W_{23,f} = {f + s² + t³} ⊂
    P(1,1,1,3,2)` with `E₂₃ = {a² + b³ + c⁶} ⊂ P(3,2,1)`, `ε = (a, b, ηc)`, differential `db/a` of character η,
    covering curves of weighted degree `e ≤ 2` (three binary e-forms, one 3e-form, one 2e-form: `8e + 5`
    coefficients, `6e + 1` equations, `6e × 8e` freeness matrix), and a μ₆ (not μ₃) translated-diagonal
    contraction. A complete solve needs this **in addition to** the K3 sector (`HasCertificate.add`).

---

## 5. Candidate attack plan (ordered by value / cost)

Sizes below are for the (U,V) frame, parity-restricted unless stated; `rows ≈ #monomials of degree 2d+6 / 2`.

**P1. Independent replication of the degree ≤ 7 exclusion in the (U,V) frame.** Idea: Koszul parametrization
`D = Σ_i a_i ∂_i f(U) + Σ_i b_i ∂_i f(V)` is immediate in the tensor ring `R_f ⊗ R_f` (`uvring.py`); build the
exact quadratic system `B_D(Âa) ∈ L` degree by degree (`koszul_excl.py` already writes Singular input for d = 5:
96 parameters, 1505 equations), then the chart-and-square projective emptiness at d = 7 (120 parameters).
Cost: d ≤ 6 minutes, d = 7 ≈ 5 min (thread: 261 s), memory < 4 GB. Success chance (replicating): high. Value:
makes "degree ≥ 8" a two-source fact; the sum principle (P2) rests on it.

**P2. Verify the sum principle and the dichotomy computationally, then post.** Compute parity-free served
spaces (`pfsolve.py`) for `U0`, `f(U)`, `f(U) − f(V)`, `f(U)·U0`, `f(U)·(U0 − V0)`, `U0 − V0` and check
`served(D₁D₂) ≡ served(D₁) + served(D₂) (mod J)` and `K3_projection_dim = 0` for each. Sizes: d = 6:
≈ 16.8k × (2376 + 226 + 462 target columns) mod p (≈ 2 min each); d = 7: ≈ 26.6k × 4.5k (≈ 10 min).
Chance: high (it is a theorem check). Value: converts the thread's "degree ≥ 8" into "degree ≥ 8 **and**
irreducible factor of degree ≥ 8 **and** (if ρ = 1) equivalent to GHC for Λ²T" — publishable progress on #115.

**P3. Degree-8 decision for irreducible D.** Two routes. (a) Universal: the 315-parameter relaxation restricted
off the reducible locus — saturate the degree-8 quadric ideal by the ideal of the family `{f_U h, (f_U − f_V) h,
P_i·cubic}` (the spurious components), or add the linear "transversality" constraints that kill them; then
msolve/Singular over F₁₀₁ on ≤ 315 variables with 1050 quadrics (the thread's unrestricted run timed out at
600 s with 10 GiB; the local machine has 40 cores / 125 GB and a 512 GB cloud box is available). (b) Fixed-D
screen of natural irreducible octics in the (U,V) frame (Fermat octics `ΣU_i⁸ ± ΣV_i⁸`, weighted Fermat octics,
cones `U0⁸+U1⁸+V0⁸+V1⁸`, cones over the Hessian-minor curves `M_ij(f)(Y)`, `M_ij(U) ± M_ij(V)`, polar octics
`Σ f_i(U)V_i³ ± swap`, `Σ f_i(U)U_i³ + swap`): 40 326 × 7 497 sparse mod-p systems, ≈ 10–20 min and ≈ 4 GB
each (`candidates_d8.json`, 24 entries, ready). Chance of a witness: ≈ 0 by §4.3.9–12; chance of completing an
exclusion: moderate for (a) with the extra constraints, high for (b). Value: sharpens the lower bound to ≥ 9.

**P4. Determine ρ(S₂).** Compute `#S(F_{p^n})` at a second good prime (`p = 11` or `13`; avoid 5 and 7) for
enough `n` to pin the Frobenius polynomial (needs p-adic/controlled-reduction point counting for `w² = f` in
`P(1,1,1,3)`; brute force only reaches `n ≤ 6` at p = 11). If some prime gives `ρ(S̄_p) = 2`, compare
discriminants (Artin–Tate at that prime vs. the lattice `⟨H, R_+⟩` of discriminant `−3·…` at 3 if
`ρ(S̄_3)` were 2 — here `ρ(S̄_3) = 4`, so one needs two primes with rank 2 or a rank-2 prime plus an
independent square-class argument). Chance: moderate (software-bound, days). Value: `ρ = 1` would make the
equivalence theorem unconditional (no loophole) — i.e. a clean statement "this instance is an open case of the
generalized Hodge conjecture".

**P5. Loophole test (only if P4 gives ρ ≥ 2).** Certified periods of `S₂` and of the class `η = [z0²z1²z2²]`
(Lairez–Sertöz `PeriodSuite`/`numperiods`, LLL for NS), then test `η ∈ NS_prim ⊗ ℂ`. If it held, `r_K3` would
lie in the `T ⊗ NS_prim` sector, which **is** supported (curve × K3), and the whole instance would reduce to
constructive work (P6). Chance that the loophole is real: very low (`≤ 2` algebraic directions among 19, and
`η` is a "random" rational-coefficient class); but it is the **only** route to a certificate of reasonable
degree, so it must be settled.

**P6. Constructive certificate for the mixed (sixth-root) part `A_mix = A − A_K3`** where
`A_K3(X,Y) = [(X0+Y0)²(X1+Y1)²(X2+Y2)² − (X0−Y0)²(X1−Y1)²(X2−Y2)²]/64` (both Y-odd). Follow
`crossed-ternary-asymmetric` / `split-ternary-q2` but with the weighted companion `W_{23,f} = {7(f + s² + t³)=0}
⊂ P(1,1,1,3,2)`: find a polynomial curve family through the two quotient points (Hom-scheme search `e ≤ 2`,
`8e+5` unknowns, `6e+1` equations, check the `6e × 8e` freeness matrix mod a small prime, lift), get companion
residue primitives for the two sources of U- resp. V-bidegree (4,·) and (2,·), convert to the characters
`η, η⁵` on `E₂₃` (`db/a`), contract along a translated diagonal of `E₂₃` with the infinity correction, trace,
transport. Chance: moderate but a months-long Lean project; **does not solve the instance** (the K3 part
remains) — but it is required by any eventual solve and would be the first explicit sixth-root-sector
construction for any polarized split sextic (reusable for other polarized-split targets with `A(X,X) = 0`).

**P7. Kuga–Satake algebraicity for this K3 (theory).** Survey and test every known sufficient condition:
Kummer / Shioda–Inose (`ρ ≥ 17`–19): excluded (`ρ ≤ 3`); Paranjape six-line double planes (`ρ ≥ 16`): excluded
(and `C` is smooth, not six lines); Bolognesi–Laterveer order-3: excluded; RM/CM: excluded (`End = ℚ`);
Voisin/Floccari-type results needing special Hodge loci: inapplicable (very general `T`). What is left is a
genuinely new theorem (e.g. algebraicity of KS for *all* degree-2 K3s, which would also settle many other
things). Chance: negligible within this project; record as the real mathematical content of the instance.

**P8. Hidden geometry of C (cheap screens, mod 3 / mod 101).** Totally tangent quartics (`f|_Q = square`,
classes `E² = 2, H·E = 4`), `Aut(C)` via Gröbner on `PGL₃` (expected trivial; implied by `ρ = 1`),
tritangent-line and 6-tangent-conic counts over `F̄₃` already done; sextic contact curves `f = g² + h·k`
with `deg h ≥ 2` (more general contact presentations à la `polarized-contact`). Chance of finding anything:
low (`ρ ≤ 3`); value only through P5.

**P9. Degree-10 structured candidates (low priority).** Polar pairings `Σ f_i(U)f_i(V)`,
`Σ f_i(U)f_{i+1}(V) + swap` (irreducible, τ-even, in the E1 kernel): systems ≈ 85k × 18.9k mod p, 30–60 min and
3–5 GB each. Expected inconsistent by §4.3; worth one run as a control of the dichotomy diagnostic.

**P10. Post to #115** (after P1–P2): the equivalence theorem, dichotomy, sum principle, the ℚ-rigorous
injectivity of `[Â]` through d = 9, the kernel dimensions 294/261 at d = 10, and the negative table; credit
`1337isnot1337`'s packet and the book. Keep a single comment, update it as P3–P5 complete.

---

## 6. Feasibility assessment

**Explicit witness.** Under the proven facts (`End_Hdg(T) = ℚ`, `1 ≤ ρ ≤ 3`, no low-degree special curves) and
the argued equivalence (§4.3.9), a witness exists iff an open case of the generalized Hodge conjecture holds
(unless the transcendental loophole P5 is real). Even if GHC is true (as expected), the only known mechanism
produces supports of astronomical degree (Kuga–Satake surfaces in `A × A`, `dim A = 2^{rank T − 2}`), and any
divisor of degree ≤ 12 would have to carry a weight-2 sub-Hodge structure of adjoint type with `h^{2,0} ≥ 17` —
no hand-written hypersurface section does that. Fixed-D searches at d = 8–12 are therefore expected to return
"inconsistent" and can only sharpen the lower bound. Probability of finding an explicit certificate by search:
essentially zero; probability that one of degree ≤ 12 exists at all: very low.

**Nonconstructive argument.** The repository accepts existential Lean proofs (`HasCertificate`), as in
`crossed-ternary-asymmetric`. For this instance that would require (i) a proof of the support of `Λ²T` — new
mathematics (an algebraic Kuga–Satake correspondence for this K3, or a new divisor on `S₂ × S₂` with the right
Hodge structure), and (ii) its formalization through the prescribed-support theorem (mixed Hodge modules,
Bott vanishing, strictness) — far beyond current Lean/mathlib infrastructure. Not feasible.

**What a Lean proof would require, concretely.** Either an explicit `(D, N, Q)` with
`A D² − D Σ N_i ∂_iF + F(D Σ ∂_iN_i − Σ N_i ∂_iD) = F² Q`, `deg N_i = d + 1`, `deg Q = 2d − 6`, `D ≠ 0`, and
complex coprimality of `D` with `F` (irreducibility of `D` over ℚ is not enough; the accepted proofs use
specialization/Eisenstein-type arguments or `AlgebraicSpecialization`) — kernel replay of such identities is
routine up to roughly the 243 KB `Witness.lean` sizes already accepted, i.e. `d ≲ 10`; or an existential
construction in the style of the elliptic-companion library, which exists only for the μ₃ companion and only
for block-sign/crossed sectors. The honest statement for the record: **`polarized-k3-remainder` is, modulo
`ρ = 1` (to be decided, P4) a reformulation of the generalized Hodge conjecture for `Λ²T` of a very general
degree-2 K3 surface; the useful contributions available now are the sharpened obstruction, the equivalence
theorem, the Picard-rank determination, and a constructive sixth-root-sector certificate for the mixed part of
the target.**

---

## Appendix A. Independent checks performed for this dossier (exact unless noted)

* `expand(F − ½(f(U) + f(V))) == 0`, `F(X,−Y) = F`, `A(X,−Y) = −A`, `A(X,X) = X0²X1²X2²` (sympy).
* Smith normal form of exponent-difference lattices: F → `(1,1,1,1,2,0)` (projective diagonal group ℤ/2 = ⟨τ⟩);
  f → `(1,1,0)` (trivial). No coordinate permutation preserves f.
* `7f ≡ 2(z0⁶ + z1⁶ + z2⁶) (mod 5)`; all mixed coefficients of F divisible by 5.
* `Hilb(R_f) = (1,3,6,10,15,18,19,18,15,10,6,3,1)` over ℚ (python-flint `fmpq_mat`); `rank J(f)_6 = 9`,
  `rank(J(f)_6 + η) = 10` ⇒ `η = [z0²z1²z2²] ≠ 0` in `R_{f,6}`.
* Middle catalecticant of f: rank 9 over ℚ, kernel spanned by `z1 z2²` (= the packet's `b c²`).
* Sector counts: 226 Y-odd sextics = 28 pure + 198 mixed antisymmetric (U,V) combinations; `dim J(F)_6 = 36`,
  antisymmetric part 18; `208 = 19 + 189`.
* Unknown counts in §1.3 from binomial coefficients; d = 8 matches the engine's `N 6006, Q 1491`.

## Appendix B. Pointers

* Thread: https://github.com/hodge-computer/sextics/discussions/115 ; packet gist
  https://gist.github.com/1337isnot1337/d504e828c4f50d8dade0687792ea9ec0 ; local copy
  `/tmp/hodge/downloads/k3-remainder-gist/packet/` (`REPORT.md`, `DEGREE_SIX.md`, `DEGREE_SEVEN.md`,
  `DEGREE_EIGHT.md`, `PICARD_BOUND.md`, `PICARD_LIFT_OBSTRUCTION.md`, `ORDER_THREE_EXCLUSION.md`,
  `FROBENIUS_ORDER_THREE.md`, `FROBENIUS_GALOIS.md`, `HODGE_ENDOMORPHISMS.md`, `SMALL_GEOMETRY.md`,
  `VERIFICATION.md`, scripts `sextics-*.py`, `count-k3-reduction-3.cpp`, Singular/msolve logs).
* External sweep: https://github.com/DiamondJM/hodgeHoldout (`sweep/RESULTS.md`,
  `sweep/logs/polarized-k3-remainder.log`); local `/tmp/hodge/ext/hodgeHoldout/sweep/`.
* Crashed-session tools: `/tmp/hodge/results/polarized-k3-remainder/` (`uvring.py` tensor Jacobian ring,
  `pfsolve.py` parity-free served space + K3-block projection, `k3proj.py`, `e1sweep*.py`, `koszul_excl.py`,
  `gen_candidates.py`, `candidates_d8.json`, `candidates_d10.json`, `runD.py`, `log.jsonl`); engine
  `~/devel/sextics-work/hodge/` (`fixed_d.py` v2 — results before 2026-11-01 17:25 are unreliable).
* Book: `docs/03a-split-and-k3.md` §§3A.4–3A.6 (decomposition, polynomial test, Kuga–Satake criterion,
  order-three family, "direct boundary to the elliptic trick"); `docs/03-geometric-support.md` §3.1
  (Gysin/N¹, (3,1) needs non-normal or non-rational singularities), §3.6 (support ⇒ every q);
  `docs/01-certificates.md` (Euler gauges, prescribed-support theorem, `D ↦ D + tF` gauge at d = 6);
  `docs/04-known-locus-and-frontier.md` §4.1 (polarized family: normal form 11, saturation 24, quotient 7),
  §4.2 (rank-1073 supported substructure; "No unconditional support theorem for that remaining piece").
* Companions: `data/split-ternary-q2/solution.md` (same f), `data/crossed-ternary-asymmetric/solution.md`,
  `data/polarized-contact-coupled-d2/{solution.md,witness.json}` (D = X2² − Y2², splitting-line mechanism).
