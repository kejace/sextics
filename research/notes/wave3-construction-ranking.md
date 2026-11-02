# Wave 3 — construction-side synthesis and ranking (16 open lanes)

STATUS: COMPLETE (2026-11-02, UTC). Synthesis only — no new witness, no new computation beyond re-reading receipts and one
literature check. Author: construction-synthesis agent (Claude) for `kejace`. Nothing posted, no git.
Log of this synthesis: `/tmp/hodge/results/wave3-construct/log.jsonl`.

Inputs read in full: every `wave2-*-construct*.md` (incl. part files), `wave2-theory-{cyclic,deformation,generic,permutation}.md`,
`wave2-trio-theory.md`, `probe-four-plus-two-jacobian.md`, `probe-polarized-k3.md`, book maps `book-02/03/03a/04`,
the accepted catalogue `accepted-group0..10.md` (cross-cutting and relation-to-open sections), `wave2-status.md` §§2, 6, 7,
the referee reports under `/tmp/hodge/posts/*/referee-report.md`, the attack plans of the dossiers whose construct lens is a
stub (klein, quadric-contact, quadratic-cubic, rank-one, symmetric-sextic, four-plus-two, one-reflection-skew),
the six shared-library READMEs under `formal/lean/HodgeComputer/Constructions/`, `formal/README.md`, and `hodge/rjob.py`.
State of the queued construct jobs (checked with `rjob status` at write time): `construct-dp-j1j2`, `construct-off-paired-served30`,
`construct-rcr-all` RUNNING; `construct-off-paired-mw` FINISHED (msolve `[-1]` for both signs: no polynomial Mordell–Weil
section of `Y_∓`, see §2.8); `construct-dp-fibre`, `construct-epc-fanoplanes`, `construct-spp-vert/k7var`, `construct-ctr-ztsmooth`,
`construct-kcc-sigmaXQ`, `construct-V3`, `construct-dense-rational-cones/univ12`, `construct-polarized-k3-remainder-w23seed`
have no output (never started or starved).

Conventions. Variables `(X0,X1,X2,Y0,Y1,Y2)`, `τ(X,Y) = (X,−Y)`, residue character of a monomial `Z^e` under diagonal weights `w`
is `Σ_j (e_j+1) w_j` (book §3.3). "Served" = in the Gysin image of a resolution of `X ∩ V(D)`; E1 = `A·D ∈ J(F)`.
Probabilities below are subjective: **P(solve)** = probability that the route, pursued for one wave (≈ one agent-week plus
remote compute), ends in an accepted Lean proof of the lane's protected theorem; **P(decide)** = probability that the named
next computation returns a definite yes/no for the route. "Verified" means a receipt exists in the lane directory or a
referee replay; "cited" means taken from a thread post without replay.

## 0. Rubric and the five mechanisms as posed

For each lane I score the five mechanism classes of the task — (i) a new cover creating a reflection/exchange (book §2.7 /
§4.6 descent), (ii) a Fano or `p_g = q = 0`-fibration companion of a subgroup quotient (book §3.2–3.4, §4.2–4.3),
(iii) a Galois-closure correspondence (§3.1 averaging ideal on intermediate quotients), (iv) `q`-raising from a served
sub-family (§4.5) or deformation/continuation, (v) companion-support modification (joining divisors, twisted slices,
`HasCertificate.add` bundles) — and keep the one or two that are not excluded by an exact computation.

Lean-feasibility scale used in the entries: **L-A** = explicit witness, checkable by the Kronecker big-integer route
(`lean-capacity.md`: 3.1 M raw products in 181–263 s / 11–13.5 GB; `ring` template dies above ≈ 30 k products);
**L-B** = existential proof assembled from the shared modules (`PolynomialCurves` lifting/normalized data on a uniruled
auxiliary variety + a new cover/quotient field layer + `CertificateOperations`), i.e. the holdout / reflection-coupled-fano /
cyclic-triple-fano architecture; **L-C** = existential proof needing a new library layer (elliptic companion for a non-Fermat
curve, higher residues on a new family); **L-D** = no formalizable route in the present library (Hodge-theoretic existence
only: Kuga–Satake, elliptic integrals, GHC-type statements).

## 1. Cross-lane facts that fix the ranking (all with receipts; see the per-lane entries for files)

F1. **Residue character / averaging ideal.** A quotient-type correspondence `X ← Z → Z/K` serves the target's representation
`ρ` only if `ρ^K ≠ 0` for a `K` whose quotient is uniruled or `p_g = q = 0`-fibred (book §3.1, §3.2). Every open lane is the
complementary character of its solved sibling or of the only available quotient; this was re-verified lane by lane
(cycle-q2 part 2: all 226 subgroups of the 216-element group; cyclic-triple part 1: all 128 subgroups of the order-108 Galois
closure; eight-plane part 2: `S4 × Z/2`; reflection-coupled §2.2: `B3 × Z/2`; polarized-k3 A5: `D6`; off-paired §G: `D6`).

F2. **ω-descent principle** (eight-plane part 2, Theorem; same argument in dense-rational part 1, symmetric-sextic). `X` is
Calabi–Yau; for any finite group acting on a cover `Z → X`, the holomorphic 4-form `ω` and the target class `α` span
isomorphic representations up to the sign of a `τ`-type element. A quotient that carries `α` therefore carries `ω` or
`ω²` unless the group contains an element acting as `−1` on `ω`, `+1` on `α` **with a codimension-one fixed locus** (a
created reflection). On `X` itself only `τ` separates `ω` from `α`, and `τ` has codimension-3 fixed curves (age 3/2), so
every honest quotient has `κ ≥ 0`. Consequence: only (i) covers with a *created* reflection, or non-quotient
correspondences (divisor geometry on `X`, twisted slices, Prym/elliptic transports) can reach any open target.

F3. **Reid–Tai audits**: in every open lane every nontrivial element that fixes the residue has age `≥ 1` at every fixed point
(polarized-k3's junior `h, h^5` serve only the joining sector); so no `X/K` and no elliptic/K3-twisted product `(X × Y)/K`
is uniruled (theory-permutation O1; off-paired §D).

F4. **§2.7-type exchange covers are excluded exactly** on cycle-q2 (60 coordinate choices + Fermat-curve automorphisms),
cyclic-triple (Singular elimination: only `X1`, `Y0` directions, both blocked by two couplings), dense-rational
(`V1 = V2 = ∅` coordinate-free, msolve leading ideals), eight-plane (Lemma A, exact ranks 6/6), off-paired (no-hub theorem),
short-power (subset-sum covers, row-space argument), reflection-coupled (six coordinate pencils), klein (no linear fibrewise
exchange; `Aut` of the generic fibre of `X → P²_Y` is exactly `Z/14`, so no birational one either), diagonal-permutation
(pure-`X`/`Y` partners excluded; mixed partner msolve still running), quadric-contact and four-plus-two (no merging Fermat
pair). Not analysed: symmetric-sextic, rank-one, polarized-k3 (irrelevant there: `Aut(C_f) = 1`), one-reflection-skew
(DiamondJM: rescaled exchange needs `α⁶ = −1`, finite base points — not a cover).

F5. **Surface-fibration companions**: minimal fibre `p_g` over all natural fibrations is 2 (cyclic-triple, cycle-q2's S3),
1 (diagonal-permutation), 3 (dense-rational over `P²_X`, reflection-coupled), 20–100 (short-power, rank-one, off-paired's
pencil surfaces `p_g = 45`). None is 0; the book's §3.2 theorem never applies to an open lane.

F6. **(3,1)-support constraint and sector-variation bounds**: any support has a non-normal or non-rationally-singular
component (book §3.1); very general supports need `Σ h^{2,0} ≥ 208` (dense lanes, §4.4), `≥ 126` (one-reflection residual
sector), `≥ 67` (cyclic-triple, exact `K = K2 = 0` computation). Low-degree denominators are therefore excluded for
structural reasons, consistent with the rigorous `deg D ≥ 7` (all `q = 1` lanes), `≥ 8` (polarized-k3), `≥ 6` (cycle-q2)
bounds of the exclusion lenses.

F7. **Lean capacity** (lean-capacity.md): explicit witnesses of degree `≥ 7` on dense lanes need the Kronecker checker;
existential proofs reuse `PolynomialCurves` (generic), `AveragingProjectors`, `CertificateOperations`, `ResidueTransport`;
`HoldoutCover` binds the holdout sextic (a new cover needs its own field layer); `EllipticCompanions` is stated for the
Fermat cubic only (`ScalarEmbedding` allows a parameterized companion cubic); `HigherResidues` is the `q = 2` five-coordinate
engine; `TransitiveCubic` is the six-coefficient paired-sign family. No module integrates an elliptic differential.

F8. **Literature check (today)**: algebraic Kuga–Satake correspondences for K3 surfaces are known for Picard rank `≥ 17`,
Paranjape's six-line double planes, Bolognesi–Laterveer's order-three family and Varesco's RM families; nothing covers a
very general degree-2 K3 with `ρ ≤ 3`, `End_Hdg(T) = Q`. The K3 part of `polarized-k3-remainder` stays an open
GHC statement (sources: arXiv 2203.09778, 2304.02519, 2210.02948).

## 2. One-line summary table (details in §3; ranking in §4)

| lane | best mechanism left | auxiliary object | decisive computation | P(decide) | P(solve) | Lean |
|---|---|---|---|---|---|---|
| cycle-q2 | M*: S3-transport, non-fibrewise rational current | RES over `P¹_T` with Kummer twist, `V_{ω²}` rank 12 | twisted Gauss–Manin over `C(c1,c2)`, class of `σ` | 0.6 | 0.04 | L-A (q=2 witness) |
| cyclic-triple-residual | j=0 elliptic exceptional surfaces / scrolls over `𝓕` | Fano variety of lines `𝓕` of `W ⊂ P(1⁵,3)` | compute `𝓕`, twisted irregularity of line surfaces | 0.5 | 0.03 | L-A only |
| dense-rational | none (all universal divisors have budget 0) | — | — | — | <0.01 | — |
| diagonal-permutation | second reflection on `(01)`-pencil fibres ⇒ §4.6 | `X_μ = X ∩ {r = μt}`, `W_t = X_μ/s` | J2/J2f msolve (running) | 0.9 | 0.02 | L-B |
| eight-plane-configuration | none (ω-descent closes plane-pair covers) | — | Fano-scheme check (formality) | 0.95 | <0.01 | — |
| four-plus-two-residual | core-side support via exact 5-variable reduction | `Y_H: U⁶ = −H ⊂ P⁴`, core surfaces `V(D) ∩ Y_H` | adjoint-budget design + reduced solver sweep | 0.7 | 0.04 | L-A |
| klein-coupled-character | cone over `B = {K = 0}` (char-3 forms) | `Z = X ∩ V(g(X))`, fibres `S_y` | Gysin test of cones, base degree 9–12 | 0.6 | 0.02 | L-A |
| off-paired | descent-through-reflection via birational `φ` of `W` | `Z = X ×_W X^φ`, `W = X/σ2` Fano | MW sections with denominators; rigidity of `W` | 0.5 | 0.03 | L-B |
| one-reflection-skew | fqlx/benediamond twisted CM companion (E4/E5) | quartic cover `R(β) = 0`, `E_β`, base cover `B''` | κ(B̄''), monodromy of `L`, period test | 0.7 | 0.04 | L-B/L-C |
| polarized-k3-remainder | Part B: `A_mix` via `W_{23,f}` (partial only) | `W_{23,f} ⊂ P(1,1,1,3,2)`, seed over `F_5` | none for `A_K3` (open GHC) | — | 0.005 | L-C (partial) |
| quadratic-cubic-pullback | non-ACM curves / imprimitive cover `C'` | `X → C` (Fermat cubic 4-fold), degree 32 | monodromy of `φ` (E7), adjoint screening at d=10 | 0.8 | 0.02 | L-A |
| quadric-contact-residual | cover `T² = δ` with lifted exchange; NL-transverse supports | `X̃: T² = X0² − Y0²` | averaging-ideal bookkeeping (E8) in minutes | 0.9 | 0.03 | L-B if positive |
| rank-one-biquartic | correlated two-pencil support; degree-12 split loci | `Γ ⊂ P²×P²`, `p ∈ S(2,0), q ∈ S(0,2)` | E1/E2 of degree-12 pull-backs (sparse solver) | 0.6 | 0.02 | L-A |
| reflection-coupled-residual | M5 twisted slices / M1 hidden fibre reflection | `W = X/r ⊂ P(2,1⁵)`, conic family `e = 2` | Hodge bundle `λ` on slices + Gysin; `construct-rcr-all` | 0.5 | 0.03 | L-D / L-B |
| short-power-presentation | vertical divisors on the `(ε_X, ε_Y, t)` base | `D_res` (deg 10), ring `C[e(X),e(Y),L,M]` | E1 of vertical families d = 7..12 (job never ran) | 0.9 | 0.02 | L-A |
| symmetric-sextic-residual | cones `Δ_X·s(X)` over symmetric plane curves | pencil `S_λ ⊂ P³` with `sgn_Y` three times | Hodge-bundle twist on `C = {Δ_X s = 0}` (E6) | 0.6 | 0.03 | L-A |

## 3. Per-lane entries

Format: target and sector → mechanisms kept (with the explicit auxiliary variety) → hypotheses verified / not verified →
decisive next computation with size → probabilities → Lean route. Excluded mechanisms are listed once with their receipt.

### 3.1 cycle-q2 (q = 2, `A = X0²X1²X2³Y0²Y1²Y2` on the cyclic three-link sextic; residue character `ω²` on `G/K`, G-sector (0,6,26,6,0), `α^{3,1} ≠ 0` by fqlx)

Excluded with receipts: (i) all 60 coordinate §2.7 covers, non-coordinate forms by Fermat-curve automorphisms (part 3);
(ii)/(iii) all 226 subgroups of `Γ = G ⋊ ⟨R⟩` (order 216): 154 reach `ρ = Ind χ_A`, all have `κ = 0` and their only surface
fibration has `p_g = 2` fibres (part 2, `w2c_gamma.json`); (iv) the q-raising identity is verified but all 61 degree-6
factorisations `A = M'·B'` either keep the character on the sibling or force the S3 base change (part 4, `w2c_factor.json`);
vertical divisors of the `(c1,c2)`-fibration cannot support `α` (part 5, Leray + injective KS).

**M-a (kept): M\* — S3-transport with base absorption, i.e. a non-fibrewise rational primitive.** Auxiliary data (part 5,
`w2c_s3.py`): cubic `U = a x³ + b y³ + c z³ + d x²y + e y²z + f z²x`, `(a,b,c) = (1+t³, 1+c1⁶t³, 1+c2⁶t³)`,
`(d,e,f) = (t/2)(c1, c1c2, c2)`, `t = v² = (Y0/X0)²`; Jacobian `η² = ξ³ + s4 ξ − s6/4` with `s4 = c1c2t³(c1⁶c2t³ + c1²c2⁶t³ + c1c2²t³ + c1² + c1c2² + c2)/8`,
`s6` the standard degree-18 polynomial in `t³`; `S3(c1,c2)` over `P¹_t` (χ = 3, p_g = 2, κ = 1), `RES(c1,c2)` over `P¹_T`, `T = t³`.
The target's `(3,1)`-part sits fibrewise in the sector of the untwisted form `ω₀ = dt∧dξ/η`; the remaining equation is
`[σ] = 0` in `H²_dR(B°, V_{ω²})`, `V_{ω²} = H¹(P¹_T, R¹_RES ⊗ L_{ω²})` of rank 12, type (1,10,1), over `B = C²_{c1,c2}` minus the discriminant.
Verified: (H1) `s4, s6 ∈ C[t³]`; (H2) `deg s4 = 6, deg s6 = 18`, minimal at `t = 0, ∞`, discriminant degree 36, `j` non-constant;
(H3) characters of `ω₀, ω₁`; cited: (H4) DiamondJM — the fibrewise (0,2)-component of the transported source lies in
`J37 + ∂_{c1}G·S19 + ∂_{c2}G·S19`; (H5) benediamond — vertical remainder nonzero (5-adic), so no fibrewise primitive exists.
Not verified: the Gauss–Manin connection of `V_{ω²}` over `C(c1,c2)`; the class `[σ]`; whether poles along a horizontal
non-normal divisor `Z` (Noether–Lefschetz jump of the `ω₀`-sector) are needed.
Decisive computation: Picard–Fuchs/Griffiths–Dwork reduction for the Kummer-twisted RES in the frame `{t^j T^k dξ/η}`:
two `12 × 12` connection matrices over `Q(c1,c2)` (entries of degree ≲ 10²), the source `σ` from the second variation
`2AΩ/F³ = ∂_{s1}∂_{s2} Ω/F_s|_0` transported through the cubic (DiamondJM's frame or the lens's `w2c_s3.py` data), then
solve `∇Ψ = σ` for `Ψ ∈ V ⊗ Ω¹_B` rational (linear algebra over `Q(c1,c2)` with pole-order bounds along the discriminant).
Size: a few CPU-hours of Singular/sympy on the remote; the only two-variable D-module computation in the programme.
Outcome A (rational `Ψ` exists): fqlx's pole-three clearing adapter turns `Ψ` into an explicit `q = 2` certificate with
denominator the polar locus of `Ψ` (degree unknown, probably large) → **L-A** via the Kronecker checker.
Outcome B (`[σ] ≠ 0` on `B°∖Z` for every tested `Z`): the support must be a horizontal non-normal divisor (dossier P5);
the lane is then a search problem in the sibling spaces.
P(decide) 0.6 (the connection is computable; the pole bookkeeping is the risk), P(solve) 0.04.
No finite-field incidence seed exists for any auxiliary surface carrying the `ω²`-sector (dominance fails: `S3` has `κ = 1`),
so `PolynomialCurves`/`HoldoutCover` cannot be the engine here.

**M-b (secondary): q-raising from a character-free sibling** — the five families with `G' = K` (`M' ∈ {X1²X2²Y0², X0X2Y0Y1²Y2,
X0²Y0²Y1², X0²X2²Y0², X0²X1X2Y1Y2}`) give `q = 1` problems `(F + sM', B')` with `B'` invariant under `Sym(F_s)`, but `F_s`
leaves the 26-coefficient Γ-family (`M'` is not `h`-invariant), so no sector theorem applies; only a fixed-degree search
(`deg D ≥ 7`) at a few `s` values is left. P(solve) 0.01.

### 3.2 cyclic-triple-residual (`A = X0³Y0³`, `ζ3`-sector exact Hodge numbers (41,119,26), `τ`-invariant)

Excluded with receipts: all 128 subgroups of the order-108 group of the Galois closure `Z̃ = X ×_{P⁴} Y` (38 serve the class,
none with negative canonical degree; `c2_abstract.json`); abelian-extension lemma (no tower of cyclic covers helps); §2.7
exchange covers for every linear form (Singular char 0: only `e_{X1}`, `e_{Y0}` directions, both blocked by
`X0²Y0Y1Y2²`/`X0X2³Y1Y2` resp. `X0X1Y1³Y2`/`X0X2³Y1Y2`); surface fibrations of `X/K`: minimal `p_g = 2`; §4.6 descent would need
`Z̃/⟨r,σ⟩` (deg K = 1) supported. Sector variation: `Def_G·ξ = 41`, `Def_G²·ξ = 26`, `K = K2 = 0` ⇒ very general supports need
`Σ h^{2,0} ≥ 67` (exact over Q by rank semicontinuity; `c5_sector_vhs.json`, `c5b_K2.json`).

**M-a (kept): isotrivial j = 0 elliptic exceptional surfaces along triple-point curves, equivalently scrolls of lines over
irregular surfaces of the Fano variety of lines.** `W = X/μ3 ⊂ P(1,1,1,1,1,3)` (coordinates `X0,X1,Z = X2³,Y`), branch divisor
`B = {Z = 0} ≅ {H = 0} ⊂ P⁴`; for a weighted line `L ⊄ B`, `E_L = π⁻¹(L) = {X2³ = c_L(s,t)}` is a j = 0 cubic with `μ3` acting as CM;
`V_1 ⊕ V_2 ↪ H³(𝓕, L_χ ⊕ L_χ̄)` with `𝓕` the 3-dimensional family of lines (rbajaj5's free line `[0:s:0:is:t:it]` is one member).
Supports: `T_S = ∪_{L ∈ S} E_L` for a surface `S ⊂ 𝓕` with twisted irregularity `q_χ(S') ≥ 41`; the triple-point-curve
denominators of part 4 (`D` of order 3 along a plane `Π ⊂ {X2 = 0}`, leading cubic `c3(n1,n2) + aX2³`) are their degenerate shadow.
Verified: Hodge bookkeeping, degree-9 families (120 `(plane, sector)` pairs, dims 29–66, four random members inconsistent at
two primes), the `≥ 67` bound. Not verified: `𝓕` itself (never computed), the twisted irregularity of any surface of lines,
smoothness of `Z̃` (job `construct-ctr-ztsmooth` never ran).
Decisive computation: (1) compute `𝓕` — incidence of five linear forms and one binary cubic in the equation of `W` modulo
`GL2`, ≈ 17 unknowns over `F_p` (Gröbner, hours on the remote, 8 cpus); (2) for the surface `S_B ⊂ 𝓕` of lines tangent to `B`
and for `S_Σ` (lines meeting a fixed surface `Σ = W ∩ {ℓ = ℓ' = 0}`) compute `q_χ` of a smooth model (point counts or
Hodge-bundle degree on a pencil of curves in `S`); (3) if `q_χ ≥ 41`, derive the equation of `T_S` (elimination) and run the
sparse sector solver. P(decide) 0.5, P(solve) 0.03. Lean: **L-A only** (an explicit witness); the elliptic-fibre integration
has no library counterpart (elliptic integrals, not polynomial sources).

### 3.3 dense-rational (`A = X0³Y0³`, 236-term even sextic; sector (0,208,866,208,0), expected simple)

Everything definable from the equation has zero (2,0)-budget: `V1 = V2 = ∅` (no exchange cover, no quasi-reflection in any
coordinates; msolve leading ideals, 4.9 s / 374 s), `V3` (linear reflections) pending; ternary cones serve no primitive class
(Lefschetz on the net of `P³`'s through `P²_Y`) and fail E1 for `k ≤ 8` exactly (9..12 pending); `R`, `R'` (residual tangency,
`D_{R'} = (F24 + F06)² + 4F60F06`, degree 12), `det Hess F`, bihomogeneous sections, products: budget 0 by Lefschetz on `W = X/τ`
(smooth cubic cover of `P² × P²`) and the equivariance reduction (only `H²((Y/τ)~)` can serve the invariant class).
No mechanism survives; a support must be a non-universal singular divisor of degree `≥ 10` whose singular locus carries a
1282-dimensional transcendental structure. No auxiliary variety, no seed. P(solve) < 0.01. (The pending jobs are consistency
checks of the structural argument, not routes.)

### 3.4 diagonal-permutation (`A = X0³Y0³ − X1³Y1³`, standard representation of the simultaneous `S3`; `B = T0 + T1 − 2T2` is the `(01)`-invariant partner)

Excluded with receipts: all quotients (ages ≥ 1, `K ≡ 0`, table of `h^{3,1}(X/K)`); `⟨s,τ⟩`-stable linear-section fibrations
have fibre `p_g ≥ 1` (exact character count, J3); Galois closure of the cubic cover with the commuting `S3^sym`: every `K` with
`(V_std ⊗ V_std)^K ≠ 0` gives a non-uniruled quotient; pure-`X`/`Y` §2.7 partners excluded by hand; 13 natural degree-≤6
divisors (Shioda-type cones, sign quadric, `τ·s` cones) serve nothing (X1/J5/J5b); fixed-degree deformation from the
three-reflection sub-families is obstructed for `d ≤ 7` (kernel bundle rank 0) and `d = 8` (the D8 line, inconsistent).

**M-a (kept): a second reflection on the `(01)`-pencil fibres after base change ⇒ §4.6 descent** (holdout architecture).
Auxiliary: pencil `{r = μt} ⊃ P³_0 = Fix(s)` (`r = X0 − X1`, `t = Y0 − Y1`), fibres `X_μ` sextic threefolds on which `s` is a
reflection with fixed surface `S`; `W_t = X_μ/s` the CY3 sextics in `P(1,1,1,1,2)` of #177; a second (hidden or created)
reflection on `X_μ` would make `Z/⟨s̃,σ⟩` fibrewise Fano in `P(1,1,1,2,2)` and give `H⁴(X)^{(01)} ⊂ N¹` (hence `A` via
`B`'s `(01)`-equivariant data). Verified: `S` is smooth on every `X_μ`; the fibre quotient needs exactly two reflections.
Not verified: J2 (mixed exchange partners on `P⁵`: `ℓ_v ∣ ∂_v³F`, `F(v) ≠ 0`, six msolve charts; job `construct-dp-j1j2`
RUNNING since 20:23 EST) and J2f (hidden reflections of `F_μ` at `μ = 3/7`; `construct-dp-fibre` has no output — resubmit
`w2c_fibre.sh`). Decisive: let J2/J2f finish (msolve, 4 cpus, hours; numerical screens found nothing). P(decide) 0.9,
P(nonempty) ≈ 0.02 — if nonempty the route is **L-B** (new cover fields + `PolynomialCurves`).
**M-b**: the served-space criterion `A ∈ Q[S3]·served(D) + J_6^-` applied to `D8` and `(01)`-stable degree-8/10 candidates
(`w2c_served_orbit.py`, `w2c_d8_served.py` ready). P(solve) 0.02.

### 3.5 eight-plane-configuration (`A = X0³Y0³`, `F = Σ Q_i K_i`, `Q_i = X_i² + Y_i²`; τ-invariant (3,1) sector 208)

Excluded with receipts: all §2.7 covers (Lemma A; the 15 coefficient quadrics of `F_{4,2}`, `F_{2,4}` span 6 dimensions, Gram
ranks 3); the four quartic plane-pair covers `X ⇢ Π_ε × Π_{−ε}` (`F(sp + tp') = st·q`, `S4` monodromy, reduced discriminant,
non-isotrivial elliptic carrier `w² = q`): every intermediate quotient of the `S4 × Z/2` closure carrying `α` has `κ ≥ 0`
(ω-descent theorem, (H5) checked up to isolated points); ramification divisors `X ∩ V(Σ ε_j R_j F)` are excluded degree-6
Jacobian denominators; `Aut(X) = ⟨τ⟩` among block-monomial maps (exact over `Q(i)`), conditional on the Fano scheme of planes
being the eight planes (job `construct-epc-fanoplanes` never ran; recipe in part 3 §6: 20 msolve charts, 9 unknowns,
28 sextics, `p = 2147483629`, 8 cpus, ≈ 1.5 h). Nothing remains except a non-linear created reflection, for which the
arrangement offers no candidate (planes rigid, type (2,2); the self-correspondence `x ↦` other residual points acts alike on
`ω` and `α`). P(solve) < 0.01. Only action: run the Fano-scheme job to make `Aut(X) = ⟨τ⟩` unconditional.

### 3.6 four-plus-two-residual (`A = X0³Y0³`, `F = H(X0,X1,Y0,Y1) + B(X2,Y2)`; sector `M' = H³(T_H)^{ι1=−1}_χ ⊗_K H¹(P)_{χ⁻¹} ⊕ conj`, rank 288, (64,160,64), level 2)

Structure (verified): `T_H: v³ = H ⊂ P(1,1,1,1,2)` smooth CY threefold, `χ`-sector (1,68,35,0); `C_B: w³ = B` genus 4,
`J(C_B) ~ E' × P` with `P` the 3-dimensional Prym; `P` is **absolutely simple with `End⁰ = Q(ζ3)`** (irreducible Frobenius
sextics at 28 of 36 split primes `p ≤ 400`, distinct Frobenius fields) ⇒ the 3A elliptic contraction is unavailable; `W_H =
{H + s³ + t³ = 0} ⊂ P(1,1,1,1,2,2)` is a terminal Gorenstein Fano fourfold (`−K = O(2)`, `h^{3,1} = 56`) supporting
`[H³(T_H) ⊗ H¹(E)]^{μ3}`, but the target needs a second supported factor on the binary side (attack1 §1.4 dimension count:
`Γ ∩ (D_H × W_B × X)` has dimension 4). Rigorous: `deg D ≥ 6`; core kernels exact (k = 6: 20/13, 7: 30/30, 8: 41/34, 9: 38/38).

**M-a (kept): core-side support through the exact 5-variable reduction.** `X̃ = Bl_{S_H} X → P¹_{(X2:Y2)}` with fibres
`Y_t ≅ Y_H: U⁶ = −H ⊂ P⁴`; a core-only `D(X0,X1,Y0,Y1)` cuts `Z_D = ∪_t V(D) ∩ Y_t`, and `Z_D` supports `[A]` iff the fibre
class `[aU] ∈ H^{2,1}(Y_H)_{ζ6²}` is supported on `S_D = V(D) ∩ Y_H` (the `μ6`-cover of the surface `V(D) ⊂ P³` branched along
`V(D) ∩ S_H`). Algebraic form (proved direction): a core field `ñ` with `ñ(H) = aD`, `ñ(D) = D·div ñ` lifts to a strict
6-variable certificate `N = ñ + mE`; the reduced solver (`threefold.py: reduced_test`, U-degree sector, ≈ 20× smaller) detects
the known core-only witnesses of `family-hesse-web` and `f3-block-x0-quartic-strict`. Verified negatives: every rref basis
element and random member of the core kernels for `k ≤ 8`, random members at `k = 9, 10`, cones with vertex on `S_H` over
plane curves of degree 2–6, six-plane unions, the split neighbour (`deg ≤ 8`). Not verified: any **designed** singular core
surface: `V(D)` must be singular so that `S̃_D` acquires irregularity in the `ζ6²`-character — the classical sources are
(a) non-rational transversal singularities along curves `Γ ⊂ S_H` (cones with vertex on `S_H` over curves of genus `≥ 1`;
double surfaces along curves of `S_H`), (b) Zariski-type branch curves `V(D) ∩ S_H` with special singularities (cusps on
conics give irregular cyclic covers), (c) unions of surfaces through curves of `S_H` (Mayer–Vietoris joining). The adjoint
criterion of the quadratic-cubic lens (`h^{2,0}(S̃)_χ = h¹(I_adj(d))` sector-wise) prices each configuration in seconds.
Decisive computation: enumerate configurations of degree 9–14 with positive `ζ6²`-adjoint budget (plane sections of `S_H`
through chosen points/lines, the genus-4 curves `C_f = {X0 = X1 = 0} ∩ Y_H`, `C_g`, tangent developables), impose them as
linear conditions on `D` (`D ∈ I_Γ^m` restricted to `Y_H`), intersect with the exact core E1 kernel, and run `reduced_test`
on every member (5-variable systems of 3–10 k unknowns: seconds to minutes each; a full sweep is hours on 8 cpus).
Complementary theory step: Frobenius on `H³(T_H)^−_χ` (surface/curve point counts of `S_H`, `C_f`, `C_g`) to see whether the
core Weil-type structure (signature (32,16)) has `End ⊋ Q(ζ3)` — if it does, an explicit Prym/Jacobian correspondence
(the missing "second supported factor") is plausible and would name the surfaces.
P(decide) 0.7, P(solve) 0.04. Lean: **L-A** — a strict 5-variable identity lifts to a strict 6-variable witness (`Q = 0`),
checkable by the Kronecker route; no new library layer.
**M-b**: reducible `c·X2`, `c·Y2` products (joining across the fibre `Y_0`): 33 candidates in the pending sum-principle job.
P(solve) 0.01.

### 3.7 klein-coupled-character (`A = x³t³`, `σ`-character 3 of `Z/7`; sector (26,123,37); `Sym = Z/14`)

Excluded with receipts: every quotient of `X` (averaging over `σ^k` kills `A`; `X/τ` has `κ = 0`); the §4.3 weighted Fano
`𝒲 ⊂ P(2,2,2,1,1,1)` serves only character 0; no linear fibrewise exchange on any `Y`-ratio fibration (the `σ⁻¹`-eigenvalues
`ζ⁶, ζ⁵, ζ³` are absent among `x, y, z`), and **no birational one**: `Aut(S_(1,2,3)) = Z/14` exactly (msolve char 0 on the
Hessian-degeneracy locus; monomial automorphisms = diagonal only), so by Matsusaka–Mumford the geometric generic fibre of
`X → P²_Y` has `Aut = ⟨σ,τ⟩` (`w2c/sigma_*.out`); E1 forces `deg D ≥ 9` unless `D ∈ J(F)`; `Jhat` and all polars excluded
exactly.

**M-a (kept): cones with vertex `P²_X` over plane curves of base degree ≥ 9, or fibrations over `X`-ratios.** For
`D = g(X)`, `Z = X ∩ V(g)` is singular along `B = {K = 0}` (the Klein plane sextic, `H^{1,0}(B) = R(K)_3 ∋ x³, yz²` of
character 3) and fibred over `C_g` by the sextic surfaces `S_y = {K(x) + λ²M(x,y) + λ⁶G(y) = 0}` whose `H^{2,0}` contains
character 3 once (`xy`); E1 shows the Gysin image misses `x³t³` for base degree `≤ 8`. Verified: the geometry above; the
E1 tables through `d = 12` by character (search lens). Not verified: the Gysin image for base degree 9–12, and any
fibration over an `X`-ratio (where the `σ`-eigenvalues do mix) — the lens's part 2 (linear classification) was never written.
Decisive: (1) sector-restricted served-space test of `D = g(X)` for `g` in the character-3-compatible part of the degree-9..12
E1 kernels (`w2lib` sector solver, 1–5 k unknowns per sector, minutes each); (2) list the `X`-ratio fibrations and repeat the
fibrewise-automorphism computation (`Aut` of the generic fibre of `X → P²_X`-type maps) — an afternoon of msolve.
P(decide) 0.6, P(solve) 0.02. Lean: **L-A** only.
**M-b**: the §4.6 cover `u⁷ = (xt²/z³)` (adjoin a 7th root of the degree-0 character-3 function) with an exchange — speculative
(dossier §5.5); no exchange symmetry of the cover equation is in sight. P(solve) 0.01.

### 3.8 off-paired (`A = X0³Y0³`; det character of `D6 = ⟨h, σ2⟩` on the `(X2,Y2)`-block; sector (22,50,22), all `σ2`-anti-invariant)

Structure (verified exactly, `c1_identities.py`, Singular): `F = H(x) + C(x)X2Y2 + X2⁶ + Y2⁶`; `X/D6 = P(1,1,1,1,2)`;
Klein cover `Ŷ = X/⟨h²⟩: t1² = H + Cm − 2m³, t2² = H + Cm + 2m³` (`m = X2Y2`); `Y_∓ = X/⟨h², σ2 h^{0/3}⟩` are **Fano**
Weierstrass fibrations `t² = ±2m³ − Cm − H` over `P³_x` (`−K = O(3)`, discriminant `∝ 2C³ ∓ 27H²`) serving the other two
characters; the target is the `(−,−)` character, i.e. the `σ2`-anti-invariant part of `H⁴(X/⟨h⟩)`; the genus-2 family
`v² = (H + Cu)² − 4u⁶` depends only on `κ = C³/H²` and its total space `Σ ⊂ P(1,2,3,3)` is rational (2A1 + 2A2, explicit
parametrisation), so all `(3,1)`-classes involve `H^{2,0}` of the pencil surfaces `S̃_κ = {C = ks², H = ks³} ⊂ P(1,1,1,1,2)`
(`p_g = 45`); base changes `κ = ρⁿ` have `p_g(Σ̃_n) = 0,0,0,1,2,1` (n = 1..6). Excluded: all `X/K` (Reid–Tai, `P_2 ≥ 1`),
twisted products, §2.7 hubs (no-hub theorem: `Q_X(Y)` has rank 3), Galois closure (det served only by `p_g = 1` rotation
quotients), vertical denominators with `κ` constant or with linearly-equivalent conjugate components; served-space scan:
`X2⁶ − 2Y2⁶`, `C² − (X0²+Y1²)⁴`, `C² − (X0Y1+X1Y0)⁴`, `C² − X0⁸` serve only `J` in all 18 Y-odd sectors (two pass E1);
**no polynomial Mordell–Weil section of `Y_∓`** (`m` quadric, `t` cubic): msolve `[-1]` for `ε = −1` at 1000003, 2147483647,
1073741827 and for `ε = +1` at 1073741827 and (job `construct-off-paired-mw`, finished) at the two further primes.

**M-a (kept): descent through a reflection (new sufficient criterion, proved in part 3 §E).** If `p: Z → X` has a deck
automorphism `e`, an involution `σ` with `σeσ = d` and `p∘d = σ2∘p`, and `R` commutes with `σ` acting trivially on `p^*ξ`,
then divisor support of `H⁴` of a smooth model of `Z/⟨σ,R⟩` gives `ξ ∈ N¹` — no supported quotient of `X` is needed.
Realisation: `Z = X ×_W X^φ` with `W = X/σ2 = {H + Cm + s⁶ − 6s⁴m + 9s²m² − 2m³ = 0} ⊂ P(1,1,1,1,1,2)` (`s = X2+Y2`, Fano,
`−K = O(1)`), `φ` a birational involution of `W` **not lifting to `X`**, `e = (1,σ2)`, `d = (σ2,1)`, `σ = swap`, `R = ⟨(τ,τ)⟩`.
Verified: every automorphism `φ` (centraliser `⟨g², h³⟩`) makes `Z` reducible; a hypothetical linear involution of
`P(1,1,1,1,1,2)` preserving `W` gives `p_g(Z/⟨σ,r⟩) ≥ 1`; the natural birational candidates — fibrewise translations by
Mordell–Weil sections of `Y_∓` — have no polynomial instance (above). Not verified: sections with denominators
(`m = M/d², t = T/d³`, `deg d ≥ 1`); any non-fibre-preserving birational involution of `W` (a degree-6 index-1 weighted Fano
4-fold is expected to be birationally rigid, so `Bir(W)` should be generated by `Aut(W)` and untwisting involutions centred
on special points/curves — none identified); uniruledness of `Z/⟨σ,r⟩` for a given `φ` (fibrewise: a degree-4 cover of
`W/φ` modulo the swap must be Fano-fibred over a curve or `p_g = q = 0`-fibred).
Decisive: (1) MW sections with `deg d = 1` (`M` quartic 35, `T` sextic 84, `d` linear 4: 123 unknowns, 455 equations of
degree 12; msolve `-g 1` over `F_p`, hours) and `deg d = 2` (314 unknowns, 1330 equations; borderline) — P(section) ≈ 0.02;
(2) a rigidity/untwisting classification of `W` (theory; Pukhlikov-style for quasi-smooth weighted hypersurfaces) deciding
whether non-automorphic birational involutions exist at all. P(decide) 0.5, P(solve) 0.03.
Lean if positive: **L-B** — the holdout architecture verbatim (polynomial curve family on the Fano quotient `Z/⟨σ,r⟩`,
polynomial normal, integration, trace, transport; `PolynomialCurves` + `CertificateOperations`; a new cover-field layer
replaces `HoldoutCover`'s `F_H`-bound modules).
**M-b**: the type-V1 loophole — a Picard local system of the pencil `S̃_κ` with finite monodromy (compute `ρ(S̃_κ)` at a few
rational `κ` by point counts / van Luijk); if `ρ = 1` generically only mixed-type (M) divisors remain (fibrewise correspondences
coupling `J(E_κ)` to `H²_tr(S_κ)`). P(solve) 0.01. Pending: `construct-off-paired-served30` (degree-12 `H² − ν⁶`, `(X2³+Y2³)(X2³−Y2³)H`).

### 3.9 one-reflection-skew (`A = X0³Y0³`, `F = F_eq + (1/7)X0X1Y0Y1(Y0² − Y1²)`; even in `Y2` only; diagonal group of order 24; sector (0,44,77,3,0))

No construct lens was written for this lane; the material is fqlx/benediamond/DiamondJM/rbajaj5's thread work (dossier §3)
and the search lens's exact `6 → 5`-variable reduction (a certificate exists iff a rational 5-variable field `u = N/D` has
`A D² − D N(f) + (6/5) f (D div N − N(D)) = 0`, `f = F − Y2⁶`; then `(D², (N_iD), Y2(D div N − N(D))/5, 0)` is a strict
6-variable certificate of degree `2d`; validated on `family-hesse-web`). Excluded: the two hidden reflections of `one-reflection`
(the skew term is odd in `P = Y0+Y1` and in `M = Y0−Y1`); the §2.7 cover (the skew pair `(αX0X1/7)·yT·(α²T² − y²)` is not
`y ↔ T` symmetric); the `Y2`-reflection serves the wrong character; no cubic-surface reduction with one active `X`
coordinate (DiamondJM, Gröbner); `deg D ≥ 7` exact, degree-7 `Y2`-odd and three pure even characters excluded, character
`1110` + mixtures open; rigorous exclusions of several bounded primitive families (R1, R2 Lean-checked, R4 rank 26/27 and
38/39, C8 1,698-parameter family).

**M-a (kept): fqlx's twisted CM elliptic companion, completed by E4/E5.** Chart `Y0 = 1`, `X0 = x`, `X1 = rx`, `X2 = wx`,
`Y1 = t`, `Y2 = z`, `U = x²`: `F = z⁶ + P(U)`, `P = aU³ + bU² + cU + d`, `(a,b,c,d) = (1 + r⁶ + w⁶, t(r + r³w), rt(1 − t²)/7, 1 + t⁶)`
(the `U = x²` reduction is legitimate because the simultaneous `X`-sign change is a symmetry). Residue
`θ = (5U²/(72z¹¹)) dU∧ν`, `ν = dr∧dw∧dt`; Hermite reduction `θ = θ₀ + d(Sν)`, `θ₀ = (H0ω0 + H1ω1)/(72δ)∧ν`, `δ = disc_U(P)`.
Quartic cover `B = K[β]/(R)`, `R = 3a²β⁴ − 4abβ³ + 6acβ² − 12adβ + 4bd − c²` (`K = C(r,w,t)`), `M = b − 3aβ`,
`N = (c − 3aβ²)/2`, `k = N − Mβ`, `P = a(U+β)³ + (MU+N)²/M`; twisted j = 0 curves `E_β: υ² + M(1 + aξ³) = 0` with degree-3 maps
`φ_β: C_B → E_β` (`ξ = (U+β)/z²`, `υ = (MU+N)/z³`), `Tr_{B/K}(φ^*φ_*) = 6·I₂`, and `θ = d(Sν) + (1/6)Tr_{B/K}(φ_β^*Θ_β)`,
`Θ_β = C_β ω_E∧ν`, `C_β = M²(3aβ²k − 2N²)/(12ak⁴)`. Remaining equation on the total space `𝓛 = B(ξ,υ)`:
`∂_E W + Σ_i (D_iV_i − γ_iV_i) = C_β`, `γ = da/(3a) + dM/(2M)` (the twist), `W, V_i ∈ 𝓛` with moving poles allowed.
Verified (reported by fqlx, replayed mod p by DiamondJM): all identities; the shortcut via a translation field `T` fails
(`div T ≠ 0`); R3's rank-17 reformulation; R4's simpler cover (degree-18 `f(r) = 0`, independent of `t`) and the explicit
current `Ξ` reducing the source to a one-row term. Not verified: **E5** — whether the remaining base 2-form class vanishes on
the complement of some divisor (if a transcendental period is nonzero, no rational primitive exists in the compact
presentation; `W = 0`, base-only `V_i`); **E4** — the monodromy of the rank-2 unitary CM local system `L` (finite subgroup of
`U(2)`) and the Kodaira dimension / rational connectedness of the trivialising cover `B'' → (r,w,t)-space (Galois closure of
`R(β) = 0` with `√M` and `a^{1/3}`; branch locus `δ = 0` (degree ≤ 24), `a = 0`, `M = 0`). If a smooth model `B̄''` is
uniruled (or merely has `H³(B̄'')_{χ'} ⊂ N¹`), the class lies in `H³(B̄'') ⊗ H¹(E_ρ)` and is supported on `Z × E_ρ` for a divisor
`Z ⊂ B̄''` (covering-curve theorem on `B̄''`; the elliptic factor rides along as a closed (1,0)-form, no contraction needed),
and the quotient-transfer lemma pushes the support to `X`.
Decisive: E4(a)+(b) — explicit construction of `B''` (field extension of degree ≤ 144 over `C(r,w,t)`), its branch data and
`h^{3,0}`/`κ` of a resolution (weighted Hurwitz on a compactification; Singular/sympy, days); E5 — numerical periods of
`C_β ν` (or R4's reduced density) over 2-cycles of the base cover (one careful high-precision computation).
Honest expectation: a cover of `P³` branched along a degree-24 surface is typically of general type, so E4 most likely
returns `κ ≥ 0` (closing the route) — but the computation has never been done, and this is the only lane where a complete
transport identity to an auxiliary variety already exists. P(decide) 0.7, P(solve) 0.04.
Lean if positive: **L-B/L-C** — `PolynomialCurves` on `B̄''` (incidence seed over `F_p`, polynomial normal, integration,
trace to `K`), the twisted elliptic substitution through `ScalarEmbedding`-type lemmas, a new cover-field layer for
`R(β)`, and the exact source identities `θ = d(Sν) + (1/6)Tr(φ^*Θ_β)` (already in closed form).
**M-b**: E3 — first-order deformation of the accepted degree-30 half-certificate `(D₊, N₊)` of `one-reflection` along the
skew term (`≈ 1.4 M × 4.1 M` sparse modular system after `τ`-parity; engine-v3 `sge` on 64 cores could attempt it); gives a
first-order family only. P(solve) 0.02.

### 3.10 polarized-k3-remainder (`A = X0X1X2Y0Y1Y2`, `F = ½(f(U) + f(V))`; `A = A_K3 + A_mix`, both nonzero mod `J`)

Verified: `Aut(X) = D6 = ⟨h, τ⟩` exactly (corrected from the dossier; six planes `V = λU`, `λ⁶ = −1`); averaging-ideal
exhaustion (Theorem A5): no quotient of `X` by any subgroup of `D6`, and no intermediate quotient of `X → P² × P²`, serves
the K3 sector; `Aut(C_f) = 1`, `ρ(S) ≤ 4` unconditional (`≤ 3` with the thread's lift obstruction), `End_Hdg(T) = Q`
(thread); hence `A_K3` is equivalent to GHC (coniveau 1) for `Λ²T` of a very general degree-2 K3 — **new mathematics, not
a construction**, and the literature check (F8) confirms no applicable algebraic Kuga–Satake theorem.
**M-a (partial, kept as a deliverable): the six mixed monomials `A_mix` (sixth-root sector) via the second elliptic
companion.** Chart reduction for every `V`-degree `b` (sympy); explicit strict certificate for the joining sector
(`D = f(V)`, `N = (A/3)E_V`, verified with `core.verify_witness`); scalar conversion identity for weights (3,2) on
`W_{23,f} = {7f(x,a,y) + 7S² + 7T³ = 0} ⊂ P(1,1,1,3,2)` (smooth, `−K = O(2)`); **finite-field incidence seed found** (profile A,
13 unknowns: `x = 4t, a = 1 + t + 2t², y = 2 + 3t², S = 3 + 4t + 4t⁵ + t⁶, T = 1 + 2t² + t⁴` mod 5, `det J = 3`, freeness ranks
12/12, dominant), Hensel lift converges mod `5⁶`. Not verified/absent: Lean modules for the (3,2)-weight scalar calculus,
the `E_23` (`μ6`) contraction with its own theta kernels and infinity recurrence, the `b = 2,4` chart identities (≈ PR #149
size, months). Deliverable: a `HasCertificate` component for `A_mix` only. P(solve of the instance) 0.005; P(component
formalised in a long wave) ≈ 0.5. Lean: **L-C**.

### 3.11 quadratic-cubic-pullback (`A = X0³Y0³`, `F = Σ Q_i³`, `φ = (Q_i): P⁵ → P⁵` finite flat of degree 32, `X = φ⁻¹(Fermat cubic 4-fold C)`)

Verified: `f^*H^{3,1}(C) = ⟨δ⟩`, `δ = det(∂Q/∂Z)`; `A, δ` independent mod `J`; group `{1,τ}`; `F` not a join; the Galois
closure's intransitive quotients are of general type; all natural degree-6 supports serve only `J ⊕ ⟨δ⟩` (joins
`Q_a³ + Q_b³ + Q_c³` serve exactly `δ`, consistent with the adjoint-cokernel criterion); the two-plane Fermat mechanism is
unrealisable in any degree (`c_{20}F_{42} = −c_{02}F_{60}` impossible); `deg D ≥ 7` exact, `≥ 10` unless polar.
**M-a (kept): non-ACM adjoint schemes at degree ≥ 10 / an imprimitive intermediate cover.** Decisive: E7 — Frobenius cycle
types of `φ` over `F_p` (factor `(Q_i(z) − u_i)` for 50 random `u`; hours of setup, minutes of compute): if the monodromy of
`X/τ → C` is imprimitive there is an intermediate cover `C' → C` whose uniruledness would open §3.2 for the part of `V_15`
it carries (the mod-5 limit is the Kummer cover with group `(Z/2)⁵`, imprimitive; generically destroyed). Then E8b/E3 —
adjoint pre-screen of singular configurations at `d = 10` (curves on the 405 K3 surfaces `{Q_i = Q_j = Q_k = 0} ∩ X`,
`C_S ∪ C_{S'}` with 96 common points, corank-1 locus of the ramification surface) intersected with the degree-10 E1 kernel
(294/261), each survivor one sparse fixed-`D` solve (18 900 unknowns, 84 903 equations; ≈ 1 h with `sge`). P(decide) 0.8,
P(solve) 0.02. Lean: **L-A**.

### 3.12 quadric-contact-residual (`A = P0³ = (X0Y0)³` on `F_b = S0³ + 2P0²S0 + 5S0Q1 + (5/7)P0Q2 + X1⁶+Y1⁶+X2⁶+Y2⁶`; residue invariant under `G_aff` (48), odd under `s0: X0↔Y0` and `r1: X1↔Y1`; sector (0,22,51,22,0))

Verified: `deg D ≥ 7` (fqlx, replayed); `(J:A)_d = J_d` for `d ≤ 7`; unique non-Jacobian degree-8 direction `D8* = 4G_S² − G_P²`
(∝ `G_uG_v`, universal for `F = G((X0−Y0)², (X0+Y0)², …)`) inconsistent; degree-7 needs support on ≥ 3 of the 8 parity
classes; every `X/H`, `H ⊆ Γ⁺`, has `κ = 0` (Reid–Tai); no certificate keeps `(Y1/X1, Y2/X2)` passive (fqlx's Gauss–Manin
argument on the weighted-K3 fibration); the companion's `(X0,Y0) ↦ (S,T)` cubic-surface reduction has an `s0`-odd Jacobian,
so `P0³` is never a quotient density; no merging Fermat cover: `(X1,Y1)` is already a Fermat pair with `X/r1 ⊂ P(1,1,1,1,2,2)`
Fano of index 2 (serves only `r1`-odd forms), `(X2,Y2)` couples through `T², T⁴` only, and the `(X0,Y0)`-block
`S0³ + 2S0P0² = (X0²+Y0²)(X0⁴ + 4X0²Y0² + Y0⁴)` is not a sum of two sixth powers (angular gaps 35.3°, 35.3°, 109.4°).
The structure is the crossed-ternary pattern without the third mechanism: the residue is anti-invariant under `s0` and `r1`
and invariant under `s0r1`, whose fixed locus has codimension 2 (age 1), and no block-sign split supplies the `ρ = +1` sector.

**M-a (kept): a cyclic cover with a lifted exchange of nontrivial deck character (E5 of the dossier).** Adjoin
`T² = δ = X0² − Y0²` (the Jacobian factor of the companion's reduction) or `T⁶ = (X0−Y0)⁶ + …`; `s0` lifts to
`(X0,Y0,T) ↦ (Y0,X0,±iT)`; the question is whether `AΩ/F²` pulled back is odd under a lifted reflection of the cover whose
quotient is Fano or `p_g = q = 0`-fibred, with the correction carrying a nontrivial deck character (killed by the trace, as in
§2.7). Verified: nothing beyond the character table of `Γ`. Not verified: the averaging ideal of any proposed cover
(E8 bookkeeping is a finite character computation, minutes); uniruledness of the candidate quotients. Decisive: write the
cover `X̃ = {T² = δ} ⊂ P(1⁶,1)`, compute the lifted group (`Γ` extended by the deck `T ↦ −T`), its characters on `H⁴(X̃)`,
and test — for every subgroup `K` containing a reflection of `X̃` — whether `ρ_A^K ≠ 0` and `K_{X̃/K} < 0` (weighted adjunction).
P(decide) 0.9, P(positive) ≈ 0.03. Lean if positive: **L-B** (new cover fields; the companion `quadric-contact-broken`
supplies the cubic-surface seed machinery for the quotient integration).
**M-b**: supports transverse to the K3 fibration `X/⟨σ1,σ2⟩ → A²_{(t1,t2)}` — preimages of Noether–Lefschetz curves (threefolds
fibred in K3s on which the restricted class becomes algebraic); their `(3,1)` Gysin image depends on their singularities
(E4). P(solve) 0.01. Also pending: the structured `D8*·(form)` solves at degree 9–12 (E3).

### 3.13 rank-one-biquartic (`A = X0³Y0³`, `F = f(X) + g(Y) + (5/7)L(X)²M(Y)`; group `{1,τ}` only)

Verified: at `t = 0` the strict degree-6 witness `D = f`, `N = AE_X/6`; the first-order continuation is obstructed only by
E2 (`B_{D_X}(AX) = 8tAL²M` for `D_X = E_X(F)`), and the degree-6 family `dF_t(G_t)` with `G_0 = E_X/6` stays E1-admissible
(verify L5 — so "no degree-6 family" is **not** proved; what is proved is that no Y-even/odd/mixed degree-6 denominator works
at `t = 5/7`, `hessrelax6`); `X/τ → P² × P²` is a non-Galois triple cover (`S3` monodromy, branch `{g = 0} ∪ Δ'`,
`Δ' = 4t³L⁶M³ + 27f²g`); the 24 contracted ruled surfaces carry no 2-forms; `Z_f`, `Z_g`, `X ∩ V(L)`, `X ∩ V(M)` serve only `J`
(the predicted 3 classes of `C_f × Σ` die in the Gysin map); the Hessian symmetrizer is 1-dimensional (no linear splitting);
`deg D ≥ 7`, three Y-even degree-7 kernel elements (aggregator) being solved.
Excluded as mechanisms: all quotients and Galois-closure quotients (book §4.7 verbatim); the 3A elliptic contraction (the
`μ3` structure of the triple cover exists only at `t = 0`; for `t ≠ 0` the cover is non-Galois, so the companion
`{f + tL²M(Y) + S³ + T³}` of accepted-group8 §7 has no deck action); §2.7 covers (not analysed in a lens, but `F` has no
Fermat pair: every coordinate appears in `L` or `M`).

**M-a (kept): a correlated two-pencil support (fqlx's programme) — concretely the pull-backs of the degree-12 split loci.**
On `LM ≠ 0`, `X/τ` is the family of products of the pencil curves `C_X(u) = {f − uL³ = 0}`, `C_Y(k) = {g − kM³ = 0}` over the
incidence surface `u + ct_s + kt_s³ = 0` (`c = 5/7`); the fibre cubic `f s³ + tL²M s²u + g u³ = 0` splits over
`Γ ⊂ {g p³ + tL²M p q² + f q³ = 0}` with `p ∈ S(2,0)`, `q ∈ S(0,2)` (smallest such `Γ` has degree 12), and `π⁻¹(Γ)` is the
first reducible pull-back whose components can carry transcendental 2-forms; the `A5` chains of the ramification divisor
`R = X ∩ V(E_X F)` over `Σ = {w² = Mg}` give `P¹`-bundles over `Σ` (3 forms each). Not verified: E1/E2 of any degree-12
split-locus pull-back (fixed-`D` at `d = 12`: 25 704 + 16 797 unknowns, 162 248 equations — needs the sparse `sge` solver,
≈ hours), the (2,0)-budget of `R`'s resolution (predicted ≤ 72, Gysin unknown). Decisive: build `D_Γ = g p³ + tL²M p q² + f q³`
for random and structured `(p, q)` (bidegree (2,0)/(0,2) forms), test E1 (seconds) and, for survivors, the served space with
`served.py`/`fixed_d.target_image` at `d = 12`. P(decide) 0.6, P(solve) 0.02. Lean: **L-A**.
**M-b**: the Galois-closure companion `W = (Z × E)/μ3` with `Z` the double cover of `X/τ` branched along `Δ'` (a fivefold):
uniruledness and a `V_std ⊗ H¹(E)` contraction are both unknown; dimension 5 is beyond every accepted construction. P(solve) 0.01.

### 3.14 reflection-coupled-residual (`A = X0³Y0³`, Klein group `{1, r, ρ, τ}`, residue `(r,ρ) = (−,−)`; residual sector (126,492,126))

Verified: the residual sector is `H¹(P¹, j_*V₁₁)` for the `(−,−)` variation of the pencil `Y0 = λY1` (fibres: Klein-symmetric
sextic threefolds `X_λ`), `V₁₁` of rank 186 with Hodge numbers (3,90,90,3) (`p = 32003`); companion sector = `H³` of the Fano
fibres `V_λ ⊂ P(1,1,1,2,2)` (rank 82, level one); **no level-one sub-variation of `V₁₁` is seen by the nodal fibres**
(proposition: every sub-variation containing the three (3,0)-forms contains all nodal `(−,−)` vanishing classes), so
Abel–Jacobi/irregular-surface supports are impossible and a support must be a horizontal divisor carrying `(−,−)` 2-forms;
Galois closure of `X → P⁴_x` is `B3 × Z/2`, every quotient reaching `V ⊗ sgn_ρ` is non-uniruled (`X/τ`, a (4,6)-CI mod an
involution with ≥ 9 invariant canonical forms, a degree-8 hypersurface mod an involution); no exchange on any coordinate
pencil; only one polynomial Piola map (`Φ = (x, Y2²)`, the companion); `deg D ≥ 7`, degree 7 closed in three Klein characters,
`K8` excluded exactly.

**M-a (kept): M5 — twisted slices of the companion's rational-curve family.** `W = X/r ⊂ P(2,1⁵)` (Fano, `−K = O(1)`),
weighted conics (`e = 2`) on `W` form a 7-dimensional family; a 2-parameter slice `S` sweeps `E_W ⊂ W`, `E = π⁻¹(E_W) ⊂ X` is
fibred over `S` in genus-one hyperelliptic lifts, and `H^{2,0}(Ẽ)^{r=−1} ↪ H⁰(S, Ω¹_S ⊗ λ)` with `λ = f_*ω_{Ẽ/S}` a line bundle;
`e = 1` (lines) gives rank 0 (companion classes only). Not verified: a slice with `H⁰(S, Ω¹_S ⊗ λ) ≠ 0` **and** nonzero Gysin
image into the residual sector; the degree of `λ` on test curves of the conic family; `E = V₁₁` (`gm_rank.py`, queued in
`construct-rcr-all`). Decisive: compute `λ` (the Hodge bundle of the genus-one family over the conic parameter space: its
first Chern class from the discriminant of `u|_{C_s}`), choose positive slices, and evaluate the Gysin pairing of `[E]`'s
2-forms with the target by a Gauss–Manin/period computation on the pencil. P(decide) 0.5. This yields at best a
Hodge-theoretic existence statement at the generic point of the one-reflection locus (158 parameters), not a transport at
this `F`; termwise integration is obstructed by `H³(S, R¹p_*L) ≠ 0` for `e ≥ 2`. P(solve) 0.02. Lean: **L-D**.
**M-b**: M1 — a non-coordinate hidden/created fibre reflection (Exp B of `construct-rcr-all`, RUNNING: 0-dimensional msolve
system for reflections of `X` and of the generic fibres of the six coordinate pencils, with a 65-reflection Fermat control).
If found: §4.6 descent with `Z/⟨r,σ⟩` fibrewise Fano in `P(1,1,1,2,2)` → **L-B**. P(nonempty) ≈ 0.02.

### 3.15 short-power-presentation (`A = X0³Y0³`, `F = e(X) + e(Y) + c((L+M)⁶ + (L−M)⁶)`, `c = 5/7⁷`; `Aut = {1,τ}`)

Verified: shared-variable Thom–Sebastiani form `F = ft(L,X1,X2) + gt(M,Y1,Y2) + 30cL²M²(L² + M²)`; `X/τ` is the cubic cover
`(ε_X + 2c) + 30cw + 30cw² + (ε_Y + 2c)w³ = 0` of `P² × P²` with `S3` Galois group and a rational cuspidal branch curve; no
product structure (commuting-subgroup lemma: no ch. 3A contraction); no exchange on any subset-sum Fermat cover; lines move in
≤ 1-parameter families; **`X ≅ W ∩ {u = L²}` with `W = {e(X) + e(Y) + 2c(M⁶ + 15M⁴u + 15M²u² + u³) = 0} ⊂ P(1⁶,2)` a quasi-smooth
Fano fivefold (`−K = O(2)`)**, the reflection quotient of the smooth split join `Ĵ = {e(X) + g(Y,λ) = 0} ⊂ P⁶`; support transfer
from `W` fails because `K_W + X = 0` (`W∖X` has log Kodaira dimension 0, curves meet `X` twice: only the tautological identity
`mξ + Tξ ∈ N¹`); all E1-admissible denominators through degree 8 are polar residuals of the two (2+3)-split hyperplane sections
`X ∩ {L = 0}`, `X ∩ {M = 0}` and serve only `J`.
**M-a (kept): vertical divisors on the `(ε_X, ε_Y, t)` base** — pull-backs of curves on the base (ring `C[e(X), e(Y), L, M]`)
and of surfaces on the `(x,t)`/`(y,t)` bases (rings `C[X0,X1,X2,M]`, `C[L,Y0,Y1,Y2]`); distinguished member
`D_res = e(X)(L⁴ − L²M² + M⁴) + 2cL⁶(L⁴ + 14L²M² + M⁴)` (degree 10, the closure of `{ε_X = ε_Y}`, non-normal of multiplicity 4
along `X ∩ {L = M = 0}`). Not verified: E1 of any vertical family at `d = 7..10` (job `construct-spp-vert` never started;
resubmit `/tmp/hodge/results/short-power-presentation/w2c_pending_jobs.sh`, minutes). Decisive: that job, then full solves
of E1-survivors (`w2c_vert2.py`). P(decide) 0.9, P(solve) 0.02. Lean: **L-A**.
**M-b**: none (the Fano fivefolds `Φ⁶ ∩ H_λ ⊃ X` fail for the same `K + X = 0` reason).

### 3.16 symmetric-sextic-residual (`A = e3(X)e3(Y)` on `p6(X) + p6(Y) + (5/7)Q_XQ_Y(Q_X + 2Q_Y)`; `G = S3 × S3`; residue character sgn⊗sgn, τ = +; residue sector = polynomial (triv,triv) Y-odd column, Hodge numbers (17,52,17) — the dossier's (1,10,1) was the wrong column)

Verified: the residue of `A` and the 4-form `Ω/F` have the **same** `S3 × S3`-character (`det = sgn⊗sgn`), separated only by
`τ` — so F2 applies: every quotient of `X` or of a cover by a group acting the same way on `ω` and `α` has `κ ≥ 0`; all 35
residue-invariant non-identity elements have age ≥ 1 (ages table); the symmetric Fano `W ⊂ P(1,2,3,1,2,3)` serves only
multiples of `Δ_XΔ_Y`; `V = X/ker χ_res` is a CY-type double cover of `W`; `deg D ≥ 7` (any parity), no semi-invariant
degree-7 denominator; every legal `D` of degree 3–6 serves nothing beyond `J^G_6`; sign sectors reduce to the invariant one
(`D = Δ_X^aΔ_Y^b S ⇒ [A] ∈ served(D) ⇔ [A] ∈ served(S)`).
**M-a (kept): cones `D = Δ_X·s(X)` over `S3`-semi-invariant plane curves with vertex `{X = 0}`** (dossier 3.3(i)): the divisor
`E_D` fibres over `C = {Δ_X s = 0} ⊂ P²_X` with fibre the sextic surface `S_x ⊂ P³_{(s:Y)}` depending on `x` only through
`λ = p6(x)/e2(x)³`, pulled back from the single `S3(Y)`-symmetric pencil `S_λ = {λs⁶ + p6(Y) + (5/7)s⁴e2(Y) + (10/7)s²e2(Y)²}`
whose `H^{2,0}` contains `sgn_Y` three times (`s², e1Y², e2Y`); the `sgn_X` character comes from `Δ_X` and the twist of the
Hodge bundle along `C`. The unique degree-7 E1 survivors `D7 = Δ_X e2(X²)`, `D7Y` are of this type and fail. Not verified:
`H⁰(C, φ^*(Hodge bundle)^{sgn_Y} ⊗ twist)^{sgn_X}` for `deg s = 4..8` (E6: Griffiths residues `q(s,Y)Ω₃/F_s` and their
`λ`-dependence, splitting degrees on `P¹_λ`; one week). Decisive: E6, then E1 + the character-twisted solver (E2) on the first
curve with a nonzero count. P(decide) 0.6, P(solve) 0.03. Lean: **L-A**.
**M-b**: theory-permutation's twisted carrier `(X × Y)/K`: killed by O1 (ages ≥ 1) and the top-form test. P(solve) 0.

## 4. Cross-lane ranking of the five best routes

Criterion: P(solve) × fit to a Lean-provable argument (L-A explicit witness via the Kronecker checker, or L-B existential via
`PolynomialCurves` + a cover/quotient layer) × decidability of the named next computation. All five are long shots
(no route exceeds 4 %); the ordering is about where a wave of work produces either a certificate or a definite closure.

| rank | lane / route | why it ranks here | P(solve) | P(decide) | Lean |
|---|---|---|---|---|---|
| 1 | **one-reflection-skew — twisted CM companion, closed by E4 (geometry of the base cover `B''`) and E5 (period test)** | the only open lane with a *complete* transport identity to an auxiliary variety already in closed form (`θ = d(Sν) + (1/6)Tr_{B/K}(φ_β^*Θ_β)`, verified); the decisive geometric question (is a smooth model of the trivialising cover uniruled / is `H³(B̄'')_{χ'} ⊂ N¹`) has never been computed; if positive the route is the book's recipe on `B̄''` (seed, normal, integration, trace) and existence follows from §3.2 + quotient transfer | 0.04 | 0.7 | L-B / L-C |
| 2 | **four-plus-two-residual — core-side support via the exact 5-variable reduction** | the cheapest decisive tests in the programme (reduced solver, seconds per `D`), an explicit design principle (adjoint budget in the `ζ6²` character of `μ6`-covers of singular core surfaces), and a payoff that is an explicit strict witness (L-A, no new library) | 0.04 | 0.7 | L-A |
| 3 | **cycle-q2 — M\*: non-fibrewise rational primitive, decided by the twisted Gauss–Manin computation** | the remaining equation is precisely stated (`[σ] = 0` in `H²_dR(B°, V_{ω²})`, rank 12) and every other mechanism is excluded by exact computation; a rational current gives an explicit `q = 2` certificate through fqlx's adapter (L-A); the cost is one two-variable D-module computation, outside the library but well within Singular/sympy | 0.04 | 0.6 | L-A |
| 4 | **off-paired — descent through a reflection via a birational involution `φ` of the Fano `W = X/σ2`** | a proved new sufficient criterion and the richest uniruled geometry on the same block family (`Y_∓` Fano Weierstrass fibrations, `W` Fano of index 1); first steps are cheap (MW sections with denominators); but the required `φ` probably does not exist (expected birational rigidity of `W`; no polynomial MW section at 3 primes) | 0.03 | 0.5 | L-B |
| 5 | **polarized-k3-remainder — Part B, the `A_mix` component on `W_{23,f}`** | the only route with a verified finite-field incidence seed and a Hensel lift; it delivers a reusable `HasCertificate` component (sixth-root sector of polarized-split sextics) but cannot close the instance, whose K3 part is an open GHC statement (F8) | 0.005 (instance) | — | L-C (partial) |

Routes that would jump the table if an in-flight check comes back positive (each ≈ 2 %): diagonal-permutation J2/J2f
(hidden/created second reflection on the `(01)`-pencil fibres → §4.6 descent, L-B; `construct-dp-j1j2` running,
`construct-dp-fibre` to resubmit), reflection-coupled M1 (non-coordinate fibre reflection, `construct-rcr-all` running),
quadric-contact E5/E8 (lifted exchange on `T² = X0² − Y0²`; minutes of character bookkeeping, never done), cyclic-triple
E5 (the Fano variety of lines `𝓕` of `W`; hours of Gröbner; payoff L-A only), short-power vertical families (job never ran).

What would change the ranking: (a) a positive E4 for the skew lane makes it the first open lane with an existence
theorem and a clear Lean path — it moves to "work it to a PR"; (b) a negative E4 (κ(B̄'') ≥ 0) plus a nonzero period in E5
closes the companion route and demotes the lane below cycle-q2; (c) any positive hidden-reflection search (J2f, rcr-all)
immediately outranks everything (the §4.6 architecture is formalized and the seed machinery exists); (d) `Frobenius` data
showing `End(H³(T_H)^−_χ) ⊋ Q(ζ3)` for four-plus-two would raise its P(solve) to ≈ 0.1 by naming the Prym correspondence.

Lanes with no construction route worth a work order: dense-rational, eight-plane-configuration (both: every universal
object has zero (2,0)-budget; only a non-linear created reflection could help and none is suggested by the geometry),
klein-coupled-character and symmetric-sextic-residual (cone-type supports are the only character-compatible geometry; both
need the Hodge-bundle computation of §3.7/§3.16 before any solve), rank-one-biquartic, quadratic-cubic-pullback,
short-power-presentation (search-type candidates only), reflection-coupled-residual (Hodge-theoretic existence at the
generic point of the one-reflection locus at best), cyclic-triple-residual (elliptic, not rational, geometry: L-A only).

Skeptic's note on the book recipe. Of the sixteen lanes, a finite-field incidence seed can be *written down* today only for
polarized-k3's companion `W_{23,f}` (verified) and, conditionally, for off-paired's `Z/⟨σ,r⟩` (if `φ` exists) and the skew
lane's `B̄''` (if uniruled). Every other lane fails the recipe at **dominance** (the auxiliary surface/fourfold carrying the
character has `κ ≥ 0`: cycle-q2's `S3`, cyclic-triple's elliptic sweeps, reflection-coupled's genus-one lifts, all honest
quotients by F2/F3), not at seed-finding or lifting. That is why the two L-A routes (four-plus-two, cycle-q2) are ranked
alongside the recipe-shaped ones: an explicit witness needs no auxiliary uniruled variety at all.
