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
