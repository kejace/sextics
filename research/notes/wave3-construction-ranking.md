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
| one-reflection-skew | fqlx/benediamond twisted CM companion; E4 now negative (base cover of general type, §3.9 addendum), E5 diagnostic | quartic cover `R(β) = 0`, `E_β`, base cover `B''` | period test E5; R3 rank-17 formulation | 0.7 | 0.02 | L-B only if a smaller uniruled carrier is found |
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

**Addendum (this synthesis, sympy-verified, `/tmp/hodge/results/wave3-construct/skew_branch.py`).** E4's optimistic branch
is closed by a Hurwitz count: `disc_β(R) = −6912·a⁴·δ²` with `δ = disc_U(P)` (so the quartic cover `B_4 → (r,w,t)`-space is
branched exactly along `{δ = 0}`, total degree 24, and `{a = 0}`, degree 6), and `N_{B/K}(M) = −9a²(3ac − b²)²` (so the `√M`
twist is branched along `{3ac − b² = 0}`, degree 10, and `{a = 0}`); the `a^{1/3}` twist is branched along `{a = 0}`. For any
Galois cover `B'' → P³` containing these, `K_{B''} = π^*(O(−4) + Σ(1 − 1/e_i)B_i) ≥ π^*O(−4 + 12 + 5 + 4) = π^*O(17)`, and the
branch surfaces have only discriminant-type (canonical) singularities, so every smooth model of `B''` is of general type —
**not uniruled**; the covering-curve theorem cannot be applied to `B̄''`, and `H³(B̄'')_{χ'}` is at best a GHC statement.
Consequently the skew route is a diagnostic (E5: a nonzero transcendental period of `C_β ν` on the base cover proves that no
rational primitive exists in the compact/companion presentation and explains the failed ansätze), not a construction, unless a
*smaller* uniruled carrier of the `χ'`-class is found (none is in sight). Revised: P(solve) 0.02, Lean only via a new carrier.

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
`PolynomialCurves` + a cover/quotient layer) × decidability of the named next computation. All five are long shots (no route
exceeds 4 %); the ordering is about where one wave of work produces either a certificate or a definite closure. The
one-reflection-skew route, which looked strongest on paper, was demoted during this synthesis by the Hurwitz count of
§3.9 (its auxiliary base cover is of general type).

| rank | lane / route | why it ranks here | P(solve) | P(decide) | Lean |
|---|---|---|---|---|---|
| 1 | **four-plus-two-residual — core-side support via the exact 5-variable reduction** | the cheapest decisive tests in the programme (reduced solver, seconds per `D`, validated on two accepted core-only witnesses); an explicit design principle (adjoint/irregularity budget of `μ6`-covers of singular core surfaces in the `ζ6²` character, Zariski/Libgober-type superabundance); the sector is level 2 of abelian type (GHC predicts support, no `(4,0)`); payoff = explicit strict witness, no new library | 0.04 | 0.7 | L-A |
| 2 | **cycle-q2 — M\*: non-fibrewise rational primitive, decided by the twisted Gauss–Manin computation** | the remaining equation is precisely stated (`[σ] = 0` in `H²_dR(B°, V_{ω²})`, rank 12 over `C(c1,c2)`), every other mechanism is excluded exactly, the frames (DiamondJM's second variation, the lens's `w2c_s3.py` Weierstrass data) exist; a rational current gives an explicit `q = 2` certificate through fqlx's pole-three adapter | 0.04 | 0.6 | L-A |
| 3 | **off-paired — descent through a reflection via a birational involution `φ` of the Fano `W = X/σ2`** | the only route with a proved *new* sufficient criterion that needs no supported quotient of `X`, on the block family where the holdout cover works; the auxiliary geometry is uniruled (`W` Fano of index 1, `Y_∓` Fano elliptic fibrations); first steps are cheap; the required `φ` probably does not exist (expected birational rigidity; no polynomial MW section at three primes), which is why it is third despite being the best recipe-shaped route | 0.03 | 0.5 | L-B |
| 4 | **one-reflection-skew — fqlx's twisted CM companion, now as a diagnostic (E5), plus R3's rank-17 formulation** | the only lane with a complete transport identity to an auxiliary variety, but the base cover trivialising the CM local system is of general type (§3.9 addendum), so the covering-curve theorem does not apply; E5 (period test) would close or revive fqlx's compact primitive problem; a positive answer would still need a smaller uniruled carrier for Lean | 0.02 | 0.7 | L-B only with a new carrier |
| 5 | **polarized-k3-remainder — Part B, the `A_mix` component on `W_{23,f}`** | the only verified finite-field incidence seed in the programme; delivers a reusable `HasCertificate` component (sixth-root sector of polarized-split sextics) but cannot close the instance (K3 part = open GHC, F8) | 0.005 (instance) | — | L-C (partial) |

Routes that would jump the table if an in-flight check comes back positive (each ≈ 2 %): diagonal-permutation J2/J2f
(hidden or created second reflection on the `(01)`-pencil fibres → §4.6 descent, L-B; `construct-dp-j1j2` running,
`construct-dp-fibre` to resubmit), reflection-coupled M1 (non-coordinate fibre reflection, `construct-rcr-all` running),
cyclic-triple E5 (the Fano variety of lines `𝓕` of `W`; hours of Gröbner; payoff L-A only), short-power vertical families
(job never ran, minutes), quadric-contact E8 bookkeeping for lifted exchanges (minutes; expected negative because any
automorphism of a cyclic cover `T^n = g(Z)` normalises the deck group and descends to `X` — a non-descending exchange needs
a coordinate that mixes with `T`, as `y ↔ T` does in holdout, and `F_b` offers none).

What would change the ranking: (a) any positive hidden-reflection search (J2f, rcr-all) immediately outranks everything (the
§4.6 architecture is formalized and the seed machinery exists); (b) Frobenius data showing `End(H³(T_H)^−_χ) ⊋ Q(ζ3)` for
four-plus-two would raise its P(solve) to ≈ 0.1 by naming the Prym correspondence; (c) a rational current in cycle-q2's
M\* turns that lane into a pure Kronecker-checking exercise; (d) an explicit non-automorphic birational involution of
off-paired's `W` with uniruled `Z/⟨σ,r⟩` turns off-paired into a holdout-style PR.

Lanes with no construction route worth a work order: dense-rational, eight-plane-configuration (every universal object has
zero (2,0)-budget; only a non-linear created reflection could help and the geometry suggests none), klein-coupled-character
and symmetric-sextic-residual (cone-type supports are the only character-compatible geometry; both need the Hodge-bundle
computations of §3.7/§3.16 before any solve), rank-one-biquartic, quadratic-cubic-pullback, short-power-presentation
(search-type candidates only), reflection-coupled-residual (Hodge-theoretic existence at the generic point of the
one-reflection locus at best), cyclic-triple-residual (elliptic, not rational, geometry: L-A only).

Skeptic's note on the book recipe. Of the sixteen lanes, a finite-field incidence seed can be *written down* today only for
polarized-k3's companion `W_{23,f}` (verified) and, conditionally, for off-paired's `Z/⟨σ,r⟩` (if `φ` exists). Every other lane
fails the recipe at **dominance** (the auxiliary variety carrying the character has `κ ≥ 0`: cycle-q2's `S3`, cyclic-triple's
elliptic sweeps, reflection-coupled's genus-one lifts, the skew lane's base cover, all honest quotients by F2/F3), not at
seed-finding or lifting. That is why the two L-A routes rank first: an explicit witness needs no auxiliary uniruled variety.

## 5. Work orders for the top three routes

Each work order is self-contained: a fresh agent needs only this section, the engine, and the named files. Common rules:
log every step to `<lane dir>/wave3-construct.log.jsonl` (append, one JSON object per step, immediately); rewrite the lane
report `/tmp/hodge/notes/wave3-<lane>-construct.md` after every result; heavy jobs go through
`cd /home/kejace/devel/sextics-work && python3 -m hodge.rjob submit --lane <lane> --name <job> --cpus N --wait 3600 -- <cmd>`
(outputs to files named after the job inside the lane dir; `python3 -m hodge.rjob pull --lane <lane>` to collect); never
`rjob up` a lane whose remote jobs are appending to local-named files; write squares as `^` in msolve inputs; use primes
`< 2^28` or `998244353` for msolve; a modular "consistent" is only a witness after the exact rational re-solve
(`fixed_d.solve_fixed_D_fast(..., want_exact=True)` + `core.verify_witness`). Never post; never run git.

### WO-1 — four-plus-two-residual: designed singular core surfaces tested in the exact 5-variable model

**Goal.** Find a core-only denominator `D(X0,X1,Y0,Y1)` of degree `k ≤ 14` such that the reduced strict system is consistent
(then a strict 6-variable witness follows by the proved lift), or establish that every core surface whose `μ6`-cover has
positive `ζ6²`-irregularity budget up to degree 14 fails.

**Inputs.** `data/four-plus-two-residual/instance.rec` (`F = H + B`, `H = f(X0,X1) + g(Y0,Y1) + K`, `f = X0⁶ + (5/7)X0⁵X1 +
(10/7)X0X1⁵ + X1⁶`, `g = Y0⁶ + (5/7)Y0⁵Y1 + (15/7)Y0Y1⁵ + Y1⁶`, `K = (5/7)X0³X1Y0Y1 + (10/7)X0X1Y0³Y1`, `B = X2⁶ + (15/7)X2⁴Y2² +
(10/7)X2²Y2⁴ + Y2⁶`; `A = X0³Y0³`). Engine `/home/kejace/devel/sextics-work/hodge` (`core.load_instance`, `fixed_d.e1_solvable`,
`fixed_d.target_image`, `fixed_d.solve_fixed_D_fast`, `core.verify_witness`). Lane dir `/tmp/hodge/results/four-plus-two-residual/`:
`threefold.py` (`reduced_test(Hcore, a, core_idx, D, p)` — fixed-`D` system of `(H + U⁶, aU, D)` in the `U`-degree sector,
returns rank/augmented rank; `e1_kernel_core(k, par, p)`; `lift_basis(k, par, nprimes)`; `run3(D, note, family)` logs to
`log.jsonl`), `corekernel.json`, `wave2-excl/excl-kernels-corereps.json` (exact rational core-kernel representatives,
`k ≤ 12`), `h.py`, `threefold_cones.py` (template for cone families), `prym.py` (point-count template). Reference for the
geometry: attack1 §1.5/§4.2 (`/tmp/hodge/notes/attack1-four-plus-two-residual.md`), cyclic-triple construct part 4 (the
`x³ = c3(n1,n2)` Ẽ6 mechanism), quadratic-cubic attack1 §1.3 (adjoint-condition criterion).

**Mathematical frame (what the scripts implement).** `Y_H: U⁶ = −H ⊂ P⁴` (sextic threefold, `μ6` acting on `U`),
`S_H = {H = 0} ⊂ P³` smooth. For a core `D` of degree `k`, `S_D = V(D) ∩ Y_H` is the `μ6`-cover of the surface `V = V(D) ⊂ P³`
branched along `Γ = V ∩ S_H ∈ |O_V(6)|`. The target's core class `[aU] ∈ H^{2,1}(Y_H)` lies in a fixed `μ6`-eigenspace; it is
supported on `S_D` only if the resolution `S̃_D` has `H^{1,0}` in the matching character (`ζ6²` in the convention of
attack1 §1.1: `aU ↔ (2,4)` sectors; verify the convention against the reduced solver on `family-hesse-web`). For smooth `V`
and smooth `Γ`, `H¹(S_D, O)_{χ^i} = H¹(V, O_V(−i)) = 0` (Lefschetz), so **irregularity can only come from singularities**:
(S1) singular points of `Γ` created by tangency of `V` with `S_H` or by singular points of `V` on `S_H` — by Esnault–Viehweg
/ Libgober, the `χ^i`-irregularity of the resolution equals the superabundance `h¹(V, O_V(−i) ⊗ J((i/6)Γ))` (`J` = multiplier
ideal); for the character paired with `ζ6²` the contributing singularities are those with log-canonical threshold `≤ 2/3`:
ordinary triple points of `Γ` (`lct = 2/3`, `J = m_p`), i.e. **ordinary triple points of `V` lying on `S_H`** (each imposes 10
linear conditions on `D`), `D4`-type and worse; cusps (`lct 5/6`) feed the `ζ6^{±1}` characters and are useless here;
(S2) non-normal `V` along a curve `C ⊂ S_H` (`D ∈ I_C²`), whose normalisation data can carry `H⁰(C̃, …)`-type classes;
(S3) curves of non-rational transversal type (`D ∈ I_C³` for a curve `C ⊄ S_H` meeting `S_H` in `6·deg C` points, each an
Ẽ6-type point of `S_D` with a j = 0 elliptic exceptional curve of character `ζ3`).
Before any search, the required budget must be known.

**Step 0 (sector-variation bound, local, minutes).** In the 4-variable Jacobian ring `R(H)` (Hilbert function
`(1,4,10,20,35,52,68,80,85,…)`), with `Def = R(H)_6^{ι1-even}` (the `(Y0,Y1)`-even sextics in 4 variables modulo `J(H)`) and
`ξ = [a]`, compute `dim Def·ξ ⊂ R(H)_{12}^{odd}`, `dim Def²·ξ ⊂ R(H)_{18}^{odd}`, `K = {v ∈ R(H)_6^{odd} : Def·v ⊂ J}`,
`K2 = {u ∈ R(H)_{12}^{odd} : Def·u ⊂ J}` mod two primes (copy `c5_sector_vhs.py` / `c5b_K2.py` from
`/tmp/hodge/results/cyclic-triple-residual/wave2-construct/`, replacing the ring and parity filter). If `K = K2 = 0` the
smallest sub-VHS containing `[a]` over the `ι1`-symmetric family is the whole `(32,80,32)` sector, and every support of the
core class on a very general member needs `h^{1,0}(S̃_D)_χ ≥ 32`: then (S1) needs `≥ 32` superabundant triple points
(`k ≥ 12` by the count `10m ≤ dim S_k − 1`), and (S3) needs curves with `Σ genus ≥ 32`. Record the numbers; they set the
degree range below. (If `K` or `K2` is nonzero the bound is smaller and low degrees come back into play — report it.)

**Step 1 (budget calculator, write `w3_budget.py`, local, seconds per candidate).** Input: a candidate configuration
(list of triple points `p_j ∈ S_H(F_p)` with chosen tangent data, or a curve `C` with multiplicity). Output: (a) the linear
space `L_k = {D ∈ S_k : D has ordinary triple points at the p_j}` (or `D ∈ I_C^m`), its dimension, and its intersection with the
exact core E1 kernel (`A·D ∈ J(H)` — use `e1_kernel_core` or the stored representatives; the intersection must leave `J(H)_k`);
(b) the superabundance `h¹(V, O_V(−i) ⊗ ∏ m_{p_j})` for the relevant `i` computed as
`dim H⁰(V, O_V(d_i) ⊗ ∏m_{p_j}) − (dim H⁰(V, O_V(d_i)) − #conditions)` with `d_i` the Serre-dual twist (`K_V = O_V(k − 4)`;
write out the exact twist for `i = 2` and `i = 4` and check both against the known Fermat join control in attack1 §1.3);
the count is a rank of an evaluation matrix of monomials at the points mod `p`. Keep only configurations with positive budget
in the target character.

**Step 2 (enumeration + reduced solve, remote).** Write `w3_design_sweep.py`: for `k = 9..14` and `m = 1..⌊(dim S_k − 1)/10⌋`
triple points sampled among `F_p`-rational points of `S_H` (random and structured: points on the curves `C_f = {X0 = X1 = 0}`,
`C_g`, on the lines `X_i = 0`, on the hyperplane sections through the coordinate points), build `L_k ∩ (core E1 kernel)`,
take random members and all basis members, and run `reduced_test` (modular, `p = 2147483647` then a second prime);
for (S3) use `D ∈ I_C³` for lines and conics `C ⊂ P³` (twisted cubics at `k ≥ 12`); for (S2) `D ∈ I_C²` with `C = S_H ∩ V(G)`,
`deg G ≤ 3`. Submit as `python3 -m hodge.rjob submit --lane four-plus-two-residual --name w3-design-k<k> --cpus 8 --wait 3600 --
python3 w3_design_sweep.py --k <k> --out w3-design-k<k>.jsonl` (one job per degree; each reduced solve is 3–10 k unknowns,
seconds to minutes; budget ≈ 2–6 h per degree). Every "consistent" hit: re-solve exactly (CRT lift of `D`, then
`solve_fixed_D_fast(want_exact=True)` on the 6-variable system with the lifted strict field) and `core.verify_witness`.

**Step 3 (theory probe, remote, optional but cheap).** Point counts of `S_H` over `F_{p^n}` (`n ≤ 4`, `p ≤ 200`) and of the
genus-4 curves `C_f`, `C_g` to get Frobenius on `H³(Y_H)`/`H²(S_H)` sectors (template `prym.py`); test whether the core
sector's Frobenius polynomials factor more than `Q(ζ3)` forces. Extra endomorphisms would point to an explicit Prym
correspondence (name the surfaces for Step 2).

**Success criterion.** A rational core `D` with a strict reduced certificate verified exactly, lifted to a strict
6-variable witness (`Q = 0`) passing `core.verify_witness` and `IsRelPrime` by specialization (`D` omits `X2, Y2`, in which
`F` is monic) → `Witness.lean` via the Kronecker route (`lean-capacity.md`), `Solution.lean` in the accepted template.
Negative result worth posting (thread #120): the budget bound from Step 0 plus "every positive-budget configuration of
degree ≤ 14 is inconsistent" (modular, two primes), with the exact E1-kernel representatives.

### WO-2 — cycle-q2: decide `[σ] = 0` in the rank-12 Kummer-twisted cohomology, and build the q = 2 certificate if it holds

**Goal.** Compute the Gauss–Manin connection of the rank-12 local system `V = V_{ω²}` over the `(c1,c2)`-plane carrying the
target's `(3,1)`-part, express DiamondJM's second-variation source `σ` in it, and decide whether `σ` is `∇`-exact with
rational coefficients (poles along the discriminant only, or along a chosen horizontal divisor). A rational primitive `Ψ`
yields an explicit `q = 2` certificate through fqlx's pole-three clearing adapter.

**Inputs.** `data/cycle-q2/instance.rec` (`F = Σ(X_i⁶ + Y_i⁶) + ½(X0³X1Y0Y1 + X1³X2Y1Y2 + X0X2³Y0Y2)`, `q = 2`,
`A = X0²X1²X2³Y0²Y1²Y2`). Lane dir `/tmp/hodge/results/cycle-q2/`: `w2c_s3.py` (the Weierstrass data `s4, s6` of the cubic `U`
over `(t; c1, c2)`), `w2c_qraise.py` (the verified q-raising identity, useful as a template for certificate-level identities),
`w2c_factor.json`. Construct report parts 1–5 (`/tmp/hodge/notes/wave2-cycle-q2-construct-part*.md`), dossier
`/tmp/hodge/notes/dossier-cycle-q2.md` (§5 P2/P5, A.7: the adapter and the pole-three clearing), DiamondJM's repository mirror
`/tmp/hodge/ext/hodgeHoldout` (second-variation frames; the facts "P ∉ J37, P ∈ J37 + ∂_{c1}G·S19 + ∂_{c2}G·S19"), the accepted
`data/cycle-three/` proof (`PhysicalPrimitive.lean`, `solution.md`: the RES transport at `q = 1`, moving-map and normal-change
corrections, and the Lean modules `ResidueTransport`, `TransitiveCubic` that already encode the cubic's covariants).

**Objects.** Paired-sign chart: `u = (X0², X1², X2²)`, `v = Y0/X0`, `c_i = (Y_i/X_i)/v`, `t = v²`; cubic
`U = a x³ + b y³ + c z³ + d x²y + e y²z + f z²x`, `(a,b,c) = (1 + t³, 1 + c1⁶t³, 1 + c2⁶t³)`, `(d,e,f) = (t/2)(c1, c1c2, c2)`;
Jacobian `η² = ξ³ + s4 ξ − s6/4` (`s4, s6 ∈ Q[c1,c2][t³]`, degrees 6 and 18 in `t`); `S3 = S3(c1,c2)` the elliptic surface over
`P¹_t` (χ = 3, `p_g = 2`), `RES` the same over `P¹_T`, `T = t³`; `μ3` acts on `S3` by `t ↦ ωt`; `V := H²_tr(S3)_{(sector of ω₀ = dt∧dξ/η)}
≅ H¹(P¹_T, R¹_RES ⊗ L_{ω²})`, rank 12, Hodge type (1,10,1), a local system over `B° = C²_{c1,c2} ∖ Δ` (Δ = locus where the
Weierstrass fibration degenerates further: vanishing of the discriminant's leading/constant coefficients, `c1c2 = 0`, …).

**Step 1 (frame and connection, write `w3_gm.py` with sympy/flint; remote, hours).** Use the Griffiths–Dwork method on the
Weierstrass surface: with `G(t, ξ, η) = η² − ξ³ − s4 ξ + s6/4` (weighted-homogenised in a `P(1,1,2,3)`-bundle chart over `P¹_t`
or simply affine in `(t, ξ, η)` with pole-order filtration in `t` at `∞` and at the singular fibres), the classes
`t^j ξ^i dt∧dξ/η` (`i ∈ {0,1}`, `j` in the residue class fixed by the `ω`-character under `t ↦ ωt`) span `V` after reduction;
pick the 12-dimensional basis by Gaussian elimination mod `p` at a random rational `(c1,c2)` first (to learn the degrees),
then over `Q(c1,c2)` by evaluation at ≈ 100 points mod `p` and rational reconstruction (the connection matrices
`Ω_1 = ∇_{c1}`, `Ω_2 = ∇_{c2}` are `12 × 12` with entries in `Q(c1,c2)`; expect denominators = factors of `Δ`). Reduction rule:
`∇_{c_i}[P dt∧dξ/η] = [(∂_{c_i}P) dt∧dξ/η − P·(∂_{c_i}G/(2η²))·… ]` — i.e. differentiate under the integral and reduce the
pole order of `η` with the Weierstrass relation and the `t`-derivative (standard Picard–Fuchs reduction for elliptic surfaces;
Schütt–Shioda §§6–8 fix the minimal model and `χ`). Check: integrability `∂_{c2}Ω_1 − ∂_{c1}Ω_2 + [Ω_1, Ω_2] = 0`, and the local
monodromy around `Δ` has the expected finite/unipotent form.

**Step 2 (the source `σ`).** Transport the second-variation source `2AΩ/F³ = ∂_{s1}∂_{s2}(Ω/F_s)|_0`, `F_s = F + s1·X0X1X2Y0Y1Y2 +
s2·X0X1X2²Y0Y1`, through the paired-sign cubic exactly as `cycle-three` transports its `q = 1` source (`source −V_α/(144r)`,
multiplier `−144 r v³`, moving-map and normal-change corrections, now applied twice — the dossier's A.7 and DiamondJM's
second-variation identity give the formulas); keep the factor `t = v²`; reduce the result to the frame of Step 1: `σ = Σ_{k}
σ_k(c1,c2) e_k ⊗ dc1∧dc2` with `σ_k ∈ Q(c1,c2)`. Check H4 (DiamondJM): the fibrewise `(0,2)`-component of `σ` is absorbed, i.e.
`σ` has no component along the `(0,2)`-direction after adding `∇`-exact terms — a consistency test of the frame.

**Step 3 (decision).** Solve `∇(Ψ_1 dc1 + Ψ_2 dc2) = σ dc1∧dc2`, i.e. `∇_{c1}Ψ_2 − ∇_{c2}Ψ_1 = σ`, for `Ψ_i ∈ V ⊗ Q(c1,c2)` with
denominators `Δ^N` and numerator degrees `≤ M`, increasing `(N, M)` (linear algebra over `Q`, sizes `12·C(M+2,2)` unknowns;
mod `p` first). If a solution appears: lift to `Q` (CRT), verify exactly, go to Step 4. If none up to `(N, M) = (8, 40)`:
compute `dim H²_dR(B°, V)` by the pole-order filtration (the complex `V → V ⊗ Ω¹ → V ⊗ Ω²` with polynomial coefficients and
poles along `Δ`; the reduction algorithm bounds the pole order needed to represent every class — regular singularities —
so exactness is decidable); if `[σ] ≠ 0`, test the candidate horizontal divisors of the dossier's P5 (Noether–Lefschetz
loci of the `ω₀`-sector: curves `Γ_c ⊂ S3` found by the lens's sector solver at low degree) by allowing poles along their
images. Record the dimension of `H²_dR(B°, V)` and the class of `σ` in it.

**Step 4 (certificate).** With rational `Ψ`, the dossier's P2/A.7 adapter (fqlx) produces the `q = 2` field: pull `Ψ` back
through the paired-sign chart (`u_i = X_i²`, `v`, `c_i`), apply the two moving-map corrections, clear poles along the polar
locus of `Ψ` (horizontal, non-normal — this is the expected support), and obtain `(D, N, Q)` with `D = (polar locus)·(chart
divisors)`; verify with `core.verify_witness` (`q = 2`, `deg N_i = deg D + 7`); legality by specialization. Expect `deg D`
in the tens; the Lean witness goes through the Kronecker checker (`lean-capacity.md`: slice if products exceed 3 M).

**Remote jobs.** `w3-gm-frame` (Step 1, `--cpus 8`, sympy/flint evaluation-interpolation, 1–4 h), `w3-gm-source` (Step 2,
`--cpus 4`), `w3-gm-solve` (Step 3, `--cpus 16`, the linear systems; hours), all under `--lane cycle-q2`, outputs
`w3-gm-*.json` in the lane dir.
**Success criterion.** A verified rational `Ψ` and a `q = 2` witness passing `verify_witness`; or the rigorous statement
`[σ] ≠ 0 ∈ H²_dR(B°, V)` with the computed dimension, posted to thread #61 as the first exact obstruction at the base level.

### WO-3 — off-paired: decide whether the descent-through-reflection criterion can be realised (`φ`), and if so seed it

**Goal.** Either exhibit a birational involution `φ` of the Fano fourfold `W = X/σ2` not lifting to `X` such that
`Z = X ×_W X^φ` modulo `⟨swap, (τ,τ)⟩` has a Fano fibration or `p_g = q = 0` surface fibration (then the holdout recipe gives
existence and a seed), or prove that no such `φ` exists (closing the only new criterion of the lane).

**Inputs.** `data/off-paired/instance.rec` (`F = H(x) + C(x)X2Y2 + X2⁶ + Y2⁶`, `x = (X0,X1,Y0,Y1)`, `H = X0⁶ + X1⁶ + Y0⁶ + Y1⁶ +
X0³X1Y0Y1/2`, `C = X1³Y1/2 + X0²X1Y0/3`; `A = X0³Y0³`). Construct report parts 1–4 (`/tmp/hodge/notes/wave2-off-paired-construct-part*.md`;
the lemma of part 3 §E with proof, the no-hub theorem, the Hodge-bundle table), lane scripts
`/tmp/hodge/results/off-paired/wave2-construct/` (`c1_identities.py`, Singular files `c2*`, `c4_Skappa.sing`, `c6*`; `run_mw.sh` and
`mw_sections_eps±1.ms` for the Mordell–Weil systems; `c5_served.py` for served-space scans with `oplib.served_targets`),
`/tmp/hodge/results/off-paired/wave2-search/oplib.py`, `w2lib.py`. Geometry: `W = {H + Cm + s⁶ − 6s⁴m + 9s²m² − 2m³ = 0} ⊂
P(1,1,1,1,1,2)` (`s = X2 + Y2` weight 1, `m = X2Y2` weight 2; the ambient singular point `[0:…:0:1]` is not on `W` since the
`m³`-coefficient is `−2`, so `W` is smooth iff quasi-smooth — check with Singular (minutes) before Step 2; `−K_W = O(1)`), `Y_∓ = {t² = ±2m³ − Cm − H} ⊂ P(1,1,1,1,2,3)` (Fano Weierstrass fibrations over `P³_x`,
`j` a function of `κ = C³/H²`), `W → Y_−` of degree 3 (`t = i·s(s² − 3m)`), fibres of `W → P³_x` genus 4 (`Aut = {1, s ↦ −s}`
generically).

**Step 1 (Mordell–Weil rank of `Y_−`, `Y_+`; local + remote, a day).** (a) Shioda–Tate on the resolution: `ρ(Y_∓) = 1`
(Lefschetz for quasi-smooth weighted hypersurfaces) and the only singular point is `p_∞ = [0:0:0:0:1:±√2]` (where all fibres
meet); resolve `p_∞` (weighted blow-up; count exceptional divisors `e`), count reducible fibres over the discriminant surface
`2C³ ∓ 27H² = 0` (generic fibre over it is nodal: no extra components except over its singular curves), and apply
`rank MW = ρ(Ỹ) − 2 − Σ(fibre components)`; alternatively compute `h^{1,1}(Ỹ_∓)` by Griffiths–Steenbrink for the weighted
hypersurface (`h^{2,2}`-type counts) and `b_2`. If `rank MW = 0`, no fibrewise-translation `φ` exists at all (torsion sections
are the 2-torsion mirrors `t = 0`, which lift to `X`). (b) Only if (a) gives rank `> 0`: search sections with denominators
`m = M/d², t = T/d³` (`deg d = 1`: `M` quartic 35, `T` sextic 84, `d` linear 4 → 123 unknowns, 455 equations of degree 12;
msolve `-g 1`, `p = 998244353`, `--cpus 8`, cap 4 h; `deg d = 2`: 314 unknowns, run only if `deg d = 1` finishes).

**Step 2 (birational rigidity of `W`; literature + one computation).** `W` is a smooth Fano fourfold hypersurface of
degree 6 and index 1 in `P(1,1,1,1,1,2)`, a triple cover of `P⁴` via the cubic in `m`. Check the literature on birational
(super)rigidity of such weighted Fano hypersurfaces (Pukhlikov; Cheltsov–Park; "triple spaces"); if rigid, `Bir(W) =
Aut(W) = ⟨g², h³⟩`-image and no `φ` exists — the route is closed rigorously. If not covered, test the maximal-singularity
candidates: compute the locus of points of `W` where the triple cover `W → P⁴_(x,s)` is totally ramified (`m`-discriminant
of the cubic has a triple root) and the curves of `W` of degree 1 in `O(1)` (lines through the cover), the natural centres
of untwisting involutions; for each centre write the candidate Geiser/Bertini-type involution explicitly (fibrewise over the
pencil of hyperplane sections through the centre) and test whether it preserves `W` (Singular, `p = 32003`).

**Step 3 (only if a `φ` exists).** Build `Z = X ×_W X^φ` as a subvariety of `X × X` (equations: `π(x1) = φ(π(x2))` cleared),
`e = (1,σ2)`, `d = (σ2,1)`, `σ = swap`, `r = (τ,τ)`; verify `σeσ = d` and `p∘d = σ2∘p` on the equations; compute the
quotient `Z/⟨σ,r⟩` fibrewise over `P³_x` (or over `W/φ`) and its canonical class by weighted adjunction; if Fano-fibred over a
curve or `p_g = q = 0`-fibred over a surface, run the book recipe: covering-curve incidence over `F_p` (polynomial curves of
weighted degree `≤ 2n` on the Fano fibre with a square incidence system and invertible Jacobian — model on
`w23_seed_search.py` of the polarized lane), polynomial normal, source identity `J·B = A − C·T^e` with a trace-zero
correction, integration, trace; Lean = `PolynomialCurves.square_system_generic_lift` + `normalized_data_from_polynomial_family`
+ a new cover-field layer + `CertificateOperations.add` for the `σ2`-averaging.

**Remote jobs.** `w3-mw-rank` (Step 1a, Singular/flint, `--cpus 4`), `w3-mw-d1` (Step 1b, msolve, `--cpus 8`, 4 h cap),
`w3-rigidity-centres` (Step 2, Singular, `--cpus 4`); lane `off-paired`; outputs `w3-*.out/.jsonl`.
**Success criterion.** Either (i) an explicit `φ` with `Z/⟨σ,r⟩` uniruled, a finite-field seed with invertible incidence
Jacobian and the source identity verified (then a holdout-style PR is in reach), or (ii) a proof that `rank MW(Y_∓) = 0`
and `Bir(W) = Aut(W)` (post to thread #57: the descent-through-reflection lemma, the Klein-cover structure, the rational
genus-2 surface, the `p_g` table of base changes, and the closure of the realisation).

## 6. Appendix

### 6.1 What this synthesis adds to the record (new, verifiable)

1. **one-reflection-skew, E4 closed (sympy exact, `/tmp/hodge/results/wave3-construct/skew_branch.py`):**
   `disc_β(R) = −6912·a⁴·δ²`, `N_{B/K}(M) = Res_β(R, b − 3aβ) = −9·a²·(3ac − b²)²`, total degrees in `(r,w,t)`:
   `a: 6, b: 5, c: 4, d: 6, 3ac − b²: 10, δ = disc_U(P): 24`. Hence the branch locus of the quartic cover is `{δ = 0} ∪ {a = 0}`,
   of the `√M` twist `{3ac − b² = 0} ∪ {a = 0}`, of `a^{1/3}` `{a = 0}`, and any Galois cover `B'' → P³` trivialising the CM local
   system has `K_{B''} ≥ π^*O(17)`: general type. The covering-curve route on `B̄''` is closed; fqlx's companion route can only
   survive on a smaller carrier.
2. **Mechanism map across all 16 lanes** (§3) with every exclusion tied to a receipt, and the honest statement that the
   book recipe fails at *dominance* everywhere except polarized-k3's `W_{23,f}` (verified seed) and, conditionally,
   off-paired's `Z/⟨σ,r⟩`.
3. **The literature status of algebraic Kuga–Satake for degree-2 K3 surfaces** (F8): nothing covers `ρ ≤ 3`,
   `End_Hdg(T) = Q`; `polarized-k3-remainder`'s K3 part remains an open GHC case.
4. The `construct-off-paired-mw` job result (no polynomial Mordell–Weil section of `Y_∓` for either sign at the further primes)
   is now recorded (§3.8) — the fibrewise-translation realisation of the descent lemma has no polynomial instance.

### 6.2 Pending jobs that would refine this ranking (collect with `rjob pull --lane <lane>`)

| job | lane | what it decides | effect on ranking |
|---|---|---|---|
| `construct-dp-j1j2` (RUNNING), `construct-dp-fibre` (resubmit `w2c_fibre.sh`) | diagonal-permutation | hidden/created second reflection on `P⁵` / `(01)`-pencil fibres | positive ⇒ §4.6 descent, jumps to #1 |
| `construct-rcr-all` (RUNNING) | reflection-coupled-residual | non-coordinate fibre reflections (Exp B), `gm_rank` (E = V₁₁?), Fermat-pencil continuation, critical points | positive Exp B ⇒ jumps to #1; `rank < 186` ⇒ a level-one sub-variation exists (revives Abel–Jacobi supports) |
| `construct-off-paired-served30` (RUNNING) | off-paired | degree-12 vertical/mirror denominators in sector (3,0) | nontrivial served class ⇒ reopens type-V1 supports |
| `construct-epc-fanoplanes` (never ran) | eight-plane | Fano scheme of planes = 8 planes ⇒ `Aut(X) = ⟨τ⟩` unconditional | none (closure only) |
| `construct-spp-vert`, `construct-spp-k7var` (never ran; `w2c_pending_jobs.sh`) | short-power | E1 of vertical families d = 7..10; generic degree-7 kernel | survivor ⇒ one fixed-`D` solve |
| `construct-ctr-ztsmooth` (never ran) | cyclic-triple | smoothness of `Z̃`, `Y`, `B` for the klt/BDPP step | none (rigor only) |
| `construct-kcc-sigmaXQ` (char-0 rerun) | klein | `Σ` on `X` over `Q̄` (currently mod p) | none (rigor only) |
| `construct-V3`, `construct-dense-rational-cones`, `-univ12` | dense-rational | linear reflections coordinate-free; cone E1 k = 9..12; `R'` E1 | none (closure only) |
| `construct-polarized-k3-remainder-w23seed` | polarized-k3 | extra `F_5` and random `F_11/F_13` seeds for `W_{23,f}` | strengthens the Part B deliverable |

### 6.3 Probability bookkeeping (for the aggregator)

Sum over lanes of P(solve) by the best route ≈ 0.30 for one wave; the three work orders carry ≈ 0.11 of it. No lane has a
route above 0.04. The dominant failure mode of the programme is structural (F2/F3/F6), not computational; the two routes
with explicit-witness payoff (WO-1, WO-2) are preferred precisely because they do not need an auxiliary uniruled variety.

### 6.4 Reproduction

* This file: `/tmp/hodge/notes/wave3-construction-ranking.md`; log `/tmp/hodge/results/wave3-construct/log.jsonl`
  (steps: start, read, read-complete, write parts 1–7, new-result, patch).
* The one computation of this synthesis: `python3 -I /tmp/hodge/results/wave3-construct/skew_branch.py` (sympy 1.12, < 1 min).
* Every other number is quoted from the lane reports and receipts named in §3; the referee reports under
  `/tmp/hodge/posts/<lane>/referee-report.md` give the VERIFIED/PLAUSIBLE status of the exclusion-side claims reused here.
