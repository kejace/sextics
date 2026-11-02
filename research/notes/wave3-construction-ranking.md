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
