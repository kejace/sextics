# External work: DiamondJM/hodgeHoldout (branch `claude/inspiring-turing-tx2pbu`)

Source: clone at `/tmp/hodge/ext/hodgeHoldout`, branch `claude/inspiring-turing-tx2pbu`
(HEAD `af0589e`, 2026-10-23; `origin/main` is at `bfb2fb9`, 2026-10-20, and does **not** contain
`sweep/`, `skew/`, `cycleq2/` — those exist only on the claude branch). All work there was
"prepared with Claude Code assistance on behalf of `DiamondJM`". Reviewed 2026-11-01.

Status legend used below: **[proven]** = exact identity / exact rational linear algebra /
short proof sketch that I read and find sound; **[modular]** = computation mod one prime,
evidence only; **[claim]** = asserted in the repo but not reproducible from committed logs;
**[failed]** = attempted, not completed.

---

## 1. The "1,300 natural degree-6 supports" sweep (`sweep/`)

### 1.1 What was swept

Files: `sweep/fixed_d.py` (solver), `sweep/run.sh` (driver), `sweep/instances.txt` (16 ids),
`sweep/instances_q1.txt` (the 15 with q=1), `sweep/logs/<id>.log`, `sweep/RESULTS.md`.
Commits `f5e2e24` … `fc57193`, all on 2026-10-22 (18:35–22:55 UTC).

The 16 instances are **exactly the 16 currently open instances** of the catalogue
(verified against `data/*/proof.json` in the sextics repo on 2026-11-01):

| instance | q | A | candidates | run | solvable |
|---|---|---|---|---|---|
| cycle-q2 | 2 | X0^2X1^2X2^3Y0^2Y1^2Y2 | 59 | **not run** (killed, >15 GB) | – |
| cyclic-triple-residual | 1 | X0^3Y0^3 | 71 | 71 | 0 |
| dense-rational | 1 | X0^3Y0^3 | 106 | 106 | 0 |
| diagonal-permutation | 1 | X0^3Y0^3 − X1^3Y1^3 | 77 | 77 | 0 |
| eight-plane-configuration | 1 | X0^3Y0^3 | 106 | 106 | 0 |
| four-plus-two-residual | 1 | X0^3Y0^3 | 86 | 86 | 0 |
| klein-coupled-character | 1 | X0^3Y0^3 | 53 | 53 | 0 |
| off-paired | 1 | X0^3Y0^3 | 57 | 57 | 0 |
| one-reflection-skew | 1 | X0^3Y0^3 | 55 | 55 | 0 |
| polarized-k3-remainder | 1 | X0X1X2Y0Y1Y2 | 100 | 100 | 0 |
| quadratic-cubic-pullback | 1 | X0^3Y0^3 | 97 | 97 | 0 |
| quadric-contact-residual | 1 | X0^3Y0^3 | 79 | 79 | 0 |
| rank-one-biquartic | 1 | X0^3Y0^3 | 101 | 101 | 0 |
| reflection-coupled-residual | 1 | X0^3Y0^3 | 61 | 61 | 0 |
| short-power-presentation | 1 | X0^3Y0^3 | 106 | 106 | 0 |
| symmetric-sextic-residual | 1 | X0X1X2Y0Y1Y2 | 65 | 65 | 0 |

Totals: 1,279 candidates enumerated, **1,220 actually tested (all q=1), 0 solvable**.
"About 1,300" in the posted comments refers to this.

### 1.2 Denominator families (all of degree exactly 6)

`candidates(F)` in `fixed_d.py` generates, per instance:

1. **`powS`** — sums of sixth powers of a k-subset S ⊂ {X0,X1,X2,Y0,Y1,Y2}, k ∈ {2,3,4}:
   C(6,2)+C(6,3)+C(6,4) = 15+20+15 = **50 candidates, identical for every instance**
   (e.g. `pow01` = X0^6+X1^6, `pow012` = X0^6+X1^6+X2^6, `pow0345` = X0^6+Y0^6+Y1^6+Y2^6).
   Note `pow` with k=2 is a product of six hyperplanes; k=3,4 are Fermat-type joins.
2. **`FsubS`** — "coordinate sub-equations": for S ⊂ {0..5}, |S| ∈ {2,3,4,5}, the sum of those
   monomials of F supported on the variables in S, kept only if it has **more than |S| terms**
   (i.e. contains at least one coupling monomial, not just the pure powers). Duplicates (same
   polynomial from different S) removed. Count varies with the coupling structure of F:
   3 (klein-coupled-character, off-paired… have few couplings) up to 56 (dense F with 236 terms:
   dense-rational, short-power-presentation, eight-plane-configuration).

No other shapes: **no products, no pencils, no generic coefficients, no degree ≠ 6, no
semi-invariant general D in a character, no D depending on the target A.**

### 1.3 Method (exact description of `fixed_d.py`)

Cleared certificate condition, taken verbatim from `formal/README.md`:

    A·D² − q·D·Σ_i N_i ∂F/∂Z_i + F·(D·Σ_i ∂N_i/∂Z_i − Σ_i N_i ∂D/∂Z_i) = F^{q+1}·Q.

For fixed D this is **linear in (N_0..N_5, Q)**. Degrees: deg N_i = deg D + 6q − 5 (= 7 for
d=6, q=1), deg Q = deg A + 2·deg D − 6(q+1) (= 6 for q=1, deg A = 6).

* **Unknown space**: the *full* monomial space for each N_i and for Q, restricted only by
  Y-parity: a monomial m is admitted in N_i iff ypar(m) ≡ 1 + ypar(D) + [i ≥ 3] (mod 2), and in
  Q iff ypar = 1. This makes every term of the identity Y-odd. This is a legitimate exact
  isotypic splitting (not a restriction) **provided F is Y-even and A is Y-odd**. I verified
  with sympy from the 16 `instance.rec` files on 2026-11-01: every F has Y-degrees in {0,2,4,6}
  and every A has Y-degree 3 (cycle-q2: 5). So the parity cut loses nothing for all 16.
  ypar(D) is read off the *first* monomial of D (fine: all candidates are Y-homogeneous-parity
  because F is Y-even).
* **Linear algebra**: dense `nmod_mat` over F_P, **P = 268435399** (single prime), via
  python-flint. Columns = images of each admitted monomial (N_i-columns:
  −q·D·m·∂F/∂Z_i + F·(D·∂m/∂Z_i − m·∂D/∂Z_i); Q-columns: −F^{q+1}·m); right-hand side A·D².
  Rows = all monomials occurring. Verdict `SOLVABLE` iff rank[M | A·D²] = rank M mod P.
  Typical sizes for q=1, dense F: rows = 16,797 (Y-odd degree-18 monomials), cols = 2,602
  (6 × ~396 degree-7 monomials of the right parity + ~226 Y-odd degree-6 Q-monomials),
  rank ≈ 2,304–2,307 (kernel ≈ 295–298 = trivial syzygies). Sparse-structured F give the same
  column count and fewer rows.
* **No symmetry/character block reduction** inside `fixed_d.py` (unlike `packet/cert.py`,
  which splits by the diagonal group of the holdout sextic). The sweep is therefore a
  *complete* fixed-D test for the given (F,q,A) and D, modulo the single prime.
* **Evidence status** [modular]: a `no` mod P means rank_P M < rank_P [M|rhs]. A rational
  solution would reduce mod P to a solution unless P divides a denominator of the reduced
  echelon form; false negatives have heuristic probability ~1/P per candidate. A `SOLVABLE`
  mod P would likewise not be a proof. No exact-over-Q rerun was done.
* **q=1 prescribed-support theorem** (docs/01): for q=1, failure of the *complete* fixed-D
  test excludes the whole class's support on V(D) with any multiplicity. The sweep's negatives
  would therefore become genuine exclusions of each V(D) if redone exactly over Q (as
  `packet/exact_support.py` does for holdout). This upgrade was **not** done for the 16 open
  instances.

### 1.4 Timing / memory

* Per candidate (q=1, degree 6): 34–81 s (dense `nmod_mat` rank of ~16.8k × 2.6k, twice).
  Per-instance totals: 1,581 s (one-reflection-skew) … 5,232 s (dense-rational).
  Sum over the 15 q=1 instances: **46,757 s ≈ 13.0 CPU-hours**, run with `xargs -P3`
  (3 processes), wall ≈ 4.3 h (consistent with commit timestamps 18:35 → 22:55 UTC).
* Memory per q=1 process: two dense 16.8k × 2.6k matrices of 64-bit entries ≈ 0.7 GB.
* **cycle-q2 (q=2)**: deg N_i = 6+12−5 = 13, deg Q = 11+12−18 = 5; rows ≈ 49k (Y-odd
  degree-23 monomials), cols ≈ 6 × 4,284 + 126 ≈ 25.8k; dense ≈ 10 GB per matrix.
  Log says only `Killed`; comment says "exceeded 15 GB". **Not run.** A sparse or
  character-split solver would make it feasible (the 72-element diagonal group of cycle-q2 splits
  the columns into ≤72 blocks).
* The host was `/home/user/hodgeHoldout` with sextics at `/home/user/sextics` (paths are
  hard-coded in `run.sh` and `skew/supp.py`). Machine specs not recorded; RAM limit implied
  ≈ 15 GB (also "exceeded 14 GB" and "9 GB" aborts elsewhere).

### 1.5 Controls

* **[claim]** `RESULTS.md`: "Controls: the Fermat join D = X0^6+X1^6+X2^6 serves X0^3Y0^3;
  X1^6+X2^6+Y1^6 does not." This positive control is **not in any committed log of
  `fixed_d.py`**; the committed Fermat controls (`packet/exact_support.log`) were produced by
  the different code path `packet/cert.py`/`exact_support.py` (exact fmpq, character-block,
  whole six-dimensional (3,1)-sector; supported-sector dims 4 and 3 as predicted). So the
  sweep script itself has no committed positive hit. I attempted to rerun that positive control
  through `fixed_d.py` on 2026-11-01 but execution of the cloned code was not permitted in my
  session; it remains to be done (≈1 min, see §4).
* **Negative calibration on a solved instance** [modular], `skew/supp.py`: on **one-reflection**
  (solved by `fqlx` with a degree-30 support D₊ = 9δZ + 4h³T²), the degree-6 candidates
  s²+T², S²+T², T²−s²/3, T²−S²/3 (s = Y0³−Y1³, S = Y0³+Y1³, T = Y2³) carry **none** of
  A, A₊, A₋ (A± = X0³(Y0³ ± Y1³)/2). Same for one-reflection-skew. Conclusion drawn by the
  author and which I endorse: **degree-6 negatives are weak evidence on this family**; the
  accepted supports can have degree 30.

### 1.6 Conclusions (author's and mine)

* No natural degree-6 support (as defined in §1.2) carries the target for any of the 15 open
  q=1 instances, mod 268435399. [modular]
* This says nothing about degree > 6, about products/unions of these candidates, about
  pencils/generic members of families, or about cycle-q2.

---

## 2. Other work in the repository (not the sweep)

### 2.1 `holdout` (solved; packet/ — posted to discussion #44, comment `DC_kwDOVLOOJc4BIp9D`)

F = ΣZ_i^6 + X0^3X1Y0Y1 + X1^3X2Y1Y2, q=1, A = X0^3Y0^3. Diagonal group G of order 1296
(216 mod scalars), characters c1=(3,1,0,1,1,0), c2=(0,3,1,0,1,1).

* **[proven, exact over Q]** `exact_support.py`: complete fixed-D systems (character block, full
  numerator space, unrestricted quotient) for the whole 6-dim (3,1) target sector
  {X0³Y0³, X1³Y1³, X2³Y2³, X1²X2Y0²Y2, X2Y1⁴Y2, X1Y0⁴Y1}: for every D listed the augmented rank
  rises by exactly 6 ⇒ no class of the sector supported on V(D) (any multiplicity, by the q=1
  prescribed-support theorem). Degree 6 (50 unknowns): D0 = X0⁶+X1⁶+X2⁶, D0+C1+C2, Y-block,
  Φ0 = X0⁶+Y0⁶+C1, R = X1⁶+X2⁶+Y1⁶+Y2⁶+C2, X0⁶+Y0⁶, X1⁶±Y1⁶, X2⁶−Y2⁶, U0²+U1² (= F−Ψ), C1+C2,
  and all 24 "3+3 block sums containing X0 with any subset of {C1,C2}". Degree 12 (586 unknowns):
  432Φ0²−(X1Y1)⁶, D0·Qy, X2¹²−X1⁶Y1⁶, (X2Y2)⁶−(X1Y1)⁶, (X0X1)⁶−(Y0Y1)⁶.
* **[modular, P = 268435399]** `bilin2.py`: all Y-even G-semi-invariant D, per character block,
  parametrise E1 (ΣN_i∂_iF = A·D) by nullspace, reduce E2 (D·divN − N(D) ∈ (F)) mod F to quadrics
  in the free parameters, dehomogenise by random linear forms vanishing on F·S_{d−6} (excludes
  D ∈ (F)), msolve Gröbner mod P ⇒ unit ideal for **every Y-even character block at degrees
  8 and 9** (95 and 108 admissible blocks). `bilin2_tau.py`: X2↔Y2 (anti)symmetric sectors
  through degree 12 (ε=−1) and 9 (ε=+1); **degree-10 ε=+1 blocks (~63 free parameters)
  exceeded 14 GB — inconclusive** [failed]. `ansatz.py`: N = A·g, D = g(F), all semi-invariant g
  through degree 11: empty. `ansatz2.py`: N = A·g + Koszul (complete for D ∈ J(F)) through
  degree 9: empty. `pencil.py`/`pencil_run*.py`: 6,200 binomial pencils D = m1^k + b·m2^k
  (deg m_i ≤ 6, all b ∈ F̄_P via projected-determinant interpolation) through degree 12: empty.
* **[modular]** `gm.py`/`gm_all.py`: Gauss–Manin along ΣZ⁶ + t(C1+C2) at t0 = 86928197 mod P:
  for each of the six (3,1) seeds and several (2,2) seeds, the derivatives of [A·Ω/F²] span the
  full rank-19 character block, Hodge numbers (6,9,4). An earlier version dropping pole-order
  components was withdrawn.
* Subsequently superseded: holdout was solved by `fqlx` (PR #128, moving conic on the Fermat
  cover (Y0,Y1,Y2) = (αT, y, βT), T⁶ = Y0⁶+Y2⁶).

### 2.2 Holdout-family extension (`tools/gen_holdout_targets.py`, `proposals/holdout-*`, `comment-holdout-extension.md`)

**[proven, sympy-checked identities]** fqlx's construction is target-independent except for the
quotient source and the deck-trace companion. For a Y-odd monomial A = X^a Y0^{b0}Y1^{b1}Y2^{b2},
with k = b0 + b2 + 1 − b1: J_c·B = A − A·(y/T)^k, J_c = 2α⁵(T²−y²)/T, B polynomial of weighted
degree 5 iff k even; companion killed by the trace iff 6 ∤ k. Table of (A, k, B) for the nine
monomial holdout targets (k ∈ {−2, 2, 4}); `holdout-bundle`'s second summand X0³Y0Y1² has k=0
and is outside this route. Lean proposal folders were generated for all nine plus
`one-reflection-y2-sector` (via `gen_oneref.py`: same cover with β=α, new conic seed = real root
of a degree-8 polynomial g in (−13,−12)). **Outcome**: all nine holdout-* targets and
one-reflection-y2-sector are now SOLVED in the catalogue; the accepted holdout-* proofs credit
"OpenAI Codex assistance for `benediamond`" (not DiamondJM), and one-reflection-y2-sector
credits fqlx's conic family; `pr-body-oneref.md` on the branch is DiamondJM's PR body for it.
`tools/seed-search/`: brute-force conic seed search for `holdout-ternary-deformation` **[failed]**
(msolve F4 reached degree 14, 77k×199k, no finish in 20–25 min over Q or mod 65521); that
instance is now solved by someone else with a "nonsingular modulo-five seed".
`scan_cover.py`: of all open instances only `one-reflection-y2-sector` had the Fermat-cover
symmetry (so none of the current 16 do).

### 2.3 `one-reflection-skew` (`skew/`, posted to discussion #127, comment `DC_kwDOVLOOJc4BI_aX`)

F = ΣZ⁶ + (X0³X1 + X1³X2)p + (1/7)X0X1·p·t, p = Y0Y1, t = Y0²−Y1², A = X0³Y0³, q=1.

* **[proven]** Diagonal symmetry group has order **24 mod scalars** (216 for one-reflection),
  generated by τ (Y→−Y), σ: (X1,Y1)→(−X1,−Y1), μ6 on Y2. Hidden reflections Y0↔±Y1 negate the
  skew term; no diagonal X-twist repairs it (needs a+b≡3, 3a+b≡0 mod 6, impossible by parity).
  σ survives ⇒ split A = A₊ + A₋ still available.
* **[proven]** Fermat cover T⁶ = Y0⁶+Y2⁶: the scaled exchange preserves the skew image only if
  α⁶ = −1 — a constant, not the generic point.
* **[proven, sketch + Gröbner]** No cubic-surface reduction F = G(u, ℓ1, ℓ2, ℓ3) with one active
  X-coordinate u = X0², passive X-ratios, weight-2 Y-rational ℓ_i, and density ∝ Y0³ (or S, s):
  (i) if two ℓ_i are Y2-free the third's μ6-orbit has size ≤ 3 ⇒ wrong Y2-character of the
  Jacobian; (ii) otherwise ℓ2ℓ3 = pt − cp², eliminating gives Y2⁶ = Φ(ℓ2) with Φ Laurent ⇒ Φ has
  one pole ⇒ ℓ2 a polynomial in Y2 of degree ≥ 2 ⇒ det ∂(p,ℓ2,P4/ℓ2)/∂Y ∝ (∂ℓ2/∂Y2)/ℓ2 depends
  on Y2. Contradiction. Gröbner checks in Singular (char 32003, Y0=2, Y1=3): Möbius ℓ2 and
  degree-1 ℓ2 give ⟨1⟩; degree-2 polynomial ℓ2 gives a nonempty variety (killed by step (ii));
  a (2,2) rational ansatz **exceeded 9 GB [failed]**. **Not excluded**: two active X-coordinates,
  or passive ratios mixing X and Y.
* **[proven]** Fibration over (r, w, x = Y0/Y1): fibre is the plane sextic g(X0², Y1²) + Y2⁶ = 0;
  the residue character (+1, −1, ζ) selects a 2-dim purely (1,0) eigenspace {X0Y1², X0³} of
  H¹(fibre) — signature (2,0) unitary piece of a μ6-cover of P¹, not an elliptic curve. Identical
  for the solved one-reflection, so neither obstructs nor serves.
* **[modular]** `fqlx_cert.py` reproduces fqlx's strict half-certificate (N₊, D₊) mod 268435009
  (a prime with √3): N₊(F0) = A₊D₊ and D₊divN₊ − N₊(D₊) = 0 exactly; D₊ has 725 terms, N₊ 352–546
  per component. First-order deformation along the skew term sized at ~1.13M + 2×0.16M unknowns
  vs ~4.1M degree-60 rows (D₊ semi-invariant only under τ) — **not attempted [failed/too big]**.
  Side remark: s²+T² does not support A±, so fqlx's support is not the tritangent-plane curves
  {u = ρ_k p, s = ±iT}.

### 2.4 `cycle-q2` (`cycleq2/`, discussion #61, comments `DC_kwDOVLOOJc4BI_1U`, `…BJC7I`, `…BJD3V`)

F = ΣZ⁶ + ½(X0³X1Y0Y1 + X0X2³Y0Y2 + X1³X2Y1Y2), q=2, A = X0²X1²X2³Y0²Y1²Y2.

* **[proven]** `structure.py`: with α=Y0/X0, β=Y1/X1, γ=Y2/X2 passive and U,V,W = X0²,X1²,X2²:
  F = (1+α⁶)U³ + (1+β⁶)V³ + (1+γ⁶)W³ + ½(αβU²V + βγV²W + αγUW²), A = α²β²γU²V²W²: a family of
  plane cubics over a 3-dim base (as in cycle-three), each pair map Jacobian 2.
  Diagonal group **72 mod scalars**; residue character of A takes values ζ6^{0,2,4}, trivial on an
  index-3 subgroup of order 24 (contains τ and all paired signs), primitive cube root on the
  quotient μ3. holdout-q2's endpoint reduction (cubic surface + 27 lines) has no analogue: each
  pair enters through both ZW and Z³W, so an active pair needs (X2², X2Y2, Y2²) ⇒ (2,3) complete
  intersection (K3) fibres, no cubic surface. (An earlier claim that cycle-three's model is the
  μ3-quotient was corrected: its passive coordinates c1, c2 pick up ζ^∓2.)
* **[proven]** `second_variation.py`: 2A·Ω/F³ = ∂_{s1}∂_{s2}(Ω/F_s)|₀, F_s = F + s1·X0X1X2Y0Y1Y2 +
  s2·X0X1X2²Y0Y1 (m1m2 = A); both directions add an xyz term to cycle-three's cubic
  (g1 = v³c1c2, g2 = v²c1 in the chart v = Y0/X0, c1 = (Y1/X1)/v, c2 = (Y2/X2)/v).
  `invariants.py`: Aronhold s4, s6 of ax³+by³+cz³+dx²y+ey²z+fz²x+m·xyz, normalised to
  cycle-three (formulas in NOTES.md / comment); point-count checked mod 101–109 [modular].
* **[proven bookkeeping]** Variations of A = s4, B = −s6/4: A1,B1 carry v³; A2,B2 carry v²;
  A12 = 0; B12 carries v⁵. With the multiplier −144Rv³ the transported source carries v² = T3/R:
  a function on the **triple cover** M3 ((r:t) = (ρ³:τ³), τ/ρ = v²), not on cycle-three's M1.
  M3 fibres over fixed (c1,c2): elliptic surfaces S3 in P(1,1,6,9), coefficients degrees (12,18),
  minimal, χ=3, p_g=2, nonconstant j. Deck τ→ωτ: invariant sector = M1 (p_g=0), each nontrivial
  sector h^{2,0}=h^{0,2}=1. ⇒ cycle-three's "every polynomial source has a primitive" is false here.
* **[proven, exact over Q at (c1,c2) = (2,3),(3,5),(1,4)]** Weighted Jacobian ring of S3:
  dim R_1 = 2, R_19 = 28, R_37 = 2, R_38 = 1 (socle), R_39 = 0. Target sector (τ-exponent 0 mod 3)
  Hodge numbers (1,10,1); invariant sector (0,8,0). The cubic-pole numerator P (weighted degree
  33 ≤ 37) is **not** in J_37 (nonzero fibrewise (0,2) part — no fibrewise primitive), but
  P ∈ J_37 + ∂_{c1}G·S_19 + ∂_{c2}G·S_19 (Kodaira–Spencer absorbs it) — no obstruction at (0,2).
  KS: R_1 → R_19 injective (image 2), R_19 → R_37 onto.
* **[modular, P = 1000003]** `gm_engine.py`/`gm_jets.py`: target-sector Gauss–Manin connection in a
  monomial normal-form frame, jets to any order; flatness ∂1Γ2 − ∂2Γ1 + [Γ1,Γ2] = 0 holds
  exactly. Span of c-derivatives of s = Res(ρΩ/G) up to orders 0..6: **1, 3, 6, 10, 12, 12** at
  (2,3), (3,5), (7,11) ⇒ holomorphic form generates all of V_ω (rank 12): **no algebraic classes
  in the target sector**. (A rank-4 result from a buggy ring missing −x³ was withdrawn.)
* **[failed]** `degree_probe.py`: connection entries along c2=3 do not fit rational functions of
  degree ≤ ~130 in c1 (140 points); spurious frame poles ⇒ naive 2-variable reconstruction
  infeasible.
* **Open (author's reformulation)**: is Σ ∈ ∇1(K¹²) + ∇2(K¹²), K = Q(c1,c2), for the rank-12
  system — a 2-variable rational-integrability / creative-telescoping problem; by Gysin, removing
  base curves only helps along Noether–Lefschetz curves of the target sector. Caveats still open
  for a certificate: s-independence of the multiplier −144Rv³ with the xyz term; second-order
  moving-map correction.

---

## 3. What was NOT covered (by the sweep and by the repo as a whole)

1. **cycle-q2 fixed-D sweep**: nothing run (memory). A character-split or sparse solver is needed.
2. **Degrees ≠ 6** for the 16 open instances: no degree-7..12 denominators, no degree-12 products
   of the degree-6 candidates, no degree-30-type supports (the kind fqlx actually used on
   one-reflection). The deeper searches (semi-invariant D per character at degrees 8–9,
   pencils, Koszul/ansatz families) exist **only for holdout**, not for any open instance.
3. **Denominator shapes**: no pencils/generic coefficients, no sums of cross-block couplings not
   already sub-sums of F, no D built from A or from invariants of the symmetry group, no
   non-Y-even D (a Y-odd D would be allowed by the formalism; the parity cut assumes D Y-even).
4. **Target sector**: only the catalogue's A was tested, not the full isotypic (3,1)-sector of A
   (for holdout this was done exactly; for open instances the prescribed-support conclusion would
   need the whole sector or at least an exact rerun).
5. **Exactness**: all 1,220 negatives are mod one prime; none re-run over Q (fmpq) — so none is a
   theorem yet, though each would be cheap (< 1 min each per the holdout timings at degree 6).
6. **Positive control of `fixed_d.py` itself** is not in the committed logs (see §1.5).
7. For one-reflection-skew: cubic-surface reductions with **two active X-coordinates** or passive
   ratios mixing X and Y; the ~1.4M-unknown first-order deformation of fqlx's certificate.
8. For cycle-q2: the global de Rham/creative-telescoping question; Noether–Lefschetz search; the
   two certificate caveats. No Lean.
9. No work at all on: dense-rational, polarized-k3-remainder, klein-coupled-character,
   symmetric-sextic-residual, rank-one-biquartic, four-plus-two-residual, quadric-contact-residual,
   reflection-coupled-residual, cyclic-triple-residual, short-power-presentation,
   quadratic-cubic-pullback, eight-plane-configuration, diagonal-permutation, off-paired beyond
   the degree-6 sweep row.

---

## 4. Reusable code (paths under `/tmp/hodge/ext/hodgeHoldout/`)

Dependencies: Python ≥ 3.11, python-flint 0.9.0; msolve 0.6.5 and Singular 4.3.2 for
`packet/bilin2*.py`, `packet/ansatz*.py`, `skew/groebner/*.sing`; sympy/mpmath for cycleq2 and
seed-search. Note: these are files from an external clone — run them from a separate scratch
directory with `python3 -I` and pass paths as arguments.

| path | what it does | notes |
|---|---|---|
| `sweep/fixed_d.py` | **Generic fixed-D solvability test for any catalogue instance** (any q). `parse_rec(path)` parses `instance.rec` (q, F, A; handles `/` as exact fmpq) into flint `fmpq_mpoly` in variables X0..Y2; `monoms(n)`; `test(F,q,A,D) -> (solvable, rows, cols, rank)` mod P=268435399 with Y-parity cut; `candidates(F)` builds the pow/Fsub list; CLI `python3 fixed_d.py <data-dir> <id> [max]`, resumable via `logs/<id>.log`. | Dense `nmod_mat`; ~40 s / candidate at d=6, q=1; infeasible for q=2 at d=6. Easy to swap P or to replace `nmod_mat` by `fmpq_mat` for an exact run. |
| `sweep/run.sh` | `xargs -P3` driver over `instances_q1.txt` | hard-coded `/home/user/...` paths |
| `packet/cert.py` | Holdout-specific but structurally reusable: diagonal-group characters (`char`, `char_add`, `poly_char`, `monoms_char`), `build_system`/`try_D` for E1/E2 split into character blocks, `PRIME = 2^31−1` default | group hard-coded as c1=(3,1,0,1,1,0), c2=(0,3,1,0,1,1); F_HOLD, A hard-coded; variables x0..y2 |
| `packet/exact_support.py` | **Exact fmpq_mat sector test**: rank(M) vs rank([M \| A_1D … A_6D]) over Q for a whole target sector; gives "supported sector dimension" | the right template for turning sweep negatives into prescribed-support theorems; depends on `cert.py` |
| `packet/koszul.py` | fixed-D test with Koszul parametrisation of the E1 kernel and E2 reduced mod F by division — scales to higher degree than `cert.try_D` | holdout-specific imports |
| `packet/bilin2.py`, `bilin2_tau.py` | unknown-D search per (degree, character): E1 nullspace → E2 quadrics → msolve GB mod P → Singular decomposition, discarding D ∈ (F) | the most elaborate solver in the repo; needs msolve + Singular |
| `packet/ansatz.py`, `ansatz2.py` | N = A·g (+Koszul), D = g(F) families | |
| `packet/pencil.py`, `pencil_run*.py` | pencil D = D1 + b·D2: special b as roots of det of projected systems | |
| `packet/span31.py` | for fixed D, the subspace of (3,1) monomial targets admitting a certificate mod p | |
| `packet/gm.py`, `gm_all.py` | one-parameter Gauss–Manin by Griffiths–Dwork over F_p[u]/(u^N) in a character block; rank of the sub-local-system generated by [AΩ/F²], Hodge numbers | holdout group hard-coded |
| `cycleq2/gm_engine.py`, `gm_jets.py` | mod-p (P=1000003) Gauss–Manin connection + arbitrary-order jets for weighted elliptic surfaces y² = x³ + a(ρ,τ;c)x + b(ρ,τ;c) in P(1,1,6,9), sector-split by deck action; usage `python3 gm_jets.py c1 c2 K` | reusable for any 2-parameter family of such surfaces |
| `cycleq2/jacobian_ring.py`, `ks_test.py`, `sectors.py` | weighted Jacobian ring dimensions, (0,2)-test P ∈ J_37, Kodaira–Spencer absorption test, sector dimensions | exact over Q (flint) |
| `cycleq2/invariants.py`, `inv_check.py` | Aronhold s4/s6 of a ternary cubic with xyz term, normalised to cycle-three; point-count checker | |
| `cycleq2/gm.py`, `step2.py`, `solve_rational.py` | fibrewise Griffiths–Dwork on η² = ξ³ + aξ + b over Q(t); rational-solution search for G' + G/(3t) + MG = H with denominators t^jΔ^k | |
| `cycleq2/structure.py` | paired-ratio cubic + full diagonal-symmetry/character enumeration (brute force over μ6⁶) | trivially adaptable to any instance: generic character-table tool |
| `skew/fqlx_cert.py` | rebuilds fqlx's one-reflection half-certificate (N₊, D₊) mod 268435009 and checks both identities | template for reproducing accepted certificates mod p |
| `skew/supp.py` | runs `fixed_d.test` on arbitrary hand-written (D, A) pairs for one-reflection(-skew) | |
| `skew/groebner/*.sing` | Singular scripts for the pole lemma (char 32003) | |
| `tools/gen_holdout_targets.py` (903 lines), `gen_oneref.py` | generate Lean proposal folders from fqlx's holdout proof for a new monomial target (replaces source/companion/trace lemmas; sympy-checks identities) | now moot (all targets solved) but documents which 6 modules are target-specific |
| `tools/lean_compile.py` | compiles a proposal folder the way the sextics verifier does (mirrors `scripts/proposals/sandbox.py`) | |
| `tools/seed-search/scan_cover.py` | tests every open instance for the holdout-type Fermat-cover symmetry (y ↔ T) in either block | result: none of the current 16 has it |
| `.github/workflows/post-discussion.yml`, `open-pr.yml`, `pr-status.yml`, `pr-diagnose.yml`, `lean-check.yml` | GitHub Actions to post/update/reply to hodge.computer discussion comments (needs `HODGE_DISCUSSION_TOKEN`), open/close PRs, report checks | |

## 5. Suggested immediate uses

* Rerun the 1,220 sweep negatives **exactly over Q** (swap `nmod_mat` → `fmpq_mat` in
  `fixed_d.test`; or use `exact_support.py`'s pattern) to turn them into genuine
  prescribed-support exclusions for the 15 q=1 instances — and first run the Fermat positive
  control through `fixed_d.py`.
* Run the cycle-q2 degree-6 sweep with a character-block split (72-element group) or a sparse
  mod-p solver.
* Port `bilin2.py`'s unknown-D per-character search (degrees 7–9) from holdout to the open
  instances; it needs the instance's own diagonal group in place of c1, c2.
