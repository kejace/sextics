# Attack 1 on `off-paired` (q = 1, A = X0^3 Y0^3) — working notes

**No witness found.** Every experiment below is negative or structural. (First session crashed on 2026-11-01 ~17:45
while launching the degree-8 msolve run; resumed and completed here. All numbers are engine v2, i.e. after the 17:25
sign fix.)

Engine: `/home/kejace/devel/sextics-work/hodge` (core, fixed_d). Own code: `/tmp/hodge/results/off-paired/{oplib.py,
bilinear.py, step0_sanity.py, step1_colon.py, step3_structural.py, step4_ages.py, step5_sectors.py, step6_served.py,
step7_scan.py, step8_disc.py, step9_fermat_disc.py}`. Log of every tested denominator (669+ rows):
`/tmp/hodge/results/off-paired/log.jsonl`; bilinear sweep records `bilinear_results.jsonl`; served-target scan
`served_scan.jsonl`; msolve inputs/outputs `msolve_d{6,7,8}_c*.{ms,out}`; `delta_t_factors.txt`;
`control_target_image.out`; `step8_{minus,plus}.out`; `step9.out`.
Dossier (other agent, read in full): `/tmp/hodge/notes/dossier-off-paired.md`; its structured-support sweep lives in
`/tmp/hodge/work/sweep_*.txt` and is **not** repeated here (listed for completeness in 2.6).

Instance: `F = X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6 + X0^3X1Y0Y1/2 + X1^3X2Y1Y2/2 + X0^2X1X2Y0Y2/3`, `A = X0^3Y0^3`, q = 1.
Notation: `f6 = X0^6+X1^6+X2^6`, `g6 = Y0^6+Y1^6+Y2^6`, `Qx = F − f6 − g6` (the mixed part: a quadratic form in Y with
quartic X-coefficients `a = X0^3X1, b = X0^2X1X2, c = X1^3X2`), `H = X0^6+X1^6+Y0^6+Y1^6+X0^3X1Y0Y1/2`,
`C = X1^3Y1/2 + X0^2X1Y0/3` (so `F = H + C·X2Y2 + X2^6 + Y2^6`), `D0 = F_{X0}/(6X0)` (the unique quartic colon element),
`F_i = ∂_iF`, `E_I(F) = Σ_{i∈I} Z_i F_i`. Diagonal group `G ≅ Z/6×Z/6`, generators `g = diag(1,ζ²,1,ζ,ζ³,ζ³)`,
`h = diag(1,1,ζ,1,1,ζ⁵)`; characters written `(g-weight, h-weight)`. `A` has character (3,0); its residue is G-invariant
and anti-invariant under `σ₂: X2 ↔ Y2`.

---

## 1. The obstruction, and the mechanism a certificate needs

**Geometric meaning of a certificate.** At q = 1 the class is `Res(A Ω/F²) ∈ H^{3,1}(X)` (pole order 2; X is a
Calabi–Yau sextic fourfold, h^{3,1} = 426). A certificate with denominator D says the class is supported on the threefold
`Z = X ∩ V(D)`: it lies in the Gysin image of `H^{2,0}(Z̃)` for a resolution Z̃. A smooth ample divisor has
`H^2(Z) = H^2(X)` (algebraic), so Z must have non-rational singularities or be non-normal, and the (2,0)-forms live on
the exceptional geometry. At q = 1 only the reduced support matters (prescribed-support theorem), so the question is:
**which reduced divisors of P^5 cut X in a threefold whose resolution carries a (2,0)-form with Gysin image `[A]`?**

**The Fermat limit is served by a cone (positive control, 2.1).** For `F₀ = f6 + g6`, `D = f6` works:
`N = (A/6)(X0,X1,X2,0,0,0)`, Q = 0 (recomputed exactly). Geometrically `Z₀ = {f6 = g6 = 0}` is the join of the two
Fermat plane sextic curves `C_X`, `C_Y` (genus 10 each); its resolution is a P¹-bundle over `C_X × C_Y`, with
`H^{2,0} = H^{1,0}(C_X) ⊗ H^{1,0}(C_Y)` (dimension 100), and `A = X0^3 · Y0^3` is a product of two canonical forms
(degree-3 forms are the canonical forms of a plane sextic). The target is a **join class**, served by the cone over
`C_X` with vertex the plane `{X = 0}` (or symmetrically by `g6`).

**What `Qx` does to this mechanism — verified exactly (2.7).** On `off-paired`, `Z = X ∩ V(f6) = {f6 = 0, g6 + Qx = 0}`
is still singular along both curves `C_Y = {X=0, g6=0}` and `C_X = {Y=0, f6=0}` (Qx is quadratic in Y and quartic in X),
but the transversal types differ: along `C_Y` it is still the cone over `C_X`; along `C_X` it is the quadric cone
`Qx(x,Y) + g6(Y) = 0` in the Y-directions — an ordinary double point wherever `det Qx(x) ∝ abc ≠ 0`, i.e. away from the
18 points `x0x1x2 = 0` of `C_X`: a **rational** singularity. Blowing up the vertex plane `{X=0}`, Z̃ is a fibration over
`C_X` by the sextic surfaces `S_x = {g6(Y) + λ⁴ Q_x(Y) = 0} ⊂ P³_{λ,Y}` (p_g = 10, q = 0); with `λ ∈ O(−1)|_{C_X}`,
`Z̃ ∈ |O_P(6)|` in `P(O(−1) ⊕ O³)`, `K_{Z̃} = O_P(2) ⊗ π*O_{C_X}(4)`, Hodge bundle `f_*ω_{Z̃/C_X} = O(3) ⊕ O(2)³ ⊕ O(1)⁶`
with 46 global sections. But `h^{2,0}(Z̃)` is only the *kernel of the infinitesimal period map* on those sections (the
monodromy-invariant (2,0)-classes of the family), and the exact computation says it contributes nothing:
`target_image(F, 1, f6)` has dimension 18 = dim `J(F)_6^{odd}` — **the cone `f6` serves no nontrivial (3,1)-class at
all on `off-paired`, in any character sector** (same for `g6` and `f6+g6`). The one mixed term `Qx` (bidegree (4,2))
moves the fibers `S_x` non-isotrivially enough to kill all 100 join forms, not only the (3,3) one. **Control:** on the
Fermat `F₀` the same computation gives `target_image(F₀, 1, f6)` of dimension 118 = 18 (Jacobian) + **100**, exactly
`h^{1,0}(C_X)·h^{1,0}(C_Y)`, and `A` lies in the span (263 s). So the cone over `C_X` serves 100 classes on the Fermat
and 0 on `off-paired`. [The intermediate point `F₀ + Qx/2`: see 2.7.]

**Consequences.** A serving divisor needs a resolution with a G-invariant, σ₂-anti-invariant (2,0)-form. Geometric
candidates: (i) cones with a plane vertex over curves of genus ≥ 1 — the Fermat type, whose two coordinate instances
are dead; (ii) cones with a line vertex over surfaces with p_g > 0 (degree ≥ 4 in four variables) — `H`, `C`, the
genus-2 discriminant components `27H² ∓ 2C³`, base divisors of the genus-2 fibration `X/<h> → P³` in general; (iii)
natural divisors of the structure — `Qx`, the cubic-cover discriminant `Δ_t = 4Qx³ + 27f6²g6`, Jacobian pieces,
Hessian minors — and products. Everything in (i)–(iii) that fits in memory was tested below. None passes even E1 except
the `(block)·D0` family and Jacobian-type products, which are inconsistent. Every mechanism that solved a neighbour
(hidden reflection + Hessian, Fermat cover with Y-exchange, paired-sign cubic reduction) is a way of restoring a *second*
non-rational factor for the (3,3) target, and `Qx` breaks each (dossier §2). That is the obstruction as I understand it.

---

## 2. Experiments

Conventions: p = 1000003 unless stated; "E1" = `A·D ∈ J(F)` tested in D's character sector (necessary condition);
"full" = the complete fixed-D linear system (six numerators of degree d+1, Q of degree 2d−6) restricted to D's character
sector — no loss for eigenvector D since solutions can be averaged over G; "nullity" = dimension of the solution space of
the homogeneous system (kernel), reported even when inconsistent.

### 2.1 Sanity and controls

| test | result |
|---|---|
| `D = 1` | A ∉ J(F): the class is nonzero (rank 19 of 20 in degree 6) |
| Fermat control `F₀ = f6+g6`, `D = f6` | E1 2214/2214; **solved** exactly in 70 s: `N0 = X0^4Y0^3/6`, Q = 0 |
| `D0 = X0^4 + X0X1Y0Y1/4 + X1X2Y0Y2/9` | E1 passes; full inconsistent (reduced rank 67, nullity 8; unreduced 706/707 over 4461×765) — reproduces fqlx |
| `f6`, `g6`, `f6+g6` on off-paired | E1 fails (rank 176, augmented 177, 216×192 reduced system). Not "one condition away": the E1 cokernel in this sector has dimension 40 and `A·f6` is non-orthogonal to 34 of its 40 functionals (same numbers for `g6`); the obstruction functional has 117-term support with non-small coefficients mod p. |
| served targets of `f6`, `g6`, `f6+g6` (full `target_image`, all characters, Y-odd) | dimension 18 each = `J(F)_6^{odd}`; in A's sector exactly `{X0·F_{Y1}/6, (3/2)Y1·F_{X0}}`: **nothing nontrivial** (108–161 s each) |
| residue character | (0,0) under g and h: G-invariant; σ₂-anti-invariant. Y-even D ⇒ g-character ∈ {0,2,4}. |

### 2.2 Bilinear sector sweep (D free in a character sector): degrees 6, 7 complete; degree 8 partial

`bilinear.py`: D = Σ t_a D_a over a basis of the colon sector `(J:A)_{d,χ}` modulo the Euler gauge `F·S_{d−6}`;
N = Σ t_a N_a + Σ s_k S_k with S_k the Koszul kernel `z^e(F_j e_i − F_i e_j)`; equations = coefficients of
`B_D(N) = D div N − N(D)` modulo F (cokernel functionals of multiplication by F). Bilinear in (t,s); one random affine
normalization of t; msolve F4 over F_1073741827.

| d | sectors with nonempty colon space | (nt, ns, neq) | msolve | time |
|---|---|---|---|---|
| 6 | 22 of 36 | nt 1–6, ns 2–23, neq 141–202 | **no solution, all 22** | < 1 s each |
| 7 | 32 of 36 | nt 1–14, ns 9–45, neq 263–346 | **no solution, all 32** | 1–42 s each |
| 8 | (0,0): nt 28, ns 94, neq 557 (`msolve_d8_c00.ms`) | run with 9 GB cap / 20 min timeout | pending, see 2.7 |

Degrees 6 and 7 re-derive, by a different method and prime, fqlx's degree-6 exclusion and deanodotdotdot's degree-7
single-character exclusion. At a random t the linear-in-s system has cokernel 128–338: the equations are far from
generically satisfiable. (The six h-invariant degree-8 sectors `(a,0)` are also being run by the other agent through
Singular `modStd`, 13 GB each, at the time of writing.)

### 2.3 Structural denominators from F (degrees 7–10)

| family | d | # | E1 pass | consistent | nullity | time |
|---|---|---|---|---|---|---|
| `F_i · m`, deg m = 2 (`partial*monomial`) | 7 | 126 | 126 | 0 | 16–53 | 194 s |
| `F_{X0} · m`, deg m = 3 | 8 | 10 | 10 | 0 | 34–83 | 106 s |
| `D0 · m`, deg m = 3 (`quartic-colon*monomial`) | 7 | 56 | 56 | 0 | 18–53 | 76 s |
| `E_I(F) · m`, 1 ≤ |I| ≤ 3, deg m = 1 (`subset-euler*monomial`; includes block radial derivatives) | 7 | 210 | 210 | 0 | 30–49 | 508 s |
| 2×2 Hessian minors `F_{ik}F_{jl} − F_{il}F_{jk}` | 8 | 60 | **0** | – | – | 568 s |

Jacobian-type products pass E1 trivially and are uniformly inconsistent; Hessian minors are not even in the colon ideal.

### 2.4 Served-target scan

`oplib.served_targets`: all degree-6 targets of character (3,0) served by D, modulo `J(F)_6`. 60 Hessian minors: served
dimension 2, nontrivial **0** every time. Random colon-space elements at d = 7, 8 in every character (first session,
`step6_served.py`): none served a nontrivial sector class.

### 2.5 Geometric divisors of the genus-2 fibration and of the cubic cover (new)

| D | deg | symmetry | E1 | full | time |
|---|---|---|---|---|---|
| `27H² − 2C³` (component of the genus-2 discriminant `Δ = 4C⁶ − 729H⁴`) | 12 | eigenvector of `ker χ_(3,0)` (index 2) | **fails** (2135/2136 in 2177×3396) | – | 47 s |
| `27H² + 2C³` | 12 | same | **fails** (2135/2136) | – | 45 s |
| `Δ_t = 4Qx³ + 27f6²g6` (discriminant of the cubic cover `f6 + Qx·t + g6·t³`) | 18 | G-invariant; **irreducible over Q** (28 terms) | **passes** (sector (0,0)) | reduced system 44107 × 17390 (N 7784, Q 9606); first solve died on a 5.7 GB allocation at 120/125 GB machine load; chunked rerun `step10_delta_t.py` — see 2.7 | 19 s build |

Neither component of the base discriminant of the genus-2 fibration supports the class (one prime). The product `Δ`
(degree 24; a different reduced support) was not tested: ≈ 24k numerator + 43k Q unknowns after the 36-fold reduction.

### 2.6 Fermat-deformation denominators (new)

| D | deg | E1 | full | nullity |
|---|---|---|---|---|
| `f6·m`, `g6·m`, `(f6+g6)·m`, `Qx·m`, m a monomial of degree 1 or 2 | 7, 8 | **all 108 fail** (e.g. 213/214, 239/240, 335/336) | – | – |
| `f6·g6`, `f6·Qx`, `g6·Qx`, `f6·H`, `g6·H`, `Qx·H`, `(f6+g6)·Qx`, `(f6+g6)·H`, `f6·(f6+g6)`, `g6·(f6+g6)` | 12 | **all fail** (1063/1064) | – | – |
| `f6·C`, `g6·C`, `(f6+g6)·C`, `Qx·C` | 10 | **all fail** (646/647) | – | – |
| `f6·D0`, `g6·D0`, `(f6+g6)·D0`, `Qx·D0` | 10 | pass | **inconsistent** (rank 1123–1126 of 1308, 5103 rows) | 182–185 |

Together with the other agent's sweep (`/tmp/hodge/work/sweep_all_7_12.txt`: all squarefree products of
`{X_i, Y_i, X2Y2, X2⁶±Y2⁶, D0, H, C, R, F_i}` at degrees 7–12, all inconsistent), **every natural semi-invariant support
of degree ≤ 12 built from the visible structure of F is excluded at the one-prime level.**

### 2.7 Results that were pending when this file was first written

State when the session had to close (jobs were terminated with the session; rerun commands are in the scripts named):

* **Fermat control — done.** `target_image(F₀, 1, f6)`: dimension 118 = 18 Jacobian + 100 nontrivial, `A` in the span
  (263 s; `control_target_image.out`). On `off-paired`: 18 + 0 (2.1). The join mechanism is verified exactly at both ends.
* **Intermediate point `F₀ + Qx/2` and the two-parameter probes** (`F₀ + chain`, `F₀ + extra term`, `F₀ + first link`;
  `deform_probe.out`): started, not finished (each ≈ 4–5 min). Rerun: the heredoc jobs recorded in this session's
  transcript / `control_target_image.out` header; expected to show whether the 100 join classes die already at first
  order in the mixed term.
* **msolve on the degree-8 sector (0,0)** (122 unknowns, 557 quadrics, 9 GB cap): still running at 13.5 min / 3.4 GB when
  the session closed (20-min timeout). No output yet. Rerun: `msolve -f msolve_d8_c00.ms -o msolve_d8_c00.out -t 6`;
  the other five `(a,0)` sectors need the `.ms` files generated by `bilinear.py 8 a0`.
* **`Δ_t` full solve** (`step10_delta_t.py`, chunked compression, needs ≈ 14 GB): was waiting for the machine to have
  22 GB free (other agents held 100+ GB); not run. E1 passes, so this is the single most interesting open fixed-D
  test of the session (degree 18, G-invariant, irreducible, geometrically meaningful).

---

## 3. Witness

None. `/tmp/hodge/results/off-paired/witness-candidate.json` does not exist.

---

## 4. Mathematical ideas developed

1. **Join-class reading of the target, and its death.** On the Fermat the class is `ω_X ⊗ ω_Y` for canonical forms of the
   two plane sextic curves, served by the cone over either curve. On `off-paired` the cone over `C_X` becomes a
   non-isotrivial sextic-surface fibration whose Hodge bundle `O(3)⊕O(2)³⊕O(1)⁶` has 46 sections but — by the exact
   `target_image` computation — **no flat (2,0)-section**: `X ∩ V(f6)` supports no (3,1)-class whatsoever. A certificate
   must build a divisor whose resolution has *two* non-rational "factors" compatible with the (3,3) target.

2. **How far the cones are from the colon ideal.** `A·f6` violates 34 of the 40 cokernel functionals of the E1 map in
   its sector (2.1), so the Fermat certificate is destroyed at *leading* order by `Qx`, not perturbatively; no correction
   `f6 + δ` can be found at degree 6 (excluded anyway) and `f6·L + δ'` at degree 7 lies in the excluded bilinear sectors
   (2.2). Any "deformed cone" denominator must have degree ≥ 8 and change the support (idea 3). A deformation
   tracking along `F_t = F₀ + t·Qx` (step 7 of §5) would show whether the 100 join classes die already at first order in
   t — the `target_image` numbers at t = 1/2 (pending, 2.7) are the first data point.

3. **Mixed cones.** The (3,3) target suggests a vertex that is neither `{X=0}` nor `{Y=0}`: cones over genus ≥ 1 curves
   in a plane `P²_{X0,Y0,Z}` (`Z ∈ {X1,X2,Y1,Y2}`), e.g. `X0^6+Y0^6+Z^6`, `X0³Y0³ + Z⁶`, or over surfaces in a P³ mixing
   X and Y other than the genus-2 base `P³_{X0X1Y0Y1}` (whose natural divisors `H`, `C`, `27H²∓2C³` are all excluded).
   G-orbit products make any such cone a legitimate semi-invariant support of degree ≤ 36·deg. Not run (memory).

4. **Discriminant divisors.** Two natural discriminants: the genus-2 base discriminant `Δ = (2C³−27H²)(2C³+27H²)`
   (components excluded) and the cubic-cover discriminant `Δ_t = 4Qx³ + 27f6²g6` (irreducible, degree 18; pending).
   Next in line: the univariate discriminants `disc_{X2}F`, `disc_{Y2}F` (degree 30) and the full Hessian (degree 24).

5. **Degree-8 closure is within reach.** The bilinear sector formulation (2.2) is complete for single characters and
   closes degree 6–7 in minutes; degree 8 (six open sectors, 110–122 unknowns, 525–557 quadrics) is a matter of msolve
   time/memory, duplicating the other agent's Singular run with a different engine.

---

## 5. Ranked next steps (estimates for this 40-core / 125 GB machine; everything one-prime unless kernels are exported)

1. **Finish degree 8** for the six `(a,0)` sectors with msolve (≤ 1 h and ≤ 10 GB each) → "no degree-8 denominator
   of any single G-character" (with deanodotdotdot's 30 sectors). Mixed-character degree 8 (436 D-coordinates) is
   hours–days.
2. **Mixed-cone served-target spaces** (idea 3): `target_image` for `X0^6+Y0^6+Z^6`, `X0³Y0³+Z⁶`, `X0^6 + Y0^6 + X1^6`,
   `H`-type quartics in `(X0,Y0,Z,W)`, and their G-orbit products: minutes each at degree 6, ≤ 1 h / ≤ 10 GB at degree
   12–18 with the compressed solver. Any D with a nontrivial served (3,0)-class would be the first such object found.
3. **Deformation tracking along `F_t = F₀ + t·Qx`** (idea 2 / step 7 merged): `target_image(F_t, 1, f6)` and the
   rank of the `f6`-system for t ∈ {1/2, 1/10, 1/100, …} and for the two-parameter family scaling the chain and the extra
   term separately (the chain alone is the solved `holdout`); seconds per point; tells whether the join classes die at
   first order and whether the holdout's serving divisor degenerates continuously.
4. **`Δ_t` full solve** if E1 passes (≈ 7k + 9k unknowns, 43k rows; compressed-light solver ≈ 1 h, ≈ 6 GB) and the
   degree-30 univariate discriminants (E1 cheap; full solve needs a sparse solver or the 512 GB box).
5. **`Δ`, `H·Δ`, `H·C·Δ`** (degrees 24–34; base-divisor supports not covered by the thread's obstructions): same solver
   requirement; alternatively settle them conceptually via the κ-fibration obstruction (dossier E5, a day of exact algebra).
6. **Nonlinear search in the fully symmetric sector at degrees 9–12** (dossier E4): (50,161,859) … (185,681,2535)
   bilinear systems; msolve or LM; minutes–hours per degree.
7. **Deformation tracking**: served-target space of `f6` and the rank of the `f6`-system along `F_t = F₀ + t·Qx`
   (seconds per t) to see whether the join classes die at first order in t or only at a finite set of t.
