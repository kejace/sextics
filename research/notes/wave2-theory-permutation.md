# Wave 2 theory report — group "permutation"

Lanes: `symmetric-sextic-residual` (SSR), `diagonal-permutation` (DP),
`polarized-k3-remainder` (PK3R), `four-plus-two-residual` (FP2R).
Author: theory-permutation agent, 2026-11-02 (on VM hodge-sextics).
Work dir `/tmp/hodge/results/_theory-permutation/`, log `log.jsonl`.
Scripts there (run with `python3 -I`). Nothing posted to GitHub; nothing in the repo modified.

Status labels: **[P]** proved / exact over Q; **[C]** exact computation, scope as stated (mod p or finite);
**[H]** heuristic / conjecture; **[X]** route shown to fail.

---

## STATUS (rewritten after every result)

Phase: framework written; equivariant sector computation in progress.

Headline findings so far:
- (framework) The book's elliptic companions (Ch. 3A.2/3A.5) **are** the abelian case of the requested
  equivariant transfer: `W_f = (T_f x E)/mu3`, a *twisted quotient* `(variety x E_chi)/K` with `K` cyclic and
  `E_chi` an elliptic curve carrying the character. The "permutation" program is the nonabelian generalization
  (`K = S3`, `S3xS3`, the block-swap `Z/2`, the Prym `Z/3`). The averaging-ideal/Maschke criterion of Ch. 3.1
  reads, for a twisted quotient `(X x Y)/K`: the K-invariants of `H(X) (x) H(Y)` contain the rho-isotypic part
  of `H(X)` whenever `H(Y)` contains `rho*`. So reaching a target sector `rho` needs (a) an auxiliary `Y` with
  `rho* ⊂ H(Y)` carrying the K-action, and (b) `(X x Y)/K` uniruled / Fano / surface-fibered with pg=q=0 (so the
  support theorem applies). (a) is always arrangeable; (b) is the real obstruction.
- (to be filled by computation) sector Hodge numbers and endomorphism/CM status per lane.

(Sections below are appended/updated as results land.)

---

## 1. The equivariant-transfer framework (precise)

### 1.1 Twisted quotients and the Maschke criterion with twisted coefficients [P]

Let a finite group `G` act on the smooth sextic fourfold `X`, `K ≤ G` a subgroup, and `Y` a smooth projective
variety with a `K`-action (the "character carrier"). On the diagonal quotient `(X x Y)/K`, rational cohomology
is the K-invariants, and by Künneth + Maschke
```
H^m((X x Y)/K, Q) = (H^*(X) (x) H^*(Y))^K = ⊕_{a+b=m} ⊕_rho  [H^a(X)_rho (x) H^b(Y)_{rho*}],
```
sum over irreducible Q-reps `rho` of `K`, where `V_rho` is the rho-isotypic projector image. Hence the
K-invariant cohomology of the twisted product *sees* the rho-isotypic part of `H(X)` exactly for those `rho`
occurring in `H(Y)`. This is the twisted-coefficient form of the averaging ideal (Ch. 3.1): an *honest* quotient
`X/K` (take `Y = point`, only `rho = triv`) serves `rho` iff `rho^K ≠ 0`; a twisted quotient with carrier `Y`
serves `rho` iff `rho* ⊂ H(Y)`. The open "residual" lanes are precisely the sectors with `rho^K = 0`, which is
why no honest quotient reaches them — and why a nontrivial carrier `Y` is the only hope.

The transfer is realized by algebraic correspondences exactly as in Ch. 3.1/3A.3: resolve `(X x Y)/K`
equivariantly, support its relevant `H` by covering curves (if uniruled/Fano) or by the surface-fibration
theorem (if it has a pg=q=0 surface fibration), then pull/push and (for the companion case) contract the
`Y`-factor against a second copy via a translated diagonal. The book's `E = Fermat cubic`/`E_23` are the
cyclic carriers; a nonabelian `rho` needs a carrier whose `H^*` contains `rho*`.

### 1.2 The two obstructions a carrier must clear [P]

**(O1) Uniruledness / Fano / fibration of the twisted quotient.** Needed so the support theorem applies.
- *Kollár–Larsen (ages).* For `k ∈ K` acting on `X x Y`, the age at a fixed point `(p,e)` is
  `age_X(k;p) + age_Y(k;e) ≥ age_X(k;p)`. So if **every** nontrivial `k` that fixes points on `X` has
  `age_X(k) ≥ 1`, then no twisted product `(X x Y)/K` is uniruled *by the age criterion* — for any carrier `Y`.
  This is a clean generalization of the per-lane age screens in the dossiers and it kills the Kollár–Larsen
  route for all four lanes at once (ages verified ≥ 1 in each dossier). It does **not** kill a *weighted-Fano*
  model (the elliptic companion `W_f` is Fano by adjunction, not by ages), so (O1) is not automatically fatal.
- *Holomorphic-top-form test.* `(X x Y)/K` (a `(4 + dim Y)`-fold) carries a holomorphic top form in the
  invariant part iff `H^{4,0}(X)_rho (x) H^{top,0}(Y)_{rho*} ≠ 0` for the twisting character. If it does, the
  smooth model has Kodaira dimension ≥ 0 and is not uniruled (this is the book's "direct boundary to the
  elliptic trick", Ch. 3A.6 last paragraph, in general form).

**(O2) Hodge type match.** Algebraic correspondences preserve Hodge type. The contraction that removes the
`Y`-factor pairs `H(Y)_{rho*}` against a second carrier; for the target's `(3,1)` piece to survive, the
carrier's relevant eigenspace must pair a `(1,0)` with a `(2,1)` (as `H^1(E)_zeta` does in 3A), never a `(2,0)`
with the target `(2,0)` (which produces a holomorphic form and blocks uniruledness). Equivalently the carrier's
`rho*`-isotypic `H^1` must be of CM/Hodge type compatible with the target's transcendental Hodge structure.

### 1.3 Reduction to the Mumford–Tate / endomorphism algebra of the sector [P/H]

Combining (O1)+(O2): the twisted transfer supports the target sector iff there is an auxiliary carrier `Y`
(uniruled or surface-fibered after the twist) whose cohomology, paired down, realizes the sector's transcendental
Hodge structure `T_sec` as a Gysin image. When `T_sec` is of **abelian/CM type** this is doable with explicit
curves/abelian carriers (the companion mechanism); when `T_sec` is of **generic/large Mumford–Tate type** (e.g.
`so(T)` for `End_Hdg(T)=Q`), only a Kuga–Satake-type carrier of astronomical dimension can carry `rho*`, and no
uniruled/low-degree model exists. So the decisive per-lane question is the **endomorphism/MT algebra of the
sector's transcendental Hodge structure** — computed below.

(Per-lane sections 2–5 appended as computations finish.)
