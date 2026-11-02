## 0. Corrections to the wave-1 record

* `scan_served.py` (wave 1, last 260 log lines, family `quot-curve6X`): the lines "A ENTERS THE SPAN" and
  "A-system consistent=True" are artefacts. The candidate `1EX(F)+1EY(F)` is `6F` (Euler), so its served
  space is everything and the cumulative span is polluted from that point on; and `S.solve()` was called on
  the system that still contained the target columns, so "consistent" is vacuous. Real content of that
  scan: every tested D of degree 3-6 serves nothing beyond `J(F)_6^G`.

## 1. Plan (wave 2)

1. E1-admissible denominator spaces per degree d = 4..14 in the character sectors triv/sgnX/sgnY/sgnsgn
   of S3xS3 (both Y-parities), plus S2xS2 / S3xS2 / S2xS3 types at d = 10 and the unrestricted kernel at d = 9, 10.
2. Single-D exact (mod p, CRT-lifted) solves of every kernel class + random combinations (never done for
   the non-invariant classes).
3. Character-reduced variable-projection (varpro) joint search over the full E1-admissible family
   (kernel + J^chi_d) per (chi, d, pD), thousands of restarts on the remote.
4. Structured families (cones over plane curves with the three G-stable vertex planes, line-vertex cones,
   discriminant/Hessian-block products, products of kernel elements) via restricted E1 kernels + solves.
5. Rational reconstruction + exact re-solve of every candidate with residual < 1e-10.

(Results sections are appended below as batches finish.)

## 2. Methods built for wave 2 (all in the lane dir; engine files untouched)

* `w2lib.CharFamily` — character-reduced *float* E-system for a linear family `D(a) = Σ a_k B_k` of
  χ-semi-invariant denominators (same orbit folding and column meaning as `charsolve.CharSystem`), with the
  variable-projection objective `f(a) = min_x |M(a)x − b(a)|²` on the unit sphere of the admissible subspace
  (E1 kernel inside the family, computed mod p for the dimension and by SVD for the float basis, minus the
  F-multiples). Gradient by the Golub–Pereyra formula; L-BFGS-B restarts; rational reconstruction with
  increasing denominators; exact confirmation by `charsolve.solve_char` (CRT lift) + `core.verify_witness`.
  Reduction factor ≈ 36 (S3×S3): d = 10 invariant system is 2.6k × 0.9k instead of 85k × 19k.
* `w2lib.JacFamily` — the Jacobian-type ("w-space") search of the task, *unrestricted by symmetry*: `D = Σ w_i F_i`
  (+ the E1-kernel classes beyond J with least-squares E1 lifts), explicit lift `N₀ = A·w/q`, unknowns = Koszul
  corrections `c_ij` (degree d−4) and `Q`; equation `T_D(c) − F·Q = −B_D(N₀)` in degree 2d (book ch. 5, full
  correction criterion; complete at q = 1 because the `F·H` freedom in E1 changes `B_D` by `−HD/q`). Sizes:
  d = 6: 3080 × 388 (K = 18), d = 7: ≈ 5.8k × 0.8k (K = 64), d = 8: ≈ 10k × 1.9k (K = 174).
* Validation (logged as `pipeline-validation`): on the Fermat sextic `p6X + p6Y` both engines recover the
  subset-Euler denominator at d = 6 (residual 1e-24) and the exact solvers lift and verify it; gradient
  checks agree to 8 digits; on the real F the invariant `J^G_6` family (dimP = 3) has a single minimum at
  rel 8.3e-3 and the unrestricted w-space at d = 6 has minima at rel ≈ 0.2.

## 3. Structural observation (reduces the semi-invariant sectors to the invariant one)

For the S3 action on three variables every sign-semi-invariant polynomial is `Δ × (symmetric)`, so
`sgnX`-type D = `Δ_X·S`, `sgnY`-type = `Δ_Y·S`, `sgnsgn`-type = `Δ_XΔ_Y·S` with `S` invariant (the orbit
counts confirm it: the sgnsgn family at (d,pD) = (10,1) has `nD = 12` = number of invariant Y-even quartic
orbits). With the sum principle `served(D₁D₂) = served(D₁) + served(D₂)` (coprime factors) and the fact that
the reflection `s = (X0 X1)` fixes the section `X ∩ {X0 = X1}` pointwise — so every class supported on it is
`s`-invariant, i.e. `served(X_i − X_j) ⊂ U = Σ_s H⁴(X)^s`, and `U` meets the (sgn,sgn) isotypic component
trivially — one gets, for the target sector,

> `[A] ∈ served(Δ_X^a Δ_Y^b S)  ⇔  [A] ∈ served(S)`.

Hence the sign sectors are redundant with the invariant sector at degree d−3 / d−6, and the wave-1 exclusions
of `D7 = Δ_X e2(X²)` and `D7Y` are *implied* by the degree-4 E1 failure of `e2(X²)`, `e2(Y²)`. Compute is
therefore concentrated on `triv` (both parities), the unrestricted w-space at d ≤ 8, and the subgroup
(S2×S2, S3×S2, S2×S3) types, which probe invariant denominators of much higher degree (orbit products).
The sign sectors are still run at low cost as consistency checks of the principle.

## 6. Concurrent agent on this lane (cross-check, not our work)

Another process submits jobs named `search-symmetric-sextic-residual-{dims,dims2,struct,vp-d8triv0,vp-d8sgnX1}` and
`excl-ssr-relax1` (its scripts `w2x_*.py`, `w2_test2.py`, outputs `w2struct.jsonl`, `w2vp_*.jsonl`, `wave2-excl/` live
only on the remote lane dir). From its outputs (pulled 19:50 EST): its own character-reduced varpro ("cvp") on
`triv d=8 pD=0` (6 restarts) reached rel 8.7e-5 and on `sgnX d=8 pD=1` (22 restarts) rel 0.094 — consistent with
our landscape (invariant minima ≈ 1e-4 at d = 7–8, no zero); its structured E1 scan tested 1014 products of
symmetric functions at d = 6–8 (26 pass E1, none full-solved in that file). Our jobs use the prefix
`search-symmetric-sextic-residual-w2b-*` and the `w2b/` directory to avoid name collisions.
