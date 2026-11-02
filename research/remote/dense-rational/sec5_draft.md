## 5. Ranked next steps (compute estimates on this 40-core machine, exact unless stated)

1. **Publish the degree-7 τ-pure exclusion** (this note §2.4) in discussion #59 with the msolve inputs
   (`polar_d7_pd{0,1}_{slice,homog}_p1073741827.ms`), the generator script and the kernel-sweep table; it raises fqlx's
   rigorous bound `deg D ≥ 7` to "`deg D ≥ 8` for τ-stable supports". Cost: 0.
2. **Degree-8 polar sector, numerically** (dossier 5.3): build the float projector onto `S_16^-/J_16^-` (SVD of
   `10152 × 13104`, ~30 min), then LM on `A·w^T Hess(F) w ⊥ J_16^-` with 156/159 gauge-fixed complex unknowns and
   `10^4` starts (~2–4 h). Any residual `< 10^{-10}` is promoted exactly by `e2norm.py`. Expected: none.
3. **Degree-10 E1 kernel: structure and E2** (§2.7 started): classify the 294 (+261) kernel directions — which have small
   support, which are products `ℓ·(degree 9)` or `(quadric)·(degree 8)` (intersect `K_full(10)` with `ℓ·S_9^±`: cheap rank
   tests) — and run `e2norm.py` (26532 × 9261, ~20 min each) on the structured ones. This is the first degree at which a
   non-polar denominator can exist; a structured kernel element (e.g. one vanishing doubly along a surface of `X`) would be
   the first genuinely new candidate.
4. **Full-E2 variable projection at d = 7..9** (dossier 5.4): outer unknowns `w` (60/156/350), inner dense least squares
   (`≤ 16797 × 4970`). ~1 s per inner solve → `10^3` starts per hour. Subsumes the Hessian relaxation and is the only
   exact-structure search available at d = 9 (where the relaxation is underdetermined: 208 equations, 350 unknowns).
5. **Gröbner basis of `J(F)` mod p** (msolve/Singular, colength 15625, minutes) → normal forms in degrees 16–24: makes the
   degree-12 block products, `det Hess(F)` and the cubic-cover discriminant (degree 24) testable (E1 is automatic for
   `d ≥ 19`; E2 needs `K` of degree 20 — only with sparse normal-form arithmetic).
6. **Hodge-theoretic simplicity** (§1.3): decide whether `H^4(X,Q)^τ_prim` is simple for this `F`. Evidence: served-target
   dimensions of every tested `D` are 18 or 226 (§2.8 + dossier 3.8). A proof (e.g. via the full monodromy of the
   even-sextic family plus a specialisation argument) would turn the budget `Σ h^{2,0}(Ỹ_i) ≥ 208` into a theorem and
   kill all "small" supports at once; it would also show that any solve must be a *general* mechanism (Ch. 4 §4.4).
7. **Deformation to all orders** (§2.6): if first order is consistent for some `I`, continue to second order (bilinear
   in the first-order kernel) to locate the obstruction order; if inconsistent at order 1 for all 12 seeds, record it as
   the formal obstruction of the Fermat degree-6 component (cheap, < 1 h).
8. **Lean**: nothing to formalise without a witness. If a witness appears at degree ≤ 8, coprimality follows the
   `two-ternary-dense` specialisation pattern and the identity is within the `ring`/`decide` budget only with structure.
