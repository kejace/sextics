STATUS: DONE 2026-11-02 (second attempt; first attempt killed for output size). NO WITNESS, NO SEED. Main results: (1) core-reduction theorem — any divisor support of the target (h-invariant WLOG) forces, in the primary Leray case, a surface in the Calabi–Yau threefold T_H: u^3 = H supporting the (2,1)-class [a]; hence existence <=> a core-only denominator works, and the binary block (Prym, K3 companion W_B, Kuga–Satake) is irrelevant; (2) every twisted-quotient/companion mechanism fails by a carrier dichotomy (curve carriers: holomorphic-top-form test on the binary side; surface carriers: contraction dimension count); (3) T_H is fibred over P^1 x P^1 by C_B-type genus-4 curves with a dominant Prym period map (exact); (4) Abel–Jacobi triviality theorem: all lifted hyperplane-section families on the contact surfaces {H = c^3} (the 4-variable form of the book's perfect-power contact / Nambu recipe) have zero cylinder map, so the explicit sweep denominators D = lam_2(l)H - lam_1(l)c^3 fail at every degree.

# Wave 2 — four-plus-two-residual — CONSTRUCTION lens (index)

## Table of contents
- part1 /tmp/hodge/notes/wave2-four-plus-two-residual-construct-part1.md — sector bookkeeping; pencil structure Xt -> P^1 with R^3 = (+) V_k (x) L^k; core-reduction theorem (primary/secondary cases, proof sketch); consequences 1–4; what any mechanism must produce.
- part2 /tmp/hodge/notes/wave2-four-plus-two-residual-construct-part2.md — mechanisms (i) covers/exchange, (ii) Fano/surface-fibration companions of quotients and twisted carriers (dichotomy), (iii) Galois closures + averaging ideal, (iv) q-raising/deformation, (v) n/a — each with its precise failure point.
- part3 /tmp/hodge/notes/wave2-four-plus-two-residual-construct-part3.md — fibration of T_H by C_B-type curves (verified), Prym period-map rank 2 (verified), contact surfaces Sigma_c, Abel–Jacobi triviality theorem, explicit sweep family D = lam_2 H - lam_1 c_hat^3 (degrees 6+r+3e) and its exact failure, remaining routes (Noether–Lefschetz loci of {V_c}, non-2-normal contact space curves, simplicity of V, secondary normal-function route), Lean/seed status.
- script /tmp/hodge/work/construct42/fibration_check.py; log /tmp/hodge/results/four-plus-two-residual/wave2-construct.log.jsonl

## Summary (20 lines)
1. Target sector: h -> zeta3, r2, r5, s -> -1; Hodge numbers (32,80,32) = [H^3(T_H)_zeta^{s=-1} (type (0,32,16,0)) (x) <dt/w> ⊂ H^1(Prym)].
2. Pencil from the plane {X2=Y2=0}: fibres Y_t ≅ Y_H = {U^6 + H = 0}; R^3 = (+)_k H^3(Y_H)_{eta^k} (x) L^k with L the rank-1 system of B^{1/6}.
3. The target is the pure tensor [aU] (x) dt/w in the Leray piece H^1(P^1, R^3); joining classes (cones over S_H) are vertical and h-invariant.
4. Theorem: a support divisor Z (h-invariant WLOG) has, in the primary case, fibre surfaces Z_t ⊂ Y_H whose Gysin image contains [aU];
   a single fibre V_{t0} ⊂ P^3 then gives a core-only certificate of degree <= deg D (converse = attack1 §4.2).
5. Secondary case (Ext^1 normal-function map from monodromy-fixed eta^2-classes of the fibres) is impossible for core-only D and would need >= 32 fixed 2-forms.
6. Hence the binary side (Prym simplicity, W_B K3, Kuga–Satake, carriers for C_B, deformations of B) cannot decide existence.
7. Existence <=> generalized Hodge conjecture (coniveau 1) for the level-1 piece V of the CY threefold T_H; if V is simple every support surface has q >= 48.
8. (i) no (6,1)-variable in F or H: no §2.7/§2.8 exchange, for coordinates or linear forms. 
9. (ii) honest quotients by K0 have kappa = 0; curve carriers fail the top-form test on the binary side; surface carriers (cyclic cubic surface, both companions Fano) fail the contraction dimension count.
10. (iii) averaging ideals of all intermediate quotients miss chi_A; (iv) deformation to the split limit loses the CM piece (period map dominant); (v) n/a.
11. T_H -> P^1 x P^1 by genus-4 curves w^3 = b_l(la^2, mu^2) (exact formula), Prym family 2-dimensional (Jacobian rank 2 at random points).
12. Sweeps by these fibres only see the fixed part of the Prym family: never V.
13. Contact surfaces Sigma_c = {u = c} ∩ T_H ≅ {H = c^3} ⊂ P^3 are regular sextic surfaces (q = 0).
14. Abel–Jacobi triviality: families of curves on the Sigma_c that are multiples of the hyperplane class (all plane contact curves, all lifted hyperplane sections) have constant AJ map, zero cylinder map.
15. The explicit family D = lam_2(l1,l2) H - lam_1(l1,l2) c_hat^3 (degree 6 + r + 3e, parameter curve mu^3 = lam of genus 2r+3e-2) therefore fails at every degree; smallest members degree 7 (g=0), 8 (g=2), 9 (g=4).
16. This is the exact 4-variable failure point of the §2.3/§3A perfect-power-contact recipe: the source surface is regular, where in 3+3 it is a genus-10 curve.
17. Remaining routes: Noether–Lefschetz loci of {V_c} ((2,b)-complete-intersection contact curves, b >= 4, 40 equations/42 parameters), non-2-normal contact space curves, the simplicity of V (degree-48 Frobenius polynomial out of reach).
18. No uniruled auxiliary variety with the target character exists in the screened classes: no curve family, no incidence Jacobian, no finite-field seed.
19. Recommendation for other lenses: restrict all searches to core-only D in the reduced 5-variable model (lossless in the primary case); do not spend effort on the binary block.
20. Nothing posted; nothing in the repo modified; no remote jobs submitted (box at load ~600; all work was light local sympy).
