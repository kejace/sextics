STATUS: NO MECHANISM FOUND — construction lens closed with rigorous coordinate-free exclusions (V1 = V2 = ∅ => no exchange cover, no quasi-reflection in any linear coordinates; V3 pending => no linear reflection), E1 exclusion of all ternary cones through degree 8 (9..12 pending on the remote), and a structural theorem (equivariance reduction to W = X/tau + Lefschetz) showing every universal divisor of the cubic-cover geometry (R, R', Δ~, det Hess, bihomogeneous sections, cones, products) has zero usable (2,0)-budget. Last update: see part 3 header. No witness. Nothing posted.

# Wave 2 — dense-rational — CONSTRUCTION lens (index)

Parts (each < 300 lines):
* part1 — landscape, rigid structure of X/tau = W (cubic cover of P^2 x P^2), triage of mechanisms (i)-(v)
  /tmp/hodge/notes/wave2-dense-rational-construct-part1.md
* part2 — coordinate-free exclusion of exchange covers and reflections (V1, V2, V3 definitions and proofs), cone test method
  /tmp/hodge/notes/wave2-dense-rational-construct-part2.md
* part3 — computational results table (msolve certificates, cone E1 tests, universal degree-12 E1 tests)
  /tmp/hodge/notes/wave2-dense-rational-construct-part3.md
* part4 — what a mechanism must produce ((2,0)-budget refinements), audit of universal divisors, R' with explicit equations,
  equivariance reduction to W
  /tmp/hodge/notes/wave2-dense-rational-construct-part4.md
Scripts: /tmp/hodge/results/dense-rational/wave2_construct_{hess,ltcheck,cones,univ}.py ; log
/tmp/hodge/results/dense-rational/wave2-construct.log.jsonl ; packet dir /tmp/hodge/results/dense-rational/wave2-construct/.

## 20-line summary
1. Target: q = 1, A = X0^3 Y0^3 on the 236-term even sextic; sector H^4(X)^tau_prim = (0,208,866,208,0), expected simple.
2. A certificate = a divisor Y = X ∩ V(D) whose resolution carries >= 208 holomorphic 2-forms with Gysin image ⊇ sector.
3. (i) Exchange covers (book §2.7, any coordinates/blocks/base form) force a point e ∉ X with rank Hess F(e) = 1.
4. V1 = {rank Hess <= 1} = ∅ and V2 = {rank Hess <= 2} = ∅, certified by msolve leading ideals (pure powers) => rigorous over C.
5. Hence no exchange cover and no quasi-reflection (orders 3,5,6) exists in ANY linear coordinates; Fermat sanity checks pass.
6. Linear reflections: F even in a linear form forces the degree-15 system V3 (cubic Taylor term vanishing on the polar
   hyperplane); msolve pending; Fermat sanity non-empty as required.
7. (ii) Aut_proj(F) = {1, tau} (diagonal/permutation exact; infinitesimal exact; reflections/quasi-reflections now excluded
   coordinate-free); X/tau not uniruled; twisted companions (X x Y)/<tau x iota> never uniruled (ages >= 3/2, top form).
8. Surface fibrations of W: fibres p_g = 3 (over P^2_X) or 20 (Ch. 4.7): theorem inapplicable.
9. (iii) Galois closure of the cubic cover: V_std^K = 2,1,0,0; the only quotients carrying the sector are W and Z~, not uniruled.
10. (iv) Deformation from Fermat is obstructed at fixed degree (modular, 4/12 seeds); no sector theorem on a sub-family transfers.
11. (v) No solved companion on this variety.
12. Ternary cones g(X), g(Y): zero budget by Lefschetz monodromy of the net of P^3's through P^2_Y, and E1 fails rigorously for
    k = 7, 8 (9..12 pending).
13. Equivariance: only H^2(Y~)^tau = H^2((Y/tau)~) can serve the invariant class, so the problem lives on the smooth fourfold W.
14. Universal divisors of W: R = {E_X F = 0} and the residual-tangency divisor R' = {(F24+F06)^2 + 4 F60 F06 = 0} are the two
    double covers of the (12,12) discriminant Δ; Δ is unibranch, Lefschetz gives H^2(Δ~) all (1,1): budget 0 for both.
15. det Hess F, Δ~, smooth bihomogeneous sections: budget 0 (Lefschetz on X or W); products: budgets add, nothing new.
16. Therefore every published mechanism and every universal divisor fails for a structural reason; the remaining space is
    non-universal singular divisors of degree >= 10 (E1 kernel 294/261 at d = 10), i.e. the search lens.
17. No finite-field incidence seed can be produced: there is no cover, quotient or curve family to seed.
18. Deliverables for a discussion post: V1/V2(/V3) certificates (new, coordinate-free), cone E1 table, the R' equation and the
    equivariance/Lefschetz argument (theory, with its hypotheses stated).
19. Caveats: Lefschetz-unibranch argument for Δ assumes the cusp-locus resolution contributes only exceptional classes
    (true by the pushforward argument); simplicity of the sector is an expectation, not a theorem.
20. Nothing was posted to GitHub; no git commands were run.
