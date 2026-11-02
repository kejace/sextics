# Wave 2 — lane `cyclic-triple-residual` — lens CONSTRUCT (new mechanism)

STATUS: finished 2026-11-02 (UTC) — NO WITNESS. All published quotient/cover/fibration mechanisms excluded for the ζ3-sector by
exact computations; the one character-compatible mechanism (isotrivial j=0 elliptic exceptional surfaces / sweeps of the
elliptic curves over lines of W) identified and tested at degree 9 (negative); new rigorous-over-Q sector-variation result:
very general supports need Σ h^{2,0} ≥ 67.  One remote job (Galois-closure smoothness, Singular) still queued behind the
capacity gate: construct-ctr-ztsmooth; output will appear on the remote in wave2-construct/ and /tmp/hodge/jobs/construct-ctr-ztsmooth.log
(pull with `python3 -m hodge.rjob pull --lane cyclic-triple-residual`).

Work dir: /tmp/hodge/results/cyclic-triple-residual/wave2-construct/ ; log: /tmp/hodge/results/cyclic-triple-residual/wave2-construct.log.jsonl
(steps resume, C1, C3, C2, notes, C4-dims, C4-test ×4, C4-summary, C5, C5b).  Engine v2 used only for the sector solver (sect.test) in C4.

## Parts
* part1 — covers/Galois closure Z̃ (order-108 symmetry group, all 128 subgroups): 38 serve the class irrep, none with negative canonical degree; abelian-extension lemma; §4.6 descent status; P^2×P^2 closure.
* part2 — §2.7 exchange covers for every linear form (Singular char 0 / msolve): only X1, Y0 directions; both blocked by two couplings; cubic variant forces a CY base.
* part3 — surface-fibration companions of X/K (exact p_g counts, 20 + 45 bases): min p_g = 2; K3 fibration of W remark.
* part4 — surviving mechanism: triple-point curves of T_W ∩ B ⇒ isotrivial j=0 elliptic exceptional surfaces with character ζ3; degree-9 families (120 nonempty, dims 29–66) and 4 mod-p E2 tests (all inconsistent).
* part5 — E5 refinement: supports = π^{-1}(scrolls of lines over irregular surfaces of 𝓕); Lean-route assessment; recommended next computation (compute 𝓕).
* part6 — sector variation over the G-family: Def·ξ = 41 of 119, Def²·ξ = 26 = all, K = K₂ = 0 ⇒ smallest flat sub-VHS ∋ ξ is all of S_A ⇒ Σ h^{2,0} ≥ 67 for very general supports (dossier E4 settled in computable form).

## Summary
1. Target ξ = [X0^3Y0^3 Ω/F^2]: μ3-character ζ3, τ-invariant, exact G-sector (41,119,26); the dossier's V_i are τ-invariant parts (the full ζ3-eigenspace also holds Ω/F).
2. Natural Galois closure Z̃ = X ×_{P^4} Y, Y = {u^3 = H} Fano; full group G108; Riemann–Hurwitz degree formula deg K_{Z̃/K} = 2 − 2[d∈K] − 2[e∈K] − #σ-conjugates in K reproduces X, W, Y, P^4.
3. Averaging-ideal test over all 128 subgroups (exact): every quotient serving ξ has deg K ≥ 0 (X, X', X/τ, X/⟨s⟩…); every uniruled quotient (P^4, W, Y, …) misses ξ. Twisted elements (g3·e^{-1} etc.) add nothing.
4. Lemma: no tower of cyclic covers over the natural structure can change this (abelian deck groups act trivially on pulled-back classes).
5. §4.6 descent would need Z̃/⟨r,σ⟩ (deg K = 1) supported — no sector theorem applies; its ρ-part contains the (4,0)-class.
6. Exchange covers: F has the y^6 + yΛ + f shape only in directions e_{X1}, e_{Y0} (Singular, char 0, 2 components); for both, the couplings X0^2Y0Y1Y2^2 / X0X1Y1^3Y2 and X0X2^3Y1Y2 break the T-block condition; the X2-cubic variant needs T^6 = H (CY base).
7. Surface fibrations of X/K, K ⊂ ⟨s⟩: fibre-quotient p_g ≥ 2 (best: base [X1:X2]×[Y0:Y1], K_b = ⟨s⟩); the (6,3) passive-ratio fibres have only τ fibrewise (p_g 20).
8. W is K3-fibred over P^2; the fibres of X are sextic surfaces with p_g = 1 + 6 + 3 by μ3-character — the twisted 9 are the fibrewise shadow of V_1 ⊕ V_2 and are never killed.
9. Only character-compatible divisor geometry: cyclic triple cover of a divisor T_W ⊂ W whose branch curve has ordinary triple points ⇒ Ẽ6 singularities with j=0 exceptional elliptic surfaces, holomorphic forms of character ζ3 (X-language: D of order 3 along planes Π ⊂ {X2=0}, leading cubic c_3 + aX2^3).
10. Degree 9: all 120 (coordinate plane, μ3-invariant sector) order-3 families nonempty (dims 29–66, mostly J-type); four random members with genuine j=0 cubics: inconsistent mod two primes (rank_aug = rank+1), same signature as all degree-8 plane families.
11. Sector variation (rigorous over Q via rank semicontinuity): Def_G·ξ spans 41 of the 119 (2,2)-classes, Def_G²·ξ all 26 (1,3)-classes, no Def-annihilated (3,1) or (2,2) classes ⇒ the smallest flat sub-VHS containing ξ is the whole sector ⇒ very general supports need ≥ 41 twisted + 26 conjugate 2-forms.
12. Hence low-degree (≤ 10) elliptic-plane denominators (≤ 10 forms per curve) cannot work generically; the lane needs either deg D ≳ 12 with ≥ 7 triple-point curves (parametric search) or the E5 geometry: π^{-1}(scroll of lines over a surface S ⊂ 𝓕 with twisted irregularity ≥ 41).
13. No Lean route in the shared library fits: explicit witness (none), quotient+integration (no serving uniruled quotient), §4.6 descent (no exchange), surface-fibration (p_g ≥ 2); elliptic-fibre integration gives elliptic integrals, outside the polynomial-source theorem.
14. Recommended next computation: the Fano variety of lines 𝓕 on W (incidence Gröbner over F_p, remote) and the twisted irregularity of the surface of lines tangent to B; then the sweep divisor's equation and a sparse sector test.
15. Pending: Z̃ smoothness/quasi-smoothness of Z̃/σ, Y, B (Singular job queued; needed only to make the klt/BDPP step of item 3 fully rigorous — the degree bookkeeping is unaffected).
