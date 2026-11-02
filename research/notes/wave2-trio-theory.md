STATUS: IN PROGRESS (2026-11-02). Parts 1-2 written (EV transfer theory, W-tower rigidity theorem, Galois closures incl. NEW qcr towers). Part 3 (sector-variation budgets rcr/qcr, klt check, verdicts) pending remote jobs trio-rcr-sectorvhs / trio-qcr-ztsing (remote saturated, waiting for the capacity gate). No witness, no new certificate on any lane.

# Wave-2: residual-companion trio — character-transfer theory (index)

Lanes: `cyclic-triple-residual` (ctr), `reflection-coupled-residual` (rcr), `quadric-contact-residual` (qcr).
Work dir: /tmp/hodge/results/_cross/ (log trio.log.jsonl, scripts trio_*.py, outputs trio-*.json).
Parts (each < 300 lines):
* part1 — Esnault–Viehweg eigen-decomposition for the three cyclic covers X → W; what a divisor support must carry (twisted (2,0)-forms); the W-tower rigidity theorem (covers Galois over W, fibre products, changed branch sections: all excluded by character theory).
* part2 — Galois closures of the natural finite maps to P^4 and the averaging-ideal test: ctr (order 108, done by construct lens), rcr (B3 × Z/2, done by construct lens), qcr (NEW: the eight reflection towers).
* part3 — decisive computations per lane, receipts, honest verdicts.

## 20-line summary
1. All three lanes are the χ ≠ 1 eigen-sector of a cyclic cover π: X → W = X/G onto a Fano quotient (ctr: μ3, W ⊂ P(1^5,3); rcr: Y2 ↦ −Y2, W ⊂ P(1^5,2); qcr: X0 ↔ Y0, W ⊂ P(1^5,2)); the companion serves H^4(W) = H^4(X)_0.
2. Esnault–Viehweg: H^{p,q}(X)_χ = H^q(W, Ω^p_W(log B) ⊗ O_W(−1)), B the branch threefold (a hyperplane section of X, smooth: Singular check); the (4,0)-form always has character χ; the lane sectors are the isotypic parts under the remaining symmetries.
3. A support must be G-invariant, E = π^{-1}(E_W), with twisted (2,0)-forms on the equivariant resolution: H^0(Ẽ_W, Ω²(log B̃) ⊗ L̃^{-1}) ≠ 0. Smooth or rational-singularity E_W never has them (Bott/Lefschetz vanishing); only cover ramification along singular loci of B ∩ E_W creates them — the hyperelliptic-slice (n = 2) and j = 0 elliptic-surface (n = 3) mechanisms are the two local models.
4. THEOREM (W-tower rigidity): for any cover Ẑ → W factoring through X, with M = Gal(Ẑ/W) ⊃ N = Gal(Ẑ/X), the averaging projector e_K kills the pulled-back class unless K ⊂ N. Hence fibre products X ×_W X' (any second branch section, any twist), diagonal covers, and all abelian towers over W see the χ-sector only through covers of X: "changing the branch section" and "a second auxiliary cover over W" transfer nothing. Only Galois closures over a different base (K ⊄ M) can.
5. Galois closures: uniform Riemann–Hurwitz formula deg K_{Z̃/K} = 4 − #(B3-reflections in K) for the (2,4,6) complete-intersection closure Z̃ ⊂ P^7 (rcr, qcr); the reaching quotients are X/H (κ = 0), the (4,6)-CI in P(1^6,2) (K = 2H) and Z̃/S3 (K = H) — never uniruled (klt + BDPP; klt check of Z̃ submitted). ctr: construct lens's 128-subgroup enumeration, explained by the theorem.
6. NEW (qcr): all four real reflections of Γ⁻ give towers with full B3 closure (sympy Galois groups at 3 random specialisations each); the lane is structurally four copies of the rcr situation, all negative.
7. NEW (qcr budget): over the 26-parameter Γ-symmetric family the Q-flat closure of ξ is the entire ε-sector (22,51,22) (Def²·ξ fills (1,3); Def·S_A^{3,1} fills (2,2); kernel 0; two primes) ⇒ any support moving with F carries ≥ 22 twisted (2,0)-forms; structured low-degree denominators cannot work generically.
8. rcr budget: same test, sector (126,492,126) — RESULT_RCR_PLACEHOLDER
9. ctr: budget 67 (construct lens), all cover mechanisms closed; decisive next computation = Fano variety of lines of W (hours).
10. Verdicts: no new route to a certificate on any of the three lanes; the cover/correspondence family of ideas is now closed by theorem rather than by enumeration; what remains is explicit twisted-form divisor geometry (slices / elliptic surfaces) with the budgets 67 / 22 / (rcr) as the minimum complexity of any generic support.
