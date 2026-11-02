STATUS: NO WITNESS, NO NEW SERVING MECHANISM FOUND (2026-11-02 ~04:00 local). Rigorous negatives (Reid–Tai, no-hub, Galois closure), one new sufficient criterion (descent through a reflection), sharpened structure (Klein cover, Fano elliptic fibrations, rational genus-2 surface, p_g of base changes); served-space scan: 4/7 structurally motivated denominators done, all trivial as predicted (two pass E1); degree-12 ones running in sector (3,0) (jobs construct-off-paired-served, -served30); Mordell–Weil screen: no polynomial section of Y± (3 primes / 1 prime).

# Wave 2 — off-paired — CONSTRUCTION lens (index)

Parts (each < 300 lines):
- part1: /tmp/hodge/notes/wave2-off-paired-construct-part1.md — target, structure summary S1–S5, mechanism inventory (i)–(v) verdicts
- part2: /tmp/hodge/notes/wave2-off-paired-construct-part2.md — Klein cover / elliptic fibrations / rational genus-2 surface Σ and its base changes (p_g table), Leray analysis of divisor supports (V1/H1/M), Reid–Tai theorem
- part3: /tmp/hodge/notes/wave2-off-paired-construct-part3.md — descent-through-reflection lemma (proof), fibre-product realisation and its p_g obstruction, no-hub theorem, Galois-closure verdict
- part4: /tmp/hodge/notes/wave2-off-paired-construct-part4.md — computation table, served-space scan results, verdicts, decisive next computations
Scripts/outputs: /tmp/hodge/results/off-paired/wave2-construct/ ; log: /tmp/hodge/results/off-paired/wave2-construct.log.jsonl

## 20-line summary
1. F = H(x) + C(x) X2Y2 + X2^6 + Y2^6 (x = X0,X1,Y0,Y1); the target carries the determinant character of the dihedral group D6 = <h, sigma2> on the (X2,Y2)-block.
2. X/D6 = P(1,1,1,1,2); Ŷ = X/<h^2> is the Klein (Z/2)^2-cover t1^2 = H+Cm−2m^3, t2^2 = H+Cm+2m^3 (m = X2Y2); the target sector is its (−,−) character.
3. The two other characters are H^4 of the FANO fourfolds Y_± = X/<h^2, reflection> = Weierstrass elliptic fibrations t^2 = ∓2m^3 − Cm − H over P^3 (disc ∝ 2C^3 ∓ 27H^2); the (−,−) part is the Prym J(E_b) of the genus-2 fibration X/<h> → P^3.
4. The h-invariant Y-odd primitive sector is (56,144,56), all sigma2-anti-invariant; g-characters 0/2/4 give (22,50,22)/(17,47,17)/(17,47,17).
5. The genus-2 family depends only on kappa = C^3/H^2; its total space Σ ⊂ P(1,2,3,3) is RATIONAL (2A1+2A2, K = O(−3)); explicit parametrisation, fibre = (3,2)-curve (1+w)(1−s^2) = 2 lambda w^3(1+s^2).
6. Hence the genus-2 variation carries no transcendental (2,0): every (3,1)-class of the sector involves H^{2,0} of the pencil surfaces S̃_kappa = {C = k s^2, H = k s^3} ⊂ P(1,1,1,1,2) (quasi-smooth, p_g = 45, q = 0).
7. Base changes kappa = rho^n have p_g(Σ̃_n) = 0,0,0,1,2,1 (n = 1..6; K3 at n = 4, 6), proved via double-plane models; twisted (2,0)-classes need characters of order ≥ 4.
8. Leray analysis: vertical denominators (pulled back from P^3) serve only V ⊗ (algebraic classes of S_kappa twisted by the component permutation); kappa-constant ones (H, C, Δ) serve nothing (= the thread's obstructions); kappa-root ones of the shape f(C,H) = mu^n have linearly equivalent conjugate components, hence serve nothing; the only loophole is a Picard local system of the pencil with monodromy.
9. Horizontal denominators pulled back from multisections of Σ (block forms, mirrors, torsion multisections) need flat twisted (2,0)-forms of the pencil: none expected, consistent with every scan so far.
10. GHC therefore demands a MIXED divisor: fibrewise a correspondence coupling J(E_kappa) to H^2_tr(S_kappa) — the kind of object the holdout cover produces on the same block family.
11. Reid–Tai: all G-elements have age ≥ 1, no pseudo-reflections, det = ±1 ⇒ every X/K (K ≤ G) has canonical singularities and P_2 ≥ 1: NOT uniruled (also for elliptic/K3-twisted companions). Item (ii) fails rigorously.
12. No-hub theorem: Q_X(Y) has rank 3, so no linear hub exists; each coordinate hub is killed by exactly one mixed monomial; the symmetrised equation does not preserve the cover relation. Item (i) fails.
13. Galois closure X → X/D6: the det character is served only by rotation quotients, all with p_g = 1. Item (iii) fails; the cubic cover is book §4.7.2.
14. Deformation/q-raising: no certificate on the chain sub-family to continue; nothing new. Item (iv) fails; (v) N/A.
15. NEW criterion (part 3 §E): a cover Z → X with deck e and involution sigma such that sigma e sigma lifts a reflection, plus Z/<sigma,R> supported ⇒ target supported; no supported quotient of X needed.
16. Realisation Z = X ×_W X^phi (W = X/sigma2 Fano) or over P(1,1,1,1,2); automorphism phi's are reducible or give p_g ≥ 1; a birational phi is required; the natural candidates (reflections through Mordell–Weil sections of Y_±) have no polynomial instance: msolve finds no section (m quadric, t cubic) mod three primes (Y_-) / one prime (Y_+) [SCREEN].
17. Served-space scan (remote): X2^6−2Y2^6, C^2−(X0^2+Y1^2)^4, C^2−(X0Y1+X1Y0)^4 serve only the Jacobian part in all 18 Y-odd sectors (the two vertical ones pass E1) — as predicted; degree-12 candidates pending (part 4 §I). Hodge-bundle table for base changes n ≤ 12 in part 2 §B'.
18. No finite-field incidence seed: no auxiliary variety with the required support exists yet.
19. Decisive next steps: Picard number/monodromy of the pencil S̃_kappa (van Luijk at a few kappa); dossier E5; search for birational involutions of W or Y_- (Mordell–Weil).
20. All exact checks reproducible with the scripts listed in part 4 §H (sympy + Singular, seconds each).
