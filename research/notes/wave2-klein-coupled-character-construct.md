# Wave 2 — klein-coupled-character — CONSTRUCTION lens (index)

STATUS: NO WITNESS; no construction mechanism found; in progress (updated 2026-11-01 ~23:40 local; char-0 replays finishing).

Log: /tmp/hodge/results/klein-coupled-character/wave2-construct.log.jsonl ; scripts/outputs: /tmp/hodge/results/klein-coupled-character/w2c/

## Parts
- part1: framework — representation-theoretic form of the book-4.6 deck+exchange criterion; a 4.6 cover is a twisted
  fibre product X x_{X/g, g-bar} X with a NON-LIFTING involution g-bar of X/g; fibrewise (dihedral) form as in holdout.
  /tmp/hodge/notes/wave2-klein-coupled-character-construct-part1.md
- part2: linear exchanges (book 2.7) for all linear forms/complements: excluded via the degenerate-polar locus (DP);
  Aut of the geometric generic fibres of every Y-ratio fibration = Z/14 (surface fibre exact; threefolds nearly exact);
  Aut(X) = Z/14 via the 21 triple points (e_x,e_y,e_z + 18 Fermat hyperflexes), exact in charts x,t.
  /tmp/hodge/notes/wave2-klein-coupled-character-construct-part2.md
- part3: quotient companions (K-L ages: only W is uniruled, via the junior element sigma^3 tau), Galois closures of
  projections (general-type intermediate quotients), deformation (no t=0 certificate), and the best candidate: the
  Geiser-twisted fibre product Z_gamma = X x_{W,gamma} X = Kummer cover w^14 = phi^5 x^3 y^9 (explicit, non-lifting proved,
  group test passes) which FAILS the support hypothesis (h^{4,0}(Z/<swap>) >= 1 rigorous; <r,swap>: ~1000 invariant forms).
  /tmp/hodge/notes/wave2-klein-coupled-character-construct-part3.md

## Summary (20 lines)
1. Target sector: zeta^3-eigenspace of H^4_prim(X)^tau, Hodge (26,123,37); every quotient of X averages it to 0.
2. The only character-changing mechanism in the book (4.6: deck d + exchange s' with s'ds' covering g) is equivalent to a
   non-lifting birational involution of X/<g>, g in {sigma^k, sigma^k tau}; the exchange must NOT descend to X.
3. Book 2.7 linear exchange (any linear form, any complement): needs a point p with Hess_F(p)|ker grad F(p) = 0, F(p) != 0;
   only p = e_t, e_r, e_s (mod p; exact replay running), and there the coupling form M (rank 3 over k(X)) blocks the
   required splitting. Rigorous exclusion of all coordinate pairings and non-coordinate linear forms.
4. Base-change (holdout-type) covers over the Y-ratio fibrations need extra automorphisms of the geometric generic fibre:
   surface fibre S_Y: Aut = Z/14 exactly (triple-point locus = {e_x,e_y,e_z} over Qbar + monomial classification +
   Matsusaka-Mumford); threefold fibres: same (charts x,y,z exact, Fermat-point charts mod p / running). Closed.
5. Aut(X) = Z/14: the triple points of X are the 3 coordinate points (multiplicity 3) and 18 Fermat hyperflexes
   (multiplicity 4); any automorphism is block-monomial, and the monomial group is Z/14. Closes dossier 1.5's gap.
6. Dossier correction: sigma^3 tau is junior (age 1/2 at e_x,e_y,e_z); that is why W = X/<sigma,tau> is Fano.
7. Quotients by subgroups of Z/14: X/sigma, X/tau have kappa = 0 (no Fano/surface-fibration companion); W serves weight 0 only.
8. Galois closures of projections from points of X: intermediate quotients of general type; constituent test fails.
9. Deformation/q-raising: no certificate at the split limit for x^3 t^3, support not deformation-invariant: nothing.
10. Best candidate: Geiser involution gamma of the cubic-surface fibre of W from e_U (explicit formula), non-lifting
    (Kummer-class divisor argument), twisted fibre product Z_gamma = Kummer cover w^14 = phi^5 x^3 y^9 of X.
11. It fails exactly at "Z_gamma/K supported": h^{4,0}(Z_gamma/<swap>) >= 1 (rigorous); with r = tau x tau the
    eigensheaf count still gives ~10^3 invariant (4,0)-forms (not a Fano fibration; no base curve).
12. Hodge facts: zeta^3-part of H^{2,0}(S_Y) = span(xy) (surface-fibration theorem inapplicable to the sector);
    H^{3,0}(V)[zeta^3] = 0 for the threefold fibres (sector of type (2,1)+(1,2)).
13. No finite-field incidence seed can be produced: no auxiliary uniruled quotient reaches weight 3.
14. Remaining unexplored: covers not of base-change or twisted-fibre-product type; Bertini-type involutions of W
    (same Kummer-cover obstruction expected); non-eigen polynomial Piola maps (dossier item 7).
