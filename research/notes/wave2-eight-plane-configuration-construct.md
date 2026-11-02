# Wave 2 — eight-plane-configuration — CONSTRUCTION lens

STATUS: parts 1-3 written (2026-11-01 22:45). NO WITNESS. Mechanisms (i),(iii),(iv),(v) excluded with verified
hypotheses; (ii): Aut(X) = <tau> proved within block-monomial maps (exact over Q(i)), conditional on the
Fano-scheme-of-planes check (remote job construct-epc-fanoplanes queued behind capacity; local probe of chart00).

Log: /tmp/hodge/results/eight-plane-configuration/wave2-construct.log.jsonl
Scripts/data: /tmp/hodge/results/eight-plane-configuration/wave2-construct/ (mech_i_ranks.py/.json,
quartic_cover.py/.json, fano_planes_gen.py, fano/chart*.ms)

## Parts
1. wave2-eight-plane-configuration-construct-part1.md — framework (what a construction must do), mechanism
   inventory table, Lemma A (exchange covers) + exact rank test, mechanisms (ii),(iv),(v).
2. wave2-eight-plane-configuration-construct-part2.md — the quartic plane-pair covers X -> Pi x Pi':
   explicit equations, computational checks (S_4 monodromy, reduced discriminant, non-isotrivial
   elliptic carrier), the omega-descent theorem and its proof, the precise failure point.
3. wave2-eight-plane-configuration-construct-part3.md — Aut(X) classification, the (H5) loophole analysis
   (Delta on the twisted diagonal = F60 F06 S^2; smooth/transverse branches), ramification divisors of the
   covers = X n V(sum eps_j R_j F) (excluded degree-6 Jacobian denominators), precise stopping statement.

## Summary (20 lines)
- The target is the tau-invariant (3,1) class; X is Calabi–Yau with tau^* omega = -omega, tau^* alpha = alpha.
- Every book construction = a uniruled 4-fold W with a correspondence to X seeing alpha. W must act
  differently on alpha and omega; only tau does so on X, and tau has codim-3 fixed locus (omega^2 descends).
- (i) Exchange covers (§2.7 type, any linear forms, any base degree m): Lemma A forces the 15 coefficient
  quadrics of F_{4,2} (resp. F_{2,4}) to span <= 3 dims; exact ranks are 6 and 6, all Gram ranks 3. Excluded.
- (iii) New instance-specific finite maps: the four quartic covers X -> Pi_eps x Pi_{-eps} (lines through
  antipodal planes; F(sp+tp') = st q). Verified: geometric monodromy S_4 (Frobenius types 4, 3+1, 2+1+1 at
  15 split primes; discriminant of bidegree (18,18) reduced and non-square), q symmetric on the twisted
  diagonal, non-isotrivial residual quartic. Theorem: every intermediate quotient Y/K of the Galois closure
  (group S_4 x Z/2) that carries alpha has kappa >= 0 — alpha and omega span isomorphic representations
  V_std x 1 and V_std x sgn, so a quotient seeing alpha sees the 4-form eta unless it contains a tau-type
  element, and tau-type elements have codim >= 2 fixed loci, so eta^2 descends. Same for point projections
  and the book's cubic cover. Only loophole: a non-canonical singularity of the Galois closure along the
  twisted diagonal (checkable; Part 3).
- (ii) Only <tau> acts; X/tau has kappa = 0; twisted quotients have wrong dimension and ages >= 1.
  Aut(X) = <tau> reduces (via "the 8 planes are the only planes") to a finite block-monomial check: running.
- (iv),(v): no q-companion, no solved companion; the Fermat sub-family's degree-6 witness cannot continue.
- Net: no construction mechanism reaches the sector; the obstruction is the omega-descent principle, and the
  only shapes not excluded are non-linear created reflections, for which the arrangement offers no candidate.
