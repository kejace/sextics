# Draft discussion comment (klein-coupled-character, #116) — wave 2, exclusion lens
(DRAFT — numbers to be finalized from /tmp/hodge/results/klein-coupled-character/wave2-excl/; do not post yet)

**Partial update (exclusions only, no certificate): every certificate denominator has degree >= 7, and (if the
degree-7 analysis closes) >= 8.** Target unchanged: q = 1, A = X0^3 Y0^3 on the coupled Klein sextic.

Setup. Write the certificate in normalized form (ch. 1): E1 N(F) = A D, E2 D div N − N(D) ∈ (F). E1 forces
A D ∈ J(F). Modulo p = 1000003 (also 32003 and 2^31−1 in earlier runs) the colon space (J:A)_d equals J(F)_d for
every d <= 8 in both Y-parities and all seven sigma-characters, and dim J_d mod p equals the complete-intersection
value in every sector; since (J:A) over Q reduces into (J:A) mod p this gives (J:A)_d = J_d over Q for d <= 8.
So every admissible denominator of degree <= 8 is D = Σ w_i F_i.

Hessian relaxation. For D = Σ w_i F_i every numerator is N = A w + K(p) with p antisymmetric (Koszul, F_i regular
sequence). Since D div K(p) − K(p)(D) ∈ J for every D, E2 modulo J reads w^T Hess(F) w ∈ (J : A)_{2d−6}: a
homogeneous quadratic system in the coefficients of w, invariant under the denominator gauge w → w + l z
(D → D + 6 l F). Projective emptiness of this system mod p implies emptiness over Q̄ (reduce a Q̄-point).

Results (all exact or modular-with-semicontinuity; scripts and receipts linked below):
* d <= 4: (J:A)_d = 0. d = 5: relaxation empty (both parities and mixed; Macaulay matrix of degree 2 has full rank).
* d = 6, Y-odd: relaxation empty (msolve/Singular: zero-dimensional cone).
* d = 6, Y-even and mixed: the relaxation locus is a single irreducible curve C, contained in the "weighted Euler"
  family D_λ = Σ λ_i Z_i F_i (over Q: every off-diagonal coordinate squares into the relaxation ideal); in the
  weighted-Euler P^4 (gauge λ_0 = 0) C is the complete intersection of three explicit quadrics. For this
  σ-invariant family the Koszul-form E2 system (59 unknowns, 294 rows) has 150 constant left-kernel vectors; they
  produce 8 exact quadrics in λ that must vanish at any consistent λ, and over Q these quadrics cut out exactly the
  gauge line λ ∝ (1,…,1) (D = 6F). Hence no D ∈ J_6 of any parity carries a certificate: **deg D >= 7**.
* Lemma (non-zero-divisors): every nonzero linear form in span(X0,X1,X2) or span(Y0,Y1,Y2) is a NZD on R_6
  (quadratic syzygies of the restricted partials: module annihilator = (1) in every chart). By the Gysin description
  of N^1 (ch. 3) served(l·D) = served(D) modulo J for such l, so such linear factors can be stripped from denominators.
* d = 7, σ-eigen pure-parity denominators: the relaxation locus is empty for Y-even characters 0,3,5,6 and Y-odd
  characters 1–6 (no such denominators), and for Y-even characters 1,2,4 it equals the lifted curves X_j·C exactly over Q
  (characteristic-0 Gröbner bases), whose points X_j D_λ are excluded by the NZD lemma + sum principle + the d=6 result.
  Y-odd character 0 (the lifted family span(Y)·C): [TO FILL]. Non-eigen degree 7: [TO FILL].
* d = 8 (σ-eigen sectors): [TO FILL]
* d = 9, 10: [TO FILL: kernel reps, random combos — modular screens; exact duals for Jhat, polars, polar products]
* Re-verification of earlier claims: Jhat exclusion (exact dual), polars, Z_j F_i, products — agree at p = 1000003;
  the earlier "served-targets" numbers were computed in character 4 (not 3) and are withdrawn/redone: [TO FILL].

What is NOT claimed: nothing about d >= 9 beyond specific denominators (modular rank evidence unless an exact dual is
listed); nothing about nonconstructive existence.

Reproduction: /tmp/hodge/results/klein-coupled-character/wave2-excl/ (w2lib.py, local_d6/, local_d7/, local_nzd/,
receipts jobA_*/jobC_*/jobD_*), log wave2-excl.log.jsonl. [gist link TO FILL]
