STATUS: RUNNING — 2026-11-01T23:40Z — setting up twisted cubic-character sums for Frobenius traces on H^3(T_H)_χ^{s=±1}
# wave4 · four-plus-two-residual · route "frobenius"

Side note for the user question "what is diamondjm": no string `diamondjm` occurs in the repo or in /tmp/hodge. The lane's
discussion-thread author is `benediamond` (character computation: binary factor = σ-anti-invariant Prym line dt/w). GitHub has an
unrelated account `DiamondJM` (Josh Diamond) with no footprint in this project.

## Plan
1. Sector: V = H^3_c(P^3 − S_H, L_χ(H))^{s=−1} ≅ H^3(T_H)_χ^{s=−1}, K-rank 48, Hodge (0,32,16,0); V^+ rank 56, (1,36,19,0).
2. Traces at k=1 from Σ_{P^3(F_p)} χ(H) and the s-twisted sum Σ χ(H_n) (n non-residue): Tr(Frob|V^±) = −(S0 ± S1)/2.
3. Rank-48 char polys are out of reach by counting; deliverable = traces for all p ≡ 1 (3) up to ~2000, exact K-cubics/quartics for
   P, E', C_f, C_g at small p, Sato–Tate second moments and cross-correlations (multiplicity estimates), plus Hodge-polygon congruence checks.
