# Wave 2 — diagonal-permutation — VERIFY lens (adversarial referee)

STATUS: in progress (started 2026-11-02 00:15; updated 00:40). Confirmed 16 / Refuted-as-stated 1 (attack1 headline "deg D <= 7 excluded" overclaims: only deg D <= 4 is excluded by E1; 5<=d<=8 forces D into J(F)_d (+Q*D8 at d=8)) / Weakened 3 (modular inconsistency stated as proof; Fermat-kernel description) / Pending: independent E1 kernels (own prime, Singular), D8 exact lift, exact dual certificates for fixed-D exclusions.

Log (one JSON line per check): /tmp/hodge/results/diagonal-permutation/wave2-verify.log.jsonl
My scripts: /tmp/hodge/verify/diagonal-permutation/ (independent of the `hodge` engine; own polynomial code, own primes, Singular where feasible)

Sources refereed: dossier-diagonal-permutation.md, attack1-diagonal-permutation.md, wave2-diagonal-permutation-construct.md (re-read every 30 min), lane logs under /tmp/hodge/results/diagonal-permutation/.

## Claim table

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
| 1 | F has 24 terms = Fermat + (5/7)O1 + (5/11)O2 + (5/13)O3; Y-degree profile 3/12/6/3 | dossier 1.1 | CONFIRMED | v_struct.py (own sympy parse) |  |
| 2 | Only sign symmetries (up to scalar): +-1, +-tau | dossier 1.1/1.4 | CONFIRMED | v_struct.py, 64 patterns |  |
| 3 | Signed-permutation symmetry group = S3 x <tau> (12 projective elements), simultaneous pair permutations; no X<->Y swap; no reflection | dossier 1.4, attack1 1.1 | CONFIRMED | v_struct.py: 46080 signed perms, 24 preserve F; no element with 5-dim fixed space | transpositions: det +1, fixdim 4; tau*s: det -1, fixdim 3; 3-cycles fixdim 2; tau*c fixdim 0 |
| 4 | Projective diagonal group = {1,tau} (Smith form diag(1,1,1,1,2,0)) | dossier 1.4 | CONFIRMED | v_struct.py sympy SNF |  |
| 5 | A o g in {+-(T_i-T_j)}; S3-average 0; A o tau = -A; A spans V_std; residue transforms like A since det(S3)=+1 | dossier 1.5 | CONFIRMED | v_struct.py table; det values in claim 3 | residue statement uses det=+1 (checked)  |
| 6 | Equivariant decomposition R_6^- = 39 triv+31 sign+69 std; R_12^- = 154+136+288 std; R_6^+ = 46+28+72 std; S_6^- character (226,10,1); G-invariant dims of R_d (d=6..13) = 46,63,94,115,144,153,167,153 | dossier 1.5, 4.1 | CONFIRMED | v_chars.py equivariant Koszul series (own code) | std multiplicities for R_13..19^- not yet re-derived (series truncated); pending |
| 7 | At Fermat, primitive Ann(A)_8^+ = 202, 'spanned by monomials in (X0^2,Y0^2) cap (X1^2,Y1^2)' | attack1 2.6, 4.4 | WEAKENED | own exact kernel at Fermat: 202 (two primes, matrix has entries 0,+-1); only 193 are such monomials, the other 9 are binomials (T0+T1)X_iY_j | number right, description incomplete; harmless |
| 8 | HEADLINE: no certificate with deg D <= 7; degree-8 certificate must use D8 | attack1 status, 4.1 | REFUTED-AS-STATED (overclaim) | logic (see log) | E1 only forces D in J_d (5<=d<=7) and D in J_8 + Q*D8 at d=8; only 18 specific Jacobian denominators were (modularly) excluded. Rigorous: deg D <= 4 impossible. |
| 9 | fixed-D inconsistency / served_dim=18 for 18 Jacobian denominators => Y_{G(F)} supports no odd (3,1) class | attack1 2.4, 4.2 | WEAKENED | one prime, compressed rows | not a proof over Q; exact dual certificates being computed (v_fixedD.py --dual) |
| 10 | D8 inconsistent (2 primes) => deg D >= 9 and all D8^m excluded | attack1 4.1 | WEAKENED | logic | two primes not a proof; J_8 and J_5..7 members unexcluded; prescribed-support step itself fine |
| 11 | E1 <=> [A][D]=0 in R_F, all parities incl. mixed | dossier 1.7, attack1 1.3 | CONFIRMED | logic |  |
| 12 | zero primitive kernel mod p => zero over Q, C | dossier 4.1, attack1 2.1 | CONFIRMED | needs dim J_{d+6} mod p = Koszul; verified at p=998244353 (own) and consistent at p=2147483647 (engine data) | over Q: K_d^Q = J_d^Q exactly for d=5,6,7 |
| 13 | Gorenstein duality and derived k_d values | attack1 1.3 | CONFIRMED | own series: R_24^+=1, R_24^-=0; arithmetic rechecked | over Q k_8^+=1 needs exact A*D8 in J (pending) |
| 14 | WLOG G-invariant D; A <=> B_2; e_std identity | dossier 1.7, 3.1 | CONFIRMED | hand check | (01)-average gives anti-equivariant N |
| 15 | age table; tau*c free on X | dossier 1.6 | CONFIRMED (wording) | own eigenvalue bookkeeping | fixed points of tau*c are (1,w^k,w^2k) eigenvectors, not coordinate points |
| 16 | p_g(X/K)=1 for K<=S3 => no p_g=q=0 fibration, no uniruled model | attack1 1.2 | CONFIRMED | logic |  |
| 17 | prescribed-support use at q=1 (D^m excluded by exact failure over D) | dossier E2, attack1 4.1 | CONFIRMED | docs/01 statement | needs EXACT inconsistency |
| 18 | Lie algebra of linear symmetries 1-dim; Hessian symmetrizer 1-dim | dossier 1.4 | CONFIRMED | own sympy rank (Q) and p=10007,1000003 |  |
| 19 | std multiplicities 288,281,248,207,156,111,69,41 in R_12..19^-; Hilbert function table | dossier 4.1, attack1 1.3 | CONFIRMED | own series |  |
| 20 | explicit (X)^3 cap J_7 family lies in (X)^3 | attack1 2.3 | CONFIRMED | own expansion |  |
