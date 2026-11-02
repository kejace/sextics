# Wave 2 — diagonal-permutation — VERIFY lens (adversarial referee)

STATUS: in progress (started 2026-11-02 00:15). Confirmed 0 / Refuted 0 / Weakened 0 / Unverified 0 — claims being extracted.

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
