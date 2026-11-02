# Draft discussion comment for #118 (`rank-one-biquartic`) — rigorous low-degree exclusion  [DRAFT — NOT POSTED]

### Low-degree denominators are excluded: no certificate with deg D <= {DMAX_PROVED} (exact, characteristic zero)

AI-assisted computation (Claude). Nothing here claims nonexistence of a certificate; only denominators of bounded
degree are excluded, and the bound is stated with its exact scope.

**Statement.** For the exact target of this instance (`F = f(X) + g(Y) + (5/7) L(X)^2 M(Y)`, `A = X0^3 Y0^3`, `q = 1`)
there is no certificate `(D, N, Q)` with `deg D <= {DMAX_PROVED}` — for *every* nonzero `D` coprime to `F`, of any
Y-parity (mixed parity included), with `Q` arbitrary (no denominator bound, no strictness assumed).
Degrees {DMAX_PROVED+1}..{DMAX_SCREENED} are excluded modulo the prime(s) listed below (screened, not proved).

**Method** (all steps are exact; scripts and receipts in the linked packet).
1. Euler gauge: WLOG `N(F) = A D`, so `A D in J(F)`; the E1-kernel `V_d = {D : A D in J}` is `J_d` for `d <= 6` and
   `J_d + K_d` for `d = 7, 8` with `dim K_7 = 3` (Y-even), `dim K_8 = 9 + 9`; the non-Jacobian generators are computed
   exactly over Q (no rational reconstruction) with exact E1 lifts.
2. Koszul parametrization `N = N_0(D) + K(c)` reduces the identity to `E2: D div N - N(D) = F Q`, a weight-2 system in
   the family coordinates `a` (of `D`), `c`, `Q`.
3. Hessian relaxation: `B_D(K(c))` and `F Q` lie in `J_{2d}`, so every functional on `(S/J)_{2d}` applied to
   `B_D(N_0(D))` is a quadratic form in `a` that must vanish; for polar `D = dF(G)` this is `[A Hess_F(G,G)] = 0`.
4. Exact certificates: rational vectors `lambda` vanishing on `J_{2d}` with `lambda(E2) = u^2` for each coordinate
   `u` of `V_d / (F S_{d-6})` (or exact quadratic forms in the ideal), verified exactly; then a Groebner computation
   over Q (Singular, char 0) shows the common zero set is the spurious locus `{D in F S_{d-6}}` (not coprime).

**Results.**
| d | E1-kernel `V_d` (dim, both parities) | weight-2 squares certified / coordinates | Groebner over Q | status |
|---|---|---|---|---|
{RESULTS_TABLE}

**Re-verification of earlier claims** (wave-1 report / dossier, different primes and code paths): {REVERIFY_SUMMARY}

**Scope and caveats.** Modular ranks are lower bounds for rational ranks; every "proved" line above rests on exact
rational certificates and a char-0 Groebner basis, not on modular results. "Screened" lines are modular only.
Nothing is claimed about denominators of degree > {DMAX_SCREENED}; the first unexcluded space is the degree-{DMAX_SCREENED+1}
E1-kernel family (`dim >= ...`).

**Reproduction.** Packet `/tmp/hodge/results/rank-one-biquartic/wave2-excl/` (to be attached as a gist): `README.md`
lists the commands (`fam_run.py 5`, `hess_pipeline.py 6/7/8`, `k7exact.py`, verification scripts) and the receipts
(`*.certs.json` exact certificates, `*.D.out` Singular char-0 transcripts, `wave2-excl.log.jsonl`).
Primes used: 1000000007, 999999937, 999999929 (CRT lifts); 65521 / 32003 (Singular quotient re-checks).
