# wave2-excl packet — `rank-one-biquartic` (q = 1, A = X0^3 Y0^3), EXCLUSION / RIGOR lens

All computations are exact (python-flint `nmod_mat` modular ranks, `fmpq_mat` exact rational solves, Fractions for
polynomials).  Instance: `F = f(X) + g(Y) + (5/7) L(X)^2 M(Y)`, `A = X0^3 Y0^3`, q = 1 (see dossier).

## Mathematical setup (what each script proves)
Let `J = (F_0..F_5)` (regular sequence; `X` smooth).  For a certificate `(D, N, Q)` of degree `d = deg D` with `D` coprime
to `F` (any Y-parity, mixed allowed):
1. **Euler gauge.** `N(F) = A D + F h` (mod-`F` reduction, `D` coprime); replacing `N` by `N - (h/6) Euler` leaves the
   cleared identity invariant, so WLOG `N(F) = A D` exactly.  Then `A D in J`: `D` lies in the E1-kernel
   `V_d = {D in S_d : A D in J_{d+6}}`.
2. **Koszul parametrization.** `N = N_0(D) + K(c)`, `K(c)_i = sum_j c_ij F_j`, `c` antisymmetric of degree `d-4`
   (syzygies of a regular sequence are Koszul).  The identity becomes `E2: D div N - N(D) = F Q`.
3. **Weight-2 ideal.** With `D = sum_k a_k D_k` over an exact Q-basis `(D_k, N_k)` of `V_d`, E2 is a system of polynomials
   of weight 2 (`a:1, c:1, Q:2`) in `(a, c, Q)` generating an ideal `I`; a certificate is a point of `V(I)` with
   `D(a) notin F * S_{d-6}` (the spurious subspace `Sp`).
4. **Hessian relaxation.** `B_D(K(c))` and `F Q` lie in `J_{2d}` (because `div K(c) = sum_{ij} (d_i c_ij) F_j`), so for
   every functional `lambda` on `(S/J)_{2d}` the quadratic form `phi_lambda(a) = lambda([B_{D(a)}(N_0(a))])` lies in `I`.
   For polar `D = dF(G)`, `[B_D(A G)] = -[A Hess_F(G,G)]` (identity (*) in `w2lib.py`).
5. **Certificates.** (a) `u^2 in I` for a coordinate `u` of `V_d / Sp`: an exact rational vector `lambda` with
   `lambda^T E' = (0 | u^2)` where `E' = [J_{2d} generators | quadratic columns]` (so `lambda` kills `J_{2d}` and hence
   all Koszul/Q columns).  (b) Exact quadratic forms in `I` (CRT-lifted rref rows, each certified as in (a)), followed by a
   Groebner computation over Q (Singular, char 0) showing `V(forms) subset Sp`.
   Either way: `V(I) subset Sp`, i.e. **no certificate of degree d**.
6. **E1-kernel over Q.** `dim_Q V_d <= dim_p V_d` whenever the `J_{d+6}` block has its regular-sequence rank mod p
   (verified in every run); `J_d subset V_d` always.  Non-Jacobian generators (`K_7`, `K_8`) are computed EXACTLY over Q
   (`k7exact.py`: pivot-submatrix exact solve, no rational reconstruction) with exact E1 lifts.

## Files
* `w2lib.py` — library (family systems, Hessian relaxation, exact certificates, sparse mod-J column representatives).
* `fam_run.py` — full weight-2 system (Koszul + Q columns) + exact certificates (used for d = 5).
* `hess_run.py` — Hessian relaxation per parity block + exact certificates `u^2 in I`.
* `hess_lift.py` — exact quadratic forms (CRT over primes + exact certification) + Singular-over-Q radical test.
* `quadforms_dump.py` — modular quadratic forms + Singular mod-p radical test (diagnostic).
* `k7exact.py` — exact non-Jacobian E1-kernel generators in the base-type subspace `S_d(X) M^a g^b`.
* `verify_e1kernel.py`, `verify_fixedD.py`, `quotJA_p*.sing` — re-verification of wave-1 / dossier claims at new primes.
* Receipts: `*.json`, `*.certs.json`, `*.log.jsonl`, `*.out`; everything is also appended to `../wave2-excl.log.jsonl`.

## Replay
```
cd /tmp/hodge/results/rank-one-biquartic/wave2-excl
python3 fam_run.py 5 1000000007 --cert                 # degree 5: all a_k^2 in I_2, exact certificates  (~10 s)
python3 hess_run.py 6 1000000007 --blocks 1 --cert     # degree 6: 25/35 coordinate squares certified   (~2 min)
python3 quadforms_dump.py 6 1000000007 && Singular -q quadforms_d6_p1000000007.sing   # mod-p radical test (~6 s)
python3 hess_lift.py 6 3 --singular --tag hesslift_d6  # exact forms + Singular over Q
python3 k7exact.py 7 && python3 hess_lift.py 7 3 --sparse --nonj nonj_d7.json --singular --tag hesslift_d7
```
