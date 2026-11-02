# Reproducibility packet — `short-power-presentation`, wave 2, exclusion/rigor lens

Instance: q = 1, A = X0^3 Y0^3, F = eX + eY + c((L+M)^6 + (L−M)^6), c = 5/7^7, L = X0+2X1+3X2, M = Y0+3Y1+5Y2
(variables X0,X1,X2,Y0,Y1,Y2; τ = global Y-sign; Y-parity of D is called even/odd; "mixed" = D = D^+ + D^-).

## Exact statements (see `receipts.jsonl` for primes, ranks, Hilbert functions, timings)

* **T1 (E1 kernel).** K_d^± := {D ∈ S_d^± : A·D ∈ J(F)}. dim_Q K_d^±/J_d^± = 0 for d ≤ 6 (both parities),
  = 1 for d = 7 (both parities), and the d = 7 classes have exact rational representatives D7^+ ∈ Q[X0,X1,X2,M],
  D7^- ∈ Q[Y0,Y1,Y2,L] with exact numerator certificates (N(F) = A·D7 over Q). Modular upper bounds use the
  Hilbert-series check (prime good in degrees d and d+6; rank_p ≤ rank_Q); the lower bound at d = 7 is the exact
  element.
* **T2 (degree bound).** Every certificate denominator D satisfies: deg D ≥ 6, and if D is Y-pure then deg D ≥ 7.
  (Proof: for d ≤ 5 and for Y-pure d = 6, K_d = J_d, so D = w(F); the E2 relaxation Ξ(w) := −A·w^T Hess(F) w ∈ J_{2d}
  has, in the Euler slice, an m-primary ideal mod a good prime (Hilbert function reaches 0), hence by
  h_{S/I_p} ≥ h_{S/I_Q} the only solutions over Q-bar are w = G·E, i.e. D = 6GF, not coprime to F.)
* **T2' (mixed parity, d = 6).** For D = D^+ + D^- (both parts nonzero) the odd part of the relaxation has no cross
  terms; linear functionals vanishing on J_12^- and on the span L^+ of the coordinate polynomials of Ξ^+ give 36 quadrics
  in the 18 odd coefficients whose ideal is m-primary mod 1073741789 (Singular dim 0, vdim 4096); dim_Q L^+ = 39 = dim_p L^+
  is certified by 114 rational dependency certificates verified with fmpq_mpoly (`exp` = mixed-annihilator-exact-certificate),
  so the functionals are reductions of rational ones. Hence w^- = 0, then w^+ = 0 by T2. **Theorem: every certificate
  denominator has degree ≥ 7.**
* **T3 (Y-pure d = 7).** D = λ·D7 + w(F). Block elimination (`excl_blocks.py`, `exp` = block-elimination-step) mod p:
  λ = 0 and the 18 coefficients w_{k,X_iY_j} (j ≠ k, Y-components) vanish; residual: 42 unknowns / 231 quadrics —
  see `receipts.jsonl` for the final verdict and the exact certificates.

## Replay
```
cd /tmp/hodge/results/short-power-presentation
python3 excl_e1kernel.py 0 1e 1o ... 10o          # E1 kernels, prime PRIMES30[0] = 1073741789, engine cross-check
python3 excl_relax.py <d> <e|o|m> [--pidx k]       # relaxation system -> msolve -g 1 (+ Singular), verdict in the log
python3 excl_kernel.py 9                           # exact D7^±, kernel generation by K_7, lambda-only sub-cases
python3 excl_slice.py <system.json> k1 k2          # exploratory: random linear slices of a relaxation variety
```
Engine: `~/devel/sextics-work/hodge` (core.load_instance, python-flint nmod_mat), msolve, Singular 4.3.2.
Every log line (`/tmp/hodge/results/short-power-presentation/wave2-excl.log.jsonl`) records exp, d, parity, prime,
ranks/Hilbert data, seconds, job name.  The derivation of the relaxation and the rigor framework are in
`/tmp/hodge/notes/wave2-short-power-presentation-excl.md` §0.

## Pending (remote) outputs at closing time
`excl_residual_d7e_p0.{ms,lead,msolve.log,sing.out}` (d = 7 Y-even residual Gröbner), `exclB-spp-certify7e.out`
(exact certificates, `excl_certify.py 7 e`), `exclB-spp-blocks2.out` (d = 7 odd/mixed, d = 8 block eliminations),
`exclB-spp-relax-mixed.out` (d = 7 mixed by R4), `exclB-spp-kernel2.out` (K_9, λ-only d = 9) — all in the remote
lane directory; their log lines are appended to the remote `wave2-excl.log.jsonl`.
Caveat on verdict strings: a `block-elimination` record with `result = residual-not-empty` and an EMPTY msolve output
(`msolve_hilbert` = [1, n, ...] with no leading monomials) means the residual Gröbner step died (memory), not that
solutions exist.
