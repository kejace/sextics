# Wave 2 — `reflection-coupled-residual` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; re-verification and rigorous lower bound being computed)

Started 2026-11-01 ~19:15 local. Work dir `/tmp/hodge/results/reflection-coupled-residual/wave2-excl/`,
log `/tmp/hodge/results/reflection-coupled-residual/wave2-excl.log.jsonl` (one JSON line per experiment).
Remote jobs named `excl-reflection-coupled-residual-*`.

## Plan
(a) re-verify wave-1/dossier rigorous claims with a different prime (p = 2147483587 / 2147483579; wave 1 used 2147483647, 2147483629; dossier used 32003), with explicit good-reduction (Hilbert function) checks;
(b) rigorous lower bound on deg D via E1 by Klein character + Hessian relaxation of E2 on the Jacobian family;
(c) non-Jacobian E1-kernel classes and structured families, degrees 7-12, character sectors, full E1+E2 (symmetric solver);
(d) packet; (e) draft post.

## Results so far

### R1. Degree 5 (Hessian relaxation, local, p = 2147483587) — DONE
For D = Σ c_i F_i (c_i constants; this is the whole E1-admissible space at d = 5 once job e1char confirms
"E1 ⇒ D ∈ J_5" in all four characters), the E2-mod-J obstruction is ob(c) = −[A Σ c_i c_j F_ij] ∈ (S/J)_10.
The Hessian has exactly one vanishing off-diagonal entry, F_{X2 Y2} = 0; the other 20 classes [A F_ij] are
linearly independent in (S/J)_10 (rank 20 in each of the four target-character blocks: sizes 6/10/2/2;
good reduction verified: rank J_10^ψ = dim S_10^ψ − Fermat count for all ψ). Hence c_i c_j = 0 for all 20 pairs,
in particular c_i² = 0, so c = 0. **No denominator of degree 5 (any character, mixed parity included).**
Receipts: `wave2-excl/hess.jsonl` (exp phi-kernel/phi-verdict, d=5).

### Jobs submitted (remote)
* `excl-rcr-e1char`: E1 kernels by Klein character, d = 1..10, p = 2147483587 (+ 2147483579 for reconstruction).
* `excl-rcr-hess`: Φ-kernel at d = 6 (two primes); Macaulay-3 rank of the quadratic system per character at d = 7;
  msolve exports of the pair/full systems at d = 7.
* `excl-rcr-k8job`: rigorous fixed-D test of K8 (restricted + full), K8 ∉ J, exact E1 lift, e2lin.

