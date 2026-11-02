# Wave 2 — lane `rank-one-biquartic` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; this lens produces exclusions and rigor audits, not witnesses). Started 2026-11-02 00:20 local.

Lane: `rank-one-biquartic` (q = 1, A = X0^3 Y0^3, F = f(X) + g(Y) + (5/7) L(X)^2 M(Y)).
Log: `/tmp/hodge/results/rank-one-biquartic/wave2-excl.log.jsonl`. Packet: `/tmp/hodge/results/rank-one-biquartic/wave2-excl/`.

## Plan
(a) adversarial re-verification of every rigorous claim of the dossier / wave-1 report at new primes (1000000007, 1000000009) and with independent code paths (engine `jacobian.mult_kernel`, own normalized-E2 solver, Singular/msolve);
(b) rigorous lower bound on deg D: E1 kernel = J for d <= 6 (both parities) + Hessian relaxation (quadratic system on the F-reduced E1 family, Groebner mod p, lifted by the homogeneous-reduction argument) at d = 5, 6 (Y-even, Y-odd, mixed);
(c) normalized E2 screening (mod p, full Koszul correction space) of all non-Jacobian E1-kernel generators and random members at d = 7, 8, 9 (both parities), structured families, degree-9 mixed sums;
(d) packet; (e) draft post.

## Results so far
(none yet)
