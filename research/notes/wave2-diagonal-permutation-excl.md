# Wave 2 — diagonal-permutation — EXCLUSION / RIGOR lens

STATUS: in progress (no witness; this lens produces exclusions only)

Started 2026-11-01 (wave 2). Packet: /tmp/hodge/results/diagonal-permutation/wave2-excl/ ; log: /tmp/hodge/results/diagonal-permutation/wave2-excl.log.jsonl

## Plan
(a) re-verify wave-1 rigorous claims at new primes / new method (E1 kernels d<=8, D8, served(partials)=0).
(b) rigorous lower bound on deg D via E1 + "E2 mod J" relaxation (Lemma below) + msolve Groebner runs, per sector
    (Y-parity, mixed parity, with/without the D8 line).
(c) non-Jacobian E1-kernel elements (D8; d=9 kernel) and sector kernels d=10..12; structured families; full E1+E2 tests.
(d) packet + README; (e) draft post.

## Lemma (E2 modulo J)  [to be proven in README]
If (D,N,Q) is a certificate at q=1 then N(F)=A D + F H, N = u + K + (H/6)E with u(F)=A D, K Koszul; and
Psi(D) := D div u - u(D)  lies in J(F) (independent of the choice of u modulo J).  For D = w(F): Psi = -A * (w^T Hess(F) w) mod J.

## Results so far
(none yet)
