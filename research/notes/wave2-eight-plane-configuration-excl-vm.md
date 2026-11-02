# Wave 2 — `eight-plane-configuration` — lens EXCLUSION / RIGOR — VM-side agent report

STATUS: in progress — no witness. (Companion report by the other excl agent: wave2-eight-plane-configuration-excl.md; this file holds the VM-side agent's independent checks and additions, kept separate to avoid sync clobbering.)

Log: /tmp/hodge/results/eight-plane-configuration/wave2-excl-vm.log.jsonl
Packet: /tmp/hodge/results/eight-plane-configuration/wave2-excl/vm/

## Findings so far
1. [exact] Six partials dF/dZ_i (d=5): Koszul-reduced E2 system (my independent implementation kred.py, ALL monomials, no parity filter, 3003 x 216): inconsistent mod 2^31-1 and mod 2147483587, and exact fmpz ranks rank_Q(M)=210, rank_Q([M|b])=211 -> inconsistent over Q (rigorous). Agrees with wave-1 and with the other agent's parity-restricted exact ranks (102/103).
2. [correction] Coefficient denominators of F: lcm = 7^4 * 11^4 * 13^2 * 17^2 (the dossier said products of 7, 11, 13 only; 17 also occurs). Irrelevant for the primes used (>2^31-100), but recorded.
3. [rigor gap flagged] The other agent's (R1) "WLOG D has pure Y-parity" is a WLOG for existence (D -> D*tau D doubles the degree), NOT for degree-by-degree exclusion: a mixed-parity D of degree d is not covered by pure-parity sector computations at degree d. E1 does split (J is parity-graded) so mixed D of degree <= 9 still lies in J, but the Hessian relaxation / E2 must be run on the full (mixed) w-space. Being closed here for d = 5, 6 (see below).

## Plan (complementary to the other agent; no redo of their results)
(a) mixed-parity closure d=5,6; independent Singular/flint cross-checks of the msolve results and of the E1 ladder.
(b) structured Jacobian subfamilies at d=7,8,9 (pure and mixed) via Hessian relaxation restricted to the subfamily (msolve) and the full Koszul-reduced E2 (modular with full-column-rank rigor check / exact).
(c) random-sample generic exclusion with rigor receipts; (d) packet README; (e) draft post (-vm suffix).
