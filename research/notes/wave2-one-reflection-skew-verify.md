# wave2 verification ledger: one-reflection-skew (adversarial referee)

STATUS: in progress — confirmed 0 / refuted 0 / weakened 0 / unverified 0 (updated continuously; see table)

Scope: every claim labelled rigorous/proved/exact/computed in `dossier-one-reflection-skew.md`, the lane logs under
`/tmp/hodge/results/one-reflection-skew/`, and any `attack1-`/`wave2-one-reflection-skew-*` notes and draft post as they appear.
Independent tooling: own polynomial code (`/tmp/hodge/verify/ors/vlib.py`, no engine import), python-flint `nmod_mat`
at primes 1000003 / 999983 (lane used 2^31-1, dossier used 32003), Singular 4.3.2 (`jac.sing`), msolve where needed.
JSON lines for every check: `/tmp/hodge/results/one-reflection-skew/wave2-verify.log.jsonl`.

Modular-logic reminders applied throughout:
* rank_p(M) <= rank_Q(M).  "Full rank mod p" => full rank over Q (OK).  "Kernel of dimension k mod p" => Q-kernel has dimension <= k; equality needs an exhibited Q-kernel of dimension k (e.g. J_d inside (J:A)_d).
* "Inconsistent mod p" alone does NOT prove Q-inconsistency (example: p x = 1).  It does when rank_p(M) = rank_Q(M), e.g. certified by an exhibited Q-kernel of the right size, or by an integer dual witness y with y^T M = 0 over Z and y^T b != 0.
* Good-prime condition for "independent modulo J over Q" from a mod-p computation: dim_p J_n must equal the char-0 value (Hilbert function of the complete intersection).

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
