STATUS: RUNNING — 2026-11-02 04:50 UTC — A = A_K3 + A_mix verified exactly; A_mix nonzero mod J; E1 kernels d≤7 = J(F)_d only; incidence-degree (msolve) and seed-lift computations launching

# Wave 4 — polarized-k3-remainder — route `amix-component`

Owner: wave-4 subagent. Log: /tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl. Scripts/data: /tmp/hodge/results/polarized-k3-remainder/wave4/
Prior work credited: wave-2 construct note (chart reduction, (3,2) scalar identity, F_5 seed), probe note (Aut(C_f)=1, ρ≤4), book ch. 3A, `1337isnot1337` (crossed-ternary companion construction, PR #149), `fqlx` (moving-conic/PolynomialCurves), thread #115 (`1337isnot1337`).

## Results table (running)
| # | step | status | result | file |
|---|------|--------|--------|------|
| 1 | exact split A = A_K3 + A_mix (physical coords) | PROVEN (exact Q) | A_mix = 13 Y-odd monomials, see amix_pieces.json | wave4/amix_pieces.json |
| 2 | Jacobian membership of A, A_K3, A_mix, A_2, A_4 | modular (p=2^31−1, 2^31−19) | none in J(F): rank 18→19 (A, A_K3, A_mix), 36→37 (A_2, A_4) | log |
| 3 | E1 kernel of [A_mix] on (R/J)_d, d=4..7, both parities | modular (p=2^31−1) | kernel 0; admissible space = J(F)_d (dims 0,3,18,63) — identical to [A] | wave4/e1_kernels_d4to7.json |
| 4 | E1 kernel d=8,9 | running | | |
| 5 | degree of the 13×13 profile-A incidence system at seed params (msolve mod 65521) | launching | | |
| 6 | 5-adic lift of the F_5 seed to high precision + exact finite-field seed data (Jacobian inverse, freeness right inverses) | launching | | |
| 7 | independent re-verification of chart identity (b=2,4) and (3,2) scalar conversion | launching | | |
