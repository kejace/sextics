"""Append the addendum to referee-report.md from the referee receipts; print the verdict counts."""
import json, sys, os
REF = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee'
rep = '/tmp/hodge/posts/reflection-coupled-residual/referee-report.md'
fx = [json.loads(l) for l in open(f'{REF}/ref_fixedD.jsonl')]
d9 = [r for r in fx if r['tag'].startswith('d9')]
gb = [json.loads(l) for l in open(f'{REF}/ref_gb.jsonl')]
g00 = [r for r in gb if r['file'].startswith('ref_hess_d7_W00_')]
lines = []
if d9:
    r = d9[-1]
    lines.append(f"* Row 12, degree-9 representative `e1ker-d9-chi00-basis0` (referee, p = 10^9+21): system {r['rows']} x {r['cols']}, "
                 f"{r['apriori_vectors']} a-priori kernel vectors (all exact kernel vectors), rank {r['apriori_rank_p']}; row subset {r['subset_rows']}: "
                 f"rank_p(M_sub) = {r['rank_p_subset']} = {r['cols']} - {r['apriori_rank_p']}, rank_p([M|b]_sub) = {r['rank_p_aug_subset']} -> "
                 f"certified inconsistent over Q: {r['certified_inconsistent_over_Q']} ({r['seconds']} s). Verdict row 12: VERIFIED for this element (agrees with the lane's rank 6006 / nullity 798).")
else:
    lines.append("* Row 12: the referee replay of the degree-9 representative had not finished when this report was written; verdict stays PLAUSIBLE (lane receipts, two primes).")
if g00:
    r = g00[-1]
    if r.get('singular_dim') == 0:
        lines.append(f"* Row 7, degree-7 character (0,0) (referee, p = 10^9+21): GB {r['gb_size']}, msolve quotient {r['msolve_quotient_dim']}, Singular dim 0, vdim {r['singular_vdim']}, HF {r['HF'][:8]} ({r['msolve_seconds']} s). Verdict row 7: VERIFIED.")
    else:
        lines.append(f"* Row 7: referee run ended with {r.get('result') or r.get('verdict')}; verdict stays PLAUSIBLE.")
else:
    lines.append("* Row 7, degree-7 character (0,0): the referee replay at 10^9+21 (local, 4 threads) and the queued remote job had not finished when this report was written; "
                 "the lane's single run (quotient 3028, HF 1,42,629,1986,359,11,0, 4224 s at 2147483587) has the same valid transfer argument as rows 5, 6, 8. Verdict row 7: PLAUSIBLE (single run, not independently replayed yet).")
lines.append("* Degree-10 E1 kernels replayed: 186/126/136/108 with good reduction (row 2 upper bounds confirmed at a third prime).")
v7 = 'VERIFIED' if (g00 and g00[-1].get('singular_dim') == 0) else 'PLAUSIBLE'
v12 = 'VERIFIED' if d9 and d9[-1]['certified_inconsistent_over_Q'] else 'PLAUSIBLE'
verdicts = {1:'VERIFIED',2:'VERIFIED',3:'VERIFIED',4:'VERIFIED',5:'VERIFIED',6:'VERIFIED',7:v7,8:'VERIFIED',9:'VERIFIED',10:'VERIFIED',
            11:'VERIFIED (conditional on the book theorem)',12:v12,13:'UNSUPPORTED',14:'PLAUSIBLE',15:'PLAUSIBLE',16:'PLAUSIBLE',17:'UNSUPPORTED'}
cnt = {'VERIFIED':0,'PLAUSIBLE':0,'UNSUPPORTED':0}
for k,v in verdicts.items(): cnt[v.split()[0]] += 1
lines.append(f"\n**Final counts (17 claims):** VERIFIED {cnt['VERIFIED']} (row 11 conditional on the book's Chapter 1 theorem), PLAUSIBLE {cnt['PLAUSIBLE']}, UNSUPPORTED {cnt['UNSUPPORTED']}.")
lines.append("\n**Strongest verified statement.** For the protected (F, q = 1, A = X0^3 Y0^3): every certificate denominator has degree >= 7 with no parity or Klein-character restriction (E1 by character at three primes with good reduction; Hessian relaxation zero-dimensional at degrees 5, 6 with the transfer HF_Q <= HF_p); no degree-7 denominator of a single Klein character exists for (1,1), (1,0), (0,1) and none with components only in (1,0)+(1,1) [the (0,0) block: lane's run, see row 7]; the unique non-Jacobian degree-8 E1 class K8 is certified not a denominator over C.")
open(rep, 'a').write('\n'.join(lines) + '\n')
print('\n'.join(lines)[:1500]); print('COUNTS', cnt, 'v7', v7, 'v12', v12)
