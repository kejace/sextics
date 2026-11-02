"""Referee finalizer: patch the D8 placeholders in referee-report.md / comment.md from ref_d8_p999999937.json,
copy the report next to the receipts, and (if the D8 replay is unfinished) say so explicitly."""
import json, os, re, shutil
REF = '/tmp/hodge/results/diagonal-permutation/wave2-referee'
POSTS = '/tmp/hodge/posts/diagonal-permutation'
rep = open(f'{POSTS}/referee-report.md').read(); com = open(f'{POSTS}/comment.md').read()
f = f'{REF}/ref_d8_p999999937.json'
rec = json.load(open(f)) if os.path.exists(f) else {}
relax = rec.get('steps', {}).get('relax')
if relax and relax.get('good_reduction_16') and relax.get('Psi_in_J_mod_p') is False:
    cell = (f"Ψ(D8;u) ∉ J_16 mod 999999937 (good reduction in degrees 14 and 16; E1 nullity 954 = Koszul; "
            f"rank J_16^- = {relax['rank_J16_odd']}; u recomputed from scratch)")
    verdict = 'VERIFIED (scope: canonical representative + F-multiples only)'
    third = 'and, in the referee replay, $999999937$'
elif relax and relax.get('Psi_in_J_mod_p') is True:
    cell = 'Ψ(D8;u) ∈ J_16 mod 999999937 — DISAGREES with the agent receipts'
    verdict = 'UNSUPPORTED (replay disagrees)'
    third = 'but the referee replay at $999999937$ found $\\Psi \\in J$, so this item is disputed'
else:
    cell = 'replay not finished at the time of writing (see ref_d8.out)'
    verdict = 'PLAUSIBLE (two primes, transfer argument checked; third-prime replay unfinished)'
    third = 'referee replay at $999999937$ unfinished at posting time'
rep = rep.replace('REPLAY_D8_PLACEHOLDER', cell).replace('REPLAY_D8_VERDICT', verdict)
com = com.replace('__D8_THIRD_PRIME__', third)
open(f'{POSTS}/referee-report.md', 'w').write(rep); open(f'{POSTS}/comment.md', 'w').write(com)
shutil.copy(f'{POSTS}/referee-report.md', f'{REF}/referee-report.md')
print('verdict row 7c:', verdict); print('comment D8 clause:', third)
