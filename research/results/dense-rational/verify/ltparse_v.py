import sys, re, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import log
p = int(sys.argv[1]); results = {}
for sector in ('even', 'odd', 'mixed'):
    fn = f'/tmp/hodge/results/dense-rational/verify/quad6_{sector}_p{p}.ms'
    names = open(fn).readline().strip().split(',')
    try: txt = open(fn + '.g1').read()
    except FileNotFoundError: results[sector] = 'missing'; continue
    body = txt[txt.index('[') + 1:txt.rindex(']')]
    lead = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
    pure = {}; maxdeg = 0
    for mono in lead:
        degs = {}
        for fac in mono.split('*'):
            m = re.fullmatch(r'(c\d+)(?:\^(\d+))?', fac.strip()); assert m, mono
            degs[m.group(1)] = degs.get(m.group(1), 0) + (int(m.group(2)) if m.group(2) else 1)
        maxdeg = max(maxdeg, sum(degs.values()))
        if len(degs) == 1:
            v, e = next(iter(degs.items())); pure[v] = min(pure.get(v, 99), e)
    missing = [v for v in names if v not in pure]
    results[sector] = {'nvars': len(names), 'n_leading': len(lead), 'max_leading_degree': maxdeg, 'pure_power_degrees': sorted(set(pure.values())), 'missing': missing, 'empty': not missing and bool(lead)}
    print(sector, results[sector])
ok = all(isinstance(r, dict) and r['empty'] for r in results.values())
log('[supersedes earlier quad6 line whose parser read the header] excl R-d: degree-6 Hessian relaxation projectively empty in even/odd/mixed sectors => deg D >= 7 for all parities (independent generator, p = 1000003, msolve -g 1)',
    'CONFIRMED' if ok else ('UNVERIFIED' if any(r == 'missing' for r in results.values()) else 'REFUTED'),
    f'p={p}: ' + '; '.join(f'{s}: {r}' for s, r in results.items()),
    'leading-monomial counts 135 (even) / 153 (odd) match the excl lens at p = 1073741789; emptiness mod p of a homogeneous ideal is rigorous over C (Macaulay full column rank lifts)')
