#!/usr/bin/env python3
"""Referee's own parser for msolve '-g 1' leading-ideal output. Variable names are read from the
'#variable order:' header of the OUTPUT (not from the .ms input). Verdict: every variable has a pure power
among the leading monomials  <=>  F_p[c]/LT(I) finite  <=>  the homogeneous ideal mod p has only the trivial zero.
Usage: python3 ref_ltparse.py file.g1.out [...]  -> one JSON line per file (also appended to referee.log.jsonl)"""
import sys, re, json, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
def check(path):
    txt = open(path).read()
    m = re.search(r'#variable order:\s*(.*)', txt)
    if not m or '[' not in txt:
        return {'file': path, 'ok': False, 'reason': 'no header / no leading-ideal block (msolve did not finish?)', 'bytes': len(txt)}
    names = [s.strip() for s in m.group(1).split(',') if s.strip()]
    p = int(re.search(r'#field characteristic:\s*(\d+)', txt).group(1))
    body = txt[txt.index('[') + 1: txt.rindex(']')]
    monos = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
    pure = {}; maxdeg = 0; bad = []
    for mono in monos:
        degs = {}
        for fac in mono.split('*'):
            mm = re.fullmatch(r'([A-Za-z_]\w*)(?:\^(\d+))?', fac.strip())
            if not mm: bad.append(mono); break
            degs[mm.group(1)] = degs.get(mm.group(1), 0) + (int(mm.group(2)) if mm.group(2) else 1)
        else:
            maxdeg = max(maxdeg, sum(degs.values()))
            if len(degs) == 1:
                v, e = next(iter(degs.items())); pure[v] = min(pure.get(v, 10**9), e)
    unknown = [v for v in pure if v not in names]
    missing = [v for v in names if v not in pure]
    ok = (not bad) and (not unknown) and (not missing) and len(monos) > 0
    return {'file': os.path.relpath(path, HERE), 'prime': p, 'nvars': len(names), 'n_leading': len(monos), 'max_lead_deg': maxdeg,
            'pure_power_exponents': sorted(set(pure.values())), 'vars_without_pure_power': missing, 'unparsed': bad[:3],
            'm_primary_mod_p': ok}
if __name__ == '__main__':
    for f in sys.argv[1:]:
        rec = check(f); rec['check'] = 'referee-ltparse'; rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S')
        print(json.dumps(rec))
        with open(os.path.join(HERE, 'referee.log.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
