#!/usr/bin/env python3
"""Check an msolve '-g 1' leading-ideal output for projective emptiness of a homogeneous ideal:
the DRL leading ideal contains a pure power of every variable  <=>  k[c]/I finite-dimensional  <=>  V_proj(I) = empty.
Usage: python3 wave2_excl_ltcheck.py file.ms file.g1.out [job-name]   -> prints + appends a JSON receipt."""
import sys, re, json, time
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
msf, outf = sys.argv[1], sys.argv[2]; job = sys.argv[3] if len(sys.argv) > 3 else ''
names = open(msf).readline().strip().split(',')
p = int(open(msf).read().split('\n')[1].strip())
txt = open(outf).read()
body = txt[txt.index('[') + 1:txt.rindex(']')] if '[' in txt else ''
lead = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
pure = {}
maxdeg = 0
for mono in lead:
    facs = mono.split('*')
    degs = {}
    for fac in facs:
        m = re.fullmatch(r'([A-Za-z_]\w*)(?:\^(\d+))?', fac.strip())
        if not m: raise SystemExit(f'cannot parse monomial {mono!r}')
        degs[m.group(1)] = degs.get(m.group(1), 0) + (int(m.group(2)) if m.group(2) else 1)
    maxdeg = max(maxdeg, sum(degs.values()))
    if len(degs) == 1:
        v, e_ = next(iter(degs.items()))
        pure[v] = min(pure.get(v, 10**9), e_)
missing = [v for v in names if v not in pure]
empty = (len(missing) == 0) and len(lead) > 0
rec = {'exp': 'msolve-leading-ideal-check', 'job': job, 'ms_file': msf, 'out_file': outf, 'prime': p, 'nvars': len(names),
       'n_leading_monomials': len(lead), 'max_leading_degree': maxdeg, 'pure_power_degrees': sorted(set(pure.values())),
       'variables_without_pure_power': missing, 'projectively_empty_mod_p': empty,
       'result': 'EMPTY (no nonzero solution over C, by the projective reduction argument)' if empty else 'NOT CERTIFIED EMPTY'}
print(json.dumps(rec), flush=True)
log2(**rec)
with open(outf.replace('.out', '_ltcheck.json'), 'w') as f: json.dump(rec, f)
