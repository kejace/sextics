#!/usr/bin/env python3
"""Parse msolve -g 1 output; the DRL leading ideal contains a pure power of every variable <=> projectively empty (rigorous over C for homogeneous input)."""
import sys, re, json, time
msf, outf, job = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else '')
names = open(msf).readline().strip().split(',')
p = int(open(msf).read().split('\n')[1].strip())
txt = open(outf).read()
body = txt[txt.index('[') + 1:txt.rindex(']')] if '[' in txt else ''
lead = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
pure, maxdeg = {}, 0
for mono in lead:
    degs = {}
    for fac in mono.split('*'):
        m = re.fullmatch(r'([A-Za-z_]\w*)(?:\^(\d+))?', fac.strip())
        if not m:
            if fac.strip() == '1': degs = {'1': 0}; break
            raise SystemExit(f'cannot parse {mono!r}')
        degs[m.group(1)] = degs.get(m.group(1), 0) + (int(m.group(2)) if m.group(2) else 1)
    maxdeg = max(maxdeg, sum(degs.values()))
    if len(degs) == 1:
        v, e_ = next(iter(degs.items())); pure[v] = min(pure.get(v, 10**9), e_)
missing = [v for v in names if v not in pure]
empty = (len(missing) == 0 and len(lead) > 0) or ('1' in pure)
rec = {'exp': 'msolve-leading-ideal-check', 'job': job, 'ms_file': msf, 'out_file': outf, 'prime': p, 'nvars': len(names),
       'n_leading_monomials': len(lead), 'max_leading_degree': maxdeg, 'pure_power_degrees': sorted(set(pure.values())),
       'variables_without_pure_power': missing, 'projectively_empty_mod_p': empty,
       'result': 'EMPTY over C (projective reduction argument)' if empty else 'NOT CERTIFIED EMPTY',
       'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
print(json.dumps(rec), flush=True)
with open('/tmp/hodge/results/dense-rational/wave2-construct.log.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
