"""Hilbert function of S/LT(I) from an msolve '-g 1' (leading ideal) output.
Usage: python3 hilb_from_msolve.py <msolve-input.ms> <msolve-output> <maxdeg>
Pure-Python standard-monomial count (numpy) per degree up to maxdeg, plus Singular hilb() cross-check when available."""
import sys, re, subprocess, json, itertools
import numpy as np
inp, out, maxdeg = sys.argv[1], sys.argv[2], int(sys.argv[3])
vars_ = open(inp).readline().strip().split(',')
n = len(vars_)
vidx = {v: i for i, v in enumerate(vars_)}
txt = open(out).read()
body = txt[txt.find('[') + 1: txt.rfind(']')]
lead = []
for t in body.replace('\n', '').split(','):
    t = t.strip()
    if not t or t == '1':
        if t == '1': lead.append(tuple([0] * n))
        continue
    e = [0] * n
    for fac in t.split('*'):
        m = re.fullmatch(r'(w\d+)(?:\^(\d+))?', fac.strip())
        if not m: raise ValueError(f'bad monomial {t!r}')
        e[vidx[m.group(1)]] += int(m.group(2) or 1)
    lead.append(tuple(e))
LT = np.array(lead, dtype=np.int16)
degs = LT.sum(axis=1)
def std_count(k):
    """number of degree-k monomials not divisible by any leading monomial"""
    cnt = 0
    L = LT[degs <= k]
    for c in itertools.combinations(range(k + n - 1), n - 1):
        e = np.zeros(n, dtype=np.int16); prev = -1
        for j, x in enumerate(c):
            e[j] = x - prev - 1; prev = x
        e[n - 1] = k + n - 2 - prev
        if not np.any(np.all(L <= e, axis=1)):
            cnt += 1
    return cnt
HF = []
for k in range(maxdeg + 1):
    from math import comb
    if comb(k + n - 1, n - 1) > 3_000_000:
        break
    HF.append(std_count(k))
res = {'nvars': n, 'n_lead_monomials': len(lead), 'max_lead_degree': int(degs.max()), 'HF': HF,
       'vanishes_from_degree': next((k for k, v in enumerate(HF) if v == 0), None)}
# Singular cross-check
sing = f"ring r = 32003, ({','.join(vars_)}), dp;\nideal I = {','.join('*'.join(f'{vars_[i]}^{e[i]}' for i in range(n) if e[i]) or '1' for e in lead)};\nI = std(I);\nintvec h = hilb(I, 1);\nprint(h);\nprint(hilb(I,2));\nquit;\n"
try:
    r = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=120)
    res['singular_stdout_tail'] = r.stdout.strip().split('\n')[-6:]
except Exception as ex:
    res['singular_error'] = repr(ex)
print(json.dumps(res))
