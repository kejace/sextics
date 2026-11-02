"""Hilbert function of a monomial (leading) ideal from an msolve '-g 1' output; usage: w2_hilb.py file.lm d parity"""
import sys, re, json, itertools, time
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from w2core import w2log
path, d, par = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
txt = open(path).read()
vars_ = re.search(r'#variable order:\s*(.*)', txt).group(1).replace(' ', '').split(',')
n = len(vars_); vidx = {v: i for i, v in enumerate(vars_)}
body = txt[txt.index('['):]
body = body[1:body.index(']')]
lms = []
for term in body.replace('\n', '').split(','):
    term = term.strip()
    if not term: continue
    e = [0] * n
    for fac in term.split('*'):
        m = re.match(r'(\w+)(?:\^(\d+))?$', fac.strip())
        e[vidx[m.group(1)]] += int(m.group(2) or 1)
    lms.append(tuple(e))
lmset = set(lms)
bydeg = {}
for e in lms: bydeg.setdefault(sum(e), 0); bydeg[sum(e)] += 1
def divisible(e):
    # is some LM a divisor of e? check all sub-monomials of degree >= 1 (small degrees only)
    deg = sum(e)
    for k in sorted(bydeg):
        if k > deg: break
        # enumerate divisors of e of degree k
        idxs = [i for i in range(n) for _ in range(e[i])]
        for comb in set(itertools.combinations(idxs, k)):
            f = [0] * n
            for i in comb: f[i] += 1
            if tuple(f) in lmset: return True
    return False
def mons(k):
    for c in itertools.combinations_with_replacement(range(n), k):
        e = [0] * n
        for i in c: e[i] += 1
        yield tuple(e)
HF = {}
t0 = time.time()
for k in range(0, 5):
    if k == 4 and n > 70: break
    cnt = sum(1 for e in mons(k) if not divisible(e))
    HF[k] = cnt
    if cnt == 0: break
res = {'batch': 'macaulay', 'd': d, 'parity': par, 'n_vars': n, 'n_leading_monomials': len(lms), 'LM_by_degree': bydeg,
       'hilbert_function': HF, 'rigorous_exclusion': any(v == 0 for v in HF.values()),
       'note': 'HF(k)=0 for some k means the relaxation quadrics span all degree-k forms mod p, hence over Q: no E1-admissible D of this degree/parity satisfies E2',
       'seconds': round(time.time() - t0, 1), 'file': path}
print(json.dumps(res), flush=True)
w2log(res)
