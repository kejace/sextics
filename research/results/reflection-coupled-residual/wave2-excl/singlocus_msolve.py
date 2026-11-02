"""Singular locus of Z = X cap V(D) via msolve (-g 1 leading ideal) + Singular dim on the monomial lead ideal.
Usage: python3 singlocus_msolve.py <tag> '<D expr>' [prime] [threads]"""
import sys, json, time, subprocess, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
from hodge.core import NV, pderiv, mul, sub
tag, expr = sys.argv[1], sys.argv[2]
p = int(sys.argv[3]) if len(sys.argv) > 3 else 2147483587
thr = sys.argv[4] if len(sys.argv) > 4 else '2'
D = core.reduce_mod(core.parse(expr), p); Fp = core.reduce_mod(L.F, p)
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
def s(P):
    return '+'.join(f'{v}*' + '*'.join(f'{V[i]}^{e[i]}' for i in range(6) if e[i]) if any(e) else str(v) for e, v in P.items()) or '0'
Fi = [pderiv(Fp, i, p) for i in range(NV)]; Di = [pderiv(D, i, p) for i in range(NV)]
polys = [Fp, D]
for i, j in itertools.combinations(range(NV), 2):
    m = sub(mul(Fi[i], Dj := Di[j], p), mul(Fi[j], Di[i], p), p)
    if m: polys.append(m)
fn = f'{L.OUT}/singlocus_{tag}.ms'
open(fn, 'w').write(','.join(V) + f'\n{p}\n' + ',\n'.join(s(P) for P in polys) + '\n')
t0 = time.time()
r = subprocess.run(['msolve', '-f', fn, '-o', fn + '.lt', '-g', '1', '-t', thr, '-v', '1'], capture_output=True, text=True, timeout=7000)
txt = open(fn + '.lt').read(); body = txt[txt.find('[') + 1: txt.rfind(']')]
mons = [t.strip() for t in body.replace('\n', '').split(',') if t.strip()]
sing = f"ring r = {p}, ({','.join(V)}), dp;\nideal I = {','.join(mons)};\nattrib(I,\"isSB\",1);\nprint(\"DIMSING\"); print(dim(I)); print(\"MULT\"); print(mult(I));\nquit;\n"
rs = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=600)
toks = rs.stdout.split()
def grab(k):
    try: return int(toks[toks.index(k) + 1])
    except Exception: return None
res = {'exp': 'singular-locus-msolve', 'family': tag, 'D': expr[:100], 'prime': p, 'n_generators': len(polys), 'lead_monomials': len(mons),
       'dim_sing_cone': grab('DIMSING'), 'mult_sing': grab('MULT'), 'seconds': round(time.time() - t0, 1)}
res['note'] = ('Z = X cap V(D) is SMOOTH mod p (singular cone = origin): no non-rational singularities, no non-normality => no (3,1)-support'
               if res['dim_sing_cone'] == 0 else f'singular locus has projective dimension {res["dim_sing_cone"] - 1 if res["dim_sing_cone"] is not None else None}')
L.receipt('singlocus', **res)
