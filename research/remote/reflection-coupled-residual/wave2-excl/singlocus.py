"""Singular locus of Z = X cap V(D) (projective) mod p via Singular: I = (F, D, 2x2 minors of [grad F; grad D]).
dim(I) is the affine dimension of the singular cone (0 => Z smooth projectively, given Z nonempty).
Usage: python3 singlocus.py <tag> '<D expr>' [prime]"""
import sys, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
tag, expr = sys.argv[1], sys.argv[2]
p = int(sys.argv[3]) if len(sys.argv) > 3 else 32003
D = core.reduce_mod(core.parse(expr), p); Fp = core.reduce_mod(L.F, p)
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
def s(P):
    return ' + '.join(f'{v}*' + '*'.join(f'{V[i]}^{e[i]}' for i in range(6) if e[i]) if any(e) else str(v) for e, v in P.items()) or '0'
sing = f"""ring r = {p}, ({','.join(V)}), dp;
poly F = {s(Fp)};
poly D = {s(D)};
ideal I = F, D, minor(jacob(ideal(F, D)), 2);
ideal J = std(I);
print("DIMSING"); print(dim(J)); print("MULT"); print(mult(J));
ideal Z = std(ideal(F, D)); print("DIMZ"); print(dim(Z)); print("DEGZ"); print(mult(Z));
quit;
"""
t0 = time.time()
r = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=3000)
toks = r.stdout.split()
def grab(k):
    try: return int(toks[toks.index(k) + 1])
    except Exception: return None
res = {'exp': 'singular-locus', 'family': tag, 'D': expr[:120] + ('...' if len(expr) > 120 else ''), 'prime': p,
       'dim_sing_cone': grab('DIMSING'), 'mult_sing': grab('MULT'), 'dim_Z_cone': grab('DIMZ'), 'deg_Z': grab('DEGZ'),
       'seconds': round(time.time() - t0, 1), 'stderr': r.stderr[-300:]}
res['note'] = ('Z = X cap V(D) is SMOOTH mod p (singular cone = origin) => smooth over Q => no (3,1)-support on Z (book ch.3 s3.1)'
               if res['dim_sing_cone'] == 0 else f'singular locus of projective dimension {res["dim_sing_cone"] - 1 if res["dim_sing_cone"] is not None else None}')
L.receipt('singlocus', **res)
