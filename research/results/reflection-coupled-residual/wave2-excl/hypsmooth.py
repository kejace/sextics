"""Smoothness of the hyperplane sections X cap {z_i = 0} mod p (Singular, dim of the Jacobian ideal of F|_{z_i=0} in 5 variables).
Smooth mod p => smooth over Q (the singular scheme is defined over Z; emptiness mod p of its fibre implies emptiness over Q by the rank
argument on the Macaulay matrices of the Jacobian ideal). Usage: python3 hypsmooth.py [prime]"""
import sys, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
p = int(sys.argv[1]) if len(sys.argv) > 1 else 32003
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
Fp = core.reduce_mod(L.F, p)
for i in range(6):
    Fi0 = core.substitute_scalar(Fp, i, 0, p)
    rest = [v for k, v in enumerate(V) if k != i]
    def s(P):
        return '+'.join(f'{v}*' + '*'.join(f'{V[k]}^{e[k]}' for k in range(6) if e[k]) if any(e) else str(v) for e, v in P.items()) or '0'
    sing = f"ring r = {p}, ({','.join(rest)}), dp;\npoly G = {s(Fi0)};\nideal I = G, jacob(G);\nideal J = std(I);\nprint(\"DIM\"); print(dim(J));\nquit;\n"
    t0 = time.time()
    r = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=1200)
    toks = r.stdout.split()
    try: dimJ = int(toks[toks.index('DIM') + 1])
    except Exception: dimJ = None
    L.receipt('hypsmooth', exp='hyperplane-section-smoothness', section=f'{V[i]}=0', prime=p, dim_singular_cone=dimJ,
              smooth=(dimJ == 0), seconds=round(time.time() - t0, 1),
              note='X cap {%s=0} is a smooth sextic threefold mod p (hence over Q): its (3,1)-served space is zero (Lefschetz), so by the sum principle a denominator divisible by %s serves exactly what its cofactor serves' % (V[i], V[i]) if dimJ == 0 else 'section singular: no reduction')
