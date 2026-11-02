"""Probe: hyperplane sections through the plane {V=0} that are TANGENT to X.
H_lambda = {lambda . U = 0}: X ∩ H_lambda is singular exactly at the points (U0:0) where the line lambda.U=0 is tangent
to the plane sextic C = {f=0}; there the local equation is t^k + f(V) = 0 (k = contact order >= 2), a non-rational
3-fold singularity whose weighted blow-up has exceptional divisor {t^2 = f(V)} = the K3 double plane S_2 itself (k=2).
We take F_p-points of C, their tangent lines, and compute the K3-class projection of served(lambda.U) mod p.
usage: tangent_probe.py [p] [ntries]"""
import sys, json, random, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from flint import nmod_poly
from hodge import core
import servedK3
p = int(sys.argv[1]) if len(sys.argv) > 1 else 2147483647
ntries = int(sys.argv[2]) if len(sys.argv) > 2 else 3
uv = json.load(open('/tmp/hodge/results/polarized-k3-remainder/uv_forms.json'))
f = core.parse(uv['f'], p)      # ternary sextic in slots 0..2
fi = [core.pderiv(f, i, p) for i in range(3)]
random.seed(7)
found = 0; tries = 0
while found < ntries and tries < 200:
    tries += 1
    z1, z2 = random.randrange(p), random.randrange(p)
    # polynomial in z0
    coeffs = [0] * 7
    for e, c in f.items():
        coeffs[e[0]] = (coeffs[e[0]] + c * pow(z1, e[1], p) * pow(z2, e[2], p)) % p
    poly = nmod_poly(coeffs, p)
    roots = poly.roots()
    if not roots: continue
    z0 = int(roots[0][0])
    pt = (z0, z1, z2)
    lam = [core.evaluate(fi[i], (pt[0], pt[1], pt[2], 0, 0, 0), p) for i in range(3)]
    assert core.evaluate(f, (pt[0], pt[1], pt[2], 0, 0, 0), p) == 0
    assert sum(l * c for l, c in zip(lam, pt)) % p == 0      # Euler: tangent line passes through the point
    Dexpr = ' + '.join(f'{lam[i]}*X{i}' for i in range(3) if lam[i])
    res = servedK3.run(f'tangentU-pt{found}', Dexpr, p, verbose=False)
    res['point'] = pt; res['lambda'] = lam
    print(json.dumps(res), flush=True)
    with open('tangent_probe.jsonl', 'a') as fh: fh.write(json.dumps(res) + '\n')
    # control: a random (non-tangent) line through the plane V=0
    lam2 = [random.randrange(1, p) for _ in range(3)]
    res2 = servedK3.run(f'randomU-{found}', ' + '.join(f'{lam2[i]}*X{i}' for i in range(3)), p, verbose=False)
    print(json.dumps(res2), flush=True)
    with open('tangent_probe.jsonl', 'a') as fh: fh.write(json.dumps(res2) + '\n')
    found += 1
print('done', found, 'tangent lines tested')
