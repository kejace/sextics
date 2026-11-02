#!/usr/bin/env python3
"""Served-space scan (mod p, sector (3,0) and all Y-odd sectors) for the structurally motivated denominators of the
construction lens (Klein cover / elliptic fibration / rational genus-2 surface picture).  Reuses the search lens's
oplib.served_targets.  Output: wave2-construct-served.jsonl in this directory.  Usage: python3 c5_served.py [--heavy]"""
import os, sys, json, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); SEARCH = os.path.join(os.path.dirname(HERE), 'wave2-search')
sys.path.insert(0, SEARCH); sys.path.insert(0, '/home/kejace/devel/sextics-work'); os.chdir(HERE)
import w2lib, oplib
from hodge import core, jacobian
F, A, q = w2lib.load()
heavy = '--heavy' in sys.argv
m = 'X2*Y2'; C = '(X1**3*Y1/2 + X0**2*X1*Y0/3)'; H = '(X0**6+X1**6+Y0**6+Y1**6+X0**3*X1*Y0*Y1/2)'
psi_m = f'(3*({m})**4 - 3*{C}*({m})**2 - 6*{H}*({m}) - {C}**2/4)'      # 3-torsion of t^2 = 2m^3 - C m - H
psi_t = f'(3*({m})**4 - 3*{C}*({m})**2 + 6*{H}*({m}) - {C}**2/4)'      # tau-image
cands = {
 'block-generic-c6=2':        'X2**6 - 2*Y2**6',                       # pullback of a rational 3-section s=const of Sigma
 'block-section-mixed':       'X2**6 + 2*X2**3*Y2**3 - 3*Y2**6',        # another curve in Sigma (v linear in H+Cm, m^3)
 'vert-twist-C2-mu4-even':    f'{C}**2 - (X0**2+Y1**2)**4',             # Lambda_mu u g Lambda_mu, mu = X0^2+Y1^2 (sqrt kappa on Lambda)
 'vert-twist-C2-mu4-odd':     f'{C}**2 - (X0*Y1+X1*Y0)**4',
 'vert-twist-C2-X0^8':        f'{C}**2 - X0**8',
 'vert-cube-H2-nu6':          f'{H}**2 - (X0**2+Y1**2)**6',             # Lambda: H = nu^3 (cube root of kappa), with g-image
 'mirror3xH':                 f'(X2**3+Y2**3)*(X2**3-Y2**3)*{H}',      # 2-torsion divisors of both Weierstrass fibrations x H
}
if heavy:
    cands['torsion3-prod'] = f'{psi_m}*{psi_t}'                      # degree 16, pullback of the 3-torsion multisections
    cands['vert-twist-x-mirror'] = f'({C}**2 - (X0**2+Y1**2)**4)*(X2**6-Y2**6)'  # degree 14
OUT = 'wave2-construct-served.jsonl'
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try: done.add(json.loads(line)['name'])
        except Exception: pass
sectors = [(a, b) for a in (1, 3, 5) for b in range(6)]
for name, expr in cands.items():
    if name in done: continue
    t0 = time.time(); D = core.parse(expr)
    rec = {'name': name, 'expr': expr, 'd': core.degree(D), 'char': oplib.char_of(D), 'parity': core.parity_of(D)}
    try:
        divF, _ = jacobian.divisible_by(F, D)
        if divF: rec['status'] = 'F-divides'
        else:
            rec['e1'] = bool(__import__('hodge.fixed_d', fromlist=['e1_solvable']).e1_solvable(F, A, q, D, oplib.P0))
            res = {}
            for cA in ([(3, 0)] if rec['d'] > 12 else sectors):
                sv = oplib.served_targets(F, q, D, p=oplib.P0, cA=cA)
                res[f'{cA[0]}{cA[1]}'] = [sv['dim_served'], sv['dim_J'], sv['dim_nontrivial']]
                if sv['dim_nontrivial'] > 0:
                    rec.setdefault('nontrivial', {})[f'{cA[0]}{cA[1]}'] = [core.to_str(oplib.lift_poly_small(B, oplib.P0) or {}) for B in sv['basis']][:5]
            rec['served'] = res; rec['status'] = 'NONTRIVIAL' if 'nontrivial' in rec else 'trivial'
    except Exception as ex:
        rec['status'] = 'error'; rec['error'] = str(ex)[:300]
    rec['time'] = round(time.time() - t0, 1); rec['ts'] = w2lib.now()
    with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)
