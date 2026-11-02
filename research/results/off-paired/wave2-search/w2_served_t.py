#!/usr/bin/env python3
"""Exact (mod p) served-target spaces along F_t = Fermat + t*Qx for the Fermat serving denominators f6, g6, f6*X0^2, f6+g6:
dim of nontrivial served classes in the target sector (3,0) as a function of t.  Cheap, deterministic."""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
import json, time
from fractions import Fraction
import w2lib, oplib
from hodge import core
from hodge.core import add, sub, scale, mul
F, A, q = w2lib.load()
f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6'); F0 = add(f6, g6); Qx = sub(F, F0)
cands = {'f6': f6, 'g6': g6, 'f6+g6': F0, 'f6*X0^2': mul(f6, core.parse('X0**2')), 'f6*X2*Y2': mul(f6, core.parse('X2*Y2')),
         'g6*X0^2': mul(g6, core.parse('X0**2'))}
ts = [Fraction(0), Fraction(1, 1000), Fraction(1, 100), Fraction(1, 10), Fraction(1, 2), Fraction(1)]
out = []
for name, D in cands.items():
    for t in ts:
        Ft = add(F0, scale(Qx, t))
        t0 = time.time()
        try:
            sv = oplib.served_targets(Ft, q, D, p=oplib.P0, cA=(3, 0))
            rec = {'D': name, 't': str(t), 'd': core.degree(D), 'dim_served': sv['dim_served'], 'dim_J': sv['dim_J'], 'dim_nontrivial': sv['dim_nontrivial'],
                   'rows': sv['rows'], 'cols': sv['cols'], 'time': round(time.time() - t0, 1)}
        except Exception as ex:
            rec = {'D': name, 't': str(t), 'error': str(ex)[:200]}
        print(json.dumps(rec), flush=True); out.append(rec)
        w2lib.append_json('wave2-served_t.jsonl', rec)
w2lib.append_json('wave2-search.remote.log.jsonl', {'ts': w2lib.now(), 'step': 'served-deformation', 'job': 'servedt',
                  'results': [(r['D'], r['t'], r.get('dim_nontrivial')) for r in out]})
