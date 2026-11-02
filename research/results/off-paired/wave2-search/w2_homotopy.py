#!/usr/bin/env python3
"""Deformation tracking F_t = F0 + t*Qx (F0 = Fermat, Qx = mixed part of off-paired) in the FULL monomial sector
family of degree d, character c (no colon restriction, which depends on t).  Start from the Fermat solution
D_0 = f6 * m (m = X0^(d-6)) or from the best random LM solution at t = 0, then warm-start LM at each t.

  python3 w2_homotopy.py --d 8 --c 0,0 --steps 40 --job hom8 [--start f6]   (writes wave2-<job>.jsonl)
"""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import argparse, json, time, math
import numpy as np
import w2lib
from hodge import core
from hodge.core import add, sub, scale, mul
from fractions import Fraction
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, default=8); ap.add_argument('--c', default='0,0'); ap.add_argument('--steps', type=int, default=40)
ap.add_argument('--job', required=True); ap.add_argument('--start', default='f6'); ap.add_argument('--tmax', type=float, default=1.0)
ap.add_argument('--maxiter', type=int, default=200); ap.add_argument('--restarts0', type=int, default=20)
args = ap.parse_args()
c = tuple(int(x) for x in args.c.split(','))
F, A, q = w2lib.load()
f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6'); F0 = add(f6, g6); Qx = sub(F, F0)
out = f'wave2-{args.job}.jsonl'
mons = w2lib.monos_char(args.d, c)
fam_basis = [{m: Fraction(1)} for m in mons]
def family_at(t):
    Ft = add(F0, scale(Qx, Fraction(t).limit_denominator(10**6)))
    return w2lib.SectorVarPro(Ft, A, q, fam_basis, c, gauge=w2lib.gauge_polys(Ft, args.d, c)), Ft
fam, Ft = family_at(0.0)
print(f'family d={args.d} c={c}: nD={fam.nD} dimP={fam.dimP} rows={fam.R} N={fam.nN} Q={fam.nQ}', flush=True)
# start point
if args.start == 'f6':
    m = tuple([args.d - 6, 0, 0, 0, 0, 0])
    D0 = {tuple(a + b for a, b in zip(e, m)): v for e, v in f6.items()}
    Dc = np.zeros(fam.nD)
    for e, v in D0.items(): Dc[mons.index(e)] = float(v)
    a = fam.Borth @ Dc; u = fam.P.T @ a
    rel0 = fam.relres(u)
    print(f't=0 start f6*X0^{args.d-6}: relres {rel0:.3e}', flush=True)
else:
    best = None
    for s in range(args.restarts0):
        r = fam.lm_from(s, maxiter=args.maxiter)
        if best is None or r['relres'] < best['relres']: best = r
    u = best['u']; print(f't=0 random start: best relres {best["relres"]:.3e}', flush=True)
w2lib.append_json(out, {'t': 0.0, 'relres': float(fam.relres(u)), 'start': args.start})
# the u-coordinates refer to the gauge complement of F_t which changes with t; carry D coefficients instead
Dc = fam.D_coef(u)[0]
ts = np.linspace(0, args.tmax, args.steps + 1)[1:]
hist = []
for t in ts:
    fam, Ft = family_at(float(t))
    # project the carried D onto the new gauge complement
    a = fam.Borth @ Dc; u0 = fam.P.T @ a
    rel_start = fam.relres(u0)
    r = fam.lm_from(0, maxiter=args.maxiter, init=u0)
    Dc = r['Dc']
    top = sorted([(abs(v), core.to_str({mons[i]: Fraction(round(v, 4))})) for i, v in enumerate(Dc)], reverse=True)[:4]
    rec = {'t': float(t), 'relres_warm_start': float(rel_start), 'relres': float(r['relres']), 'iters': r['iters'], 'status': r['status'],
           'top_terms': [x[1] for x in top]}
    hist.append(rec); w2lib.append_json(out, rec)
    print(json.dumps(rec), flush=True)
w2lib.append_json('wave2-search.remote.log.jsonl', {'ts': w2lib.now(), 'step': 'homotopy', 'job': args.job, 'd': args.d, 'sector': list(c),
                  'start': args.start, 'steps': args.steps, 'relres_path': [(round(h['t'], 4), float(f"{h['relres']:.3e}")) for h in hist]})
