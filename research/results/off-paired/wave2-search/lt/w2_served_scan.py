#!/usr/bin/env python3
"""Exact (mod p) served-space scan: for structured linear families of denominators (forms in 3/4/5-variable subsets of degree
4-6 in each G-character; random colon elements at d = 9..12), compute the sector-(3,0) served target space of random members
and report any nontrivial (mod J) served class — the first such class on this variety would be new, and by the sum principle
products of such denominators serve the span.   python3 w2_served_scan.py --workers 4 --job servedscan"""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
import argparse, json, time, itertools, random
import multiprocessing as mp
from fractions import Fraction
import w2lib, oplib
from hodge import core, jacobian
from hodge.core import NV, add, scale, mul
ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=4); ap.add_argument('--job', default='servedscan')
ap.add_argument('--members', type=int, default=3); ap.add_argument('--all-sectors', action='store_true'); ap.add_argument('--skip-colon', action='store_true')
args = ap.parse_args()
F, A, q = w2lib.load()
rng = random.Random(12345)
tasks = []
# (a) variable-subset forms
for k in (3, 4, 5):
    for V in itertools.combinations(range(NV), k):
        for deg in (4, 5, 6):
            fams = {}
            for m in core.monomials(deg):
                if all(m[i] == 0 for i in range(NV) if i not in V):
                    fams.setdefault(oplib.char(m), []).append(m)
            for c, mons in fams.items():
                if len(mons) < 2:
                    continue
                for j in range(args.members):
                    P = {m: Fraction(rng.randint(-3, 3)) for m in mons}
                    P = {m: v for m, v in P.items() if v}
                    if len(P) < 2: continue
                    tasks.append((f'vars{"".join(core.VARS[i] for i in V)}-deg{deg}-c{c[0]}{c[1]}-r{j}', P))
                if deg == 6 and k == 3:
                    P = {tuple(6 if t == i else 0 for t in range(NV)): Fraction(1) for i in V}
                    if oplib.char_of(P) is not None:
                        tasks.append((f'fermat-cone-{"".join(core.VARS[i] for i in V)}', P))
# (b) named
named = {'f6': 'X0**6+X1**6+X2**6', 'g6': 'Y0**6+Y1**6+Y2**6', 'H': 'X0**6+X1**6+Y0**6+Y1**6+X0**3*X1*Y0*Y1/2',
         'C': 'X1**3*Y1/2 + X0**2*X1*Y0/3', 'D0': 'X0**4 + X0*X1*Y0*Y1/4 + X1*X2*Y0*Y2/9', 'X2^6+Y2^6': 'X2**6+Y2**6', 'X2^6-Y2^6': 'X2**6-Y2**6',
         'Qx': 'X0**3*X1*Y0*Y1/2 + X1**3*X2*Y1*Y2/2 + X0**2*X1*X2*Y0*Y2/3', 'R': 'X0**3*Y0*Y1/2 + X0**2*X2*Y0*Y2/3 + X1**2*X2*Y1*Y2/2'}
for n, e in named.items():
    tasks.append((n, core.parse(e)))
# (c) random colon elements d = 9..12 in (0,0) and (3,0)
for d in ((9, 10, 11, 12) if not args.skip_colon else ()):
    for c in ((0, 0), (3, 0)):
        try:
            basis, gauge, info = w2lib.build_family(F, A, q, d, c, 'colon')
        except Exception as ex:
            print('colon basis failed', d, c, ex, flush=True); continue
        for j in range(2):
            P = {}
            for B in basis:
                P = add(P, scale(B, Fraction(rng.randint(-2, 2))))
            if P: tasks.append((f'colon-d{d}-c{c[0]}{c[1]}-r{j}', P))
print(f'{len(tasks)} tasks', flush=True)
OUT = f'wave2-{args.job}.jsonl'
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try: done.add(json.loads(line)['name'])
        except Exception: pass
sectors = [(a, b) for a in (1, 3, 5) for b in range(6)] if args.all_sectors else [(3, 0)]

def run(item):
    name, D = item
    t0 = time.time()
    rec = {'name': name, 'd': core.degree(D), 'char': oplib.char_of(D), 'terms': len(D)}
    try:
        divF, _ = jacobian.divisible_by(F, D)
        if divF:
            rec['status'] = 'F-divides'; return rec
        res = {}
        for cA in sectors:
            sv = oplib.served_targets(F, q, D, p=oplib.P0, cA=cA)
            res[f'{cA[0]}{cA[1]}'] = [sv['dim_served'], sv['dim_J'], sv['dim_nontrivial']]
            if sv['dim_nontrivial'] > 0:
                rec.setdefault('nontrivial', {})[f'{cA[0]}{cA[1]}'] = [core.to_str(oplib.lift_poly_small(B, oplib.P0) or {}) for B in sv['basis']][:5]
        rec['served'] = res
        rec['status'] = 'NONTRIVIAL' if 'nontrivial' in rec else 'trivial'
    except Exception as ex:
        rec['status'] = 'error'; rec['error'] = str(ex)[:200]
    rec['time'] = round(time.time() - t0, 1)
    return rec

todo = [t for t in tasks if t[0] not in done]
T0 = time.time(); n_non = 0
with mp.get_context('fork').Pool(args.workers) as pool:
    for rec in pool.imap_unordered(run, todo, chunksize=4):
        w2lib.append_json(OUT, rec)
        if rec['status'] == 'NONTRIVIAL':
            n_non += 1; print('!!!! NONTRIVIAL SERVED CLASS:', json.dumps(rec), flush=True)
recs = [json.loads(l) for l in open(OUT)]
summ = {'ts': w2lib.now(), 'step': 'served-scan', 'job': args.job, 'n': len(recs), 'nontrivial': sum(1 for r in recs if r.get('status') == 'NONTRIVIAL'),
        'by_status': {s: sum(1 for r in recs if r.get('status') == s) for s in set(r.get('status') for r in recs)},
        'sectors': [f'{a}{b}' for a, b in sectors], 'wall': round(time.time() - T0, 1)}
w2lib.append_json('wave2-search.remote.log.jsonl', summ); print(json.dumps(summ), flush=True)
