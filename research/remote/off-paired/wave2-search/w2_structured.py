#!/usr/bin/env python3
"""Structured fixed candidates not covered by wave 1: Hessian minors (3x3 deg 12, 4x4 deg 16, 5x5 deg 20),
mixed Fermat-type cones times blocks (deg 8-12).  E1 (sector-reduced, mod p) then full sector-reduced solve.
  python3 w2_structured.py --set hess3 --workers 8 --job struct-hess3
"""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import argparse, json, time, itertools
import multiprocessing as mp
from fractions import Fraction
import w2lib, oplib
from hodge import core, jacobian
from hodge.core import NV, mul, add, sub, scale, pderiv, shift, power
ap = argparse.ArgumentParser(); ap.add_argument('--set', required=True); ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--job', required=True); ap.add_argument('--dmax', type=int, default=20); ap.add_argument('--limit', type=int, default=100000)
args = ap.parse_args()
F, A, q = w2lib.load()
Fi = [pderiv(F, i) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]

def det(M):
    n = len(M)
    if n == 1: return M[0][0]
    out = {}
    for j in range(n):
        if not M[0][j]: continue
        minor = [row[:j] + row[j+1:] for row in M[1:]]
        term = mul(M[0][j], det(minor))
        out = add(out, term if j % 2 == 0 else scale(term, -1))
    return out

cands = []
if args.set.startswith('hess'):
    k = int(args.set[4])
    seen = set()
    for I in itertools.combinations(range(NV), k):
        for J in itertools.combinations(range(NV), k):
            key = tuple(sorted([I, J]))
            if key in seen: continue
            seen.add(key)
            M = [[Fij[i][j] for j in J] for i in I]
            D = det(M)
            if D: cands.append((f'hess{k}-minor I={I} J={J}', D))
elif args.set == 'cones':
    # mixed Fermat-type cones (G-eigenvectors) times monomials / blocks, degrees 8..12
    sixth = [core.parse(f'{v}**6') for v in core.VARS]
    cones = {}
    for trip in itertools.combinations(range(NV), 3):
        if not (any(i < 3 for i in trip) and any(i >= 3 for i in trip)): continue   # mixed X/Y only (pure ones excluded)
        name = '+'.join(core.VARS[i] + '^6' for i in trip)
        cones[name] = add(add(sixth[trip[0]], sixth[trip[1]]), sixth[trip[2]])
    # cubic-type curves X_a^3 Y_b^3 + Z^6 are not eigenvectors in general; add the (3,0)-type ones that are
    for a in range(3):
        for b in range(3, 6):
            for z in range(NV):
                if z in (a, b): continue
                P = add(mul(power({tuple(1 if t == a else 0 for t in range(NV)): Fraction(1)}, 3), power({tuple(1 if t == b else 0 for t in range(NV)): Fraction(1)}, 3)), sixth[z])
                if oplib.char_of(P) is not None: cones[f'{core.VARS[a]}^3{core.VARS[b]}^3+{core.VARS[z]}^6'] = P
    blocks = {'1': {(0,)*NV: Fraction(1)}}
    for deg in (2, 3, 4):
        for m in core.monomials(deg):
            blocks['m' + ''.join(str(x) for x in m)] = {m: Fraction(1)}
    D0 = core.parse('X0**4 + X0*X1*Y0*Y1/4 + X1*X2*Y0*Y2/9'); H = core.parse('X0**6+X1**6+Y0**6+Y1**6+X0**3*X1*Y0*Y1/2')
    C = core.parse('X1**3*Y1/2 + X0**2*X1*Y0/3'); XY = core.parse('X2*Y2'); Mp = core.parse('X2**6+Y2**6'); Mm = core.parse('X2**6-Y2**6')
    blocks.update({'D0': D0, 'H': H, 'C': C, 'X2Y2': XY, 'X2^6+Y2^6': Mp, 'X2^6-Y2^6': Mm})
    for cn, P in cones.items():
        for bn, B in blocks.items():
            D = mul(P, B)
            d = core.degree(D)
            if d < 8 or d > 12: continue
            if oplib.char_of(D) is None: continue
            if d == 8 and oplib.char_of(D)[1] != 0: continue     # degree 8: only (a,0) open
            cands.append((f'cone({cn})*{bn}', D))
elif args.set == 'file':
    for line in open(args.job + '.cands'):
        line = line.strip()
        if line: cands.append((line, core.parse(line)))
cands = [c for c in cands if core.degree(c[1]) <= args.dmax][:args.limit]
print(f'{len(cands)} candidates (set {args.set})', flush=True)
OUT = f'wave2-{args.job}.jsonl'
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try: done.add(json.loads(line)['name'])
        except Exception: pass

def run(item):
    name, D = item
    t0 = time.time(); d = core.degree(D)
    rec = {'name': name, 'd': d, 'char': oplib.char_of(D), 'terms': len(D)}
    try:
        divF, _ = jacobian.divisible_by(F, D)
        rec['F_divides'] = bool(divF)
        if divF: rec['status'] = 'spurious'; return rec
        if d <= 18:
            ok, info = oplib.e1_reduced(F, A, q, D); rec['e1'] = bool(ok); rec['e1_info'] = info
            if not ok: rec['status'] = 'e1-fail'; rec['time'] = round(time.time() - t0, 1); return rec
        else:
            rec['e1'] = True
        S = oplib.ReducedSystem(F, A, q, D, oplib.P0)
        rec.update({'rows': S.nrows, 'cols': S.ncols})
        if S.ncols > 9000:
            r = __import__('hodge.fixed_d', fromlist=['x']).solve_compressed(S)
        else:
            r = S.solve()
        rec.update({'consistent': bool(r['consistent']), 'rank': r['rank'], 'nullity': r['nullity']})
        rec['status'] = 'CONSISTENT' if r['consistent'] else 'inconsistent'
        if r['consistent']:
            st, data = oplib.lift_solution(F, A, q, D)
            rec['lift'] = st
            if st == 'solved':
                rec['witness'] = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
                json.dump(rec['witness'], open(f'witness-candidate-{args.job}.json', 'w')); json.dump(rec['witness'], open('witness-candidate.json', 'w'))
    except Exception as ex:
        rec['status'] = 'error'; rec['error'] = str(ex)[:300]
    rec['time'] = round(time.time() - t0, 1)
    return rec

todo = [c for c in cands if c[0] not in done]
T0 = time.time()
with mp.get_context('fork').Pool(args.workers) as pool:
    for rec in pool.imap_unordered(run, todo):
        w2lib.append_json(OUT, rec)
        oplib.log(rec['name'][:200], rec['d'], f'wave2-{args.set}', rec.get('e1'), rec.get('consistent', False), rec.get('rank', -1),
                  rec.get('nullity', -1), rec.get('time', 0), f"status={rec.get('status')} char={rec.get('char')}")
        print(json.dumps({k: v for k, v in rec.items() if k not in ('witness', 'e1_info')}), flush=True)
recs = [json.loads(l) for l in open(OUT)]
summ = {'ts': w2lib.now(), 'step': 'structured', 'job': args.job, 'set': args.set, 'n': len(recs),
        'e1_pass': sum(1 for r in recs if r.get('e1')), 'consistent': sum(1 for r in recs if r.get('consistent')),
        'by_status': {s: sum(1 for r in recs if r.get('status') == s) for s in set(r.get('status') for r in recs)},
        'degrees': sorted({r['d'] for r in recs}), 'wall': round(time.time() - T0, 1)}
w2lib.append_json(w2lib.REMOTE_LOG, summ); print(json.dumps(summ), flush=True)
