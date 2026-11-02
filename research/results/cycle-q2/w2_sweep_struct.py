"""Structured semi-invariant candidates tested with the fast modular sector solver (eq.test_D_fast, p = 4194301 and a
second prime), D and D^2 where affordable: Hessian 2x2 minors (degree 8), 3x3 minors (degree 12), subset-Euler forms
E_I(F) = sum_{i in I} Z_i F_i (degree 6) and their squares (degree 12), block binary cubics in (X_i^2, Y_i^2) with the
block-q2 coefficients, products of partials F_i F_j not yet tested.
Usage: python3 w2_sweep_struct.py [nproc] [which: minors2,minors3,euler,blocks,all]
Skips candidates already present (normalized) in log.jsonl / wave2-search.log.jsonl."""
import sys, time, json, itertools, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import multiprocessing as mp
from fractions import Fraction
import w2lib as W
import eq, modlu
from hodge import core, jacobian

LANE = '/tmp/hodge/results/cycle-q2'
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']


def det_poly(M):
    """exact determinant of a square matrix of dict polys (Laplace, small sizes)"""
    n = len(M)
    if n == 1:
        return M[0][0]
    out = {}
    for j in range(n):
        if not M[0][j]:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        sub = det_poly(minor)
        term = core.mul(M[0][j], sub)
        out = core.add(out, term) if j % 2 == 0 else core.sub(out, term)
    return out


def candidates(which):
    F, A, q, rec = core.load_instance('cycle-q2')
    Fi = [core.pderiv(F, i) for i in range(6)]
    H = [[core.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
    Z = [core.parse(v) for v in V]
    out = []
    if which in ('minors2', 'all'):
        for I in itertools.combinations(range(6), 2):
            for J in itertools.combinations(range(6), 2):
                if J < I:
                    continue
                M = [[H[i][j] for j in J] for i in I]
                D = det_poly(M)
                if D:
                    out.append((f'hess-minor2', f'det H[{",".join(V[i] for i in I)};{",".join(V[j] for j in J)}]', D, False))
    if which in ('minors3', 'all'):
        for I in itertools.combinations(range(6), 3):
            for J in itertools.combinations(range(6), 3):
                if J < I:
                    continue
                M = [[H[i][j] for j in J] for i in I]
                D = det_poly(M)
                if D:
                    out.append((f'hess-minor3', f'det H[{",".join(V[i] for i in I)};{",".join(V[j] for j in J)}]', D, False))
    if which in ('euler', 'all'):
        for r in range(1, 7):
            for I in itertools.combinations(range(6), r):
                D = {}
                for i in I:
                    D = core.add(D, core.mul(Z[i], Fi[i]))
                out.append(('subset-euler', 'E_{' + ','.join(V[i] for i in I) + '}(F)', D, True))
    if which in ('blocks', 'all'):
        # block-q2 binary cubic P(u,w) = u^3 + 2u^2w - uw^2 + w^3 in (X_i^2, Y_i^2), and all 3 pairs; also u^3 - w^3, u^3 + w^3 variants
        for i in range(3):
            u = core.parse(f'X{i}**2'); w = core.parse(f'Y{i}**2')
            P = core.add(core.add(core.power(u, 3), core.scale(core.mul(core.power(u, 2), w), Fraction(2))),
                         core.add(core.scale(core.mul(u, core.power(w, 2)), Fraction(-1)), core.power(w, 3)))
            out.append(('block-cubic', f'P_blockq2(X{i}^2,Y{i}^2)', P, True))
            out.append(('block-cubic', f'X{i}^6 - Y{i}^6', core.sub(core.power(u, 3), core.power(w, 3)), True))
        # pair products of two different pair-Fermat forms (degree 12)
        for i, j in itertools.combinations(range(3), 2):
            D = core.mul(core.parse(f'X{i}**6 + Y{i}**6'), core.parse(f'X{j}**6 + Y{j}**6'))
            out.append(('pair-fermat-product', f'(X{i}^6+Y{i}^6)(X{j}^6+Y{j}^6)', D, False))
    if which in ('fifj', 'all'):
        for i in range(6):
            for j in range(i, 6):
                out.append(('F_i*F_j', f'F_{V[i]}*F_{V[j]}', core.mul(Fi[i], Fi[j]), False))
    return out


def done_set():
    done = set()
    for fn in ('log.jsonl', 'wave2-search.log.jsonl'):
        p = os.path.join(LANE, fn)
        if not os.path.exists(p):
            continue
        for line in open(p):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if 'D_norm' in r:
                done.add(r['D_norm'])
            elif 'D' in r and isinstance(r['D'], str):
                done.add(r['D'])
    return done


def work(args):
    fam, note, Dstr, square = args
    F, A, q, rec = core.load_instance('cycle-q2')
    D = core.parse(Dstr)
    res = []
    t0 = time.time()
    try:
        r = eq.test_D_fast(F, A, q, D, fam, note=note + ' [D] wave2', verbose=False, do_e1=True)
        rec1 = {'step': 'struct-sweep', 'family': fam, 'note': note, 'D': core.to_str(D), 'd': core.degree(D), 'e1': r.get('e1'),
                'consistent': r.get('consistent'), 'rank': r.get('rank'), 'nullity': r.get('nullity'), 'cols': r.get('cols'), 'order': r.get('order'),
                'p': r.get('p'), 'seconds': r.get('seconds')}
        if r.get('e1') is not False and r.get('consistent') is not None:
            rb = eq.test_D_fast(F, A, q, D, fam, note=note + ' [D] wave2 2nd prime', verbose=False, do_e1=False, p=4194287)
            rec1['consistent_p2'] = rb.get('consistent'); rec1['rank_p2'] = rb.get('rank')
        W.log(**rec1)
        res.append(rec1)
        if square and core.degree(D) <= 6:
            D2 = core.mul(D, D)
            r2 = eq.test_D_fast(F, A, q, D2, fam + '^2', note=note + ' [D^2] wave2', verbose=False, do_e1=True, Dstr='(' + core.to_str(D) + ')**2')
            rec2 = {'step': 'struct-sweep', 'family': fam + '^2', 'note': note + ' squared', 'D': '(' + core.to_str(D) + ')**2', 'd': core.degree(D2),
                    'e1': r2.get('e1'), 'consistent': r2.get('consistent'), 'rank': r2.get('rank'), 'nullity': r2.get('nullity'), 'cols': r2.get('cols'),
                    'order': r2.get('order'), 'p': r2.get('p'), 'seconds': r2.get('seconds')}
            W.log(**rec2)
            res.append(rec2)
    except Exception as ex:
        W.log(step='struct-sweep', family=fam, note=note + f' ERROR {ex}', D=Dstr)
    return res


if __name__ == '__main__':
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    which = sys.argv[2] if len(sys.argv) > 2 else 'all'
    cands = candidates(which)
    done = done_set()
    jobs = []
    for fam, note, D, square in cands:
        Ds = core.to_str(D)
        if Ds in done:
            continue
        jobs.append((fam, note, Ds, square))
    print(f'[struct-sweep] {len(cands)} candidates, {len(jobs)} new', flush=True)
    W.log(step='struct-sweep-start', which=which, n_candidates=len(cands), n_new=len(jobs))
    t0 = time.time()
    n_cons = 0
    with mp.get_context('fork').Pool(nproc) as pool:
        for res in pool.imap_unordered(work, jobs):
            for r in res:
                flag = 'CONSISTENT!!!' if r.get('consistent') else ''
                print(f"  {r.get('family')} d={r.get('d')} e1={r.get('e1')} cons={r.get('consistent')} rk={r.get('rank')} {r.get('seconds')}s {r.get('note')[:60]} {flag}", flush=True)
                if r.get('consistent'):
                    n_cons += 1
    W.log(step='struct-sweep-end', which=which, n_new=len(jobs), n_consistent=n_cons, seconds=round(time.time() - t0, 1))
    print(f'[struct-sweep] done {time.time()-t0:.0f}s consistent={n_cons}', flush=True)
