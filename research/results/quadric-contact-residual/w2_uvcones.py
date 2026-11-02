#!/usr/bin/env python3
"""Cone families adapted to the two reflections: denominators in C[u, <subset of X1,X2,Y1,Y2>] or
C[v, <subset>] with u = X0 - Y0, v = X0 + Y0 (cones with vertex containing p+ = [1:0:0:1:0:0] resp.
p- = [1:0:0:-1:0:0], the isolated fixed points of the reflections s0*sigma0 resp. s0).  These subspaces are
not tau_Y-stable, so sectors are taken for the index-2 subgroup H = <sigma0, sigma2, h1> (3-component
characters, no Y-parity).  For each (subspace, H-sector, degree): exact E1 kernel on the span, then the
full-E2 Koszul varpro search (H-sector restricted unknowns, mixed parity) when the kernel is nonzero.

  python3 w2_uvcones.py --dmin 8 --dmax 12 --restarts 40 --procs 8
"""
import sys, os, json, time, argparse, itertools
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ.setdefault('OMP_NUM_THREADS', '1')
from fractions import Fraction
import numpy as np
import multiprocessing as mp
import w2lib as W
import w2_varpro as V

W.set_sector_mode('H')
ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=8); ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--restarts', type=int, default=40); ap.add_argument('--procs', type=int, default=8)
ap.add_argument('--seed', type=int, default=31000); ap.add_argument('--tag', default='uvcone')
ap.add_argument('--time-limit', type=float, default=1500)
ap.add_argument('--min-vars', type=int, default=3)
a = ap.parse_args()

U = {(1, 0, 0, 0, 0, 0): Fraction(1), (0, 0, 0, 1, 0, 0): Fraction(-1)}
Vv = {(1, 0, 0, 0, 0, 0): Fraction(1), (0, 0, 0, 1, 0, 0): Fraction(1)}
OTHER = [1, 2, 4, 5]   # X1, X2, Y1, Y2


def subspace_polys(d, w, subset):
    """Degree-d monomials in (w, Z_j for j in subset) expanded; grouped by H-character."""
    groups = {}
    wp = [{(0,) * 6: Fraction(1)}]
    for k in range(1, d + 1):
        wp.append(W.mul(wp[-1], w))
    for e in W.monomials(d):
        if e[3] != 0 or any(e[j] and j not in subset for j in OTHER) or (e[0] and False):
            continue
        # e[0] = exponent of w
        rest = {(0, e[1], e[2], 0, e[4], e[5]): Fraction(1)}
        P = W.mul(wp[e[0]], rest)
        c = W.chi_of(P)
        groups.setdefault(c, []).append(P)
    return groups


n_fam = 0
for d in range(a.dmin, a.dmax + 1):
    t0 = time.time()
    dims = {}
    for wname, w in (('u', U), ('v', Vv)):
        for nv in range(a.min_vars - 1, 5):
            for subset in itertools.combinations(OTHER, nv):
                name = wname + ''.join(W.core.VARS[j] for j in subset)
                groups = subspace_polys(d, w, set(subset))
                for c, polys in sorted(groups.items()):
                    ker, info = W.e1_family_span(d, polys, c)
                    key = f'{name}:{W.chi_str(c)}'
                    dims[key] = {'K': len(polys), 'kernel': len(ker)}
                    if not ker:
                        continue
                    n_fam += 1
                    print(f'=== {key} d={d}: span {len(polys)} kernel {len(ker)}', flush=True)
                    # members with least-squares lifts (H-sector, mixed parity)
                    LS = W.LiftSolver(d + 6, c, None)
                    members = []
                    maxres = 0.0
                    for k, E in enumerate(ker):
                        Nk, res = LS.lift(W.mul(W.Af, W.fpoly(E)))
                        maxres = max(maxres, res)
                        members.append((W.fpoly(E), Nk, f'E:{k}'))
                    tb = time.time()
                    try:
                        vp = W.KoszulVarPro(members, d, c, parity=None, verbose=True)
                    except Exception as ex:
                        import traceback; traceback.print_exc()
                        W.log_line({'task': 'varpro', 'd': d, 'chi': W.chi_str(c), 'vars': name, 'tag': a.tag, 'status': 'family-error', 'error': str(ex)[:300]})
                        continue
                    rec = {'task': 'varpro', 'd': d, 'chi': 'H' + W.chi_str(c), 'vars': name, 'mode': 'span', 'tag': a.tag,
                           'dimS': len(polys), 'dim_kernel': len(ker), 'n_members': len(members), 'max_lift_res': maxres,
                           'rows': vp.rows, 'nK': vp.nK, 'nQ': vp.nQ, 'K': vp.K, 'dimP': vp.dimP, 'excluded': vp.excluded,
                           't_build': round(vp.build_time, 1)}
                    if vp.dimP == 0:
                        rec['status'] = 'only-gauge'; W.log_line(rec); continue
                    V.VP = vp
                    nres = max(12, min(a.restarts, 8 * vp.dimP))
                    jobs = [(a.seed + 7919 * i + 101 * n_fam, 1, 300, 1e-12) for i in range(nres)]
                    cands = []
                    ctx = mp.get_context('fork')
                    pool = ctx.Pool(processes=min(a.procs, nres))
                    ts = time.time()
                    try:
                        for res in pool.imap_unordered(V._worker, jobs):
                            cands.extend(res)
                            if time.time() - ts > a.time_limit:
                                break
                        pool.close()
                    finally:
                        pool.terminate(); pool.join()
                    cands.sort(key=lambda x: (x['spurious'], x['rel']))
                    rels = np.array([x['rel'] for x in cands if not x['spurious']])
                    rec.update({'restarts_done': len(cands), 'n_spurious': int(sum(1 for x in cands if x['spurious'])),
                                'best_rel': float(rels.min()) if rels.size else None, 'median_rel': float(np.median(rels)) if rels.size else None,
                                'n_rel_lt_1e-6': int((rels < 1e-6).sum()), 'n_rel_lt_1e-10': int((rels < 1e-10).sum()),
                                'hist_log10': {str(k): int(v) for k, v in zip(*np.unique(np.floor(np.log10(np.maximum(rels, 1e-300))).astype(int), return_counts=True))} if rels.size else {},
                                'seconds': round(time.time() - tb, 1), 'status': 'done'})
                    outp = f'{W.RES}/varpro-{a.tag}-d{d}-{name}-c{W.chi_str(c)}.json'
                    json.dump({'rec': rec, 'members': [m[2] for m in members], 'best': cands[:10]}, open(outp, 'w'))
                    rec['out'] = os.path.basename(outp)
                    recon = []
                    for x in cands:
                        if x['spurious'] or x['rel'] > 1e-10:
                            continue
                        att = W.try_reconstruct_and_confirm(vp, np.array(x['a']), verbose=True, tag=f'{a.tag}-d{d}-{name}')
                        recon.append({'rel': x['rel'], 'attempts': att})
                        if any(t.get('witness_path') for t in att):
                            rec['WITNESS'] = True; break
                    rec['recon'] = recon
                    W.log_line(rec)
                    print(json.dumps({k: v for k, v in rec.items() if k not in ('hist_log10', 'recon')}), flush=True)
    W.log_line({'task': 'e1-uvcones', 'd': d, 'n_subspaces': len(dims), 'n_nontrivial': sum(1 for v in dims.values() if v['kernel']),
                'nontrivial': {k: v for k, v in dims.items() if v['kernel']}, 'seconds': round(time.time() - t0, 1)})
print('uvcones done', n_fam)
