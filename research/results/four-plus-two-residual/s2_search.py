"""VM search driver (agent vm-search) for four-plus-two-residual: random-restart sector VarPro over linear
denominator families, LBFGS (+LM polish) or LM from random starts; candidates with small residual are rationally
reconstructed and re-solved exactly.

usage: python3 s2_search.py --job NAME --family SPEC --degrees 6,7 [--sectors irr|all|"000;303"] [--restarts 20]
          [--batches 50] [--procs 4] [--threads 1] [--maxiter 300] [--seed0 1] [--method lbfgs|lm|lbfgs+lm] [--instance IID]
family SPEC: kernelJ | J | vars:X0,X1,Y0,Y1 | prod:<expr> | kernelJ-minus-vars:X0,X1,Y0,Y1 (admissible D with at least one
             X2/Y2 term: complement of the core-only subfamily, numerically orthogonal complement inside the family)
Outputs: wave2/s2-<job>.jsonl (one line per batch, resumable), wave2/s2-<job>-cands.jsonl, lane log wave2-search.log.jsonl."""
import os, sys, argparse, json, time, math
ap = argparse.ArgumentParser()
ap.add_argument('--job', required=True); ap.add_argument('--family', required=True); ap.add_argument('--degrees', required=True)
ap.add_argument('--sectors', default='irr'); ap.add_argument('--restarts', type=int, default=20); ap.add_argument('--batches', type=int, default=50)
ap.add_argument('--procs', type=int, default=4); ap.add_argument('--threads', type=int, default=1); ap.add_argument('--maxiter', type=int, default=300)
ap.add_argument('--seed0', type=int, default=1); ap.add_argument('--instance', default='four-plus-two-residual')
ap.add_argument('--method', default='lbfgs+lm'); ap.add_argument('--cand-rel', type=float, default=1e-8)
ap.add_argument('--lm-iter', type=int, default=60); ap.add_argument('--min-adm', type=int, default=1)
ap.add_argument('--max-dimP', type=int, default=100000)
args = ap.parse_args()
os.environ.setdefault('OPENBLAS_NUM_THREADS', str(args.threads)); os.environ.setdefault('OMP_NUM_THREADS', str(args.threads))
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
import numpy as np
import fplib
from hodge import core
from multiprocessing import get_context

LANE = fplib.LANE
F, A, q, G = fplib.load(args.instance); ch = fplib.Chars(G)
JOBF = f'{LANE}/wave2/s2-{args.job}.jsonl'; CANDF = f'{LANE}/wave2/s2-{args.job}-cands.jsonl'
done = set()
if os.path.exists(JOBF):
    for l in open(JOBF):
        try:
            r = json.loads(l); done.add((r['d'], tuple(r['sector']), r['batch']))
        except Exception:
            pass


def parse_sectors(d):
    if args.sectors == 'irr':
        return [c for c in sorted(ch.all_chars(d)) if fplib.is_irreducible_sector(ch, d, c)]
    if args.sectors == 'all':
        return sorted(ch.all_chars(d))
    out = []
    for s in args.sectors.split(';'):
        out.append(tuple(int(x) for x in s.strip().split(',')) if ',' in s else tuple(int(x) for x in s.strip()))
    return out


def build_family(d, c):
    B, Dm, info = fplib.sector_basis_numeric(F, A, q, d, c, ch)
    idx = {m: i for i, m in enumerate(Dm)}
    spec = args.family
    note = ''
    if B.shape[0] == 0:
        return None, Dm, info, 'no admissible D'
    if spec == 'kernelJ':
        fam = B
    elif spec == 'J':
        Jv = fplib.jacobian_vecs(F, d, c, ch, Dm)
        fam = fplib.orth_rows(Jv)
    elif spec.startswith('vars:'):
        vs = [core.VARS.index(v) for v in spec.split(':', 1)[1].split(',')]
        sel = np.array([all(m[i] == 0 for i in range(6) if i not in vs) for m in Dm])
        sub = np.eye(len(Dm))[sel]
        fam = fplib.restrict_rows(B, sub)
    elif spec.startswith('kernelJ-minus-vars:'):
        vs = [core.VARS.index(v) for v in spec.split(':', 1)[1].split(',')]
        sel = np.array([all(m[i] == 0 for i in range(6) if i not in vs) for m in Dm])
        sub = np.eye(len(Dm))[sel]
        inside = fplib.restrict_rows(B, sub)
        # orthogonal complement of `inside` within span(B)
        if inside.shape[0]:
            Bc = B - (B @ inside.T) @ inside
            fam = fplib.orth_rows(Bc)
        else:
            fam = B
    elif spec.startswith('prod:'):
        P = core.parse(spec.split(':', 1)[1]); dp = core.degree(P); cp = ch.char_of(P)
        if cp is None or dp > d:
            return None, Dm, info, 'factor not an eigenvector or too big'
        cm = ch.cadd(c, ch.cneg(cp))
        rowsv = [fplib.poly_to_vec(core.mul(P, {m: 1}), idx, len(Dm)) for m in ch.monos_char(d - dp, cm)]
        if not rowsv:
            return None, Dm, info, 'no multiples in sector'
        fam = fplib.restrict_rows(B, np.array(rowsv))
    else:
        raise ValueError(spec)
    if fam.shape[0] == 0:
        return None, Dm, info, 'empty family'
    return fam, Dm, info, note


FAM = None


def run_restart(seed):
    fam = FAM
    t0 = time.time()
    out = {'seed': seed}
    if args.method in ('lbfgs', 'lbfgs+lm'):
        r = fam.minimize_from(seed, maxiter=args.maxiter)
        out.update({'relres_lbfgs': r['relres'], 'iters': r['iters']})
        u = r['u']; rel = r['relres']; Dc = r['Dc']
        if args.method == 'lbfgs+lm' and rel < 1e-2:
            r2 = fam.lm_from(seed, init=u, maxiter=args.lm_iter)
            if r2['relres'] < rel:
                rel, Dc, u = r2['relres'], r2['Dc'], r2['u']
            out.update({'lm_iters': r2['iters'], 'lm_status': r2['status']})
    else:
        r = fam.lm_from(seed, maxiter=args.maxiter)
        out.update({'iters': r['iters'], 'lm_status': r['status']})
        rel, Dc, u = r['relres'], r['Dc'], r['u']
    out['relres'] = rel; out['time'] = round(time.time() - t0, 2)
    if rel < args.cand_rel:
        out['Dc'] = Dc.tolist()
    return out


def hist(rels):
    h = {}
    for r in rels:
        k = f'<1e{int(math.ceil(math.log10(max(r, 1e-300))))}'
        h[k] = h.get(k, 0) + 1
    return h


for d in [int(x) for x in args.degrees.split(',')]:
    for c in parse_sectors(d):
        key = (d, tuple(c))
        if all((d, tuple(c), b) in done for b in range(args.batches)):
            continue
        t0 = time.time()
        fam_rows, Dm, info, note = build_family(d, c)
        if fam_rows is None or info['adm'] < args.min_adm:
            rec = {'batch': 'search', 'agent': 'vm-search', 'job': args.job, 'family': args.family, 'd': d, 'sector': list(c),
                   'adm': info.get('adm', 0), 'note': note or 'skipped', 'restarts': 0}
            fplib.append_json(JOBF, {**rec, 'batch_no': -1}); fplib.wlog(rec)
            print(rec, flush=True)
            continue
        FAM = fplib.SectorVarPro(F, A, q, fam_rows, Dm, c, ch, tag=args.job)
        if FAM.dimP == 0 or FAM.dimP > args.max_dimP:
            rec = {'batch': 'search', 'agent': 'vm-search', 'job': args.job, 'family': args.family, 'd': d, 'sector': list(c),
                   'adm': info['adm'], 'K': int(FAM.rankB), 'gauge': FAM.ng, 'dimP': FAM.dimP, 'note': 'dimP zero or too large', 'restarts': 0}
            fplib.append_json(JOBF, {**rec, 'batch_no': -1}); fplib.wlog(rec)
            print(rec, flush=True)
            continue
        base = {'batch': 'search', 'agent': 'vm-search', 'job': args.job, 'family': args.family, 'd': d, 'sector': list(c),
                'adm': info['adm'], 'sv_gap': info.get('sv_gap'), 'K': int(FAM.rankB), 'gauge': FAM.ng, 'dimP': FAM.dimP,
                'rows': FAM.R, 'nN': FAM.nN, 'nQ': FAM.nQ, 'method': args.method, 'maxiter': args.maxiter,
                'build_seconds': round(time.time() - t0, 1)}
        print('FAMILY', base, flush=True)
        ctx = get_context('fork')
        for b in range(args.batches):
            if (d, tuple(c), b) in done:
                continue
            tb = time.time()
            seeds = [args.seed0 + 100000 * b + i for i in range(args.restarts)]
            if args.procs > 1:
                with ctx.Pool(args.procs) as pool:
                    res = pool.map(run_restart, seeds)
            else:
                res = [run_restart(s) for s in seeds]
            rels = [r['relres'] for r in res]
            rels_s = sorted(rels)
            rec = dict(base); rec.update({'batch_no': b, 'seed0': seeds[0], 'restarts': len(res), 'best_rel': rels_s[0],
                                           'median_rel': rels_s[len(rels_s) // 2], 'q10_rel': rels_s[max(0, len(rels_s) // 10)],
                                           'hist': hist(rels), 'mean_iters': float(np.mean([r['iters'] for r in res])),
                                           'n_cands': sum(1 for r in res if 'Dc' in r), 'seconds': round(time.time() - tb, 1),
                                           'distinct_minima': len({round(math.log10(max(r, 1e-300)), 2) for r in rels})})
            fplib.append_json(JOBF, rec); fplib.wlog(rec)
            print({k: rec[k] for k in ('d', 'sector', 'batch_no', 'restarts', 'best_rel', 'median_rel', 'hist', 'n_cands', 'seconds')}, flush=True)
            for r in res:
                if 'Dc' not in r:
                    continue
                Dc = np.array(r['Dc'])
                crec = {'job': args.job, 'd': d, 'sector': list(c), 'seed': r['seed'], 'relres': r['relres'],
                        'D_float': FAM.D_float_str(Dc), 'ts': fplib.now()}
                print('CANDIDATE', crec, flush=True)
                recs = fplib.try_reconstruct_and_confirm(FAM, Dc, G=G)
                crec['recon'] = [{k: v for k, v in rr.items() if k not in ('N', 'Q')} for rr in recs]
                fplib.append_json(CANDF, crec)
                for rr in recs:
                    if rr.get('status') == 'solved' and rr.get('verified'):
                        Dx = core.parse(rr['D']); N = [core.parse(n) for n in rr['N']]; Q = core.parse(rr['Q'])
                        ok = fplib.save_witness(args.instance, Dx, N, Q, f'{LANE}/witness-candidate.json')
                        fplib.wlog({'batch': 'WITNESS', 'agent': 'vm-search', 'job': args.job, 'd': d, 'sector': list(c), 'D': rr['D'], 'verified': ok})
                        print('!!!!! WITNESS FOUND', rr['D'], 'verified', ok, flush=True)
                        sys.exit(0)
print('done', args.job, flush=True)
