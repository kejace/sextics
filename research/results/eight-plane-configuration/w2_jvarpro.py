"""Wave-2 job: full-E2 variable-projection search over a denominator family (multiprocess restarts).

usage: python3 w2_jvarpro.py --family jacobian|kernel --d D --par P --restarts R --seed S --nproc N
                             [--max-nfev 200] [--stage1 0|1] [--stage1-only 0|1] [--npz PATH]
                             [--threads 1] [--tag NAME] [--recon-rel 1e-9]
Writes <tag>.jsonl (one line per restart) and <tag>.best.json in the lane dir; logs reconstruction attempts.
"""
import sys, os, json, time, math, argparse
ap = argparse.ArgumentParser()
ap.add_argument('--family', default='jacobian')
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--par', default='0', help='0, 1 or mixed')
ap.add_argument('--restarts', type=int, default=100)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--max-nfev', type=int, default=200)
ap.add_argument('--stage1', type=int, default=0)
ap.add_argument('--stage1-only', type=int, default=0)
ap.add_argument('--stage1-nfev', type=int, default=100)
ap.add_argument('--npz', default=None)
ap.add_argument('--threads', type=int, default=1)
ap.add_argument('--tag', default=None)
ap.add_argument('--recon-rel', type=float, default=1e-9)
ap.add_argument('--cond', type=float, default=1e-10)
ap.add_argument('--precompute', type=int, default=0, help='only build the family (and cache the R-projector), then exit')
args = ap.parse_args()
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = str(args.threads)
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from w2lib import _col, _dense
import multiprocessing as mp

PAR = None if str(args.par) == 'mixed' else int(args.par)
TAG = args.tag or f'search-{args.family}-d{args.d}-p{args.par}-s{args.seed}'
OUT = f'{WD}/{TAG}.jsonl'
BEST = f'{WD}/{TAG}.best.json'
need_C = bool(args.stage1 or args.stage1_only)

_fam = None


def get_family():
    global _fam
    if _fam is None:
        if args.family == 'jacobian':
            _fam = JacobianFamily(args.d, PAR, verbose=True, want_C=need_C, cond=args.cond)
        else:
            _fam = KernelFamily(args.d, PAR, args.npz, verbose=True, want_C=need_C, cond=args.cond)
    return _fam


def work(seed):
    fam = get_family()
    t0 = time.time()
    try:
        if args.stage1_only:
            rng = np.random.default_rng(seed)
            u0 = rng.standard_normal(fam.dimP); u0 /= np.linalg.norm(u0)
            res = fam.lm(u0, max_nfev=args.stage1_nfev, stage1=True)
            a, _ = fam.a_of_u(res.x)
            D = _col(fam.Dfull @ a); Om = fam.pat.omega(D); y0 = _col(fam.Nfull @ a); b = Om @ y0
            rp = float(np.linalg.norm(fam.C.T @ b)); bn = float(np.linalg.norm(b))
            info = {'seed': int(seed), 'stage1_only': True, 'r_part': rp, 'bnorm': bn, 'r_part_rel': rp / max(bn, 1e-300),
                    'nfev': int(res.nfev), 'status': int(res.status), 'seconds': round(time.time() - t0, 2)}
            rel = info['r_part_rel']
        else:
            info, a, o = fam.restart(seed, max_nfev=args.max_nfev, stage1_first=bool(args.stage1),
                                     stage1_nfev=args.stage1_nfev)
            rel = info['rel']
        info['family'] = fam.name; info['d'] = fam.d; info['par'] = fam.pD; info['tag'] = TAG
        info['a_top'] = [float(x) for x in np.sort(np.abs(a))[::-1][:5]]
        if rel < 1e-6 or (args.stage1_only and fam.K <= 400):
            info['a'] = [float(x) for x in a]
        if rel < args.recon_rel and not args.stage1_only:
            info['reconstruction'] = fam.confirm_rational(a, log_path=OUT, note=f'{TAG} seed {seed}')
        return info
    except Exception as ex:
        import traceback
        return {'seed': int(seed), 'error': f'{ex.__class__.__name__}: {ex}', 'trace': traceback.format_exc()[-1500:],
                'tag': TAG, 'seconds': round(time.time() - t0, 2)}


if __name__ == '__main__':
    t0 = time.time()
    fam = get_family()
    if args.precompute:
        print(f'[precompute] {fam.name} built; C cached={fam.C is not None} ({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)
        sys.exit(0)
    jlog(OUT, {'event': 'start', 'tag': TAG, 'family': fam.name, 'd': fam.d, 'par': fam.pD, 'K': fam.K, 'dimZ': fam.dimZ,
               'dimP': fam.dimP, 'rows': fam.rows, 'nN': fam.nN, 'nc': fam.nc, 'nQ': fam.nQ, 'nnz_omega': fam.pat.nnz,
               'dimC': (fam.C.shape[1] if fam.C is not None else None), 'args': vars(args),
               'ts': time.strftime('%Y-%m-%dT%H:%M:%S')})
    seeds = [args.seed * 100000 + k for k in range(args.restarts)]
    best = None
    done = 0
    ctx = mp.get_context('fork')
    with ctx.Pool(args.nproc) as pool:
        for info in pool.imap_unordered(work, seeds):
            done += 1
            jlog(OUT, info)
            key = info.get('r_part_rel') if args.stage1_only else info.get('rel')
            if key is not None and (best is None or key < best.get('_key', 1e300)):
                best = dict(info); best['_key'] = key
                with open(BEST, 'w') as f:
                    json.dump(best, f, indent=1)
            if done % 10 == 0 or done == len(seeds):
                print(f'[{TAG}] {done}/{len(seeds)} restarts, best key {best.get("_key") if best else None} '
                      f'({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    jlog(OUT, {'event': 'end', 'tag': TAG, 'restarts': len(seeds), 'best': best, 'seconds': round(time.time() - t0, 1),
               'ts': time.strftime('%Y-%m-%dT%H:%M:%S')})
