"""Driver: multi-restart fastvp search on diagonal-permutation.
usage: python3 run_fastvp.py --d 7 --family jac|jac+ker --parity 0 --sector all|inv --restarts 100 --procs 8 --seed 0 --maxiter 400 [--tag T]
Kernel elements (family jac+ker) are read from w2search/kernel-d<d>-par<par>-exact.json (list of D strings) if present,
else from ../D8_e1kernel.json for d=8.
Appends per-restart lines to w2search/fastvp-<tag>.jsonl and one batch summary line to ../wave2-search.log.jsonl.
Candidates with rel < RECON_TOL are rationally reconstructed and re-solved exactly (fixed_d); witnesses -> ../witness-candidate.json.
"""
import sys, os, time, json, argparse, math
sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/w2search'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from hodge import core, jacobian, fixed_d
import fastvp
LANE = '/tmp/hodge/results/diagonal-permutation'
RECON_TOL = 1e-10
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--family', default='jac'); ap.add_argument('--parity', type=int, default=0)
ap.add_argument('--sector', default='all'); ap.add_argument('--restarts', type=int, default=20); ap.add_argument('--procs', type=int, default=1)
ap.add_argument('--seed', type=int, default=0); ap.add_argument('--maxiter', type=int, default=400); ap.add_argument('--tag', default=None)
ap.add_argument('--ptseed', type=int, default=11); ap.add_argument('--relax', action='store_true', help='K=0 (no Koszul unknowns) relaxation')
ap.add_argument('--oversample', type=float, default=1.12)
args = ap.parse_args()
tag = args.tag or f'd{args.d}-{args.family}-p{args.parity}-{args.sector}{"-relax" if args.relax else ""}-s{args.seed}'
OUT = f'{LANE}/w2search/fastvp-{tag}.jsonl'
def log_batch(e):
    e.setdefault('engine', 2); e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); e.setdefault('lens', 'search')
    with open(f'{LANE}/wave2-search.log.jsonl', 'a') as f: f.write(json.dumps(e, default=str) + '\n')
def log_line(e):
    with open(OUT, 'a') as f: f.write(json.dumps(e, default=str) + '\n')

F, A, q, rec = core.load_instance('diagonal-permutation')
t0 = time.time()
labels, Bs, N0s = fastvp.jacobian_family(F, A, args.d, args.parity)
nj = len(Bs)
nk = 0
if 'ker' in args.family:
    kf = f'{LANE}/w2search/kernel-d{args.d}-par{args.parity}-exact.json'
    if os.path.exists(kf):
        Ks = [core.parse(s) for s in json.load(open(kf))['D']]
    elif args.d == 8 and args.parity == 0:
        Ks = [core.parse(json.load(open(f'{LANE}/D8_e1kernel.json'))['D'])]
    else:
        Ks = []
    if Ks:
        KN0, res = fastvp.e1_particular(F, A, Ks, verbose=True)
        print(f'[run] kernel part: {len(Ks)} elements, e1 particular residual {res:.2e}', file=sys.stderr)
        for i, (K_, N0_) in enumerate(zip(Ks, KN0)):
            labels.append(('ker', i)); Bs.append(K_); N0s.append(N0_)
        nk = len(Ks)
if args.relax:
    fam = None
fam = fastvp.FastFamily(F, A, Bs, N0s, labels=labels, seed=args.ptseed, sector=args.sector, verbose=True, oversample=args.oversample)
if args.relax:
    # K = 0 relaxation: override lstsq to return the plain Phi residual
    def _lsq0(Psi, Phi, _self=fam):
        if _self.sector == 'inv': Phi = Phi[0::2]
        return np.zeros(Psi.shape[1]), float(np.vdot(Phi, Phi).real)
    fam._lsq = _lsq0
sector_check = None
if args.sector == 'inv':
    full, out = fam.check_sector()
    sector_check = {'full': full, **{k: v for k, v in out.items()}}
    print(f'[run] sector check (lossless Koszul eigenspace): full={full:.6e} minus={out["-"]} plus={out["+"]}', file=sys.stderr)
    # pick the lossless one
    best = min(((abs(v - full) if v is not None else 1e300), k) for k, v in out.items())[1]
    fam.Kcols_basis = fam._sector_bases[best]
    print(f'[run] using Koszul eigenspace {best} dim {fam.Kcols_basis.shape[1]} of {fam.nK}', file=sys.stderr)
fd = fam.fd_check()
print(f'[run] fd check: f={fd[0]:.4e} g.v={fd[1]:.6e} fd={fd[2]:.6e}', file=sys.stderr)
info = {'d': args.d, 'parity': args.parity, 'family': args.family, 'sector': args.sector, 'relax': args.relax, 'K': fam.K, 'nJ': nj, 'nKer': nk,
        'dimP': fam.dimP, 'excluded': fam.excluded, 'P': fam.P, 'Rprime': fam.Rprime, 'nK': fam.nK, 'nQ_elim': fam.nQ,
        'FresMax': fam.Fres, 'build_s': round(fam.build_time, 1), 'fd_check': fd, 'sector_check': sector_check}
log_line({'exp': 'fastvp-family', 'tag': tag, **info})
print('[run] family', json.dumps(info, default=str)[:500], file=sys.stderr)

def worker(w):
    os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
    seed = args.seed * 1000 + w
    n = args.restarts // args.procs + (1 if w < args.restarts % args.procs else 0)
    def lg(e):
        e['worker'] = w; log_line(e)
    cands = fam.search(n_restarts=n, seed=seed, maxiter=args.maxiter, verbose=False, log=lg, tag=tag)
    return [{k: v for k, v in c.items()} for c in cands[:5]]

t1 = time.time()
if args.procs > 1:
    import multiprocessing as mp
    ctx = mp.get_context('fork')
    with ctx.Pool(args.procs) as pool:
        res = pool.map(worker, range(args.procs))
    cands = [c for r in res for c in r]
else:
    cands = worker(0)
cands.sort(key=lambda c: c['rel'])
rels = np.array([c['rel'] for c in cands])
summary = {'exp': 'fastvp-batch', 'tag': tag, **{k: info[k] for k in ('d', 'parity', 'family', 'sector', 'relax', 'K', 'dimP', 'P', 'nK')},
           'restarts': args.restarts, 'best_rel': float(rels.min()), 'median_rel_top': float(np.median(rels)),
           'n_below_1e-6': int((rels < 1e-6).sum()), 'n_below_1e-10': int((rels < 1e-10).sum()), 'seconds': round(time.time() - t1, 1),
           'best_a_top5': [np.round(c['a'], 5).tolist() for c in cands[:1]]}
# rational reconstruction + exact re-solve for the good ones
confirmed = []
for c in cands:
    if c['rel'] > RECON_TOL: break
    for fr in fam.reconstruct(np.array(c['a'])):
        st, data = fam.confirm(fr, verbose=True)
        log_line({'exp': 'fastvp-confirm', 'tag': tag, 'rel': c['rel'], 'frac': [str(x) for x in fr], 'status': st,
                  'data': {k: v for k, v in data.items() if k in ('rank', 'nullity', 'note', 'verify')}})
        if st == 'solved':
            D = data['D']; N = data['N']; Q = data['Q']
            w = {'d': args.d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q), 'tag': tag}
            json.dump(w, open(f'{LANE}/witness-candidate.json', 'w'), indent=1)
            confirmed.append(w); summary['WITNESS'] = True
            print('WITNESS FOUND', file=sys.stderr)
            break
    if confirmed: break
summary['confirmed'] = len(confirmed)
log_batch(summary)
print('[run] SUMMARY', json.dumps(summary, default=str)[:800], file=sys.stderr)
