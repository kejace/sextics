"""Single-D solves (mod p) of the E1-kernel classes of the subgroup types S2xS2 / S3xS2 / S2xS3 (all four sign
characters) at d = 8, 9 (+ random combos with J), 8 workers.  These probe non-S3xS3-symmetric denominators
(orbit products would be invariant of degree 9d / 3d)."""
import sys, os, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from w2lib import *
from charsolve import e1_kernel_char, J_char_basis
import multiprocessing as mp
ap = argparse.ArgumentParser(); ap.add_argument('--job', default='singlesub'); ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--types', default='S22_++,S22_-+,S22_+-,S22_--,S3S2_++,S3S2_-+,S3S2_+-,S3S2_--,S2S3_++,S2S3_-+,S2S3_+-,S2S3_--')
ap.add_argument('--dmin', type=int, default=8); ap.add_argument('--dmax', type=int, default=9); ap.add_argument('--combos', type=int, default=2)
args = ap.parse_args()
JOBLOG = f'{RES}/w2/{args.job}.log.jsonl'; os.makedirs(f'{RES}/w2', exist_ok=True)
F, A, q = load(); p = P0

def run_task(task):
    kind, label, chi, D, meta = task
    G = GROUPS[chi](); t0 = time.time()
    try:
        r = solve_single(F, A, q, D, G, p, lift=True)
        rec = dict(meta, job=args.job, family=kind, label=label, chi=chi, d=core.degree(D), pD=core.parity_of(D),
                   e1=r['e1'], consistent=r['consistent'], rank=r['rank'], nullity=r['nullity'], rows=r['rows'], cols=r['cols'],
                   seconds=round(time.time() - t0, 1))
        if r['consistent']:
            rec['lift_status'] = r.get('lift_status')
            if 'N' in r:
                ok, info = core.verify_witness(F, A, q, D, r['N'], r['Q']); rec['verified'] = bool(ok)
                if ok:
                    w = witness_json(D, r['N'], r['Q']); w.update({'instance': LANE, 'family': f'{kind}:{label}', 'engine': ENGINE})
                    path = f'{RES}/witness-candidate.json'; json.dump(w, open(path, 'w'), indent=1); rec['WITNESS'] = path
                    print('WITNESS FOUND', path, flush=True)
        return rec
    except Exception as ex:
        return dict(meta, job=args.job, family=kind, label=label, chi=chi, error=repr(ex)[:300], seconds=round(time.time() - t0, 1))

tasks = []
for chi in args.types.split(','):
    G = GROUPS[chi]()
    for d in range(args.dmin, args.dmax + 1):
        for pD in (0, 1):
            t0 = time.time()
            info, polys = e1_kernel_char(F, A, q, d, pD, G, p)
            rec = {'job': args.job, 'family': 'e1-kernel-dims', 'chi': chi, 'd': d, 'pD': pD, 'nD': info['nD'], 'dimK': info['dimK'],
                   'dimJ': info['dimJ'], 'ker_beyond_J': info['ker'], 'seconds': round(time.time() - t0, 1)}
            wlog(rec, JOBLOG); print(json.dumps(rec), flush=True)
            if info['ker'] == 0: continue
            Jp = J_char_basis(F, d, pD, G)
            for i, P in enumerate(polys):
                tasks.append(('kernel-class', f'{chi}/d{d}/p{pD}/class{i}', chi, P, {'d': d, 'pD': pD, 'class': i, 'n_classes': info['ker']}))
            rng = np.random.default_rng(100 * d + pD + 7)
            for trial in range(args.combos):
                D = {}
                for P in polys: D = core.add(D, core.scale(P, Fraction(int(rng.integers(-3, 4)))))
                if Jp:
                    for j in rng.choice(len(Jp), size=min(4, len(Jp)), replace=False): D = core.add(D, core.scale(Jp[j], Fraction(int(rng.integers(-3, 4)))))
                if D: tasks.append(('kernel-random-combo', f'{chi}/d{d}/p{pD}/combo{trial}', chi, D, {'d': d, 'pD': pD, 'trial': trial}))
print(f'{len(tasks)} tasks', flush=True)
tasks.sort(key=lambda t: -core.degree(t[3]))
ctx = mp.get_context('fork')
with ctx.Pool(args.workers) as pool:
    for rec in pool.imap_unordered(run_task, tasks):
        wlog(rec, JOBLOG)
        print(json.dumps({k: v for k, v in rec.items() if k in ('family', 'label', 'chi', 'd', 'pD', 'e1', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds', 'error', 'WITNESS')}), flush=True)
print('ALL DONE', flush=True)
