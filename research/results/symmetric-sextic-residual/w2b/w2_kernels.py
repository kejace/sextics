"""E1-admissible denominator spaces per character sector and degree (mod p), served-beyond-J at random
points, Fermat-deformation scan, and the unrestricted kernels at d = 9, 10."""
import os, sys, os, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2lib import *
from charsolve import e1_kernel_char, CharSystem, J_char_basis
from hodge import jacobian
ap = argparse.ArgumentParser(); ap.add_argument('--job', default='kernels'); ap.add_argument('--part', default='chars,sub,fermat,full')
args = ap.parse_args()
JOBLOG = f'{RES}/w2b/{args.job}.log.jsonl'; os.makedirs(f'{RES}/w2b', exist_ok=True)
F, A, q = load(); p = P0
parts = args.part.split(',')
if 'chars' in parts:
    for chi in ('triv', 'sgnX', 'sgnY', 'sgnsgn'):
        G = GROUPS[chi]()
        for d in range(4, 15):
            for pD in (0, 1):
                t0 = time.time()
                info, polys = e1_kernel_char(F, A, q, d, pD, G, p)
                rec = {'job': args.job, 'family': 'e1-kernel-dims', 'chi': chi, 'd': d, 'pD': pD, 'nD': info['nD'],
                       'dimK': info['dimK'], 'dimJ': info['dimJ'], 'ker_beyond_J': info['ker'], 'e1_rows': info['rows'],
                       'e1_nN': info['nN'], 'seconds': round(time.time() - t0, 1)}
                # served targets beyond J at up to 2 random admissible points (kernel class + random J element), d <= 10
                if d <= 10 and info['ker'] > 0:
                    rng = np.random.default_rng(d * 10 + pD)
                    Jp = J_char_basis(F, d, pD, G)
                    served = []
                    for trial in range(2):
                        D = {}
                        for P in polys:
                            D = core.add(D, core.scale(P, Fraction(int(rng.integers(-3, 4)))))
                        for P in (Jp[i] for i in rng.choice(len(Jp), size=min(3, len(Jp)), replace=False)) if Jp else []:
                            D = core.add(D, core.scale(P, Fraction(int(rng.integers(-3, 4)))))
                        if not D: continue
                        try:
                            S = CharSystem(F, A, q, D, G, p=p, targets=True)
                            tg = S.served_targets()
                            # beyond J(F)_6: rank([J6 | tg]) - rank(J6), in coordinates of the 19 invariant odd sextic orbits
                            import invsolve
                            Treps = invsolve.orbit_reps(6, 1, invsolve.canon_codes)
                            vec = lambda P_: [P_.get(r, 0) % p for r in Treps]
                            Jv = [vec(core.reduce_mod(P_, p)) for P_ in J_char_basis(F, 6, 1, GROUPS['triv'](), p)]
                            def rank_of(vs):
                                return nmod_mat(len(vs), len(Treps), [v[i] for v in vs for i in range(len(Treps))], p).rank() if vs else 0
                            vs = [vec(P_) for P_ in tg]
                            beyond = rank_of(Jv + vs) - rank_of(Jv)
                            served.append({'served_total': len(vs), 'beyond_J': beyond})
                        except Exception as ex:
                            served.append({'error': repr(ex)[:200]})
                    rec['served_random_points'] = served
                wlog(rec, JOBLOG); print(json.dumps(rec), flush=True)
if 'sub' in parts:
    for chi in ('S22_++', 'S22_-+', 'S22_+-', 'S22_--', 'S3S2_++', 'S3S2_-+', 'S3S2_+-', 'S3S2_--', 'S2S3_++', 'S2S3_-+', 'S2S3_+-', 'S2S3_--'):
        G = GROUPS[chi]()
        for d in (10,):
            for pD in (0, 1):
                t0 = time.time()
                try:
                    info, polys = e1_kernel_char(F, A, q, d, pD, G, p)
                    rec = {'job': args.job, 'family': 'e1-kernel-dims', 'chi': chi, 'd': d, 'pD': pD, 'nD': info['nD'],
                           'dimK': info['dimK'], 'dimJ': info['dimJ'], 'ker_beyond_J': info['ker'], 'e1_rows': info['rows'],
                           'e1_nN': info['nN'], 'seconds': round(time.time() - t0, 1)}
                except Exception as ex:
                    rec = {'job': args.job, 'family': 'e1-kernel-dims', 'chi': chi, 'd': d, 'pD': pD, 'error': repr(ex)[:300]}
                wlog(rec, JOBLOG); print(json.dumps(rec), flush=True)
if 'fermat' in parts:
    # F_t = Fermat + t * e2X e2Y (e2X + 2 e2Y): invariant and sgnX kernels for several t
    QX, QY = e_k(X, 2), e_k(Y, 2)
    M = core.mul(core.mul(QX, QY), core.add(QX, core.scale(QY, 2)))
    Fermat = core.add(p_k(X, 6), p_k(Y, 6))
    for t in (Fraction(0), Fraction(5, 7), Fraction(1), Fraction(2), Fraction(1, 3), Fraction(-1), Fraction(7, 5), Fraction(3, 11), Fraction(-5, 7)):
        Ft = core.add(Fermat, core.scale(M, t))
        for chi in ('triv', 'sgnX'):
            G = GROUPS[chi]()
            for d in range(6, 11):
                for pD in (0, 1):
                    t0 = time.time()
                    info, polys = e1_kernel_char(Ft, A, q, d, pD, G, p)
                    rec = {'job': args.job, 'family': 'fermat-deformation-kernel', 't': str(t), 'chi': chi, 'd': d, 'pD': pD,
                           'dimK': info['dimK'], 'dimJ': info['dimJ'], 'ker_beyond_J': info['ker'], 'seconds': round(time.time() - t0, 1)}
                    wlog(rec, JOBLOG); print(json.dumps(rec), flush=True)
if 'full' in parts:
    for d in (9, 10, 11):
        for pD in (0, 1):
            t0 = time.time()
            try:
                out = jacobian.mult_kernel(F, A, d, p, parity=pD, verbose=True)
                rec = {'job': args.job, 'family': 'e1-kernel-dims-unrestricted', 'chi': 'all', 'd': d, 'pD': pD,
                       'ker_beyond_J': out['dim'], 'kernel_dim': out['kernel_dim'], 'jacobian_dim': out['jacobian_dim'],
                       'ring_dim': out['ring_dim'], 'rows': out['rows'], 'cols': out['cols'], 'seconds': round(time.time() - t0, 1)}
                import pickle
                pickle.dump(out['basis'], open(f'{RES}/w2b/{args.job}_fullkernel_d{d}_p{pD}.pkl', 'wb'))
            except Exception as ex:
                rec = {'job': args.job, 'family': 'e1-kernel-dims-unrestricted', 'd': d, 'pD': pD, 'error': repr(ex)[:300]}
            wlog(rec, JOBLOG); print(json.dumps(rec), flush=True)
print('ALL DONE', flush=True)
