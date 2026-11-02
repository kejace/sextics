"""Goal (c): individual exact tests of every non-Jacobian E1-kernel basis direction.
For each d in argv: (1) Y2-free family (J_f:A)_d (exact basis + exact lifts, saved to d<d>_5v_fam.pkl) and the strict
5-variable fixed-D test (exact dual) of every 'extra' (non-J_f) direction plus a few random in-sector combinations;
(2) the 6-variable gauge-fixed kernel bases from e1basis_d<d>.pkl: full sector-restricted E1+E2 test (exact dual)
of every direction not in J.  Usage: python3 job_extra_dirs.py <nproc> <d1> <d2> ..."""
import sys, os, json, time, pickle, random, multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x, w2x5
from w2x import ctx
random.seed(7)
P2 = 2147483629

def t5(args):
    nm, Dd, fam = args
    D = ctx.from_dict(Dd)
    out = w2x5.test_fixed_D5(D, nm, fam, want_exact_dual=True)
    return {k: v for k, v in out.items() if k != 'witness'} | ({'witness': out['witness']} if 'witness' in out else {})

def t6(args):
    nm, Dd, fam = args
    D = ctx.from_dict(Dd)
    out = w2x.test_fixed_D(D, nm, fam, want_exact_dual=True)
    return {k: v for k, v in out.items() if k != 'witness'} | ({'witness': out['witness']} if 'witness' in out else {})

if __name__ == '__main__':
    nproc = int(sys.argv[1]); ds = [int(a) for a in sys.argv[2:]]
    for d in ds:
        t0 = time.time()
        famf = f'd{d}_5v_fam.pkl'
        if os.path.exists(famf):
            obj = pickle.load(open(famf, 'rb'))
            items = [(nm, ctx.from_dict(D), [ctx.from_dict(n) for n in N], den, tuple(chi), kind) for (nm, D, N, den, chi, kind) in obj['items']]
            fam5 = w2x5.Family5(d, items)
        else:
            fam5 = w2x5.Family5.from_e1(d, verbose=True)
            pickle.dump({'d': d, 'items': [(nm, D.to_dict(), [n.to_dict() for n in N], den, chi, kind) for (nm, D, N, den, chi, kind) in fam5.items]}, open(famf, 'wb'))
        extra5 = [(it[0], it[1].to_dict(), f'5v-extra-dirs-d{d}') for it in fam5.items if it[5] == 'extra']
        # random in-sector combinations of extra directions (3 per sector)
        bysec = {}
        for it in fam5.items:
            if it[5] == 'extra': bysec.setdefault(it[4], []).append(it[1])
        combos = []
        for chi, Ds in bysec.items():
            if len(Ds) < 2: continue
            for k in range(3):
                Pp = ctx.from_dict({})
                for D in Ds: Pp += D * random.randint(-3, 3)
                if not Pp.is_zero(): combos.append((f'5v:d{d}:{chi}:randcombo{k}', Pp.to_dict(), f'5v-extra-combos-d{d}'))
        print(f'd={d}: 5v family {fam5.n}, extra {len(extra5)}, combos {len(combos)}', round(time.time()-t0, 1), flush=True)
        w2x.logexp(exp=f'5v-d{d}-family', d=d, prime=None, family='(J_f:A)_d exact basis with exact lifts (Y2-free)', result=f'n={fam5.n} extra={len(extra5)}', rank=fam5.n, seconds=round(time.time()-t0, 1), job='excl-ors-extradirs')
        with mp.Pool(nproc) as pool, open(f'extra_dirs_d{d}.jsonl', 'a') as fh:
            for rec in pool.imap_unordered(t5, extra5 + combos):
                fh.write(json.dumps(rec, default=w2x._jsdefault) + '\n'); fh.flush()
                print(rec['name'], rec['status'], rec.get('rank'), rec.get('rank_aug'), rec.get('seconds'), flush=True)
                if rec['status'] == 'WITNESS':
                    json.dump(rec, open('WITNESS_FOUND_5v.json', 'w'), default=w2x._jsdefault); print('!!!!!!!! WITNESS', rec, flush=True)
        # 6-variable extra directions
        bf = f'e1basis_d{d}.pkl'
        if os.path.exists(bf):
            bases = pickle.load(open(bf, 'rb'))
            tasks6 = []
            for chi, blist in bases.items():
                for k, Dd in enumerate(blist):
                    D = ctx.from_dict(Dd)
                    if not w2x.in_jacobian_exact(D):
                        tasks6.append((f'6v:d{d}:{chi}:{k}', Dd, f'6v-extra-dirs-d{d}'))
            print(f'd={d}: 6v non-Jacobian basis directions {len(tasks6)}', flush=True)
            with mp.Pool(nproc) as pool, open(f'extra_dirs6_d{d}.jsonl', 'a') as fh:
                for rec in pool.imap_unordered(t6, tasks6):
                    fh.write(json.dumps(rec, default=w2x._jsdefault) + '\n'); fh.flush()
                    print(rec['name'], rec['status'], rec.get('rank'), rec.get('rank_aug'), rec.get('seconds'), flush=True)
                    if rec['status'] == 'WITNESS':
                        json.dump(rec, open('WITNESS_FOUND_6v.json', 'w'), default=w2x._jsdefault); print('!!!!!!!! WITNESS', rec, flush=True)
        print('d done', d, round(time.time()-t0, 1), flush=True)
    print('DONE', flush=True)
