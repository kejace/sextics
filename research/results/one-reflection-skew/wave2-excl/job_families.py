"""Structured-family exact sweep (goal c): products of named diagonal-eigenvector invariants, degrees dmin..dmax,
at most `maxf` factors.  Stage 1: E1 test mod p in batch per (degree, sector), rigorous because rank_p(J_{d+6}[chi])
is certified equal to the exact rank (exact J-duals).  Stage 2: survivors get the sector-restricted E1+E2 fixed-D
test with an exact dual certificate (EXCLUDED-exact) or an exact witness.
Usage: python3 job_families.py <dmin> <dmax> <nproc> <maxf>"""
import sys, os, json, time, itertools, math, multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x
from w2x import ctx, Z
from flint import nmod_mat
P = 2147483647
X0, X1, X2, Y0, Y1, Y2 = Z
LET = {
 'X0': X0, 'X1': X1, 'X2': X2, 'Y0': Y0, 'Y1': Y1, 'Y2': Y2,
 'cq': X0**3 + X1**2 * X2, 'c1': X0**3 + 7 * X1**2 * X2,
 't': Y0**2 - Y1**2, 'm': Y0**2 + Y1**2, 'u': X0**2 + X2**2, 'v': X0**2 - X2**2,
 'G0': X0**6 + X1**6 + X2**6, 'Y6': Y0**6 + Y1**6, 'f': w2x.F7 - 7 * Y2**6,
 'h1': X0**3 * X1 * Y0 * Y1 + Y0**6, 'h2': X1**3 * X2 * Y0 * Y1 + Y1**6,
 'h3': X0 * X1 * Y0**3 * Y1 - X0 * X1 * Y0 * Y1**3 + 7 * X2**6,
 'g1': X0**6 + Y0**6, 'g2': X1**6 + Y1**6, 'g3': X2**6 + Y2**6, 'g4': X0**6 + Y2**6,
 'q2': X0 * Y1**2 + X2 * Y0**2, 'r1': X0**2 * Y0 * Y1 + X1 * X2 * Y1**2 * 0 + X1**4, 'r2': X0 * X1 * Y0 * Y1 + X2**4,
 'w1': X0**3 * Y0**3 + X1**6, 'w2': X0**3 * Y0 * Y1**2 + X2**6,
}
LET = {k: v for k, v in LET.items() if w2x.poly_chr(v) is not None}
DEG = {k: w2x.poly_deg(v) for k, v in LET.items()}

def candidates(dmin, dmax, maxf):
    names = sorted(LET, key=lambda k: (DEG[k], k))
    seen = {}
    for r in range(1, maxf + 1):
        for combo in itertools.combinations_with_replacement(names, r):
            dd = sum(DEG[k] for k in combo)
            if dd < dmin or dd > dmax: continue
            if combo.count('f') > 1: continue
            Pp = ctx.from_dict({(0,)*6: 1})
            for k in combo: Pp = Pp * LET[k]
            key = tuple(sorted(Pp.to_dict().items()))
            if key in seen: continue
            seen[key] = '*'.join(combo)
    return [(nm, ctx.from_dict(dict(key))) for key, nm in seen.items()]

def e1_batch(args):
    """All candidates of one (d, chi): E1 mod p with certified J rank.  Returns list of (name, Ddict, e1pass)."""
    d, chi, group = args
    n = d + 6; tchi = w2x.chr_add(chi, w2x.CHI_A)
    rows = w2x.secmons(n, tchi)
    Y, rankJ_exact = w2x.jac_dual(n, tchi)
    gens = w2x.jgens(n, tchi)
    cols = [w2x.vec_of(g[2], n, tchi) for g in gens]
    Jp = w2x.nmod_from_cols(cols, len(rows), P)
    rankJ_p = Jp.rank()
    assert rankJ_p == rankJ_exact, ('J rank drop mod p', d, chi, rankJ_p, rankJ_exact)
    # left kernel mod p of J, applied to A*D
    Lk, k = Jp.transpose().nullspace()           # len(rows) x k
    out = []
    for nm, Dd in group:
        D = ctx.from_dict(Dd)
        v = w2x.vec_of(w2x.A * D, n, tchi)
        ok = True
        for j in range(k):
            s = 0
            for r, c in v.items():
                s = (s + int(Lk[r, j]) * c) % P
            if s: ok = False; break
        out.append((nm, Dd, ok, d, list(chi)))
    return out

def e2_work(args):
    nm, Dd = args
    D = ctx.from_dict(Dd); d = w2x.poly_deg(D)
    out = w2x.test_fixed_D(D, nm, 'structured-products', p=P, want_exact_dual=True)
    rec = {'name': nm, 'd': d, 'chi': list(out['chi']), 'e1': True, 'status': out['status'], 'rank': out['rank'], 'rank_aug': out['rank_aug'],
           'nN': out['nN'], 'nQ': out['nQ'], 'rows': out['rows'], 'seconds': out['seconds']}
    if 'witness' in out: rec['witness'] = out['witness']
    return rec

if __name__ == '__main__':
    dmin, dmax = int(sys.argv[1]), int(sys.argv[2]); nproc = int(sys.argv[3]); maxf = int(sys.argv[4])
    t0 = time.time()
    cands = candidates(dmin, dmax, maxf)
    print('letters', {k: DEG[k] for k in LET}, 'candidates', len(cands), flush=True)
    groups = {}
    for nm, D in cands:
        groups.setdefault((w2x.poly_deg(D), w2x.poly_chr(D)), []).append((nm, D.to_dict()))
    tasks = [(d, chi, g) for (d, chi), g in sorted(groups.items())]
    outp = f'families_d{dmin}-{dmax}_f{maxf}.jsonl'
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            try: done.add(json.loads(line)['name'])
            except Exception: pass
    e1pass = []
    nfail = 0
    with mp.Pool(nproc) as pool, open(outp, 'a') as fh:
        for res in pool.imap_unordered(e1_batch, tasks):
            for nm, Dd, ok, d, chi in res:
                if nm in done: continue
                if ok: e1pass.append((nm, Dd))
                else:
                    nfail += 1
                    fh.write(json.dumps({'name': nm, 'd': d, 'chi': chi, 'e1': False, 'status': 'E1-fails-rigorous'}) + '\n')
        fh.flush()
        print(f'E1 stage: {nfail} fail (rigorous: J rank certified), {len(e1pass)} pass -> E2', round(time.time()-t0, 1), flush=True)
        w2x.logexp(exp='family-E1-batch', family='structured-products', d=f'{dmin}-{dmax}', prime=P, result=f'{nfail} E1-fail (rigorous), {len(e1pass)} E1-pass',
                   rank=None, seconds=round(time.time()-t0, 1), job='excl-ors-families', note=f'{len(cands)} candidates, <= {maxf} factors; J_{{d+6}}[chi] rank mod p certified equal to exact rank in every sector used')
        e1pass.sort(key=lambda x: sum(next(iter(x[1]))))
        for rec in pool.imap_unordered(e2_work, e1pass):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print(rec['name'], rec['d'], rec['status'], rec.get('seconds'), flush=True)
            if rec['status'] == 'WITNESS':
                print('!!!!!!!! WITNESS FOUND !!!!!!!!', rec, flush=True)
                with open('WITNESS_FOUND.json', 'w') as wf: json.dump(rec, wf)
    print('DONE', round(time.time()-t0, 1), flush=True)
