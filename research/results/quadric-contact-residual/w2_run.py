"""Wave-2 drivers.  Usage:
  w2_run.py dims DMIN DMAX                                  -> E1 colon dims per (d, chi), numeric + mod p
  w2_run.py hess D SECTORS N WORKERS SEED [MAXITER]         -> Hessian relaxation over Jacobian-type family
  w2_run.py koszul D SECTORS N WORKERS SEED [MAXITER]       -> full E2 (Koszul correction) over Jacobian-type family
  w2_run.py varpro D SECTORS N WORKERS SEED [MAXITER] [jac|full] -> full-N sector varpro over colon family
SECTORS: all | each | even | odd | alleven | allodd | 0000,0101,...
"""
import sys, os, json, time
import numpy as np
import w2lib as W
from w2lib import *

def main():
    cmd = sys.argv[1]
    if cmd == 'dims':
        dmin, dmax = int(sys.argv[2]), int(sys.argv[3])
        out_path = RES + '/w2_dims.json'
        allres = json.load(open(out_path)) if os.path.exists(out_path) else {}
        for d in range(dmin, dmax + 1):
            t0 = time.time()
            per = {}
            tot_extra = 0; tot_colon = 0; tot_J = 0; mism = []
            for c in chars_of_degree(d):
                key = f'{d}:{cstr(c)}'
                if key in allres:
                    per[cstr(c)] = allres[key]
                else:
                    mp_ = colon_dims_modp(d, c)
                    cf = colon_family(d, c, want_extras=(mp_['extra'] > 0))
                    rec = dict(mp_); rec['num_colon'] = cf['dim_colon']; rec['num_extra'] = cf['n_extra']
                    rec['mismatch'] = (cf['dim_colon'] != mp_['colon']) if mp_['extra'] > 0 else False
                    allres[key] = rec; per[cstr(c)] = rec
                    json.dump(allres, open(out_path, 'w'))
                r = per[cstr(c)]
                tot_extra += r['extra']; tot_colon += r['colon']; tot_J += r['dimJ']
                if r.get('mismatch'): mism.append(cstr(c))
                print(f'd={d} c={cstr(c)} {r}', file=sys.stderr, flush=True)
            wlog(kind='dims', d=d, n_sectors=len(per), dimJ_total=tot_J, colon_total=tot_colon, extra_total=tot_extra,
                 extra_sectors={k: v['extra'] for k, v in per.items() if v['extra']}, numeric_mismatch=mism,
                 seconds=round(time.time() - t0, 1))
        return
    d = int(sys.argv[2]); spec = sys.argv[3]; n = int(sys.argv[4]); workers = int(sys.argv[5]); seed = int(sys.argv[6])
    maxiter = int(sys.argv[7]) if len(sys.argv) > 7 else 300
    mode = sys.argv[8] if len(sys.argv) > 8 else 'full'
    for CD in sector_list(spec, d):
        tag = spec if spec in ('all', 'alleven', 'allodd') else ','.join(cstr(c) for c in CD)
        t0 = time.time()
        try:
            if cmd in ('hess', 'koszul') and not jac_C_basis(d, CD):
                wlog(kind='batch', family=cmd, d=d, sector=tag, note='empty Jacobian-type family (K=0)'); continue
            if cmd == 'hess':
                fam = HessRelax(d, CD, verbose=True)
                family = 'hess-relax'
                info = {'K': fam.K, 'dimR': fam.dimR, 'gauge': fam.excluded_dim, 'dimP': fam.dimP}
            elif cmd == 'koszul':
                fam = KoszulVarPro(d, CD, verbose=True)
                family = 'koszul-E2'
                info = {'K': fam.K, 'rows': fam.rows, 'cols': fam.ncols, 'nz': fam.nA, 'nQ': fam.nC, 'gauge': fam.excluded_dim, 'dimP': fam.dimP}
            elif cmd == 'varpro':
                basis = []; dims = []
                for c in CD:
                    cf = colon_family(d, c, verbose=True, want_extras=(mode != 'jac'))
                    basis += cf['jac'] + (cf['extras'] if mode != 'jac' else [])
                    dims.append({'c': cstr(c), 'dimS': cf['dimS'], 'dimJ': cf['dimJ'], 'colon': cf['dim_colon'], 'extra': cf['n_extra']})
                if not basis:
                    wlog(kind='batch', family='sector-varpro', d=d, sector=tag, note='empty family'); continue
                fam = SectorVarPro(basis, verbose=True, name=f'varpro d={d} {tag}')
                family = 'sector-varpro-' + mode
                info = {'K': fam.K, 'rows': fam.rows, 'cols': fam.ncols, 'nN': fam.nA, 'nQ': fam.nC, 'gauge': fam.excluded_dim, 'dimP': fam.dimP, 'dims': dims}
            else:
                raise SystemExit('unknown command')
        except Exception as ex:
            wlog(kind='error', family=cmd, d=d, sector=tag, error=repr(ex)[:500]); print('ERROR', ex, file=sys.stderr); continue
        if fam.dimP == 0:
            wlog(kind='batch', family=family, d=d, sector=tag, note='family is pure gauge', **info); continue
        res = search_parallel(fam, n, workers, seed_base=seed, maxiter=maxiter)
        rec = {k: v for k, v in res.items() if k not in ('low', 'best_a')}
        rec['n_low'] = len(res['low'])
        wlog(kind='batch', family=family, d=d, sector=tag, restarts=n, seed_base=seed, maxiter=maxiter, build_s=round(fam.build_time, 1),
             **info, **rec)
        print(f'[batch] {family} d={d} {tag}: best {res["best_rel"]:.3e} hist {res["hist"]} ({res["seconds"]}s)', file=sys.stderr, flush=True)
        # best candidate record (float D) even if not low
        Db = fam.D_float(np.asarray(res['best_a']))
        with open(CAND, 'a') as f:
            f.write(json.dumps({'family': family, 'd': d, 'sector': tag, 'rel': res['best_rel'], 'a': res['best_a'], 'best': True,
                                'Dfloat_top': sorted(((str(e), v) for e, v in Db.items()), key=lambda t: -abs(t[1]))[:40]}) + '\n')
        if res['low']:
            _, found = process_candidates(res['low'], fam, family, d, tag)
            if found:
                print('WITNESS FOUND -- stopping', file=sys.stderr); return

if __name__ == '__main__':
    main()
