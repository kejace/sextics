"""Wave-2 re-verification sweep (lens excl, lane cycle-q2).

For every candidate denominator D: engine-v2 sector-restricted fixed-D solve at prime P1 (and P2),
and, when inconsistent, an exact rational Farkas certificate (prime-independent proof of
inconsistency over C).  Candidates: dossier A.5/A.6 and wave-1 structural families.

usage: python3 wave2-excl/reverify.py --workers 32 [--only TAGPREFIX] [--no-farkas] [--skip-d12]
"""
import sys, os, json, time, argparse, itertools, traceback
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hodge import core
from hodge.core import mul, add, scale, pderiv, parse, to_str, power

OUT = f'{w2lib.PACKET}/reverify.results.jsonl'
CERTDIR = f'{w2lib.PACKET}/certs'

def candidates(F):
    X = ['X0','X1','X2','Y0','Y1','Y2']
    Fi = [pderiv(F, i) for i in range(6)]
    C = []
    def addc(tag, fam, D, note=''):
        C.append({'tag': tag, 'family': fam, 'D': D, 'note': note})
    # --- dossier A.5 ---
    for i in range(6):
        addc(f'F_{X[i]}', 'deg5-J', Fi[i], 'partial dF/d'+X[i]); 
    for i in range(6):
        addc(f'F_{X[i]}^2', 'deg5-J-sq', mul(Fi[i], Fi[i]), 'square of partial')
    B6 = {}
    for i in range(3):
        B6 = add(B6, mul(parse(X[i]), Fi[i]))
    addc('B6', 'fqlx-support', B6, 'sum_i X_i F_{X_i} (fqlx B6)')
    addc('B6^2', 'fqlx-support-sq', mul(B6, B6))
    for i in range(3):
        Gr = add(mul(parse(X[i]), Fi[i]), mul(parse(X[i+3]), Fi[i+3]))
        addc(f'Grad_{i}', 'block-radial', Gr, '(X_i d/dX_i + Y_i d/dY_i) F'); addc(f'Grad_{i}^2', 'block-radial-sq', mul(Gr, Gr))
    for j in range(6):
        for i in range(6):
            addc(f'{X[j]}*F_{X[i]}', 'ZjFi', mul(parse(X[j]), Fi[i]))
    addc('X0^2X1^2X2^2', 'monomial', parse('X0**2*X1**2*X2**2')); addc('Y0^2Y1^2Y2^2', 'monomial', parse('Y0**2*Y1**2*Y2**2'))
    K10 = add(mul(parse('2*Y2**5'), Fi[4]), scale(mul(parse('2*Y1**5'), Fi[5]), -1))
    addc('K10', 'fqlx-support', K10, '2(Y2^5 d_{Y1} - Y1^5 d_{Y2}) F')
    for i in range(6):
        for j in range(i+1, 6):
            addc(f'F_{X[i]}*F_{X[j]}', 'FiFj', mul(Fi[i], Fi[j]))
    # --- wave-1 structural families (D and D^2) ---
    P = parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
    M = parse('X0**3*X1*Y0*Y1 + X1**3*X2*Y1*Y2 + X0*X2**3*Y0*Y2')
    fam = {'pure-power': P, 'coupling-sum': M, 'X-block': parse('X0**6 + X1**6 + X2**6'), 'Y-block': parse('Y0**6 + Y1**6 + Y2**6'),
           'pair-fermat-0': parse('X0**6 + Y0**6'), 'pair-fermat-1': parse('X1**6 + Y1**6'), 'pair-fermat-2': parse('X2**6 + Y2**6'),
           'edge-01': parse('X0**6 + X1**6 + Y0**6 + Y1**6 + X0**3*X1*Y0*Y1/2'), 'edge-12': parse('X1**6 + X2**6 + Y1**6 + Y2**6 + X1**3*X2*Y1*Y2/2'),
           'edge-20': parse('X0**6 + X2**6 + Y0**6 + Y2**6 + X0*X2**3*Y0*Y2/2'),
           'coupling-01': parse('X0**3*X1*Y0*Y1'), 'coupling-12': parse('X1**3*X2*Y1*Y2'), 'coupling-20': parse('X0*X2**3*Y0*Y2')}
    for k, D in fam.items():
        addc(k, 'wave1-structural', D); addc(k+'^2', 'wave1-structural-sq', mul(D, D))
    for sub4 in itertools.combinations(range(6), 4):
        D = {}
        for i in sub4: D = add(D, parse(X[i]+'**6'))
        nm = 'cone4-' + ''.join(X[i] for i in sub4)
        addc(nm, 'fermat-cone4', D); addc(nm+'^2', 'fermat-cone4-sq', mul(D, D))
    for sub3 in itertools.combinations(range(6), 3):
        D = {}
        for i in sub3: D = add(D, parse(X[i]+'**3'))
        nm = 'cone3-' + ''.join(X[i] for i in sub3)
        addc(nm, 'fermat-cone3', D); addc(nm+'^2', 'fermat-cone3-sq', mul(D, D))
        D4 = {}
        for i in sub3: D4 = add(D4, parse(X[i]+'**4'))
        addc('cone3q-' + ''.join(X[i] for i in sub3), 'fermat-cone3-d4', D4)
    return C

def run_one(args):
    cand, do_farkas, p1, p2 = args
    F, A, q, G = w2lib.load()
    D = cand['D']
    t0 = time.time()
    rec = {'exp': 'reverify', 'tag': cand['tag'], 'family': cand['family'], 'd': core.degree(D), 'note': cand.get('note', ''),
           'D': to_str(D) if len(to_str(D)) < 300 else to_str(D)[:300] + '...'}
    try:
        r1, S, _ = w2lib.solve_modp(F, A, q, D, p1, G)
        rec.update({'prime': p1, 'consistent_p1': r1['consistent'], 'rank_aug_p1': r1['rank'], 'cols': r1['cols'], 'rows': r1['rows'], 'stab': r1['stab']})
        r2, _, _ = w2lib.solve_modp(F, A, q, D, p2, G)
        rec.update({'prime2': p2, 'consistent_p2': r2['consistent'], 'rank_aug_p2': r2['rank']})
        if r1['consistent']:
            rec['result'] = 'CONSISTENT mod p1 -- needs exact lift (potential witness!)'
        else:
            rec['result'] = 'inconsistent mod p1' + ('' if not r2['consistent'] else ' but CONSISTENT mod p2 (bad prime?)')
            rec['rank'] = f"{r1['rank']-1} (M) / {r1['rank']} (aug) of {r1['cols']} cols"
            if do_farkas:
                os.makedirs(CERTDIR, exist_ok=True)
                fr = w2lib.farkas(F, A, q, D, p1, G, save=f"{CERTDIR}/{cand['tag'].replace('*','_').replace('^','pow').replace('/','_')}.json", Dstr=rec['D'])
                rec.update({'farkas': fr['status'], 'farkas_rows': fr.get('farkas_rows'), 'bits_max': fr.get('bits_max'), 'sha256_y': fr.get('sha256_y'), 'yb': fr.get('yb')})
                if fr['status'] == 'certified':
                    rec['result'] += '; exact Farkas certificate verified (rigorous over C)'
    except Exception as e:
        rec['result'] = 'ERROR ' + repr(e); rec['trace'] = traceback.format_exc()[-800:]
    rec['seconds'] = round(time.time() - t0, 1); rec['job'] = os.environ.get('W2JOB', 'reverify')
    w2lib.log(**rec)
    with open(OUT, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=8); ap.add_argument('--only', default=None)
    ap.add_argument('--no-farkas', action='store_true'); ap.add_argument('--skip-d12', action='store_true'); ap.add_argument('--max-d', type=int, default=99)
    a = ap.parse_args()
    F, A, q, G = w2lib.load()
    C = candidates(F)
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try: done.add(json.loads(line)['tag'])
            except Exception: pass
    todo = [c for c in C if c['tag'] not in done and (a.only is None or c['tag'].startswith(a.only)) and core.degree(c['D']) <= a.max_d
            and not (a.skip_d12 and core.degree(c['D']) >= 12)]
    # big systems first for load balance
    todo.sort(key=lambda c: -core.degree(c['D']))
    print(f'{len(C)} candidates, {len(todo)} to do', flush=True)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run_one, (c, not a.no_farkas, w2lib.P1, w2lib.P2)) for c in todo]
        for f in as_completed(futs):
            r = f.result()
            print(r['tag'], r['d'], r['result'], r.get('rank'), r.get('farkas'), r['seconds'], flush=True)
    print('ALL DONE', flush=True)
