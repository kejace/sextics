"""Faster bilinear zero-pair test (see bilinear_pairs.py for the statement).  Sparse matrix fill, random row subsampling with
exact verification (full column rank of a subset certifies full column rank), skips pairs already logged.
Usage: bilinear_pairs2.py quads.json a b [pairs like 000-011,010-101 | all]"""
import sys, json, itertools, time, os, random
from flint import nmod_mat
fn = sys.argv[1]; A = int(sys.argv[2]); B = int(sys.argv[3]); want = sys.argv[4] if len(sys.argv) > 4 else 'all'
LOG = os.environ.get('HODGE_W2LOG', '/tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl')
done = set()
if os.path.exists(LOG):
    for l in open(LOG):
        if '"bilinear-pair"' in l:
            d = json.loads(l)
            if d.get('bidegree') == [A, B] and d['result'].startswith('NO ZERO PAIR'): done.add(tuple(d['pair']))
J = json.load(open(fn)); p = J['p']; chi = [tuple(c) for c in J['chi']]
quads = [{(a, b): v for a, b, v in cf} for cf in J['quads']]
n = len(chi)
cls = lambda c: (c[0], c[1], c[3])
classes = sorted({cls(c) for c in chi})
def monos(vars_, k): return list(itertools.combinations_with_replacement(vars_, k))
def hchar(mono): return sum(chi[a][2] for a in mono) % 6
random.seed(1)
for s, t in itertools.combinations(classes, 2):
    key = (''.join(map(str, s)), ''.join(map(str, t)))
    if key in done or (want != 'all' and '-'.join(key) not in want.split(',')): continue
    t0 = time.time()
    Vs = [a for a in range(n) if cls(chi[a]) == s]; Vt = [a for a in range(n) if cls(chi[a]) == t]
    Ss, St = set(Vs), set(Vt)
    bil = []
    for cf in quads:
        b = {}
        for (a, c), v in cf.items():
            if a in Ss and c in St: b[(a, c)] = v
            elif a in St and c in Ss: b[(c, a)] = v
        if b: bil.append(b)
    cross = [(a, c) for a in Vs for c in Vt]
    M1 = nmod_mat(len(bil), len(cross), p)
    cidx = {m: i for i, m in enumerate(cross)}
    for r, b in enumerate(bil):
        for m, v in b.items(): M1[r, cidx[m]] = v
    rk1 = M1.rank()
    cols = {}
    for ms in monos(Vs, A):
        for mt in monos(Vt, B):
            k2 = ms + mt; cols.setdefault(hchar(k2), []).append(k2)
    gens = [(ms, mt, r) for ms in monos(Vs, A - 1) for mt in monos(Vt, B - 1) for r in range(len(bil))]
    random.shuffle(gens)
    target = {h: len(cl) + 40 for h, cl in cols.items()}
    blocks = []; full = True
    rows = {h: [] for h in cols}
    pos = [0]
    def fill(limit):
        while pos[0] < len(gens) and any(len(rows[h]) < limit[h] for h in rows):
            ms, mt, r = gens[pos[0]]; pos[0] += 1
            row = {}
            for (a, c), v in bil[r].items():
                k2 = tuple(sorted(ms + (a,))) + tuple(sorted(mt + (c,)))
                row[k2] = (row.get(k2, 0) + v) % p
            h = hchar(next(iter(row)))
            if len(rows[h]) < limit[h]: rows[h].append(row)
    fill(target)
    for h, cl in cols.items():
        ci = {m: i for i, m in enumerate(cl)}
        rk = 0
        for attempt in range(2):
            rw = rows[h]
            if not rw: break
            M = nmod_mat(len(rw), len(cl), p)
            for r, row in enumerate(rw):
                for m, v in row.items(): M[r, ci[m]] = v
            rk = M.rank()
            if rk == len(cl) or attempt == 1: break
            fill({hh: (10**9 if hh == h else len(rows[hh])) for hh in rows})
        blocks.append((h, len(rows[h]), len(cl), rk))
        if rk < len(cl): full = False
    rec = {'exp': 'bilinear-pair', 'prime': p, 'time': time.strftime('%Y-%m-%dT%H:%M:%S'), 'd': 7, 'pair': list(key),
           'n_s': len(Vs), 'n_t': len(Vt), 'n_bilinear_forms': len(bil), 'bilinear_rank': rk1, 'n_cross_monomials': len(cross), 'bidegree': [A, B],
           'blocks_h_rows_cols_rank': blocks, 'result': 'NO ZERO PAIR (empty in P x P)' if full else 'not closed at this bidegree', 'seconds': round(time.time() - t0, 1),
           'note': 'cross term H(c_s,c_t) of the degree-7 Hessian relaxation; emptiness of {H(c_s,c_t)=0} in P(V_s) x P(V_t) via bidegree Macaulay rank per h1-block'}
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    print(key, len(Vs), len(Vt), 'forms', len(bil), 'rank', rk1, '/', len(cross), rec['result'], f"{rec['seconds']}s", flush=True)
print('BILINEAR DONE', flush=True)
