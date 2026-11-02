"""Bilinear zero-pair test for pairs of (Z/2)^3 parity classes (s,t) of the degree-7 system.
Claim used: if c in V has character support S (set of parity classes) and a pair {s,t} in S has a difference s+t realised by
no other pair inside S, then the cross term H(c_s, c_t) must vanish on its own.  We test: is there (c_s != 0, c_t != 0) with
H(c_s, c_t) = 0 ?  The bihomogeneous ideal of the bilinear forms is tested by the bidegree-(a,b) Macaulay matrix
(full column rank in every h1-character block  =>  the variety in P(V_s) x P(V_t) is empty).  Exact mod p (python-flint).
Usage: bilinear_pairs.py quads.json a b"""
import sys, json, itertools, time, os
from flint import nmod_mat
fn = sys.argv[1]; A = int(sys.argv[2]); B = int(sys.argv[3])
LOG = os.environ.get('HODGE_W2LOG', '/tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl')
J = json.load(open(fn)); p = J['p']; chi = [tuple(c) for c in J['chi']]
quads = [{(a, b): v for a, b, v in cf} for cf in J['quads']]
n = len(chi)
cls = lambda c: (c[0], c[1], c[3])
classes = sorted({cls(c) for c in chi})
def monos(vars_, k): return list(itertools.combinations_with_replacement(vars_, k))
def hchar(mono): return sum(chi[a][2] for a in mono) % 6
results = {}
for s, t in itertools.combinations(classes, 2):
    t0 = time.time()
    Vs = [a for a in range(n) if cls(chi[a]) == s]; Vt = [a for a in range(n) if cls(chi[a]) == t]
    # bilinear forms: cross parts of the quadrics (c_a c_b with a in Vs, b in Vt)
    Ss, St = set(Vs), set(Vt)
    bil = []
    for cf in quads:
        b = {}
        for (a, c), v in cf.items():
            if (a in Ss and c in St): b[(a, c)] = v
            elif (a in St and c in Ss): b[(c, a)] = v
        if b: bil.append(b)
    # rank of the bilinear forms (as linear forms in the |Vs|*|Vt| cross monomials)
    cross = [(a, c) for a in Vs for c in Vt]; cidx = {m: i for i, m in enumerate(cross)}
    M1 = nmod_mat(len(bil), len(cross), [bil[r].get(m, 0) for r in range(len(bil)) for m in cross], p)
    rk1 = M1.rank()
    # bidegree (A,B) Macaulay: rows (m_s * m_t) * bil, m_s in S_{A-1}(Vs), m_t in S_{B-1}(Vt); columns bidegree (A,B) monomials; blocks by h-char
    cols = {}
    for ms in monos(Vs, A):
        for mt in monos(Vt, B):
            key = ms + mt; cols.setdefault(hchar(key), []).append(key)
    rows = {}
    for ms in monos(Vs, A - 1):
        for mt in monos(Vt, B - 1):
            for b in bil:
                row = {}
                for (a, c), v in b.items():
                    key = tuple(sorted(ms + (a,))) + tuple(sorted(mt + (c,)))
                    row[key] = (row.get(key, 0) + v) % p
                h = hchar(next(iter(row)))
                rows.setdefault(h, []).append(row)
    full = True; blocks = []
    for h, cl in cols.items():
        ci = {m: i for i, m in enumerate(cl)}
        rw = rows.get(h, [])
        if not rw:
            full = False; blocks.append((h, 0, len(cl), 0)); continue
        flat = [0] * (len(rw) * len(cl))
        for r, row in enumerate(rw):
            for m, v in row.items(): flat[r * len(cl) + ci[m]] = v
        rk = nmod_mat(len(rw), len(cl), flat, p).rank()
        blocks.append((h, len(rw), len(cl), rk))
        if rk < len(cl): full = False
    rec = {'exp': 'bilinear-pair', 'prime': p, 'time': time.strftime('%Y-%m-%dT%H:%M:%S'), 'd': 7, 'pair': [''.join(map(str, s)), ''.join(map(str, t))],
           'n_s': len(Vs), 'n_t': len(Vt), 'n_bilinear_forms': len(bil), 'bilinear_rank': rk1, 'n_cross_monomials': len(cross), 'bidegree': [A, B],
           'blocks_h_rows_cols_rank': blocks, 'result': 'NO ZERO PAIR (empty in P x P)' if full else 'not closed at this bidegree', 'seconds': round(time.time() - t0, 1),
           'note': 'cross term H(c_s,c_t) of the degree-7 Hessian relaxation; emptiness of {H(c_s,c_t)=0} in P(V_s) x P(V_t) via bidegree Macaulay rank per h1-block'}
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    print(rec['pair'], len(Vs), len(Vt), 'forms', len(bil), 'rank', rk1, '/', len(cross), rec['result'], f"{rec['seconds']}s", flush=True)
