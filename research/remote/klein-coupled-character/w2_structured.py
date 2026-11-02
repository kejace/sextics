"""Structured-family tests (task item 4): evaluate geometrically motivated E1-admissible candidates of degree 6..12
inside their character sector: numerical fixed-D residual rel(D) via the Koszul varpro functional (one evaluation),
float E1 membership, exact e1_solvable (mod p), and exact kcc2 Koszul test for canonical candidates / any rel < 1e-8.
Families: polar x monomial (pure-block and squarefree), polar x {K,G,M,radX,radY}, polar products x variable,
radial products, principal Hessian minors (2x2: d=8, 3x3: d=12), Jhat x monomials (d=10..12).
Usage: python3 w2_structured.py [dmin dmax]"""
import sys, os, json, time, itertools
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
from fractions import Fraction
import w2fast as wf, w2lib, kcc2
from hodge import core, fixed_d
from hodge.core import NV, monomials, pderiv, mul, add, sub, scale, shift
F, A, q = wf.load()
dmin, dmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (6, 12)
Fi = [pderiv(F, i) for i in range(NV)]
X = [(0,)*i + (1,) + (0,)*(5-i) for i in range(NV)]
def mono(e): return {tuple(e): Fraction(1)}
Kb = {e: v for e, v in F.items() if sum(e[3:]) == 0}
Gb = {e: v for e, v in F.items() if sum(e[:3]) == 0}
Mb = sub(sub(F, Kb), Gb)
radX = {}; radY = {}
for i in range(3): radX = add(radX, shift(Fi[i], X[i], 1))
for i in range(3, 6): radY = add(radY, shift(Fi[i], X[i], 1))
H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
def det3(I):
    (a, b, c) = I
    def m(i, j): return H[i][j]
    t = mul(m(a,a), sub(mul(m(b,b), m(c,c)), mul(m(b,c), m(c,b))))
    t = sub(t, mul(m(a,b), sub(mul(m(b,a), m(c,c)), mul(m(b,c), m(c,a)))))
    t = add(t, mul(m(a,c), sub(mul(m(b,a), m(c,b)), mul(m(b,b), m(c,a)))))
    return t
def det2(I):
    a, b = I
    return sub(mul(H[a][a], H[b][b]), mul(H[a][b], H[b][a]))
rec_src = open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0]
Jhat = core.parse('7*X0*X2**2*(9*(Y0**6 - Y1**6 + Y2**6) - 9*(' + rec_src + ')) + 15*X1*Y2*Y0*(2*Y0**6 - 4*Y1**6 + 8*Y2**6 - 2*(' + rec_src + '))')
names = ['X0','X1','X2','Y0','Y1','Y2']
def mname(e): return '*'.join(f'{names[i]}^{e[i]}' if e[i] > 1 else names[i] for i in range(6) if e[i]) or '1'
cands = []   # (D, tag, family, canonical)
def addc(D, tag, fam, canon=False):
    if not D: return
    d = core.degree(D)
    if d < dmin or d > dmax: return
    if core.parity_of(D) is None or wf.char_of(D) is None: return
    cands.append((D, tag, fam, canon))
# polar x monomials
for deg in range(1, 8):
    for e in monomials(deg):
        pureX = sum(e[3:]) == 0; pureY = sum(e[:3]) == 0; sqfree = max(e) == 1
        if not (pureX or pureY or sqfree or (deg == 7 and min(e) >= 1)):
            continue
        for i in range(NV):
            addc(shift(Fi[i], e, 1), f'F_{i}*{mname(e)}', 'polar-x-monomial')
# polar x forms
for nm, P in (('K', Kb), ('G', Gb), ('M', Mb), ('radX', radX), ('radY', radY)):
    for i in range(NV):
        addc(mul(Fi[i], P), f'F_{i}*{nm}', 'polar-x-form', True)
# polar products x variable / quadratic monomials
for i in range(NV):
    for j in range(i, NV):
        P = mul(Fi[i], Fi[j])
        for e in monomials(1): addc(shift(P, e, 1), f'F_{i}*F_{j}*{mname(e)}', 'polar-product-x-var')
        for e in monomials(2): addc(shift(P, e, 1), f'F_{i}*F_{j}*{mname(e)}', 'polar-product-x-quad')
# radial products
for (n1, P1), (n2, P2) in itertools.combinations_with_replacement((('radX', radX), ('radY', radY), ('K', Kb), ('G', Gb), ('M', Mb)), 2):
    if 'rad' in n1 or 'rad' in n2:
        addc(mul(P1, P2), f'{n1}*{n2}', 'radial-product', True)
# Hessian minors
for I in itertools.combinations(range(NV), 2): addc(det2(I), f'Hess2{I}', 'hessian-2x2-principal', True)
for I in itertools.combinations(range(NV), 3): addc(det3(I), f'Hess3{I}', 'hessian-3x3-principal', True)
# Jhat multiples
for deg in (1, 2, 3):
    for e in monomials(deg): addc(shift(Jhat, e, 1), f'Jhat*{mname(e)}', 'jhat-x-monomial', deg == 1)
print('candidates', len(cands), flush=True)
groups = {}
for D, tag, fam, canon in cands:
    groups.setdefault((core.degree(D), core.parity_of(D), wf.char_of(D)), []).append((D, tag, fam, canon))
out = open('wave2-structured.jsonl', 'a')
done = set()
if os.path.exists('wave2-structured.jsonl'):
    for line in open('wave2-structured.jsonl'):
        try: done.add(json.loads(line)['tag'])
        except Exception: pass
for (d, par, chi), lst in sorted(groups.items()):
    lst = [c for c in lst if c[1] not in done]
    if not lst: continue
    t0 = time.time()
    sec = wf.Sector(F, A, d, par, chi)
    Bc, Bn, info = sec.e1_space()
    kvp = wf.KVP(sec, Bc, Bn) if info['K'] else None
    summary = {'n': len(lst), 'e1_fail': 0, 'spurious': 0, 'rels': [], 'exact': []}
    for D, tag, fam, canon in lst:
        rec = {'tag': tag, 'family': fam, 'd': d, 'parD': par, 'chi': chi, 'terms': len(D), 'canonical': canon}
        try:
            c = wf.poly_to_vec(D, sec.Didx)
        except KeyError:
            rec['note'] = 'monomial outside sector list'; out.write(json.dumps(rec) + '\n'); continue
        if kvp is None:
            rec['e1_float'] = False; rec['e1'] = False
        else:
            a = Bc.T @ c; resid = float(np.linalg.norm(Bc @ a - c) / max(np.linalg.norm(c), 1e-300))
            rec['e1_float'] = resid < 1e-8; rec['e1_resid'] = resid
            if rec['e1_float']:
                spur, srel = kvp.is_spurious(a); rec['F_divides'] = bool(spur)
                if spur:
                    summary['spurious'] += 1
                else:
                    f, bb, xP, r, rank = kvp.value(a, exact=True)
                    rec['rel'] = f / max(bb, 1e-300); rec['f'] = f; rec['rankT'] = rank; rec['nullityT'] = sec.nP - rank
                    summary['rels'].append(rec['rel'])
        if not rec.get('e1_float'):
            ok, einfo = fixed_d.e1_solvable(F, A, 1, D)
            rec['e1'] = bool(ok)
            if not ok: summary['e1_fail'] += 1
        need_exact = (rec.get('rel', 1.0) < 1e-8) or (canon and d <= 11 and rec.get('e1_float') and not rec.get('F_divides'))
        if need_exact:
            e = kcc2.full_test(F, A, D, tag, 'wave2-structured-' + fam, note='structured candidate (wave 2 search lens)', want_log=True)
            rec['exact'] = {'e1': e['e1'], 'consistent': e['consistent'], 'rank': e['rank'], 'nullity': e['nullity'], 'seconds': e['seconds']}
            summary['exact'].append((tag, e['consistent']))
            if e['consistent']:
                r = kcc2.lift_certificate(F, A, D, weights=True, nprimes=8)
                rec['lift'] = None if r is None else r['ok']
                if r and r['ok']:
                    wf.save_witness_candidate(F, A, D, r['N'], r['Q'])
                    print('WITNESS FOUND', tag, flush=True)
                    wf.batch_log({'lens': 'structured', 'WITNESS': True, 'd': d, 'parD': par, 'chi': chi, 'D': core.to_str(D)})
        out.write(json.dumps(rec, default=str) + '\n'); out.flush()
    rels = summary['rels']
    wf.batch_log({'lens': 'structured', 'd': d, 'parD': par, 'chi': chi, 'K': info['K'], 'n_candidates': len(lst),
                  'e1_fail': summary['e1_fail'], 'F_multiples': summary['spurious'], 'n_in_family': len(rels),
                  'best_rel': min(rels) if rels else None, 'median_rel': float(np.median(rels)) if rels else None,
                  'exact_tests': summary['exact'], 'families': sorted({c[2] for c in lst}), 'seconds': round(time.time() - t0, 1)})
    print(f'sector d={d} par={par} chi={chi}: {len(lst)} cands, e1_fail {summary["e1_fail"]}, in family {len(rels)}, best rel {min(rels) if rels else None}, exact {summary["exact"]} [{time.time()-t0:.0f}s]', flush=True)
    os.system('python3 /tmp/hodge/results/klein-coupled-character/w2_report.py > /dev/null 2>&1')
print('STRUCTURED_DONE', flush=True)
