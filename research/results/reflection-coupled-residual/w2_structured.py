"""Structured denominator candidates (degrees 5..12, Klein-sector pure) + exact E1 screen (mod p, cached annihilators).
Usage: python3 w2_structured.py gen   -> structured-cands.jsonl
       python3 w2_structured.py screen [p_index] -> structured-e1.jsonl (one line per candidate), batch lines in the lane log."""
import sys, os, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
import w2rcr
from hodge import core, jacobian
from hodge.core import NV, mul, add, sub, scale, pderiv, power, shift
from fractions import Fraction
F, A, q = w2rcr.load()
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
P = core.parse
Fi = [pderiv(F, i) for i in range(NV)]
LDIR = w2rcr.LDIR


def restrict(Pol, keep):
    return {e: c for e, c in Pol.items() if all(e[i] == 0 for i in range(NV) if i not in keep)}


def bideg(Pol, xdeg):
    return {e: c for e, c in Pol.items() if e[0] + e[1] + e[2] == xdeg}


def det3(M):
    return sub(add(add(mul(mul(M[0][0], M[1][1]), M[2][2]), mul(mul(M[0][1], M[1][2]), M[2][0])), mul(mul(M[0][2], M[1][0]), M[2][1])),
               add(add(mul(mul(M[0][2], M[1][1]), M[2][0]), mul(mul(M[0][0], M[1][2]), M[2][1])), mul(mul(M[0][1], M[1][0]), M[2][2])))


def gen():
    C = []   # (family, poly)
    def addc(fam, Pol):
        if Pol and core.is_homogeneous(Pol) and w2rcr.sector_of(Pol) is not None and 5 <= core.degree(Pol) <= 12:
            C.append((fam, Pol))
    # A/B/C Fermat cones
    for k in (3, 4, 5, 6):
        for S in itertools.combinations(range(NV), k):
            base = [f'{V[i]}**6' for i in S]
            addc(f'fermat{k}', P(' + '.join(base)))
            if k <= 4:
                for j in range(k):
                    terms = [('-' if t == j else '+') + ' ' + b for t, b in enumerate(base)]
                    addc(f'fermat{k}-sign', P(' '.join(terms).lstrip('+ ')))
            if k == 3:
                u, v, w = [V[i] for i in S]
                for lam in (-3, 1, 3, -1):
                    addc('hesse3', P(f'{u}**6 + {v}**6 + {w}**6 + ({lam})*{u}**2*{v}**2*{w}**2'))
                addc('cubic3-sq', P(f'({u}**3 + {v}**3 + {w}**3)**2'))
                addc('cubic3-hesse-sq', P(f'({u}**3 + {v}**3 + {w}**3 - 3*{u}*{v}*{w})**2'))
    # D restrictions of F
    for k in (3, 4, 5):
        for S in itertools.combinations(range(NV), k):
            addc(f'F-restrict{k}', restrict(F, S))
    # E bidegree pieces
    parts = {xd: bideg(F, xd) for xd in (6, 4, 2, 0)}
    for r in range(1, 4):
        for S in itertools.combinations((6, 4, 2, 0), r):
            Pol = {}
            for xd in S:
                Pol = add(Pol, parts[xd])
            addc('F-bideg-' + ''.join(str(x) for x in S), Pol)
    # F minus one or two sixth powers (keeps couplings)
    for k in (1, 2):
        for S in itertools.combinations(range(NV), k):
            Pol = dict(F)
            for i in S:
                e = tuple(6 if j == i else 0 for j in range(NV)); Pol.pop(e, None)
            addc(f'F-minus{k}', Pol)
    # F radial blocks
    for k in (1, 2, 3):
        for S in itertools.combinations(range(NV), k):
            Pol = {}
            for i in S:
                Pol = add(Pol, mul(P(V[i]), Fi[i]))
            addc(f'radial{k}', Pol)
    # G Hessian minors
    H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
    for i, j in itertools.combinations(range(NV), 2):
        addc('hess2', sub(mul(H[i][i], H[j][j]), mul(H[i][j], H[i][j])))
    for S in itertools.combinations(range(NV), 3):
        M = [[H[a][b] for b in S] for a in S]
        addc('hess3', det3(M))
    # H discriminant in u = Y2^2 and H_u powers
    pu = P('5*X0*X1*(Y0**2 + 2*Y1**2)/7'); r0 = restrict(F, (0, 1, 2, 3, 4))
    disc = sub(scale(power(pu, 3), -4), scale(power(r0, 2), 27))
    addc('disc-u', disc)
    for lam in (1, -1, 4, -4):
        addc('disc-u-var', add(power(pu, 3), scale(power(r0, 2), lam)))
    Hu = P('3*Y2**4 + 5*X0*X1*(Y0**2 + 2*Y1**2)/7')
    addc('Hu2', power(Hu, 2)); addc('Hu3', power(Hu, 3))
    for k in (3, 4):
        for S in itertools.combinations(range(NV), k):
            addc(f'Hu-fermat{k}', mul(Hu, P(' + '.join(f'{V[i]}**6' for i in S))))
    addc('Hu-F60', mul(Hu, parts[6])); addc('Hu-F06', mul(Hu, parts[0])); addc('Hu-r0', mul(Hu, r0))
    # I partials and products
    for i in range(NV):
        addc('partial', Fi[i])
    for i, j in itertools.combinations_with_replacement(range(NV), 2):
        addc('partial-prod', mul(Fi[i], Fi[j]))
    # J companion denominators from Klein-compatible accepted witnesses
    for iid in ['cofactor-eight', 'cofactor-mixed', 'family-block-cube', 'family-hesse-web', 'family-shear', 'family-trace-zero',
                'twoblock-cubic', 'xtype-sparse-rational-1m142', 'xtype-strict-rational-113', 'block-q2', 'fermat-q2',
                'hesse-triangle', 'tie-singleton-3odd-strict', 'm3-coupled', 'family-sylvester-norm']:
        try:
            d, D, N, Q = core.load_witness(iid)
            addc('companion-' + iid, D)
        except Exception as ex:
            print('skip', iid, ex)
    # K linear-factor times structured (Y2-odd / Y0,Y1-odd sectors)
    for lin in ('Y2', 'Y0', 'Y1', 'X0'):
        L = P(lin)
        for i in range(NV):
            addc(f'{lin}-partial', mul(L, Fi[i]))
        for S in itertools.combinations(range(NV), 3):
            addc(f'{lin}-fermat3', mul(L, P(' + '.join(f'{V[i]}**6' for i in S))))
        addc(f'{lin}-Hu', mul(L, Hu)); addc(f'{lin}-Hu2', mul(L, power(Hu, 2)))
        addc(f'{lin}-F60', mul(L, parts[6])); addc(f'{lin}-F06', mul(L, parts[0])); addc(f'{lin}-r0', mul(L, r0))
        for k in (4,):
            for S in itertools.combinations(range(NV), k):
                addc(f'{lin}-F-restrict{k}', mul(L, restrict(F, S)))
    for quad in ('Y0*Y1', 'Y0*Y2', 'Y1*Y2', 'Y2**2', 'Y0**2 + Y1**2', 'X0*Y0', 'X0*Y2'):
        Lq = P(quad)
        for S in itertools.combinations(range(NV), 3):
            addc(f'({quad})-fermat3', mul(Lq, P(' + '.join(f'{V[i]}**6' for i in S))))
        addc(f'({quad})-Hu', mul(Lq, Hu)); addc(f'({quad})-r0', mul(Lq, r0)); addc(f'({quad})-F60', mul(Lq, parts[6]))
    # M K8-derived
    try:
        K8 = P(json.load(open(f'{LDIR}/kernel_d8_X03Y03.json'))[0]) if isinstance(json.load(open(f'{LDIR}/kernel_d8_X03Y03.json')), list) else None
    except Exception:
        K8 = None
    if K8 is None:
        try:
            info = json.load(open(f'{LDIR}/e1space-d8-chi00.json')); K8 = P(info['reps'][0])
        except Exception as ex:
            print('no K8', ex)
    if K8:
        addc('K8', K8)
        for lin in ('Y2', 'Y0', 'Y1', 'X0', 'X1', 'X2'):
            addc(f'{lin}-K8', mul(P(lin), K8))
        addc('Hu-K8', mul(Hu, K8))
    # dedupe
    seen = set(); out = []
    for fam, Pol in C:
        s = core.to_str(Pol)
        if s in seen:
            continue
        seen.add(s); out.append({'family': fam, 'd': core.degree(Pol), 'sector': list(w2rcr.sector_of(Pol)), 'D': s})
    with open(f'{LDIR}/structured-cands.jsonl', 'w') as f:
        for rec in out:
            f.write(json.dumps(rec) + '\n')
    from collections import Counter
    print(len(out), 'candidates;', Counter((r['d'], tuple(r['sector'])) for r in out))


def screen(pidx=0):
    p = w2rcr.PRIMES[pidx]
    cands = [json.loads(l) for l in open(f'{LDIR}/structured-cands.jsonl')]
    done = set()
    outp = f'{LDIR}/structured-e1.jsonl'
    if os.path.exists(outp):
        done = {json.loads(l)['D'] for l in open(outp)}
    cands = [c for c in cands if c['D'] not in done]
    cands.sort(key=lambda c: (c['d'], tuple(c['sector'])))
    cache = {}
    fh = open(outp, 'a')
    t0 = time.time(); npass = 0
    for key, group in itertools.groupby(cands, key=lambda c: (c['d'], tuple(c['sector']))):
        group = list(group); tg = time.time()
        res = []
        for c in group:
            D = P(c['D'])
            divF, _ = jacobian.divisible_by(F, D, p) if core.degree(D) >= 6 else (False, None)
            ok = w2rcr.e1_test_modp(F, A, D, p, cache)
            rec = dict(c); rec.update({'e1': bool(ok), 'F_divides': bool(divF), 'p': p, 'engine': 3})
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            res.append(rec); npass += int(ok and not divF)
        w2rcr.log_batch({'step': 'structured-e1', 'd': key[0], 'sector': list(key[1]), 'n': len(group),
                         'e1_pass': sum(int(r['e1'] and not r['F_divides']) for r in res),
                         'pass_families': [r['family'] for r in res if r['e1'] and not r['F_divides']],
                         'seconds': round(time.time() - tg, 1), 'p': p})
        print(key, len(group), 'pass', sum(int(r['e1'] and not r['F_divides']) for r in res), f'({time.time()-tg:.0f}s)', flush=True)
    print('total pass', npass, f'({time.time()-t0:.0f}s)')


if __name__ == '__main__':
    if sys.argv[1] == 'gen':
        gen()
    else:
        screen(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
