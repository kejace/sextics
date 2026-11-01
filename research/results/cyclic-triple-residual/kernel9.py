"""Extract degree-9 kernel elements in the small sectors (2,5),(4,0),(4,1),(4,2),(4,4) and test them individually (exact)."""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
def log(e):
    with open(f'{WD}/log.jsonl', 'a') as f: f.write(json.dumps(e) + '\n')
F, A, q, rec = core.load_instance('cyclic-triple-residual')
d = 9
primes = PRIMES[:4]
want = [(2,5), (4,0), (4,1), (4,2), (4,4), (0,0)]
res_by_p = {}
for p in primes:
    res, gens = e1_kernel(F, A, d, p=p)
    res_by_p[p] = res
chi = make_chi(gens)
Fz = integral(F); dF = [ideriv(Fz, i) for i in range(6)]
def canon(p, c, vecs):
    Dm = [m for m in mons(d) if chi(m) == c]; idx = {m: i for i, m in enumerate(Dm)}
    Jrows = []
    for i in range(6):
        ci = cadd(c, chi(E6[i]))
        for m in mons(d - 5):
            if chi(m) != ci: continue
            row = [0]*len(Dm)
            for e, v in dF[i].items():
                k = tuple(a+b for a, b in zip(m, e)); row[idx[k]] = (row[idx[k]] + v) % p
            Jrows.append(row)
    J = nmod_mat(len(Jrows), len(Dm), [x for r in Jrows for x in r], p)
    Jr, rk = J.rref()
    piv = []
    for i in range(rk):
        j = 0
        while int(Jr[i, j]) == 0: j += 1
        piv.append(j)
    out = []
    for vec in vecs:
        v = [0]*len(Dm)
        for m, x in vec.items(): v[idx[m]] = x % p
        for i, j in enumerate(piv):
            if v[j]:
                f = v[j] * pow(int(Jr[i, j]), -1, p) % p
                for k in range(len(Dm)):
                    if int(Jr[i, k]): v[k] = (v[k] - f*int(Jr[i, k])) % p
        if any(v): out.append(v)
    Mo = nmod_mat(len(out), len(Dm), [x for r in out for x in r], p)
    Mr, rk2 = Mo.rref()
    return Dm, [[int(Mr[i, j]) for j in range(len(Dm))] for i in range(rk2)]
D8 = core.parse(json.load(open(f'{WD}/kernel_d8_0_4_0.json'))['D'])
for c in want:
    cans = {p: canon(p, c, res_by_p[p][c]['vectors']) for p in primes}
    Dm = cans[primes[0]][0]; k = len(cans[primes[0]][1])
    print('sector', c, 'kernel mod J dim', k, flush=True)
    elems_found = []
    for j in range(k):
        acc, m = cans[primes[0]][1][j], primes[0]
        for p in primes[1:]:
            acc = [core.crt_pair(a, m, b, p)[0] for a, b in zip(acc, cans[p][1][j])]; m *= p
        rat = [core.ratrecon(x, m) for x in acc]
        if any(x is None for x in rat):
            print('   ratrecon failed for element', j, flush=True); continue
        Dk = {Dm[i]: rat[i] for i in range(len(Dm)) if rat[i]}
        elems_found.append(Dk)
        Dexpr = core.to_str(Dk)
        print(f'   element {j}: parity={core.parity_of(Dk)} terms={len(Dk)}: {Dexpr[:300]}', flush=True)
        with open(f'{WD}/kernel_d9_{c[0]}_{c[1]}_{j}.json', 'w') as f: json.dump({'d': 9, 'sector': c, 'D': Dexpr}, f)
        t0 = time.time()
        out = test(F, A, 1, Dk, exact=True, verbose=False)
        ent = {'engine': 2, 'D': Dexpr[:400], 'd': 9, 'family': 'e1-kernel-element', 'e1': True, 'consistent': out['consistent'], 'rank': out['rank'],
               'nullity': out['nullity'], 'seconds': round(time.time()-t0, 1),
               'note': f"degree-9 E1-kernel element (sector {c}, canonical rref representative {j} of {k}); exact rows={out['rows']} cols={out['cols']} rank_aug={out['rank_aug']}"}
        if out.get('verified'):
            ent['witness'] = {'d': 9, 'D': Dexpr, 'N': out['N'], 'Q': out['Q']}; ent['note'] += ' | *** SOLVED ***'
            with open(f'{WD}/witness-candidate-{int(time.time())}.json', 'w') as f:
                json.dump({'instanceId': 'cyclic-triple-residual', 'family': 'e1-kernel-element', 'verified_by_core': True, **ent['witness']}, f, indent=1)
            print('*** SOLVED ***', flush=True)
        log(ent); print('      ', out['consistent'], out['rank'], out['rank_aug'], out['cols'], ent['seconds'], flush=True)
    # random combinations inside the sector kernel (mod p test), if k >= 2
    if k >= 2 and elems_found:
        random.seed(3)
        for trial in range(4):
            D = {}
            for Dk in elems_found:
                D = core.add(D, core.scale(Dk, Fraction(random.randint(-5, 5))))
            if not D: continue
            t0 = time.time(); out = test(F, A, 1, D, exact=False, verbose=False)
            ent = {'engine': 2, 'D': f'random combination in degree-9 kernel sector {c} (trial {trial})', 'd': 9, 'family': 'e1-kernel-combo', 'e1': True,
                   'consistent': out['consistent'], 'rank': out['rank'], 'nullity': out['nullity'], 'seconds': round(time.time()-t0, 1),
                   'note': f"mod p rows={out['rows']} cols={out['cols']} rank_aug={out['rank_aug']} consistent2={out.get('consistent2')}"}
            log(ent); print('   combo', trial, out['consistent'], out['rank'], out['rank_aug'], flush=True)
