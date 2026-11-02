"""E1-kernel by character sector for d = 8..12; extract and test every kernel element (mod J) as a denominator."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import sect
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
def log(e):
    with open(f'{WD}/log.jsonl', 'a') as f: f.write(json.dumps(e) + '\n')
F, A, q, rec = core.load_instance('cyclic-triple-residual')
p0, p1 = PRIMES[0], PRIMES[1]
dmin = int(sys.argv[1]) if len(sys.argv) > 1 else 8
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
for d in range(dmin, dmax + 1):
    t0 = time.time()
    res, gens = e1_kernel(F, A, d, p=p0)
    tot = {'kerA': sum(r['kerA_in_S'] for r in res.values()), 'dimJ': sum(r['dimJ'] for r in res.values()),
           'kermodJ': sum(r['ker_mod_J'] for r in res.values()), 'dimS': sum(r['dim_S'] for r in res.values())}
    secs = round(time.time() - t0, 1)
    print(f'd={d} {tot} {secs}s', flush=True)
    for c, r in res.items():
        if r['ker_mod_J']:
            print('  sector', c, {k: v for k, v in r.items() if k != 'vectors'}, flush=True)
    log({'engine': 2, 'D': f'E1-kernel sweep (sector) degree {d}', 'd': d, 'family': 'e1-kernel', 'e1': None, 'consistent': None,
         'rank': tot['kerA'], 'nullity': tot['kermodJ'], 'seconds': secs,
         'note': f"dim ker([A]: S_d -> R_(d+6)) = {tot['kerA']}; dim J_d = {tot['dimJ']}; kernel mod J = {tot['kermodJ']} (all parities, by G-sector; p={p0})",
         'sectors': {str(c): {k: v for k, v in r.items() if k != 'vectors'} for c, r in res.items() if r['ker_mod_J']}})
    # extract kernel elements modulo J: for each sector with ker_mod_J > 0, reduce kernel vectors modulo J_d (mod p),
    # then rational-reconstruct with a second prime.
    for c, r in res.items():
        if not r['ker_mod_J']: continue
        # recompute at p1 for reconstruction
        res1, _ = e1_kernel(F, A, d, p=p1, G=None)
        r1 = res1[c]
        # reduce mod J: build J_d^c basis (mod p), then rref [J | kernel vectors] to find vectors independent of J,
        # canonical representative: project onto complement of J's pivot monomials.
        for (pp, rr) in ((p0, r), (p1, r1)):
            pass
        # Simpler canonical form: rref of [J-basis ; kernel-basis] rows (as polynomials in the monomial basis),
        # the rows of the rref beyond rank(J) give representatives; same pivot structure at both primes -> CRT.
        def canon(pp, rr):
            Fz = integral(F); dF = [ideriv(Fz, i) for i in range(6)]
            chi = make_chi(gens)
            Dm = [m for m in mons(d) if chi(m) == c]; idx = {m: i for i, m in enumerate(Dm)}
            Jrows = []
            if d >= 5:
                for i in range(6):
                    ci = cadd(c, chi(E6[i]))
                    for m in mons(d - 5):
                        if chi(m) != ci: continue
                        row = [0]*len(Dm)
                        for e, v in dF[i].items():
                            k = tuple(a+b for a, b in zip(m, e)); row[idx[k]] = (row[idx[k]] + v) % pp
                        Jrows.append(row)
            J = nmod_mat(len(Jrows), len(Dm), [x for row in Jrows for x in row], pp) if Jrows else None
            rJ = J.rank() if J is not None else 0
            # rref of J to get pivot columns
            piv = set()
            if J is not None:
                Jr, rk = J.rref()
                for i in range(rk):
                    j = 0
                    while j < len(Dm) and int(Jr[i, j]) == 0: j += 1
                    piv.add(j)
            out = []
            for vec in rr['vectors']:
                v = [0]*len(Dm)
                for m, x in vec.items(): v[idx[m]] = x % pp
                # reduce v modulo row space of J: eliminate pivot columns
                if J is not None:
                    Jr, rk = J.rref()
                    for i in range(rk):
                        j = 0
                        while int(Jr[i, j]) == 0: j += 1
                        if v[j]:
                            f = v[j] * pow(int(Jr[i, j]), -1, pp) % pp
                            for k in range(len(Dm)):
                                if int(Jr[i, k]): v[k] = (v[k] - f*int(Jr[i, k])) % pp
                if any(v):
                    # normalize: first nonzero = 1
                    j0 = next(k for k in range(len(Dm)) if v[k]); inv = pow(v[j0], -1, pp)
                    v = [(x*inv) % pp for x in v]
                    out.append(v)
            # independent set
            if out:
                Mo = nmod_mat(len(out), len(Dm), [x for row in out for x in row], pp)
                Mr, rk = Mo.rref()
                out = [[int(Mr[i, j]) for j in range(len(Dm))] for i in range(rk)]
            return Dm, out
        Dm0, v0 = canon(p0, r); Dm1, v1 = canon(p1, r1)
        assert Dm0 == Dm1 and len(v0) == len(v1)
        for k in range(len(v0)):
            m = p0 * p1
            xs = [core.crt_pair(a, p0, b, p1)[0] for a, b in zip(v0[k], v1[k])]
            rat = [core.ratrecon(x, m) for x in xs]
            if any(x is None for x in rat):
                print('  rational reconstruction failed in sector', c, flush=True); continue
            Dk = {Dm0[j]: rat[j] for j in range(len(Dm0)) if rat[j]}
            Dexpr = core.to_str(Dk)
            par = core.parity_of(Dk)
            print(f'  kernel element d={d} sector={c} parity={par} terms={len(Dk)}:\n    {Dexpr}', flush=True)
            with open(f'{WD}/kernel_d{d}_{"_".join(map(str,c))}_{k}.json', 'w') as f:
                json.dump({'d': d, 'sector': c, 'D': Dexpr, 'parity': par}, f)
            # test as denominator (exact)
            try:
                out = test(F, A, q, Dk, exact=True, verbose=True)
            except AssertionError as ex:
                out = {'error': str(ex)}
            ent = {'engine': 2, 'D': Dexpr[:400], 'd': d, 'family': 'e1-kernel-element', 'e1': True, 'consistent': out.get('consistent'),
                   'rank': out.get('rank', -1), 'nullity': out.get('nullity', -1), 'seconds': out.get('t_rank', 0),
                   'note': f"degree-{d} E1-kernel element (sector {c}, parity {par}); sector-reduced exact test rows={out.get('rows')} cols={out.get('cols')} rank_aug={out.get('rank_aug')}" + (' | *** SOLVED ***' if out.get('verified') else '')}
            if out.get('verified'):
                ent['witness'] = {'d': d, 'D': Dexpr, 'N': out['N'], 'Q': out['Q']}
                with open(f'{WD}/witness-candidate-{int(time.time())}.json', 'w') as f:
                    json.dump({'instanceId': 'cyclic-triple-residual', 'family': 'e1-kernel-element', 'verified_by_core': True, **ent['witness']}, f, indent=1)
                print('\n*** SOLVED ***', Dexpr, flush=True)
            log(ent)
