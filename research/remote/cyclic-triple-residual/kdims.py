"""(A) D8-coset perturbations D8 + random J_8 element (mod p ranks), then kernel dims d=10..12 (all sectors)."""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
def log(e):
    with open(f'{WD}/log.jsonl', 'a') as f: f.write(json.dumps(e) + '\n')
F, A, q, rec = core.load_instance('cyclic-triple-residual')
p0 = PRIMES[0]
D8 = core.parse(json.load(open(f'{WD}/kernel_d8_0_4_0.json'))['D'])
Fz = integral(F); G = semi_invariant_subgroup(elems(Fz), [integral(A), integral(D8)]); gens = generators(G); chi = make_chi(gens)
c8 = chi_of(integral(D8), chi)
Fi = [core.pderiv(F, i) for i in range(6)]
Jbasis = []
for i in range(6):
    for m in mons(3):
        if cadd(chi(m), tuple((-x) % 6 for x in chi(E6[i]))) == c8:   # chi(m F_i) = chi(m) - e_i
            Jbasis.append(core.mul({m: Fraction(1)}, Fi[i]))
print('J_8 basis in sector of D8:', len(Jbasis), flush=True)
random.seed(11)
for trial in range(6):
    if trial < 3: coeffs = [Fraction(random.randint(-3, 3)) for _ in Jbasis]
    else: coeffs = [Fraction(random.randint(-20, 20), random.randint(1, 5)) for _ in Jbasis]
    D = D8
    for c, Jb in zip(coeffs, Jbasis):
        if c: D = core.add(D, core.scale(Jb, c))
    t0 = time.time()
    out = test(F, A, 1, D, exact=False, verbose=False)
    ent = {'engine': 2, 'D': f'D8 + random J_8 element (trial {trial})', 'd': 8, 'family': 'D8-coset', 'e1': True, 'consistent': out['consistent'],
           'rank': out['rank'], 'nullity': out['nullity'], 'seconds': round(time.time()-t0, 1), 'coeffs': [str(c) for c in coeffs],
           'note': f"sector solver mod p={p0} rows={out['rows']} cols={out['cols']} rank_aug={out['rank_aug']}; second prime consistent2={out.get('consistent2')}"}
    log(ent); print(ent['D'], out['consistent'], out['rank'], out['rank_aug'], ent['seconds'], flush=True)
for d in range(10, 13):
    t0 = time.time()
    res, gens = e1_kernel(F, A, d, p=p0, want_vectors=False)
    tot = {'kerA': sum(r['kerA_in_S'] for r in res.values()), 'dimJ': sum(r['dimJ'] for r in res.values()),
           'kermodJ': sum(r['ker_mod_J'] for r in res.values()), 'dimS': sum(r['dim_S'] for r in res.values())}
    secs = round(time.time() - t0, 1)
    print(f'd={d} {tot} {secs}s', flush=True)
    log({'engine': 2, 'D': f'E1-kernel sweep (sector) degree {d}', 'd': d, 'family': 'e1-kernel', 'e1': None, 'consistent': None,
         'rank': tot['kerA'], 'nullity': tot['kermodJ'], 'seconds': secs,
         'note': f"dim ker([A]: S_d -> R_(d+6)) = {tot['kerA']}; dim J_d = {tot['dimJ']}; kernel mod J = {tot['kermodJ']} (all parities, by G-sector; p={p0})",
         'sectors': {str(c): {k: v for k, v in r.items() if k != 'vectors'} for c, r in res.items() if r['ker_mod_J']}})
