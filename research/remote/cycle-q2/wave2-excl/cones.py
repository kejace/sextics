"""Structured families at degree d = 7..12 (lens excl, lane cycle-q2): cones with vertex a coordinate line
(forms in 4 of the 6 variables) and vertex a coordinate plane (forms in 3 variables).
For each coordinate subset V and each G-character sector: E1 kernel {D in S_d(V)^kappa : A D in J}.
A zero kernel excludes the whole family in that sector (rigorous at a smooth prime); nonzero kernels are
tested element-wise (Farkas) and by the sector Hessian relaxation (Singular dim).
usage: python3 wave2-excl/cones.py --deg 8 [--prime P] [--workers 8] [--subsets 4,3] [--no-farkas]
"""
import sys, json, time, argparse, os, itertools, random, re, subprocess
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib, hess6, sectors
from hess6 import key_of, poly_key, SectorReducer
from sectors import solve_E1, Hpart, singular_dim
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials, reduce_mod
from flint import nmod_mat

X = ['X0','X1','X2','Y0','Y1','Y2']

def analyze(args):
    subset, d, p, do_farkas, nrandom = args
    t0 = time.time()
    F, A, q, G = w2lib.load(); gens = G.generators
    F2 = scale(F, 2); Fi = [reduce_mod(pderiv(F2, i), p) for i in range(6)]
    Ap = reduce_mod(A, p); kA = poly_key(A, gens)
    degE = d + 12; degH = 2 * d + 6
    mb = {}
    for k in sorted({d, degE, degE - 5, degH, degH - 5, max(d - 5, 0)}):
        mb[k] = {}
        for m in monomials(k):
            mb[k].setdefault(key_of(m, gens), []).append(m)
    # monomials supported on the subset, grouped by sector
    sub_mons = {}
    for m in monomials(d):
        if all(m[i] == 0 for i in range(6) if i not in subset):
            sub_mons.setdefault(key_of(m, gens), []).append(m)
    name = ''.join(X[i] for i in subset)
    red18 = {}
    out = {'exp': f'cones-d{d}', 'd': d, 'prime': p, 'subset': name, 'nsectors': len(sub_mons), 'dimV': sum(len(v) for v in sub_mons.values()), 'sectors': {}}
    total_K = 0; elems = []
    for kD, Dmons in sorted(sub_mons.items()):
        kervecs, Rd, kAD = sectors.sector_kernel(Fi, Ap, kA, kD, Dmons, gens, p, mb, d, degE, red18)
        k = len(kervecs)
        total_K += k
        srec = {'dimS': len(Dmons), 'dimK': k}
        if k == 0:
            out['sectors'][str(kD)] = srec; continue
        # J-part within the cone family: those kernel vectors in J_d
        inJ = sum(1 for v in kervecs if Rd.normal_form(v) == {})
        srec['nJ_basis_in_J'] = inJ
        tag = f'cone{len(subset)}_{name}_d{d}_p{p}_k' + ''.join(map(str, kD))
        tests = [(f'{tag}_b{i}', v, 'basis') for i, v in enumerate(kervecs)]
        rng = random.Random(hash((name, kD)) & 0xffff)
        if k >= 2:
            for r_ in range(nrandom):
                Dr = {}
                for v in kervecs: Dr = add(Dr, scale(v, rng.randrange(1, p), p), p)
                tests.append((f'{tag}_r{r_}', Dr, 'random'))
        for tname, Dt, kind in tests:
            t1 = time.time()
            Dt = {e: (v if v <= p // 2 else v - p) for e, v in Dt.items()}
            Dstr = core.to_str(Dt)
            try:
                if do_farkas:
                    fr = w2lib.farkas(F, A, q, Dt, p, G, save=f'{w2lib.PACKET}/certs/{tname}.json', Dstr=Dstr[:300])
                    st = fr['status']
                    e = {'exp': f'cones-d{d}-element', 'd': d, 'prime': p, 'subset': name, 'sector': list(kD), 'kind': kind, 'tag': tname, 'D': Dstr[:400],
                         'result': 'CONSISTENT mod p (lift needed!)' if st == 'consistent-mod-p' else ('excluded; exact Farkas certificate' if st == 'certified' else 'inconsistent mod p; Farkas FAILED'),
                         'rank': f"{fr.get('rank_M')}/{fr.get('rank_aug')} of {fr['cols']}", 'bits_max': fr.get('bits_max'), 'sha256_y': fr.get('sha256_y'), 'seconds': round(time.time() - t1, 1)}
                else:
                    o, S, r = w2lib.solve_modp(F, A, q, Dt, p, G)
                    e = {'exp': f'cones-d{d}-element', 'd': d, 'prime': p, 'subset': name, 'sector': list(kD), 'kind': kind, 'tag': tname, 'D': Dstr[:400],
                         'result': 'CONSISTENT mod p (lift needed!)' if o['consistent'] else 'inconsistent mod p', 'rank': f"{o['rank']}/{o['cols']}", 'seconds': round(time.time() - t1, 1)}
            except Exception as ex:
                e = {'exp': f'cones-d{d}-element', 'd': d, 'subset': name, 'sector': list(kD), 'tag': tname, 'result': 'ERROR ' + repr(ex)[:300]}
            e['job'] = os.environ.get('W2JOB', 'cones'); w2lib.log(**e); elems.append(e)
        if k >= 2:
            Ns = []
            for v in kervecs:
                N = solve_E1(Fi, reduce_mod(mul(Ap, v, p), p), kAD, gens, p, mb, degE); assert N is not None; Ns.append(N)
            redH = {}; quad = {}; inv2 = pow(2, -1, p)
            for a_ in range(k):
                for b_ in range(a_, k):
                    Hab = scale(Hpart(kervecs[a_], Ns[b_], kervecs[b_], Ns[a_], p), inv2, p)
                    Hab = {e_: v % p for e_, v in Hab.items() if v % p}
                    if not Hab: continue
                    kH = poly_key(Hab, gens)
                    if kH not in redH: redH[kH] = SectorReducer(Fi, degH, kH, gens, p, mb[degH], mb[degH - 5])
                    nf = redH[kH].normal_form(Hab)
                    for e_, v in nf.items():
                        quad.setdefault((kH, e_), {})[(a_, b_)] = v if a_ == b_ else (2 * v) % p
            vars_ = [f'c{i+1}' for i in range(k)]
            polys = ['+'.join(f'{c}*{vars_[a_]}*{vars_[b_]}' if a_ != b_ else f'{c}*{vars_[a_]}^2' for (a_, b_), c in sorted(t.items())) for t in quad.values()]
            os.makedirs(f'{w2lib.PACKET}/sectors', exist_ok=True)
            try:
                dm, vd, txt = singular_dim(polys, k, p, tag)
            except Exception as ex:
                dm, vd = None, None
            srec.update({'relax_nquadrics': len(polys), 'relax_dim_affine': dm, 'relax_vdim': vd})
        out['sectors'][str(kD)] = srec
    out['total_dimK'] = total_K
    out['result'] = (f'E1-admissible cone family of dim {total_K} in {sum(1 for s in out["sectors"].values() if s["dimK"])} sectors; ' +
                     ('all sector relaxations: ' + ','.join(str(s.get('relax_dim_affine')) for s in out['sectors'].values() if s['dimK'] >= 2) if any(s['dimK'] >= 2 for s in out['sectors'].values()) else '') +
                     f"; element tests: {sum(1 for e in elems if 'excluded' in e['result'])} excluded(Farkas), {sum(1 for e in elems if 'CONSISTENT' in e['result'])} consistent, {sum(1 for e in elems if 'ERROR' in e['result'] or 'FAILED' in e['result'])} errors")
    if total_K == 0:
        out['result'] = f'EXCLUDED: no E1-admissible form of degree {d} in variables {name} (zero kernel in every sector)'
    out['seconds'] = round(time.time() - t0, 1); out['job'] = os.environ.get('W2JOB', 'cones')
    w2lib.log(**out)
    return out, elems

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--deg', type=int, required=True); ap.add_argument('--prime', type=int, default=w2lib.P1)
    ap.add_argument('--workers', type=int, default=8); ap.add_argument('--subsets', default='4,3'); ap.add_argument('--no-farkas', action='store_true'); ap.add_argument('--random', type=int, default=1)
    a = ap.parse_args()
    sizes = [int(s) for s in a.subsets.split(',')]
    jobs = [(sub, a.deg, a.prime, not a.no_farkas, a.random) for sz in sizes for sub in itertools.combinations(range(6), sz)]
    print(f'd={a.deg}: {len(jobs)} coordinate subsets', flush=True)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    allrec = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(analyze, j) for j in jobs]
        for f in as_completed(futs):
            out, elems = f.result(); allrec.append(out)
            print(out['subset'], 'dimV', out['dimV'], '|', out['result'][:160], round(out['seconds'], 1), flush=True)
    os.makedirs(f'{w2lib.PACKET}/sectors', exist_ok=True)
    json.dump(allrec, open(f'{w2lib.PACKET}/sectors/cones_d{a.deg}_p{a.prime}.json', 'w'), indent=1)
    print('ALL DONE', flush=True)
