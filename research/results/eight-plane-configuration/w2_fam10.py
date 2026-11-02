"""Wave-2 job: structured denominator families at degree d (default 10) intersected with the E1-admissible
space (kernel mod J + J_d), exactly mod p; extracts a basis of fam ∩ admissible in family coordinates and
lifts it to Q by rational reconstruction when possible; candidates go to <tag>-cands.jsonl for the sweep.

Families (from wave-1 families.py, parity-restricted): G_k(Q0,Q1,Q2) net pull-backs; C[Q,K,H,diag] subring;
(Q)^k * S_{d-2k}; products of conjugate plane-pair cones I_Pi^k ∩ I_Pibar^k (symmetric/antisymmetric);
plus here: products of two Jacobian partials, cones with a plane vertex over plane curves (X-only and Y-only
forms), Fermat-part powers, and sum-of-sixth-powers products.

usage: python3 w2_fam10.py --d 10 --pars 0,1 [--tag fam10]
"""
import sys, os, json, time, argparse, pickle, itertools
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, default=10)
ap.add_argument('--pars', default='0,1')
ap.add_argument('--tag', default='fam10')
args = ap.parse_args()
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from lib import *            # F, A, q, WD, core, fixed_d, P, Q0..Q2
import families as fam_mod
from hodge import jacobian
from hodge.core import monomials, mul, add, sub, scale, power, reduce_mod, parity_of, shift, pderiv
from flint import nmod_mat
from fractions import Fraction

p = core.PRIMES[1]
assert p % 4 == 1
OUT = f'{WD}/{args.tag}.jsonl'
CANDS = f'{WD}/{args.tag}-cands.jsonl'


def jlog(obj):
    with open(OUT, 'a') as f:
        f.write(json.dumps(obj, default=str) + '\n')


def extra_families(d, par):
    """Additional structured families (exact rational polys, reduced mod p), degree d, parity par."""
    Fi = [pderiv(F, i) for i in range(6)]
    fams = {}
    if d == 10:
        # products of two partials F_i F_j (in J, so admissible; sum principle says they serve only J-classes)
        lst = []
        for i in range(6):
            for j in range(i, 6):
                Pq = mul(Fi[i], Fi[j])
                if parity_of(Pq) == par:
                    lst.append(reduce_mod(Pq, p))
        fams['F_i*F_j'] = lst
    # cones with plane vertex: forms in X only (even) or products X-form * Y-form
    if par == 0:
        lst = [{m: 1} for m in monomials(d, 0) if all(e == 0 for e in m[3:])]
        fams['X-only (cone over plane curve, vertex Y-plane)'] = lst
        lst = [{m: 1} for m in monomials(d, 0) if all(e == 0 for e in m[:3])]
        fams['Y-only (cone, vertex X-plane)'] = lst
    # Fermat-part powers and sixth-power products
    if par == 0 and d % 2 == 0:
        fer = [core.parse(f'X{i}**6+Y{i}**6') for i in range(3)]
        lst = []
        if d == 12:
            for i in range(3):
                for j in range(i, 3):
                    lst.append(reduce_mod(mul(fer[i], fer[j]), p))
        fams['sixth-power products'] = lst
    return fams


def rank_of(polys, monos):
    if not polys:
        return 0
    R, C = len(monos), len(polys)
    flat = [0] * (R * C)
    for j, P_ in enumerate(polys):
        for e, v in P_.items():
            flat[monos[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p).rank()


def intersection_basis(famL, admL, monos):
    """Vectors x (family coordinates) with sum x_j fam_j in span(adm): nullspace of [fam | -adm]."""
    R = len(monos); C1, C2 = len(famL), len(admL)
    flat = [0] * (R * (C1 + C2))
    for j, P_ in enumerate(famL):
        for e, v in P_.items():
            flat[monos[e] * (C1 + C2) + j] = v % p
    for j, P_ in enumerate(admL):
        for e, v in P_.items():
            flat[monos[e] * (C1 + C2) + C1 + j] = (-v) % p
    M = nmod_mat(R, C1 + C2, flat, p)
    X, nul = M.nullspace()
    vecs = []
    if nul:
        proj = nmod_mat(nul, C1, [int(X[j, k]) for k in range(nul) for j in range(C1)], p)
        Pr, rk = proj.rref()
        for i in range(rk):
            vecs.append([int(Pr[i, j]) for j in range(C1)])
    return vecs


d = args.d
for par in [int(x) for x in args.pars.split(',')]:
    t0 = time.time()
    pk = f'{WD}/e1kernel_basis_d{d}_par{par}.pkl'
    kernel = pickle.load(open(pk, 'rb'))
    J = fam_mod.J_basis(d, par, p)
    adm = kernel + J
    fams = fam_mod.families(d, par, p)
    fams.update(extra_families(d, par))
    monos = {}
    for lst in [adm] + list(fams.values()):
        for P_ in lst:
            for e in P_:
                monos.setdefault(e, len(monos))
    rA = rank_of(adm, monos); rJ = rank_of(J, monos)
    jlog({'event': 'adm', 'd': d, 'par': par, 'p': p, 'dim_admissible': rA, 'dim_J': rJ, 'dim_kernel_mod_J': rA - rJ,
          'n_families': len(fams)})
    for name, lst in fams.items():
        t1 = time.time()
        rF = rank_of(lst, monos)
        if rF == 0:
            continue
        rAF = rank_of(adm + lst, monos); rJF = rank_of(J + lst, monos)
        capA = rF + rA - rAF; capJ = rF + rJ - rJF
        rec = {'event': 'family', 'd': d, 'par': par, 'family': name, 'n_gens': len(lst), 'dim_fam': rF,
               'fam_cap_admissible': capA, 'fam_cap_J': capJ, 'new_directions': capA - capJ, 'seconds': round(time.time() - t1, 1)}
        print(json.dumps(rec), flush=True)
        if capA > 0:
            vecs = intersection_basis(lst, adm, monos)
            rec['n_intersection_vectors'] = len(vecs)
            lifted = 0
            for vi, x in enumerate(vecs):
                # rational reconstruction of the family coordinates
                fr = [core.ratrecon(v, p) for v in x]
                if all(f is not None for f in fr) and max(abs(f.denominator) for f in fr) <= 10**6:
                    lifted += 1
                    rec.setdefault('lifted_coords', []).append([str(f) for f in fr if f][:12])
                    # the family generators are known only mod p here; record the coordinate vector and the
                    # generator names/indices for an exact rebuild in the sweep step
                    with open(CANDS, 'a') as f:
                        f.write(json.dumps({'d': d, 'par': par, 'family': name, 'coords': [str(f) for f in fr],
                                            'gen_index': list(range(len(lst)))}) + '\n')
                else:
                    rec.setdefault('unlifted', 0)
                    rec['unlifted'] += 1
            rec['lifted'] = lifted
        jlog(rec)
    jlog({'event': 'done', 'd': d, 'par': par, 'seconds': round(time.time() - t0, 1)})
