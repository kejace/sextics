"""Non-Jacobian E1-kernel elements of the K3-sector target, tested with the full E2 system in the (U,V) frame using the
engine's character-restricted (mu6) solver.  Everything is done with 6-variable dict polys in UV coordinates
(slots 0..2 = U, 3..5 = V); the certificate identity is GL-equivariant so a UV-frame certificate is a certificate.
F-hat = (f(U)+f(V))/2,  A_K3-hat = (U0U1U2)^2/64 - (V0V1V2)^2/64 (both mu6-eigenvectors), D mu6-homogeneous.
(1) exact rational basis of K4 := ker(eta: S_4 -> R_{f,10}) (quartics g with eta*g in J(f)); structured D's:
    D_sym = sum_k g_k(U) g_k(V) (char 4, swap-symmetric), D_alt, and rank-1 controls;
(2) kernel of [A_K3] on (R/J)_8 per character (UV frame, engine mult_kernel + character split) and random
    homogeneous elements;  (3) symmetry.solve_fixed_D_symmetric for each candidate (mod p, CRT lift if consistent).
usage: k3kernel_uv.py <tag> [nrandom]
"""
import sys, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from fractions import Fraction as Fr
from flint import fmpq_mat, fmpz_mat
from hodge import core, jacobian, fixed_d, symmetry
tag = sys.argv[1]; nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 2
only_kernel = '--only-kernel' in sys.argv
OUT = f'k3kernel-uv-{tag}.jsonl'
uv = json.load(open('/tmp/hodge/results/polarized-k3-remainder/uv_forms.json'))
Fh = core.parse(uv['Fhat']); f = core.parse(uv['f'])
AK3 = {(2,2,2,0,0,0): Fr(1,64), (0,0,0,2,2,2): Fr(-1,64)}
p = core.PRIMES[0]
def vdeg(e): return e[3]+e[4]+e[5]
def emit(rec):
    rec['time'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['tag'] = tag
    with open(OUT, 'a') as fh: fh.write(json.dumps(rec) + '\n')
    print(json.dumps(rec)[:500], flush=True)
group = symmetry.diagonal_group(Fh)
emit({'exp': 'uv-diagonal-group', 'order': group.order, 'size': len(group), 'generators': group.generators})

def test_D(label, D, family):
    t0 = time.time()
    chs = {vdeg(e) % 6 for e in D}
    if len(chs) != 1 or core.parity_of(D) is None:
        emit({'exp': 'fixedD-AK3-uv', 'label': label, 'family': family, 'result': 'skipped: not mu6-homogeneous / mixed V-parity', 'chars': sorted(chs)}); return
    inJ, _ = jacobian.in_jacobian(Fh, D, p=p)
    okK, infoK = fixed_d.e1_solvable(Fh, AK3, 1, D, p=p)
    rec = {'exp': 'fixedD-AK3-uv', 'label': label, 'family': family, 'd': core.degree(D), 'char': chs.pop(), 'prime': p,
           'D_in_J': bool(inJ), 'e1_AK3': bool(okK), 'terms': len(D)}
    if not okK:
        rec.update({'result': 'E1 fails for A_K3', 'seconds': round(time.time()-t0, 1)}); emit(rec); return
    status, data = symmetry.solve_fixed_D_symmetric(Fh, AK3, 1, D, verbose=False, group=group, compressed=True, max_primes=4)
    rec.update({'result': status, 'rank': data.get('rank'), 'nullity': data.get('nullity'), 'rows': data.get('rows'), 'cols': data.get('cols'),
                'full_cols': data.get('full_cols'), 'subgroup_size': data.get('subgroup_size'), 'seconds': round(time.time()-t0, 1)})
    if status == 'solved':
        json.dump({'frame': 'UV', 'd': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q']), 'target': 'A_K3-hat'},
                  open(f'WITNESS-AK3-UV-{tag}-{label}.json', 'w'))
        rec['result'] = 'SOLVED -- K3-sector target served (check witness file!)'
    emit(rec)

# (1) exact K4 = ker(eta: S_4 -> R_10)  via J(f)_10 membership: eta*g in span of (z^m f_i, |m|=5)
tern = lambda P: {(e[0], e[1], e[2]): v for e, v in P.items()}
f3 = tern(f); fi3 = [tern(core.pderiv(f, i)) for i in range(3)]
def tmons(k): return [(a, b, k-a-b) for a in range(k+1) for b in range(k+1-a)]
m10 = tmons(10); i10 = {m: i for i, m in enumerate(m10)}
def tmul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2]); out[e] = out.get(e, 0) + c1*c2
    return {e: v for e, v in out.items() if v}
gens = [tmul({m: Fr(1)}, fi3[i]) for i in range(3) for m in tmons(5)]
m4 = tmons(4)
cols = gens + [tmul({m: Fr(1)}, {(2,2,2): Fr(1)}) for m in m4]
den = 1
import math
Mz = fmpz_mat(len(m10), len(cols))
for j, c in enumerate(cols):
    for e, v in c.items(): Mz[i10[e], j] = int(Fr(v) * 7)   # 7f integral
Xz, nulz = Mz.nullspace()
K4 = []
ng = len(gens)
for k in range(nulz):
    g = {m4[i]: int(Xz[ng + i, k]) for i in range(len(m4)) if int(Xz[ng + i, k])}
    if g: K4.append(g)
# independent quartic parts
Gm = fmpz_mat(len(K4), len(m4)) if K4 else None
for r, g in enumerate(K4):
    for i, m in enumerate(m4): Gm[r, i] = g.get(m, 0)
rk = Gm.rank() if K4 else 0
emit({'exp': 'K4-exact', 'dim_ker_eta_S4_to_R10': rk, 'nullspace_total': nulz, 'note': 'quartics g (exact integer coefficients) with eta*g in J(f)_10; J(f)_4 = 0 so this is ker(eta) on R_4',
      'sample': [core.to_str({(m[0],m[1],m[2],0,0,0): v for m, v in g.items()}) for g in K4[:3]]})
# reduce to an independent set
basis = []
Bm = []
for g in K4:
    row = [g.get(m, 0) for m in m4]
    if fmpz_mat(Bm + [row]).rank() > len(Bm): Bm.append(row); basis.append(g)
def toU(g): return {(m[0], m[1], m[2], 0, 0, 0): Fr(v) for m, v in g.items()}
def toV(g): return {(0, 0, 0, m[0], m[1], m[2]): Fr(v) for m, v in g.items()}
Dsym = {}; Dalt = {}
for k, g in enumerate(basis):
    t = core.mul(toU(g), toV(g)); Dsym = core.add(Dsym, t); Dalt = core.add(Dalt, core.scale(t, (-1)**k))
if not only_kernel: test_D('K4-sym-pairing', Dsym, 'sum_k g_k(U) g_k(V), g_k basis of ker(eta|R_4) (char 4)')
if not only_kernel: test_D('K4-alt-pairing', Dalt, 'alternating-sign pairing')
if len(basis) >= 2 and not only_kernel:
    g0, g1 = basis[0], basis[1]
    test_D('K4-g0g1+g1g0', core.add(core.mul(toU(g0), toV(g1)), core.mul(toU(g1), toV(g0))), 'symmetrized rank-2 tensor')
    test_D('K4-g0g0', core.mul(toU(g0), toV(g0)), 'rank-1 g0(U)g0(V) (reducible control)')
# random swap-symmetric combinations within K4 (x) K4
random.seed(3)
for t in (range(nrand) if not only_kernel else []):
    D = {}
    for i, gi in enumerate(basis):
        for j, gj in enumerate(basis):
            if i <= j:
                c = random.randrange(-2, 3)
                if c:
                    term = core.add(core.mul(toU(gi), toV(gj)), core.mul(toU(gj), toV(gi))) if i != j else core.mul(toU(gi), toV(gi))
                    D = core.add(D, core.scale(term, c))
    if D: test_D(f'K4-random-sym-{t}', D, 'random symmetric tensor in ker(eta)_4 (x) ker(eta)_4')
# (2) full kernel of [A_K3] on (R/J)_8 in the UV frame, split by character; random homogeneous elements
for par in (0, 1):
    k = jacobian.mult_kernel(Fh, AK3, 8, p=p, parity=par)
    bych = {}
    for b in k['basis']:
        for ch in range(6):
            comp = {e: v for e, v in b.items() if vdeg(e) % 6 == ch}
            if comp: bych.setdefault(ch, []).append(comp)
    dims = {}
    for ch, vecs in bych.items():
        mons = sorted({e for v in vecs for e in v}); idx = {m: i for i, m in enumerate(mons)}
        from flint import nmod_mat
        Mm = nmod_mat(len(vecs), len(mons), [vecs[r].get(m, 0) for r in range(len(vecs)) for m in mons], p)
        dims[ch] = Mm.rank()
    emit({'exp': 'e1kernel-AK3-uv-d8', 'V-parity': par, 'kernel_mod_J': k['dim'], 'ring_dim': k['ring_dim'], 'J_dim': k['jacobian_dim'], 'char_dims': dims})
    for ch, vecs in sorted(bych.items()):
        for t in range(min(nrand, 2)):
            D = {}
            for v in vecs:
                c = random.randrange(-2, 3)
                if c: D = core.add(D, core.scale({e: Fr(x if x <= p//2 else x - p) for e, x in v.items()}, c))
            if D: test_D(f'ker-par{par}-char{ch}-random{t}', D, f'random combination of character-{ch} components of E1-kernel reps (symmetric residues mod p)')
emit({'exp': 'done'})
