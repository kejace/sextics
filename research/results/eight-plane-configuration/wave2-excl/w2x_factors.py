"""Sum-principle bookkeeping: for low-degree natural factors D_k, compute
  (i)  Ann_R(D_k)_6^{odd} = ker(.D_k : R_6^- -> R_{6+k})  (E1-level served targets modulo J; superset of served(D_k) mod J),
  (ii) for d <= 4: the full served target space fixed_d.target_image(F, q, D_k) and its dimension modulo J_6^-.
If every factor of a product has served(D_k) subset J, the product is excluded (served(D1 D2) = served(D1) + served(D2) mod J).
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
p = P_NEW
Fp, Fi = partials(p)
Kf = [reduce_mod(core.parse(s), p) for s in json.load(open(f'{WD}/K_decomp.json'))['K']]
Hf = [reduce_mod(core.parse(s), p) for s in json.load(open(f'{WD}/K_decomp.json'))['H']]
Q = [reduce_mod(core.parse(f'X{i}**2+Y{i}**2'), p) for i in range(3)]
diag = [reduce_mod(core.parse(f'X{i}**4-X{i}**2*Y{i}**2+Y{i}**4'), p) for i in range(3)]
facs = {}
for i in range(3):
    facs[f'Q{i}'] = Q[i]; facs[f'K{i}'] = Kf[i]; facs[f'H{i}'] = Hf[i]; facs[f'diag{i}'] = diag[i]; facs[f'Q{i}^2'] = mul(Q[i], Q[i], p)
    facs[f'X{i}Y{i}'] = reduce_mod(core.parse(f'X{i}*Y{i}'), p)
    facs[f'X{i}^2-Y{i}^2'] = reduce_mod(core.parse(f'X{i}**2-Y{i}**2'), p)
for i in range(3):
    for j in range(i + 1, 3):
        facs[f'Q{i}Q{j}'] = mul(Q[i], Q[j], p); facs[f'Q{i}+Q{j}'] = add(Q[i], Q[j], p); facs[f'Q{i}-Q{j}'] = sub(Q[i], Q[j], p)
facs['Q0+Q1+Q2'] = add(add(Q[0], Q[1], p), Q[2], p)
facs['Q0Q1Q2'] = mul(mul(Q[0], Q[1], p), Q[2], p)
facs['X0^2'] = reduce_mod(core.parse('X0**2'), p); facs['X0X1'] = reduce_mod(core.parse('X0*X1'), p); facs['X0Y1'] = reduce_mod(core.parse('X0*Y1'), p)
facs['Fermat'] = reduce_mod(core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6'), p)
# Ann_R(D)_6^{odd} via mult_kernel with "A" := D_k acting on degree-6 odd forms
done = set()
for l in open(f'{WD}/wave2-excl.log.jsonl'):
    try:
        e = json.loads(l)
        if e.get('exp') == 'c-factor-annihilator': done.add(e['D'])
    except Exception: pass
rows = []
for name, Dk in facs.items():
    if name in done: continue
    t0 = time.time()
    k = degree(Dk)
    out = jacobian.mult_kernel(F, {e: int(v) for e, v in Dk.items()}, 6, p=p, parity=1)   # kernel of .D_k on R_6^{odd}
    e = {'exp': 'c-factor-annihilator', 'D': name, 'd': k, 'prime': p, 'result': f"dim Ann_R({name})_6^odd = {out['dim']}", 'rank': None,
         'ann_dim': out['dim'], 'seconds': round(time.time() - t0, 1), 'job': 'local', 'note': 'E1-level served targets mod J for this factor (sum principle bookkeeping)'}
    if k <= 4:
        t1 = time.time()
        dim, basis, S = fixed_d.target_image(F, q, Dk, p=p)
        # dimension modulo J_6^{odd} (18): rank of [basis ; J_6^- gens] - 18
        gens = J_generators(Fi, 6, 1, p)
        mons = {}
        for P_ in basis + gens:
            for ee in P_: mons.setdefault(ee, len(mons))
        def rk(L):
            if not L: return 0
            flat = [0] * (len(L) * len(mons))
            for r, P_ in enumerate(L):
                for ee, v in P_.items(): flat[r * len(mons) + mons[ee]] = v % p
            return nmod_mat(len(L), len(mons), flat, p).rank()
        rJ = rk(gens); rBJ = rk(basis + gens)
        e['served_dim'] = int(dim); e['served_mod_J'] = int(rBJ - rJ); e['t_target_image'] = round(time.time() - t1, 1)
        e['result'] += f"; served(D) dim {dim}, mod J_6^- : {rBJ - rJ}"
    print(json.dumps(e), flush=True)
    wlog(e)
