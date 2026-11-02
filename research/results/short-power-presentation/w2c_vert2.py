"""Job: for given (d, parity, family) with a nonzero restricted E1 kernel, recompute the kernel mod two primes,
rationally reconstruct the family-coefficient vectors, verify A*D in J exactly mod a third prime, and run the full
fixed-D solve (compressed, mod p; exact lift if consistent) on every non-F-multiple candidate.
Usage: python3 w2c_vert2.py <out.jsonl> <family> <d> <par> [<family> <d> <par> ...]"""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2c_lib import *
from hodge import fixed_d
import flint
try: flint.ctx.threads = int(os.environ.get('FLINT_THREADS', '8'))
except Exception: pass
out_path = sys.argv[1]
specs = [(sys.argv[i], int(sys.argv[i+1]), int(sys.argv[i+2])) for i in range(2, len(sys.argv), 3)]
F, A, q, rec = core.load_instance(IID)
p1, p2, p3 = core.PRIMES[0], core.PRIMES[1], core.PRIMES[2]
fams = {'vert': fam_vert, 'coneXM': fam_coneXM, 'coneLY': fam_coneLY}
def recon(v1, v2):
    keys = set(v1) | set(v2); out = {}
    for k in keys:
        x, mod = core.crt_pair(v1.get(k, 0), p1, v2.get(k, 0), p2)
        r = core.ratrecon(x, mod)
        if r is None: return None
        if r: out[k] = r
    return out
for fname, d, par in specs:
    fam = fams[fname](d, par)
    t0 = time.time()
    vecs1, info1 = restricted_kernel(F, A, fam, d, p1)
    vecs2, info2 = restricted_kernel(F, A, fam, d, p2)
    print('kernels', info1['kernel_dim'], info2['kernel_dim'], flush=True)
    cands = []
    for v1, v2 in zip(vecs1, vecs2):
        r = recon(v1, v2)
        if r is None:
            print('ratrecon failed', flush=True); continue
        name2P = {n: P for n, P in fam}
        D = {}
        for n, cf in r.items(): D = core.add(D, core.scale(name2P[n], cf))
        ok, info = fixed_d.e1_solvable(F, A, q, D, p=p3)
        divF, _ = jacobian.divisible_by(core.reduce_mod(F, p3), core.reduce_mod(D, p3), p3)
        cands.append((r, D, ok, bool(divF)))
        print('candidate', {k: str(x) for k, x in r.items()}, 'E1(p3)=', ok, 'F|D=', bool(divF), flush=True)
    for r, D, ok, divF in cands:
        recd = {'exp': 'vertical-candidate', 'family': fname, 'd': d, 'parD': par, 'coeffs': {k: str(x) for k, x in r.items()},
                'E1_p3': ok, 'F_divides': divF, 'job': 'construct-spp-vert2', 'prime': p3}
        if ok and not divF:
            t1 = time.time()
            try:
                out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True)
                recd.update({'consistent': out['consistent'], 'rank': out['rank'], 'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols'],
                             'result': 'solved' if out.get('exact') else ('inconsistent' if not out['consistent'] else 'lift-failed'), 'seconds': round(time.time() - t1, 1)})
                if out.get('exact'):
                    with open(f'{RES}/witness-candidate-construct.json', 'w') as f:
                        json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])}, f)
                    print('!!! WITNESS FOUND', flush=True)
            except Exception as e:
                recd.update({'result': 'error: %s' % e})
        else:
            recd['result'] = 'skipped (E1 fails or F|D)'
        with open(out_path, 'a') as f: f.write(json.dumps(recd) + '\n')
        print('DONE', recd.get('result'), flush=True)
print('ALL DONE', flush=True)
