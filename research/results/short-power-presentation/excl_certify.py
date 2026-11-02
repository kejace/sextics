"""Exact certificates for the block-elimination steps of the degree-d relaxation (excl lens).

Replays the block elimination (scheme comp, iterated passes) on the cached mod-p pair table and, for every step that
eliminated unknowns, certifies dim_Q span(other pair polynomials mod J_{2d}) = dim_p by exact rational dependency
certificates (excl_lib.exact_dependencies, verified with fmpq_mpoly).  Needs excl_kernel_D7.json for d = 7.
Usage: python3 excl_certify.py <d> <e|o> [--pidx 0] [--passes 3]
"""
import sys, os, json, time, math
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import excl_lib as W
import excl_blocks as B
from hodge import core
from flint import nmod_mat

d = int(sys.argv[1]); parsarg = sys.argv[2]
par = {'e': 0, 'o': 1}[parsarg]
pidx = int(sys.argv[sys.argv.index('--pidx') + 1]) if '--pidx' in sys.argv else 0
passes = int(sys.argv[sys.argv.index('--passes') + 1]) if '--passes' in sys.argv else 3
p = W.PRIMES30[pidx]
job = os.environ.get('W2JOB', f'excl_certify_d{d}{parsarg}')
F, A, q, rec = W.load()
t00 = time.time()
cache = os.path.join(W.RES, f'excl_blocks_d{d}{parsarg}_p{pidx}.pairs.json')
C = json.load(open(cache))
unk = [tuple(u) if u[0] == 'lam' else (u[0], u[1], tuple(u[2])) for u in C['unknowns']]
neq = C['neq']
pairvec = {(int(k.split(',')[0]), int(k.split(',')[1])): {int(c): v for c, v in vec.items()} for k, vec in C['pairvec'].items()}
nu = len(unk)
kern_q = []
if d >= 7:
    ex = json.load(open(os.path.join(W.RES, 'excl_kernel_D7.json')))[str(par)]
    kern_q.append((core.parse(ex['D']), [core.parse(s) for s in ex['N']]))
exact_pairs = B.build_exact_pairs(F, A, d, kern_q, [par])
assert set(exact_pairs) >= set(pairvec), 'exact pair table does not cover the cached pairs'
print(f'loaded: {nu} unknowns, {len(pairvec)} pairs, exact pairs {len(exact_pairs)} [{time.time()-t00:.0f}s]', flush=True)
blocks = {}
for i, u in enumerate(unk):
    b = 'lam' if u[0] == 'lam' else f'w{u[1]}'
    blocks.setdefault(b, []).append(i)
order = sorted(blocks, key=lambda b: (b != 'lam', b))   # lam FIRST: its certificate needs only w-pair dependencies
alive = set(range(nu))
certs = []
for ps in range(passes):
    progress = False
    for b in order:
        Bl = [i for i in blocks[b] if i in alive]
        if not Bl:
            continue
        Bset = set(Bl)
        inB = [k for k in pairvec if k[0] in Bset and k[1] in Bset]
        others = [k for k in pairvec if (k[0] in alive and k[1] in alive) and not (k[0] in Bset and k[1] in Bset)]
        rowsO = [pairvec[k] for k in others]
        flat = [0] * (len(rowsO) * neq)
        for r, vec in enumerate(rowsO):
            for c, v in vec.items():
                flat[r * neq + c] = v
        MO = nmod_mat(len(rowsO), neq, flat, p) if rowsO else None
        if MO is not None:
            rO = MO.rank(); Psi, npsi = MO.nullspace()
            Psi = nmod_mat(neq, npsi, [int(Psi[i, k]) for i in range(neq) for k in range(npsi)], p)
        else:
            rO, npsi = 0, neq
            Psi = nmod_mat(neq, neq, [1 if i == j else 0 for i in range(neq) for j in range(neq)], p)
        flatB = [0] * (len(inB) * neq)
        for r, k in enumerate(inB):
            for c, v in pairvec[k].items():
                flatB[r * neq + c] = v
        MB = nmod_mat(len(inB), neq, flatB, p)
        Qr, rq = (MB * Psi).transpose().rref()
        local = {}
        eqs = []
        for i in range(rq):
            eq = {}
            for c, k in enumerate(inB):
                v = int(Qr[i, c])
                if v:
                    a = local.setdefault(k[0], len(local)); bb = local.setdefault(k[1], len(local))
                    eq[(min(a, bb), max(a, bb))] = v
            eqs.append(eq)
        visible = sorted(local, key=lambda i: local[i])
        if rq == 0:
            continue
        empty, hinfo = B.hilbert_verdict({'unknowns': [None] * len(local), 'eqs': eqs}, os.path.join(W.RES, f'excl_certify_d{d}{parsarg}_{b}_p{ps}'), p, '4')
        print(f'[pass {ps} block {b}] quadrics {rq} in {len(visible)} visible, m-primary {empty}, others span {rO}', flush=True)
        if not empty:
            continue
        # exact certificate: dim_Q span(others mod J) == rO
        t1 = time.time()
        polys = [exact_pairs[k] for k in others]
        dep = W.exact_dependencies(F, polys, 2 * d, 1, primes=W.PRIMES31 + W.PRIMES30 + W.PRIMES29, max_primes=12)
        ok = dep.get('exact_equals_modp', False) and dep.get('rank_modp_span_mod_J') == rO
        certs.append({'pass': ps, 'block': b, 'eliminated': len(visible), 'others': len(others), 'rank_modp': rO, 'dep': {k: v for k, v in dep.items() if k != 'primes'}, 'ok': ok})
        W.wlog(exp='block-elimination-exact-certificate', d=d, pars=parsarg, scheme='comp', block=b, pass_=ps, prime=dep.get('primes'),
               eliminated=len(visible), others=len(others), rank_modp_others=rO, result='exact' if ok else 'INCOMPLETE',
               n_deps_verified=dep.get('n_deps_verified'), n_deps_modp=dep.get('n_deps_modp'), seconds=round(time.time() - t1, 1), job=job,
               note='dim_Q of the span of the other pair polynomials modulo J_{2d} equals the mod-p dimension => the annihilator functionals of this step are reductions of rational ones')
        print(f'   exact certificate: {ok} ({dep.get("n_deps_verified")}/{dep.get("n_deps_modp")} dependencies) [{time.time()-t1:.0f}s]', flush=True)
        for i in visible:
            alive.discard(i)
        progress = True
    if not progress:
        break
print('remaining', len(alive), 'all certificates ok:', all(c['ok'] for c in certs), f'[{time.time()-t00:.0f}s]', flush=True)
W.wlog(exp='block-elimination-certified', d=d, pars=parsarg, prime=p, remaining=len(alive), steps=len(certs), result='all-exact' if certs and all(c['ok'] for c in certs) else 'incomplete', seconds=round(time.time() - t00, 1), job=job)
