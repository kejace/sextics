"""Wave-2 (excl lens): block elimination of the E2 relaxation by annihilator functionals.

The relaxed E2 condition for D = sum_a lam_a D_a + w(F) is  Xi(u) := sum_{k=(u1,u2)} u1*u2*P_k  in J_{2d},
with P_k the pair polynomials (quadratic coefficients).  Working with the normal forms v_k of the P_k in
R_{2d} = S_{2d}/J_{2d} (both parities), for a block B of unknowns put
    O_B = span{ v_k : k not contained in B },   Psi_B = annihilator of O_B.
For every solution u, Xi(u) = 0 in R_{2d}, so  psi(Xi(u)) = psi( sum_{k subset B} u1 u2 v_k ) = 0  for psi in Psi_B:
a system of quadrics in the B-unknowns only.  If its ideal is m-primary (in the unknowns that actually occur), every
solution has those unknowns = 0; they are dropped and the next block is processed on the reduced pair table.
The residual system (remaining unknowns) is solved directly at the end.

Rigor (mod p -> Q): each step needs dim_Q(O_B) = dim_p(O_B), certified by exact rational dependency certificates
among the corresponding pair polynomials modulo J (excl_lib.exact_dependencies), and the Hilbert-function
inequality for the quadric ideals.  --exact runs those certificates (needs exact kernel data for d >= 7).

Usage: python3 excl_blocks.py <d> <e|o|m> [--pidx k] [--scheme comp|xy|parity] [--threads T] [--exact] [--tag T]
"""
import sys, time, json, os, math, re
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import excl_lib as W
import excl_relax as R
from hodge import core
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, reduce_mod
from flint import nmod_mat
from fractions import Fraction


def hilbert_verdict(sysd, base, p, threads):
    """msolve leading ideal + Singular on a small quadric system; returns (empty, info)."""
    nu = len(sysd['unknowns'])
    W.write_msolve(sysd, base + '.ms', p)
    W.write_singular(sysd, base + '.sing', p)
    ms = R.run_msolve(base + '.ms', base + '.ms.out', threads)
    mons = R.parse_msolve_leading(base + '.ms.out')
    degs = [sum(k for _, k in m) for m in mons]
    maxdeg = max(degs) if degs else 0
    tmax = 0
    for t in range(1, maxdeg + 4):
        cost = sum(math.comb(t - dm + nu - 1, nu - 1) for dm in degs if dm <= t)
        if cost > 3e7:
            break
        tmax = t
    h = R.hilbert_from_leading(mons, nu, tmax)[0] if tmax >= 1 else None
    sg = R.run_singular(base + '.sing', base + '.sing.out', nu) if nu <= 40 else {}
    empty = (h is not None and 0 in h[1:]) or str(sg.get('dim', '')) == '0'
    return bool(empty), {'msolve_hilbert': h, 'msolve_seconds': ms['seconds'], 'singular_dim': sg.get('dim'),
                         'singular_vdim': sg.get('vdim'), 'singular_hilbert': sg.get('hilbert_function'),
                         'singular_seconds': sg.get('seconds'), 'unknowns': nu, 'quadrics': len(sysd['eqs'])}


def main():
    d = int(sys.argv[1]); parsarg = sys.argv[2]
    pars = {'e': [0], 'o': [1], 'm': [0, 1]}[parsarg]
    pidx = int(sys.argv[sys.argv.index('--pidx') + 1]) if '--pidx' in sys.argv else 0
    scheme = sys.argv[sys.argv.index('--scheme') + 1] if '--scheme' in sys.argv else 'comp'
    threads = sys.argv[sys.argv.index('--threads') + 1] if '--threads' in sys.argv else '4'
    tag = sys.argv[sys.argv.index('--tag') + 1] if '--tag' in sys.argv else ''
    do_exact = '--exact' in sys.argv
    p = W.PRIMES30[pidx]
    job = os.environ.get('W2JOB', f'excl_blocks_d{d}{parsarg}')
    base = os.path.join(W.RES, f'excl_blocks_d{d}{parsarg}_{scheme}_p{pidx}{tag}')
    F, A, q, rec = W.load()
    Fp, Ap = reduce_mod(F, p), reduce_mod(A, p)
    t00 = time.time()
    # ---- kernel elements (exact if available for d = 7) ---------------------------------------------------
    kern = []
    kern_exact = []
    exact_ok = True
    path7 = os.path.join(W.RES, 'excl_kernel_D7.json')
    for par in pars:
        if d < 7:
            continue
        if d == 7 and os.path.exists(path7):
            ex = json.load(open(path7))[str(par)]
            Dq = core.parse(ex['D']); Nq = [core.parse(s) for s in ex['N']]
            kern_exact.append((Dq, Nq))
            kern.append((reduce_mod(Dq, p), [reduce_mod(n, p) for n in Nq]))
        else:
            exact_ok = False
            info, basis = W.e1_kernel(F, A, d, par, p, want_basis=True)
            for Da in basis or []:
                Na, _ = W.membership_mod_p(Fp, mul(Ap, Da, p), p)
                kern.append((Da, Na))
    print(f'kernel elements: {len(kern)} (exact available: {exact_ok}) [{time.time()-t00:.0f}s]', flush=True)
    # ---- full pair table in normal-form coordinates (cached) ----------------------------------------------------
    cache = os.path.join(W.RES, f'excl_blocks_d{d}{parsarg}_p{pidx}.pairs.json')
    if os.path.exists(cache):
        C = json.load(open(cache))
        unk = [tuple(u) if u[0] == 'lam' else (u[0], u[1], tuple(u[2])) for u in C['unknowns']]
        neq = C['neq']
        pairvec = {(int(k.split(',')[0]), int(k.split(',')[1])): {int(c): v for c, v in vec.items()} for k, vec in C['pairvec'].items()}
        print(f'pair table loaded from cache {cache}', flush=True)
    else:
        full = W.build_relaxation(F, A, d, kern, pars, p, slice_euler=True, verbose=False)
        unk = full['unknowns']
        neq = len(full['eqs'])
        pairvec = {}
        for e_idx, eq in enumerate(full['eqs']):
            for k, v in eq.items():
                pairvec.setdefault(k, {})[e_idx] = v
        with open(cache, 'w') as f:
            json.dump({'unknowns': [list(u) if u[0] == 'lam' else [u[0], u[1], list(u[2])] for u in unk], 'neq': neq,
                       'pairvec': {f'{k[0]},{k[1]}': vec for k, vec in pairvec.items()}}, f)
    nu = len(unk)
    print(f'unknowns {nu}, pairs {len(pairvec)}, equation coordinates {neq} [{time.time()-t00:.0f}s]', flush=True)

    def blockname(u):
        if u[0] == 'lam':
            return 'lam'
        j = u[1]
        if scheme == 'comp':
            return f'w{j}'
        if scheme == 'xy':
            return 'wX' if j < 3 else 'wY'
        if scheme == 'parity':
            pj = 0 if j < 3 else 1
            return 'w+' if (core.yparity(u[2]) + pj) % 2 == 0 else 'w-'
        if scheme == 'mono':
            return 'm' + ''.join(str(e) for e in u[2])
        if scheme == 'fine':
            return 'u' + ','.join(str(e) for e in u[2]) + f'_{j}'
        if scheme == 'xtype':   # component block split by the X-degree of the monomial
            return f'w{j}x{sum(u[2][:3])}'
        raise ValueError(scheme)

    blocks = {}
    for i, u in enumerate(unk):
        blocks.setdefault(blockname(u), []).append(i)
    order = sorted(blocks, key=lambda b: (b == 'lam', b))   # lam last
    if scheme == 'fine':
        order = sorted(blocks, key=lambda b: (b == 'lam', unk[blocks[b][0]][1] < 3, b))   # Y-components first, then X, lam last
    print('blocks:', {b: len(v) for b, v in blocks.items()}, 'order', order, flush=True)
    alive = set(range(nu))
    steps = []
    maxpass = int(sys.argv[sys.argv.index('--passes') + 1]) if '--passes' in sys.argv else 4
    sched = []
    for ps in range(maxpass):
        sched += [(ps, b) for b in order]
    done_blocks = set()
    last_progress_pass = 0
    for ps, b in sched:
        if ps > last_progress_pass + 1:
            break
        if b in done_blocks:
            continue
        t0 = time.time()
        B = [i for i in blocks[b] if i in alive]
        if not B:
            done_blocks.add(b)
            continue
        Bset = set(B)
        inB = [k for k in pairvec if k[0] in Bset and k[1] in Bset]
        others = [k for k in pairvec if (k[0] in alive and k[1] in alive) and not (k[0] in Bset and k[1] in Bset)]
        # annihilator of span(others)
        rowsO = [pairvec[k] for k in others]
        if rowsO:
            flat = [0] * (len(rowsO) * neq)
            for r, vec in enumerate(rowsO):
                for c, v in vec.items():
                    flat[r * neq + c] = v
            MO = nmod_mat(len(rowsO), neq, flat, p)
            rO = MO.rank()
            Psi, npsi = MO.nullspace()
            Psi = nmod_mat(neq, npsi, [int(Psi[i, k]) for i in range(neq) for k in range(npsi)], p)
        else:
            rO, npsi = 0, neq
            Psi = nmod_mat(neq, neq, [1 if i == j else 0 for i in range(neq) for j in range(neq)], p)
        # quadric coefficients for in-block pairs
        flatB = [0] * (len(inB) * neq)
        for r, k in enumerate(inB):
            for c, v in pairvec[k].items():
                flatB[r * neq + c] = v
        MB = nmod_mat(len(inB), neq, flatB, p)
        Qm = (MB * Psi).transpose()      # npsi x #inB
        Qr, rq = Qm.rref()
        rB = MB.rank()
        # independent quadrics -> system in block unknowns that occur
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
        invisible = [i for i in B if i not in local]
        info = {'block': b, 'block_size': len(B), 'pairs_in_block': len(inB), 'dim_span_in_block': rB,
                'dim_span_others': rO, 'functionals': npsi, 'quadrics': rq, 'visible': len(visible),
                'invisible': [unk[i] for i in invisible]}
        if rq == 0:
            empty = False
            hinfo = {}
        else:
            sysd = {'unknowns': [None] * len(local), 'eqs': eqs}
            empty, hinfo = hilbert_verdict(sysd, base + f'.block_{b}', p, threads)
        info.update(hinfo)
        info['proj_empty_in_visible'] = empty
        info['seconds'] = round(time.time() - t0, 1)
        info['pass'] = ps
        if empty:
            for i in visible:
                alive.discard(i)
            info['eliminated'] = len(visible)
            last_progress_pass = ps
            if not invisible:
                done_blocks.add(b)
        steps.append(info)
        print(f'[block {b}] size {len(B)} | pairs in block {len(inB)} (span {rB}) | others span {rO} | functionals {npsi} | '
              f'quadrics {rq} in {len(visible)} visible unknowns | m-primary: {empty} | invisible {len(invisible)} [{time.time()-t0:.0f}s]', flush=True)
        W.wlog(exp='block-elimination-step', d=d, pars=parsarg, scheme=scheme, prime=p, **{k: (v if k != 'invisible' else len(v)) for k, v in info.items()},
               result='eliminated' if empty else 'not-eliminated', job=job)
    # ---- residual system -----------------------------------------------------------------------------------
    rem = sorted(alive)
    print(f'remaining unknowns after block elimination: {len(rem)}: {[unk[i] for i in rem][:40]}', flush=True)
    res_info = {'remaining': len(rem)}
    verdict = 'unknown'
    if not rem:
        verdict = 'proj-empty'
    else:
        remset = set(rem)
        kk = [k for k in pairvec if k[0] in remset and k[1] in remset]
        rowsR = [pairvec[k] for k in kk]
        flat = [0] * (len(rowsR) * neq)
        for r, vec in enumerate(rowsR):
            for c, v in vec.items():
                flat[r * neq + c] = v
        MR = nmod_mat(len(rowsR), neq, flat, p)
        # equations = columns -> basis of the span of the quadrics: rref of MR^T
        Qr, rq = MR.transpose().rref()
        local = {i: n for n, i in enumerate(rem)}
        eqs = []
        for i in range(rq):
            eq = {}
            for c, k in enumerate(kk):
                v = int(Qr[i, c])
                if v:
                    a, bb = local[k[0]], local[k[1]]
                    eq[(min(a, bb), max(a, bb))] = v
            eqs.append(eq)
        sysd = {'unknowns': [None] * len(rem), 'eqs': eqs}
        if len(rem) <= 40 or ('--residual-gb' in sys.argv and len(rem) <= 70):
            empty, hinfo = hilbert_verdict(sysd, base + '.residual', p, threads)
            res_info.update(hinfo)
            verdict = 'proj-empty' if empty else 'residual-not-empty'
        else:
            W.write_msolve(sysd, base + '.residual.ms', p)
            W.write_singular(sysd, base + '.residual.sing', p)
            res_info['quadrics'] = rq
            verdict = 'residual-too-large-written'
    print('verdict:', verdict, res_info, flush=True)
    rec = {'exp': 'block-elimination', 'd': d, 'pars': parsarg, 'scheme': scheme, 'prime': p, 'unknowns': nu,
           'steps': [{k: (v if k != 'invisible' else len(v)) for k, v in s.items()} for s in steps],
           'residual': res_info, 'result': verdict, 'seconds': round(time.time() - t00, 1), 'job': job,
           'note': 'proj-empty => every solution of the relaxed E2 system is trivial (mod p; exact certificates needed for the functionals of each step)'}
    # ---- exact certificates for the functionals (dependencies among the "others" pair polynomials mod J) -------
    if do_exact and verdict == 'proj-empty' and (d < 7 or exact_ok):
        # rebuild exact pair polynomials for all pairs, grouped by parity, and certify dim_Q span(others) per step
        print('exact certification of the annihilator steps ...', flush=True)
        kern_q = kern_exact if d >= 7 else []
        exact_pairs = build_exact_pairs(F, A, d, kern_q, pars)
        certs = []
        alive2 = set(range(nu))
        for s in steps:
            b = s['block']
            B = [i for i in blocks[b] if i in alive2]
            Bset = set(B)
            others = [k for k in exact_pairs if (k[0] in alive2 and k[1] in alive2) and not (k[0] in Bset and k[1] in Bset)]
            cert = {}
            for par in sorted({core.yparity(mo) for k in others for mo in exact_pairs[k]}):
                polys = [{mo: v for mo, v in exact_pairs[k].items() if core.yparity(mo) == par} for k in others]
                polys = [P for P in polys if P]
                dep = W.exact_dependencies(F, polys, 2 * d, par, primes=W.PRIMES31, max_primes=5)
                cert[par] = dep
            certs.append({'block': b, 'cert': cert})
            W.wlog(exp='block-elimination-exact-certificate', d=d, pars=parsarg, scheme=scheme, block=b,
                   result='exact' if all(c.get('exact_equals_modp') for c in cert.values()) else 'INCOMPLETE',
                   detail={k: {kk: vv for kk, vv in v.items() if kk != 'primes'} for k, v in cert.items()}, job=job)
            if s.get('proj_empty_in_visible'):
                for i in B:
                    if i not in [j for j in B if unk[j] in s['invisible']]:
                        alive2.discard(i)
        rec['exact_certificates'] = certs
    W.wlog(**rec)


def build_exact_pairs(F, A, d, kern_q, pars):
    """Exact rational pair polynomials (unreduced) keyed like excl_lib.build_relaxation's unknown indices."""
    Fi = [pderiv(F, i) for i in range(NV)]
    H = [[pderiv(Fi[j], l) for l in range(NV)] for j in range(NV)]
    Ai = [pderiv(A, i) for i in range(NV)]
    unknowns = [('lam', a) for a in range(len(kern_q))]
    for parD in pars:
        for (j, m) in W.field_basis(d, parD):
            if j == 0 and m[0] > 0:
                continue
            unknowns.append(('w', j, m))
    terms = {}

    def addterm(i, j, P):
        if not P:
            return
        if i > j:
            i, j = j, i
        terms[(i, j)] = add(terms.get((i, j), {}), P)

    divN = []
    for Da, Na in kern_q:
        dv = {}
        for i in range(NV):
            dv = add(dv, pderiv(Na[i], i))
        divN.append(dv)
    Dg = [[pderiv(Da, i) for i in range(NV)] for Da, _ in kern_q]
    for a, (Da, Na) in enumerate(kern_q):
        for b, (Db, Nb) in enumerate(kern_q):
            B = mul(Da, divN[b])
            for i in range(NV):
                B = sub(B, mul(Nb[i], Dg[a][i]))
            addterm(a, b, B)
    widx = [(k, u) for k, u in enumerate(unknowns) if u[0] == 'w']
    for a, (Da, Na) in enumerate(kern_q):
        HN = []
        for j in range(NV):
            s = {}
            for i in range(NV):
                s = add(s, mul(Na[i], H[i][j]))
            HN.append(s)
        for k, (_, j, m) in widx:
            P = shift(sub(sub(mul(Da, Ai[j]), mul(A, Dg[a][j])), HN[j]), m, 1)
            addterm(a, k, P)
    AH = [[mul(A, H[j][l]) for l in range(NV)] for j in range(NV)]
    for x in range(len(widx)):
        k1, (_, j, m) = widx[x]
        for y in range(x, len(widx)):
            k2, (_, l, m2) = widx[y]
            mm = tuple(u + v for u, v in zip(m, m2))
            addterm(k1, k2, shift(AH[j][l], mm, -1 if x == y else -2))
    return terms


if __name__ == '__main__':
    main()
