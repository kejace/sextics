"""Wave-2 step (1): E1-admissible denominator spaces per degree d = 4..12 and character sector for four-plus-two-residual.
For each (d, sector): #monomials, dim J(F)_d, dim {D : A D in J(F)} (mod p, two primes) and numerically (SVD),
kernel mod J, F-multiples dim; plus the Thom-Sebastiani tensor prediction from the 4-variable core kernels
ker(a : R(H)_k -> R(H)_{k+6}) (k = 0..12, both (Y0,Y1)-parities) and the binary ring R(B).  Output: wave2-dims.json(l)."""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '4'; os.environ['OMP_NUM_THREADS'] = '4'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
from w2lib import *
from w2lib import _nmod_from_polys, _dense
import multiprocessing as mp

OUT = f'{LANE}/wave2-dims.jsonl'
inst = Instance()
F, A = inst.F, inst.A
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}
Bb = {e: c for e, c in F.items() if e[2] or e[5]}
CORE = [0, 1, 3, 4]

def cmon(deg, par=None):
    return [e for e in monomials(deg, par) if e[2] == 0 and e[5] == 0]

def core_kernel_dim(k, par, p):
    """dim ker(a: R(H)_k -> R(H)_{k+6}) on (Y0,Y1)-parity par, plus dim J(H)_k and dim R(H)_k of that parity (mod p)."""
    Hp = reduce_mod(H, p); ap = reduce_mod(A, p)
    Hi = [pderiv(Hp, i, p) for i in CORE]
    mons = cmon(k, par)
    tpar = (1 + par) % 2
    Jcols = []
    for idx, i in enumerate(CORE):
        for m2 in cmon(k + 1, (tpar - (0 if i < 3 else 1)) % 2):
            Jcols.append(shift(Hi[idx], m2, 1, p))
    Acols = [shift(ap, m, 1, p) for m in mons]
    rJ = _nmod_from_polys(Jcols, p)[0].rank() if Jcols else 0
    rAJ = _nmod_from_polys(Acols + Jcols, p)[0].rank()
    dimKer = len(mons) - (rAJ - rJ)
    Jk = []
    if k >= 5:
        for idx, i in enumerate(CORE):
            for m3 in cmon(k - 5, (par - (0 if i < 3 else 1)) % 2):
                Jk.append(shift(Hi[idx], m3, 1, p))
    dimJk = _nmod_from_polys(Jk, p)[0].rank() if Jk else 0
    return {'k': k, 'par': par, 'nmons': len(mons), 'dimKer': dimKer, 'dimJk': dimJk, 'dimKer_modJ': dimKer - dimJk, 'ring_dim': len(mons) - dimJk}

def binary_ring_dims():
    """dim R(B)_j by (X2-parity, Y2-parity): standard monomials X2^a Y2^b, a,b <= 4."""
    out = {}
    for a in range(5):
        for b in range(5):
            out.setdefault((a + b, a % 2, b % 2), 0)
            out[(a + b, a % 2, b % 2)] += 1
    return out

def sector_job(args):
    d, lab = args
    t0 = time.time()
    p1, p2 = core.PRIMES[0], core.PRIMES[1]
    nD = len(inst.sector_monomials(d, lab))
    rec = {'kind': 'dims', 'd': d, 'sector': list(lab), 'nmons': nD}
    if nD == 0:
        rec['dimFam'] = 0; rec['seconds'] = 0; return rec
    m1 = e1_family_dim_modp(inst, d, lab, p1)
    m2 = e1_family_dim_modp(inst, d, lab, p2)
    rec['dimFam_modp'] = [m1['dim'], m2['dim']]
    rec['dimJ'] = jacobian_dim_sector(inst, d, lab, p1)
    rec['dimFmult'] = fmult_dim_sector(inst, d, lab)
    rec['rank_gens'] = m1['rank_gens']; rec['rows'] = m1['rows']
    fam = e1_family(inst, d, lab)
    rec['dimFam_num'] = fam['K']; rec['sv_kept_max'] = fam['sv_kept_max']; rec['sv_dropped_min'] = fam['sv_dropped_min']
    rec['dimFam'] = m1['dim']
    rec['dimKernel_modJ'] = m1['dim'] - rec['dimJ']
    rec['parityY'] = core.parity_of({inst.sector_monomials(d, lab)[0]: 1})
    rec['seconds'] = round(time.time() - t0, 1)
    return rec

if __name__ == '__main__':
    p = core.PRIMES[0]
    # 1. core kernels k = 0..12
    core_dims = {}
    for k in range(0, 13):
        for par in (0, 1):
            r = core_kernel_dim(k, par, p)
            r2 = core_kernel_dim(k, par, core.PRIMES[1])
            r['agree_two_primes'] = (r['dimKer'] == r2['dimKer'])
            core_dims[(k, par)] = r
            log_line(OUT, dict(kind='core-kernel', **r))
            print('[core]', r, flush=True)
    # 2. tensor prediction: dim ker([A]) on (R/J)_d^{Y-parity} = sum_j sum K_core(d-j, pc) * dim R(B)_j^{(x2par, y2par)} with pc + y2par = parY
    rb = binary_ring_dims()
    pred = {}
    for d in range(0, 13):
        for parY in (0, 1):
            tot = 0
            for (j, xp, yp), nb in rb.items():
                k = d - j
                if k < 0: continue
                pc = (parY - yp) % 2
                tot += core_dims[(k, pc)]['dimKer_modJ'] * nb
            pred[(d, parY)] = tot
            log_line(OUT, {'kind': 'tensor-pred', 'd': d, 'parityY': parY, 'dimKernel_modJ_pred': tot})
            print(f'[tensor] d={d} parY={parY}: predicted dim ker[A] mod J = {tot}', flush=True)
    # 3. per sector d = 4..12
    tasks = []
    for d in range(4, 13):
        for lab in sorted(inst.labels_in_degree(d)):
            tasks.append((d, lab))
    tot = {}
    with mp.get_context('fork').Pool(8) as pool:
        for rec in pool.imap_unordered(sector_job, tasks):
            log_line(OUT, rec)
            print('[sector]', {k: rec[k] for k in rec if k not in ('kind',)}, flush=True)
            key = (rec['d'], rec.get('parityY'))
            if rec['nmons']:
                t = tot.setdefault(key, {'dimFam': 0, 'dimJ': 0, 'dimKernel_modJ': 0, 'nmons': 0, 'nsectors': 0})
                t['dimFam'] += rec['dimFam']; t['dimJ'] += rec['dimJ']; t['dimKernel_modJ'] += rec['dimKernel_modJ']; t['nmons'] += rec['nmons']; t['nsectors'] += 1
    for (d, parY), t in sorted(tot.items()):
        t.update({'kind': 'dims-total', 'd': d, 'parityY': parY, 'tensor_pred_kernel_modJ': pred.get((d, parY))})
        log_line(OUT, t)
        print('[total]', t, flush=True)
    print('[dims] done', flush=True)
