"""Sector-wise analysis of E1-admissible denominators of degree d (lens excl, lane cycle-q2).

For every G-character sector kappa of S_d:
  K_kappa = {D in S_d^kappa : A D in J}  (E1 kernel), J_kappa = J_d cap S_d^kappa, complement = K/J (non-Jacobian classes).
  (1) dims; (2) full E1+E2 sector solve (engine v2) + exact Farkas certificate for every complement basis element
      (reduced normal form w.r.t. J) and, if dim K >= 2, for a few random points of the family (screening);
  (3) Hessian relaxation restricted to the sector family: quadrics in dim K_kappa variables, Singular dim/vdim
      (dim 0 => the whole semi-invariant family of this sector is excluded; dim 1 => finitely many candidate rays).
usage: python3 wave2-excl/sectors.py --deg 7 [--prime P] [--max-sectors N] [--no-farkas] [--random 2] [--workers W]
"""
import sys, json, time, argparse, os, random, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib, hess6
from hess6 import key_of, poly_key, SectorReducer
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials, reduce_mod
from flint import nmod_mat

def sector_kernel(Fi, Ap, kA, kD, Dmons, gens, p, mb, d, degE, red18):
    kAD = tuple((x + y) % 6 for x, y in zip(kA, kD))
    if kAD not in red18:
        red18[kAD] = SectorReducer(Fi, degE, kAD, gens, p, mb[degE], mb[degE - 5])
    R = red18[kAD]
    std = R.standard; sidx = {m: i for i, m in enumerate(std)}
    K = len(Dmons); Cn = len(std)
    if Cn == 0:
        kervecs = [{m: 1} for m in Dmons]
    else:
        flat = [0] * (Cn * K)
        for j, m in enumerate(Dmons):
            nf = R.normal_form(shift(Ap, m, 1))
            for e, v in nf.items():
                flat[sidx[e] * K + j] = v
        M = nmod_mat(Cn, K, flat, p)
        X, nullity = M.nullspace()
        kervecs = [{Dmons[j]: int(X[j, c]) for j in range(K) if int(X[j, c])} for c in range(nullity)]
    Rd = SectorReducer(Fi, d, kD, gens, p, mb[d], mb[max(d - 5, 0)] if d >= 5 else {})
    return kervecs, Rd, kAD

def solve_E1(Fi, AD, kAD, gens, p, mb, degE):
    cols = []; colp = []
    for k in range(6):
        kk = key_of(next(iter(Fi[k])), gens)
        need = tuple((x - y) % 6 for x, y in zip(kAD, kk))
        for m in mb[degE - 5].get(need, []):
            cols.append((k, m)); colp.append(shift(Fi[k], m, 1))
    rows = {}
    for P in colp:
        for e in P: rows.setdefault(e, len(rows))
    for e in AD: rows.setdefault(e, len(rows))
    Rn, Cc = len(rows), len(cols)
    flat = [0] * (Rn * (Cc + 1))
    for j, P in enumerate(colp):
        for e, val in P.items(): flat[rows[e] * (Cc + 1) + j] = val % p
    for e, val in AD.items(): flat[rows[e] * (Cc + 1) + Cc] = val % p
    M = nmod_mat(Rn, Cc + 1, flat, p)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j <= Cc and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    if piv and piv[-1] == Cc:
        return None
    N = [{} for _ in range(6)]
    for i, j in enumerate(piv):
        val = int(Mr[i, Cc])
        if val:
            k, m = cols[j]; N[k][m] = val
    return N

def Hpart(Da, Nb, Db, Na, p):
    out = {}
    divNb = {}; divNa = {}
    for k in range(6):
        divNb = add(divNb, pderiv(Nb[k], k, p), p); divNa = add(divNa, pderiv(Na[k], k, p), p)
    out = add(out, mul(Da, divNb, p), p); out = add(out, mul(Db, divNa, p), p)
    for k in range(6):
        out = core.sub(out, mul(Nb[k], pderiv(Da, k, p), p), p)
        out = core.sub(out, mul(Na[k], pderiv(Db, k, p), p), p)
    return out

def singular_dim(polys, nvars, p, tag):
    vars_ = [f'c{i+1}' for i in range(nvars)]
    path = f'{w2lib.PACKET}/sectors/{tag}.sing'
    with open(path, 'w') as f:
        f.write(f'ring r = {p}, ({",".join(vars_)}), dp;\nideal I = ' + (',\n'.join(polys) if polys else '0') + ';\n')
        f.write('ideal G = std(I); "DIM", dim(G); if (dim(G) <= 1) { "VDIM", vdim(G); } "LEAD"; lead(G); quit;\n')
    r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=600)
    m = re.search(r'DIM (\d+)', r.stdout)
    v = re.search(r'VDIM (-?\d+)', r.stdout)
    return (int(m.group(1)) if m else None), (int(v.group(1)) if v else None), r.stdout[-400:]

def analyze_sector(args):
    (kD, Dmons, d, p, do_farkas, nrandom) = args
    t0 = time.time()
    F, A, q, G = w2lib.load()
    gens = G.generators
    F2 = scale(F, 2)
    Fi = [reduce_mod(pderiv(F2, i), p) for i in range(6)]
    Ap = reduce_mod(A, p); kA = poly_key(A, gens)
    degE = d + 12; degH = 2 * d + 6
    mb = {}
    for k in sorted({d, degE, degE - 5, degH, degH - 5, max(d - 5, 0)}):
        mb[k] = {}
        for m in monomials(k):
            mb[k].setdefault(key_of(m, gens), []).append(m)
    red18 = {}
    kervecs, Rd, kAD = sector_kernel(Fi, Ap, kA, kD, Dmons, gens, p, mb, d, degE, red18)
    dimK = len(kervecs); dimJ = Rd.rank
    rec = {'exp': f'sector-d{d}', 'd': d, 'prime': p, 'sector': list(kD), 'dimS': len(Dmons), 'dimK': dimK, 'dimJ': dimJ, 'complement': dimK - dimJ}
    os.makedirs(f'{w2lib.PACKET}/sectors', exist_ok=True)
    tag = f'd{d}_p{p}_k' + ''.join(map(str, kD))
    elems = []
    if dimK == 0:
        rec['result'] = 'no E1-admissible D in this sector'; rec['seconds'] = round(time.time() - t0, 1)
        w2lib.log(**rec); return rec, elems
    # complement reps: normal forms mod J, row reduce
    nfs = [Rd.normal_form(v) for v in kervecs]
    std = Rd.standard; sidx = {m: i for i, m in enumerate(std)}
    comp = []
    if std and any(nfs):
        flat = [0] * (len(nfs) * len(std))
        for r_, nf in enumerate(nfs):
            for e, v in nf.items(): flat[r_ * len(std) + sidx[e]] = v
        M = nmod_mat(len(nfs), len(std), flat, p)
        Mr, rk = M.rref()
        for i in range(rk):
            comp.append({std[j]: int(Mr[i, j]) for j in range(len(std)) if int(Mr[i, j])})
    assert len(comp) == dimK - dimJ, (len(comp), dimK, dimJ)
    # full E1+E2 tests of complement reps (+ random points)
    tests = []
    for i, Dc in enumerate(comp):
        tests.append((f'{tag}_c{i}', Dc, 'complement-basis'))
    rng = random.Random(hash(kD) & 0xffff)
    if dimK >= 2:
        for r_ in range(nrandom):
            Dr = {}
            for v in kervecs:
                Dr = add(Dr, scale(v, rng.randrange(1, p), p), p)
            tests.append((f'{tag}_r{r_}', Dr, 'random-point'))
    for name, Dt, kind in tests:
        t1 = time.time()
        Dt = {e: (v if v <= p // 2 else v - p) for e, v in Dt.items()}   # symmetric lift to Z for the exact machinery
        Dstr = core.to_str(Dt)
        try:
            if do_farkas:
                fr = w2lib.farkas(F, A, q, Dt, p, G, save=f'{w2lib.PACKET}/certs/{name}.json', Dstr=Dstr[:300])
                e = {'exp': f'sector-d{d}-element', 'd': d, 'prime': p, 'sector': list(kD), 'kind': kind, 'tag': name, 'D': Dstr[:400],
                     'status': fr['status'], 'rank': f"{fr.get('rank_M')} (M) / {fr.get('rank_aug')} (aug) of {fr['cols']} cols" if fr['status'] != 'consistent-mod-p' else f"consistent rank {fr['rank']}",
                     'farkas_rows': fr.get('farkas_rows'), 'bits_max': fr.get('bits_max'), 'sha256_y': fr.get('sha256_y'), 'seconds': round(time.time() - t1, 1)}
                e['result'] = 'CONSISTENT mod p (lift needed!)' if fr['status'] == 'consistent-mod-p' else ('excluded; exact Farkas certificate' if fr['status'] == 'certified' else 'inconsistent mod p; Farkas FAILED')
            else:
                out, S, r = w2lib.solve_modp(F, A, q, Dt, p, G)
                e = {'exp': f'sector-d{d}-element', 'd': d, 'prime': p, 'sector': list(kD), 'kind': kind, 'tag': name, 'D': Dstr[:400],
                     'result': 'CONSISTENT mod p (lift needed!)' if out['consistent'] else 'inconsistent mod p', 'rank': f"{out['rank']-1 if not out['consistent'] else out['rank']} / {out['cols']} cols", 'seconds': round(time.time() - t1, 1)}
        except Exception as ex:
            e = {'exp': f'sector-d{d}-element', 'd': d, 'sector': list(kD), 'tag': name, 'result': 'ERROR ' + repr(ex)[:300], 'seconds': round(time.time() - t1, 1)}
        e['job'] = os.environ.get('W2JOB', 'sectors')
        w2lib.log(**e); elems.append(e)
    # Hessian relaxation on the sector family (all of K_kappa)
    if dimK >= 2:
        Ns = []
        for v in kervecs:
            AD = reduce_mod(mul(Ap, v, p), p)
            N = solve_E1(Fi, AD, kAD, gens, p, mb, degE)
            assert N is not None
            Ns.append(N)
        redH = {}; quad = {}
        inv2 = pow(2, -1, p)
        for a_ in range(dimK):
            for b_ in range(a_, dimK):
                Hab = scale(Hpart(kervecs[a_], Ns[b_], kervecs[b_], Ns[a_], p), inv2, p)
                Hab = {e: v % p for e, v in Hab.items() if v % p}
                if not Hab: continue
                kH = poly_key(Hab, gens)
                if kH not in redH:
                    redH[kH] = SectorReducer(Fi, degH, kH, gens, p, mb[degH], mb[degH - 5])
                nf = redH[kH].normal_form(Hab)
                for e, v in nf.items():
                    quad.setdefault((kH, e), {})[(a_, b_)] = v if a_ == b_ else (2 * v) % p
        vars_ = [f'c{i+1}' for i in range(dimK)]
        polys = ['+'.join(f'{c}*{vars_[a_]}*{vars_[b_]}' if a_ != b_ else f'{c}*{vars_[a_]}^2' for (a_, b_), c in sorted(t.items())) for t in quad.values()]
        try:
            dm, vd, txt = singular_dim(polys, dimK, p, tag)
        except Exception as ex:
            dm, vd, txt = None, None, repr(ex)
        rec.update({'relax_nquadrics': len(polys), 'relax_dim_affine': dm, 'relax_vdim': vd})
        rec['result'] = ('sector family EXCLUDED by Hessian relaxation (affine dim 0)' if dm == 0 else
                         (f'relaxation leaves finitely many rays (affine dim 1, vdim {vd})' if dm == 1 else f'relaxation variety of affine dim {dm} (not decided)'))
    else:
        rec['result'] = 'dim K = 1: single element tested above (complement %d)' % (dimK - dimJ)
    rec['seconds'] = round(time.time() - t0, 1); rec['job'] = os.environ.get('W2JOB', 'sectors')
    w2lib.log(**rec)
    return rec, elems

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--deg', type=int, required=True); ap.add_argument('--prime', type=int, default=w2lib.P1)
    ap.add_argument('--max-sectors', type=int, default=0); ap.add_argument('--no-farkas', action='store_true'); ap.add_argument('--random', type=int, default=2)
    ap.add_argument('--workers', type=int, default=8); ap.add_argument('--only-nonJ', action='store_true')
    a = ap.parse_args()
    F, A, q, G = w2lib.load(); gens = G.generators; d = a.deg; p = a.prime
    byk = {}
    for m in monomials(d):
        byk.setdefault(key_of(m, gens), []).append(m)
    sectors = sorted(byk.items())
    if a.max_sectors: sectors = sectors[:a.max_sectors]
    print(f'd={d}: {len(sectors)} sectors, dimS_d={len(monomials(d))}', flush=True)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    allrec = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(analyze_sector, (kD, Dmons, d, p, not a.no_farkas, a.random)) for kD, Dmons in sectors]
        for f in as_completed(futs):
            rec, elems = f.result(); allrec.append(rec)
            print(rec['sector'], 'dimS', rec['dimS'], 'K', rec['dimK'], 'J', rec['dimJ'], '|', rec['result'], '|', [e['result'][:40] for e in elems], round(rec['seconds'], 1), flush=True)
    os.makedirs(f'{w2lib.PACKET}/sectors', exist_ok=True)
    json.dump(allrec, open(f'{w2lib.PACKET}/sectors/summary_d{d}_p{p}.json', 'w'), indent=1)
    tot = sum(r['dimK'] for r in allrec); totJ = sum(r['dimJ'] for r in allrec)
    print(f'TOTAL dim K_{d} = {tot}, dim J_{d} = {totJ}, non-Jacobian = {tot - totJ}', flush=True)
    w2lib.log(exp=f'sector-d{d}-summary', d=d, prime=p, result=f'dim K_{d} = {tot} (J {totJ}, non-J {tot-totJ}); sectors: ' + '; '.join(f"{r['sector']}:{r['result'][:50]}" for r in allrec if r['dimK']), job=os.environ.get('W2JOB','sectors'))
