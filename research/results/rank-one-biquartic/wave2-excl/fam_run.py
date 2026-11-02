"""Family exclusion driver: python3 fam_run.py <d> <prime> [--parity 0|1|all] [--cert] [--tag TAG] [--nonj FILE]
Builds the E1-kernel family basis for degree d (polar part exact; non-J part from exact JSON lifts if given),
builds the weight-2 E2 Macaulay matrix, computes its rank mod p, the quadratic forms in a lying in I_2, and checks
(mod p, then with exact rational certificates if --cert) that u^2 in I_2 for every coordinate u of V_d / Sp,
Sp = F * S_{d-6} (the non-coprime denominators)."""
import sys, json, time, os, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
import w2lib as W
from w2lib import *
from w2log import log

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('prime', type=int)
ap.add_argument('--parity', default='all'); ap.add_argument('--cert', action='store_true'); ap.add_argument('--tag', default='')
ap.add_argument('--nonj', default='', help='JSON file with exact non-J kernel generators [{"D":str,"N":[str*6]}]')
a = ap.parse_args()
d, p = a.d, a.prime
par = None if a.parity == 'all' else int(a.parity)
tag = a.tag or f'fam_d{d}_p{p}_{a.parity}'
LOGF = f'{PK}/{tag}.log.jsonl'
t0 = time.time()
basis = polar_basis(d, par)
nonj = []
if a.nonj:
    for k, rec_ in enumerate(json.load(open(a.nonj))):
        D = core.parse(rec_['D']); N = [core.parse(s) for s in rec_['N']]
        if par is not None and core.parity_of(D) != par: continue
        basis.append((('K', k), D, N)); nonj.append(k)
K = len(basis)
labidx = {b[0]: k for k, b in enumerate(basis)}
print(f'd={d} p={p} parity={a.parity}: basis size {K} (polar {K-len(nonj)}, non-J {len(nonj)})', flush=True)
FS = FamilySystem(d, basis)
info, pivots, unit_in, quad_forms, Mr = FS.rowspace_analysis(p)
# spurious subspace Sp = F*S_{d-6}; F*m = (1/6) sum_i (m z_i) F_i  <-> labels ('J', i, m*z_i) with coefficient 1/6.
sp_groups = []
if d >= 6:
    for m in monomials(d - 6):
        grp = []
        for i in range(NV):
            mm = list(m); mm[i] += 1; mm = tuple(mm)
            lab = ('J', i, mm)
            if lab in labidx: grp.append(labidx[lab])
        if len(grp) == NV: sp_groups.append(grp)
        else: print('WARNING spurious group incomplete (parity restriction?)', m, grp, flush=True)
# coordinates on V/Sp: all labels not in any spurious group (as a_k), and within each group the 5 differences a_g0 - a_gi
in_group = {k for grp in sp_groups for k in grp}
coords = []   # list of (name, linear form as dict k->coef)
for k in range(K):
    if k not in in_group:
        coords.append((f'a[{basis[k][0]}]', {k: Fraction(1)}))
for gi, grp in enumerate(sp_groups):
    for i in range(1, NV):
        coords.append((f'a[{basis[grp[0]][0]}]-a[{basis[grp[i]][0]}]', {grp[0]: Fraction(1), grp[i]: Fraction(-1)}))
def square_vec(lin):
    """column vector (dict col->Fraction) of the quadratic form (sum_k lin_k a_k)^2 in the aa-columns."""
    out = {}
    ks = sorted(lin)
    for x in range(len(ks)):
        for y in range(x, len(ks)):
            k, l = ks[x], ks[y]
            coef = lin[k] * lin[l] * (1 if k == l else 2)
            out[FS.colidx[('aa', k, l)]] = out.get(FS.colidx[('aa', k, l)], 0) + coef
    return out
C = FS.ncols
modp_ok = []
for name, lin in coords:
    modp_ok.append(reduce_by_rref(Mr, pivots, square_vec(lin), p, C))
n_ok = sum(modp_ok)
res = dict(exp='family-weight2', d=d, prime=p, parity=a.parity, K=K, nonJ=len(nonj), rank=info['rank'], rows=info['nrows'], cols=info['ncols'],
           n_aa=info['n_aa'], n_unit_aa_in_I2=info['n_unit_aa_in_I2'], n_quad_forms_in_I2=info['n_quad_forms_in_I2'],
           n_spurious_dirs=len(sp_groups), n_coords=len(coords), coord_squares_in_I2_modp=n_ok,
           failing_coords=[name for (name, _), ok in zip(coords, modp_ok) if not ok][:40],
           t_build=round(FS.build_time), t_rref=info['t_rref'], seconds=round(time.time() - t0),
           result=('all coordinate squares in I_2 mod p' if n_ok == len(coords) else f'{len(coords)-n_ok} coordinate squares NOT in I_2 mod p'))
print(json.dumps(res, default=str), flush=True)
json.dump(dict(res=res, quad_forms_modp=[{str(k): v for k, v in qf.items()} for qf in quad_forms], basis_labels=[str(b[0]) for b in basis],
               spurious_groups=[[str(basis[k][0]) for k in grp] for grp in sp_groups], coords=[(n, {str(k): str(v) for k, v in lin.items()}) for n, lin in coords]),
          open(f'{PK}/{tag}.json', 'w'), default=str)
log(dict(res, job=tag, note='weight-2 Macaulay matrix of E2 over the E1-kernel family (Euler gauge + Koszul); coordinates u on V_d/(F*S_{d-6}); u^2 in I_2 => no certificate of degree d'), path=LOGF)
if a.cert and n_ok > 0:
    targets = [(name, square_vec(lin)) for (name, lin), ok in zip(coords, modp_ok) if ok]
    certs = exact_certificate(FS, pivots, None, targets, p)
    okc = {name: (certs[name] is not None) for name, _ in targets} if certs else {}
    nok = sum(1 for v in okc.values() if v)
    # save the certificates (sparse rational row combinations) for the packet
    if certs:
        json.dump({name: ({str(r): str(v) for r, v in lam.items()} if lam else None) for name, lam in certs.items()},
                  open(f'{PK}/{tag}.certs.json', 'w'))
        json.dump(okc, open(f'{PK}/{tag}.certs_ok.json', 'w'))
    log(dict(exp='family-weight2-exact-certificates', d=d, prime=p, parity=a.parity, certified=nok, of=len(coords), job=tag,
             seconds=round(time.time() - t0), result=('PROVED: no certificate of degree %d' % d if nok == len(coords) else 'incomplete'),
             note='exact rational combinations lambda^T E = u^2 verified over Q for each coordinate u of V_d/(F*S_{d-6})'), path=LOGF)
    print('certified', nok, 'of', len(coords), flush=True)
print('done', round(time.time() - t0), 's', flush=True)
