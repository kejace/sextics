"""Helper: run fixed-D tests for four-plus-two-residual and append to log.jsonl."""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

IID = 'four-plus-two-residual'
RES = f'/tmp/hodge/results/{IID}'
LOG = f'{RES}/log.jsonl'
F, A, q, rec = core.load_instance(IID)

def log(entry):
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')

def run(D, family, note='', expr=None, fast=None, e1_only=False, p=None):
    """D: poly dict or expr string. Returns dict."""
    if isinstance(D, str):
        expr = expr or D
        D = core.parse(D)
    expr = expr or core.to_str(D)
    d = core.degree(D)
    t0 = time.time()
    if d <= 4 or e1_only:
        ok, info = fixed_d.e1_solvable(F, A, q, D, p=p)
        ent = {'D': expr, 'd': d, 'family': family, 'e1': bool(ok), 'consistent': None,
               'rank': info['rank'], 'nullity': info['cols'] - info['rank'], 'seconds': 0.0, 'note': note, 'engine': 2}
    else:
        ok = True
        ent = {'D': expr, 'd': d, 'family': family, 'e1': None, 'consistent': None,
               'rank': None, 'nullity': None, 'seconds': 0.0, 'note': note, 'engine': 2}
    if not ok or e1_only:
        ent['seconds'] = round(time.time() - t0, 2)
        ent['note'] = (note + ' | E1 fails' if not ok else note + ' | E1 only')
        log(ent); return ent
    if fast is None:
        fast = d >= 6
    if fast:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, p=p, verbose=False)
        ent.update({'consistent': bool(out['consistent']), 'rank': out['rank'], 'nullity': out['nullity'],
                    'rows': out['rows'], 'cols': out['cols']})
        if out['consistent']:
            ent['note'] = note + f" | CONSISTENT exact={out.get('exact')}"
            if out.get('exact'):
                save_witness(D, out['N'], out['Q'], expr)
    else:
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=False)
        ent['consistent'] = (status != 'inconsistent')
        ent['rank'] = data.get('rank'); ent['nullity'] = data.get('nullity')
        ent['note'] = note + f' | {status}'
        if status == 'solved':
            save_witness(D, data['N'], data['Q'], expr)
    ent['seconds'] = round(time.time() - t0, 2)
    log(ent)
    return ent

def save_witness(D, N, Q, expr):
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    w = {'instance': IID, 'd': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in N],
         'Q': core.to_str(Q), 'verified': bool(ok), 'info': str(info)}
    with open(f'{RES}/witness-candidate.json', 'w') as f:
        json.dump(w, f, indent=1)
    print('!!! WITNESS SAVED verified=', ok, info, file=sys.stderr)

def P(s):
    return core.parse(s)

def run_modp(D, p, family, expr, note=''):
    """Consistency test mod p only (D may have mod-p coefficients). Logs and returns entry."""
    d = core.degree(D)
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, D, p, verbose=False)
    r = S.solve(verbose=False)
    ent = {'D': expr, 'd': d, 'family': family, 'e1': None, 'consistent': bool(r['consistent']),
           'rank': r['rank'], 'nullity': r['nullity'], 'rows': r['rows'], 'cols': r['cols'],
           'seconds': round(time.time() - t0, 2), 'note': note + f' | mod p={p} only', 'engine': 2, 'p': p}
    log(ent)
    return ent

def lin(coeffs):
    """linear form from dict var->coeff"""
    return {tuple(1 if j == i else 0 for j in range(6)): Fraction(c) for i, c in coeffs.items() if c}

# ---- sector diagnostics -------------------------------------------------
from flint import nmod_mat
SECTOR = [e for e in core.monomials(6, 1) if e[2] == 0 and e[5] == 0]   # core sextics, (Y0,Y1)-odd: 40 monomials
ODD6 = core.monomials(6, 1)                                               # 226 Y-odd sextic monomials
IDX = {e: i for i, e in enumerate(ODD6)}

def _rank_of_polys(polys, p):
    if not polys: return 0
    M = nmod_mat(len(polys), len(ODD6), [int(P_.get(e, 0)) % p for P_ in polys for e in ODD6], p)
    return M.rank()

def jacobian6(p):
    Fp = core.reduce_mod(F, p)
    J = []
    for i in range(6):
        Fi = core.pderiv(Fp, i, p)
        for e in core.monomials(1, None):
            Q = core.shift(Fi, e, 1, p)
            if core.parity_of(Q) == 1: J.append(Q)
    return J

def image_sector(D, p=None, verbose=False):
    """Targets served by D (mod p): total dim, dim in the core-odd sector, and modulo J."""
    p = p or core.PRIMES[0]
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D, p=p, verbose=verbose)
    J = jacobian6(p)
    rJ = _rank_of_polys(J, p)
    sector_polys = [{e: 1} for e in SECTOR]
    r_img = _rank_of_polys(basis, p)
    r_img_sec = _rank_of_polys(basis + sector_polys, p)
    dim_img_cap_sec = r_img + len(SECTOR) - r_img_sec
    # image ∩ sector modulo J: dim((img∩sec)+J) - dim J ; compute via rank of [basis ; sector] restricted... use:
    # dim(img ∩ sec mod J) = dim((img∩sec) + J) - dim J = (dim(img∩sec) - dim(img∩sec∩J)) ... simpler: rank(basis+J)+40-rank(basis+J+sector) - (rank(J)+40-rank(J+sector)) = dim((img+J)∩sec) - dim(J∩sec)
    r_imgJ = _rank_of_polys(basis + J, p)
    r_imgJ_sec = _rank_of_polys(basis + J + sector_polys, p)
    r_J_sec = _rank_of_polys(J + sector_polys, p)
    dim_imgJ_cap_sec = r_imgJ + len(SECTOR) - r_imgJ_sec
    dim_J_cap_sec = rJ + len(SECTOR) - r_J_sec
    Ap = core.reduce_mod(A, p)
    r_imgA = _rank_of_polys(basis + [Ap], p)
    out = {'d': core.degree(D), 'dim_image': r_img, 'dim_image_mod_J': r_imgJ - rJ, 'dim_J6': rJ,
           'dim_image_cap_sector': dim_img_cap_sec, 'dim_(image+J)_cap_sector': dim_imgJ_cap_sec,
           'dim_J_cap_sector': dim_J_cap_sec, 'sector_served_mod_J': dim_imgJ_cap_sec - dim_J_cap_sec,
           'A_in_image': r_imgA == r_img, 'seconds': round(time.time() - t0, 1), 'p': p}
    return out, basis
