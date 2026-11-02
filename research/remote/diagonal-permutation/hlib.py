"""Helper: run fixed-D tests on diagonal-permutation and append to log.jsonl."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction
IID = 'diagonal-permutation'
LOG = f'/tmp/hodge/results/{IID}/log.jsonl'
F, A, q, rec = core.load_instance(IID)

def log(entry):
    entry.setdefault('engine', 2)
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')

def test_D(Dexpr, family, note='', full=True, fast=True, verbose=False, D=None, want_exact=True):
    """E1 test then (optionally) full fixed-D solve.  Returns dict."""
    if D is None:
        D = core.parse(Dexpr)
    d = core.degree(D)
    t0 = time.time()
    entry = {'D': Dexpr, 'd': d, 'family': family, 'e1': None, 'consistent': None,
             'rank': -1, 'nullity': -1, 'seconds': 0.0, 'note': note}
    if core.parity_of(D) is None:
        entry['note'] += ' [mixed parity D, skipped]'
        log(entry); return entry
    if not core.is_homogeneous(D):
        entry['note'] += ' [non-homogeneous D, skipped]'
        log(entry); return entry
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    entry['e1'] = bool(ok)
    entry['e1_info'] = {k: int(v) for k, v in info.items()}
    if not ok or not full:
        entry['seconds'] = round(time.time() - t0, 2)
        log(entry); return entry
    if fast:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=want_exact)
        entry['consistent'] = bool(out['consistent']); entry['rank'] = int(out['rank']); entry['nullity'] = int(out['nullity'])
        entry['rows'] = int(out['rows']); entry['cols'] = int(out['cols'])
        if out['consistent'] and out.get('exact'):
            entry['note'] += ' SOLVED-EXACT'
            entry['witness'] = {'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])}
            with open(f'/tmp/hodge/results/{IID}/witness-candidate.json', 'w') as f:
                json.dump({'d': d, **entry['witness']}, f, indent=1)
        elif out['consistent']:
            entry['note'] += ' CONSISTENT-MODP (exact lift not done)'
    else:
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
        entry['consistent'] = (status != 'inconsistent')
        entry['rank'] = int(data.get('rank', -1)); entry['nullity'] = int(data.get('nullity', -1))
        entry['status'] = status
        if status == 'solved':
            entry['note'] += ' SOLVED-EXACT'
            entry['witness'] = {'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
            with open(f'/tmp/hodge/results/{IID}/witness-candidate.json', 'w') as f:
                json.dump({'d': d, **entry['witness']}, f, indent=1)
    entry['seconds'] = round(time.time() - t0, 2)
    log(entry)
    return entry

# ---------------------------------------------------------------------------
# compressed target_image: served target space (mod p) + consistency for A
# ---------------------------------------------------------------------------
import numpy as np
from flint import nmod_mat
def target_image_fast(D, p=None, k_extra=40, seed=1, verbose=False):
    """Returns (dimA, served_basis_modp, contains_A, info).  dimA = dim of Y-odd degree-6q targets served by D."""
    p = p or core.PRIMES[0]
    t0 = time.time()
    A0 = {core.monomials(6 * q, 1)[0]: 1}
    S = fixed_d.FixedDSystem(F, A0, q, D, p, with_target_space=True, verbose=verbose)
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    C = S.ncols
    k = min(S.nrows, C + k_extra)
    Cm = fixed_d._compress(M, p, k, seed)
    t1 = time.time()
    Mn = nmod_mat(k, C, Cm.ravel().tolist(), p)
    Xn, nullity = Mn.nullspace()
    t2 = time.time()
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    if nullity == 0:
        return 0, [], False, {'rows': S.nrows, 'cols': C, 'k': k}
    proj = nmod_mat(len(a_idx), nullity, [int(Xn[j, kk]) for j in a_idx for kk in range(nullity)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v: poly[S.cols[j][2]] = v
        basis.append(poly)
    Ap = core.reduce_mod(A, p)
    mons = {m: i for i, m in enumerate(core.monomials(6 * q, 1))}
    flat = [0] * ((rk + 1) * len(mons))
    for i, poly in enumerate(basis):
        for e, v in poly.items(): flat[i * len(mons) + mons[e]] = v
    for e, v in Ap.items(): flat[rk * len(mons) + mons[e]] = v
    rk2 = nmod_mat(rk + 1, len(mons), flat, p).rank()
    contains = (rk2 == rk)
    info = {'rows': S.nrows, 'cols': C, 'k': k, 't_build': round(t1 - t0, 1), 't_null': round(t2 - t1, 1)}
    if verbose:
        print(f'[timage p={p}] served dim={rk} containsA={contains} nullity={nullity} ({info})', file=sys.stderr)
    return rk, basis, contains, info

def test_D_image(Dexpr, family, note='', D=None):
    """E1 test, then compressed target-image test (consistency for A + served dimension).  Logs and returns entry."""
    if D is None:
        D = core.parse(Dexpr)
    d = core.degree(D)
    t0 = time.time()
    entry = {'D': Dexpr, 'd': d, 'family': family, 'e1': None, 'consistent': None, 'rank': -1, 'nullity': -1,
             'seconds': 0.0, 'note': note}
    if core.parity_of(D) is None or not core.is_homogeneous(D):
        entry['note'] += ' [bad D, skipped]'; log(entry); return entry
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    entry['e1'] = bool(ok)
    if not ok:
        entry['seconds'] = round(time.time() - t0, 2); log(entry); return entry
    dimA, basis, contains, inf = target_image_fast(D)
    entry['consistent'] = bool(contains); entry['served_dim'] = int(dimA); entry['rows'] = inf['rows']; entry['cols'] = inf['cols']
    entry['note'] += f' served_dim={dimA}'
    if contains:
        entry['note'] += ' CONSISTENT-MODP -> exact solve'
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True, want_exact=True)
        entry['rank'] = int(out['rank']); entry['nullity'] = int(out['nullity'])
        if out.get('exact'):
            entry['note'] += ' SOLVED-EXACT'
            entry['witness'] = {'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])}
            with open(f'/tmp/hodge/results/{IID}/witness-candidate.json', 'w') as f:
                json.dump({'d': d, **entry['witness']}, f, indent=1)
    entry['seconds'] = round(time.time() - t0, 2)
    log(entry)
    return entry
