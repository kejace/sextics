"""Served-space profiling: for structured denominators D, compute the whole space of Y-odd degree-6 targets served by D
(fixed_d.FixedDSystem with with_target_space=True, compressed rows, nullspace mod p), reduce it modulo J_6 and split it by
Klein residue sector.  A denominator whose served space meets the target's sector (A = X0^3 Y0^3, Klein (0,1)) is a
building block even if it does not serve A itself (served(D1 D2) = served(D1) + served(D2)).
Usage: python3 job_served.py <candidate file: tag<TAB>expr per line> [max_degree]"""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from flint import nmod_mat
from hodge import core, fixed_d, jacobian
from hodge.core import monomials
import rcr2
from rcr2 import chi, sec_str, log_event, P
F, A, q, rec = core.load_instance('reflection-coupled-residual')
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 8
# reducer for degree-6 Y-odd targets modulo J_6
red_std, red = jacobian.jacobian_ring_basis(F, 6, P, parity=1)
std_idx = {m: i for i, m in enumerate(red.standard)}
def profile(D, tag):
    t0 = time.time()
    d = core.degree(D)
    A0 = {monomials(6, 1)[0]: 1}
    S = fixed_d.FixedDSystem(F, A0, q, D, P, with_target_space=True)
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    k = min(S.nrows, S.ncols + 60)
    Cm = fixed_d._compress(M, P, k, seed=3)
    Mn = nmod_mat(k, S.ncols, Cm.ravel().tolist(), P)
    X, nullity = Mn.nullspace()
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    out = {'event': 'served', 'tag': tag, 'd': d, 'rows': S.nrows, 'cols': S.ncols, 'k': k, 'nullity': nullity}
    if nullity == 0:
        out['served_raw'] = 0; out['served_modJ'] = 0; out['by_sector'] = {}
    else:
        proj = nmod_mat(len(a_idx), nullity, [int(X[j, kk]) for j in a_idx for kk in range(nullity)], P)
        T = proj.transpose()
        Tr, rk = T.rref()
        polys = []
        for i in range(rk):
            poly = {}
            for jj, j in enumerate(a_idx):
                v = int(Tr[i, jj])
                if v:
                    poly[S.cols[j][2]] = v
            polys.append(poly)
        out['served_raw'] = rk
        # reduce mod J_6 and split by sector
        reduced = [red.normal_form(pp) for pp in polys]
        bysec = {}
        for pp in reduced:
            for e, v in pp.items():
                bysec.setdefault(sec_str(chi(e)), []).append(None)
        # rank per sector of the reduced vectors
        sec_rank = {}
        for s in ('00', '01', '10', '11'):
            vecs = []
            for pp in reduced:
                row = [0] * len(red.standard)
                nz = False
                for e, v in pp.items():
                    if sec_str(chi(e)) == s:
                        row[std_idx[e]] = v; nz = True
                if nz:
                    vecs.append(row)
            if vecs:
                sec_rank[s] = nmod_mat(len(vecs), len(red.standard), [x for r in vecs for x in r], P).rank()
            else:
                sec_rank[s] = 0
        allrows = [[0] * len(red.standard) for _ in reduced]
        for r, pp in enumerate(reduced):
            for e, v in pp.items():
                allrows[r][std_idx[e]] = v
        out['served_modJ'] = nmod_mat(len(reduced), len(red.standard), [x for r in allrows for x in r], P).rank() if reduced else 0
        out['by_sector'] = sec_rank
        # does the served space contain A's class?  (A reduced mod J, membership)
        Ared = red.normal_form(core.reduce_mod(A, P))
        arow = [0] * len(red.standard)
        for e, v in Ared.items():
            arow[std_idx[e]] = v
        out['A_served'] = bool(reduced) and nmod_mat(len(reduced) + 1, len(red.standard), [x for r in allrows for x in r] + arow, P).rank() == out['served_modJ']
    out['seconds'] = round(time.time() - t0, 1)
    print('[served]', out, flush=True)
    log_event(out)
    return out
for line in open(sys.argv[1]):
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    tag, expr = line.split('\t', 1) if '\t' in line else ('cand', line)
    try:
        D = core.parse(expr)
    except Exception as ex:
        print('parse error', tag, ex); continue
    if not D or not core.is_homogeneous(D) or core.parity_of(D) is None or core.degree(D) > dmax:
        continue
    try:
        profile(D, tag)
    except Exception as ex:
        log_event({'event': 'served', 'tag': tag, 'error': repr(ex)[:300]})
        print('error', tag, ex, flush=True)
print('SERVED_DONE')
