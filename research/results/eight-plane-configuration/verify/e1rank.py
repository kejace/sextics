"""Independent E1-kernel computation by ranks (no nullspace/projection):
V_p(d,par) = {D in S_d^par : A*D in J_{d+6}} ;  dim V_p = nD - rank_p[M_A | M_J] + rank_p(M_J).
Also rank_p(M_J) is compared with the exact Koszul dimension dim_Q J_{d+6}^{par'} (no-drop check),
and rank_p(J_d^par) with dim_Q J_d^par.  Then dim_Q V <= dim V_p (rigorous when no-drop holds)."""
import sys, json, time; sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/verify')
from v_lib import *
from flint import nmod_mat
P_DEFAULT = 4294967291  # 2^32 - 5, not in the engine's prime list

def build(d, par, p):
    F, A, q, _ = load_instance()
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    parT = (1 + par) % 2            # parity of A*D (A odd)
    cols = []
    Dm = monomials(d, par)
    for m in Dm: cols.append(shift(Ap, m, 1, p))
    nD = len(cols)
    for i in range(6):
        pm = (parT - (0 if i < 3 else 1)) % 2   # m*F_i must have parity parT
        for m in monomials(d+1, pm): cols.append(shift(Fi[i], m, 1, p))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e]*C + j] = v
    M = nmod_mat(R, C, flat, p)
    # J_d^par generators
    jc = []
    for i in range(6):
        pm = (par - (0 if i < 3 else 1)) % 2
        for m in monomials(d-5, pm) if d >= 5 else []: jc.append(shift(Fi[i], m, 1, p))
    return M, nD, jc, R, C

def run(d, par, p=P_DEFAULT):
    t0 = time.time()
    M, nD, jc, R, C = build(d, par, p)
    t1 = time.time()
    rank_full = M.rank()
    t2 = time.time()
    # rank of M_J = columns nD..C-1
    MJ = nmod_mat(R, C-nD, [int(M[r, c]) for r in range(R) for c in range(nD, C)], p)
    rank_J = MJ.rank()
    t3 = time.time()
    parT = (1 + par) % 2
    exact_J_target = dimJ_exact(d+6, parT)
    dimV = nD - rank_full + rank_J
    rank_Jd = 0
    if jc:
        rws = {}
        for c in jc:
            for e in c: rws.setdefault(e, len(rws))
        fl = [0]*(len(rws)*len(jc))
        for j, c in enumerate(jc):
            for e, v in c.items(): fl[rws[e]*len(jc)+j] = v
        rank_Jd = nmod_mat(len(rws), len(jc), fl, p).rank()
    exact_Jd = dimJ_exact(d, par) if d >= 5 else 0
    forced = max(0, dimR_exact(d, par) - dimR_exact(d+6, parT))
    out = {'d': d, 'parD': par, 'p': p, 'rows': R, 'cols': C, 'nD': nD, 'rank_full': rank_full, 'rank_J_target_p': rank_J,
           'exact_J_target': exact_J_target, 'nodrop_target': rank_J == exact_J_target, 'dimV_p': dimV,
           'rank_Jd_p': rank_Jd, 'exact_Jd': exact_Jd, 'nodrop_Jd': rank_Jd == exact_Jd, 'kermodJ_p': dimV - rank_Jd,
           'forced_lower_bound': forced, 'AnnR_A_d_exact_Q': (dimV - exact_Jd) if (rank_J == exact_J_target and (dimV - exact_Jd) == forced) else None,
           't_build': round(t1-t0,1), 't_rank_full': round(t2-t1,1), 't_rank_J': round(t3-t2,1), 'method': 'own code e1rank.py (ranks)'}
    return out

if __name__ == '__main__':
    p = P_DEFAULT
    args = sys.argv[1:]
    if args and args[0].startswith('p='): p = int(args[0][2:]); args = args[1:]
    for a in args:
        d, par = (int(a.split(':')[0]), int(a.split(':')[1])) if ':' in a else (int(a), None)
        for pr in ([par] if par is not None else [0, 1]):
            out = run(d, pr, p)
            print(json.dumps(out), flush=True)
            log({'check': f'e1kernel d={d} par={pr}', **out})
            with open('/tmp/hodge/results/eight-plane-configuration/verify/e1rank.jsonl', 'a') as f: f.write(json.dumps(out)+'\n')
