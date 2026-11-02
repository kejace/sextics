"""Own-code E1 kernel sweep: K_d^par = {D in S_d^par : A D in J_{d+6}}, mod p (default 1000000007).
Certifies exactness: rank of the J-block must equal dim J_{d+6}^{par'} = S - R from the parity Hilbert function
(then dim K computed mod p is an upper bound for dim K over Q, and J_d gives the lower bound).
usage: python3 e1_kernel_own.py dmin dmax [p]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
p = int(sys.argv[3]) if len(sys.argv) > 3 else 1000000007
F, A, q = V.load_instance_F_A()
Fp, Ap = V.modp(F, p), V.modp(A, p)
Fi = [V.pderiv(Fp, i, p) for i in range(6)]
HR = V.hilbert_R_parity()
for d in range(dmin, dmax + 1):
    for par in (0, 1):
        t0 = time.time()
        tpar = (par + 1) % 2                          # parity of A*D
        cols = V.monomials(d + 6, tpar)
        Jrows = []
        for i in range(6):
            fpar = 0 if i < 3 else 1                  # parity of F_i
            for m in V.monomials(d + 1, (tpar - fpar) % 2):
                Jrows.append(V.shift(Fi[i], m, 1, p))
        Dm = V.monomials(d, par)
        Arows = [V.shift(Ap, m, 1, p) for m in Dm]
        MJ = V.to_matrix_rows(Jrows, cols, p)
        rJ = MJ.rank()
        M = V.to_matrix_rows(Jrows + Arows, cols, p)
        rM = M.rank()
        dimK = len(Dm) - (rM - rJ)
        # dim J_d^par = sum_i #monomials(d-5, par - fpar)
        dimJd = sum(len(V.monomials(d - 5, (par - (0 if i < 3 else 1)) % 2)) for i in range(6)) if d >= 5 else 0
        expected_rJ = len(cols) - HR[d + 6][tpar]
        rec = dict(d=d, parity=par, p=p, cols=len(cols), rowsJ=len(Jrows), rankJ=rJ, rankJ_expected=expected_rJ,
                   rank_stack=rM, nD=len(Dm), dimK=dimK, dimJd=dimJd, dimK_mod_J=dimK - dimJd,
                   ringdim_target=len(cols) - rJ, ringdim_expected=HR[d + 6][tpar], seconds=round(time.time() - t0, 1))
        exact_ok = (rJ == expected_rJ)
        verdict = 'CONFIRMED' if exact_ok else 'WEAKENED'
        log(f'C13 E1 kernel d={d} parity={par}: dim K_d/J_d (own code, p={p})', verdict,
            f"dimK={dimK} dimJd={dimJd} dimK_mod_J={dimK - dimJd}; rankJ={rJ} (Hilbert expects {expected_rJ}, exactness certified: {exact_ok}); {rec['seconds']}s",
            **rec)
        with open('/tmp/hodge/results/quadratic-cubic-pullback/verify/e1_kernel_own.jsonl', 'a') as f:
            f.write(json.dumps(rec) + '\n')
