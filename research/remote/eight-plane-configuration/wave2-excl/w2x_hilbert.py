"""Exact parity-refined Hilbert function of S and of the Jacobian ring R = S/J(F) for a smooth sextic
in 6 variables with the Y-parity grading (tau: Y -> -Y).  Koszul complex of the regular sequence
(F_X0,F_X1,F_X2 : degree 5, even ; F_Y0,F_Y1,F_Y2 : degree 5, odd).
Hilbert series in t (degree) and s (parity, s^2 = 1):
   H_S = 1/((1-t)^3 (1-st)^3),   H_R = H_S * (1-t^5)^3 (1-s t^5)^3,   H_J = H_S - H_R.
Also: check wave-1 e1kernel.jsonl rows: implied rank_p J_{d+6}^{par'} = kerdim + nN - nullity vs exact dim.
"""
import sys, json
from math import comb
sys.path.insert(0, '/home/kejace/devel/sextics-work')

def S_dim(k, par=None):
    """# monomials of degree k in X0,X1,X2,Y0,Y1,Y2 with Y-degree parity par (None = all)."""
    tot = 0
    for j in range(k + 1):           # j = Y-degree
        if par is not None and j % 2 != par: continue
        tot += comb(j + 2, 2) * comb(k - j + 2, 2)
    return tot

def hilbert_R(kmax=24):
    """dict (k, par) -> dim R_k^par via Koszul: R_k = sum_{a,b} (-1)^{a+b} C(3,a) C(3,b) S_{k-5(a+b)}^{par - b}."""
    out = {}
    for k in range(kmax + 1):
        for par in (0, 1):
            v = 0
            for a in range(4):
                for b in range(4):
                    kk = k - 5 * (a + b)
                    if kk < 0: continue
                    v += (-1) ** (a + b) * comb(3, a) * comb(3, b) * S_dim(kk, (par - b) % 2)
            out[(k, par)] = v
    return out

if __name__ == '__main__':
    HR = hilbert_R(26)
    print('k  S_k^+  S_k^-  | R_k^+ R_k^- R_k | J_k^+ J_k^-')
    for k in range(0, 27):
        Rp, Rm = HR[(k, 0)], HR[(k, 1)]
        print(f'{k:2d} {S_dim(k,0):6d} {S_dim(k,1):6d} | {Rp:5d} {Rm:5d} {Rp+Rm:5d} | {S_dim(k,0)-Rp:6d} {S_dim(k,1)-Rm:6d}')
    # total R_k check against dossier numbers
    tot = [HR[(k,0)] + HR[(k,1)] for k in range(25)]
    print('R_k totals:', tot)
    # wave-1 e1kernel.jsonl consistency: implied rank_p(J_{d+6}^{par'}) = kerdim + nN - nullity ; par' = (1 + parD) % 2 (A odd)
    print('\nwave-1 e1kernel rows: implied rank_p J_{d+6} vs exact')
    for fn in ['/tmp/hodge/results/eight-plane-configuration/e1kernel.jsonl', '/tmp/hodge/remote/eight-plane-configuration/e1kernel_remote.jsonl']:
        try:
            for line in open(fn):
                line = line.strip()
                if not line: continue
                r = json.loads(line)
                d, parD = r['d'], r['parD']
                nN = r['cols'] - r['nD']
                implied = r['kerdim'] + nN - r['nullity']
                parp = (1 + parD) % 2
                exact = S_dim(d + 6, parp) - HR[(d + 6, parp)]
                print(f"{fn.split('/')[-1]} d={d} parD={parD} p={r.get('p','2147483647')} implied rank_p J_{d+6}^{parp}={implied} exact={exact} {'OK' if implied==exact else 'MISMATCH'}  kermodJ={r['kermodJ']} exact dimJ_d={S_dim(d,parD)-HR[(d,parD)]} logged dimJ_d={r['dimJ_d']}")
        except FileNotFoundError:
            pass
