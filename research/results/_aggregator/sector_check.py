# Molien-Koszul count of the isotypic Hodge numbers of symmetric-sextic-residual (independent check of the verify lens)
# dim R_n^{chi} = (1/|G|) sum_g conj(chi(g)) [t^n] det(1 - t^5 g)/det(1 - t g),  G = S3 x S3 x <tau>, acting on (X0,X1,X2,Y0,Y1,Y2)
import itertools, numpy as np
from fractions import Fraction
N = 19
def perm_matrix(p):
    M = np.zeros((3,3), dtype=int)
    for i,j in enumerate(p): M[j,i] = 1
    return M
def sgn(p):
    s = 1
    for i in range(3):
        for j in range(i+1,3):
            if p[i] > p[j]: s = -s
    return s
def series_det_ratio(g):
    # coefficients of prod_i (1 - t^5 lam_i)/(1 - t lam_i) up to t^N, via eigenvalues (complex)
    lam = np.linalg.eigvals(g.astype(complex))
    num = np.zeros(N+1, dtype=complex); num[0] = 1
    for l in lam:
        new = num.copy()
        new[5:] -= l**5 * num[:-5] if N >= 5 else 0
        num = new
    den = np.zeros(N+1, dtype=complex); den[0] = 1
    for l in lam:
        # multiply by 1/(1 - l t) = sum l^k t^k
        new = np.zeros(N+1, dtype=complex)
        for k in range(N+1):
            new[k] = sum(den[j] * l**(k-j) for j in range(k+1))
        den = new
    return den * 0 + np.convolve(num, den)[:N+1]
perms = list(itertools.permutations(range(3)))
rows = {}
chars = {'triv,triv': lambda a,b: 1, 'sgn,sgn': lambda a,b: sgn(a)*sgn(b), 'sgn,triv': lambda a,b: sgn(a), 'triv,sgn': lambda a,b: sgn(b)}
for cname, chi in chars.items():
    for ypar, label in ((1,'Y-even'),(-1,'Y-odd')):
        tot = np.zeros(N+1, dtype=complex)
        for a in perms:
            for b in perms:
                for s in (1,-1):
                    g = np.zeros((6,6), dtype=int)
                    g[:3,:3] = perm_matrix(a); g[3:,3:] = s*perm_matrix(b)
                    psi = chi(a,b) * (1 if ypar == 1 else s)   # tau-character of the sector
                    tot += psi * series_det_ratio(g)
        dims = np.round(tot.real/72).astype(int)
        rows[(cname,label)] = dims
        print(f'{cname:10s} {label:7s} dims n=0..18: {list(dims)}   (3,1),(2,2),(1,3) = {dims[6]},{dims[12]},{dims[18]}')
