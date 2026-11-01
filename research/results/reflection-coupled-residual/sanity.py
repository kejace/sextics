import sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction
IID='reflection-coupled-residual'
F, A, q, rec = core.load_instance(IID)
print('q=',q,'A=',core.to_str(A)); print('F terms',len(F)); print(core.to_str(F))
# sanity: D = 1
D1 = core.parse('1')
t0=time.time(); ok, info = fixed_d.e1_solvable(F, A, q, D1); print('E1 with D=1 (A in J(F)?):', ok, info, f'{time.time()-t0:.1f}s')
t0=time.time(); st, data = fixed_d.solve_fixed_D(F, A, q, D1); print('solve D=1:', st, {k:v for k,v in data.items() if k not in ('N','Q')}, f'{time.time()-t0:.1f}s')
# diagonal symmetry group: find all (a0..a5) in Z/6 (scaling by 6th roots... general: solve for weights w in Z/n with <w,e> const for all monomials)
# Projective diagonal symmetries: diag(z^w0..z^w5) with z primitive n-th root; F(z^w Z) = z^k F(Z) for all monomials e in F: <w,e> = k mod n.
exps = list(F.keys())
def diag_syms(n):
    out=[]
    for w in itertools.product(range(n), repeat=6):
        vals = {sum(wi*ei for wi,ei in zip(w,e)) % n for e in exps}
        if len(vals)==1:
            out.append((w, vals.pop()))
    return out
for n in [2,3,4,6,12]:
    S = diag_syms(n)
    # quotient by scalar (1,1,1,1,1,1)*c
    print(f'n={n}: {len(S)} diagonal symmetries (incl. scalars) ')
    if n in (2,6,12):
        for w,k in S:
            # character of residue of A*Omega/F^2: A has weight <w,eA>, Omega weight sum(w), F^2 weight 2k
            eA = next(iter(A))
            chi = (sum(wi*ei for wi,ei in zip(w,eA)) + sum(w) - 2*k) % n
            print('  w=',w,'F-weight',k,'residue char exponent',chi)
