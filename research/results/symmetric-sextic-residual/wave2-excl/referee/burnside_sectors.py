# Independent recount: S3xS3-equivariant Hilbert function of the Jacobian ring in the Fermat monomial model
# (basis: monomials with all exponents <= 4; exact for every smooth S3xS3-invariant F by the equivariant Koszul argument).
# Polynomial character multiplicities by Burnside; the geometric (residue) character is chi_poly * sgn(piX)sgn(piY).
from itertools import permutations, product
def sgn(p):
    s=1
    for i in range(3):
        for j in range(i+1,3):
            if p[i]>p[j]: s=-s
    return s
G=[(a,b) for a in permutations(range(3)) for b in permutations(range(3))]
def fixed_count(n, ypar, g):
    a,b=g; cnt=0
    for e in product(range(5),repeat=6):
        if sum(e)!=n or (e[3]+e[4]+e[5])%2!=ypar: continue
        f=(e[a[0]],e[a[1]],e[a[2]],e[3+b[0]],e[3+b[1]],e[3+b[2]])
        if f==e: cnt+=1
    return cnt
chars={'triv,triv':lambda g:1,'sgn,sgn':lambda g:sgn(g[0])*sgn(g[1]),'sgn,triv':lambda g:sgn(g[0]),'triv,sgn':lambda g:sgn(g[1])}
for ypar in (1,0):
    for name,chi in chars.items():
        dims=[]
        for n in (6,12,18):
            tot=sum(chi(g)*fixed_count(n,ypar,g) for g in G)
            assert tot%36==0
            dims.append(tot//36)
        geo={'triv,triv':'sgn,sgn','sgn,sgn':'triv,triv','sgn,triv':'triv,sgn','triv,sgn':'sgn,triv'}[name]
        print(f"Y-{'odd' if ypar else 'even'} polynomial char {name:10s} -> residue char {geo:10s}: (R6,R12,R18) = {tuple(dims)}")
