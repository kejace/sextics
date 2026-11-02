# dim of the polynomially S3xS3-invariant, Y-odd part of the Jacobian ring in degrees 14 and 16 (Fermat model, exact by the equivariant Koszul argument)
from itertools import permutations, product
G=[(a,b) for a in permutations(range(3)) for b in permutations(range(3))]
def dim_inv(n, ypar):
    tot=0
    for e in product(range(5),repeat=6):
        if sum(e)!=n or (e[3]+e[4]+e[5])%2!=ypar: continue
        for a,b in G:
            f=(e[a[0]],e[a[1]],e[a[2]],e[3+b[0]],e[3+b[1]],e[3+b[2]])
            if f==e: tot+=1
    assert tot%36==0; return tot//36
for n in (14,16,18):
    print(f"degree {n}: invariant Y-odd dim R = {dim_inv(n,1)}, invariant Y-even = {dim_inv(n,0)}")
