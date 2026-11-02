# Hodge numbers of character sectors of H^4(X)_prim via the Fermat-ring count (valid along the smooth
# equivariant family, as in dossier 1.4): residue of monomial Z^e (deg 6q, q=1,2,3 for (3,1),(2,2),(1,3))
# has characters  g: sum (e_j+1) w_j,  h: ..., sigma2: swaps e2<->e5 with sign -1 on the volume.
from itertools import product
wg=(0,2,0,1,3,3); wh=(0,0,1,0,0,5)
def chars(e):
    return (sum((e[j]+1)*wg[j] for j in range(6))%6, sum((e[j]+1)*wh[j] for j in range(6))%6)
box=[e for e in product(range(5),repeat=6)]
def count(pred):
    return tuple(sum(1 for e in box if sum(e)==6*q and pred(e)) for q in (1,2,3))
# tau-sector (Y-odd numerators): all residues of Y-odd A are tau-invariant
yodd=lambda e: (e[3]+e[4]+e[5])%2==1
print("tau-sector (Y-odd)            :", count(yodd))
print("h-invariant (any g), Y-odd    :", count(lambda e: yodd(e) and chars(e)[1]==0))
print("G-invariant                   :", count(lambda e: chars(e)==(0,0)))
# sigma2 acts on monomial residues: e -> swap(e2,e5), with sign -1 (volume). Self-paired monomials (e2==e5)
# are anti-invariant; pairs {e, swap e} split 1 inv + 1 anti.
def sig_split(pred):
    inv=[0,0,0]; anti=[0,0,0]
    for qi,q in enumerate((1,2,3)):
        S=[e for e in box if sum(e)==6*q and pred(e)]
        for e in S:
            if e[2]==e[5]: anti[qi]+=1
            else: inv[qi]+=0.5; anti[qi]+=0.5
    return tuple(inv), tuple(anti)
inv,anti=sig_split(lambda e: yodd(e) and chars(e)[1]==0)
print("h-inv Y-odd: sigma2-inv / anti:", inv, anti)
inv,anti=sig_split(lambda e: chars(e)==(0,0))
print("G-inv      : sigma2-inv / anti:", inv, anti)
# h-invariant sectors by g-character
for a in range(6):
    inv,anti=sig_split(lambda e: yodd(e) and chars(e)==(a,0))
    print(f"  (g={a},h=0) Y-odd: inv {inv} anti {anti}")
