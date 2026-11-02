"""Generate msolve inputs for the Fano scheme of 2-planes in X=V(F): 20 affine charts of G(3,6),
plane = row space of 3x6 matrix with identity in columns (i,j,k) and unknowns m1..m9 elsewhere.
Equations: the 28 coefficients of F(s*row1+t*row2+r*row3) as a ternary sextic. Mod p=2147483629 (1 mod 4)."""
import sys, itertools, json, os
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
p=2147483629
F,A,q,rec=core.load_instance('eight-plane-configuration')
Fp={e:(c.numerator*pow(c.denominator,-1,p))%p for e,c in F.items()}
# polynomials in 12 vars: (s,t,r,m1..m9) as dict exp-tuple -> coeff
def pmul(a,b):
    c={}
    for ea,ca in a.items():
        for eb,cb in b.items():
            e=tuple(x+y for x,y in zip(ea,eb)); c[e]=(c.get(e,0)+ca*cb)%p
    return {e:v for e,v in c.items() if v}
def padd(a,b,scale=1):
    c=dict(a)
    for e,v in b.items(): c[e]=(c.get(e,0)+scale*v)%p
    return {e:v for e,v in c.items() if v}
outdir='/tmp/hodge/results/eight-plane-configuration/wave2-construct/fano'
os.makedirs(outdir,exist_ok=True)
charts=list(itertools.combinations(range(6),3))
for ci,chart in enumerate(charts):
    others=[c for c in range(6) if c not in chart]
    # linear forms Z_c in 12 vars
    Z=[]
    for c in range(6):
        lin={}
        if c in chart:
            k=chart.index(c); e=[0]*12; e[k]=1; lin[tuple(e)]=1
        else:
            kk=others.index(c)
            for row in range(3):
                e=[0]*12; e[row]=1; e[3+3*row+kk]=1; lin[tuple(e)]=1
        Z.append(lin)
    pows=[]
    for c in range(6):
        pw=[{tuple([0]*12):1}]
        for k in range(6): pw.append(pmul(pw[-1],Z[c]))
        pows.append(pw)
    tot={}
    for e,cf in Fp.items():
        term={tuple([0]*12):cf}
        for c in range(6):
            if e[c]: term=pmul(term,pows[c][e[c]])
        tot=padd(tot,term)
    eqs={}
    for e,v in tot.items():
        key=e[:3]; eqs.setdefault(key,{})[e[3:]]=v
    assert len(eqs)==28
    names=[f"m{i}" for i in range(1,10)]
    with open(f"{outdir}/chart{ci:02d}.ms","w") as f:
        f.write(",".join(names)+"\n"+str(p)+"\n")
        lines=[]
        for key,poly in eqs.items():
            terms=[]
            for e,v in poly.items():
                mon="*".join(f"{names[i]}^{e[i]}" if e[i]>1 else names[i] for i in range(9) if e[i])
                terms.append(f"{v}*{mon}" if mon else f"{v}")
            lines.append("+".join(terms))
        f.write(",\n".join(lines)+"\n")
    print(ci,chart,"nterms",sum(len(v) for v in eqs.values()),flush=True)
json.dump({'p':p,'charts':charts},open(f"{outdir}/charts.json","w"))
