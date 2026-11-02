"""Mechanism (i) test: a section-2.7-type exchange cover (any linear form, any base degree m)
forces the 15 quadrics q_j(Y) (coefficients of the X-quartic monomials in F) to span <=3 dims
(and <=2 dims for m>=3); dually for the (2,4)-part. Exact ranks over Q."""
import sys, json, itertools, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction
F,A,q,rec = core.load_instance('eight-plane-configuration')
quadpairs = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
def rank(rows):
    M=[r[:] for r in rows]; rk=0; ncol=len(M[0])
    for c in range(ncol):
        piv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]/M[rk][c]; M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
out={}
for name,(blockA,blockB) in {'F42':((0,1,2),(3,4,5)),'F24':((3,4,5),(0,1,2))}.items():
    # rows: monomials of degree 4 in blockA; columns: quadratic monomials in blockB
    rows=[]; gram_ranks=[]
    for e in core.monomials(4):
        if any(e[i] for i in blockB): continue
        row=[]
        for (i,j) in quadpairs:
            ee=list(e); ee[blockB[i]]+=1; ee[blockB[j]]+=1
            row.append(F.get(tuple(ee),Fraction(0)))
        rows.append(row)
        # 3x3 Gram matrix rank of this quadric
        G=[[Fraction(0)]*3 for _ in range(3)]
        for k,(i,j) in enumerate(quadpairs):
            if i==j: G[i][i]=row[k]
            else: G[i][j]=G[j][i]=row[k]/2
        gram_ranks.append(rank(G))
    out[name]={'n_quadrics':len(rows),'span_rank':rank(rows),'gram_ranks':gram_ranks}
print(json.dumps(out))
json.dump(out,open('/tmp/hodge/results/eight-plane-configuration/wave2-construct/mech_i_ranks.json','w'))
