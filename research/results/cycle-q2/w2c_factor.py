#!/usr/bin/env python3
"""(iv) classify every factorization A = M'*B' with M' Y-even (the q-raising families F_s = F + s M').
Columns: sibling residue character of (F, q=1, B') on G (exponent values), order of the diagonal symmetry
group of F + sM' (w in G with chi(M') = chi_F), whether M' is paired-sign invariant (cubic descent survives),
and the set of v-degree residues mod 6 of invariant terms (torus-weight bookkeeping) that are linear in s:
the elliptic surface over the v-line stays the RES (polynomial in v^6) iff the set is {0}."""
import itertools, json
A=(2,2,3,2,2,1)
Fm=[(6,0,0,0,0,0),(0,6,0,0,0,0),(0,0,6,0,0,0),(0,0,0,6,0,0),(0,0,0,0,6,0),(0,0,0,0,0,6),(3,1,0,1,1,0),(0,3,1,0,1,1),(1,0,3,1,0,1)]
def wt(m,w): return sum(e*x for e,x in zip(m,w))%6
G=[]
for w in itertools.product(range(6),repeat=5):
    w=w+(0,)
    if len({wt(m,w) for m in Fm})==1: G.append(w)
chiF=lambda w: wt(Fm[0],w)
def res_char(P,q,w): return (wt(P,w)+sum(w)-(q+1)*chiF(w))%6
P=[(3,0,0,3,0,0),(0,3,0,0,3,0),(0,0,3,0,0,3)]
rows=[]
for Mp in itertools.product(*[range(a+1) for a in A]):
    if sum(Mp[3:])%2: continue
    if sum(Mp)!=6: continue   # F + sM' homogeneous: deg M' = 6
    Bp=tuple(a-m for a,m in zip(A,Mp))
    chiB={res_char(Bp,1,w) for w in G}
    sym=[w for w in G if wt(Mp,w)==chiF(w)]
    Pinv=all(wt(Mp,p)==chiF(p) for p in P)
    symfib=[w for w in sym if (w[3]-w[0])%6==(w[4]-w[1])%6==(w[5]-w[2])%6]
    chiB_fib=sorted({res_char(Bp,1,w) for w in symfib})
    symK=all(res_char(A,2,w)==0 for w in sym)
    # torus bookkeeping in the cubic chart: u-exponent m=(a+b)/2 (needs a_i+b_i even), v-degree |b|
    vres=None
    if Pinv:
        m=tuple((Mp[i]+Mp[3+i])//2 for i in range(3)); vb=sum(Mp[3:])
        coup=[(2,1,0),(0,2,1),(1,0,2)]
        vres=set()
        for nd,ne,nf in itertools.product(range(3),repeat=3):
            for npr in (1,2):
                tot=[nd*coup[0][i]+ne*coup[1][i]+nf*coup[2][i]+npr*m[i] for i in range(3)]
                if len({t%3 for t in tot})==1:
                    vres.add((2*(nd+ne+nf)+npr*vb)%6)
    rows.append({"M":Mp,"B":Bp,"deg_M":sum(Mp),"Ydeg_M":sum(Mp[3:]),"sibling_char":sorted(chiB),"sym_order":len(sym),"P_invariant":Pinv,
                 "v_residues":sorted(vres) if vres is not None else None,"symfib_order":len(symfib),"sibling_char_on_symfib":chiB_fib,"sym_in_K":symK})
print("factorizations deg M'=6, Y-even:",len(rows))
triv=[r for r in rows if r["sibling_char"]==[0]]
print("sibling with trivial G-character:",len(triv))
for r in triv: print("  M'",r["M"],"B'",r["B"],"sym",r["sym_order"],"P",r["P_invariant"],"vres",r["v_residues"])
full=[r for r in rows if r["sym_order"]==72]
print("families keeping all of G (M' G-invariant):",len(full))
for r in full: print("  M'",r["M"],"B'",r["B"],"char",r["sibling_char"],"P",r["P_invariant"],"vres",r["v_residues"])
res=[r for r in rows if r["P_invariant"] and r["v_residues"]=={0} or (r["P_invariant"] and r["v_residues"]==[0])]
print("P-invariant M' with RES surviving (v-residues {0}):",len(res))
for r in res: print("  M'",r["M"],"B'",r["B"],"char",r["sibling_char"],"sym",r["sym_order"])
json.dump(rows,open("w2c_factor.json","w"))

print("--- RES-surviving (P-invariant, v-residues {0} or empty) ---")
for r in rows:
    if r["P_invariant"] and (r["v_residues"] in ([0],[])):
        print("  M'",r["M"],"B'",r["B"],"sym",r["sym_order"],"symfib",r["symfib_order"],"char(B') on G",r["sibling_char"],"on G'_fib",r["sibling_char_on_symfib"],"G'<=K",r["sym_in_K"],"vres",r["v_residues"])
print("--- all P-invariant with G' <= K (sibling trivial on its own symmetry group) ---")
for r in rows:
    if r["P_invariant"] and r["sym_in_K"]:
        print("  M'",r["M"],"B'",r["B"],"sym",r["sym_order"],"vres",r["v_residues"])
