# Do the pickled Farkas denominators equal the named polynomials? (own constructions; engine used only to load F)
import pickle, sys
from fractions import Fraction
P='/tmp/hodge/results/symmetric-sextic-residual/wave2-excl/'
def var(i): e=[0]*6; e[i]=1; return {tuple(e):Fraction(1)}
def mul(A,B):
    o={}
    for e,c in A.items():
        for f,d in B.items():
            g=tuple(a+b for a,b in zip(e,f)); o[g]=o.get(g,0)+c*d
    return {e:c for e,c in o.items() if c}
def add(A,B,s=1):
    o=dict(A)
    for e,c in B.items(): o[e]=o.get(e,0)+s*c
    return {e:c for e,c in o.items() if c}
X=[var(i) for i in range(3)]; Y=[var(3+i) for i in range(3)]
def delta(V): return mul(mul(add(V[0],V[1],-1),add(V[0],V[2],-1)),add(V[1],V[2],-1))
def e2(V): return add(add(mul(V[0],V[1]),mul(V[0],V[2])),mul(V[1],V[2]))
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
F,A,q,_=core.load_instance('symmetric-sextic-residual')
Fq={tuple(int(x) for x in e):Fraction(c) for e,c in F.items()}
def pd(Pp,i):
    o={}
    for e,c in Pp.items():
        if e[i]:
            f=list(e); f[i]-=1; o[tuple(f)]=o.get(tuple(f),0)+c*e[i]
    return {e:c for e,c in o.items() if c}
def euler(V, off):
    out={}
    for i in range(3): out=add(out,mul(V[i],pd(Fq,off+i)))
    return out
EXF=euler(X,0); EYF=euler(Y,3)
named={'farkas_D7.pkl':mul(delta(X),e2([mul(x,x) for x in X])),'farkas_D7Y.pkl':mul(delta(Y),e2([mul(y,y) for y in Y])),
       'farkas_DXxDY.pkl':mul(delta(X),delta(Y)),'farkas_EXF.pkl':EXF,'farkas_EYF.pkl':EYF,'farkas_EXFmEYF.pkl':add(EXF,EYF,-1)}
def normalize(Pp):
    k=next(iter(sorted(Pp))); c=Pp[k]; return {e:v/c for e,v in Pp.items()}
for name,mine in named.items():
    rec=pickle.load(open(P+name,'rb'))
    D={tuple(int(x) for x in k.strip('()').split(',')):Fraction(v) for k,v in rec['D'].items()}
    print(f"{name}: label={rec['label']} chi={rec['chi']} terms={len(D)} equals own construction up to scalar: {normalize(D)==normalize(mine)}")
