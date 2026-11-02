"""Adjoint-condition count for V = X ∩ V(D), D|_X ∈ I_{Γ_Y}^{m1} (X-directions) and I_{Γ_X}^{m2} (Y-directions),
Γ_Y = X∩{X=0}, Γ_X = X∩{Y=0} smooth plane sextics (g=10, K=O(3)).  Adjoint orders a_i = m_i - 2.
Lower bound h^{2,0}(V~)^{sector} >= cond_sector - (S_d^sector - ker_sector), ker = ((X)^{a1}+(F)) ∩ ((Y)^{a2}+(F)) in S_d.
Realisability: D = (cF)_{X-deg >= m1}, c in S_{d-6}^even, needs (cF) to have no Y-even bidegree (k,l) with k<m1, l<m2.
"""
import sys; sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core
from flint import nmod_mat
p = core.PRIMES[0]
F = core.reduce_mod(qcp.F, p)

def h0(k):  # plane sextic curve, K = O(3)
    return {0:1, 1:3, 2:6, 3:10}.get(k, 6*k-9) if k < 4 else 6*k-9

def bideg(e): return (e[0]+e[1]+e[2], e[3]+e[4]+e[5])

def rank_span(polys, mons):
    if not polys: return 0
    col = {e:j for j,e in enumerate(mons)}
    flat = [0]*(len(polys)*len(mons))
    for r,P in enumerate(polys):
        for e,v in P.items(): flat[r*len(mons)+col[e]] = v % p
    return nmod_mat(len(polys), len(mons), flat, p).rank()

def count(d, a1, a2, parity):
    mons = [e for e in core.monomials(d, parity)]
    # conditions
    cond = sum((j+2)*(j+1)//2 * h0(d-j) for j in range(a1) if (d-j)%2 == parity)
    cond += sum((j+2)*(j+1)//2 * h0(d-j) for j in range(a2) if j%2 == parity)
    # kernel: U = (X)^{a1}_d + F*S_{d-6},  W = (Y)^{a2}_d + F*S_{d-6}  (parity sector)
    Fm = [core.shift(F, m, 1, p) for m in core.monomials(d-6, parity)] if d >= 6 else []
    U = [{e:1} for e in mons if bideg(e)[0] >= a1] + Fm
    W = [{e:1} for e in mons if bideg(e)[1] >= a2] + Fm
    rU, rW, rUW = rank_span(U, mons), rank_span(W, mons), rank_span(U+W, mons)
    ker = rU + rW - rUW
    img = len(mons) - ker
    return dict(d=d, a1=a1, a2=a2, parity=parity, S=len(mons), cond=cond, ker=ker, img=img, bound=cond-img)

def realisable(d, m1, m2):
    # does there exist Y-even c in S_{d-6} with (cF) having no Y-even bidegree (k,l), k<m1, l<m2 ?  (dimension of solution space)
    cm = core.monomials(d-6, 0)
    rows = []
    bad = lambda e: bideg(e)[0] < m1 and bideg(e)[1] < m2
    cols = [core.shift(F, m, 1, p) for m in cm]
    badmons = sorted(set(e for c in cols for e in c if bad(e)))
    if not badmons: return len(cm)
    col = {e:j for j,e in enumerate(badmons)}
    flat = [0]*(len(cm)*len(badmons))
    for r,c in enumerate(cols):
        for e,v in c.items():
            if e in col: flat[r*len(badmons)+col[e]] = v
    return len(cm) - nmod_mat(len(cm), len(badmons), flat, p).rank()

print("Fermat-type check on OUR F is not applicable; counts below are for V singular along both coordinate plane sextics.")
for d, m1, m2 in [(6,6,6),(6,6,4),(6,4,4),(8,6,4),(8,6,6),(8,4,4),(10,6,6),(10,8,4),(10,8,2),(10,6,4),(12,8,6),(12,10,4),(12,8,8),(12,10,2),(12,6,6)]:
    r = realisable(d, m1, m2)
    o = count(d, m1-2, m2-2, 1); e = count(d, m1-2, m2-2, 0)
    print(f"d={d} orders (m1,m2)=({m1},{m2}) realisable-dim(c)={r:4d} | ODD: cond={o['cond']:5d} ker={o['ker']:5d} img={o['img']:5d} bound={o['bound']:5d} | EVEN: cond={e['cond']:5d} ker={e['ker']:5d} img={e['img']:5d} bound={e['bound']:5d}", flush=True)
