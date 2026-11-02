"""Own-code check of the attack's two-plane realisability table (adjcount.py): for D of degree d (Y-even) with
D|_X in (X)^m1 + (F) and in (Y)^m2 + (F): dimension of {c in S_{d-6}^even : cF has no monomial with Xdeg<m1 and Ydeg<m2}
plus dim ((X)^m1 cap (Y)^m2)_d. Also recomputes the 'bound = cond - img' numbers of the table (odd sector)."""
import sys, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
p = 1000000007
F, A, q = V.load_instance_F_A()
Fp = V.modp(F, p)
def bideg(e): return (e[0]+e[1]+e[2], e[3]+e[4]+e[5])
def realisable(d, m1, m2):
    cm = V.monomials(d-6, 0) if d >= 6 else []
    direct = sum(1 for e in V.monomials(d, 0) if bideg(e)[0] >= m1 and bideg(e)[1] >= m2)
    if not cm: return 0, direct
    cols = [V.shift(Fp, m, 1, p) for m in cm]
    bad = sorted({e for c in cols for e in c if bideg(e)[0] < m1 and bideg(e)[1] < m2})
    if not bad: return len(cm), direct
    r = V.rank_mod_p(cols, bad, p)
    return len(cm) - r, direct
def h0_curve(k):  # smooth plane sextic, K = O(3): h^0(O(k)) = 6k - 9 for k >= 4 (Riemann-Roch; deg = 6k, g = 10), small k by hand
    return {0: 1, 1: 3, 2: 6, 3: 10}.get(k, 6*k - 9) if k < 4 else 6*k - 9
def count(d, a1, a2, parity):
    mons = V.monomials(d, parity)
    cond = sum((j+2)*(j+1)//2 * h0_curve(d-j) for j in range(a1) if (d-j) % 2 == parity)
    cond += sum((j+2)*(j+1)//2 * h0_curve(d-j) for j in range(a2) if j % 2 == parity)
    Fm = [V.shift(Fp, m, 1, p) for m in V.monomials(d-6, parity)] if d >= 6 else []
    U = [{e: 1} for e in mons if bideg(e)[0] >= a1] + Fm
    W = [{e: 1} for e in mons if bideg(e)[1] >= a2] + Fm
    rU, rW, rUW = V.rank_mod_p(U, mons, p), V.rank_mod_p(W, mons, p), V.rank_mod_p(U + W, mons, p)
    ker = rU + rW - rUW
    img = len(mons) - ker
    return dict(S=len(mons), cond=cond, ker=ker, img=img, bound=cond - img)
rows = []
for d, m1, m2 in [(6,6,6),(6,6,4),(6,4,4),(8,6,4),(8,6,6),(8,4,4),(10,6,6),(10,8,8),(10,10,6),(10,12,4),(10,6,10),(10,4,12),(12,10,8),(12,8,10),(12,12,6),(12,6,12),(12,14,4),(12,9,9),(14,10,10),(14,12,8),(14,14,6)]:
    rc, direct = realisable(d, m1, m2)
    o = count(d, m1-2, m2-2, 1); e = count(d, m1-2, m2-2, 0)
    rows.append(dict(d=d, m1=m1, m2=m2, realisable_c=rc, direct=direct, odd=o, even=e))
    print(json.dumps(rows[-1]), flush=True)
json.dump(rows, open('/tmp/hodge/results/quadratic-cubic-pullback/verify/realisability_own.json', 'w'))
# verdicts
tab = {(r['d'], r['m1'], r['m2']): r for r in rows}
ok_table = (tab[(6,6,6)]['odd']['bound'] == 100 and tab[(6,6,6)]['realisable_c'] == 0 and tab[(8,6,4)]['realisable_c'] == 12
            and tab[(8,6,4)]['odd']['bound'] == 0 and tab[(10,6,6)]['realisable_c'] == 66 and tab[(10,6,6)]['odd']['bound'] == 0)
log('C24 attack §1.4 table rows (6,6,6),(8,6,4),(10,6,6): realisable c-dim and odd-sector adjoint bound (own code)',
    'CONFIRMED' if ok_table else 'REFUTED', json.dumps({k: (v['realisable_c'], v['odd']['bound']) for k, v in tab.items() if k in [(6,6,6),(8,6,4),(10,6,6)]}))
big = {k: v for k, v in tab.items() if k[1] + k[2] >= k[0] + 6 and k[0] >= 8}
all_zero = all(v['realisable_c'] == 0 and v['direct'] == 0 for v in big.values())
log('C25 attack §1.4 "Fermat two-plane mechanism not realisable in any degree": for d=8..14 every (m1,m2) with m1+m2>=d+6 tested has realisable dimension 0 (own code); beyond d=14 unverified',
    'CONFIRMED' if all_zero else 'REFUTED', json.dumps({str(k): (v['realisable_c'], v['direct'], v['odd']['bound']) for k, v in big.items()}))
