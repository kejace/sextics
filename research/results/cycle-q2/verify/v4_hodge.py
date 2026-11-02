"""Equivariant Hodge numbers by direct linear algebra (independent of the Koszul/character-formula method):
dim (R_k)_sector for k = 0,6,12,18,24 and the relevant sectors, mod p (good reduction => exact)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/verify')
import vcheck as V
import numpy as np
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 999983
F, A, q = V.load_instance()
F2 = {m: int(v) for m, v in V.scale(F, 2).items()}
G = V.diag_group(F)
# use a generating set to define sector keys (characters of G are determined by generators)
def gens_of(Gm):
    gens = []; S = {(0,)*6}
    for w in Gm:
        tw = tuple(int(x) for x in w)
        if tw in S: continue
        gens.append(tw)
        frontier = list(S)
        while frontier:
            nxt = []
            for s in frontier:
                t = tuple((a+b) % 6 for a, b in zip(s, tw))
                if t not in S: S.add(t); nxt.append(t)
            frontier = nxt
        if len(S) == len(Gm): break
    return np.array(gens, dtype=np.int64)
gens = gens_of(G)
def key(m): return tuple(int(x) for x in (np.array(m) @ gens.T) % 6)
eA = list(A.keys())[0]
sumw = lambda w: int(sum(w) % 6)
def sector_dims(k):
    """returns dict key -> (dim S_k^key, dim J_k^key)"""
    mk = V.monomials(k)
    bykey = {}
    for m in mk: bykey.setdefault(key(m), []).append(m)
    out = {}
    F2i = [V.pderiv(F2, i) for i in range(6)]
    if k >= 5:
        mk5 = V.monomials(k - 5)
        gcols = {}
        for i in range(6):
            kFi = key(list(F2i[i].keys())[0])
            for m in mk5:
                kk = tuple((a + b) % 6 for a, b in zip(key(m), kFi))
                gcols.setdefault(kk, []).append(V.mul(F2i[i], {m: 1}))
    for kk, ms in bykey.items():
        cols = gcols.get(kk, []) if k >= 5 else []
        if cols:
            rows = {m: j for j, m in enumerate(ms)}
            flat = [0] * (len(ms) * len(cols))
            for j, c in enumerate(cols):
                for m, v in c.items():
                    flat[rows[m] * len(cols) + j] = v % p
            rk = nmod_mat(len(ms), len(cols), flat, p).rank()
        else:
            rk = 0
        out[kk] = (len(ms), rk)
    return out
t0 = time.time()
res = {}
for k in (0, 6, 12, 18, 24):
    res[k] = sector_dims(k); print('degree', k, 'done', round(time.time() - t0, 1), 's', flush=True)
# sector characters: numerator weight key for the target sector = key(A) (chi_F = 0 for all w since F has all 6th powers)
kA = key(eA)
def dims_for(keyfun):
    return [res[k].get(keyfun(k), (0, 0))[0] - res[k].get(keyfun(k), (0, 0))[1] for k in (0, 6, 12, 18, 24)]
target = dims_for(lambda k: kA)
weight0 = dims_for(lambda k: (0,) * len(gens))
# residue-invariant sector: weight(P) + sum(w) == 0 for all w  -> key(P) = -sum(gens) componentwise
kinv = tuple((-sumw(g)) % 6 for g in gens)
inv = dims_for(lambda k: kinv)
total = [sum(v[0] - v[1] for v in res[k].values()) for k in (0, 6, 12, 18, 24)]
# K-sector (classes with trivial character on K = ker chi_A): residue character rho(P,w) = key(P).w + sum(w) ; need rho = 0 on K
# equivalently sum over G-characters chi with chi|_K trivial, i.e. chi in {rho_A^j : j=0,1,2}... compute directly by enumerating G:
Gp = np.array(sorted({tuple((w - w[2]) % 6) for w in G}))
chiA = lambda w: int(((np.array(eA) + 1) @ np.array(w)) % 6)
K = [w for w in Gp if chiA(w) == 0]
def rho_on(Pkeyvec, w):
    # need weight of P under w from key: reconstruct weight via gens? simpler: recompute from a representative monomial
    return None
# For the K-sector use representative monomials: for each key, take one monomial and compute rho on K
def dims_subgroup_invariant(sub):
    out = []
    for k in (0, 6, 12, 18, 24):
        tot = 0
        for kk, (ds, dj) in res[k].items():
            # representative monomial with this key
            rep = next(m for m in V.monomials(k) if key(m) == kk)
            if all((int(np.dot(rep, w)) + sumw(w)) % 6 == 0 for w in sub):
                tot += ds - dj
        out.append(tot)
    return out
Ksec = dims_subgroup_invariant(K)
# paired signs P and Gamma_book
P = [(0,)*6, (3,0,0,3,0,0), (0,3,0,0,3,0), (3,3,0,3,3,0)]
def span(gs):
    S = {(0,)*6}; fr = [(0,)*6]
    while fr:
        nx = []
        for s in fr:
            for g in gs:
                t = tuple((a+b) % 6 for a, b in zip(s, g))
                if t not in S: S.add(t); nx.append(t)
        fr = nx
    return [np.array(x) for x in S]
Gb = span([(3,0,0,3,0,0), (0,3,0,0,3,0), (2,1,0,3,2,1)])
Psec = dims_subgroup_invariant([np.array(x) for x in P]); Gbsec = dims_subgroup_invariant(Gb)
tau_sec = dims_subgroup_invariant([np.array((0,0,0,0,0,0)), np.array((0,0,0,3,3,3))])
summary = {'p': p, 'total': total, 'target_sector_weightA': target, 'weight0_sector': weight0, 'residue_invariant_sector': inv,
           'K_invariant': Ksec, 'pairedsigns_invariant': Psec, 'Gamma_book_invariant': Gbsec, 'tau_invariant': tau_sec,
           'R6_weightA': res[6].get(kA), 'R12_weightA': res[12].get(kA), 'seconds': round(time.time() - t0, 1)}
print(json.dumps(summary, default=str))
V.log(claim='equivariant Hodge numbers: total (1,426,1751,426,1); target sector (0,6,26,6,0); residue-invariant (0,6,26,6,0); weight-0 sector (1,6,30,8,0); K-sector (0,18,78,18,0); paired-signs (1,126,479,126,1); Gamma_book (0,22,74,22,0); tau (0,208,866,208,0)',
      source='dossier §1.6/A.2 (Koszul character formula)', verdict='CONFIRMED' if (total == [1,426,1751,426,1] and target == [0,6,26,6,0] and inv == [0,6,26,6,0] and Ksec == [0,18,78,18,0] and Psec == [1,126,479,126,1] and Gbsec == [0,22,74,22,0] and tau_sec == [0,208,866,208,0] and weight0 == [1,6,30,8,0]) else 'CHECK', evidence=summary, prime=p,
      notes='direct sector-wise rank of J_k mod p (independent of the character formula); exact because J has good reduction (smooth mod p)')
