"""Krull dimension and Hilbert series data of a Groebner basis from msolve output (-g 2 format), independent of Singular.
usage: python3 krull.py gb.ms nvars  -> prints DIM (affine Krull dimension of F_p[u]/I) from leading monomials (degrevlex)."""
import sys, re, itertools
txt = open(sys.argv[1]).read()
n = int(sys.argv[2])
# msolve -g 2 output: "#Reduced Groebner basis ... [p1, p2, ...]:"  polynomials separated by commas; variables u0..u{n-1}
body = txt[txt.index('['):txt.rindex(']')]
polys = [s.strip() for s in body.split(',\n')]
lead = []
for s in polys:
    # first term is the leading term in msolve's degrevlex output
    t = s.split('+')[0].split('-')[0] if not s.startswith('-') else s[1:].split('+')[0].split('-')[0]
    t = s.strip().lstrip('-').split('+')[0]
    # monomials like 3*u0^2*u5 ; strip coefficient
    exps = [0] * n
    for part in t.split('*'):
        m = re.match(r'u(\d+)(?:\^(\d+))?$', part.strip())
        if m: exps[int(m.group(1))] = int(m.group(2) or 1)
    lead.append(tuple(exps))
# Krull dimension of a monomial ideal = max size of a variable subset S such that no leading monomial is supported in S
# (greedy exact search via complement cover): dim = n - (min vertex cover of the hypergraph of supports)
supports = [frozenset(i for i, e in enumerate(L) if e) for L in lead]
supports = [s for s in supports if s]
# min hitting set by brute force on sizes (ok for small n)
best = None
for k in range(0, n + 1):
    for comb in itertools.combinations(range(n), k):
        c = set(comb)
        if all(s & c for s in supports):
            best = k; break
    if best is not None: break
print('NGB', len(polys), 'MINLEADDEG', min(sum(L) for L in lead), 'MAXLEADDEG', max(sum(L) for L in lead))
print('DIM', n - best)
