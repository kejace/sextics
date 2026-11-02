"""Step 4a: structured fixed candidates (principal Hessian minors of F; degree 8 and 12), written for hodge.sweep."""
import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import w2s_lib as L
from hodge import core
from fractions import Fraction
I = L.Inst()
F = I.F
H = [[core.pderiv(I.Fi[a], b) for b in range(6)] for a in range(6)]
names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']

def det(M):
    n = len(M)
    if n == 1:
        return dict(M[0][0])
    out = {}
    for j in range(n):
        if not M[0][j]:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = core.mul(M[0][j], det(minor))
        out = core.add(out, core.scale(term, (-1) ** j))
    return out

out = open('/tmp/hodge/results/cycle-q2/w2s_structured_cands.jsonl', 'w')
n = 0
for size in (2, 3):
    for S in itertools.combinations(range(6), size):
        M = [[H[a][b] for b in S] for a in S]
        Dm = det(M)
        Dm = {e: v for e, v in Dm.items() if v != 0}
        if not Dm:
            print('zero minor', S); continue
        par = core.parity_of(Dm)
        keys = {I.key(e) for e in Dm}
        rec = {'D': core.to_str(Dm), 'family': f'hess-minor-{"".join(names[a] for a in S)}', 'parity': par, 'nkeys': len(keys), 'd': core.degree(Dm)}
        if par == 0 and len(keys) == 1:
            out.write(json.dumps({'D': rec['D'], 'family': rec['family']}) + '\n'); n += 1
        else:
            print('skip', rec['family'], 'parity', par, 'keys', len(keys))
out.close()
print('wrote', n, 'candidates')
