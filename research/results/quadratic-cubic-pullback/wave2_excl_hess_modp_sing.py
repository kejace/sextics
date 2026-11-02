"""Candidate locus of a Hessian relaxation kernel in characteristic p via Singular (screening when the kernel is too
large to lift exactly).  usage: python3 wave2_excl_hess_modp_sing.py <KH json file> [timeout]"""
import sys, os, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import logline
from flint import nmod_mat
path = sys.argv[1]; TO = int(sys.argv[2]) if len(sys.argv) > 2 else 900
K = json.load(open(path))
p = K['p']; Wkeys = [tuple(x) for x in K['Wkeys']]; nw = len(K['wb']); nW = len(Wkeys)
basis = [{int(a): v % p for a, v in row.items()} for row in K['basis']]
rk = len(basis)
t0 = time.time()
Bm = nmod_mat(rk, nW, [basis[i].get(a, 0) for i in range(rk) for a in range(nW)], p)
Lam, nlam = Bm.nullspace()
quadrics = []
for k in range(nlam):
    terms = []
    for idx, (a, b) in enumerate(Wkeys):
        c = int(Lam[idx, k])
        if c:
            c = c if c <= p // 2 else c - p
            terms.append(f'{c}*w{a}*w{b}' if a != b else f'{c}*w{a}^2')
    quadrics.append('+'.join(terms).replace('+-', '-'))
sing = path.replace('.json', '.modp.sing')
with open(sing, 'w') as f:
    f.write(f'ring r = {p}, (' + ','.join(f'w{a}' for a in range(nw)) + '), dp;\n')
    f.write('ideal I = ' + ',\n'.join(quadrics) + ';\n')
    f.write('ideal G = std(I);\nprint("DIM " + string(dim(G)));\nprint("DEG " + string(mult(G)));\n')
    f.write('LIB "primdec.lib";\nlist L = minAssGTZ(I);\nprint("NCOMP " + string(size(L)));\n')
    f.write('int k; for (k = 1; k <= size(L); k++) { print("COMPONENT " + string(k) + " dim " + string(dim(std(L[k]))) + " ngens " + string(size(L[k]))); if (size(L[k]) < 80) { print(L[k]); } }\nquit;\n')
try:
    r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=TO)
    out = r.stdout
except subprocess.TimeoutExpired:
    out = 'TIMEOUT'
rec = dict(exp='hess-modp-singular', file=os.path.basename(path), d=K['d'], parity=K['parity'], prime=p, dim_KH=rk, nw=nw, n_quadrics=nlam,
           singular_tail=out[-2500:], seconds=round(time.time() - t0, 1), note='characteristic-p screening of the candidate locus {w : w w^T in K_H(p)}')
logline(os.path.join(w2x.LANE, 'wave2-excl.hessmodp.jsonl'), rec)
print(out[-2500:])
print(json.dumps({k: v for k, v in rec.items() if k != 'singular_tail'}))
