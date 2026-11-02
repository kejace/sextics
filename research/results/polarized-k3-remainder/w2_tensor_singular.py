"""Export the exact (mod p) tensor-frame Hessian relaxation at Jacobian-type degree D as quadrics for Singular.
usage: python3 w2_tensor_singular.py D P
Variables: cofactor coordinates c_(i,m) minus the 21 Euler coordinates (0, U_0*h) (set to 0).  J8 = ideal of the
pure-U component P_8 (all cofactor coordinates of positive V-degree), so that  sat(I, I + J8)  detects components
other than the P_8 one.  Output: w2-tensor-relax-d{D}-{P}.sing / .json"""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np, w2tensor, w2lib as L
from hodge import core
d = int(sys.argv[1]); p = int(sys.argv[2])
t0 = time.time()
TR = w2tensor.TensorRelaxation(d, verbose=True, modp=p)
K = TR.K
pos = {(i, m): k for k, (i, m) in enumerate(TR.basis)}
drop = set()
for h in TR.quadrics_uv:
    m = list(h); m[0] += 1; drop.add(pos[(0, tuple(m))])
keep = [k for k in range(K) if k not in drop]; kpos = {k: i for i, k in enumerate(keep)}
nk = len(keep)
Wm = TR.Wmod                      # pairs x tensor-dim (mod p), coordinate t of pair (k,l) = coefficient of c_k c_l (k<l) or c_k^2
pairidx = {pr: i for i, pr in enumerate(TR.pairs)}
coords = TR.coord_rows            # independent antisymmetric coordinates
vdeg = [m[3] + m[4] + m[5] for (_, m) in TR.basis]
lines = []
for t in coords:
    terms = []
    for i, (k, l) in enumerate(TR.pairs):
        c = int(Wm[i, t])
        if c and k in kpos and l in kpos:
            if k != l: c = (2 * c) % p
            if c > p // 2: c -= p
            terms.append(f'({c})*x({kpos[k]+1})*x({kpos[l]+1})')
    if terms: lines.append('+'.join(terms))
name = f'w2-tensor-relax-d{d}-{p}'
with open(name + '.sing', 'w') as fh:
    fh.write('LIB "elim.lib";\n')
    fh.write(f'ring r = {p}, (x(1..{nk})), dp;\noption(redSB);\nideal I;\n')
    for s, q in enumerate(lines):
        fh.write(f'I[{s+1}] = {q};\n')
    fh.write('ideal J8;\n')
    for i, k in enumerate(keep):
        if vdeg[k] > 0:
            fh.write(f'J8[size(J8)+1] = x({i+1});\n')
    fh.write('''
int t = timer; system("--ticks-per-sec", 1);
"nvars =", nvars(r), " nquadrics =", size(I), " nP8coords =", nk - size(J8);
ideal G = std(I);
"std(I): dim =", dim(G), " mult =", mult(G), " size =", size(G), " time =", timer - t;
intvec hv = hilb(G, 1); "hilbert numerator:"; hv;
t = timer;
ideal GJ = std(I + J8);
"std(I + J8) [the P8 component]: dim =", dim(GJ), " mult =", mult(GJ), " time =", timer - t;
t = timer;
ideal Isat = sat(I, I + J8);
ideal Gs = std(Isat);
"sat(I, I+J8): dim =", dim(Gs), " mult =", mult(Gs), " size =", size(Gs), " time =", timer - t;
"dim -1 (or 0) means: the relaxation variety is exactly the P8 component (mod p)";
if (dim(Gs) >= 1) { "OTHER COMPONENTS PRESENT: dim", dim(Gs), "mult", mult(Gs); intvec hv2 = hilb(Gs, 1); hv2; }
quit;
'''.replace('nk - size(J8)', str(nk) + ' - size(J8)'))
meta = {'experiment': 'tensor-relaxation-singular-export', 'd': d, 'p': p, 'K': K, 'nvars': nk, 'nquadrics': len(lines), 'coords': int(len(coords)),
        'P8_vars': sum(1 for k in keep if vdeg[k] == 0), 'seconds': round(time.time() - t0, 1)}
json.dump(meta, open(name + '.json', 'w'), indent=1); L.log_line(f'w2-tensor-relax-export.jsonl', meta)
print(json.dumps(meta), flush=True)
