"""Independent Groebner engine check: convert an msolve input (.ms) to Singular, compute std with redSB,
report dim (Krull dimension of the quotient; 0 for a homogeneous ideal means V(I) = {0}) and vdim (colength)."""
import sys, subprocess, json, re, time
ms = sys.argv[1]
lines = [l for l in open(ms).read().splitlines() if l.strip()]
names = lines[0].strip(); p = int(lines[1].strip())
polys = [l.strip().rstrip(',') for l in lines[2:]]
polys = [q for q in polys if q]
sing = f'ring r = {p}, ({names}), dp;\nideal I = ' + ',\n'.join(polys) + ';\noption(redSB);\nideal G = std(I);\nprint("NGB=" + string(size(G)));\nprint("DIM=" + string(dim(G)));\nprint("VDIM=" + string(vdim(G)));\nquit;\n'
t0 = time.time()
r = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=3600)
out = r.stdout
rec = {'ms': ms, 'prime': p, 'n_vars': len(names.split(',')), 'n_eqs': len(polys),
       'singular_NGB': int(re.search(r'NGB=(\d+)', out).group(1)), 'singular_dim': int(re.search(r'DIM=(-?\d+)', out).group(1)),
       'singular_vdim': int(re.search(r'VDIM=(-?\d+)', out).group(1)), 'seconds': round(time.time() - t0, 2),
       'verdict': None}
rec['verdict'] = 'dim 0 for a homogeneous ideal => V(I) = {0} over the algebraic closure of GF(p)' if rec['singular_dim'] == 0 else f'positive dimensional (dim {rec["singular_dim"]})'
print(json.dumps(rec)); open('ref_singular_check.jsonl', 'a').write(json.dumps(rec) + '\n')
