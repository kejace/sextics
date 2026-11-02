"""Number of distinct points (over the algebraic closure of F_p) of a zero-dimensional projective scheme given by
quadrics in an msolve file: dehomogenize at a random affine chart containing the Euler point (sum of the
Euler coordinates = 1), take the radical, count vdim."""
import sys, subprocess, random
path = sys.argv[1]; eul = sys.argv[2]   # eul: comma list of variable indices with Euler coordinate 1
lines = [l.strip() for l in open(path) if l.strip()]
vars_ = lines[0]; p = lines[1]
polys = [s.rstrip(',') for s in lines[2:] if s.rstrip(',')]
eidx = [int(x) for x in eul.split(',')]
random.seed(1)
nv = len(vars_.split(','))
# chart: ell = sum_{i in eul} w_i + random combination of the others (so E is in the chart and the chart is generic)
ell = ' + '.join([f'w{i}' for i in eidx] + [f'{random.randint(1, 1000)}*w{i}' for i in range(nv) if i not in eidx])
sing = f"""LIB "primdec.lib";
ring r = {p}, ({vars_}), dp;
option(redSB);
ideal I = {','.join(polys)};
ideal Ja = I, ({ell}) - 1;
ideal G = std(Ja);
"VDIM affine =", vdim(G);
ideal R = radical(G);
ideal GR = std(R);
"VDIM radical =", vdim(GR);
"RADICAL GENS", size(GR);
GR;
quit;
"""
open(path + '.rad.sing', 'w').write(sing)
r = subprocess.run(['Singular', '-q', path + '.rad.sing'], capture_output=True, text=True, timeout=int(sys.argv[3]) if len(sys.argv) > 3 else 900)
open(path + '.rad.out', 'w').write(r.stdout + '\n' + r.stderr)
print('\n'.join(l for l in r.stdout.splitlines() if l.startswith('VDIM') or l.startswith('RADICAL'))); print(r.stdout[-1500:] if 'VDIM radical' in r.stdout else r.stderr[-300:])
