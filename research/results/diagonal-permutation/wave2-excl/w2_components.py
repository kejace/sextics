"""Analyse a positive-dimensional relaxation variety: Singular dim/degree/minimal associated primes of the GB ideal,
then test rational points of each component (if linear) with the reduced E2' system.

usage: python3 w2_components.py <tag> [--timeout S]     (reads relax_<tag>.gb written by msolve -g 2, and relax_<tag>.ms)
Writes relax_<tag>.components.json; logs exp='relax_components'.  Only feasible for moderate systems (<= ~70 variables)."""
import sys, os, re, json, time, subprocess
import w2lib as W

tag = sys.argv[1]
timeout = int(sys.argv[sys.argv.index('--timeout') + 1]) if '--timeout' in sys.argv else 2 * 3600
ms = open(f'{W.PKT}/relax_{tag}.ms').read().splitlines()
names = ms[0].split(','); p = int(ms[1])
gb = open(f'{W.PKT}/relax_{tag}.gb').read()
body = '\n'.join(l for l in gb.splitlines() if not l.startswith('#')).strip().lstrip('[').rstrip(':').rstrip(']')
polys = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
sing = f'{W.PKT}/relax_{tag}.components.sing'
with open(sing, 'w') as f:
    f.write('LIB "primdec.lib";\n')
    f.write(f'ring r = {p}, ({",".join(names)}), dp;\n')
    f.write('ideal G = ' + ',\n'.join(polys) + ';\n')
    f.write('attrib(G, "isSB", 1);\n')
    f.write('print("DIM " + string(dim(G)));\n')
    f.write('print("MULT " + string(mult(G)));\n')
    f.write('list L = minAssGTZ(G);\n')
    f.write('print("NCOMP " + string(size(L)));\n')
    f.write('int i; for (i = 1; i <= size(L); i++) { ideal P = std(L[i]); print("COMP " + string(i) + " dim " + string(dim(P)) + " mult " + string(mult(P)) + " ngens " + string(size(L[i])) + " maxdeg " + string(deg(L[i][size(L[i])]))); print(L[i]); }\n')
    f.write('quit;\n')
t0 = time.time()
try:
    pr = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=timeout)
    out = pr.stdout + pr.stderr; status = pr.returncode
except subprocess.TimeoutExpired:
    out = ''; status = 'timeout'
secs = round(time.time() - t0, 1)
open(f'{W.PKT}/relax_{tag}.components.out', 'w').write(out)
dim_ = re.search(r'DIM (\S+)', out); ncomp = re.search(r'NCOMP (\S+)', out)
comps = re.findall(r'COMP (\d+) dim (\d+) mult (\d+) ngens (\d+) maxdeg (\d+)', out)
rec = {'tag': tag, 'status': status, 'dim': dim_.group(1) if dim_ else None, 'ncomp': ncomp.group(1) if ncomp else None,
       'components': comps, 'seconds': secs}
json.dump(rec, open(f'{W.PKT}/relax_{tag}.components.json', 'w'), indent=1)
W.log('relax_components', family=tag, prime=p, result=f"dim {rec['dim']}, {rec['ncomp']} components: {comps}", seconds=secs)
print('DONE', rec, flush=True)
