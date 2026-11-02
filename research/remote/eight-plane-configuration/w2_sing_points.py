"""Primary decomposition (char p) of a zero-dimensional homogeneous quadric ideal from an msolve file; prints
each component (degree, dimension) and, for F_p-rational points, the normalized coordinates."""
import sys, subprocess, re, json
path = sys.argv[1]
lines = [l.strip() for l in open(path) if l.strip()]
vars_ = lines[0]; p = lines[1]
polys = [s.rstrip(',') for s in lines[2:] if s.rstrip(',')]
sing = f"""LIB "primdec.lib";
ring r = {p}, ({vars_}), dp;
option(redSB);
ideal I = {','.join(polys)};
ideal G = std(I);
list L = primdecGTZ(G);
"NCOMP", size(L);
int i; int j;
for (i = 1; i <= size(L); i++) {{
  ideal P = std(L[i][2]); ideal Qc = std(L[i][1]);
  "COMP", i, "primedim", dim(P), "primedeg", degree(P), "primarydeg", degree(Qc), "ngens", size(P);
  "GENS"; P;
}}
quit;
"""
open(path + '.pd.sing', 'w').write(sing)
r = subprocess.run(['Singular', '-q', path + '.pd.sing'], capture_output=True, text=True, timeout=1800)
open(path + '.pd.out', 'w').write(r.stdout + '\n' + r.stderr)
print(r.stdout[:6000]); print(r.stderr[-500:])
