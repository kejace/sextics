"""Convert an msolve quadric file to Singular input and compute dim/degree/Hilbert data (char p)."""
import sys, subprocess, re
path = sys.argv[1]
lines = [l.strip() for l in open(path) if l.strip()]
vars_ = lines[0]; p = lines[1]
polys = [s.rstrip(',') for s in lines[2:]]
polys = [s for s in polys if s]
sing = f"""ring r = {p}, ({vars_}), dp;
option(redSB);
ideal I = {','.join(polys)};
ideal G = std(I);
int d = dim(G); int dg = degree(G);
"dim(affine cone) =", d, " projective dim =", d-1, " degree =", dg;
intvec hv = hilb(G, 1);
"first Hilbert series numerator =", hv;
"Hilbert function values (deg 0..8):";
int k; for (k = 0; k <= 8; k++) {{ kbase(G, k); "deg", k, ":", size(kbase(G, k)); }}
quit;
"""
open(path + '.sing', 'w').write(sing)
r = subprocess.run(['Singular', '-q', path + '.sing'], capture_output=True, text=True, timeout=600)
out = r.stdout
# keep only summary lines
keep = [l for l in out.splitlines() if l.startswith('dim') or l.startswith('first') or l.startswith('deg') or l.startswith('Hilbert')]
print('\n'.join(keep)); print(r.stderr[-300:])
