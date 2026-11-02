"""Table of all characteristic-0 msolve results for the degenerate-polar loci."""
import glob, os, sympy as sp
from parse_q import parse_q
rows = []
for name in sorted(glob.glob('sigmaQP_*.out')):
    if os.path.getsize(name) == 0: rows.append((name, 'EMPTY FILE (not finished)')); continue
    txt = open(name).read().strip()
    if txt.startswith('[-1]'): rows.append((name, 'no solutions over Qbar')); continue
    try:
        deg, d, fac = parse_q(name)
        rows.append((name, f'ideal degree {deg}, #points over Qbar = {d}, factor degrees {fac}'))
    except Exception as ex:
        rows.append((name, 'PARSE ERROR ' + str(ex)[:60]))
for r in rows: print(r[0].ljust(32), r[1])
