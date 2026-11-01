"""Verdict on a homogeneous quadric system file (msolve format): msolve (with dummy var) + Singular dim cross-check.
Usage: verdict.py file.ms  -> prints JSON {msolve_dim, msolve_nsol, singular_dim, singular_vdim}"""
import sys, subprocess, re, json, os
f = sys.argv[1]
lines = open(f).read().strip().split('\n')
vars_ = lines[0].split(','); p = lines[1]; eqs = [l.rstrip(',') for l in lines[2:] if l.strip()]
# msolve with dummy variable
ms2 = f + '.dummy.ms'
open(ms2, 'w').write(','.join(vars_ + ['zz']) + '\n' + p + '\n' + ',\n'.join(eqs) + '\n')
r = subprocess.run(['nice', '-n', '5', 'timeout', '900', 'msolve', '-v', '1', '-t', '2', '-f', ms2, '-o', ms2 + '.out'], capture_output=True, text=True)
log = r.stdout + r.stderr
out = open(ms2 + '.out').read() if os.path.exists(ms2 + '.out') else ''
msd = 'positive' if 'positive dimension' in log.lower() or out.startswith('[1,') or out.startswith('[2,') else ('zero' if out.startswith('[0,') else 'unknown')
m = re.search(r'Dimension of quotient:\s*(\d+)', log)
quot = int(m.group(1)) if m else None
# Singular cross-check
sg = f + '.sing'
with open(sg, 'w') as g:
    g.write(f'ring r = {p}, ({",".join(vars_)}), dp;\nideal I = {",".join(eqs)};\nideal G = std(I);\nprint("SDIM " + string(dim(G)));\nif (dim(G) == 0) {{ print("SVDIM " + string(vdim(G))); }}\nquit;\n')
r2 = subprocess.run(['nice', '-n', '5', 'timeout', '900', 'Singular', '-q', sg], capture_output=True, text=True)
sd = re.search(r'SDIM (-?\d+)', r2.stdout); sv = re.search(r'SVDIM (\d+)', r2.stdout)
res = {'file': os.path.basename(f), 'nvars': len(vars_), 'neqs': len(eqs), 'msolve_dim': msd, 'msolve_quotient_dim': quot,
       'singular_dim': int(sd.group(1)) if sd else None, 'singular_vdim': int(sv.group(1)) if sv else None, 'singular_rc': r2.returncode,
       'msolve_out_head': out[:120].replace('\n', ' ')}
print(json.dumps(res))
