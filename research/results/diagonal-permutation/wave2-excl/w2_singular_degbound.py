"""Degree-truncated Groebner computation (Singular, degBound) for a relaxation receipt relax_*.ms.

usage: python3 w2_singular_degbound.py <tag> <degBound> [--timeout S]
Reads relax_<tag>.ms (msolve format), writes relax_<tag>.sing and runs Singular with option(redSB), degBound = k.
The truncated standard basis gives leading terms L_k that are a subset of in(I); if every variable has a pure power
among them, then S/I is finite-dimensional, i.e. V(I) = {0} (homogeneous ideal) -- rigorous as in w2_relax.py.
Also prints the number of degree-k standard monomials of the truncated leading ideal (0 means I ⊇ m^k, also rigorous).
Logs exp='relax_singular'."""
import sys, os, time, json, subprocess, re
import w2lib as W

tag = sys.argv[1]; k = int(sys.argv[2])
timeout = int(sys.argv[sys.argv.index('--timeout') + 1]) if '--timeout' in sys.argv else 4 * 3600
ms = open(f'{W.PKT}/relax_{tag}.ms').read().splitlines()
names = ms[0].split(','); p = int(ms[1]); polys = [l.rstrip(',') for l in ms[2:] if l.strip()]
n = len(names)
sing = f'{W.PKT}/relax_{tag}_db{k}.sing'
with open(sing, 'w') as f:
    f.write(f'ring r = {p}, ({",".join(names)}), dp;\n')
    f.write('option(redSB); option(redTail);\n')
    f.write(f'degBound = {k};\n')
    f.write('ideal I = ' + ',\n'.join(polys) + ';\n')
    f.write('int t0 = timer; ideal G = std(I); int t1 = timer;\n')
    f.write('print("TIME " + string(t1 - t0));\n')
    f.write('ideal L = lead(G);\n')
    f.write('print("NGB " + string(size(G)));\n')
    # pure powers per variable and degree-k standard monomials count
    f.write('int i; int cnt = 0;\n')
    f.write(f'for (i = 1; i <= {n}; i++) {{ int j; int found = 0; for (j = 1; j <= size(L); j++) {{ if (leadmonom(L[j]) == var(i)^deg(L[j])) {{ found = 1; break; }} }} cnt = cnt + found; }}\n')
    f.write('print("PUREPOWERS " + string(cnt));\n')
    f.write(f'print("HILB_k " + string(size(kbase(L, {k}))));\n')
    f.write('quit;\n')
t0 = time.time()
try:
    pr = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=timeout)
    out = pr.stdout + pr.stderr; status = pr.returncode
except subprocess.TimeoutExpired:
    out = ''; status = 'timeout'
secs = round(time.time() - t0, 1)
open(f'{W.PKT}/relax_{tag}_db{k}.out', 'w').write(out[-5000:])
m = {key: (re.search(key + r' (\S+)', out).group(1) if re.search(key + r' (\S+)', out) else None) for key in ('TIME', 'NGB', 'PUREPOWERS', 'HILB_k')}
verdict = 'unknown'
if m['PUREPOWERS'] is not None and int(m['PUREPOWERS']) == n:
    verdict = 'V(I) = {0} (every variable has a pure power among truncated leading terms)'
elif m['HILB_k'] is not None and int(m['HILB_k']) == 0:
    verdict = f'V(I) = {{0}} (truncated leading ideal contains all monomials of degree {k})'
elif status == 0:
    verdict = f'inconclusive at degBound {k}: pure powers {m["PUREPOWERS"]}/{n}, degree-{k} standard monomials {m["HILB_k"]}'
else:
    verdict = f'Singular status {status}'
W.log('relax_singular', family=tag, d=int(re.search(r'd(\d+)', tag).group(1)), prime=p, result=verdict, seconds=secs,
      note=f'degBound {k}; GB size {m["NGB"]}; singular time {m["TIME"]}')
print('DONE', verdict, m, flush=True)
