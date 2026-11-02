"""m-primary verdict for a homogeneous quadric system in msolve format via `msolve -g 1` (leading ideal, DRL).
Usage: python3 msverdict.py file.ms [threads] [timeout_s]  -> prints JSON {mprimary, pure_powers, n_lead, rc, seconds}
mprimary=True  <=> every variable has a pure power among the leading monomials  <=> only c=0 over the algebraic closure of F_p.
(Default msolve solving mode segfaults for primes >= 2^31; squares must be written ck^2, never ck*ck.)"""
import sys, subprocess, re, json, time, os
f = sys.argv[1]; threads = sys.argv[2] if len(sys.argv) > 2 else '2'; to = sys.argv[3] if len(sys.argv) > 3 else '3600'
names = open(f).readline().strip().split(',')
out = f + '.g1.out'
t0 = time.time()
r = subprocess.run(['timeout', to, 'msolve', '-g', '1', '-t', threads, '-f', f, '-o', out], capture_output=True, text=True)
txt = open(out).read() if os.path.exists(out) else ''
body = txt[txt.find('['):] if '[' in txt else ''
monos = [m.strip() for m in body.strip().lstrip('[').rstrip(':').rstrip(']').split(',') if m.strip()]
pure = {}
for m in monos:
    parts = m.split('*')
    if len(parts) == 1:
        mm = re.match(r'(\w+)\^(\d+)$', parts[0]) or re.match(r'(\w+)$', parts[0])
        if mm:
            v = mm.group(1); e = int(mm.group(2)) if mm.lastindex and mm.lastindex >= 2 else 1
            pure[v] = min(pure.get(v, 10**9), e)
res = {'file': os.path.basename(f), 'nvars': len(names), 'rc': r.returncode, 'seconds': round(time.time() - t0, 1), 'n_lead': len(monos),
       'pure_powers': pure, 'mprimary': (r.returncode == 0 and all(v in pure for v in names))}
print(json.dumps(res))
