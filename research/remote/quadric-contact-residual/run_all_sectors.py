"""Gauge-fixed quadratic E2 closure for D in every G_aff character sector at degree d, then msolve. Usage: run_all_sectors.py d"""
import sys, json, time, subprocess, os
from inv import *
d = int(sys.argv[1])
chars = sorted({chi(e) for e in monomials(d)})
results = {}
outf = f'/tmp/hodge/results/quadric-contact-residual/sectors_d{d}.json'
if os.path.exists(outf):
    results = json.load(open(outf))
for c in chars:
    key = str(c)
    if key in results and 'msolve' in results[key]:
        continue
    t0 = time.time()
    try:
        res, Ds, Ns, quads = quadratic_closure(d, verbose=False, out_prefix='/tmp/hodge/results/quadric-contact-residual/closure', c=c, sym=False)
    except Exception as ex:
        results[key] = {'error': repr(ex)}; json.dump(results, open(outf, 'w')); continue
    res = dict(res)
    if res.get('n_params', 0) == 0:
        res['msolve'] = 'no parameters'
    else:
        tag = f'_d{d}_c' + ''.join(str(x) for x in c)
        msf = f'/tmp/hodge/results/quadric-contact-residual/closure_quads{tag}.ms'
        nvars = res['n_params'] - res['n_gauge']
        nq = sum(1 for _ in open(msf)) - 2
        if nvars <= 0:
            res['msolve'] = 'only gauge'
        elif nq == 0:
            res['msolve'] = 'NO EQUATIONS (cannot close)'
        else:
            try:
                out = subprocess.run(['nice', '-n', '5', 'timeout', '600', 'msolve', '-v', '0', '-t', '2', '-f', msf, '-o', msf + '.out'], capture_output=True, text=True)
                txt = open(msf + '.out').read() if os.path.exists(msf + '.out') else ''
                res['msolve'] = txt[:400].replace('\n', ' ')
                res['msolve_rc'] = out.returncode
            except Exception as ex:
                res['msolve'] = 'ERR ' + repr(ex)
        res['nvars_after_gauge'] = nvars; res['n_quadrics_after_gauge'] = nq
    res['seconds_total'] = round(time.time() - t0, 1)
    results[key] = res
    json.dump(results, open(outf, 'w'), indent=0)
    print(key, {k: res.get(k) for k in ('n_params', 'n_gauge', 'n_extra', 'dimR_2d_sector', 'n_quadrics', 'quadric_rank', 'msolve', 'seconds_total')}, flush=True)
