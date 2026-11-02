"""msolve -g 1 verdict for a homogeneous ideal: V(I) = {0} over the algebraic closure of F_p  <=>  the leading ideal
(grevlex) contains a pure power of every variable (monomial ideal is m-primary).  Usage: msverdict.py file.ms threads timeout"""
import sys, subprocess, re, os, time, json
def verdict(ms_file, threads=4, timeout=3600, verbose=2, logfile=None):
    t0 = time.time()
    out = ms_file + '.lead'
    cmd = ['timeout', str(timeout), 'msolve', '-v', str(verbose), '-t', str(threads), '-g', '1', '-f', ms_file, '-o', out]
    if logfile:
        with open(logfile, 'w') as lf:
            r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)
        tail = open(logfile).read()[-600:]
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
        tail = (r.stdout + r.stderr)[-600:]
    res = {'rc': r.returncode, 'seconds': round(time.time() - t0, 1), 'log_tail': tail}
    if not os.path.exists(out) or r.returncode != 0:
        res['verdict'] = 'msolve-failed-or-timeout'
        return res
    txt = open(out).read()
    vars_ = open(ms_file).readline().strip().split(',')
    m = re.search(r'\[(.*)\]:', txt, re.S)
    if not m:
        res['verdict'] = 'unparsed'; return res
    lead = [s.strip() for s in m.group(1).replace('\n', '').split(',') if s.strip()]
    res['n_lead'] = len(lead)
    if lead == ['1']:
        res['verdict'] = 'UNIT-IDEAL (impossible for homogeneous)'; return res
    powers = {}
    for mono in lead:
        mm = re.fullmatch(r'([A-Za-z_]\w*)(?:\^(\d+))?', mono)
        if mm:
            v, k = mm.group(1), int(mm.group(2) or 1)
            powers[v] = min(powers.get(v, 10**9), k)
    missing = [v for v in vars_ if v not in powers]
    res['pure_powers'] = powers
    res['missing_variables'] = missing
    res['max_pure_power'] = max(powers.values()) if powers else None
    res['verdict'] = 'ONLY-ORIGIN (m-primary leading ideal)' if not missing else 'POSITIVE-DIMENSIONAL'
    return res
if __name__ == '__main__':
    f = sys.argv[1]; th = int(sys.argv[2]) if len(sys.argv) > 2 else 4; to = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
    v = verdict(f, th, to)
    print(json.dumps(v))
