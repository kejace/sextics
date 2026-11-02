"""Stage 2 of the Hessian relaxation: Groebner basis of the quadric ideal (msolve, char p) and the
m-primary test (every variable has a pure power among the leading monomials => affine variety = {0}
=> projective relaxation variety empty).  Also optionally an independent Singular dim() check.
usage: python3 wave2-excl/gb_check.py <file.ms> [--threads T] [--singular]
"""
import sys, os, re, json, time, subprocess, argparse
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib

def parse_ms(path):
    lines = open(path).read().split('\n')
    vars_ = lines[0].split(','); p = int(lines[1])
    polys = [l.strip().rstrip(',') for l in lines[2:] if l.strip()]
    return vars_, p, polys

def leading_pure_powers(gbtext, vars_):
    """return set of variables having a pure power as the leading (first-printed) monomial of some GB element."""
    body = gbtext
    body = '\n'.join(l for l in body.split('\n') if not l.startswith('#'))
    body = body.strip().strip(':').strip()
    if body.startswith('['): body = body[1:]
    if body.endswith(']'): body = body[:-1]
    polys = [s.strip() for s in body.split(',') if s.strip()]
    pure = set(); lead = []
    for P in polys:
        first = re.split(r'(?<!\^)[+-]', P.strip('+'), maxsplit=1)[0] if not P.startswith('-') else re.split(r'(?<!\^)[+-]', P[1:], maxsplit=1)[0]
        first = re.sub(r'^\s*\d+\s*\*?', '', first).strip()
        lead.append(first)
        m = re.fullmatch(r'([A-Za-z_]\w*)(\^\d+)?', first)
        if m and m.group(1) in vars_:
            pure.add(m.group(1))
    return pure, polys, lead

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('msfile'); ap.add_argument('--threads', type=int, default=16); ap.add_argument('--singular', action='store_true')
    ap.add_argument('--tag', default=None)
    a = ap.parse_args()
    vars_, p, polys = parse_ms(a.msfile)
    tag = a.tag or os.path.basename(a.msfile).replace('.ms', '')
    out = a.msfile.replace('.ms', '.gb')
    t0 = time.time()
    cmd = ['msolve', '-v', '1', '-t', str(a.threads), '-g', '2', '-f', a.msfile, '-o', out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    t1 = time.time()
    res = {'exp': 'groebner', 'file': a.msfile, 'nvars': len(vars_), 'npolys': len(polys), 'prime': p, 'msolve_rc': r.returncode,
           'msolve_seconds': round(t1 - t0, 1), 'msolve_stderr_tail': r.stderr[-600:], 'msolve_stdout_tail': r.stdout[-600:]}
    verdict = 'unknown'
    if r.returncode == 0 and os.path.exists(out):
        txt = open(out).read()
        pure, gbpolys, lead = leading_pure_powers(txt, vars_)
        res['gb_size'] = len(gbpolys); res['pure_power_vars'] = len(pure); res['lead_sample'] = lead[:10]
        if gbpolys == ['1'] or (len(gbpolys) == 1 and gbpolys[0].strip() == '1'):
            verdict = 'unit ideal (impossible for homogeneous input; check)'
        elif len(pure) == len(vars_):
            verdict = 'EMPTY projective variety over F_p-bar (ideal is m-primary: pure powers of all variables lead)'
        else:
            verdict = f'NOT m-primary by leading-term test: {len(vars_) - len(pure)} variables without pure-power leading term (positive-dimensional or parse issue)'
    res['verdict'] = verdict
    if a.singular:
        sfile = a.msfile.replace('.ms', '.sing')
        if os.path.exists(sfile):
            t2 = time.time()
            rs = subprocess.run(['Singular', '-q', sfile], capture_output=True, text=True)
            res['singular_seconds'] = round(time.time() - t2, 1); res['singular_out'] = rs.stdout[-1500:]
            m = re.search(r'dim\(G\) = (\d+)', rs.stdout)
            if m: res['singular_dim'] = int(m.group(1))
    res['seconds'] = round(time.time() - t0, 1); res['job'] = os.environ.get('W2JOB', 'gb_check')
    json.dump(res, open(a.msfile.replace('.ms', '.gb.receipt.json'), 'w'), indent=1)
    w2lib.log(exp=f'groebner-{tag}', family=tag, prime=p, result=verdict, rank={'nvars': len(vars_), 'npolys': len(polys), 'gb_size': res.get('gb_size')},
              seconds=res['seconds'], job=res['job'], note=f"msolve rc {r.returncode}; singular_dim={res.get('singular_dim')}")
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main()
