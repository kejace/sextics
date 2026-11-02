import pickle, subprocess, re, time, json, sys
P2 = 1000000007
def write_ms(pkl, out, p):
    d = pickle.load(open(pkl, "rb")); alive = sorted(d['alive']); quads = d['quads']
    names = [f'l{a}' for a in alive]; lines = []
    for q in quads:
        terms = []
        for (a,b), v in q['coef'].items():
            if a not in alive or b not in alive: continue
            vv = int(v) % p
            if vv: terms.append(f'{vv}*l{a}*l{b}' if a != b else f'{vv}*l{a}^2')
        if terms: lines.append('+'.join(terms))
    lines = list(dict.fromkeys(lines))
    open(out, 'w').write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(lines) + '\n')
    return len(lines), len(names), len(quads)
def pure_powers(leadfile):
    txt = open(leadfile).read()
    vars_ = re.search(r'variable order:\s*(.*)', txt).group(1).replace(' ','').split(',')
    body = txt[txt.index('['):]
    monos = [m.strip() for m in body.strip('[]:\n ').split(',') if m.strip()]
    pure = {}
    for m in monos:
        fac = m.split('*')
        if len(fac) == 1:
            v, _, k = fac[0].partition('^'); pure[v] = int(k or 1)
    return len(monos), pure, [v for v in vars_ if v not in pure]
R = '/tmp/hodge/remote/one-reflection-skew/wave2-excl/'
jobs = [('5v','0010', R+'d7_5v_tor_0010.pkl'), ('5v','1000', R+'d7_5v_tor_1000.pkl'),
        ('6v','0100', R+'d7_sub_tor_0100.pkl'), ('6v','0010', R+'d7_sub_tor_0010.pkl'), ('6v','1000', R+'d7_sub_tor_1000.pkl')]
for kind, case, pkl in jobs:
    tag = f'ref{kind}_tor_{case}_p2'
    nl, nv, nq = write_ms(pkl, tag+'.ms', P2)
    t0 = time.time()
    try:
        subprocess.run(['nice','-n','10','msolve','-t','1','-g','1','-f',tag+'.ms','-o',tag+'.lead'], capture_output=True, text=True, timeout=1500)
        n, pure, missing = pure_powers(tag+'.lead')
        res = dict(exp=f'msolve-second-prime-{kind}', case=case, prime=P2, n_leading=n, missing=missing, pure_degrees=sorted(set(pure.values())), nquads=nl, nvars=nv, seconds=round(time.time()-t0,1), note='inputs regenerated from exact integer quadrics in the lane pkl (remote mirror); lane prime 1073741827')
    except Exception as e:
        res = dict(exp=f'msolve-second-prime-{kind}', case=case, prime=P2, error=repr(e)[:200], seconds=round(time.time()-t0,1))
    print(json.dumps(res), flush=True)
    with open('referee.log.jsonl','a') as fh: fh.write(json.dumps(res)+'\n')
