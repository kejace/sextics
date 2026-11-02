"""msolve (char p) screening of the cone {w : w w^T in K_H(p)}: random affine charts (sum r_a w_a = 1).
usage: python3 wave2_excl_hess_modp_msolve.py <KH json> [timeout] [threads]"""
import sys, os, json, time, subprocess, random
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import logline
from flint import nmod_mat
path = sys.argv[1]; TO = int(sys.argv[2]) if len(sys.argv) > 2 else 600; TH = sys.argv[3] if len(sys.argv) > 3 else '2'
K = json.load(open(path))
p = K['p']; Wkeys = [tuple(x) for x in K['Wkeys']]; nw = len(K['wb']); nW = len(Wkeys)
basis = [{int(a): v % p for a, v in row.items()} for row in K['basis']]
rk = len(basis)
t0 = time.time()
Bm = nmod_mat(rk, nW, [basis[i].get(a, 0) for i in range(rk) for a in range(nW)], p)
Lam, nlam = Bm.nullspace()
quadrics = []
for k in range(nlam):
    terms = []
    for idx, (a, b) in enumerate(Wkeys):
        c = int(Lam[idx, k])
        if c:
            terms.append(f'{c}*w{a}*w{b}' if a != b else f'{c}*w{a}^2')
    quadrics.append('+'.join(terms))
out = {}
for chart in (1, 2):
    random.seed(chart + 10)
    r = [random.randint(1, 1000) for _ in range(nw)]
    lin = '+'.join(f'{c}*w{a}' for a, c in enumerate(r)) + f'+{p-1}'
    f = path.replace('.json', f'.chart{chart}.ms')
    with open(f, 'w') as fh:
        fh.write(','.join(f'w{a}' for a in range(nw)) + f'\n{p}\n' + ',\n'.join(quadrics + [lin]) + '\n')
    o = f.replace('.ms', '.out')
    try:
        rr = subprocess.run(['msolve', '-v', '1', '-t', TH, '-P', '1', '-f', f, '-o', o], capture_output=True, text=True, timeout=TO)
        res = open(o).read() if os.path.exists(o) else ''
        msg = [l for l in rr.stdout.split('\n') if 'solution' in l.lower() or 'dimension' in l.lower() or 'basis' in l.lower()]
        out[f'chart{chart}'] = {'rc': rr.returncode, 'out_head': res[:200], 'out_len': len(res), 'msg': msg[-4:], 'seconds': round(time.time() - t0, 1)}
    except subprocess.TimeoutExpired:
        out[f'chart{chart}'] = {'timeout': TO}
        break
rec = dict(exp='hess-modp-msolve', file=os.path.basename(path), d=K['d'], parity=K['parity'], prime=p, dim_KH=rk, nw=nw, n_quadrics=nlam, msolve=out,
           seconds=round(time.time() - t0, 1), note='char-p screening; msolve syntax uses w^2 for squares')
logline(os.path.join(w2x.LANE, 'wave2-excl.hessmodp.jsonl'), rec)
print(json.dumps(rec))
