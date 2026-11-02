"""Rigorous per-denominator test (q = 1) in Singular mod p:  E1 (A D in J) and E2-mod-J (B_D(N) in J_{2d} for one
E1-solution N obtained by lift).  Both are membership tests of fixed polynomials in J, so with the Hilbert check
(std(J) has the smooth Hilbert function) a failure mod p is a failure over Q, and for a fixed rational D the linear
certificate system has a complex solution iff it has a rational one.  Pure-parity D only (mixed D: test parts).
Input: JSON list of {name, family, D} (D as a polynomial string, rational or already reduced mod p).
Usage: python3 w2x_e2fixed.py <cands.json> <tag> [--run]"""
import sys, os, time, json, subprocess
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave2-excl')
from w2x_lib import *
from hodge.core import parse, to_str, degree, parity_of

p = P
def write_sing(path, cands):
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    with open(path, 'w') as f:
        f.write(f'ring R = {p}, (X0,X1,X2,Y0,Y1,Y2), dp;\n')
        f.write('poly F = ' + to_str(Fp).replace('**', '^') + ';\n')
        f.write('poly A = ' + to_str(Ap).replace('**', '^') + ';\n')
        f.write('ideal J = jacob(F);\nint t0 = timer;\nideal G = std(J);\nprint("std seconds = " + string(timer - t0));\n')
        f.write('intvec hv = hilb(G, 2);\nprint("hilb2 = " + string(hv));\n')
        f.write('poly D; poly AD; matrix T; poly divN; poly ND; poly B; int i; int e1; int e2;\n')
        for it in cands:
            D = core.reduce_mod(parse(it['D']), p)
            assert parity_of(D) is not None, it['name']
            f.write('D = ' + to_str(D).replace('**', '^') + ';\n')
            f.write('AD = A*D; t0 = timer; e1 = (reduce(AD, G) == 0); e2 = -1;\n')
            f.write('if (e1) { T = lift(J, AD); divN = 0; ND = 0; for (i = 1; i <= 6; i++) { divN = divN + diff(T[i,1], var(i)); ND = ND + T[i,1]*diff(D, var(i)); }\n')
            f.write('  B = D*divN - ND; e2 = (reduce(B, G) == 0); }\n')
            f.write(f'print("RES|{it["name"]}|{it.get("family","")}|{degree(D)}|{parity_of(D)}|" + string(e1) + "|" + string(e2) + "|" + string(timer - t0));\n')
        f.write('quit;\n')

if __name__ == '__main__':
    cands = json.load(open(sys.argv[1])); tag = sys.argv[2]
    path = f'{PK}/e2fixed_{tag}.sing'
    write_sing(path, cands)
    print(f'{len(cands)} candidates -> {path}', flush=True)
    if '--run' in sys.argv:
        t0 = time.time()
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=72000)
        out = r.stdout + r.stderr
        with open(f'{PK}/e2fixed_{tag}.out', 'w') as f: f.write(out)
        hilb = next((l.split('=', 1)[1].strip() for l in out.splitlines() if l.startswith('hilb2 =')), None)
        expected = ','.join(str(hilb_R(k)) for k in range(25))
        res = []
        for l in out.splitlines():
            if l.startswith('RES|'):
                _, nm, fam, d, par, e1, e2, sec = l.split('|')
                res.append({'name': nm, 'family': fam, 'd': int(d), 'parD': int(par), 'E1': e1 == '1', 'E2modJ': (e2 == '1') if e2 != '-1' else None, 'seconds': int(sec)})
        summary = {'exp': 'e2fixed', 'tag': tag, 'p': p, 'n': len(res), 'hilbert_ok': hilb == expected,
                   'E1_pass': [r_['name'] for r_ in res if r_['E1']], 'E2modJ_pass': [r_['name'] for r_ in res if r_['E2modJ']],
                   'seconds': round(time.time() - t0, 1)}
        with open(f'{PK}/e2fixed_{tag}.json', 'w') as f: json.dump({'summary': summary, 'results': res}, f, indent=1)
        for r_ in res: log(dict(r_, exp='e2fixed-item', tag=tag, p=p, hilbert_ok=summary['hilbert_ok']))
        log(summary); print(json.dumps(summary), flush=True)
