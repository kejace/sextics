"""E1 screening of structured denominators by Groebner normal form modulo J(F) in Singular (mod p).
Independent of the linear-algebra E1 tests of wave 1.  A D in J is necessary for a certificate (q=1):
for mixed-parity D = D+ + D- it is equivalent to A D+ in J and A D- in J (J is tau-stable, A is Y-odd).
Rigor: the script prints the Hilbert function of S/J mod p; if it equals (1+t+..+t^4)^6 in all degrees used,
J mod p is the reduction of J over Q and a nonzero normal form proves A D not in J over Q.
Usage: python3 w2x_e1screen.py <family-set> [--run]   ->  e1screen_<set>.sing / .out / .json"""
import sys, os, time, json, itertools, subprocess
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave2-excl')
from w2x_lib import *
from hodge.core import parse, mul, add, sub, scale, power, degree, parity_of, to_str, monomials

p = P
Z = [parse(v) for v in ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']]
L = parse('X0 + 2*X1 + 3*X2'); M = parse('Y0 + 3*Y1 + 5*Y2')
c = Fraction(5, 7 ** 7)
l7 = add(L, M); l8 = sub(L, M)            # kappa dropped (scalar)
ell = Z + [l7, l8]                         # indices 0..7  (paper's l_1..l_8)
eX = parse('X0**6 + X1**6 + X2**6'); eY = parse('Y0**6 + Y1**6 + Y2**6')

def pw(P_, n): return power(P_, n)
def D_I(I):
    out = {}
    for k in I:
        t = pw(ell[k], 6)
        if k >= 6: t = scale(t, c)
        out = add(out, t)
    return out
def Fsub(I):
    F, A, q, rec = load()
    out = {}
    for e, v in F.items():
        if all(e[i] == 0 for i in range(6) if i not in I): out[e] = v
    return out

def parity_parts(D):
    Dp = {e: v for e, v in D.items() if core.yparity(e) == 0}
    Dm = {e: v for e, v in D.items() if core.yparity(e) == 1}
    return Dp, Dm

def families(which):
    cands = []  # (name, family, D)
    idx6 = range(6)
    if which in ('deg6', 'all'):
        cands.append(('1', 'constant', parse('1')))
        for r in (2, 3, 4):
            for I in itertools.combinations(idx6, r):
                cands.append(('pow' + ''.join(map(str, I)), 'coordinate power sums (sweep)', D_I(I)))
        for r in (2, 3, 4, 5):
            for I in itertools.combinations(idx6, r):
                cands.append(('Fsub' + ''.join(map(str, I)), 'sub-equations of F (sweep)', Fsub(I)))
        for r in (2, 3, 4):
            for I in itertools.combinations(range(8), r):
                if 6 in I or 7 in I:
                    cands.append(('DI' + ''.join(map(str, I)), 'eight-form subset sums D_I with l7/l8', D_I(I)))
        for nm, D in [('eX', eX), ('eY', eY), ('eX+eY', add(eX, eY)), ('B=(L+M)^6+(L-M)^6', add(pw(l7, 6), pw(l8, 6))),
                      ('eX+2cL6', add(eX, scale(pw(L, 6), 2 * c))), ('eY+2cM6', add(eY, scale(pw(M, 6), 2 * c))),
                      ('eX-L6', sub(eX, pw(L, 6))), ('eX+L6', add(eX, pw(L, 6))), ('eY-M6', sub(eY, pw(M, 6))), ('eY+M6', add(eY, pw(M, 6))),
                      ('eX+M6', add(eX, pw(M, 6))), ('eY+L6', add(eY, pw(L, 6))), ('L6', pw(L, 6)), ('M6', pw(M, 6)),
                      ('L4M2', mul(pw(L, 4), pw(M, 2))), ('L2M4', mul(pw(L, 2), pw(M, 4))), ('eX+L4M2', add(eX, mul(pw(L, 4), pw(M, 2)))),
                      ('(L2+M2)^3', pw(add(pw(L, 2), pw(M, 2)), 3)), ('(L2-M2)^3', pw(sub(pw(L, 2), pw(M, 2)), 3)),
                      ('L2M2(L2+M2)', mul(mul(pw(L, 2), pw(M, 2)), add(pw(L, 2), pw(M, 2)))),
                      ('X0X1X2Y0Y1Y2', mul(mul(mul(Z[0], Z[1]), mul(Z[2], Z[3])), mul(Z[4], Z[5]))),
                      ('prod_l1..l6', mul(mul(mul(Z[0], Z[1]), mul(Z[2], Z[3])), mul(Z[4], Z[5]))),
                      ('(L+M)^5(L-M)', mul(pw(l7, 5), l8)), ('(L+M)^3(L-M)^3', mul(pw(l7, 3), pw(l8, 3)))]:
            cands.append((nm, 'degree-6 special', D))
    if which in ('deg7', 'all'):
        # l_j * D_I, |I| = 2..4, all eight l_j ;  products of seven of the eight l's
        for r in (2, 3, 4):
            for I in itertools.combinations(range(8), r):
                DI = D_I(I)
                for j in range(8):
                    cands.append((f'l{j}*DI' + ''.join(map(str, I)), 'l_j * D_I (deg 7)', mul(ell[j], DI)))
        for j in range(8):
            Pm = parse('1')
            for k in range(8):
                if k != j: Pm = mul(Pm, ell[k])
            cands.append((f'prod_l_except{j}', 'seven-form product (deg 7)', Pm))
        for nm, D in [('X0*eX', mul(Z[0], eX)), ('Y0*eY', mul(Z[3], eY)), ('L*(eX+2cL6)', mul(L, add(eX, scale(pw(L, 6), 2 * c)))),
                      ('M*(eY+2cM6)', mul(M, add(eY, scale(pw(M, 6), 2 * c)))), ('L*eY', mul(L, eY)), ('M*eX', mul(M, eX)),
                      ('L*B', mul(L, add(pw(l7, 6), pw(l8, 6)))), ('M*B', mul(M, add(pw(l7, 6), pw(l8, 6)))),
                      ('L*(eX+eY)', mul(L, add(eX, eY))), ('M*(eX+eY)', mul(M, add(eX, eY))),
                      ('X0*(eX+eY)', mul(Z[0], add(eX, eY))), ('Y0*(eX+eY)', mul(Z[3], add(eX, eY))),
                      ('L^7', pw(L, 7)), ('M^7', pw(M, 7)), ('(L+M)^7', pw(l7, 7)), ('X0^7', pw(Z[0], 7)), ('Y0^7', pw(Z[3], 7))]:
            cands.append((nm, 'degree-7 special', D))
    if which in ('deg8', 'all'):
        allprod = parse('1')
        for k in range(8): allprod = mul(allprod, ell[k])
        cands.append(('prod_all8', 'eight-form product (deg 8)', allprod))
        B2 = add(pw(L, 2), pw(M, 2)); B2m = sub(pw(L, 2), pw(M, 2))
        for r in (2, 3, 4):
            for I in itertools.combinations(range(8), r):
                DI = D_I(I)
                for nm, Qd in [('(L2+M2)', B2), ('(L2-M2)', B2m), ('L2', pw(L, 2)), ('M2', pw(M, 2)), ('LM', mul(L, M)),
                               ('X0^2', pw(Z[0], 2)), ('Y0^2', pw(Z[3], 2)), ('X0Y0', mul(Z[0], Z[3]))]:
                    cands.append((f'{nm}*DI' + ''.join(map(str, I)), 'quadric * D_I (deg 8)', mul(Qd, DI)))
        for j, k in itertools.combinations(range(8), 2):
            for I in [(0, 1), (0, 3), (3, 4), (0, 1, 2), (3, 4, 5), (0, 1, 3, 4), (6, 7), (0, 6, 7), (0, 1, 6, 7), (0, 3, 6, 7)]:
                cands.append((f'l{j}l{k}*DI' + ''.join(map(str, I)), 'l_j l_k * D_I (deg 8, selected I)', mul(mul(ell[j], ell[k]), D_I(I))))
        for nm, D in [('eX*X0^2', mul(eX, pw(Z[0], 2))), ('eX*X0X1', mul(eX, mul(Z[0], Z[1]))), ('eY*Y0^2', mul(eY, pw(Z[3], 2))),
                      ('eX*Y0^2', mul(eX, pw(Z[3], 2))), ('eY*X0^2', mul(eY, pw(Z[0], 2))), ('(eX+eY)(L2+M2)', mul(add(eX, eY), B2)),
                      ('B*(L2+M2)', mul(add(pw(l7, 6), pw(l8, 6)), B2)), ('(eX+2cL6)(L2+M2)', mul(add(eX, scale(pw(L, 6), 2 * c)), B2)),
                      ('(eX+2cL6)*M2', mul(add(eX, scale(pw(L, 6), 2 * c)), pw(M, 2))), ('(eY+2cM6)*L2', mul(add(eY, scale(pw(M, 6), 2 * c)), pw(L, 2))),
                      ('(L2+M2)^4', pw(B2, 4)), ('(L2-M2)^4', pw(B2m, 4)), ('L^8', pw(L, 8)), ('(L+M)^8', pw(l7, 8)), ('X0^8', pw(Z[0], 8))]:
            cands.append((nm, 'degree-8 special', D))
    if which in ('deg10', 'all'):
        B4 = add(add(pw(L, 4), scale(mul(pw(L, 2), pw(M, 2)), 14)), pw(M, 4))
        for r in (2, 3, 4):
            for I in itertools.combinations(range(8), r):
                cands.append((f'B4*DI' + ''.join(map(str, I)), '(L4+14L2M2+M4) * D_I (deg 10)', mul(B4, D_I(I))))
        for nm, D in [('B4*eX', mul(B4, eX)), ('B4*eY', mul(B4, eY)), ('B4*(eX+eY)', mul(B4, add(eX, eY))), ('B4*B', mul(B4, add(pw(l7, 6), pw(l8, 6)))),
                      ('L4*(eX+2cL6)', mul(pw(L, 4), add(eX, scale(pw(L, 6), 2 * c)))), ('M4*(eY+2cM6)', mul(pw(M, 4), add(eY, scale(pw(M, 6), 2 * c)))),
                      ('L2M2*(eX+eY)', mul(mul(pw(L, 2), pw(M, 2)), add(eX, eY)))]:
            cands.append((nm, 'degree-10 special', D))
    if which in ('deg12', 'all'):
        pieces = [('eX', eX), ('eY', eY), ('eX+eY', add(eX, eY)), ('B', add(pw(l7, 6), pw(l8, 6))), ('eX+2cL6', add(eX, scale(pw(L, 6), 2 * c))),
                  ('eY+2cM6', add(eY, scale(pw(M, 6), 2 * c))), ('eX-L6', sub(eX, pw(L, 6))), ('eY-M6', sub(eY, pw(M, 6))),
                  ('L6', pw(L, 6)), ('M6', pw(M, 6)), ('(L+M)^6', pw(l7, 6)), ('(L-M)^6', pw(l8, 6))]
        for (n1, P1), (n2, P2) in itertools.combinations_with_replacement(pieces, 2):
            cands.append((f'{n1}*{n2}', 'products of degree-6 blocks (deg 12)', mul(P1, P2)))
        subs = [I for r in (2, 3, 4) for I in itertools.combinations(range(8), r)]
        for I in subs:
            # D_I * tau(D_I) for mixed I (exactly one of 6,7) ; D_I * D_J for disjoint J (coordinate subsets)
            if (6 in I) != (7 in I):
                It = tuple(sorted((7 if k == 6 else 6 if k == 7 else k) for k in I))
                cands.append((f'DI{"".join(map(str,I))}*DI{"".join(map(str,It))}', 'D_I * tau(D_I), mixed I (deg 12)', mul(D_I(I), D_I(It))))
        for I in subs:
            for Jj in subs:
                if I < Jj and not set(I) & set(Jj) and len(I) + len(Jj) <= 6:
                    cands.append((f'DI{"".join(map(str,I))}*DI{"".join(map(str,Jj))}', 'D_I * D_J, disjoint (deg 12)', mul(D_I(I), D_I(Jj))))
    return cands

def write_sing(path, cands):
    F, A, q, rec = load()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    with open(path, 'w') as f:
        f.write(f'ring R = {p}, (X0,X1,X2,Y0,Y1,Y2), dp;\n')
        f.write('poly F = ' + to_str(Fp).replace('**', '^') + ';\n')
        f.write('poly A = ' + to_str(Ap).replace('**', '^') + ';\n')
        f.write('ideal J = jacob(F);\nint t0 = timer;\nideal G = std(J);\nprint("std seconds = " + string(timer - t0));\n')
        f.write('print("dim = " + string(dim(G)));\nprint("vdim = " + string(vdim(G)));\nintvec hv = hilb(G, 2);\nprint("hilb2 = " + string(hv));\n')
        f.write('poly Dp; poly Dm; int okp; int okm;\n')
        for nm, fam, D in cands:
            Dp, Dm = parity_parts(core.reduce_mod(D, p))
            sp = to_str(Dp).replace('**', '^') if Dp else '0'
            sm = to_str(Dm).replace('**', '^') if Dm else '0'
            f.write(f'Dp = {sp}; Dm = {sm};\n')
            f.write('okp = (reduce(A*Dp, G) == 0); okm = (reduce(A*Dm, G) == 0);\n')
            f.write(f'print("RES|{nm}|{fam}|{degree(D)}|" + string(okp) + "|" + string(okm) + "|" + string(size(Dp)) + "|" + string(size(Dm)));\n')
        f.write('quit;\n')

if __name__ == '__main__':
    which = sys.argv[1]
    cands = families(which)
    path = f'{PK}/e1screen_{which}.sing'
    write_sing(path, cands)
    print(f'{len(cands)} candidates -> {path}', flush=True)
    if '--run' in sys.argv:
        t0 = time.time()
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=72000)
        out = r.stdout + r.stderr
        with open(f'{PK}/e1screen_{which}.out', 'w') as f: f.write(out)
        hilb = next((l.split('=', 1)[1].strip() for l in out.splitlines() if l.startswith('hilb2 =')), None)
        expected = ','.join(str(hilb_R(k)) for k in range(25))
        hilb_ok = (hilb == expected)
        res = []
        for l in out.splitlines():
            if l.startswith('RES|'):
                _, nm, fam, d, okp, okm, np_, nm_ = l.split('|')
                res.append({'name': nm, 'family': fam, 'd': int(d), 'E1_even_part': okp == '1', 'E1_odd_part': okm == '1',
                            'E1': okp == '1' and okm == '1', 'terms+': int(np_), 'terms-': int(nm_)})
        passing = [r_['name'] for r_ in res if r_['E1']]
        summary = {'exp': 'e1screen', 'set': which, 'p': p, 'n_candidates': len(cands), 'n_results': len(res), 'hilbert_ok': hilb_ok,
                   'hilb2': hilb, 'E1_passing': passing, 'seconds': round(time.time() - t0, 1)}
        with open(f'{PK}/e1screen_{which}.json', 'w') as f: json.dump({'summary': summary, 'results': res}, f, indent=1)
        log(summary)
        print(json.dumps(summary), flush=True)
