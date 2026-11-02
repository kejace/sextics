"""Referee replay of the Hessian-relaxation verdicts with an independent builder (ref_lib.hess_relaxation).
usage: python3 ref_relax.py <mode>   modes: d5 | d6odd | d6even | d6mixed | d7empty | d7lift | d7odd0build"""
import sys, json, time, subprocess, os, re
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
from fractions import Fraction
mode = sys.argv[1]; p = L.P1
Fp, Fip, Hp = L.Fpolys(p); FQ, FiQ, HQ = L.Fpolys(None)
OUT = 'relax'; os.makedirs(OUT, exist_ok=True)
def sing(path, timeout=3000):
    t = time.time(); r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
    return r.stdout.strip() + (('\nSTDERR: ' + r.stderr.strip()[:300]) if r.stderr.strip() else ''), time.time() - t
def vname(u): i, m = u; return f'({i},{"".join(f"{n}{e}" for n, e in zip("XXXYYY", m) )})'
def run_dim(tag, d, pars, chars):
    s = L.hess_relaxation(d, pars, chars, Fip, Hp, p)
    path = f'{OUT}/{tag}_p{p}.sing'
    L.write_singular(s, path, p, 'ideal SI = std(I);\n"gb_size", size(SI), "dim", dim(SI), "vdim", vdim(SI);\nquit;\n')
    out, t = sing(path)
    print(f'{tag}: n_vars {s["n"]} (gauge removed {s["n_gauge"]}) n_eqs {s["n_eqs"]} | {out} [{t:.1f}s]', flush=True)
    return s, out
if mode == 'd5':
    for tag, pars in [('d5_even', [0]), ('d5_odd', [1]), ('d5_mixed', [0, 1])]: run_dim(tag, 5, pars, range(7))
elif mode == 'd6odd':
    run_dim('d6_odd_allchars', 6, [1], range(7))
elif mode == 'd7empty':
    for k in (0, 3, 5, 6): run_dim(f'd7_even_k{k}', 7, [0], [k])
    for k in range(1, 7): run_dim(f'd7_odd_k{k}', 7, [1], [k])
elif mode in ('d6even', 'd6mixed'):
    pars = [0] if mode == 'd6even' else [0, 1]
    s = L.hess_relaxation(6, pars, range(7), FiQ, HQ, None)
    vars_ = s['vars']; n = s['n']
    diag = [a for a, (i, m) in enumerate(vars_) if par(m) == PI[i] and sum(m) == 1 and m[i] == 1] if False else None
    diag = [a for a, (i, m) in enumerate(vars_) if m[i] == 1 and (L.par(m) + L.PI[i]) % 2 == 0]   # w_i = lambda_i Z_i (Y-even D)
    off = [a for a in range(n) if a not in diag]
    print(f'{mode}: n_vars {n} n_eqs {s["n_eqs"]}; diagonal unknowns {[(a, vname(vars_[a])) for a in diag]}', flush=True)
    extra = 'ideal SI = std(I);\n"gb_size", size(SI), "dim", dim(SI);\nint k; int bad = 0;\n'
    for a in off:
        extra += f'k = 1; while (k <= 8) {{ if (reduce(u{a}^k, SI) == 0) {{ break; }} k = k + 1; }} if (k > 8) {{ bad = bad + 1; "FAIL u{a}"; }} else {{ "u{a} power", k, "in I"; }}\n'
    # restriction to the diagonal subspace: substitute off-diagonals = 0, collect the quadrics (ideal of C in the diagonal coords)
    extra += 'ideal Irest = I;\n' + ''.join(f'Irest = subst(Irest, u{a}, 0);\n' for a in off) + 'ideal SC = std(Irest);\n"C_ideal_gb"; print(SC); "dim", dim(SC), "mult", mult(SC);\n'
    # compare with the lens's 3 exact C quadrics (l1..l5 -> diagonal unknowns in order i = 1..5)
    Cq = json.load(open('../local_d6/C_quadrics_exact.json'))
    dmap = {i: a for a, (i, m) in enumerate(vars_) if a in diag}   # i -> unknown index
    lens = []
    for qd in Cq:
        lens.append('+'.join(f'({v})*u{dmap[int(i)]}*u{dmap[int(j)]}' for ij, v in qd.items() for i, j in [ij.split(',')]))
    extra += 'ideal CL = ' + ',\n'.join(lens) + ';\nideal SCL = std(CL);\n"lens C quadrics reduce to 0 mod mine:", size(reduce(CL, SC)) == 0, "| mine reduce to 0 mod lens:", size(reduce(SC, SCL)) == 0;\n'
    extra += '"total failures", bad;\nquit;\n'
    path = f'{OUT}/{mode}_exactQ.sing'; L.write_singular(s, path, None, extra)
    out, t = sing(path, 7200); print(out[-3000:], f'[{t:.1f}s]', flush=True)
    json.dump({'vars': vars_, 'diag': diag, 'off': off, 'n_eqs': s['n_eqs']}, open(f'{OUT}/{mode}_vars.json', 'w'))
elif mode == 'd7lift':
    # exact over Q: V(I) subset of the lifted family X_j . C for k = w_j (j = 0,1,2), using MY C quadrics from d6even (read back from Singular)
    # C quadrics: take from the lens file but verified equal to mine in mode d6even (ideal equality printed there)
    Cq = [{tuple(int(x) for x in kk.split(',')): Fraction(v) for kk, v in qd.items()} for qd in json.load(open('../local_d6/C_quadrics_exact.json'))]
    for j, k in [(0, 1), (1, 2), (2, 4)]:
        s = L.hess_relaxation(7, [0], [k], FiQ, HQ, None); vars_ = s['vars']; n = s['n']
        Xj = tuple(1 if t == j else 0 for t in range(6))
        lam = {}   # i -> unknown index of (i, X_j Z_i)
        for a, (i, m) in enumerate(vars_):
            Zi = tuple(1 if t == i else 0 for t in range(6))
            if m == tuple(x + y for x, y in zip(Xj, Zi)): lam[i] = a
        vanish = [a for a in range(n) if a not in lam.values()]
        assert set(lam) == {1, 2, 3, 4, 5}, lam
        images = ['+'.join(f'({v})*u{lam[i]}*u{lam[ip]}' for (i, ip), v in sorted(qd.items())) for qd in Cq]
        extra = 'ideal SI = std(I);\n"gb_size", size(SI), "dim", dim(SI), "mult", mult(SI);\nint k; int bad = 0;\n'
        for a in vanish:
            extra += f'k = 1; while (k <= 8) {{ if (reduce(u{a}^k, SI) == 0) {{ break; }} k = k + 1; }} if (k > 8) {{ bad = bad + 1; "FAIL u{a}"; }} else {{ "u{a} power", k, "in I"; }}\n'
        for q in images:
            extra += f'k = 1; while (k <= 8) {{ if (reduce(({q})^k, SI) == 0) {{ break; }} k = k + 1; }} if (k > 8) {{ bad = bad + 1; "FAIL Cquadric"; }} else {{ "C-quadric image power", k, "in I"; }}\n'
        # lifting check: substitute the parametrization (vanishing coords = 0) into I and reduce modulo the C ideal
        extra += 'ideal Isub = I;\n' + ''.join(f'Isub = subst(Isub, u{a}, 0);\n' for a in vanish) + 'ideal CI = ' + ','.join(images) + ';\nideal SCI = std(CI);\n"lifted family inside V(I):", size(reduce(Isub, SCI)) == 0;\n'
        extra += '"total failures", bad;\nquit;\n'
        path = f'{OUT}/d7_even_k{k}_exactQ.sing'; L.write_singular(s, path, None, extra)
        out, t = sing(path, 7200)
        print(f'd7 even k={k} (lift X{j}.C): n_vars {n} n_eqs {s["n_eqs"]} | ' + out.replace('\n', ' | ')[-1500:], f'[{t:.1f}s]', flush=True)
elif mode == 'd7odd0build':
    s = L.hess_relaxation(7, [1], [0], Fip, Hp, p)
    L.write_msolve(s, f'{OUT}/d7_odd_k0_p{p}.ms', p)
    json.dump({'vars': s['vars'], 'n': s['n'], 'n_eqs': s['n_eqs']}, open(f'{OUT}/d7_odd_k0_vars.json', 'w'))
    print(f'd7 odd k0: n_vars {s["n"]} n_eqs {s["n_eqs"]} written to {OUT}/d7_odd_k0_p{p}.ms')
print('RELAX_DONE', mode)
