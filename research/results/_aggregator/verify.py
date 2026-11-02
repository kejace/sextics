#!/usr/bin/env python3
"""Aggregator spot-check runner (wave 2).

usage: python3 verify.py <checks.json> <out.jsonl> [--workers N] [--timeout S] [--only n1,n2]

Each check is a JSON object {"name", "lane", "kind", ...}; one JSON result line per check is
appended to <out.jsonl> as soon as it finishes (also on error/timeout).  Kinds:

  e1kernel   {"d", "parity", "p"}            -> jacobian.mult_kernel: dim (kernel mod J), kernel_dim, jacobian_dim, ring_dim
  hilbert    {"d", "p"}                      -> (R/J)_d and J_d dims for both parities
  a_in_j     {"p"}                           -> jacobian.in_jacobian(F, A)
  e1         {"D", "p", ["target"]}          -> fixed_d.e1_solvable
  fixedD     {"D", "p", "method": dense|compressed|exact, ["target"], ["F"]} -> consistency/rank/nullity
  diag_group {}                              -> symmetry.diagonal_group(F): order, size
  poly_eq    {"lhs", "rhs"}                  -> exact equality of two polynomial specs

Polynomial specs: a string expression, {"expr": s}, {"partial": i}, {"euler": [i,..]} (sum Z_i F_i),
{"file": path, "key": k, "index": j}, {"product": [spec, spec]}, {"kernel_rep": {"d","parity","index"}}
(a basis element of the E1 kernel mod J computed mod p; only meaningful with the same p).
"""
import sys, os, json, time, traceback
os.environ.setdefault('FLINT_NUM_THREADS', '2')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import multiprocessing as mp

NV = 6


def build_poly(spec, F, A, q, p):
    from hodge import core, jacobian
    if isinstance(spec, str):
        return core.parse(spec)
    if 'expr' in spec:
        return core.parse(spec['expr'])
    if 'partial' in spec:
        return core.pderiv(F, spec['partial'])
    if 'euler' in spec:
        D = {}
        for i in spec['euler']:
            e = [0] * NV; e[i] = 1
            D = core.add(D, core.shift(core.pderiv(F, i), tuple(e)))
        return D
    if 'file' in spec:
        path = spec['file']
        if path.endswith('.json'):
            obj = json.load(open(path))
            if isinstance(obj, dict):
                obj = obj[spec.get('key', 'D')]
            if isinstance(obj, list):
                obj = obj[spec.get('index', 0)]
        else:
            obj = open(path).read()
        return core.parse(obj)
    if 'product' in spec:
        P = build_poly(spec['product'][0], F, A, q, p)
        for s in spec['product'][1:]:
            P = core.mul(P, build_poly(s, F, A, q, p))
        return P
    if 'lincomb' in spec:
        from fractions import Fraction
        P = {}
        for c, s in spec['lincomb']:
            P = core.add(P, core.scale(build_poly(s, F, A, q, p), Fraction(c)))
        return P
    if 'kernel_rep' in spec:
        kr = spec['kernel_rep']
        out = jacobian.mult_kernel(F, A, kr['d'], p, parity=kr.get('parity', 0))
        if kr['index'] >= len(out['basis']):
            raise IndexError(f"kernel basis has only {len(out['basis'])} elements (dim {out['dim']})")
        return out['basis'][kr['index']]
    raise ValueError(f'bad poly spec {spec}')


def run_check(chk):
    from hodge import core, fixed_d, jacobian, symmetry
    t0 = time.time()
    lane = chk['lane']
    p = chk.get('p', 1000000007)
    F, A, q, rec = core.load_instance(lane)
    if 'F' in chk:
        F = build_poly(chk['F'], F, A, q, p)
    if 'target' in chk:
        A = build_poly(chk['target'], F, A, q, p)
    if 'q' in chk:
        q = chk['q']
    kind = chk['kind']
    res = {}
    if kind == 'e1kernel':
        out = jacobian.mult_kernel(F, A, chk['d'], p, parity=chk.get('parity', 0))
        res = {k: out[k] for k in ('dim', 'kernel_dim', 'jacobian_dim', 'ring_dim', 'rows', 'cols')}
        if chk.get('save_basis'):
            res['basis'] = [core.to_str(b) for b in out['basis']]
    elif kind == 'hilbert':
        for par in (0, 1):
            std, red = jacobian.jacobian_ring_basis(F, chk['d'], p, par)
            res[f'ring_dim_par{par}'] = red.ring_dim
            res[f'jacobian_dim_par{par}'] = red.jacobian_dim
    elif kind == 'a_in_j':
        ok, info = jacobian.in_jacobian(F, A, p)
        res = {'in_jacobian': ok, **info}
    elif kind == 'e1':
        D = build_poly(chk['D'], F, A, q, p)
        ok, info = fixed_d.e1_solvable(F, A, q, D, p)
        res = {'e1': ok, 'd': core.degree(D), 'terms': len(D), **info}
    elif kind == 'fixedD':
        D = build_poly(chk['D'], F, A, q, p)
        method = chk.get('method', 'dense')
        d = core.degree(D)
        res['d'] = d; res['terms'] = len(D)
        if method == 'dense':
            S = fixed_d.FixedDSystem(F, A, q, D, p)
            r = S.solve()
            res.update({k: r[k] for k in ('consistent', 'rank', 'nullity', 'rows', 'cols', 't_build', 't_rref')})
        elif method == 'compressed':
            S = fixed_d.FixedDSystem(F, A, q, D, p)
            r = fixed_d.solve_compressed(S, k_extra=chk.get('k_extra', 40), seed=chk.get('seed', 1))
            res.update({k: r[k] for k in ('consistent', 'rank', 'nullity', 'rows', 'cols', 'k', 't_compress', 't_rref')})
        elif method == 'exact':
            status, data = fixed_d.solve_fixed_D(F, A, q, D, primes=[p] + [pp for pp in core.PRIMES if pp != p], max_primes=chk.get('max_primes', 6))
            res['status'] = status
            res.update({k: v for k, v in data.items() if k not in ('N', 'Q')})
            if status == 'solved':
                res['N'] = [core.to_str(n) for n in data['N']]
                res['Q'] = core.to_str(data['Q'])
        else:
            raise ValueError(method)
    elif kind == 'diag_group':
        G = symmetry.diagonal_group(F)
        res = {'order': G.order, 'size': len(G), 'generators': [list(g) for g in G.generators]}
    elif kind == 'poly_eq':
        L = build_poly(chk['lhs'], F, A, q, p); R = build_poly(chk['rhs'], F, A, q, p)
        diff = core.sub(L, R)
        res = {'equal': not diff, 'diff_terms': len(diff)}
    else:
        raise ValueError(kind)
    res['seconds'] = round(time.time() - t0, 1)
    return res


def worker(chk, out_path):
    line = {'name': chk['name'], 'lane': chk['lane'], 'kind': chk['kind'], 'p': chk.get('p', 1000000007),
            'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'spec': {k: v for k, v in chk.items() if k not in ('name', 'lane', 'kind')}}
    try:
        line['result'] = run_check(chk)
        line['ok'] = True
        if 'expect' in chk:
            exp = chk['expect']
            got = {k: line['result'].get(k) for k in exp}
            line['agree'] = all(line['result'].get(k) == v for k, v in exp.items())
            line['expect'] = exp; line['got'] = got
    except Exception as ex:
        line['ok'] = False
        line['error'] = f'{type(ex).__name__}: {ex}'
        line['trace'] = traceback.format_exc()[-1500:]
    with open(out_path, 'a') as f:
        f.write(json.dumps(line, default=str) + '\n')
        f.flush()


def main():
    args = sys.argv[1:]
    checks_path, out_path = args[0], args[1]
    workers = 8; timeout = 3600.0; only = None
    i = 2
    while i < len(args):
        if args[i] == '--workers': workers = int(args[i + 1]); i += 2
        elif args[i] == '--timeout': timeout = float(args[i + 1]); i += 2
        elif args[i] == '--only': only = set(args[i + 1].split(',')); i += 2
        else: i += 1
    checks = json.load(open(checks_path))
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            try: done.add(json.loads(line)['name'])
            except Exception: pass
    todo = [c for c in checks if c['name'] not in done and (only is None or c['name'] in only)]
    print(f'{len(todo)} checks to run ({len(done)} already done), {workers} workers', flush=True)
    running = []  # (proc, chk, t0)
    while todo or running:
        while todo and len(running) < workers:
            chk = todo.pop(0)
            pr = mp.Process(target=worker, args=(chk, out_path))
            pr.start()
            running.append((pr, chk, time.time()))
            print(f'[start] {chk["name"]}', flush=True)
        time.sleep(2)
        still = []
        for pr, chk, t0 in running:
            if pr.is_alive():
                if time.time() - t0 > timeout:
                    pr.terminate(); pr.join(5)
                    with open(out_path, 'a') as f:
                        f.write(json.dumps({'name': chk['name'], 'lane': chk['lane'], 'kind': chk['kind'], 'ok': False,
                                            'error': f'timeout after {timeout}s', 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}) + '\n')
                    print(f'[timeout] {chk["name"]}', flush=True)
                else:
                    still.append((pr, chk, t0))
            else:
                print(f'[done] {chk["name"]} ({time.time()-t0:.0f}s, exit {pr.exitcode})', flush=True)
                if pr.exitcode not in (0, None):
                    with open(out_path, 'a') as f:
                        f.write(json.dumps({'name': chk['name'], 'lane': chk['lane'], 'kind': chk['kind'], 'ok': False,
                                            'error': f'worker exit code {pr.exitcode}', 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}) + '\n')
        running = still
    print('all done', flush=True)


if __name__ == '__main__':
    main()
