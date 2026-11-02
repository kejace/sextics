"""Parse msolve output (-P 1, char p) for a zero-dimensional system and recover F_p-rational solutions, validated
against the equations in the input .ms file (so msolve's sign/denominator conventions do not matter)."""
import ast, re, sys, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from flint import nmod_poly

def parse_out(path):
    txt = open(path).read().strip().rstrip(':')
    data = ast.literal_eval(txt)
    if data == [-1]:
        return {'status': 'no-solutions'}
    if data[0] == 1:
        return {'status': 'positive-dimensional', 'nvars': data[1]}
    if data[0] != 0:
        return {'status': f'unknown {data[0]}'}
    body = data[1]
    p, nvars, deg, vars_, lin = body[0], body[1], body[2], body[3], body[4]
    param = body[5]
    elim = param[1][0]; den = param[1][1]; pars = param[1][2]
    def unwrap(x):
        while isinstance(x, list) and len(x) == 1 and isinstance(x[0], list):
            x = x[0]
        return x
    plist = []
    for pp in pars:
        pp = unwrap(pp)
        if isinstance(pp, list) and len(pp) == 2 and isinstance(pp[1], list):
            plist.append(pp[1])
        elif isinstance(pp, list) and len(pp) == 2 and isinstance(pp[0], list):
            plist.append(unwrap(pp[0])[1])
        else:
            plist.append(pp)
    return {'status': 'zero-dimensional', 'p': p, 'nvars': nvars, 'degree': deg, 'vars': vars_, 'linform': lin,
            'elim': elim[1], 'den': den[1], 'params': plist}

def parse_ms(path):
    """Equations of an msolve input file as list of dict {(i,j)|(i,)|(): coeff}; variables a1..an -> 0..n-1."""
    lines = open(path).read().split('\n')
    vars_ = lines[0].split(','); p = int(lines[1])
    vidx = {v: i for i, v in enumerate(vars_)}
    eqs = []
    for line in lines[2:]:
        line = line.strip().rstrip(',')
        if not line: continue
        eq = {}
        for term in re.findall(r'[+-]?[^+-]+', line.replace(' ', '')):
            parts = term.split('*')
            coef = 1; vs = []
            for part in parts:
                if re.fullmatch(r'[+-]?\d+', part): coef = int(part)
                else:
                    m = re.fullmatch(r'([+-]?)([A-Za-z_]\w*)(?:\^(\d+))?', part)
                    sign, v, e = m.group(1), m.group(2), m.group(3)
                    if sign == '-': coef = -coef
                    vs += [vidx[v]] * (int(e) if e else 1)
            key = tuple(sorted(vs))
            eq[key] = (eq.get(key, 0) + coef) % p
        eqs.append(eq)
    return vars_, p, eqs

def evaluate(eq, x, p):
    s = 0
    for key, c in eq.items():
        t = c
        for i in key: t = t * x[i] % p
        s = (s + t) % p
    return s

def solutions(ms_path, out_path):
    vars_, p, eqs = parse_ms(ms_path)
    info = parse_out(out_path)
    if info['status'] != 'zero-dimensional':
        return info, []
    el = nmod_poly(info['elim'], p)
    roots = [int(r) for r, mult in el.roots()]
    den = nmod_poly(info['den'], p)
    dden = den.derivative() if False else None
    elp = el.derivative()
    vidx = {v: i for i, v in enumerate(vars_)}
    # msolve: variables listed; the eliminating linear form is the last variable (or a random linear form)
    lin = info['linform']
    pts = []
    for t in roots:
        cands = []
        for denom in (den(t), elp(t)):
            dv = int(denom)
            if dv == 0: continue
            inv = pow(dv, -1, p)
            for sgn in (1, p - 1):
                vals = [int(nmod_poly(par, p)(t)) * inv * sgn % p for par in info['params']]
                cands.append(vals)
        # assemble full vectors by variable NAME: msolve lists its own variable order; the eliminating variable is
        # the last one of that order (when the linear form is a coordinate), params give the others in order.
        ovars = info['vars']
        ok = False
        for vals in cands:
            byname = {}
            if len(vals) == len(ovars) - 1:
                for v, val in zip(ovars[:-1], vals): byname[v] = val
                byname[ovars[-1]] = t
            elif len(vals) == len(ovars):
                for v, val in zip(ovars, vals): byname[v] = val
            else:
                continue
            x = [byname[v] for v in vars_]
            if all(evaluate(eq, x, p) == 0 for eq in eqs):
                pts.append({'t': t, 'x': x}); ok = True; break
        if not ok:
            pts.append({'t': t, 'x': None, 'note': 'could not validate parametrization conventions'})
    return info, pts

if __name__ == '__main__':
    info, pts = solutions(sys.argv[1], sys.argv[2])
    print({k: v for k, v in info.items() if k not in ('elim', 'den', 'params')})
    print('F_p-rational points:', len(pts))
    for pt in pts[:10]: print(pt)
