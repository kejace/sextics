"""(b) Groebner analysis of a relaxation system: msolve (grevlex, mod p) -> leading ideal -> Hilbert series via Singular.
Rigor argument (see README): if the Hilbert polynomial of F_p[v]/I_p equals that of the trivial family L (a linear
subspace P^r, HP = C(k+r, r)), then over Q (hence C) the relaxation variety is exactly L, so no denominator of this
degree/parity-mode is coprime to F.   Usage: python3 b2_gb.py <tag> [threads] [expected_r]"""
import sys, json, subprocess, time, re, os; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import log
tag = sys.argv[1]; thr = int(sys.argv[2]) if len(sys.argv) > 2 else 4
meta = json.load(open(tag + '.meta.json'))
n = meta['n_unknowns']; d = meta['d']; mode = meta['mode']; p = meta['p']
# expected trivial family: Y-even part w_i = l z_i, l in R_{d-6} => projective dim = dim R_{d-6} - 1 (only if an even part exists)
from math import comb
r_expected = comb(d - 6 + 5, 5) - 1 if (mode in ('even', 'mixed') and d >= 6) else -1
if len(sys.argv) > 3: r_expected = int(sys.argv[3])
t0 = time.time()
env = dict(os.environ, OMP_NUM_THREADS=str(thr))
res = subprocess.run(['msolve', '-g', '2', '-t', str(thr), '-f', tag + '.ms', '-o', tag + '.gb'], capture_output=True, text=True, env=env, timeout=36000)
t_gb = round(time.time() - t0, 1)
txt = open(tag + '.gb').read()
body = txt[txt.index('['):txt.rindex(']')].strip('[]')
polys = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
vars_ = open(tag + '.ms').readline().strip().split(',')
vidx = {v: i for i, v in enumerate(vars_)}
def parse_poly(s):
    terms = []
    for t in re.split(r'(?=[+-])', s):
        t = t.strip()
        if not t: continue
        sign = 1
        if t[0] == '+': t = t[1:]
        elif t[0] == '-': sign = -1; t = t[1:]
        factors = t.split('*')
        coef = 1; exp = [0] * n
        for f in factors:
            if '^' in f:
                v, e = f.split('^'); exp[vidx[v]] += int(e)
            elif f in vidx: exp[vidx[f]] += 1
            else: coef = coef * int(f)
        terms.append((tuple(exp), sign * coef % p))
    return terms
def lead_grevlex(terms):
    def key(e):
        return (sum(e), tuple(-x for x in reversed(e)))
    return max((e for e, c in terms), key=key)
GB = [parse_poly(s) for s in polys]
leads = [lead_grevlex(t) for t in GB]
# minimalize the leading monomial ideal
def divides(a, b): return all(x <= y for x, y in zip(a, b))
mins = [m for m in leads if not any(o != m and divides(o, m) for o in leads)]
deg_hist = {}
for t in GB: deg_hist[sum(lead_grevlex(t))] = deg_hist.get(sum(lead_grevlex(t)), 0) + 1
n_linear = sum(1 for t in GB if sum(lead_grevlex(t)) == 1)
# Singular: Hilbert series of the leading ideal
def mono(e): return '*'.join(f'{vars_[i]}^{x}' for i, x in enumerate(e) if x) or '1'
sing = f'ring r = {p}, ({",".join(vars_)}), dp;\nideal L = {",".join(mono(m) for m in mins)};\nattrib(L, "isSB", 1);\n'
sing += 'int dm = dim(L); int mu = mult(L); intvec h2 = hilb(L, 2); intvec h1 = hilb(L, 1);\nprint("DIM " + string(dm)); print("MULT " + string(mu)); print("H2 " + string(h2)); print("H1 " + string(h1));\n$\n'
open(tag + '.hilb.sing', 'w').write(sing)
sres = subprocess.run(['Singular', '-q', tag + '.hilb.sing'], capture_output=True, text=True, timeout=3600)
out = sres.stdout
dm = int(re.search(r'DIM (-?\d+)', out).group(1)); mu = int(re.search(r'MULT (-?\d+)', out).group(1))
h2 = re.search(r'H2 (.*)', out).group(1).strip(); h1 = re.search(r'H1 (.*)', out).group(1).strip()
# HP equals that of P^r iff Krull dim = r+1 and second Hilbert series numerator is "1,0" (i.e. 1)
h2_list = [int(x) for x in h2.split(',')]
hp_is_linear_space = (dm == r_expected + 1) and (h2_list[:-1] == [1] or h2_list == [1, 0] or h2_list == [1])
ent = dict(exp='b2-groebner', tag=tag, d=d, mode=mode, p=p, n_unknowns=n, n_input_quadrics=meta['n_quadrics'], gb_size=len(GB),
           gb_lead_degree_hist=deg_hist, n_linear_in_gb=n_linear, lead_ideal_min_gens=len(mins), krull_dim_affine=dm, degree=mu,
           hilb_second=h2, hilb_first_len=len(h1.split(',')), r_expected_projective=r_expected,
           HP_equals_trivial_family=bool(hp_is_linear_space), msolve_seconds=t_gb, threads=thr,
           conclusion=('RIGOROUS: relaxation variety = trivial family over C => no coprime denominator of this degree in mode ' + mode) if hp_is_linear_space else 'NOT CONCLUSIVE: extra components, inspect GB')
log(ent)
json.dump(dict(ent, leads=[mono(m) for m in mins]), open(tag + '.gbinfo.json', 'w'), indent=1)
