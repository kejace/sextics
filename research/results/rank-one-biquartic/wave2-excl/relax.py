"""Hessian relaxation on the F-reduced E1 family of degree d / parity (0, 1, or m = mixed):
quadratic system Phi(a) = [B_{D(a)}(n(a))] in R_{2d}; Groebner basis (Singular std) mod p; dim S/I = 0 means V = {0}.
Usage: python3 relax.py d parity prime [family]   family in {all, bihom:a:b, base, xshape, kY2, kX2}"""
import sys, json, time, subprocess, os
sys.path.insert(0, '.')
import w2
d = int(sys.argv[1]); par = sys.argv[2]; p = int(sys.argv[3]); famname = sys.argv[4] if len(sys.argv) > 4 else 'all'
pD = None if par == 'm' else int(par)
R, Fp, Ap, Fi, F, A = w2.load(p)
t0 = time.time()
pred = w2.not_div_X06
if famname.startswith('bihom'):
    _, a_, b_ = famname.split(':'); a_, b_ = int(a_), int(b_)
    pred = lambda e: (e[0] + e[1] + e[2] == a_) and (e[3] + e[4] + e[5] == b_)
fam = w2.e1_family(R, Fp, Ap, Fi, d, pD, dpred=pred)
m = fam['dimK']
rJ, rJV = w2.jacobian_dim(R, Fi, d, pD, dpred=pred)
print(f'family {famname} d={d} par={par} p={p}: dim K cap V = {m}, dim J_d = {rJ}, dim J_d cap V = {rJV}', flush=True)
rec = dict(exp='relaxation', d=d, parity=par, p=p, family=famname, m=m, dimJ_d=rJ, dimJ_cap_V=rJV,
           e1_rows=fam['rows'], e1_cols=fam['cols'], nD=fam['nD'])
if m == 0:
    rec.update(result='empty family', seconds=round(time.time() - t0, 1)); w2.log(rec); print(rec); sys.exit(0)
quad, red = w2.relaxation_quadrics(R, Fi, fam, d, pD)
tag = f'relax_d{d}_{par}_{famname.replace(":", "")}_{p}'
nq = w2.write_singular(tag + '.sing', quad, m, p)
w2.write_msolve(tag + '.ms', quad, m, p)
json.dump({'p': p, 'd': d, 'parity': par, 'family': famname, 'basis': [[str(D), [str(x) for x in n]] for D, n in fam['basis']],
           'quadrics': [{f'{k},{l}': v for (k, l), v in q.items()} for q in quad]}, open(tag + '.json', 'w'))
rec.update(n_quadrics=nq, ring_dim_2d=red.ring_dim, rank_J_2d=red.rank, t_setup=round(time.time() - t0, 1))
print('running Singular on', tag, flush=True)
t1 = time.time()
try:
    r = subprocess.run(['Singular', '-q', tag + '.sing'], capture_output=True, text=True, timeout=int(os.environ.get('SING_TIMEOUT', '7200')))
    outp = r.stdout
except subprocess.TimeoutExpired:
    outp = 'TIMEOUT'
open(tag + '.sing.out', 'w').write(outp)
print(outp, flush=True)
dim = None
for line in outp.splitlines():
    if line.startswith('dimension of S/I'):
        dim = int(line.split('=')[1])
rec.update(singular_dim=dim, singular_out=outp.strip()[:500], t_groebner=round(time.time() - t1, 1),
           result=('V(Phi) = {0}: every nonzero F-reduced E1-admissible D of this family fails E2 (rigorous mod-p argument)' if dim == 0
                   else 'V(Phi) nontrivial or not computed'))
w2.log(rec)
print(json.dumps({k: v for k, v in rec.items() if k != 'singular_out'}), flush=True)
