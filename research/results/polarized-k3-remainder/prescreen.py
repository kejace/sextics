"""For remaining low-degree candidates: instant test D in J(F)_d (mod p), equivalent to E1 for d<=9 by the
maximal-rank theorem; direct engine E1 as well when d<=4 (cheap); full lean solve only for members of J(F)."""
import sys, json, time; sys.path.insert(0,'/home/kejace/devel/sextics-work'); sys.path.insert(0,'/tmp/hodge/results/polarized-k3-remainder')
from hodge import core, fixed_d
from flint import nmod_mat
from logutil import log
from runD2 import uv_to_xy
from fastsolve import solve_lean
F,A,q,rec = core.load_instance('polarized-k3-remainder'); p = core.PRIMES[0]
Fp = core.reduce_mod(F,p); Fi = [core.pderiv(Fp,i,p) for i in range(6)]
done = {json.loads(l)['D'] for l in open('log.jsonl') if l.strip()}
def in_J(Dp, d):
    if d < 5: return False, 0
    par = core.parity_of(Dp)
    cols = []
    for i in range(6):
        pi = (par + (1 if i >= 3 else 0)) % 2  # d_i F has parity of F (even) shifted: x-partials even, y-partials odd
        for m in core.monomials(d-5, (par + (0 if i < 3 else 1)) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    rows = {}
    for c in cols + [Dp]:
        for e in c: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols); flat = [0]*(R*(C+1))
    for j,c in enumerate(cols):
        for e,v in c.items(): flat[rows[e]*(C+1)+j] = v
    for e,v in Dp.items(): flat[rows[e]*(C+1)+C] = v
    Ma = nmod_mat(R, C+1, flat, p); M0 = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    r0 = M0.rank(); return Ma.rank() == r0, r0
for fam, label, expr, frame in json.load(open('candidates_low.json')):
    if label in done: continue
    D = core.parse(expr); D = uv_to_xy(D) if frame == 'uv' else D
    d = core.degree(D); par = core.parity_of(D)
    if par is None:
        log(label, d, fam, None, None, -1, -1, 0.0, f'mixed parity; skipped; frame={frame}; expr={expr}'); continue
    t0 = time.time(); Dp = core.reduce_mod(D, p)
    member, rJ = in_J(Dp, d)
    note = f'J-membership test: D {"in" if member else "NOT in"} J(F)_{d} (rank {rJ}); frame={frame}; expr={expr}'
    e1 = None
    if d <= 4:
        ok, info = fixed_d.e1_solvable(F, A, q, D); e1 = ok; note += f'; direct E1={ok} (rank {info["rank"]}, aug {info["rank_aug"]})'
    if not member:
        log(label, d, fam, False if e1 is None else e1, False, -1 if e1 is None else info['rank'], -1, time.time()-t0, 'E1 fails; ' + note)
        print(label, 'E1-dead', note[:80], flush=True); continue
    out = solve_lean(F, A, q, D, verbose=True)
    log(label, d, fam, True, out['consistent'], out['rank'], out['nullity'], time.time()-t0, f'lean; rows={out["rows"]} cols={out["cols"]} k={out["k"]}; ' + note)
    print(label, out['consistent'], out['rank'], out['nullity'], flush=True)
print('PRESCREEN-DONE', flush=True)
