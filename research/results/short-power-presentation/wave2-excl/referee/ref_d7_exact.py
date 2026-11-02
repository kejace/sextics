"""Referee replay: exact verification of the degree-7 E1-kernel certificates D7^+ (even), D7^- (odd).
Checks (1) A*D = sum_i N_i dF/dZ_i exactly over Q (fmpq_mpoly); (2) D7^+ in Q[X0,X1,X2,M] via the derivations
(d/dY1 - 3 d/dY0), (d/dY2 - 5 d/dY0); mirror for D7^-; (3) mod p = 1073741741 (fresh): D7 not in J_7 (Hilbert-checked
rank), L*D7^+ in J_8 and M*D7^- in J_8; (4) degrees/parities.  Appends one JSON line per check to referee.jsonl."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import excl_lib as W
from hodge import core
from hodge.core import NV, pderiv, mul, add, sub, scale, reduce_mod, monomials, degree, parity_of
from fractions import Fraction
OUT = '/tmp/hodge/results/short-power-presentation/wave2-excl/referee/referee.jsonl'
def log(**r):
    r['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); r['lens'] = 'referee'
    open(OUT, 'a').write(json.dumps(r) + '\n'); print(json.dumps(r)[:400], flush=True)
F, A, q, rec = W.load()
ex = json.load(open('/tmp/hodge/results/short-power-presentation/wave2-excl/referee/excl_kernel_D7.json'))
ctx = W._mpoly_ctx()
Fi = [pderiv(F, i) for i in range(NV)]
Fm = [W._to_mpoly(f, ctx) for f in Fi]
Am = W._to_mpoly(A, ctx)
p = 1073741741
Fp = reduce_mod(F, p)
Lpoly = core.parse('X0+2*X1+3*X2'); Mpoly = core.parse('Y0+3*Y1+5*Y2')
for par in ['0', '1']:
    D = core.parse(ex[par]['D']); N = [core.parse(s) for s in ex[par]['N']]
    t0 = time.time()
    Dm = W._to_mpoly(D, ctx); Nm = [W._to_mpoly(n, ctx) for n in N]
    lhs = Am * Dm; rhs = sum((Nm[i] * Fm[i] for i in range(NV)), ctx.from_dict({}))
    ok_id = (lhs == rhs)
    degs = {'deg D': degree(D), 'par D': parity_of(D), 'deg N': [degree(n) for n in N if n], 'nterms D': len(D)}
    # cone property
    if par == '0':
        d1 = sub(pderiv(D, 4), scale(pderiv(D, 3), 3)); d2 = sub(pderiv(D, 5), scale(pderiv(D, 3), 5)); ring = 'Q[X0,X1,X2,M]'
    else:
        d1 = sub(pderiv(D, 1), scale(pderiv(D, 0), 2)); d2 = sub(pderiv(D, 2), scale(pderiv(D, 0), 3)); ring = 'Q[Y0,Y1,Y2,L]'
    cone = (d1 == {} and d2 == {})
    # mod p: D not in J_7 (with Hilbert check), L*D7+ in J_8 / M*D7- in J_8
    Dp = reduce_mod(D, p)
    g, info7 = W.membership_mod_p(Fp, Dp, p)
    gens7 = W.jac_gens(Fp, 7, int(par), p)
    mons7 = monomials(7, int(par)); idx7 = {m: i for i, m in enumerate(mons7)}
    rJ7 = W.cols_matrix([g_[2] for g_ in gens7], idx7, p).rank()
    hil7 = (rJ7 == W.dimJ(7, int(par)))
    mult = mul(Lpoly if par == '0' else Mpoly, D)
    multp = reduce_mod(mult, p)
    g8, info8 = W.membership_mod_p(Fp, multp, p)
    gens8 = W.jac_gens(Fp, 8, int(par), p)
    mons8 = monomials(8, int(par)); idx8 = {m: i for i, m in enumerate(mons8)}
    rJ8 = W.cols_matrix([g_[2] for g_ in gens8], idx8, p).rank()
    hil8 = (rJ8 == W.dimJ(8, int(par)))
    log(exp='referee-D7-exact', parD=int(par), identity_AD_eq_NF_over_Q=bool(ok_id), cone_ring=ring, cone=bool(cone),
        D_in_J7_modp=bool(info7['consistent']), rankJ7_modp=rJ7, rankJ7_hilbert=W.dimJ(7, int(par)), hilbert7_ok=bool(hil7),
        linear_times_D_in_J8_modp=bool(info8['consistent']), rankJ8_modp=rJ8, rankJ8_hilbert=W.dimJ(8, int(par)), hilbert8_ok=bool(hil8),
        prime=p, seconds=round(time.time() - t0, 1), **degs,
        note='identity exact => dim_Q K_7/J_7 >= 1 given D not in J_7 (Hilbert-checked, see referee README); cone via derivations')
