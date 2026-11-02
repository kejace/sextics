"""Wave-2 job: E1-admissible denominator spaces per degree (kernel of mult. by [A] on R_d, both parities).

usage: python3 w2_kernel.py --degrees 11,12 --pars 0,1 [--tag kernel]
Records for each (d, par): dim of the kernel mod J_d, dim J_d, dim (R/J)_d, and saves the basis (mod p)
as a pickle; also records the parity-split Hilbert function of R predicted by the Koszul complex and the
duality prediction  dim Ann(A)_d^{par} = dim R_d^{par} - rank(A) with rank(A) = min(dim R_d^par, dim R_{d+6}^{par'}).
"""
import sys, os, json, time, argparse, pickle
from math import comb
ap = argparse.ArgumentParser()
ap.add_argument('--degrees', default='11,12')
ap.add_argument('--pars', default='0,1')
ap.add_argument('--tag', default='kernel')
ap.add_argument('--no-basis', type=int, default=0)
args = ap.parse_args()
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.core import monomials
IID = 'eight-plane-configuration'
WD = f'/tmp/hodge/results/{IID}'
F, A, q, rec = core.load_instance(IID)
PRIME = core.PRIMES[1]
OUT = f'{WD}/{args.tag}.jsonl'


def jlog(obj):
    with open(OUT, 'a') as f:
        f.write(json.dumps(obj) + '\n')


def nS(d, par):
    return len(monomials(d, par)) if d >= 0 else 0


def koszul_dims(dmax=24):
    """Parity-split Hilbert function of R = S/J via the Koszul complex on (F_0,F_1,F_2 even; F_3,F_4,F_5 odd)."""
    out = {}
    for d in range(dmax + 1):
        for par in (0, 1):
            tot = 0
            for k in range(0, 7):
                for a in range(0, 4):
                    b = k - a
                    if b < 0 or b > 3:
                        continue
                    sign = (-1) ** k
                    # Lambda^k generator with a even and b odd partials: shift degree 5k, parity b
                    tot += sign * comb(3, a) * comb(3, b) * nS(d - 5 * k, (par - b) % 2)
            out[(d, par)] = tot
    return out


Rdims = koszul_dims()
pred = {}
for d in range(0, 19):
    for par in (0, 1):
        rd = Rdims[(d, par)]; rt = Rdims[(d + 6, 1 - par)]
        pred[f'{d},{par}'] = {'R_d': rd, 'R_d+6_target': rt, 'min_kernel': max(0, rd - rt)}
jlog({'event': 'hilbert', 'R_dims_by_parity': {f'{d},{p}': v for (d, p), v in Rdims.items()}, 'kernel_prediction': pred,
      'note': 'kernel prediction assumes maximal rank of multiplication by [A] (verified mod p for d<=10)'})
print(json.dumps(pred), flush=True)

for d in [int(x) for x in args.degrees.split(',')]:
    for par in [int(x) for x in args.pars.split(',')]:
        t0 = time.time()
        res = jacobian.mult_kernel(F, A, d, p=PRIME, parity=par, verbose=True)
        rec_ = {'event': 'mult_kernel', 'd': d, 'par': par, 'p': PRIME, 'kernel_mod_J': res['dim'],
                'kernel_dim_total': res['kernel_dim'], 'dim_J_d': res['jacobian_dim'], 'dim_R_d': res['ring_dim'],
                'rows': res['rows'], 'cols': res['cols'], 'predicted_min_kernel': pred[f'{d},{par}']['min_kernel'],
                'predicted_R_d': Rdims[(d, par)], 't_build': res['t_build'], 't_null': res['t_null'],
                'seconds': round(time.time() - t0, 1), 'engine': 2}
        jlog(rec_)
        print(json.dumps(rec_), flush=True)
        if res['basis'] and not args.no_basis:
            with open(f'{WD}/e1kernel_basis_d{d}_par{par}.pkl', 'wb') as f:
                pickle.dump(res['basis'], f)
