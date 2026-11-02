"""Exactly-true synthetic certificate identities (template shape) with controllable size.

Construction (q = 1 or 2): v = divergence-free field from antisymmetric M01=(X0 X1)^(3q-2), M34=(Y0 Y1)^(3q-2):
  v = (d1 M01, -d0 M01, 0, d4 M34, -d3 M34, 0),  div v = 0.
  D = random homogeneous of degree d with nD terms;  N_i = D * v_i;  A = q * sum_i dF_i v_i;  Q = 0.
Then  q * sum dF_i N_i = A D  and  D * div N - sum N_i dD_i = D^2 div v = 0 = F^q * Q  hold exactly.

Usage: python3 synth.py --F m3-coupled|dense-rational --d 4 --nD 40 [--q 1] [--seed 0] [--which e1,e2] [--yeven]
Writes src/Cap/S_<label>.lean and src/Cap/S_<label>.json (meta), prints the label.
"""
import sys, json, argparse, random, itertools, os
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV
from fractions import Fraction
import ident
CAP = '/tmp/hodge/lean/capacity'
FSHORT = {'m3-coupled': 'm3', 'dense-rational': 'dense', 'contact46006': 'c46'}

def monos(deg, yeven=False):
    out = []
    for e in itertools.product(range(deg + 1), repeat=NV):
        if sum(e) == deg and (not yeven or (e[3] + e[4] + e[5]) % 2 == 0):
            out.append(e)
    return out

def mono(e, c=1):
    return {tuple(e): Fraction(c)}

def build(Fid, d, nD, q=1, seed=0, yeven=False):
    F, A0, q0, rec = core.load_instance(Fid)
    rng = random.Random(seed)
    ms = monos(d, yeven)
    if nD >= len(ms):
        chosen = ms
    else:
        chosen = rng.sample(ms, nD)
    D = {}
    for e in chosen:
        c = rng.choice([c for c in range(-9, 10) if c])
        D[e] = Fraction(c)
    k = 3 * q - 2
    M01 = mono([k, k, 0, 0, 0, 0]); M34 = mono([0, 0, 0, k, k, 0])
    v = [core.pderiv(M01, 1), core.scale(core.pderiv(M01, 0), -1), {},
         core.pderiv(M34, 4), core.scale(core.pderiv(M34, 3), -1), {}]
    N = [core.mul(D, v[i]) for i in range(NV)]
    A = {}
    for i in range(NV):
        A = core.add(A, core.mul(core.pderiv(F, i), v[i]))
    A = core.scale(A, q)
    Q = {}
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    assert ok, info
    return F, A, q, D, N, Q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--F', default='m3-coupled'); ap.add_argument('--d', type=int, default=4)
    ap.add_argument('--nD', type=int, default=40); ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--which', default='e1,e2')
    ap.add_argument('--yeven', action='store_true'); ap.add_argument('--label', default=None)
    ap.add_argument('--hb', type=int, default=2000000); ap.add_argument('--rd', type=int, default=4000)
    a = ap.parse_args()
    F, A, q, D, N, Q = build(a.F, a.d, a.nD, a.q, a.seed, a.yeven)
    which = tuple(a.which.split(','))
    label = a.label or f'{FSHORT.get(a.F, a.F)}_d{a.d}_n{len(D)}' + (f'_q{q}' if q != 1 else '') + ('' if which == ('e1', 'e2') else '_' + ''.join(which))
    ns = f'HodgeComputer.Instances.S_{label}'
    st = ident.stats(F, A, q, D, N, Q)
    st.update({'label': label, 'Fid': a.F, 'which': list(which), 'seed': a.seed, 'yeven': a.yeven})
    txt = ident.module_text(ns, F, A, q, D, N, Q, which, a.hb, a.rd)
    path = f'{CAP}/src/Cap/S_{label}.lean'
    open(path, 'w').write(txt)
    st['bytes'] = len(txt.encode())
    json.dump(st, open(f'{CAP}/src/Cap/S_{label}.json', 'w'), indent=1)
    print(label, json.dumps(st))

if __name__ == '__main__':
    main()
