"""Monomial-count statistics of the certificate identities for accepted lanes + dense-rational F."""
import sys, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV
ROOT = '/home/kejace/devel/sextics'
def stats(F, A, q, D, N, Q):
    dF = [core.pderiv(F, i) for i in range(NV)]
    dD = [core.pderiv(D, i) for i in range(NV)]
    # e1: q * sum dF_i * N_i = A * D
    e1_pairs = sum(len(dF[i]) * len(N[i]) for i in range(NV)) + len(A) * len(D)
    lhs1 = {}
    for i in range(NV):
        lhs1 = core.add(lhs1, core.mul(dF[i], N[i]))
    rhs1 = core.mul(A, D)
    lhs1 = core.scale(lhs1, q)
    # e2: D * sum d_i N_i - sum dD_i N_i = F^q Q
    divN = {}
    for i in range(NV):
        divN = core.add(divN, core.pderiv(N[i], i))
    e2_pairs = len(D) * len(divN) + sum(len(dD[i]) * len(N[i]) for i in range(NV))
    Fq = core.power(F, q)
    e2_pairs += len(Fq) * len(Q)
    lhs2 = core.sub(core.mul(D, divN), _sumprod(dD, N))
    rhs2 = core.mul(Fq, Q)
    return {
        'F': len(F), 'A': len(A), 'q': q, 'D': len(D), 'N': [len(n) for n in N], 'Q': len(Q), 'd': core.degree(D),
        'e1_pairs': e1_pairs, 'e1_monos': len(lhs1), 'e2_pairs': e2_pairs, 'e2_monos': len(lhs2),
        'e1_ok': core.sub(lhs1, rhs1) == {}, 'e2_ok': core.sub(lhs2, rhs2) == {},
        'pairs_total': e1_pairs + e2_pairs, 'monos_total': len(lhs1) + len(lhs2),
    }
def _sumprod(dD, N):
    s = {}
    for i in range(NV):
        s = core.add(s, core.mul(dD[i], N[i]))
    return s
out = {}
for iid in sorted(os.listdir(f'{ROOT}/data')):
    if not os.path.exists(f'{ROOT}/data/{iid}/witness.json'):
        continue
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    st = stats(F, A, q, D, N, Q)
    st['witness_lean_bytes'] = os.path.getsize(f'{ROOT}/data/{iid}/Witness.lean') if os.path.exists(f'{ROOT}/data/{iid}/Witness.lean') else None
    out[iid] = st
    print(iid, json.dumps(st))
F, A, q, rec = core.load_instance('dense-rational')
out['dense-rational(F,A only)'] = {'F': len(F), 'A': len(A), 'q': q}
print('dense-rational', len(F), len(A), q)
json.dump(out, open('/tmp/hodge/lean/capacity/results/sizes.json', 'w'), indent=1)
