"""Wave-2 EXCLUSION lens helpers (independent of wave-1 scripts except hodge.core / hodge.jacobian).
Lane short-power-presentation: q=1, A=X0^3 Y0^3.
All ranks are modulo a prime p (default core.PRIMES[1] = 2147483629, different from wave 1's PRIMES[0]).
Rigor convention: rank_p(M) <= rank_Q(M); for the Jacobian ideal of a smooth tau-invariant sextic the
parity-graded Hilbert function of R = S/J is forced, so rank_p J_k^par == expected  <=>  rank_p = rank_Q."""
import sys, os, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from flint import nmod_mat
from fractions import Fraction

LANE = 'short-power-presentation'
RES = '/tmp/hodge/results/short-power-presentation'
PK = os.path.join(RES, 'wave2-excl')
LOGF = os.path.join(RES, 'wave2-excl.log.jsonl')
P = core.PRIMES[1]

def log(rec):
    rec = dict(rec); rec.setdefault('lens', 'excl'); rec.setdefault('engine', 2); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOGF, 'a') as f: f.write(json.dumps(rec) + '\n')
    return rec

def load():
    return core.load_instance(LANE)  # F, A, q, rec

# ---- expected parity-graded Hilbert function of R = S/J(F), F smooth sextic, tau-invariant ----
def _polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out
def _polypow(a, n):
    out = [1]
    for _ in range(n): out = _polymul(out, a)
    return out
_H = _polypow([1, 1, 1, 1, 1], 6)                       # (1+t+..+t^4)^6
_Hchi = _polymul(_polypow([1, 1, 1, 1, 1], 3), _polypow([1, -1, 1, -1, 1], 3))  # equivariant trace of tau
def hilb_R(k, par=None):
    """dim R_k (par None) or dim R_k^{par} (par 0 even / 1 odd under Y -> -Y)."""
    h = _H[k] if k < len(_H) else 0
    c = _Hchi[k] if k < len(_Hchi) else 0
    if par is None: return h
    return (h + c) // 2 if par == 0 else (h - c) // 2
def dim_S(k, par=None):
    return len(core.monomials(k, par))
def expected_rank_J(k, par=None):
    return dim_S(k, par) - hilb_R(k, par)

def jac_rank(F, k, par, p=P):
    """rank mod p of J(F)_k^{par} via the engine's jacobian_matrix (dense nmod_mat)."""
    M, mons = jacobian.jacobian_matrix(F, k, p, parity=par)
    return M.rank(), len(mons)
