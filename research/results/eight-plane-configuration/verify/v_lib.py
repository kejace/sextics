"""Independent (engine-free) polynomial toolkit for refereeing the eight-plane-configuration lane.
Polynomials: dict {exp tuple (6 ints): Fraction or int}.  Vars X0,X1,X2,Y0,Y1,Y2 = indices 0..5.
Written from scratch; does not import hodge.*"""
import re, sys, json, time
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb

NV = 6
VARS = ['X0','X1','X2','Y0','Y1','Y2']
REC = '/home/kejace/devel/sextics/data/eight-plane-configuration/instance.rec'
LOG = '/tmp/hodge/results/eight-plane-configuration/wave2-verify.log.jsonl'

def log(entry):
    entry = dict(entry); entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); entry.setdefault('agent', 'verify')
    with open(LOG, 'a') as f: f.write(json.dumps(entry) + '\n')

_term_re = re.compile(r'\s*([+-])?\s*([^+-]+)')
def parse_poly(s):
    """Parse sums of terms  c*X0**a*Y1**b / den  (rec / witness syntax)."""
    s = s.replace(' ', '')
    # split on +/- at top level (no parentheses expected)
    terms = re.findall(r'[+-]?[^+-]+', s)
    P = {}
    for t in terms:
        sign = 1
        if t[0] == '-': sign = -1; t = t[1:]
        elif t[0] == '+': t = t[1:]
        num, _, den = t.partition('/')
        c = Fraction(1); e = [0]*NV
        for f in re.split(r'(?<!\*)\*(?!\*)', num):
            if f == '': continue
            if f[0] in 'XY':
                m = re.fullmatch(r'([XY][012])(?:\*\*(\d+))?', f)
                if not m: raise ValueError(f)
                e[VARS.index(m.group(1))] += int(m.group(2) or 1)
            else:
                c *= Fraction(int(f))
        if den: c /= Fraction(int(den))
        e = tuple(e)
        P[e] = P.get(e, 0) + sign * c
    return {e: c for e, c in P.items() if c != 0}

def load_instance():
    txt = open(REC).read()
    lines = txt.split('\n')
    fields = {}; cur = None
    for ln in lines:
        if ln.startswith('#'): continue
        m = re.match(r'^([A-Za-z_]+):\s*(.*)$', ln)
        if m and not ln.startswith(' '):
            cur = m.group(1); fields[cur] = m.group(2)
        elif cur and ln.strip():
            fields[cur] += ' ' + ln.strip()
    F = parse_poly(fields['F']); A = parse_poly(fields['A']); q = int(fields['q'])
    return F, A, q, fields

def deg(P): return sum(next(iter(P))) if P else -1
def ypar(e): return (e[3]+e[4]+e[5]) % 2
def parity_of(P):
    ps = {ypar(e) for e in P}
    return ps.pop() if len(ps) == 1 else None

def monomials(d, par=None):
    out = []
    def rec(i, rem, cur):
        if i == NV-1:
            e = tuple(cur + [rem])
            if par is None or ypar(e) == par: out.append(e)
            return
        for k in range(rem, -1, -1):
            rec(i+1, rem-k, cur+[k])
    rec(0, d, [])
    return out

def add(P, Q, p=None):
    R = dict(P)
    for e, c in Q.items():
        v = R.get(e, 0) + c
        if p: v %= p
        if v: R[e] = v
        elif e in R: del R[e]
    return R
def scale(P, c, p=None):
    R = {}
    for e, v in P.items():
        w = v*c
        if p: w %= p
        if w: R[e] = w
    return R
def mul(P, Q, p=None):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a+b for a, b in zip(e1, e2))
            v = R.get(e, 0) + c1*c2
            if p: v %= p
            R[e] = v
    return {e: v for e, v in R.items() if v}
def shift(P, m, c=1, p=None):
    R = {}
    for e, v in P.items():
        w = v*c
        if p: w %= p
        if w: R[tuple(a+b for a, b in zip(e, m))] = w
    return R
def pderiv(P, i, p=None):
    R = {}
    for e, v in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            w = v*e[i]
            if p: w %= p
            if w: R[tuple(e2)] = w
    return R
def reduce_mod(P, p):
    R = {}
    for e, v in P.items():
        if isinstance(v, Fraction):
            if v.denominator % p == 0: raise ZeroDivisionError('p divides a denominator')
            w = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
        else:
            w = v % p
        if w: R[e] = w
    return R
def to_str(P):
    def mono(e):
        return '*'.join(f'{VARS[i]}**{k}' if k > 1 else VARS[i] for i, k in enumerate(e) if k) or '1'
    return ' + '.join(f'{c}*{mono(e)}' for e, c in sorted(P.items(), reverse=True))

# ---- parity-refined dimensions for S and for the Jacobian ring of a smooth sextic in 6 vars
def dimS(d, par=None):
    if d < 0: return 0
    if par is None: return comb(d+5, 5)
    return sum(comb(k+2, 2)*comb(d-k+2, 2) for k in range(d+1) if k % 2 == par)

def dimJ_exact(k, par):
    """dim of parity-par part of J_k = (F_0..F_5), F_i of degree 5, F even in Y (F_X even, F_Y odd),
    via the Koszul complex (regular sequence): sum_j (-1)^j dim(Lambda^j-generators part)."""
    # generators: 3 even (F_X) + 3 odd (F_Y), degree 5 each.  Lambda^j of the 6 generators in degree 5j:
    # a wedge of a X-gens and b Y-gens (a+b=j) has parity b%2 and multiplicity C(3,a)C(3,b).
    tot = 0
    for j in range(1, 7):
        for a in range(0, j+1):
            b = j - a
            if a > 3 or b > 3: continue
            mult = comb(3, a)*comb(3, b)
            # need monomial m of degree k-5j with parity par - b
            pm = (par - b) % 2
            tot += (-1)**(j+1) * mult * dimS(k-5*j, pm)
    return tot

def dimR_exact(k, par): return dimS(k, par) - dimJ_exact(k, par)
