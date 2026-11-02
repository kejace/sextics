#!/usr/bin/env python3
"""referee_replay.py -- independent replay of the wave2-excl exact certificates (lane rank-one-biquartic).

Own polynomial arithmetic (dict exponent-tuple -> Fraction), own parser (sympy) for the instance record and the
K_7 generators, python-flint nmod_mat for modular ranks.  The ONLY things taken from the packet's w2lib / hodge.core
are (a) the row-index -> monomial map and (b) the polar-basis label order, both needed to read the packet's
certificate JSON files.  Everything that is *checked* (E2 columns, J generators, pairings) is rebuilt here.

Parts (all run by default; pass a subset of  d5 d6 e1 k7  to select):
  d5 : the 6 exact certificates lambda^T E = a_k^2 of fam_d5_p1000000007_all.certs.json against MY E2 columns
       (a_k a_l, a_k c, Q) -> proves a_k = 0 on the degree-5 solution variety, i.e. no degree-5 denominator.
  d6 : the 58 exact certificates of hpipe_d6.blk1.certs.json (25 coordinate squares + 33 quadratic forms):
       lambda kills every Y-odd generator m*F_i of J_12 and pairs with the TRUE E2 columns P_kl + P_lk
       (dense B_D(N), not the sparse Hessian representatives) to the claimed quadratic form.
  e1 : E1-kernel dimensions dim{D in S_d : A D in J_{d+6}} (all parities together) at the unused prime
       1000000009 for d = 4..7, plus the regular-sequence rank check of the J block.
  k7 : the 3 exact non-Jacobian degree-7 generators of nonj_d7.json: N(F) = A D exactly over Q, and the three
       classes are independent modulo J_7 at the new prime (non-membership mod p is rigorous over Q).
Receipt: referee_replay.json / referee_replay.out in the packet directory.
"""
import sys, json, time, itertools, re
from fractions import Fraction
from flint import nmod_mat
import sympy

PK = '/tmp/hodge/results/rank-one-biquartic/wave2-excl'
REC = '/home/kejace/devel/sextics/data/rank-one-biquartic/instance.rec'
NV = 6
VARS = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')
SYMS = sympy.symbols(VARS)
PNEW = 1000000009
T0 = time.time()
OUT = {}
LOG = open(f'{PK}/referee_replay.out', 'a')


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


# ----------------------------------------------------------------------------- polynomials
def parse(s):
    expr = sympy.sympify(s, locals=dict(zip(VARS, SYMS)))
    P = sympy.Poly(sympy.expand(expr), *SYMS)
    out = {}
    for e, c in P.terms():
        c = sympy.Rational(c)
        out[tuple(int(x) for x in e)] = Fraction(int(c.p), int(c.q))
    return out


def add(P, Q):
    R = dict(P)
    for e, v in Q.items():
        w = R.get(e, 0) + v
        if w: R[e] = w
        elif e in R: del R[e]
    return R


def sub(P, Q):
    return add(P, {e: -v for e, v in Q.items()})


def scale(P, c):
    c = Fraction(c)
    return {e: v * c for e, v in P.items()} if c else {}


def mul(P, Q):
    R = {}
    for e1, v1 in P.items():
        for e2, v2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            w = R.get(e, 0) + v1 * v2
            if w: R[e] = w
            elif e in R: del R[e]
    return R


def pderiv(P, i):
    R = {}
    for e, v in P.items():
        if e[i]:
            f = list(e); f[i] -= 1
            R[tuple(f)] = v * e[i]
    return R


def mono(e, c=1):
    return {tuple(e): Fraction(c)}


def degree(P):
    return max(sum(e) for e in P) if P else -1


def ypar(e):
    return (e[3] + e[4] + e[5]) % 2


def monomials(deg, par=None):
    out = []
    for e in itertools.product(range(deg + 1), repeat=NV):
        if sum(e) == deg and (par is None or ypar(e) == par):
            out.append(e)
    return out


def divergence(N):
    out = {}
    for i in range(NV):
        if N[i]: out = add(out, pderiv(N[i], i))
    return out


def apply_field(N, D):
    out = {}
    for i in range(NV):
        if N[i]: out = add(out, mul(N[i], pderiv(D, i)))
    return out


def B_op(D, N):
    """E2 operator  B_D(N) = D div N - N(D)."""
    return sub(mul(D, divergence(N)), apply_field(N, D))


def pair(lam, P):
    """<lambda, P> for lambda a dict monomial -> Fraction."""
    s = Fraction(0)
    for e, v in P.items():
        c = lam.get(e)
        if c is not None: s += c * v
    return s


def fermat_hf(k):
    """dim (S/(z_i^5))_k = regular-sequence Hilbert function of any smooth sextic's Jacobian ring."""
    return sum(1 for e in monomials(k) if max(e) <= 4)


def modp(v, p):
    v = Fraction(v)
    return (v.numerator * pow(v.denominator, -1, p)) % p


# ----------------------------------------------------------------------------- instance
def load_rec():
    txt = open(REC).read().splitlines()
    fields = {}; key = None
    for line in txt:
        if line.startswith('#'): continue
        m = re.match(r'^([A-Za-z]+):\s*(.*)$', line)
        if m and not line.startswith(' '):
            key = m.group(1); fields[key] = m.group(2)
        elif key and line.startswith(' '):
            fields[key] += ' ' + line.strip()
    return fields


REC_ = load_rec()
assert REC_['q'].strip() == '1'
F = parse(REC_['F']); A = parse(REC_['A'])
Fi = [pderiv(F, i) for i in range(NV)]
assert degree(F) == 6 and degree(A) == 6 and len(A) == 1 and A == mono((3, 0, 0, 3, 0, 0))
say(f'[setup] F parsed: {len(F)} terms, degree 6; A = X0^3 Y0^3; q = 1')

# packet helpers (row maps / label order only)
sys.path.insert(0, PK)
import w2lib as W  # noqa: E402  (loads hodge.core; used only for FS.rows / monomial order / labels)


def wpoly_to_mine(P):
    return {tuple(e): Fraction(v) for e, v in P.items()}


# ----------------------------------------------------------------------------- part d5
def part_d5():
    t0 = time.time()
    FS = W.FamilySystem(5, W.polar_basis(5), verbose=False)
    inv = {i: tuple(e) for e, i in FS.rows.items()}
    assert all(sum(e) == 10 for e in inv.values())
    certs = json.load(open(f'{PK}/fam_d5_p1000000007_all.certs.json'))
    # my E2 columns for the family D = sum a_k F_k, N = sum a_k A e_k + K(c), Q
    cols = []
    def P_kl(k, l):  # B_{D_k}(N_l) with D_k = F_k, N_l = A e_l
        return sub(mul(Fi[k], pderiv(A, l)), mul(A, pderiv(Fi[k], l)))
    for k in range(NV):
        for l in range(k, NV):
            cols.append((('aa', k, l), P_kl(k, k) if k == l else add(P_kl(k, l), P_kl(l, k))))
    kos = []
    for i in range(NV):
        for j in range(i + 1, NV):
            for m in monomials(1):
                K = [{} for _ in range(NV)]
                K[i] = mul(mono(m), Fi[j]); K[j] = scale(mul(mono(m), Fi[i]), -1)
                kos.append(((i, j, m), K))
    for k in range(NV):
        for lab, K in kos:
            cols.append((('ac', k, lab), B_op(Fi[k], K)))
    for m in monomials(4):
        cols.append((('Q', m), mul(mono(m), F)))
    say(f'[d5] my E2 system: {len(cols)} columns ({NV*(NV+1)//2} aa, {NV*len(kos)} ac, {len(monomials(4))} Q); '
        f'packet rows {FS.nrows} cols {FS.ncols}')
    assert len(cols) == FS.ncols
    res = {}
    for k in range(NV):
        name = f"a[('J', {k}, (0, 0, 0, 0, 0, 0))]"
        lam = {inv[int(r)]: Fraction(v) for r, v in certs[name].items()}
        bad = []
        for lab, P in cols:
            want = Fraction(1) if lab == ('aa', k, k) else Fraction(0)
            if pair(lam, P) != want: bad.append(lab)
        h = max(max(abs(x.numerator), x.denominator) for x in lam.values())
        res[f'a{k}'] = dict(verified=(not bad), support=len(lam), height_digits=len(str(h)), n_bad_columns=len(bad))
        say(f'[d5] certificate a_{k}^2: verified={not bad} support={len(lam)} height_digits={len(str(h))}'
            + (f' FIRST BAD {bad[:3]}' if bad else ''))
    OUT['d5'] = dict(result=('VERIFIED: a_k^2 in I for all k -> no degree-5 denominator (any parity)' if all(r['verified'] for r in res.values()) else 'FAILED'),
                     certs=res, seconds=round(time.time() - t0, 1))
    say(f'[d5] {OUT["d5"]["result"]} ({OUT["d5"]["seconds"]}s)')


# ----------------------------------------------------------------------------- part d6
def part_d6(limit=None):
    t0 = time.time()
    labels = [b[0] for b in W.polar_basis(6)]
    assert len(labels) == 36
    rowinv = {i: tuple(e) for i, e in enumerate(W.core.monomials(12, 1))}
    certs = json.load(open(f'{PK}/hpipe_d6.blk1.certs.json'))
    Cj = json.load(open(f'{PK}/hpipe_d6.C.json'))
    forms = {n: {tuple(int(x) for x in re.findall(r'\d+', kk)): Fraction(v) for kk, v in fm.items()} for n, fm in Cj['exact_forms']['1'].items()}
    # my basis (D_k, N_k) from the labels ('J', i, m): D = m F_i, N = A m e_i
    basis = []
    for lab in labels:
        _, i, m = lab
        D = mul(mono(m), Fi[i]); N = [{} for _ in range(NV)]; N[i] = mul(A, mono(m))
        assert sub(apply_field(N, F), mul(A, D)) == {}, 'E1 lift check failed'
        basis.append((D, N))
    # true E2 quadratic columns (dense), all 666 pairs
    cols = {}
    for k in range(36):
        for l in range(k, 36):
            Dk, Nk = basis[k]; Dl, Nl = basis[l]
            cols[(k, l)] = B_op(Dk, Nk) if k == l else add(B_op(Dk, Nl), B_op(Dl, Nk))
    # Y-odd generators of J_12
    Jgens = [mul(mono(m), Fi[i]) for i in range(NV) for m in monomials(7) if ypar(tuple(m[j] + (1 if j == i else 0) for j in range(NV))) == 0]
    # (m F_i has the parity of m z_i; odd rows <-> parity 1 of m*F_i, i.e. parity of monomial m*z_i == 1)
    Jgens = [mul(mono(m), Fi[i]) for i in range(NV) for m in monomials(7)
             if ypar(tuple(m[j] + (1 if j == i else 0) for j in range(NV))) == 1]
    say(f'[d6] my basis 36, true quadratic columns {len(cols)}, Y-odd J_12 generators {len(Jgens)} (packet: 2376)')
    assert len(Jgens) == 2376
    names = [n for n, v in certs.items() if v]
    if limit: names = names[:limit]
    res = {}; nok = 0
    for name in names:
        lam = {rowinv[int(r)]: Fraction(v) for r, v in certs[name].items()}
        assert all(sum(e) == 12 and ypar(e) == 1 for e in lam)
        if name.startswith('sq:'):
            k = int(name[4:]); target = {(k, k): Fraction(1)}
        else:
            target = forms[name]
        badJ = sum(1 for g in Jgens if pair(lam, g) != 0)
        badC = [kl for kl, P in cols.items() if pair(lam, P) != target.get(kl, Fraction(0))]
        ok = (badJ == 0 and not badC); nok += ok
        h = max(max(abs(x.numerator), x.denominator) for x in lam.values())
        res[name] = dict(verified=ok, support=len(lam), height_digits=len(str(h)), bad_J=badJ, bad_cols=len(badC), form_terms=len(target))
        say(f'[d6] {name}: verified={ok} support={len(lam)} height_digits={len(str(h))} form_terms={len(target)}'
            + ('' if ok else f' badJ={badJ} badcols={badC[:3]}'))
    OUT['d6'] = dict(result=(f'VERIFIED: {nok}/{len(names)} exact certificates (lambda kills J_12^- and reproduces the quadratic form on the TRUE E2 columns)' if nok == len(names) else f'FAILED {len(names)-nok}'),
                     n_verified=nok, n_total=len(names), certs=res, seconds=round(time.time() - t0, 1))
    say(f'[d6] {OUT["d6"]["result"]} ({OUT["d6"]["seconds"]}s)')


# ----------------------------------------------------------------------------- part e1
def part_e1(p=PNEW, ds=(4, 5, 6, 7)):
    t0 = time.time()
    res = {}
    for d in ds:
        t1 = time.time()
        rows = monomials(d + 6); ridx = {e: i for i, e in enumerate(rows)}
        Sd = monomials(d); Sd1 = monomials(d + 1)
        colsA = [mul(A, mono(m)) for m in Sd]
        colsJ = [mul(mono(m), Fi[i]) for i in range(NV) for m in Sd1]
        def build(colist):
            R, C = len(rows), len(colist); flat = [0] * (R * C)
            for j, P in enumerate(colist):
                for e, v in P.items(): flat[ridx[e] * C + j] = modp(v, p)
            return nmod_mat(R, C, flat, p)
        rJ = build(colsJ).rank()
        rAll = build(colsA + colsJ).rank()
        dimV = len(Sd) - (rAll - rJ)
        expJ = len(rows) - fermat_hf(d + 6)
        dimJd = len(Sd) - fermat_hf(d)
        res[d] = dict(p=p, dim_S_d=len(Sd), rank_J=rJ, rank_J_expected=expJ, rankJ_ok=(rJ == expJ), dim_V_d_modp=dimV, dim_J_d=dimJd,
                      kernel_mod_J=dimV - dimJd, seconds=round(time.time() - t1, 1))
        say(f'[e1] d={d} p={p}: rank J_{d+6} = {rJ} (regular-sequence value {expJ}, ok={rJ==expJ}); dim V_d mod p = {dimV}, dim J_d = {dimJd}, '
            f'kernel mod J = {dimV-dimJd} ({res[d]["seconds"]}s)')
    OUT['e1'] = dict(result='dim_p V_d >= dim_Q V_d whenever rankJ_ok; kernel_mod_J = 0 => V_d = J_d over Q', table=res, seconds=round(time.time() - t0, 1))


# ----------------------------------------------------------------------------- part k7
def part_k7(p=PNEW):
    t0 = time.time()
    gens = json.load(open(f'{PK}/nonj_d7.json'))
    rows = monomials(7); ridx = {e: i for i, e in enumerate(rows)}
    J7 = [mul(mono(m), Fi[i]) for i in range(NV) for m in monomials(2)]
    def build(colist):
        R, C = len(rows), len(colist); flat = [0] * (R * C)
        for j, P in enumerate(colist):
            for e, v in P.items(): flat[ridx[e] * C + j] = modp(v, p)
        return nmod_mat(R, C, flat, p)
    rJ = build(J7).rank()
    Ds = []; res = {}
    for g in gens:
        D = parse(g['D']); N = [parse(s) for s in g['N']]
        e1 = (sub(apply_field(N, F), mul(A, D)) == {})
        pars = {ypar(e) for e in D}
        basetype = all(e[3] + e[4] + e[5] in (0, 2, 4, 6) for e in D)
        res[g['label']] = dict(E1_exact=e1, deg_D=degree(D), terms_D=len(D), deg_N=[degree(x) for x in N], y_parity=sorted(pars))
        say(f"[k7] {g['label']}: N(F) = A D exactly over Q: {e1}; deg D {degree(D)}, {len(D)} terms, Y-parity {sorted(pars)}")
        Ds.append(D)
    rAug = build(J7 + Ds).rank()
    say(f'[k7] p={p}: rank J_7 = {rJ} (expected {len(rows)-fermat_hf(7)}), rank [J_7 | D1 D2 D3] = {rAug} -> independent mod J_7: {rAug - rJ} of 3 (rigorous over Q)')
    OUT['k7'] = dict(gens=res, rank_J7=rJ, rank_J7_expected=len(rows) - fermat_hf(7), rank_aug=rAug, independent_mod_J=rAug - rJ,
                     result=('VERIFIED: 3 exact non-Jacobian E1 directions in degree 7' if rAug - rJ == 3 and all(r['E1_exact'] for r in res.values()) else 'FAILED'),
                     seconds=round(time.time() - t0, 1))


if __name__ == '__main__':
    parts = sys.argv[1:] or ['d5', 'd6', 'e1', 'k7']
    say(f'=== referee_replay.py start {time.strftime("%Y-%m-%d %H:%M:%S")} parts={parts} ===')
    for part in parts:
        {'d5': part_d5, 'd6': part_d6, 'e1': part_e1, 'k7': part_k7}[part]()
        json.dump({part: OUT[part]}, open(f'{PK}/referee_replay_{part}.json', 'w'), indent=1, default=str)
    say(f'=== done {round(time.time()-T0)}s ===')
