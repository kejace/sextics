"""Degree-d (default 10) E1-kernel analysis: structured-subspace intersections + certified fixed-D tests.

usage: python3 wave2_excl_k10.py <jobname> <workers> <parity> <primeindex> [d] [nrandom]
Needs wave2-excl-K_d<d>_par<par>_p<p>.json from wave2_excl_kernel.py.
Rigor: intersections computed mod p are upper bounds for the rational intersection dimensions when the generator
ranks and dim K are certified (dim K is certified when the parity dims sum to the Gorenstein bound).  Candidates that
rationally reconstruct (small heights) are exact polynomials and get the certified fixed-D test; mod-p-only
candidates (random kernel elements) are screening only.
"""
import sys, os, json, time, traceback, random, itertools
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, q, Qs, monomials, mul, add, sub, scale, shift, pderiv, power, logline, rank_of, var
from hodge import jacobian, fixed_d
from flint import nmod_mat
from fractions import Fraction
import multiprocessing as mp

JOB = sys.argv[1]; WORKERS = int(sys.argv[2]); PAR = int(sys.argv[3]); PI = int(sys.argv[4])
DD = int(sys.argv[5]) if len(sys.argv) > 5 else 10
NRAND = int(sys.argv[6]) if len(sys.argv) > 6 else 6
P = core.PRIMES[PI]
LOG = os.path.join(w2x.LANE, f'wave2-excl.{JOB}.jsonl')
KFILE = os.path.join(w2x.LANE, f'wave2-excl-K_d{DD}_par{PAR}_p{P}.json')


def key2mon(s):
    return tuple(int(x) for x in s.split())


def load_K():
    K = json.load(open(KFILE))
    Kb = [{key2mon(k): v for k, v in vec.items()} for vec in K['K_basis']]
    KmJ = [{key2mon(k): v for k, v in vec.items()} for vec in K['K_mod_J_basis']]
    return K, Kb, KmJ


def modp(v, p):
    if isinstance(v, Fraction):
        return (v.numerator * pow(v.denominator, -1, p)) % p
    return int(v) % p


def mat_of(polys, mons, p):
    col = {e: j for j, e in enumerate(mons)}
    flat = [0] * (len(polys) * len(mons))
    for r, Pl in enumerate(polys):
        for e, v in Pl.items():
            flat[r * len(mons) + col[e]] = modp(v, p)
    return nmod_mat(len(polys), len(mons), flat, p)


def intersection(Wgens, Kb, mons, p):
    """Returns (dim W, dim K, dim W cap K, list of a-vectors (coefficients in W-generators) spanning W cap K)."""
    m, r = len(Wgens), len(Kb)
    Mw = mat_of(Wgens, mons, p); rW = Mw.rank()
    Mk = mat_of(Kb, mons, p); rK = Mk.rank()
    Ms = mat_of(Wgens + Kb, mons, p); rS = Ms.rank()
    dimI = rW + rK - rS
    avecs = []
    if dimI > 0:
        X, nul = Ms.transpose().nullspace()
        rows = []
        for kk in range(nul):
            a = [int(X[i, kk]) for i in range(m)]
            if any(a):
                rows.append(a)
        if rows:
            Ar = nmod_mat(len(rows), m, [x for a in rows for x in a], p)
            Arr, rk = Ar.rref()
            avecs = [[int(Arr[i, j]) for j in range(m)] for i in range(rk)]
    return rW, rK, dimI, avecs


def ratlift(vec, p, maxden=10**6):
    out = []
    for v in vec:
        r = core.ratrecon(v % p, p)
        if r is None or r.denominator > maxden or abs(r.numerator) > maxden:
            return None
        out.append(r)
    return out


def structured_subspaces(par):
    """name -> (generator polys (exact), description)."""
    W = {}
    X0, X1, X2, Y0, Y1, Y2 = [var(i) for i in range(6)]
    delta = w2x.delta_poly()
    if par == 0:
        # quintics in Q
        gens = []
        for c in itertools.combinations_with_replacement(range(6), DD // 2):
            Pl = {(0,) * 6: Fraction(1)}
            for i in c:
                Pl = mul(Pl, Qs[i])
            gens.append(Pl)
        W['CQ_half'] = (gens, f'degree-{DD//2} forms in Q_0..Q_5 (pullbacks of hypersurface sections of the cubic)')
        W['delta_S_odd'] = ([shift(delta, m, 1) for m in monomials(DD - 6, 1)], 'delta times Y-odd forms of degree d-6')
        for S in [(0, 1, 2), (3, 4, 5), (0, 1, 3), (0, 3, 4)]:
            gens = []
            for c in itertools.combinations_with_replacement(S, DD // 2):
                Pl = {(0,) * 6: Fraction(1)}
                for i in c:
                    Pl = mul(Pl, Qs[i])
                gens.append(Pl)
            W[f'CQ_{"".join(map(str, S))}'] = (gens, f'forms of degree {DD//2} in Q_S, S={S}')
    else:
        W['delta_S_even'] = ([shift(delta, m, 1) for m in monomials(DD - 6, 0)], 'delta times Y-even forms of degree d-6')
        gens = []
        for c in itertools.combinations_with_replacement(range(6), (DD - 6) // 2):
            Pl = dict(delta)
            for i in c:
                Pl = mul(Pl, Qs[i])
            gens.append(Pl)
        W['delta_CQ'] = (gens, 'delta times forms in Q')
    return W


def plane_windows(par, mons):
    """Monomial subspaces: X-degree >= a and Y-degree >= b (vanishing orders along the coordinate planes)."""
    out = {}
    for a, b in [(2, 2), (3, 3), (4, 4), (5, 5), (2, 4), (4, 2), (3, 5), (5, 3), (2, 6), (6, 2), (4, 6), (6, 4), (0, 4), (4, 0), (0, 6), (6, 0), (0, 8), (8, 0)]:
        sel = [e for e in mons if (e[0] + e[1] + e[2]) >= a and (e[3] + e[4] + e[5]) >= b]
        if sel:
            out[f'planes_X{a}_Y{b}'] = ([{e: 1} for e in sel], f'monomials with X-degree>={a} and Y-degree>={b} (orders along X=0 / Y=0 planes)')
    return out


def k3_windows(par, p):
    """(Q_S)^2 cap (Q_S')^2 in degree DD, for the join split S=(0,1,2), S'=(3,4,5), and (Q_S)^2 alone."""
    out = {}
    def ideal_power_gens(S, k):
        gens = []
        for c in itertools.combinations_with_replacement(S, k):
            Pl = {(0,) * 6: Fraction(1)}
            for i in c:
                Pl = mul(Pl, Qs[i])
            for m in monomials(DD - 2 * k, par):
                gens.append(shift(Pl, m, 1))
        return gens
    out['QS012_sq'] = (ideal_power_gens((0, 1, 2), 2), '(Q0,Q1,Q2)^2 in degree d (double along the K3 curve C_012)')
    out['QS345_sq'] = (ideal_power_gens((3, 4, 5), 2), '(Q3,Q4,Q5)^2 in degree d')
    return out


def k3_surface_windows(par, p):
    """Denominators of multiplicity >= 2 along the K3 surfaces S = {Q_a + z Q_b = Q_c + z' Q_d = Q_e + z'' Q_f = 0}
    (pullbacks of the 405 planes of the Fermat cubic; z^3 = 1).  Rational for z = z' = z'' = 1 (15 pairings).
    For the complex ones we use a primitive cube root of unity in F_p (p = 1 mod 3): mod-p screening only."""
    out = {}
    pairings = []
    idx = list(range(6))
    def rec(rem, cur):
        if not rem:
            pairings.append(list(cur)); return
        a = rem[0]
        for b in rem[1:]:
            rec([x for x in rem if x not in (a, b)], cur + [(a, b)])
    rec(idx, [])
    # cube root of unity mod p
    w3 = None
    if p % 3 == 1:
        for g in range(2, 200):
            c = pow(g, (p - 1) // 3, p)
            if c != 1:
                w3 = c; break
    import itertools as it
    for pr in pairings:
        tag = ''.join(f'{a}{b}' for a, b in pr)
        for zs in it.product(range(3), repeat=3):
            exact = (zs == (0, 0, 0))
            if not exact and w3 is None:
                continue
            Ls = []
            for (a, b), zk in zip(pr, zs):
                if exact:
                    Ls.append(add(Qs[a], Qs[b]))
                else:
                    z = pow(w3, zk, p)
                    Ls.append(add({e: modp(v, p) for e, v in Qs[a].items()}, scale({e: modp(v, p) for e, v in Qs[b].items()}, z, p), p))
            gens = []
            for i in range(3):
                for j in range(i, 3):
                    LL = mul(Ls[i], Ls[j]) if exact else mul(Ls[i], Ls[j], p)
                    for m in monomials(DD - 4, par):
                        gens.append(shift(LL, m, 1) if exact else shift(LL, m, 1, p))
            name = f'K3double_{tag}_z{zs[0]}{zs[1]}{zs[2]}'
            out[name] = (gens, f'I(S)^2 in degree d, S = K3 surface over the plane pairing {pr} with zeta exponents {zs}' + ('' if exact else ' (mod-p, omega in F_p)'))
    return out


def qS_cube_windows(par, p):
    out = {}
    for S in [(0, 1, 2), (3, 4, 5)]:
        gens = []
        for c in itertools.combinations_with_replacement(S, 3):
            Pl = {(0,) * 6: Fraction(1)}
            for i in c:
                Pl = mul(Pl, Qs[i])
            for m in monomials(DD - 6, par):
                gens.append(shift(Pl, m, 1))
        out[f'QS{"".join(map(str,S))}_cube'] = (gens, f'(Q_S)^3 in degree d, S={S} (triple along the K3 curve)')
    return out


def span_intersection_gens(G1, G2, mons, p):
    """Generators (mod p) of span(G1) cap span(G2)."""
    rW, rK, dimI, avecs = intersection(G1, G2, mons, p)
    gens = []
    for a in avecs:
        Pl = {}
        for coef, g in zip(a, G1):
            if coef:
                Pl = add(Pl, scale({e: modp(v, p) for e, v in g.items()}, coef, p), p)
        gens.append(Pl)
    return gens, dimI


def test_candidate(args):
    name, fam, D, note, exact = args
    try:
        rec = w2x.full_test(name, fam, D, P, LOG, note=note + (' [EXACT candidate]' if exact else ' [mod-p candidate: screening only]'))
        rec['exact_candidate'] = exact
    except Exception:
        rec = logline(LOG, dict(exp='fixedD', name=name, family=fam, prime=P, result='ERROR', error=traceback.format_exc()[-800:], job=JOB))
    print(json.dumps({k: v for k, v in rec.items() if k not in ('D_head',)}), flush=True)
    return rec


def main():
    t0 = time.time()
    K, Kb, KmJ = load_K()
    mons = monomials(DD, PAR)
    rec0 = dict(exp='k10-load', d=DD, parity=PAR, prime=P, dimK=K['dimK'], dimK_mod_J=K['dimK_mod_J'], dimJd_modp=K['dimJd_modp'], job=JOB)
    logline(LOG, rec0); print(json.dumps(rec0), flush=True)
    Jgens = [g[2] for g in jacobian.jacobian_generators(F, DD, P, PAR)]
    cands = []
    allW = {}
    allW.update(structured_subspaces(PAR))
    allW.update(plane_windows(PAR, mons))
    allW.update(k3_windows(PAR, P))
    allW.update(qS_cube_windows(PAR, P))
    allW.update(k3_surface_windows(PAR, P))
    # join window: (Q_S)^2 cap (Q_S')^2
    g1, _ = allW['QS012_sq']; g2, _ = allW['QS345_sq']
    gI, dI = span_intersection_gens(g1, g2, mons, P)
    allW['QS012sq_cap_QS345sq'] = (gI, f'(Q0,Q1,Q2)^2 cap (Q3,Q4,Q5)^2 in degree d (dim {dI}); deformed-join-type double curves')
    for name, (gens, desc) in allW.items():
        if not gens:
            continue
        t1 = time.time()
        gens_p = [{e: modp(v, P) for e, v in g.items()} for g in gens]
        rW, rK, dimI, avecs = intersection(gens_p, Kb, mons, P)
        rWJ, rJ, dimIJ, _ = intersection(gens_p, Jgens, mons, P)
        rec = dict(exp='k10-intersection', d=DD, parity=PAR, prime=P, subspace=name, desc=desc, n_gens=len(gens), dimW=rW, dimK=rK,
                   dim_W_cap_K=dimI, dim_W_cap_J=dimIJ, new_mod_J=dimI - dimIJ, seconds=round(time.time() - t1, 1), job=JOB)
        logline(LOG, rec); print(json.dumps(rec), flush=True)
        if dimI - dimIJ > 0 and dimI <= 12:
            # candidates: elements of W cap K not in J.  try exact lift of the a-vectors
            for ii, a in enumerate(avecs):
                Dp = {}
                for coef, g in zip(a, gens_p):
                    if coef:
                        Dp = add(Dp, scale(g, coef, P), P)
                inJ, _ = jacobian.in_jacobian(F, Dp, P)
                if inJ:
                    continue
                lift = ratlift(a, P) if all(isinstance(v, Fraction) for g in gens for v in g.values()) else None
                if lift is not None:
                    Dx = {}
                    for coef, g in zip(lift, gens):
                        if coef:
                            Dx = add(Dx, scale(g, coef))
                    # verify E1 at another prime
                    ok2, _ = fixed_d.e1_solvable(F, A, q, Dx, core.PRIMES[PI + 1])
                    cands.append((f'{name}_cap_K_{ii}', f'k10-{name}', Dx, f'exact lift; E1 at second prime: {ok2}', True))
                else:
                    cands.append((f'{name}_cap_K_{ii}_modp', f'k10-{name}', Dp, 'no small-height lift', False))
        elif dimI - dimIJ > 0:
            # large intersection: test two random (mod p) elements not in J
            random.seed(11)
            for ii in range(2):
                Dp = {}
                for a in avecs:
                    c = random.randrange(1, P)
                    for coef, g in zip(a, gens_p):
                        if coef:
                            Dp = add(Dp, scale(g, coef * c % P, P), P)
                cands.append((f'{name}_cap_K_rand{ii}_modp', f'k10-{name}', Dp, f'random element of a {dimI}-dim intersection (mod p)', False))
    # random elements of K mod J (mod p)
    random.seed(7)
    for ii in range(NRAND):
        Dp = {}
        for v in KmJ:
            c = random.randrange(1, P)
            Dp = add(Dp, scale(v, c, P), P)
        cands.append((f'K{DD}_random{ii}_modp', f'k10-random', Dp, 'random element of K_d/J_d (mod p)', False))
    rec = dict(exp='k10-candidates', d=DD, parity=PAR, n_candidates=len(cands), names=[c[0] for c in cands], job=JOB, seconds=round(time.time() - t0, 1))
    logline(LOG, rec); print(json.dumps(rec), flush=True)
    with mp.get_context('fork').Pool(WORKERS) as pool:
        for _ in pool.imap_unordered(test_candidate, cands):
            pass
    print('ALL DONE', flush=True)


if __name__ == '__main__':
    main()
