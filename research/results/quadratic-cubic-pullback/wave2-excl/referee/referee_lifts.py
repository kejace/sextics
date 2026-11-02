#!/usr/bin/env python3
"""Exact check of the lens's lifted K_H bases (degree 6): each lifted vector v (entries in {0,+-1,+-2}, symmetric lift of
the mod-268435399 rref basis) must satisfy  sum_{a<=b} v_ab (2-delta_ab) m_a m_b F_{i_a i_b}  in J_6  (both parity parts),
which (K_6 = J_6 => [A] injective on R_6) is equivalent to the packet's A(...) in J_12.  Own exact J_6 rref (fmpq)."""
import sys, os, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from hodge import core
import referee_replay as R
from flint import fmpq_mat, fmpq

red = {}
for tp in (0, 1):
    cols = R.monos(6, tp)
    gens = [core.shift(R.Fi[i], m, 1) for i in range(6) for m in R.monos(1, (tp - R.PARZ[i]) % 2)]
    M = R.qmat(gens, cols); Rr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while Rr[i, j] == 0: j += 1
        piv.append(j)
    rows = [{cols[k]: Fraction(int(Rr[i, k].p), int(Rr[i, k].q)) for k in range(len(cols)) if Rr[i, k] != 0} for i in range(rank)]
    red[tp] = (cols, piv, rows)

def in_J6(T):
    for tp in (0, 1):
        part = {e: v for e, v in T.items() if (e[3] + e[4] + e[5]) % 2 == tp}
        cols, piv, rows = red[tp]
        for pj, row in zip(piv, rows):
            c = part.get(cols[pj], 0)
            if c: part = core.sub(part, core.scale(row, c))
        if any(v != 0 for v in part.values()): return False
    return True

base = '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/receipts/'
for fn in ['wave2-excl-hess-hessL6-d6-par1-p268435399-KH.json', 'wave2-excl-hess-hessL6-d6-par0-p268435399-KH.json',
           'wave2-excl-hess-hessL6-d6-parall-p268435399-KH.json']:
    t0 = time.time()
    K = json.load(open(base + fn))
    wb = [(i, tuple(m)) for i, m in K['wb']]; Wkeys = [tuple(x) for x in K['Wkeys']]
    basis = [{int(a): v for a, v in row.items()} for row in K['basis']]
    bad = 0; indep_rows = []
    for row in basis:
        T = {}
        for a, v in row.items():
            (ia, ma), (ib, mb) = wb[Wkeys[a][0]], wb[Wkeys[a][1]]
            t = core.shift(R.Fij[ia][ib], tuple(x + y for x, y in zip(ma, mb)), 1)
            if Wkeys[a][0] != Wkeys[a][1]: t = core.scale(t, 2)
            T = core.add(T, core.scale(t, Fraction(v)))
        if not in_J6(T): bad += 1
    # independence of the lifted vectors over Q
    nW = len(Wkeys)
    Bm = fmpq_mat(len(basis), nW, [fmpq(basis[i].get(a, 0)) for i in range(len(basis)) for a in range(nW)])
    R.log(dict(check='lens_KH_lift_exact_membership_J6', file=fn, dim=len(basis), nW=nW, failures=bad, all_in_J6=(bad == 0),
               lifted_rank_exact=Bm.rank(), max_abs=max(abs(v) for row in basis for v in row.values()), seconds=round(time.time() - t0, 1),
               note='exact over Q via own J_6 rref; equivalent to A(w^T Hess F w) in J_12 by injectivity of [A] on R_6 (K_6 = J_6)'))
