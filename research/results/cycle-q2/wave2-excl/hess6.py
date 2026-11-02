"""Hessian relaxation at deg D = 6 (lens excl, lane cycle-q2), stage 1: build the quadrics.

K_6 = {D in S_6 : A D in J_18}  (E1-admissible denominators), computed sector-wise mod p (G-characters).
For a basis D_a of K_6 choose N_a with N_a(F2) = A D_a (F2 = 2F, i.e. 2 N_a(F) = A D_a).
E2 modulo J:  H(D) = D div N - N.grad D  must lie in J_18 (Koszul/Euler gauge terms are in J or vanish).
Bilinear form  H(D_a,D_b) = (D_a div N_b - N_b.grad D_a + D_b div N_a - N_a.grad D_b)/2,
H(sum c_a D_a) = sum_a c_a^2 H(D_a,D_a) + sum_{a<b} 2 c_a c_b H(D_a,D_b).
Normal forms modulo J_18 (sector-wise) give 426 quadratic forms in the dim K_6 (expected 62) coordinates.
Output: msolve / Singular input files; if the ideal is m-primary over F_p-bar then no deg-6 certificate exists
over C (reduction argument: E1 and H in J transfer to any good prime because coker of the Jacobian matrices is
torsion-free at primes of smoothness).
usage: python3 wave2-excl/hess6.py [--prime P] [--deg 6]
"""
import sys, json, time, argparse, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
from hodge import core
from hodge.core import mul, add, scale, pderiv, shift, monomials, reduce_mod
from flint import nmod_mat

def key_of(m, gens):
    return tuple(sum(a*b for a, b in zip(m, g)) % 6 for g in gens)

def poly_key(P, gens):
    ks = {key_of(m, gens) for m in P}
    assert len(ks) == 1, ks
    return ks.pop()

class SectorReducer:
    """row-reduced J_deg in one character sector mod p; normal forms over standard monomials."""
    def __init__(self, Fi, deg, key, gens, p, mons_by_key_deg, mons_by_key_degm5):
        self.p = p; self.deg = deg; self.key = key
        mons = mons_by_key_deg.get(key, [])
        self.mons = mons; self.idx = {m: i for i, m in enumerate(mons)}
        gl = []
        for k in range(6):
            kk = key_of(next(iter(Fi[k])), gens)
            need = tuple((a - b) % 6 for a, b in zip(key, kk))   # key(m) + key(F_k) = key
            for m in mons_by_key_degm5.get(need, []):
                gl.append(shift(Fi[k], m, 1))
        C = len(mons); R = len(gl)
        if R == 0 or C == 0:
            self.rank = 0; self.pivots = []; self.rows = []; self.standard = list(mons); return
        flat = [0] * (R * C)
        for r, P in enumerate(gl):
            base = r * C
            for e, v in P.items():
                flat[base + self.idx[e]] = int(v) % p
        J = nmod_mat(R, C, flat, p)
        Jr, rank = J.rref()
        self.rank = rank
        piv = []; rowsl = []
        for i in range(rank):
            j = 0
            while j < C and int(Jr[i, j]) == 0: j += 1
            piv.append(j)
            rowsl.append({mons[k]: int(Jr[i, k]) for k in range(j, C) if int(Jr[i, k])})
        self.pivots = piv; self.rows = rowsl
        ps = set(piv)
        self.standard = [m for j, m in enumerate(mons) if j not in ps]
    def normal_form(self, P):
        p = self.p
        P = {e: int(v) % p for e, v in P.items() if int(v) % p}
        for m_piv, row in zip([self.mons[j] for j in self.pivots], self.rows):
            c = P.get(m_piv, 0)
            if c:
                for e, v in row.items():
                    nv = (P.get(e, 0) - c * v) % p
                    if nv: P[e] = nv
                    else: P.pop(e, None)
        return P

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--prime', type=int, default=w2lib.P1); ap.add_argument('--deg', type=int, default=6)
    a = ap.parse_args(); p = a.prime; d = a.deg
    t0 = time.time()
    F, A, q, G = w2lib.load()
    gens = G.generators
    F2 = scale(F, 2)
    Fi = [reduce_mod(pderiv(F2, i), p) for i in range(6)]
    Ap = reduce_mod(A, p)
    kA = poly_key(A, gens)
    degN = d + 7            # deg N = d + 6q - 5 = d + 7 ; N(F2) = A D has degree d + 12
    degE = d + 12
    mb = {}
    for k in sorted({d, degN, degE, degE - 5, 2*d + 6, 2*d + 1, max(d - 5, 0)}):
        mb[k] = {}
        for m in monomials(k):
            mb[k].setdefault(key_of(m, gens), []).append(m)
    # ---- E1 kernel per sector ----
    red18 = {}
    basis = []    # list of dict(key, D, N, inJ)
    sector_info = {}
    for kD, Dmons in sorted(mb[d].items()):
        kAD = tuple((x + y) % 6 for x, y in zip(kA, kD))
        if kAD not in red18:
            red18[kAD] = SectorReducer(Fi, degE, kAD, gens, p, mb[degE], mb[degE - 5])
        R = red18[kAD]
        # images [A m] in R_degE sector
        std = R.standard; sidx = {m: i for i, m in enumerate(std)}
        imgs = []
        for m in Dmons:
            nf = R.normal_form(shift(Ap, m, 1))
            imgs.append(nf)
        # kernel of the map (coefficients on Dmons) -> R_sector
        K = len(Dmons); Cn = len(std)
        if Cn == 0:
            kervecs = [{m: 1} for m in Dmons]
        else:
            flat = [0] * (Cn * K)
            for j, nf in enumerate(imgs):
                for e, v in nf.items():
                    flat[sidx[e] * K + j] = v
            M = nmod_mat(Cn, K, flat, p)
            X, nullity = M.nullspace()
            kervecs = []
            for c in range(nullity):
                v = {Dmons[j]: int(X[j, c]) for j in range(K) if int(X[j, c])}
                kervecs.append(v)
        # J_d in this sector (for bookkeeping)
        Rd = SectorReducer(Fi, d, kD, gens, p, mb[d], mb[d - 5] if d >= 5 else {})
        sector_info[str(kD)] = {'dimS': K, 'dimK': len(kervecs), 'dimJ': Rd.rank}
        # solve N(F2) = A D for each kernel vector: N_i in sector key(A)+kD-key(F2_i)... build sector-restricted E1 system
        for v in kervecs:
            AD = reduce_mod(mul(Ap, v, p), p)
            # columns: m F2_k with key(m F2_k) = kAD
            cols = []; colp = []
            for k in range(6):
                kk = key_of(next(iter(Fi[k])), gens)
                need = tuple((x - y) % 6 for x, y in zip(kAD, kk))
                for m in mb[degE - 5].get(need, []):
                    cols.append((k, m)); colp.append(shift(Fi[k], m, 1))
            rows = {}
            for P in colp:
                for e in P: rows.setdefault(e, len(rows))
            for e in AD: rows.setdefault(e, len(rows))
            Rn, Cc = len(rows), len(cols)
            flat = [0] * (Rn * (Cc + 1))
            for j, P in enumerate(colp):
                for e, val in P.items(): flat[rows[e] * (Cc + 1) + j] = val % p
            for e, val in AD.items(): flat[rows[e] * (Cc + 1) + Cc] = val % p
            M = nmod_mat(Rn, Cc + 1, flat, p)
            Mr, rk = M.rref()
            piv = []
            for i in range(rk):
                j = 0
                while j <= Cc and int(Mr[i, j]) == 0: j += 1
                piv.append(j)
            assert piv[-1] != Cc, 'E1 kernel vector not solvable?!'
            N = [{} for _ in range(6)]
            for i, j in enumerate(piv):
                val = int(Mr[i, Cc])
                if val:
                    k, m = cols[j]; N[k][m] = val
            # check
            chk = {}
            for k in range(6): chk = add(chk, mul(N[k], Fi[k], p), p)
            assert core.sub(chk, AD, p) == {}, 'E1 solve check failed'
            inJ = (Rd.normal_form(v) == {})
            basis.append({'key': kD, 'D': v, 'N': N, 'inJ': inJ})
    nb = len(basis)
    print(f'E1 kernel K_{d}: dim {nb} (in J: {sum(b["inJ"] for b in basis)}), sectors with kernel: {sum(1 for s in sector_info.values() if s["dimK"])}  [{time.time()-t0:.1f}s]', flush=True)
    # ---- bilinear H and normal forms mod J_{2d+6} ----
    degH = 2 * d + 6
    redH = {}
    def Hpart(Da, Nb, Db, Na):
        # D_a div N_b - N_b.grad D_a + D_b div N_a - N_a.grad D_b
        out = {}
        divNb = {}; divNa = {}
        for k in range(6):
            divNb = add(divNb, pderiv(Nb[k], k, p), p); divNa = add(divNa, pderiv(Na[k], k, p), p)
        out = add(out, mul(Da, divNb, p), p); out = add(out, mul(Db, divNa, p), p)
        for k in range(6):
            out = core.sub(out, mul(Nb[k], pderiv(Da, k, p), p), p)
            out = core.sub(out, mul(Na[k], pderiv(Db, k, p), p), p)
        return out
    inv2 = pow(2, -1, p)
    quad = {}   # (key, std monomial) -> {(a,b): coeff}
    std_index = {}
    npairs = 0
    for a_ in range(nb):
        for b_ in range(a_, nb):
            Da, Na, Db, Nb = basis[a_]['D'], basis[a_]['N'], basis[b_]['D'], basis[b_]['N']
            Hab = Hpart(Da, Nb, Db, Na)
            if a_ == b_:
                Hab = scale(Hab, inv2, p)      # H(D_a,D_a) = D_a div N_a - N_a.grad D_a
            else:
                Hab = scale(Hab, inv2, p)      # symmetric bilinear value; coefficient of c_a c_b is 2*Hab
            Hab = {e: v % p for e, v in Hab.items() if v % p}
            npairs += 1
            if not Hab: continue
            kH = poly_key(Hab, gens)
            if kH not in redH:
                redH[kH] = SectorReducer(Fi, degH, kH, gens, p, mb[degH], mb[degH - 5])
                for s in redH[kH].standard:
                    std_index[(kH, s)] = len(std_index)
            nf = redH[kH].normal_form(Hab)
            for e, v in nf.items():
                coeff = v if a_ == b_ else (2 * v) % p
                quad.setdefault((kH, e), {})[(a_, b_)] = coeff
        if a_ % 10 == 0: print(f'  row {a_}/{nb} done, {len(quad)} nonzero quadrics so far [{time.time()-t0:.1f}s]', flush=True)
    # total dim of R_degH sectors touched
    dimR = sum(len(r.standard) for r in redH.values())
    print(f'quadrics: {len(quad)} nonzero (of {dimR} standard monomials in touched sectors), pairs {npairs} [{time.time()-t0:.1f}s]', flush=True)
    # ---- write msolve and Singular inputs ----
    vars_ = [f'c{i+1}' for i in range(nb)]
    polys = []
    for (kH, e), terms in quad.items():
        s = '+'.join(f'{c}*{vars_[a_]}*{vars_[b_]}' if a_ != b_ else f'{c}*{vars_[a_]}^2' for (a_, b_), c in sorted(terms.items()))
        polys.append(s)
    tag = f'hess{d}_p{p}'
    with open(f'{w2lib.PACKET}/{tag}.ms', 'w') as f:
        f.write(','.join(vars_) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')
    with open(f'{w2lib.PACKET}/{tag}.sing', 'w') as f:
        f.write(f'ring r = {p}, ({",".join(vars_)}), dp;\nideal I = ' + ',\n'.join(polys) + ';\n')
        f.write('int t = timer; ideal G = std(I); "std done in", timer - t, "sec"; "dim(G) =", dim(G); "size(G) =", size(G);\n')
        f.write('if (dim(G) == 0) { "vdim =", vdim(G); }\n')
        f.write('"leading monomials:"; lead(G);\nquit;\n')
    json.dump({'prime': p, 'deg': d, 'variables': vars_, 'basis_keys': [list(b['key']) for b in basis],
               'quadrics': [{'key': list(kH), 'std': list(e), 'terms': {f'{a_},{b_}': c for (a_, b_), c in terms.items()}} for (kH, e), terms in quad.items()]},
              open(f'{w2lib.PACKET}/{tag}.quadrics.json', 'w'))
    receipt = {'prime': p, 'deg': d, 'dimK': nb, 'dimK_inJ': sum(b['inJ'] for b in basis), 'sectors': sector_info,
               'n_quadrics': len(quad), 'dimR_touched': dimR, 'variables': nb,
               'basis': [{'key': list(b['key']), 'inJ': b['inJ'], 'D': core.to_str(b['D'])} for b in basis],
               'seconds': round(time.time() - t0, 1)}
    json.dump(receipt, open(f'{w2lib.PACKET}/{tag}.receipt.json', 'w'), indent=1)
    w2lib.log(exp=f'hess{d}-stage1', family=f'all E1-admissible D of degree {d} (dim K_{d} = {nb} mod p)', d=d, prime=p,
              result=f'{len(quad)} quadrics in {nb} variables written to {tag}.ms/.sing (stage 2 = Groebner)', rank={'dimK': nb, 'inJ': receipt['dimK_inJ']},
              seconds=receipt['seconds'], job=f'excl-cycle-q2-hess{d}', note='Hessian relaxation H(D) in J; see hess6.py docstring')
    print('done', receipt['seconds'])

if __name__ == '__main__':
    main()
