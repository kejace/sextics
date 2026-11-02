"""(a) Adversarial re-verification of the rigorous claims of the dossier / wave-1 report.
Second prime(s), independent code paths, exact rational membership, age table, automorphisms.
Writes receipts a_*.json and log lines (experiment='a:*')."""
import sys, time, json, itertools, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual/wave2-excl')
from xlib import *
from hodge import jacobian
import numpy as np
from flint import nmod_mat

F, A, q = load()
t00 = time.time()
def say(*a):
    print(f'[{time.time()-t00:7.1f}s]', *a, flush=True)

# ---------- (1) structure of F ----------
e2X, e2Y = e_k(X, 2), e_k(Y, 2)
Fmodel = core.add(core.add(p_k(X, 6), p_k(Y, 6)),
                  core.scale(core.mul(core.mul(e2X, e2Y), core.add(e2X, core.scale(e2Y, Fraction(2)))), Fraction(5, 7)))
ok1 = (core.sub(F, Fmodel) == {})
say('F == p6X+p6Y+(5/7)e2X e2Y (e2X+2e2Y):', ok1)
log(experiment='a:structure', ok=ok1, statement='F = p6(X)+p6(Y)+(5/7) e2X e2Y (e2X+2 e2Y) exactly over Q')
# ---------- (2) smooth reduction mod 5 ----------
F5 = core.reduce_mod(F, 5)
ferm = {e: c % 5 for e, c in core.add(p_k(X, 6), p_k(Y, 6)).items()}
ok2 = (F5 == ferm)
say('F mod 5 == Fermat:', ok2)
log(experiment='a:smooth-mod5', ok=ok2, statement='F mod 5 is the Fermat sextic; J mod 5 = (z_i^5) so X is smooth over Q')
# ---------- (3) sector dimensions in R_6^{odd} ----------
sec = {}
for p in (P2, P3):
    mons6 = core.monomials(6, 1)
    # invariant orbit sums (polynomial character triv x triv)
    reps = {}
    for m in mons6:
        key = (tuple(sorted(m[:3], reverse=True)), tuple(sorted(m[3:], reverse=True)))
        reps.setdefault(key, []).append(m)
    inv = [{m: 1 for m in orb} for orb in reps.values()]
    # sgn x sgn polynomial part: Delta_X Delta_Y
    DX = core.mul(core.mul(core.sub(X[0], X[1]), core.sub(X[0], X[2])), core.sub(X[1], X[2]))
    DY = core.mul(core.mul(core.sub(Y[0], Y[1]), core.sub(Y[0], Y[2])), core.sub(Y[1], Y[2]))
    DXDY = core.mul(DX, DY)
    MJ, _, nJ = jac_matrix(F, 6, 1, p)
    rJ = MJ.rank()
    Minv, _, _ = jac_matrix(F, 6, 1, p, extra_cols=inv)
    Msg, _, _ = jac_matrix(F, 6, 1, p, extra_cols=[DXDY])
    MA, _, _ = jac_matrix(F, 6, 1, p, extra_cols=[A])
    sec[p] = dict(rankJ6odd=rJ, ci=ci_jac_rank(6, 1), n_inv_orbits=len(inv), dim_inv_sector=Minv.rank() - rJ,
                  dim_sgnsgn_poly_sector=Msg.rank() - rJ, A_in_J=(MA.rank() == rJ))
    say('sector dims mod', p, sec[p])
log(experiment='a:sector-dims', result=sec,
    statement='polynomial-invariant Y-odd sextics mod J (= residue sector sgn x sgn containing [A]) has dim 17; polynomial sgn x sgn part (Delta_X Delta_Y, residue-invariant) has dim 1')
receipt('a_sector_dims', sec)

# ---------- (4) E1 kernel d=1..8 via engine jacobian.mult_kernel, primes P2, P3 ----------
ker = {}
for p in (P2, P3):
    for d in range(1, 9):
        for par in (0, 1):
            r = jacobian.mult_kernel(F, A, d, p=p, parity=par)
            rec = dict(p=p, d=d, par=par, dim_ker_mod_J=r['dim'], kernel_dim=r['kernel_dim'], jacobian_dim=r['jacobian_dim'],
                       ring_dim=r['ring_dim'], ci_ring_dim=ci_ring_dim(d, par), rows=r['rows'], cols=r['cols'], seconds=r['t_total'])
            # rank of J_{d+6} implied: cols = K + nJ; nullity = K + nJ - rank(M); rank(M) = rankJ + K - kernel_dim
            ker[(p, d, par)] = rec
            say('mult_kernel', rec)
            log(experiment='a:e1-kernel', **rec)
            if r['dim'] and d == 7:
                with open(f'{PK}/a_e1_kernel_d7_p{par}_{p}.pkl', 'wb') as f:
                    pickle.dump(r['basis'], f)
receipt('a_e1_kernel', {f'{k[0]}_{k[1]}_{k[2]}': v for k, v in ker.items()})

# ---------- (5) rank bookkeeping: rank J_k^{par} mod p == CI value, k = 5..14 ----------
rk = {}
for p in (P2, P3):
    for k in range(5, 15):
        for par in (0, 1):
            t0 = time.time()
            M, mons, nJ = jac_matrix(F, k, par, p)
            r = M.rank()
            rec = dict(p=p, k=k, par=par, rows=len(mons), cols=nJ, rank=r, ci=ci_jac_rank(k, par), ok=(r == ci_jac_rank(k, par)), seconds=round(time.time() - t0, 1))
            rk[(p, k, par)] = rec
            say('rankJ', rec)
            log(experiment='a:rankJ', **rec)
        receipt('a_rankJ', {f'{k[0]}_{k[1]}_{k[2]}': v for k, v in rk.items()})

# ---------- (6) exact rational membership A*D7 in J_13, A*D7Y in J_13; D7 not in J_7 ----------
DX = core.mul(core.mul(core.sub(X[0], X[1]), core.sub(X[0], X[2])), core.sub(X[1], X[2]))
DY = core.mul(core.mul(core.sub(Y[0], Y[1]), core.sub(Y[0], Y[2])), core.sub(Y[1], Y[2]))
X2 = [core.power(v, 2) for v in X]; Y2 = [core.power(v, 2) for v in Y]
D7 = core.mul(DX, e_k(X2, 2)); D7Y = core.mul(DY, e_k(Y2, 2))
def exact_membership(target, k, par, label):
    """CRT-lift a particular solution of J g = target and verify exactly."""
    Fi = [core.pderiv(F, i) for i in range(6)]
    # column bookkeeping identical to jac_matrix
    colinfo = []
    for i in range(6):
        par_m = par if i < 3 else (par + 1) % 2
        for m in core.monomials(k - 5, par_m):
            colinfo.append((i, m))
    sols = []; mods = []
    pivcols_ref = None
    for p in core.PRIMES[:12]:
        M, mons, nJ = jac_matrix(F, k, par, p, extra_cols=[target])
        R, rk_ = M.rref()
        piv = []
        for i in range(rk_):
            j = 0
            while j < M.ncols() and int(R[i, j]) == 0: j += 1
            piv.append(j)
        if piv and piv[-1] == nJ:
            return dict(ok=False, note='inconsistent mod %d' % p)
        if pivcols_ref is None: pivcols_ref = piv
        elif piv != pivcols_ref:
            say('pivot set differs mod', p, '- skipping prime'); continue
        g = {}
        for i, j in enumerate(piv):
            g[j] = int(R[i, nJ])
        sols.append(g); mods.append(p)
        # try reconstruction
        if len(sols) >= 2:
            Mm = 1; comb = {}
            for g, p_ in zip(sols, mods):
                for j, v in g.items():
                    comb[j] = core.crt_pair(comb.get(j, 0), Mm, v, p_) if Mm > 1 else v
                Mm *= p_
            ratio = {}
            okr = True
            for j, v in comb.items():
                fr = core.ratrecon(v % Mm, Mm)
                if fr is None: okr = False; break
                ratio[j] = fr
            if okr:
                # exact verification
                tot = {}
                for j, fr in ratio.items():
                    i, m = colinfo[j]
                    tot = core.add(tot, core.shift(Fi[i], m, fr))
                if core.sub(tot, target) == {}:
                    hmax = max(max(abs(fr.numerator), fr.denominator) for fr in ratio.values())
                    gpolys = {}
                    for j, fr in ratio.items():
                        i, m = colinfo[j]
                        gpolys.setdefault(i, {})
                        gpolys[i] = core.add(gpolys[i], {m: fr})
                    return dict(ok=True, primes_used=mods, n_nonzero=len(ratio), height_max=str(hmax),
                                g=[core.to_str(gpolys.get(i, {})) for i in range(6)])
    return dict(ok=False, note='reconstruction failed within 12 primes')
for target, k, par, label in ((core.mul(A, D7), 13, 0, 'A*D7'), (core.mul(A, D7Y), 13, 1, 'A*D7Y')):
    t0 = time.time()
    res = exact_membership(target, k, par, label)
    res['seconds'] = round(time.time() - t0, 1)
    say(label, 'in J exactly:', res.get('ok'), res.get('note', ''), 'height', res.get('height_max'))
    log(experiment='a:exact-membership', target=label, ok=res['ok'], primes=res.get('primes_used'), height_max=res.get('height_max'), seconds=res['seconds'])
    receipt(f'a_exact_{label.replace("*","_")}', res)
# D7, D7Y not in J_7 (rank argument), and compare with wave-1 kernel representative
for p in (P2, P3):
    for Dk, par, label in ((D7, 0, 'D7'), (D7Y, 1, 'D7Y')):
        MJ, _, nJ = jac_matrix(F, 7, par, p); M2, _, _ = jac_matrix(F, 7, par, p, extra_cols=[Dk])
        rec = dict(p=p, D=label, rankJ7=MJ.rank(), rankJ7_plus_D=M2.rank(), ci=ci_jac_rank(7, par))
        say(rec); log(experiment='a:D7-not-in-J', **rec)

# ---------- (7) age table over the 72-element projective group ----------
def ages():
    out = []
    Fc = {e: float(c) for e, c in F.items()}
    def Fev(v):
        return sum(c * np.prod([v[i] ** e[i] for i in range(6)]) for e, c in Fc.items())
    for perm, sx, sy in G_elements():
        for eps in (1, -1):
            g = np.zeros((6, 6), dtype=complex)
            for i in range(6):
                g[perm[i], i] = 1 if i < 3 else eps
            chi_res = sx * sy
            w, V = np.linalg.eig(g)
            ang = [Fraction(int(round(np.angle(l) / (2 * np.pi) * 6)) % 6, 6) for l in w]
            comps = {}
            for a in set(ang):
                idx = [i for i in range(6) if ang[i] == a]
                comps[a] = idx
            rows = []
            for a, idx in comps.items():
                Vl = V[:, idx]
                # orthonormalise
                Qm, _ = np.linalg.qr(Vl)
                dim = Qm.shape[1]
                rng = np.random.default_rng(1)
                if dim >= 2:
                    vals = [abs(Fev(Qm @ (rng.normal(size=dim) + 1j * rng.normal(size=dim)))) for _ in range(3)]
                    meets = True
                    contained = max(vals) < 1e-9
                else:
                    v = Qm[:, 0]
                    meets = abs(Fev(v)) < 1e-9
                    contained = meets
                age = sum(Fraction(((b - a) % 1)) for b in ang)   # sum over all 6 eigenvalues of frac(angle(mu/lambda)); normal eigenvalue lambda^-6 = 1
                rows.append(dict(eig_angle=str(a), dim=dim, meets_X=bool(meets), contained=bool(contained), age=str(age)))
            out.append(dict(perm=perm, eps=eps, sgnX=sx, sgnY=sy, chi_res=chi_res, is_identity=(perm == tuple(range(6)) and eps == 1), comps=rows))
    return out
ag = ages()
minage = None
for g in ag:
    if g['is_identity'] or g['chi_res'] != 1: continue
    for c in g['comps']:
        if c['meets_X']:
            a = Fraction(c['age'])
            if minage is None or a < minage: minage = a
n_res_inv = sum(1 for g in ag if g['chi_res'] == 1)
say('age screen: elements with chi_res=+1:', n_res_inv, 'min age at a fixed point on X among nontrivial ones:', minage)
log(experiment='a:ages', n_elements=len(ag), n_residue_invariant=n_res_inv, min_age_nontrivial=str(minage),
    statement='every nontrivial element of ker(chi_res) (order 36) has age >= 1 at every fixed point on X (Reid-Tai: canonical quotients, Kollar-Larsen: not uniruled)')
receipt('a_ages', ag)

# ---------- (8) signed-permutation automorphisms of F ----------
auts = []
for perm in itertools.permutations(range(6)):
    for sg in itertools.product((1, -1), repeat=6):
        if sg[0] == -1: continue   # projective: normalise first sign
        G_ = signs(permute(F, perm), sg)
        if set(G_.keys()) != set(F.keys()): continue
        ratios = {G_[e] / F[e] for e in F}
        if len(ratios) == 1:
            auts.append((perm, sg, str(ratios.pop())))
block_pres = sum(1 for perm, sg, c in auts if set(perm[:3]) == {0, 1, 2})
say('signed-permutation automorphisms (projective):', len(auts), 'block-preserving:', block_pres)
log(experiment='a:monomial-automorphisms', n_projective=len(auts), n_block_preserving=block_pres,
    statement='signed permutation matrices g with F(gz) = c F(z): exactly 72 projectively = S3 x S3 x <sigma>, no block swap')
receipt('a_automorphisms', auts)
say('DONE')
