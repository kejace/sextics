"""Generate Witness.lean / Solution.lean in the accepted explicit-witness template.

Usage (python):
    from hodge import core, leangen
    files = leangen.generate('block-q2', witness_dict, name='BlockQ2')
    # files['Witness.lean'], files['Solution.lean']

The template follows the accepted proofs (e.g. data/block-q2/Witness.lean, data/m3-coupled/Witness.lean):
  * explicit MvPolynomial definitions of F, A, D, Q, N0..N5,
  * the 36 pderiv lemmas, e1 / e2 identities via MvPolynomial.funext + norm_num + ring,
  * homogeneity by isHomogeneous chains,
  * coprimality via toUni/fromUni retract + univariate specialization + Bezout identity,
  * denominator_ne_zero by evaluation at an integer point.
"""
from __future__ import annotations
import itertools, math, sys
from fractions import Fraction
from . import core
from .core import NV

VARS = core.VARS


def _coef_str(c: Fraction) -> str:
    c = Fraction(c)
    if c.denominator == 1:
        return f'{c.numerator}'
    return f'{c.numerator} / {c.denominator}'


def _term_lean(e, c):
    parts = [f'(MvPolynomial.C (σ := Fin 6) ({_coef_str(c)} : ℂ))']
    for i, k in enumerate(e):
        if k == 1:
            parts.append(f'(MvPolynomial.X (R := ℂ) ({i} : Fin 6))')
        elif k > 1:
            parts.append(f'(MvPolynomial.X (R := ℂ) ({i} : Fin 6))^({k} : ℕ)')
    return ' * '.join(parts)


def _poly_lean(P):
    if not P:
        return '0'
    terms = [_term_lean(e, c) for e, c in sorted(P.items(), reverse=True)]
    return ' +\n    '.join(terms)


def _term_homog(e, c):
    s = f'(isHomogeneous_C (σ := Fin 6) ({_coef_str(c)} : ℂ))'
    for i, k in enumerate(e):
        if k == 1:
            s = f'({s}.mul (isHomogeneous_X ℂ ({i} : Fin 6)))'
        elif k > 1:
            s = f'({s}.mul ((isHomogeneous_X ℂ ({i} : Fin 6)).pow {k}))'
    return s


def _homog_proof(P, name, deg):
    if not P:
        return f'theorem {name}_homogeneous : {name}.IsHomogeneous {deg} := by\n  simpa only [{name}] using isHomogeneous_zero (Fin 6) ℂ _\n'
    items = sorted(P.items(), reverse=True)
    s = _term_homog(*items[0])
    for e, c in items[1:]:
        s = f'({s}.add {_term_homog(e, c)})'
    return f'theorem {name}_homogeneous : {name}.IsHomogeneous {deg} := by\n  simpa only [{name}] using {s}\n'


def _deriv_lemmas():
    out = []
    for i in range(NV):
        for j in range(NV):
            val = '1' if i == j else '0'
            out.append(f'private theorem derivative_{i}_{j} : pderiv ({i} : Fin 6) (X {j} : Poly) = {val} := by norm_num [pderiv_X, Pi.single_apply, Fin.ext_iff]')
    return '\n'.join(out) + '\n'


DERIV_NAMES = ', '.join(f'derivative_{i}_{j}' for i in range(NV) for j in range(NV))


def _univariate_data(F, D, j, point):
    """Specialize all variables except j at `point` (list of 5 rationals in the order of the
    remaining variables) and return univariate polys (dict deg->Fraction) of F and D."""
    others = [i for i in range(NV) if i != j]
    def spec(P):
        out = {}
        for e, c in P.items():
            v = Fraction(c)
            for t, i in enumerate(others):
                if e[i]:
                    v *= Fraction(point[t]) ** e[i]
            if v:
                out[e[j]] = out.get(e[j], 0) + v
        return {k: v for k, v in out.items() if v}
    return spec(F), spec(D)


def _uni_to_sympy(P):
    import sympy
    x = sympy.Symbol('x')
    return sum(sympy.Rational(c.numerator, c.denominator) * x**k for k, c in P.items()), x


def _bezout(fF, fD):
    """Return (u, v) with u*fF + v*fD = 1 over Q, or None if not coprime."""
    import sympy
    f, x = _uni_to_sympy(fF)
    g, _ = _uni_to_sympy(fD)
    if g == 0:
        return None
    u, v, h = sympy.gcdex(f, g, x)
    if sympy.degree(h, x) != 0:
        return None
    h = sympy.Rational(h)  # constant
    u = sympy.Poly(sympy.expand(u / h), x)
    v = sympy.Poly(sympy.expand(v / h), x)
    return u, v


def _sympy_uni_lean(poly):
    """Polynomial ℂ expression from a sympy Poly in x."""
    import sympy
    terms = []
    for (k,), c in sorted(poly.terms(), reverse=True):
        c = sympy.Rational(c)
        cs = f'{c.p}' if c.q == 1 else f'{c.p} / {c.q}'
        if k == 0:
            terms.append(f'(Polynomial.C ({cs} : ℂ))')
        elif k == 1:
            terms.append(f'(Polynomial.C ({cs} : ℂ)) * Polynomial.X')
        else:
            terms.append(f'(Polynomial.C ({cs} : ℂ)) * Polynomial.X^{k}')
    return ' +\n    '.join(terms) if terms else '0'


def _uni_lean(P):
    """Polynomial ℂ expression from dict deg->Fraction."""
    terms = []
    for k, c in sorted(P.items(), reverse=True):
        cs = _coef_str(c)
        if k == 0:
            terms.append(f'(Polynomial.C ({cs} : ℂ))')
        elif k == 1:
            terms.append(f'(Polynomial.C ({cs} : ℂ)) * Polynomial.X')
        else:
            terms.append(f'(Polynomial.C ({cs} : ℂ)) * Polynomial.X^{k}')
    return ' +\n    '.join(terms) if terms else '0'


def _coeff_poly_lean(P, j):
    """Write P (a poly in the 5 variables other than j, as a 6-var dict with e[j]=0) as an
    MvPolynomial (Fin 5) ℂ expression using the renumbering of toUni."""
    others = [i for i in range(NV) if i != j]
    if not P:
        return '0'
    terms = []
    for e, c in sorted(P.items(), reverse=True):
        parts = [f'(MvPolynomial.C ({_coef_str(c)} : ℂ))']
        for t, i in enumerate(others):
            k = e[i]
            if k == 1:
                parts.append(f'(MvPolynomial.X {t})')
            elif k > 1:
                parts.append(f'(MvPolynomial.X {t})^{k}')
        terms.append(' * '.join(parts))
    return ' +\n    '.join(terms)


def choose_specialization(F, D, j=None, max_tries=2000, seed=0):
    """Pick variable j (F contains Z_j^6) and an integer point for the other 5 variables such
    that the univariate specializations of F and D are coprime; return (j, point, bezout)."""
    import random
    rng = random.Random(seed)
    candidates = [j] if j is not None else [i for i in range(NV)
                  if any(sum(e) == 6 and e[i] == 6 for e in F)]
    for jj in candidates:
        for trial in range(max_tries):
            if trial == 0:
                pt = [1, 1, 1, 1, 1]
            elif trial < 50:
                pt = [rng.randint(-2, 2) for _ in range(5)]
            else:
                pt = [rng.randint(-5, 5) for _ in range(5)]
            fF, fD = _univariate_data(F, D, jj, pt)
            if not fD or 6 not in fF:
                continue
            bz = _bezout(fF, fD)
            if bz is not None:
                return jj, pt, bz, fF, fD
    raise RuntimeError('no coprime specialization found')


def choose_nonzero_point(D, seed=0):
    import random
    rng = random.Random(seed)
    for trial in range(5000):
        pt = [1] * NV if trial == 0 else [rng.randint(-3, 3) for _ in range(NV)]
        if core.evaluate(D, pt) != 0:
            return pt
    raise RuntimeError('D vanishes at all tried points?')


def generate(instance_id, witness, name, heartbeats=2000000, recdepth=4000, j=None):
    """witness: dict with keys d, D, N (list of 6 strings), Q (strings in instance syntax)."""
    F, A, q, rec = core.load_instance(instance_id)
    D = core.parse(witness['D'])
    N = [core.parse(s) for s in witness['N']]
    Q = core.parse(witness['Q'])
    d = int(witness['d'])
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    if not ok:
        raise ValueError(f'witness does not verify: {info}')
    degN = d + 6 * q - 5
    jj, pt, (bu, bv), fF, fD = choose_specialization(F, D, j=j)
    others = [i for i in range(NV) if i != jj]
    nzpt = choose_nonzero_point(D)
    lc = fF[6]
    ns = f'HodgeComputer.Instances.{name}'
    L = []
    L.append('import HodgeComputer.Verify\n')
    L.append('set_option Elab.async false\n')
    L.append(f'set_option maxHeartbeats {heartbeats}\n')
    L.append(f'set_option maxRecDepth {recdepth}\n')
    L.append('set_option linter.unusedTactic false\n')
    L.append('set_option linter.unusedSimpArgs false\n')
    L.append('set_option linter.unnecessarySeqFocus false\n')
    L.append('set_option linter.unreachableTactic false\n')
    L.append('noncomputable section\n')
    L.append(f'namespace {ns}\n')
    L.append('open MvPolynomial HodgeTestbed\n')
    L.append(f'def F : Poly := {_poly_lean(F)}\n')
    L.append(f'def A : Poly := {_poly_lean(A)}\n')
    L.append(f'def D : Poly := {_poly_lean(D)}\n')
    L.append(f'def Q : Poly := {_poly_lean(Q)}\n')
    for i in range(NV):
        L.append(f'def N{i} : Poly := {_poly_lean(N[i])}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    L.append(_deriv_lemmas())
    # e1
    chain = ' + '.join(f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV))
    # nested right-assoc as produced by Fin.sum_univ_succ: a0 + (a1 + (a2 + ...))
    def nest(terms):
        return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'
    e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
    L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
             f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
             f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
             f'  simp only [F,A,D,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,{DERIV_NAMES}]\n'
             f'  apply MvPolynomial.funext; intro x\n'
             f'  norm_num <;> ring\n')
    div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
    dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
    L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
             f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
             f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
             f'  simp only [F,D,Q,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,{DERIV_NAMES}]\n'
             f'  apply MvPolynomial.funext; intro x\n'
             f'  norm_num <;> ring\n')
    L.append(_homog_proof(D, 'D', d).replace('theorem D_homogeneous', 'theorem denominator_homogeneous'))
    for i in range(NV):
        L.append(_homog_proof(N[i], f'N{i}', degN))
    L.append(f'theorem numerator_homogeneous : ∀ i, (N i).IsHomogeneous {degN} := by\n  intro i; fin_cases i\n' +
             ''.join(f'  · exact N{i}_homogeneous\n' for i in range(NV)))
    # toUni / fromUni with variable jj -> Polynomial.X, others -> C (X t)
    to_branches = []
    for i in range(NV):
        if i == jj:
            to_branches.append('(Polynomial.X)')
        else:
            to_branches.append(f'(Polynomial.C (MvPolynomial.X {others.index(i)}))')
    to_fun = ''.join(f'if i.val = {i} then {to_branches[i]} else ' for i in range(5)) + to_branches[5]
    from_branches = [f'(MvPolynomial.X {i})' for i in others]
    from_fun = ''.join(f'if i.val = {t} then {from_branches[t]} else ' for t in range(4)) + from_branches[4]
    L.append(f'def toUni : Poly →+* Polynomial Coeff := MvPolynomial.eval₂Hom (Polynomial.C.comp MvPolynomial.C) (fun i => {to_fun})\n')
    L.append(f'def fromUni : Polynomial Coeff →+* Poly := Polynomial.eval₂RingHom (MvPolynomial.eval₂Hom MvPolynomial.C (fun i => {from_fun})) (MvPolynomial.X {jj})\n')
    for i in range(NV):
        L.append(f'@[simp] theorem toUni_X{i} : toUni (MvPolynomial.X {i}) = {to_branches[i]} := by simp only [toUni,MvPolynomial.eval₂Hom_X\']; rfl\n')
    L.append('@[simp] theorem toUni_C (c : ℂ) : toUni (MvPolynomial.C c) = Polynomial.C (MvPolynomial.C c) := by simp [toUni]\n')
    for t, i in enumerate(others):
        L.append(f'@[simp] theorem fromUni_w{t} : fromUni (Polynomial.C (MvPolynomial.X {t})) = (MvPolynomial.X {i}) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X\']; rfl\n')
    L.append(f'@[simp] theorem fromUni_t : fromUni Polynomial.X = MvPolynomial.X {jj} := by simp [fromUni]\n')
    simp_from = ','.join(['fromUni_t'] + [f'fromUni_w{t}' for t in range(5)])
    L.append('theorem retract : fromUni.comp toUni = RingHom.id Poly := by\n'
             '  apply MvPolynomial.ringHom_ext\n'
             '  · intro c; simp [toUni,fromUni]\n'
             '  · intro i; fin_cases i <;> simp only [RingHom.comp_apply,RingHom.id_apply]\n' +
             ''.join(f'    · change fromUni (toUni (X {i})) = X {i}; simp only [toUni_X{i},map_add,{simp_from}] <;> ring\n' for i in range(NV)))
    # fUni: F as polynomial in variable jj with coefficients in Coeff
    by_deg = {}
    for e, c in F.items():
        k = e[jj]
        e2 = list(e); e2[jj] = 0
        by_deg.setdefault(k, {})[tuple(e2)] = c
    funi_terms = []
    for k in sorted(by_deg, reverse=True):
        funi_terms.append(f'Polynomial.C ({_coeff_poly_lean(by_deg[k], jj)}) * Polynomial.X^{k}')
    L.append('def fUni : Polynomial Coeff := ' + ' +\n    '.join(funi_terms) + '\n')
    tos = ','.join(f'toUni_X{i}' for i in range(NV))
    L.append(f'theorem toUni_F : toUni F = fUni := by\n  simp only [F, fUni, map_add, map_mul, map_pow, toUni_C, {tos}]; polynomial\n')
    L.append('theorem fUni_degree : fUni.natDegree = 6 := by\n  unfold fUni; compute_degree!\n')
    L.append(f'theorem fUni_lc : fUni.leadingCoeff = MvPolynomial.C ({_coef_str(lc)} : ℂ) := by\n'
             '  rw [Polynomial.leadingCoeff, fUni_degree]\n'
             '  simp only [fUni, Polynomial.coeff_add, Polynomial.coeff_C_mul_X_pow]\n'
             '  norm_num\n')
    L.append('def point : Fin 5 → ℂ := fun i => ' + ''.join(f'if i.val = {t} then ({pt[t]} : ℂ) else ' for t in range(4)) + f'({pt[4]} : ℂ)\n')
    L.append('def specialize : Coeff →+* ℂ := MvPolynomial.eval point\n')
    L.append(f'def fFiber : Polynomial ℂ := {_uni_lean(fF)}\n')
    L.append(f'def dFiber : Polynomial ℂ := {_uni_lean(fD)}\n')
    L.append(f'def bezoutF : Polynomial ℂ := {_sympy_uni_lean(bu)}\n')
    L.append(f'def bezoutD : Polynomial ℂ := {_sympy_uni_lean(bv)}\n')
    L.append('theorem fiber_F : (toUni F).map specialize = fFiber := by\n  rw [toUni_F]\n  norm_num [fUni,specialize,point,fFiber] <;> polynomial\n')
    L.append(f'theorem fiber_D : (toUni D).map specialize = dFiber := by\n  simp only [D,map_add,map_mul,map_pow,toUni_C,{tos}]\n  norm_num [specialize,point,dFiber] <;> polynomial\n')
    L.append('theorem fiber_coprime : IsCoprime fFiber dFiber := by\n  refine ⟨bezoutF, bezoutD, ?_⟩\n  unfold bezoutF bezoutD fFiber dFiber\n  polynomial\n')
    L.append('theorem denominator_legal : IsRelPrime F D := by\n'
             '  apply relPrime_of_retract toUni fromUni retract\n'
             '  apply relPrime_of_specialization specialize\n'
             '  · rw [toUni_F, fUni_lc]\n'
             f'    exact (isUnit_iff_ne_zero.mpr (by norm_num : ({_coef_str(lc)} : ℂ) ≠ 0)).map MvPolynomial.C\n'
             '  · rw [fiber_F, fiber_D]\n'
             '    exact fiber_coprime.isRelPrime\n')
    vec = '![' + ','.join(str(v) for v in nzpt) + ']'
    shows = ', '.join(f'show ({vec} : Fin 6 → ℂ) {i} = {nzpt[i]} from rfl' for i in range(NV))
    L.append('theorem denominator_ne_zero : D ≠ 0 := by\n'
             '  intro h\n'
             f'  have he := congrArg (MvPolynomial.eval ({vec} : Fin 6 → ℂ)) h\n'
             f'  norm_num [D, {shows}] at he\n')
    L.append(f'theorem certificate : Verify F A {q} {d} D N Q where\n'
             '  denominator_ne_zero := denominator_ne_zero\n'
             '  denominator_homogeneous := denominator_homogeneous\n'
             '  numerator_homogeneous := numerator_homogeneous\n'
             '  denominator_legal := denominator_legal\n'
             '  identity := normalized_identity (by simpa using e1) e2\n')
    L.append(f'theorem solved : HasCertificate F A {q} := ⟨{d}, D, N, Q, certificate⟩\n')
    L.append(f'end {ns}\n')
    witness_lean = '\n'.join(L)

    S = []
    S.append('import HodgeWitness\nimport Proposal.Witness\n')
    S.append('set_option maxHeartbeats 8000000\nset_option maxRecDepth 8000\nnoncomputable section\nnamespace Submission\n')
    S.append('theorem witness : HodgeComputer.Verify HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q HodgeWitness.d HodgeWitness.D HodgeWitness.N HodgeWitness.Q := by\n')
    for nm, src in (('F', 'HodgeSubmission.F'), ('A', 'HodgeSubmission.A'), ('D', 'HodgeWitness.D'), ('Q', 'HodgeWitness.Q')):
        S.append(f'  have h{nm} : {src} = {ns}.{nm} := by\n    apply MvPolynomial.funext; intro x\n    norm_num [{src}, {ns}.{nm}] <;> ring\n')
    S.append(f'  have hN : HodgeWitness.N = {ns}.N := by\n    funext i; fin_cases i\n')
    for i in range(NV):
        S.append(f'    · change HodgeWitness.N{i} = {ns}.N{i}\n      apply MvPolynomial.funext; intro x\n      norm_num [HodgeWitness.N{i}, {ns}.N{i}] <;> ring\n')
    S.append(f'  simpa only [hF,hA,hD,hQ,hN,HodgeSubmission.q,HodgeWitness.d] using {ns}.certificate\n')
    S.append('theorem certificate : HodgeComputer.HasCertificate HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q := ⟨HodgeWitness.d,HodgeWitness.D,HodgeWitness.N,HodgeWitness.Q,witness⟩\n')
    S.append('end Submission\n')
    solution_lean = ''.join(S)
    return {'Witness.lean': witness_lean, 'Solution.lean': solution_lean,
            'meta': {'j': jj, 'point': pt, 'nzpoint': nzpt}}


def witness_json(d, D, N, Q):
    return {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}
