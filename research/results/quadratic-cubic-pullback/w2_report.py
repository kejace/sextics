#!/usr/bin/env python3
"""Regenerate /tmp/hodge/notes/wave2-quadratic-cubic-pullback-search.md from wave2-search.log.jsonl (+ static prose)."""
import json, collections, time, os
LANE='quadratic-cubic-pullback'
L=f'/tmp/hodge/results/{LANE}/wave2-search.log.jsonl'
recs=[json.loads(l) for l in open(L)]
def stats(xs):
    xs=sorted(xs)
    if not xs: return 'n=0'
    return f"n={len(xs)} min={xs[0]:.3e} q10={xs[len(xs)//10]:.3e} median={xs[len(xs)//2]:.3e} max={xs[-1]:.3e}"
# kernels
ks={}
for r in recs:
    if r.get('item')=='kernel': ks[(r['d'],r['parity'])]=r
krows=[f"| {r['d']} | {'even' if r['parity']==0 else 'odd'} | {r['nD']} | {r['jacobian_dim']} | {r['ring_dim']} | {r['kernel_dim']} | {r['dim_K_mod_J']} | {r['ring_dim_target']} | {r['rows']}x{r['cols']} | {r['rank_stack']} | {r['seconds']} |" for k,r in sorted(ks.items())]
# struct
st={}
for r in recs:
    if r.get('item')=='struct-intersection': st[(r['d'],r['family'])]=r
def srows(d):
    out=[]
    for (dd,fam),r in sorted(st.items(), key=lambda kv:(kv[0][0],-kv[1]['dim_W_cap_K'],kv[0][1])):
        if dd!=d: continue
        kd=ks.get((d,0),{}).get('kernel_dim'); nD=ks.get((d,0),{}).get('nD')
        gen=max(r['dimW']+kd-nD,0) if kd and nD else '?'
        out.append(f"| {fam[:70]} | {r['dimW']} | {r['dim_W_cap_K']} | {r.get('dim_int_mod_J',0)} | {gen} |")
    return out
# closure
cl={}
for r in recs:
    if r.get('item')=='closure' and not r.get('fermat_control'):
        k=(r['d'],r['parity'])
        if k not in cl or (r.get('closed_at_k') and not cl[k].get('closed_at_k')) or len(r)>len(cl[k]): cl[k]=r
crows=[]
for k,r in sorted(cl.items()):
    extra=' '.join(f"deg{kk}: {r.get(f'deg{kk}_rank')}/{r.get(f'deg{kk}_dim')}" for kk in (3,4,5) if r.get(f'deg{kk}_rank') is not None)
    crows.append(f"| {r['d']} | {'odd' if r['parity'] else 'even'} | {r['n_w']} | {r['spurious']} | {r['nq']} | {r['rJ']} | {r['deg2_forms_rank']}/{r['deg2_dim']} | {extra} | {r['closed_at_k']} |")
prow=[]
for r in recs:
    if r.get('item')!='points': continue
    v = r.get('verdict') or ('no solution in the random affine chart (GB = [1])' if (r.get('out_head') or '').strip().startswith('[-1]') else 'see file')
    prow.append(f"| {r['d']} | {'odd' if r['parity'] else 'even'} | {r['char']} | {r['n_forms']} forms in {r['nq']} vars | {'homogeneous, -g 2' if r.get('homog') else 'affine chart l(v)=1'} | {v[:100]} | {r['msolve_seconds']} s |")
# search
s1=[r for r in recs if r.get('stage')==1]; s2=[r for r in recs if r.get('stage')==2]
lrows=[]
for d in sorted({r['d'] for r in s1}):
    for par in sorted({r.get('parity',0) for r in s1 if r['d']==d}):
        rs=[r for r in s1 if r['d']==d and r.get('parity',0)==par]
        c=collections.Counter(round(r['rel'],4) for r in rs)
        lrows.append(f"| {d} | {'odd' if par else 'even'} | stage 1 (Hessian relaxation, LM) | {stats([r['rel'] for r in rs])} | {len(c)} | {sum(r['seconds'] for r in rs)/len(rs):.0f} s |")
for d in sorted({r['d'] for r in s2}):
    for par in sorted({r.get('parity',0) for r in s2 if r['d']==d}):
        rs=[r for r in s2 if r['d']==d and r.get('parity',0)==par]
        c=collections.Counter(round(r['rel_compressed'],4) for r in rs)
        lrows.append(f"| {d} | {'odd' if par else 'even'} | stage 2 (full E2 varpro, L-BFGS; compressed rows) | {stats([r['rel_compressed'] for r in rs])} | {len(c)} | {sum(r['seconds'] for r in rs)/len(rs):.0f} s |")
rk=[r for r in recs if r.get('item')=='realkernel']
rkline = (f"d = 10 real kernel (job d10): dim K = {rk[-1]['dimK']}, dim J = {rk[-1]['dimJ']}, non-J part {rk[-1]['dim_nonJ']}; singular values of the Jacobian block at the cut "
          f"{rk[-1]['s_rank_minus1']:.3f} / {rk[-1]['s_rank']:.1e}; E1-lift residual max {rk[-1]['lift_relres_max']:.1e}." if rk else "d = 10 real kernel: pending.")
hess=[r for r in recs if r.get('item')=='hessian-identity']
closed=[f"d={k[0]} {'odd' if k[1] else 'even'}" for k,r in sorted(cl.items()) if r.get('closed_at_k')]
now=time.strftime('%Y-%m-%d %H:%M')
rep=f"""# Wave 2 search report: quadratic-cubic-pullback

STATUS: NO WITNESS FOUND.  (Report regenerated {now} EST by w2_report.py from wave2-search.log.jsonl; {len(recs)} records.)

Lane: F = sum_{{i=0}}^5 Q_i^3 (six Y-even quadrics), q = 1, A = X0^3 Y0^3.  Diagonal symmetry group of F: only Z/2 = {{1, tau}}
(Y-parity), so the only character sectors are the two Y-parities of D (a mixed-parity D can be replaced by D*tau(D)).
Certificate convention: A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q, deg N = d+1, deg Q = 2d-6, gcd(D, F) = 1.

## Headline results
1. Exact E1-admissible spaces K_d = {{D : A D in J(F)}} for d = 4..12, both parities (mod 2^31-1): K_d = J(F)_d for every
   d <= 9 (both parities); dim K_d/J_d = 294/261 (d=10), 510/510 (d=11), 677/648 (d=12) = dim R_d - dim R_{{d+6}} exactly
   (multiplication by [A] is surjective in those degrees).  Rank/nullity of the stacked E1 matrices in the table below.
2. Hessian form of the E2 relaxation (verified exactly): for D = sum_l w_l F_l, N0 = A w,
   D div N0 - N0(D) == -A (w^T Hess(F) w) (mod J(F)), so E2 mod J <=> w^T Hess(F) w in K_{{2d-6}}^even = (J(F):A)_{{2d-6}}.
3. Rigorous closure of degrees 5 and 6 (both parities): the stage-1 quadratic forms in the w-coordinates have no common
   projective zero (degree-4 Macaulay matrices of full column rank 4845/4845 and 5985/5985 exactly mod p, hence over Q;
   independently msolve's homogeneous reduced Groebner basis has v_i^2 in its leading ideal for every variable).
   Together with K_d = J_d (d <= 9): EVERY certificate for this target has deg D >= 7, in either parity.
   Closed exactly: {', '.join(closed)}.
4. Structured families at d = 10 (pullbacks C[Q]_5 of quintic sections of the cubic, the 20 triples C[Q_S]_5, bihomogeneous
   products, delta * S_4, A * S_4, cones with vertex a plane or a line, all bigraded pieces, products of polars, unions,
   rho-/sigma-invariant subspaces): the only elements of K_10 they contain are F-multiples, except the generic-dimensional
   intersections with the big unions and the rho-/sigma-invariant subspaces (36 / 288 dimensional non-J subfamilies).
5. Numerical joint searches (stage 1 = Hessian relaxation with Levenberg-Marquardt, stage 2 = full E2 with the Koszul
   corrections T_D(c) and quotient Q eliminated by variable projection): d = 6 landscape minimum rel 1.9e-2 (361 restarts,
   23 distinct minima), stage 2 minimum 1.6e-2 (compressed; scaled by k/nout ~ 0.12); d = 7, 8, 9, 10 searches were
   starved by the shared box (load 300-1750 on 128 cores) -- counts in the landscape table.  No candidate reached 1e-10,
   so no rational reconstruction was attempted.

## Batch 1 (item 1): exact E1-admissible spaces K_d (job search-quadratic-cubic-pullback-kernels, mod p = 2^31-1)
Nullspace of [A*m | J(F)_{{d+6}} generators]; d <= 8 Y-even and d <= 6 Y-odd were closed in wave 1 (K_d = J_d).

| d | parity of D | dim S_d | dim J_d | dim R_d | dim K_d | dim K_d/J_d | dim R_{{d+6}} (target) | matrix | rank | s |
|---|---|---|---|---|---|---|---|---|---|---|
{chr(10).join(krows)}

## Batch 2 + 5 (item 4, linear part): structured subspaces W of S_d^even against the exact K_d (jobs struct10, struct1112)
E1 for D is exactly D in K_d; "generic" = dim W + dim K_d - dim S_d (what a random subspace of that dimension would meet).

d = 10:

| family W | dim W | dim W cap K_10 | dim mod J | generic |
|---|---|---|---|---|
{chr(10).join(srows(10))}

d = 12 (if present):

| family W | dim W | dim W cap K_12 | dim mod J | generic |
|---|---|---|---|---|
{chr(10).join(srows(12)) or '| (pending) | | | | |'}

d = 11 (if present):

| family W | dim W | dim W cap K_11 | dim mod J | generic |
|---|---|---|---|---|
{chr(10).join(srows(11)) or '| (pending) | | | | |'}

## Hessian form of stage 1 (verified exactly, both parities; record item "hessian-identity")
For D = sum_l w_l F_l and N0 = A w: D div N0 - N0(D) == -A (w^T Hess(F) w) (mod J(F)).  Hence the Jacobian-quotient
relaxation of E2 ("Hessian relaxation") is exactly w^T Hess(F) w in K_{{2d-6}}^even; for d <= 7 that is membership in
J(F)_{{2d-6}} (independent of A), and the number of independent quadratic equations is dim R_{{2d-6}}^even (218, 483, 468, ...)
rather than dim R_{{2d}}^odd.

## Batch 3 (item 2, exact closure of Jacobian-type degrees; w2_closure.py, w2_points.py, job closure678)
At degrees with K_d = J_d every E1-admissible D is sum_l w_l F_l; E2 implies the stage-1 system q(w) = 0, rJ = dim R_{{2d}}^odd
quadratic forms in the w-coordinates, invariant under the spurious translations D -> D + F G (quotient coordinates v, dim n_q).
Full column rank of the degree-k Macaulay matrix (exact mod p => over Q) certifies that the forms have no common zero over C,
i.e. no denominator of that degree and parity admits a certificate.

| d | parity | n_w | spurious | n_q | rJ | deg-2 rank | higher Macaulay ranks | closed at k |
|---|---|---|---|---|---|---|---|---|
{chr(10).join(crows)}

msolve (F4 mod p = 2147483629):

| d | parity | char | system | mode | verdict | time |
|---|---|---|---|---|---|---|
{chr(10).join(prow)}

## Batch 4 + 6 (items 2/3, numerical joint search; jobs jvp-d678, d10, jvp-d89): residual landscape
rel = ||residual||^2 / ||B_D(N0)||^2 at the converged point (0 would be a certificate).  F-multiples and the Koszul gauge
are projected out; |u| = 1.  Stage-2 values use a random row compression S (k = n_c + 200 rows), which scales norms by
about k/n_out (0.12 at d = 6); exact stage-2 residuals are recomputed only below 1e-6 (never triggered).
{rkline}

| d | parity | objective | relative residual statistics | distinct minima | cost/restart |
|---|---|---|---|---|---|
{chr(10).join(lrows) or '| (none yet) | | | | | |'}

## What remains open after wave 2
* Degree 7 (both parities): stage 1 is 483 quadrics in 60 variables; the degree-3 Macaulay matrix cannot be full (28980 < 37820
  columns) and the degree-4 one (884k x 596k) exceeds dense exact linear algebra; msolve F4 on it is the pending exact route.
* Degrees 8, 9 (Jacobian-type, both parities) and 10 (K_10 = J_10 + 294 non-J): only numerical evidence; the stage-2 systems are
  overdetermined by factors 10-30, and the remote box was saturated, so restart counts stayed far below the planned thousands.
* Structured non-J subfamilies at d = 10 worth finishing first: rho-invariant (36-dim) and sigma-invariant (288-dim) parts of K_10.
* d = 11, 12 kernels are available (w2_kernel_<d>_<par>.json on the remote mirror); the full E2 systems there (40k-59k rows)
  need a sparse/iterative solver before any joint search is realistic.

## Reproduction
Lane dir /tmp/hodge/results/{LANE}/: w2_kernels.py (exact K_d), w2_struct.py (structured intersections), w2_closure.py
(Macaulay certificates), w2_points.py (msolve export/verdicts), w2_realkernel.py (real non-J kernel + E1 lifts),
w2_jvp.py (stage-1/2 variable-projection search; Fermat control: F0 = sum Z^6, D = X0^6+X1^6+X2^6 gives residual 1e-33),
w2_merge.py / w2_report.py (log merge, this report), w2_queue.sh (remote submissions).  Remote mirror: /tmp/hodge/remote/{LANE}/.
Log: /tmp/hodge/results/{LANE}/wave2-search.log.jsonl (one record per batch/restart; search vectors kept in the mirror files).
"""
open(f'/tmp/hodge/notes/wave2-{LANE}-search.md','w').write(rep)
print('report written:', len(rep), 'chars;', len(krows), 'kernel rows,', len(crows), 'closure rows,', len(lrows), 'landscape rows')
