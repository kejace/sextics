"""Rewrite /tmp/hodge/notes/wave2-klein-coupled-character-search.md from wave2-search.log.jsonl (STATUS first)."""
import json, os, sys, time, collections
LANE = '/tmp/hodge/results/klein-coupled-character'
LOG = LANE + '/wave2-search.log.jsonl'
OUT = '/tmp/hodge/notes/wave2-klein-coupled-character-search.md'


def load():
    rows = []
    if os.path.exists(LOG):
        for line in open(LOG):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except Exception:
                pass
    return rows


def fmt(x, nd=2):
    if x is None:
        return '-'
    try:
        return f'{x:.{nd}e}'
    except Exception:
        return str(x)


def main():
    rows = load()
    witness = [r for r in rows if r.get('family') == 'WITNESS']
    status = 'WITNESS FOUND -- see /tmp/hodge/results/klein-coupled-character/witness-candidate.json' if witness else 'NO WITNESS FOUND (search in progress)'
    if len(sys.argv) > 1 and sys.argv[1] == 'final' and not witness:
        status = 'NO WITNESS FOUND (wave 2 complete)'
    L = []
    L.append('# Wave 2 search report: klein-coupled-character (q = 1, A = X0^3 Y0^3)\n')
    L.append(f'STATUS: {status}\n')
    L.append(f'Updated {time.strftime("%Y-%m-%d %H:%M:%S")} from {len(rows)} log lines in `{LOG}`.\n')
    L.append('''
Code (lane dir `/tmp/hodge/results/klein-coupled-character/`): `w2lib.py` (sector-restricted Koszul variable projection,
sector E1 spaces with rational lifts, mod-J relaxation), drivers `w2_kernels.py`, `w2_hess.py`, `w2_varpro.py`,
`w2_struct.py`, control `w2_control.py`, queue `w2_queue.py` (remote jobs `w2-*`), this generator `w2_report.py`.
Conventions: (x,y,z,t,r,s) = (X0..Y2); sigma-weights W = (1,2,4,0,0,0) mod 7; "sector (d, parD, k)" = Y-parity parD
denominators of degree d and sigma-weight k; A has weight 3, Y-odd.  "rel" = squared residual / squared right-hand side.

## Method
* **Sector E1 spaces.** (J(F):A)_d per sector = J(F)_d ⊕ (kernel of [A] on (R/J)_d) computed mod p = 2147483647 with
  canonical bases, lifted to Q by 12-prime CRT + LLL (small-height representatives), exact E1 lifts N_k by CRT.
* **Joint E2 search (varpro, Koszul form).** D(a) = sum a_k D_k, N0(a) = sum a_k N_k, unknowns p_ij (deg d-4) and Q
  (deg 2d-6) character-restricted; Q eliminated by the t-monic remainder mod F.  f(a) = min_p |R(T_{D(a)} p + B_{D(a)}(N0(a)))|^2
  on the unit sphere orthogonal to the F-multiples.  Optimizers: L-BFGS (engine style), LM/TRF with the Kaufman
  Jacobian, joint LM in (a,x), combined (E2 + relaxation) LM.  Candidates with rel < 1e-10 are reconstructed in the
  canonical rem_F chart with denominators up to 1e6 and re-solved exactly (kcc2 Koszul test mod p, CRT lift, verify_witness).
* **Mod-J (Hessian) relaxation.** Necessary condition B_{D(a)}(N0(a)) in J(F)_{2d} (for Jacobian-type D: A (w^T Hess F w) in J).
  Real zeros searched numerically (L-BFGS, hundreds-thousands of restarts; dimension <= 45 for d <= 8 so the global
  minimum is reliably found); the exact quadric system mod 32003 in a rational chart complementary to the F-multiples
  is passed to Singular: `dim = 0` (only the origin) is a **rigorous exclusion** of every sigma-eigen denominator of that
  (d, parity, weight) -- mod-p emptiness of a homogeneous system implies emptiness over Q (closure of the generic fibre).
  `dim > 0` mod p does not imply real/rational zeros (the (6,even,0) component is an octic curve with no real points found).
''')
    # ---------------- sector dims
    dims = [r for r in rows if r.get('family') == 'sector-dims']
    if dims:
        L.append('\n## (1) E1-admissible sector dimensions (kernel of [A] on (R/J)_d^{parity}[k], mod p; upper bounds over Q)\n')
        L.append('| d | parity | nD per k=0..6 | J_d per k | kernel mod J per k (nontrivial) | rational reps (max height) | lifts failed |')
        L.append('|---|---|---|---|---|---|---|')
        by = collections.defaultdict(dict)
        for r in dims:
            by[(r['d'], r['parD'])][r['k']] = r
        for (d, parD) in sorted(by):
            ks = by[(d, parD)]
            def col(key):
                return ','.join(str(ks[k].get(key, '-')) if k in ks else '-' for k in range(7))
            heights = ','.join((str(max(ks[k]['rep_heights'])) if ks[k].get('rep_heights') else '-') if k in ks else '-' for k in range(7))
            lf = sum(ks[k].get('lift_fail', 0) or 0 for k in ks)
            L.append(f'| {d} | {"even" if parD == 0 else "odd"} | {col("nD")} | {col("J_d")} | {col("nontrivial")} | {heights} | {lf} |')
    # ---------------- relaxation
    hess = [r for r in rows if str(r.get('family', '')).startswith('hess-relax') and 'error' not in r]
    if hess:
        L.append('\n## (2) Mod-J relaxation per sector (real zeros numerically; Singular exact dimension mod 32003)\n')
        L.append('| d | parity | k | family | unknowns (dimP) | equations | restarts | best rel | median rel | #rel<1e-10 | Singular dim (vars, eqs, s) |')
        L.append('|---|---|---|---|---|---|---|---|---|---|---|')
        best = {}
        for r in hess:
            key = (r['d'], r['parD'], r['k'], r['family'])
            if key not in best or (r.get('restarts', 0) or 0) >= (best[key].get('restarts', 0) or 0):
                # merge singular info if missing
                sing = r.get('singular') or (best.get(key, {}).get('singular'))
                best[key] = dict(r); best[key]['singular'] = sing
            elif r.get('singular') and not best[key].get('singular'):
                best[key]['singular'] = r['singular']
        for key in sorted(best):
            r = best[key]
            s = r.get('singular') or {}
            sd = s.get('dim')
            stxt = '-' if not s else (f'**{sd}**' if sd == 0 else (str(sd) if sd is not None else s.get('singular', 'timeout') if isinstance(s.get('singular'), str) else 'timeout/err'))
            if s:
                stxt += f' ({s.get("nvars", "?")}, {s.get("neq", "?")}, {s.get("seconds", "?")}s)'
            L.append(f'| {r["d"]} | {"even" if r["parD"] == 0 else "odd"} | {r["k"]} | {r["family"].replace("hess-relax-", "")} | {r.get("dimP")} | {r.get("n_equations")} | {r.get("restarts")} | {fmt(r.get("best_rel"))} | {fmt(r.get("median_rel"))} | {r.get("n_lt_1e-10")} | {stxt} |')
        L.append('\nEmpty sectors (no E1-admissible eigen denominator): ' + ', '.join(f'({r["d"]},{"even" if r["parD"]==0 else "odd"},{r["k"]})' for r in rows if str(r.get('family','')).startswith('hess-relax') and r.get('K') == 0) + '\n')
    # ---------------- varpro
    vp = [r for r in rows if str(r.get('family', '')).endswith('-varpro') and 'error' not in r and r.get('K')]
    if vp:
        L.append('\n## (3) Full-E2 joint searches (sector-restricted Koszul varpro)\n')
        L.append('| d | parity | k | family | K | dimP | rows x cols (p,Q elim.) | restarts | best rel | q25 | median | #<1e-6 | #<1e-10 | methods | candidates |')
        L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
        agg = {}
        for r in vp:
            key = (r['d'], r['parD'], r['k'], r['family'])
            a = agg.setdefault(key, {'restarts': 0, 'best': 1e9, 'q25': [], 'med': [], 'n6': 0, 'n10': 0, 'methods': set(), 'cands': 0, 'r': r})
            a['restarts'] += r.get('restarts', 0) or 0
            a['best'] = min(a['best'], r.get('best_rel', 1e9) if r.get('best_rel') is not None else 1e9)
            a['q25'].append(r.get('q25_rel')); a['med'].append(r.get('median_rel'))
            a['n6'] += r.get('n_lt_1e-6', 0) or 0; a['n10'] += r.get('n_lt_1e-10', 0) or 0
            a['methods'].add(r.get('method')); a['cands'] += len(r.get('candidates') or [])
        for key in sorted(agg):
            a = agg[key]; r = a['r']
            q25 = min(x for x in a['q25'] if x is not None) if any(x is not None for x in a['q25']) else None
            med = sorted(x for x in a['med'] if x is not None)
            med = med[len(med) // 2] if med else None
            L.append(f'| {r["d"]} | {"even" if r["parD"] == 0 else "odd"} | {r["k"]} | {r["family"].replace("-varpro", "")} | {r.get("K")} | {r.get("dimP")} | {r.get("rows")}x{r.get("cols")} | {a["restarts"]} | {fmt(a["best"])} | {fmt(q25)} | {fmt(med)} | {a["n6"]} | {a["n10"]} | {",".join(sorted(m for m in a["methods"] if m))} | {a["cands"]} |')
    # ---------------- candidates
    cands = [r for r in rows if str(r.get('family', '')).startswith('candidate-confirm')]
    if cands:
        L.append('\n## (5) Exact re-solves of numerical candidates (rel < 1e-10)\n')
        L.append('| d | parity | k | rel | max_den | status | rank | nullity | terms | D (truncated) |')
        L.append('|---|---|---|---|---|---|---|---|---|---|')
        for r in cands:
            L.append(f'| {r["d"]} | {r["parD"]} | {r["k"]} | {fmt(r.get("rel"))} | {r.get("max_den")} | {r.get("status")} | {r.get("rank")} | {r.get("nullity")} | {r.get("terms")} | `{str(r.get("D", ""))[:80]}` |')
    # ---------------- structured
    st0 = [r for r in rows if str(r.get('family', '')).startswith('struct:') and 'error' not in r]
    dedup = {}
    for r in st0:
        dedup[(r['family'], r['d'], r['parD'], r['k'])] = r
    st = list(dedup.values())
    if st:
        L.append('\n## (4) Structured families (E1-admissible dimension per sector; varpro where nonzero)\n')
        nz = [r for r in st if r.get('dim_E1')]
        zero_by_fam = collections.Counter(r['family'].replace('struct:', '') for r in st if not r.get('dim_E1'))
        L.append('Cells with zero E1-admissible part (excluded as supports), by family: ' + ', '.join(f'{f}:{n}' for f, n in sorted(zero_by_fam.items())) + '\n')
        L.append(f'{len(st)} (family, d, parity, k) cells tested; {len(nz)} with a nonzero E1-admissible part (listed), the rest are excluded as supports (mod-p kernel zero).\n')
        if nz:
            L.append('| family | d | parity | k | dim V | dim E1 | nontrivial mod J | dimP | restarts | best rel | exact tests |')
            L.append('|---|---|---|---|---|---|---|---|---|---|---|')
            for r in sorted(nz, key=lambda r: (r['d'], r['parD'], r['k'], r['family'])):
                ex = r.get('reps_exact')
                extxt = ','.join(e.get('status', '?') for e in ex) if ex else '-'
                L.append(f'| {r["family"].replace("struct:", "")} | {r["d"]} | {r["parD"]} | {r["k"]} | {r.get("dim_V")} | {r.get("dim_E1")} | {r.get("nontrivial_mod_J")} | {r.get("dimP", "-")} | {r.get("restarts", "-")} | {fmt(r.get("best_rel"))} | {extxt} |')
        exps = [r for r in rows if r.get('family') == 'struct-explicit']
        if exps:
            L.append('\nExplicit structured candidates (Hessian minors etc.): ' + '; '.join(f'{r.get("name")} (d={r.get("d")}): ' + ','.join(str(x.get('status')) for x in r.get('results', [])) for r in exps) + '\n')
    errs = [r for r in rows if 'error' in r]
    if errs:
        L.append(f'\n## Errors ({len(errs)})\n')
        for r in errs[:20]:
            L.append(f'- {r.get("family")} {r.get("d")},{r.get("parD")},{r.get("k")}: {str(r.get("error"))[:200]}')
    # queue / jobs
    qlog = LANE + '/w2_queue.log'
    if os.path.exists(qlog):
        L.append('\n## Remote job queue (tail)\n```')
        L.extend(open(qlog).read().strip().splitlines()[-8:])
        L.append('```')
    notes = LANE + '/w2_notes.md'
    if os.path.exists(notes):
        L.append('\n## Notes and interpretation\n')
        L.append(open(notes).read())
    open(OUT, 'w').write('\n'.join(L) + '\n')
    print('report written', OUT, len(rows), 'rows')


if __name__ == '__main__':
    main()
