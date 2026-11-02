"""Wave-2 collector: pull the remote lane dir, append new batch summaries to wave2-search.log.jsonl, rewrite the report.

usage: python3 w2_collect.py [--no-pull] [--status-note TEXT]
"""
import sys, os, json, time, glob, argparse, subprocess, statistics
ap = argparse.ArgumentParser()
ap.add_argument('--no-pull', action='store_true')
ap.add_argument('--status-note', default='')
args = ap.parse_args()
L = 'eight-plane-configuration'
WD = f'/tmp/hodge/results/{L}'
MIR = f'/tmp/hodge/remote/{L}'
LOG = f'{WD}/wave2-search.log.jsonl'
STATE = f'{WD}/wave2-collect-state.json'
REPORT = f'/tmp/hodge/notes/wave2-{L}-search.md'
NOW = time.strftime('%Y-%m-%dT%H:%M:%S')

if not args.no_pull:
    r = subprocess.run(['python3', '-m', 'hodge.rjob', 'pull', '--lane', L], cwd='/home/kejace/devel/sextics-work',
                       capture_output=True, text=True, timeout=600)
    print(r.stdout.strip()[-200:])

state = json.load(open(STATE)) if os.path.exists(STATE) else {'keys': [], 'search_counts': {}}
keys = set(state['keys'])


def readj(path):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except Exception:
            pass
    return out


def emit(key, obj):
    if key in keys:
        return False
    keys.add(key)
    obj = dict(obj); obj['ts'] = NOW; obj['wave'] = 2
    with open(LOG, 'a') as f:
        f.write(json.dumps(obj, default=str) + '\n')
    return True


dims, hess, prep, fams, searches, recons, kern_pred = [], [], [], [], {}, [], None
witness = False
# --- kernel dims
for rec in readj(f'{MIR}/kernel.jsonl'):
    if rec.get('event') == 'mult_kernel':
        dims.append(rec)
        emit(f"mult_kernel-{rec['d']}-{rec['par']}", {'family': 'E1-kernel', **{k: rec[k] for k in rec if k != 'event'}})
    elif rec.get('event') == 'hilbert':
        kern_pred = rec
# --- hessian relaxation
for tag in ('hessA', 'hessB'):
    for rec in readj(f'{MIR}/{tag}.jsonl'):
        if rec.get('event') == 'hess':
            hess.append(rec)
            emit(f"hess-{rec['d']}-{rec['par']}", {'family': 'hessian-relaxation-quadrics', **{k: rec[k] for k in rec if k != 'event'}})
# --- k10 prep
for rec in readj(f'{MIR}/k10prep.jsonl'):
    prep.append(rec)
    emit(f"k10prep-{rec['d']}-{rec['par']}", {'family': 'kernel-float-basis', **{k: rec[k] for k in rec if k != 'event'}})
# --- fam10
for rec in readj(f'{MIR}/fam10.jsonl'):
    if rec.get('event') == 'family':
        fams.append(rec)
        emit(f"fam-{rec['d']}-{rec['par']}-{rec['family']}", {'family': 'structured-family-intersection', **{k: rec[k] for k in rec if k != 'event'}})
    elif rec.get('event') == 'adm':
        emit(f"adm-{rec['d']}-{rec['par']}", {'family': 'admissible-space', **{k: rec[k] for k in rec if k != 'event'}})
# --- searches
for path in sorted(glob.glob(f'{MIR}/search-*.jsonl') + glob.glob(f'{MIR}/stage1-*.jsonl')):
    tag = os.path.basename(path)[:-6]
    recs = readj(path)
    start = next((r for r in recs if r.get('event') == 'start'), None)
    rr = [r for r in recs if 'seed' in r and 'error' not in r and r.get('event') is None]
    errs = [r for r in recs if 'error' in r]
    rc = [r for r in recs if r.get('event') == 'reconstruct']
    for r in rc:
        recons.append((tag, r))
        if r.get('WITNESS'):
            witness = True
    if not rr:
        if start:
            searches[tag] = {'start': start, 'n': 0, 'errors': len(errs)}
        continue
    stage1 = any(r.get('stage1_only') for r in rr)
    key = 'r_part_rel' if stage1 else 'rel'
    vals = sorted(float(r[key]) for r in rr if r.get(key) is not None)
    rp = sorted(float(r['r_part_rel']) for r in rr if r.get('r_part_rel') is not None)
    summ = {'tag': tag, 'family': (start or {}).get('family'), 'd': (start or {}).get('d'), 'par': (start or {}).get('par'),
            'K': (start or {}).get('K'), 'dimP': (start or {}).get('dimP'), 'rows': (start or {}).get('rows'),
            'nc': (start or {}).get('nc'), 'nQ': (start or {}).get('nQ'), 'stage1_only': stage1,
            'n': len(rr), 'errors': len(errs), 'best': vals[0] if vals else None,
            'q10': vals[len(vals) // 10] if vals else None, 'median': statistics.median(vals) if vals else None,
            'worst': vals[-1] if vals else None, 'best_r_part_rel': rp[0] if rp else None,
            'n_below_1e-6': sum(1 for v in vals if v < 1e-6), 'n_below_1e-10': sum(1 for v in vals if v < 1e-10),
            'mean_seconds': round(statistics.mean(float(r.get('seconds', 0)) for r in rr), 1),
            'mean_nfev': round(statistics.mean(float(r.get('nfev', 0)) for r in rr), 1),
            'reconstructions': len(rc), 'start': start, 'done': any(r.get('event') == 'end' for r in recs)}
    searches[tag] = summ
    last = state['search_counts'].get(tag, 0)
    if len(rr) > last:
        state['search_counts'][tag] = len(rr)
        emit(f'search-{tag}-{len(rr)}', {'family': 'varpro-search-batch', **{k: v for k, v in summ.items() if k != 'start'}})
for tag, r in recons:
    emit(f"recon-{tag}-{json.dumps(r, sort_keys=True)[:200]}", {'family': 'rational-reconstruction', 'tag': tag, **r})

state['keys'] = sorted(keys)
json.dump(state, open(STATE, 'w'))

# ----------------------------------------------------------------------------- report
def fmt(v):
    if v is None: return '-'
    if isinstance(v, float): return f'{v:.3g}'
    return str(v)

lines = []
status = 'WITNESS FOUND' if witness else 'NO WITNESS'
nsearch = sum(s.get('n', 0) for s in searches.values())
lines.append(f'# Wave 2 search report — `{L}` (q = 1, A = X0^3 Y0^3)\n')
lines.append(f'STATUS: {status}. Updated {NOW}. {nsearch} search restarts collected over {len(searches)} batches; '
             f'{len(hess)} Hessian-relaxation systems, {len(dims)} new kernel dimensions, {len(fams)} family intersections, '
             f'{len(recons)} reconstruction attempts. {args.status_note}\n')
lines.append(open(f'{WD}/w2_report_static.md').read() if os.path.exists(f'{WD}/w2_report_static.md') else '')
lines.append('\n## 1. E1-admissible denominator spaces per degree (mod p = 2147483629)\n')
lines.append('| d | parity | dim S_d | dim J_d | dim R_d | kernel mod J (E1-admissible non-Jacobian) | admissible total | source |')
lines.append('|---:|:-:|---:|---:|---:|---:|---:|---|')
known = [(4, '0/1', '66/60', '0', '66/60', '0', '0', 'wave 1'), (5, '0/1', '126/126', '3/3', '123/123', '0', '3/3', 'wave 1'),
         (6, '0/1', '236/226', '18/18', '218/208', '0', '18/18', 'wave 1'), (7, '0/1', '396/396', '63/63', '333/333', '0', '63/63', 'wave 1'),
         (8, '0/1', '651/636', '168/168', '483/468', '0', '168/168', 'wave 1'), (9, '0/1', '1001/1001', '378/378', '623/623', '0', '378/378', 'wave 1'),
         (10, '0/1', '1512/1491', '750/747', '762/744', '294/261', '1044/1008', 'wave 1 (remote)')]
for row in known:
    lines.append('| ' + ' | '.join(str(x) for x in row) + ' |')
for rec in sorted(dims, key=lambda r: (r['d'], r['par'])):
    lines.append(f"| {rec['d']} | {rec['par']} | {rec['dim_J_d'] + rec['dim_R_d']} | {rec['dim_J_d']} | {rec['dim_R_d']} | {rec['kernel_mod_J']} (predicted min {rec['predicted_min_kernel']}) | {rec['kernel_dim_total']} | wave 2 ({rec['seconds']}s) |")
if kern_pred:
    pr = kern_pred['kernel_prediction']
    lines.append('\nKoszul-complex prediction of dim R_d^par and of the forced kernel (maximal rank of [A]): ' +
                 ', '.join(f"d={k.split(',')[0]} par={k.split(',')[1]}: R={v['R_d']}, min ker={v['min_kernel']}" for k, v in pr.items() if int(k.split(',')[0]) >= 10))
lines.append('\n## 2. Hessian relaxation (necessary condition A*Hess_F(w,w) in J(F)_{2d}), exact mod p\n')
if hess:
    lines.append('| d | parity | dim W | #pairs | dim R_2d^- | rank of quadrics | Macaulay ranks (deg: rank/cols) | projective variety | msolve |')
    lines.append('|---:|:-:|---:|---:|---:|---:|---|---|---|')
    for rec in sorted(hess, key=lambda r: (r['d'], r['par'])):
        mac = '; '.join(f"{k}: {v.get('rank', 'skip')}/{v.get('cols')}" for k, v in rec.get('macaulay', {}).items())
        var = 'EMPTY (mod p => over Q)' if rec.get('projective_variety_empty_modp') else ('single point' if rec.get('hilbert_1_at_degree') else 'undecided')
        ms = rec.get('msolve', {})
        mss = '-' if not ms else (f"timeout {ms.get('timeout')}" if 'timeout' in ms else f"rc={ms.get('returncode')} pure powers {ms.get('pure_power_vars')}/{rec['nW']} all={ms.get('all_pure')} ({ms.get('seconds')}s)")
        lines.append(f"| {rec['d']} | {rec['par']} | {rec['nW']} | {rec['npairs']} | {rec['dimR_2d_odd']} | {rec['rank_quadrics']} | {mac} | {var} | {mss} |")
else:
    lines.append('(pending)')
lines.append('\n## 3. Full-E2 variable-projection searches (residual landscape)\n')
lines.append('rel = ||E2 residual|| / ||Omega(D) N_0|| after optimal Koszul corrections and Q; r_part = component outside J(F)_{2d}.\n')
lines.append('| batch | family | d | par | params (free) | rows | restarts | best rel | 10% | median | worst | best r_part | <1e-6 | <1e-10 | s/restart | nfev | recon |')
lines.append('|---|---|---:|:-:|---|---:|---:|---|---|---|---|---|---:|---:|---:|---:|---:|')
for tag, s in sorted(searches.items()):
    if s.get('n', 0) == 0:
        lines.append(f"| {tag} | {s.get('start', {}).get('family')} | | | | | 0 (started) | | | | | | | | | | |")
        continue
    lines.append(f"| {tag}{' (stage 1 only)' if s['stage1_only'] else ''}{'' if s['done'] else ' (running)'} | {s['family']} | {s['d']} | {s['par']} | {s['K']} ({s['dimP']}) | {s['rows']} | {s['n']} | {fmt(s['best'])} | {fmt(s['q10'])} | {fmt(s['median'])} | {fmt(s['worst'])} | {fmt(s['best_r_part_rel'])} | {s['n_below_1e-6']} | {s['n_below_1e-10']} | {s['mean_seconds']} | {s['mean_nfev']} | {s['reconstructions']} |")
lines.append('\n## 4. Structured families at d = 10 intersected with the admissible space (mod p)\n')
if fams:
    lines.append('| d | par | family | #gens | dim fam | fam ∩ admissible | fam ∩ J | new directions | lifted |')
    lines.append('|---:|:-:|---|---:|---:|---:|---:|---:|---:|')
    for rec in sorted(fams, key=lambda r: (r['d'], r['par'], r['family'])):
        lines.append(f"| {rec['d']} | {rec['par']} | {rec['family']} | {rec['n_gens']} | {rec['dim_fam']} | {rec['fam_cap_admissible']} | {rec['fam_cap_J']} | {rec['new_directions']} | {rec.get('lifted', '-')} |")
else:
    lines.append('(pending)')
lines.append('\n## 5. Float admissible bases (d = 10)\n')
for rec in prep:
    lines.append(f"- d={rec['d']} par={rec['par']}: E1 matrix {rec['rows']}x{rec['nD']+rec['nN']}, rank {rec['rank_E1']}, nullity {rec['nullity']}, K={rec['K']}, singular-value gap {fmt(rec.get('gap_ratio'))}, E1 max err {fmt(rec['e1_max_err'])}, J-part in span rel err {fmt(rec['J_in_span_rel_err'])} ({rec['seconds']}s)")
lines.append('\n## 6. Rational reconstruction attempts\n')
if recons:
    for tag, r in recons:
        lines.append(f"- {tag}: {json.dumps(r, default=str)[:400]}")
else:
    lines.append('None yet (no restart reached the reconstruction threshold rel < 1e-9).')
open(REPORT, 'w').write('\n'.join(lines) + '\n')
print(f'report rewritten: {REPORT}; searches={len(searches)} restarts={nsearch} hess={len(hess)} dims={len(dims)} fams={len(fams)} witness={witness}')
