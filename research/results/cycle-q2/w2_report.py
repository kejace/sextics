"""Summarize my wave-2 search log lines into markdown tables (per degree / family): sectors searched, restarts,
system sizes, best and median residuals, candidates, structured sweep results."""
import os as _os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "FLINT_NUM_THREADS"):
    _os.environ[_v] = "1"
import json, sys, glob, os, time
import numpy as np

LANE = '/tmp/hodge/results/cycle-q2'


def load_rows():
    rows = []
    seen = set()
    files = [os.path.join(LANE, 'wave2-search.log.jsonl')] + glob.glob('/tmp/hodge/remote/cycle-q2/wave2-search*.log.jsonl')
    for f in files:
        if not os.path.exists(f):
            continue
        for line in open(f):
            line = line.strip()
            if not line or line in seen:
                continue
            seen.add(line)
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('engine') != 2 or r.get('lens') == 'excl':
                continue
            if 'step' not in r:
                continue
            rows.append(r)
    return rows


def search_tables(rows):
    out = []
    srch = [r for r in rows if r.get('step') == 'search' and r.get('best_rel') is not None]
    # keep the last line per (tag, family, d, key, part)
    last = {}
    for r in srch:
        k = (r.get('tag'), r.get('family'), r.get('d'), r.get('key'), r.get('part', 0))
        if k not in last or r.get('restarts', 0) >= last[k].get('restarts', 0):
            last[k] = r
    groups = {}
    for r in last.values():
        fam = r['family'] if not str(r['family']).startswith(('cone4:', 'cone3:', 'mult4:', 'mult5:', 'sing2:')) else str(r['family']).split(':')[0]
        g = (r.get('q', 2), str(r.get('target', 'A'))[:30], r['d'], fam, r.get('tag'))
        groups.setdefault(g, []).append(r)
    out.append('| q | target | d | family | tag | sectors (families) | restarts total | rows x cols (max) | best rel | median of sector-best rel | #cands (<1e-10) | hist <1e-4 / <1e-6 / <1e-8 |')
    out.append('|---|---|---:|---|---|---:|---:|---|---:|---:|---:|---|')
    for g in sorted(groups, key=lambda g: (g[0], g[1], g[2], g[3])):
        rs = groups[g]
        bests = np.array([r['best_rel'] for r in rs])
        tot_restarts = sum(r.get('restarts', 0) for r in rs)
        rows_max = max(r.get('rows', 0) for r in rs); cols_max = max(r.get('cols', 0) for r in rs)
        ncand = sum(r.get('n_candidates', 0) for r in rs)
        h4 = sum(r.get('hist', {}).get('lt1e-4', 0) for r in rs); h6 = sum(r.get('hist', {}).get('lt1e-6', 0) for r in rs); h8 = sum(r.get('hist', {}).get('lt1e-8', 0) for r in rs)
        out.append(f'| {g[0]} | {g[1]} | {g[2]} | {g[3]} | {g[4]} | {len(rs)} | {tot_restarts} | {rows_max} x {cols_max} | {bests.min():.2e} | {np.median(bests):.2e} | {ncand} | {h4} / {h6} / {h8} |')
    return out


def candidates(rows):
    out = []
    for r in rows:
        if r.get('step') == 'candidate':
            out.append(f"* d={r.get('d')} key={r.get('key')} family={r.get('family')} rel={r.get('rel'):.2e} -> {r.get('recon_status')} (imag_rel {r.get('imag_rel')})")
    return out


def struct_table(rows):
    out = []
    st = [r for r in rows if r.get('step') == 'struct-sweep' and r.get('d') is not None]
    fams = {}
    for r in st:
        fams.setdefault((r['family'], r['d']), []).append(r)
    out.append('| family | d | tested | E1 fails | consistent (p1/p2) | rank range (sector system) |')
    out.append('|---|---:|---:|---:|---:|---|')
    for (fam, d), rs in sorted(fams.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        e1f = sum(1 for r in rs if r.get('e1') is False)
        cons = sum(1 for r in rs if r.get('consistent'))
        ranks = [r['rank'] for r in rs if isinstance(r.get('rank'), int) and r['rank'] > 0]
        rr = f'{min(ranks)}-{max(ranks)}' if ranks else '-'
        out.append(f'| {fam} | {d} | {len(rs)} | {e1f} | {cons} | {rr} |')
    return out


if __name__ == '__main__':
    rows = load_rows()
    print(f'{len(rows)} log rows (mine)')
    print('\n'.join(search_tables(rows)))
    print()
    print('\n'.join(struct_table(rows)))
    print()
    print('\n'.join(candidates(rows)) or 'no candidates below threshold')
    wit = [r for r in rows if r.get('step') == 'WITNESS']
    print('WITNESS lines:', len(wit))
