"""Render log.jsonl / images.jsonl / pf-*.json into markdown tables for the attack note."""
import json, glob, os
os.chdir('/tmp/hodge/results/polarized-k3-remainder')
recs = [json.loads(l) for l in open('log.jsonl') if l.strip()]
recs = [r for r in recs if r.get('engine') == 2 and 'pretest' not in r['D']]
def frame_of(r): return 'uv' if 'frame=uv' in r.get('note','') else ('xy' if 'frame=xy' in r.get('note','') else '-')
def size_of(r):
    import re
    m = re.search(r'rows=(\d+) cols=(\d+)', r.get('note',''))
    return f'{m.group(1)}x{m.group(2)}' if m else '-'
print('### Fixed-denominator solves (engine v2, p = 2^31-1; "frame" = coordinates in which D is written: uv means U=X+Y, V=X-Y)\n')
print('| family | D | d | frame | E1 | consistent | rank | nullity | system | s |')
print('|---|---|---|---|---|---|---|---|---|---|')
for r in sorted(recs, key=lambda r: (r['d'], r['family'], r['D'])):
    e1 = {True:'yes', False:'**no**', None:'-'}[r['e1']]
    cons = {True:'**YES**', False:'no', None:'skipped'}[r['consistent']]
    note = ''
    if 'ILLEGAL' in r.get('note',''): cons = 'illegal (F | D)'
    if 'mixed' in r.get('note',''): cons = 'skipped (mixed parity)'
    print(f"| {r['family']} | `{r['D']}` | {r['d']} | {frame_of(r)} | {e1} | {cons} | {r['rank']} | {r['nullity']} | {size_of(r)} | {r['seconds']} |")
if os.path.exists('images.jsonl'):
    print('\n### Served-target spaces (all Y-odd sextic targets served by D) and their sector projections in R_6 = R_f (x) R_f\n')
    print('| D | served dim (XY, incl. J) | class dim mod J | blocks (6,0),(5,1),(4,2),(3,3),(2,4),(1,5),(0,6) | K3 proj | target in span |')
    print('|---|---|---|---|---|---|')
    for l in open('images.jsonl'):
        a = json.loads(l)
        b = a.get('block_projection_dims', {})
        print(f"| `{a['label']}` | {a.get('served_dim_XY')} | {a.get('class_dim_in_R6', a.get('class_dim'))} | {[b.get(k,0) for k in ['(6,0)','(5,1)','(4,2)','(3,3)','(2,4)','(1,5)','(0,6)']]} | {a.get('K3_projection_dim')} | {a.get('target_in_span_mod_J')} |")
pfs = sorted(glob.glob('pf-*.json'))
if pfs:
    print('\n### Parity-free served spaces in the (U,V) frame (all 462 sextic targets, no parity restriction)\n')
    print('| D (U=X slots, V=Y slots) | d | rows x cols | nullity | served dim | class dim mod J | swap-anti (tau-inv) | swap-sym | blocks | K3 proj | target in span |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for f in pfs:
        J = json.load(open(f)); i = J['info']; a = J['analysis']; b = a.get('block_projection_dims', {})
        print(f"| `{J['label']}` = `{J['expr'][:60]}{'...' if len(J['expr'])>60 else ''}` | {i['d']} | {i['rows']}x{i['cols']} | {i['nullity']} | {i['served_dim']} | {a.get('class_dim_mod_J',0)} | {a.get('class_dim_tau_invariant(anti-swap)','-')} | {a.get('class_dim_tau_anti(sym-swap)','-')} | {[b.get(k,0) for k in ['(6,0)','(5,1)','(4,2)','(3,3)','(2,4)','(1,5)','(0,6)']]} | {a.get('K3_projection_dim','-')} | {a.get('target_in_span_mod_J','-')} |")
