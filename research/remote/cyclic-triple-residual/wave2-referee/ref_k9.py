"""R5: exact spot-check of the agent's exact degree-9 E1-kernel elements (mirror/kernel_exact_d9.json):
A*K in J_15 and K not in J_9, exactly over Q, for the elements the agent sealed exactly and for the smallest-height
element of every sector (up to a height cap).  Each element is sector-pure, so one G-sector rank suffices."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-referee')
import ref_lib as L
OUT = '/tmp/hodge/results/cyclic-triple-residual/wave2-referee'
def log(rec):
    rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['script'] = 'ref_k9.py'
    with open(OUT + '/referee.log.jsonl', 'a') as f: f.write(json.dumps(rec, default=str) + '\n')
import os
src = OUT + '/mirror/kernel_exact_d9.remote.json' if os.path.exists(OUT + '/mirror/kernel_exact_d9.remote.json') else OUT + '/mirror/kernel_exact_d9.json'
_txt = open(src).read(); dat = json.loads(_txt[_txt.find('{'):])
els = dat['elements']
print('elements:', len(els), 'verified flags:', sum(1 for e in els if e['verified']), flush=True)
# the four elements the remote log reports as SEALED (exact fmpz fallback): (0,0)#?, (0,3)#0, (2,5)#0, (4,1)#0
sealed = []
for line in open(OUT + '/mirror/wave2-excl.log.remote.jsonl'):
    line = line.strip()
    if not line.startswith('{'): continue
    d = json.loads(line)
    if d.get('exp') == 'e1-kernel-element-sealed-test' and d['result'].startswith('SEALED'):
        s = d['D'].split('sector')[1]; sec = s.split('#')[0].strip(); idx = int(s.split('#')[1])
        sealed.append((json.loads(sec), idx))
print('sealed by agent (exact fallback):', sealed, flush=True)
def height(e): return max(abs(v) for v in e['Kz'].values())
pick = {}
for e in els:
    key = (tuple(e['sector']), e['index'])
    if (e['sector'], e['index']) in sealed or [e['sector'], e['index']] in [list(x) for x in [(list(s), i) for s, i in sealed]]:
        pick[key] = e
by_sec = {}
for e in els: by_sec.setdefault(tuple(e['sector']), []).append(e)
for sec, lst in by_sec.items():
    best = min(lst, key=height)
    pick[(sec, best['index'])] = best
n_ok = 0
for (sec, idx), e in sorted(pick.items()):
    t0 = time.time()
    K = {tuple(map(int, k.split(','))): v for k, v in e['Kz'].items()}
    h = height(e)
    if h > 10**60:
        print(f'  sector {sec} #{idx}: height {len(str(h))} digits, skipped', flush=True); continue
    cK = {L.chi(m) for m in K}
    inJ9 = L.in_J_exact(K)
    AK = L.pmul(L.AZ, K)
    inJ15 = L.in_J_exact(AK)
    ok = (not inJ9[0]) and inJ15[0] and len(cK) == 1 and next(iter(cK)) == tuple(sec)
    n_ok += ok
    print(f'  sector {sec} #{idx}: terms {len(K)}, height {len(str(h))} digits, char {cK}, K in J_9: {inJ9}, A*K in J_15: {inJ15}, genuine kernel element: {ok} ({time.time()-t0:.1f}s)', flush=True)
    log({'check': 'k9_exact_element', 'sector': str(sec), 'index': idx, 'terms': len(K), 'height_digits': len(str(h)), 'K_in_J9': inJ9, 'AK_in_J15': inJ15,
         'genuine': ok, 'agent_sealed_exact': (list(sec), idx) in [(list(s), i) for s, i in sealed], 'seconds': round(time.time() - t0, 1)})
print('genuine:', n_ok, 'of', len(pick), flush=True)
