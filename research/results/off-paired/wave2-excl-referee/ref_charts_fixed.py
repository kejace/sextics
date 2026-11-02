#!/usr/bin/env python3
"""Referee rerun of the packet's msolve chart inputs with squares written as powers.
msolve 0.6.5 parses `t5*t5` as `t5` (verified: x*x-4 -> GB [x-4]); the packet's w2_bilinear.write_ms emits squares as `v*v`,
so the original chart verdicts were computed on the wrong systems.  This script rewrites `v*v` -> `v^2` into OUT and reruns
msolve at the same prime; a `[-1]:` verdict is then a genuine modular screen for that chart.  Logs incrementally (jsonl).
usage: python3 ref_charts_fixed.py SRC_DIR OUT_DIR THREADS TIMEOUT_S"""
import glob, re, subprocess, os, sys, time, json
src, out, threads, tmo = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
os.makedirs(out, exist_ok=True)
log = open(os.path.join(out, 'results.jsonl'), 'a')
done = set()
if os.path.exists(os.path.join(out, 'results.jsonl')):
    for line in open(os.path.join(out, 'results.jsonl')):
        try: done.add(json.loads(line)['file'])
        except Exception: pass
hist = {}
for f in sorted(glob.glob(os.path.join(src, '*.ms')), key=os.path.getsize):   # smallest charts first
    b = os.path.basename(f)
    if b in done: continue
    txt = open(f).read()
    fixed, nsub = re.subn(r'\b(t\d+|s\d+)\*\1\b', r'\1^2', txt)
    assert not re.search(r'\b(t\d+|s\d+)\*\1\b', fixed)
    dst = os.path.join(out, b); open(dst, 'w').write(fixed)
    t0 = time.time()
    try:
        subprocess.run(['msolve', '-f', dst, '-o', dst[:-3] + '.out', '-t', threads, '-v', '0'], capture_output=True, text=True, timeout=tmo)
        r = open(dst[:-3] + '.out').read().strip()[:12] if os.path.exists(dst[:-3] + '.out') else 'no-output'
    except subprocess.TimeoutExpired:
        r = 'timeout'
    rec = {'file': b, 'squares_fixed': nsub, 'verdict': r, 'seconds': round(time.time() - t0, 1), 'prime': txt.split('\n')[1].strip()}
    log.write(json.dumps(rec) + '\n'); log.flush()
    hist[r] = hist.get(r, 0) + 1
print('done', sum(hist.values()), hist)
