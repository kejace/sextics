#!/usr/bin/env python3
"""List wave-2 report files under /tmp/hodge/notes and print those changed since the last scan.
usage: python3 scan.py [--all] [--list]   (state in /tmp/hodge/results/_aggregator/scan-state.json)"""
import os, sys, json, time, glob
ST = '/tmp/hodge/results/_aggregator/scan-state.json'
state = json.load(open(ST)) if os.path.exists(ST) else {}
files = sorted(glob.glob('/tmp/hodge/notes/wave2-*.md') + glob.glob('/tmp/hodge/notes/lean-capacity.md') + glob.glob('/tmp/hodge/notes/engine-bugs.md'))
files = [f for f in files if not f.endswith('wave2-status.md')]
lst = '--list' in sys.argv
print(f'{len(files)} files at {time.strftime("%H:%M:%S")} local')
changed = []
for f in files:
    st = os.stat(f)
    key = f'{st.st_mtime:.0f}:{st.st_size}'
    tag = 'NEW' if f not in state else ('CHANGED' if state[f] != key else '')
    print(f'  {time.strftime("%H:%M", time.localtime(st.st_mtime))} {st.st_size:7d} {tag:8s} {os.path.basename(f)}')
    if tag or '--all' in sys.argv:
        changed.append(f)
    state[f] = key
if not lst:
    for f in changed:
        print(f'\n########## {os.path.basename(f)} ##########')
        print(open(f).read())
    json.dump(state, open(ST, 'w'))
