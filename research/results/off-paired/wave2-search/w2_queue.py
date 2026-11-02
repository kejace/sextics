#!/usr/bin/env python3
"""Local queue runner (light): submit remote jobs from queue.txt one by one, keeping at most MAXRUN of my jobs running.
queue.txt lines:  <name>|<cpus>|<command>     (name without the search-off-paired-s2- prefix). Processed lines are moved to queue.done.
Run:  nohup python3 w2_queue.py > queue.log 2>&1 &
"""
import os, sys, time, subprocess, json
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
PREFIX = 'search-off-paired-s2-'
MAXRUN = int(os.environ.get('W2_MAXRUN', '4'))
ENGINE = '/home/kejace/devel/sextics-work'
def rjob(*args, timeout=600):
    r = subprocess.run([sys.executable, '-m', 'hodge.rjob', *args], cwd=ENGINE, capture_output=True, text=True, timeout=timeout)
    return r.stdout + r.stderr
def my_running():
    out = rjob('status', timeout=120)
    return [l.split()[1] for l in out.splitlines() if l.startswith('RUNNING ') and PREFIX in l]
def load():
    for l in out_lines():
        pass
def out_lines():
    return []
while True:
    q = [l.strip() for l in open('queue.txt')] if os.path.exists('queue.txt') else []
    q = [l for l in q if l and not l.startswith('#')]
    if not q:
        print(time.strftime('%H:%M:%S'), 'queue empty; sleeping', flush=True); time.sleep(180); continue
    try:
        running = my_running()
    except Exception as ex:
        print('status failed', ex, flush=True); time.sleep(120); continue
    if len(running) >= MAXRUN:
        print(time.strftime('%H:%M:%S'), f'{len(running)} running: {running}; waiting', flush=True); time.sleep(120); continue
    name, cpus, cmd = q[0].split('|', 2)
    print(time.strftime('%H:%M:%S'), f'submitting {name} (cpus {cpus}); running now: {running}', flush=True)
    try:
        out = rjob('submit', '--lane', 'off-paired', '--name', PREFIX + name, '--cpus', cpus, '--wait', '3600', '--', cmd, timeout=3900)
    except Exception as ex:
        print('submit failed', ex, flush=True); time.sleep(120); continue
    print(out.strip()[-400:], flush=True)
    if 'submitted' in out:
        rest = [l.strip() for l in open('queue.txt')]
        rest = [l for l in rest if l.strip() != q[0]]
        open('queue.txt', 'w').write('\n'.join(rest) + ('\n' if rest else ''))
        open('queue.done', 'a').write(time.strftime('%Y-%m-%dT%H:%M:%S') + ' ' + q[0] + '\n')
        time.sleep(30)
    else:
        time.sleep(120)
