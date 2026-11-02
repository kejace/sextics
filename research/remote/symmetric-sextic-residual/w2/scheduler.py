"""Keep at most MAXRUN of this lane's remote jobs running; submit queued jobs (w2/queue.txt, one shell command per
line, '#' comments; a line is '<name>|<cpus>|<command>') in order as slots free up.  Processed lines are moved to
w2/queue.done.txt.  Logs to w2/scheduler.log.  Runs until the queue is empty."""
import sys, os, time, subprocess, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import rjob
LANE = 'symmetric-sextic-residual'; W2 = f'/tmp/hodge/results/{LANE}/w2'
MAXRUN = int(sys.argv[1]) if len(sys.argv) > 1 else 4
QUEUE, DONE, LOG = f'{W2}/queue.txt', f'{W2}/queue.done.txt', f'{W2}/scheduler.log'
def log(msg):
    with open(LOG, 'a') as f: f.write(time.strftime('%H:%M:%S ') + msg + '\n')
def my_running():
    names = [l.split('|')[0] for l in open(DONE)] if os.path.exists(DONE) else []
    if not names: return 0, []
    cmd = '; '.join(f'if [ -f /tmp/hodge/jobs/{n}.pid ] && kill -0 $(cat /tmp/hodge/jobs/{n}.pid) 2>/dev/null; then echo RUN {n}; fi' for n in names)
    r = rjob.ssh(cmd)
    run = [l.split()[1] for l in r.stdout.splitlines() if l.startswith('RUN ')]
    return len(run), run
while True:
    lines = [l.rstrip('\n') for l in open(QUEUE)] if os.path.exists(QUEUE) else []
    pending = [l for l in lines if l.strip() and not l.startswith('#')]
    if not pending:
        log('queue empty; exiting'); break
    try:
        n, run = my_running()
    except Exception as ex:
        log(f'status error {ex!r}'); time.sleep(120); continue
    if n < MAXRUN:
        line = pending[0]; name, cpus, cmd = line.split('|', 2)
        try:
            rjob.submit(LANE, name, cmd, cpus=int(cpus), force=True)
            log(f'submitted {name} (running before: {run})')
        except SystemExit:
            log(f'submit refused {name}'); time.sleep(300); continue
        except Exception as ex:
            log(f'submit error {name}: {ex!r}'); time.sleep(300); continue
        rest = [l for l in lines if l != line]
        open(QUEUE, 'w').write('\n'.join(rest) + ('\n' if rest else ''))
        open(DONE, 'a').write(line + '\n')
        time.sleep(90)
    else:
        time.sleep(180)
