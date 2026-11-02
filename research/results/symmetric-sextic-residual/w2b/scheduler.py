"""Keep at most MAXRUN of this agent's remote jobs running; submit queued jobs (w2b/queue.txt lines
'<name>|<cpus>|<command>') in order as slots free up.  Submits run in a subprocess with a timeout (the runner's
ssh sometimes hangs after launching the job); the pid file on the remote is checked afterwards."""
import sys, os, time, subprocess, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import rjob
LANE = 'symmetric-sextic-residual'; W2 = f'/tmp/hodge/results/{LANE}/w2b'
MAXRUN = int(sys.argv[1]) if len(sys.argv) > 1 else 4
QUEUE, DONE, LOG = f'{W2}/queue.txt', f'{W2}/queue.done.txt', f'{W2}/scheduler.log'
def log(msg):
    with open(LOG, 'a') as f: f.write(time.strftime('%H:%M:%S ') + msg + '\n')
def alive(names):
    if not names: return []
    cmd = '; '.join(f'if [ -f /tmp/hodge/jobs/{n}.pid ] && kill -0 $(cat /tmp/hodge/jobs/{n}.pid) 2>/dev/null; then echo RUN {n}; fi' for n in names)
    r = rjob.ssh(cmd, timeout=120)
    return [l.split()[1] for l in r.stdout.splitlines() if l.startswith('RUN ')]
def started(name):
    r = rjob.ssh(f'test -f /tmp/hodge/jobs/{name}.pid && echo YES', timeout=120)
    return 'YES' in r.stdout
while True:
    lines = [l.rstrip('\n') for l in open(QUEUE)] if os.path.exists(QUEUE) else []
    pending = [l for l in lines if l.strip() and not l.startswith('#')]
    if not pending:
        log('queue empty; exiting'); break
    names = [l.split('|')[0] for l in open(DONE)] if os.path.exists(DONE) else []
    try:
        run = alive(names)
    except Exception as ex:
        log(f'status error {ex!r}'); time.sleep(120); continue
    if len(run) < MAXRUN:
        line = pending[0]; name, cpus, cmd = line.split('|', 2)
        with open(f'{W2}/submit_{name}.out', 'a') as fo:
            pr = subprocess.Popen(['python3', '-m', 'hodge.rjob', 'submit', '--lane', LANE, '--name', name, '--cpus', cpus, '--force', '--', *cmd.split()],
                                  cwd='/home/kejace/devel/sextics-work', stdout=fo, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                pr.wait(timeout=420)
            except subprocess.TimeoutExpired:
                log(f'submit timed out for {name}; killing its process group and checking the pid file')
                try: os.killpg(pr.pid, 9)
                except Exception: pass
        ok = False
        for _ in range(3):
            try:
                ok = started(name); break
            except Exception as ex:
                log(f'pid check error {ex!r}'); time.sleep(60)
        if ok:
            log(f'submitted {name} (running before: {run})')
            rest = [l for l in lines if l != line]
            open(QUEUE, 'w').write('\n'.join(rest) + ('\n' if rest else ''))
            open(DONE, 'a').write(line + '\n')
            time.sleep(120)
        else:
            log(f'submit FAILED for {name}; retry in 5 min'); time.sleep(300)
    else:
        time.sleep(180)
